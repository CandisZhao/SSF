#!/usr/bin/env python
# vim:ts=4 sw=4 sts=4 et:
"""\
%prog [options] reference_complexes predicted_complexes

Calculates matching scores between a set of reference and predicted complexes.
The input files must contain the reference and the predicted complexes, one
complex per line.
"""

from __future__ import division

import optparse
import sys
import os
from mwmatching import maxWeightMatching
from textwrap import dedent

__author__ = "Tamas Nepusz <nepusz@hal.elte.hu>"

###########################################################################

def canonical_protein_name(name):
    """Returns the canonical name of a protein by performing a few simple
    transformations on the name."""
    return name.strip().upper()

def is_numeric(x):
    """Returns whether the given string can be interpreted as a number."""
    try:
        float(x)
        return True
    except:
        return False

def matching_score(set1, set2):
    """Calculates the matching score between two sets (e.g., a cluster and a complex)
    using the approach of Bader et al, 2001"""
    return len(set1.intersection(set2))**2 / (float(len(set1)) * len(set2))

###########################################################################

def accuracy(reference, predicted):
    return (clusteringwise_sensitivity(reference, predicted) * \
            positive_predictive_value(reference, predicted)) ** 0.5

def clusteringwise_sensitivity(reference, predicted):
    num, den = 0., 0.
    for complex in reference:
        den += len(complex)
        num += max(len(complex.intersection(cluster)) for cluster in predicted)
    if den == 0.:
        return 0.
    return num / den

def clusteringwise_separation(reference, predicted):
    intersections = {}
    marginal_sums = [0.] * len(predicted), [0.] * len(reference)
    for i, cluster in enumerate(predicted):
        for j, complex in enumerate(reference):
            isect = len(cluster.intersection(complex))
            if isect > 0:
                intersections[i, j] = isect
            marginal_sums[0][i] += isect
            marginal_sums[1][j] += isect

    separations_complex = [0.] * len(reference)
    separations_cluster = [0.] * len(predicted)
    for i, cluster in enumerate(predicted):
        s = marginal_sums[0][i]
        for j, complex in enumerate(reference):
            isect = intersections.get((i, j), 0)
            if isect == 0:
                continue
            val = float(isect * isect) / (s * marginal_sums[1][j])
            separations_complex[j] += val
            separations_cluster[i] += val

    avg_sep_complex = sum(separations_complex) / len(separations_complex)
    avg_sep_cluster = sum(separations_cluster) / len(separations_cluster)
    return (avg_sep_complex * avg_sep_cluster) ** 0.5

def fraction_matched(reference, predicted, score_threshold=0.25):
    result = 0

    for id1, c1 in enumerate(reference):
        for id2, c2 in enumerate(predicted):
            score = matching_score(c1, c2)
            if score > score_threshold:
                result += 1
                break
    # print "matched = %d" %result
    # print "predicted = ", len(predicted)
    return result / len(reference)

def maximum_matching_ratio(reference, predicted, score_threshold=0.25):
    scores = {}

    n = len(reference)
    for id1, c1 in enumerate(reference):
        for id2, c2 in enumerate(predicted):
            score = matching_score(c1, c2)
            if score <= score_threshold:
                continue

            scores[id1, id2+n] = score

    input = [(v1, v2, w) for (v1, v2), w in scores.iteritems()]
    mates = maxWeightMatching(input)
    score = sum(scores[i, mate] for i, mate in enumerate(mates) if i < mate)
    return score / n

def positive_predictive_value(reference, predicted):
    num, den = 0., 0.
    for cluster in predicted:
        isects = [len(cluster.intersection(complex)) for complex in reference]
        isects.append(0.)
        num += max(isects)
        den += sum(isects)
    if den == 0.:
        return 0.
    return num / den

def precision(reference, predicted, score_threshold=0.25):
    result = 0
    for id1, c1 in enumerate(predicted):
        for id2, c2 in enumerate(reference):
            score = matching_score(c1, c2)
            if score > score_threshold:
                result += 1
                break
    # print "pre_matched = ", result
    # print "predicted= ", len(predicted)
    return result / len(predicted)

def recall(reference, predicted, score_threshold=0.25):
    result = 0.0
    for id1, c1 in enumerate(reference):
        for id2, c2 in enumerate(predicted):
            score = matching_score(c1, c2)
            if score > score_threshold:
                result += 1.0
                break
    return result / len(reference)

def f_measure(reference, predicted, score_threshold=0.25):
    p = precision(reference, predicted, score_threshold)
    r = recall(reference, predicted, score_threshold)
    return 2 * p * r / (p + r)


###########################################################################

class MatchApplication(object):
    def __init__(self):
        self.measures = dict(
                precision=precision,
                recall=recall,
                f_measure=f_measure,
                frac=fraction_matched,
                acc=accuracy,
                mmr=maximum_matching_ratio
        )
        self.parser = self.create_parser()

    def create_parser(self):
        parser = optparse.OptionParser(usage=dedent(sys.modules[__name__].__doc__).strip())

        parser.add_option("-m", "--measure", action="append", dest="measures", default=[],
                metavar="MEASURE", help="calculate the quality measure given by MEASURE. "
                "Possible values are: %s. May be given multiple times." %
                ", ".join(sorted(self.measures.keys())))
        parser.add_option("-n", "--network", metavar="FILE", dest="network",
                help="read the PPI network from FILE and assume that only these complexes "
                "were known to the clustering algorithm")
        parser.add_option("-q", "--quiet", action="store_true", dest="quiet",
                default=False, help="be quiet")

        return parser

    def log(self, msg):
        if self.options.quiet:
            return
        print >>sys.stderr, msg

    def read_complexes(self, fname, known_proteins=None, strictness=0.5,
            min_size=3, max_size=100):
        result = []
        for line in open(fname):
            ps = set(canonical_protein_name(x) for x in line.strip().split() if x)
            if known_proteins is not None:
                isect = ps.intersection(known_proteins)
                if len(isect) < max(min_size, len(ps) * strictness):
                    continue
                if len(isect) > max_size:
                    continue
                ps = isect
            result.append(ps)

        to_delete = set()
        for idx, cluster in enumerate(result):
            for idx2, cluster2 in enumerate(result):
                if idx == idx2 or idx2 in to_delete:
                    continue
                if cluster == cluster2:
                    to_delete.add(idx2)

        result = [r for i, r in enumerate(result) if i not in to_delete]
        return result

    def read_network(self, fname):
        known_proteins = set()
        for line in open(fname):
            parts = [canonical_protein_name(part) for part in line.strip().split()
                    if not is_numeric(part)]
            known_proteins.update(parts)
        return known_proteins

    def run(self):
        self.options, self.args = self.parser.parse_args()
        if len(self.args) != 2:
            self.parser.print_help()
            return 1

        if not self.options.measures:
            self.options.measures = sorted(self.measures.keys())

        if self.options.network:
            known_proteins = self.read_network(self.options.network)
            self.log("%d known proteins found in network" % len(known_proteins))
        else:
            known_proteins = None
        reference_complexes = self.read_complexes(self.args[0], known_proteins)
        predicted_complexes = self.read_complexes(self.args[1])
        self.log("%d reference complexes, %d predicted complexes" %
                (len(reference_complexes), len(predicted_complexes)))
        

        fout = open("temp.txt","w")
        for complex in reference_complexes:
            string = " ".join(complex) + "\n"
            fout.write(string)
        fout.close()

        command = "./Overlapping-NMI-master/onmi temp.txt  "+ self.args[1]
        os.system(command)

        for measure in self.options.measures:
            if measure not in self.measures:
                self.log("Ignoring unknown measure: %s" % measure)
                continue

            result = self.measures[measure](reference_complexes, predicted_complexes)
            print "%s = %.4f" %(measure,result)
            # print result

        return 0


def main():
    return MatchApplication().run()


if __name__ == "__main__":
    sys.exit(main())
