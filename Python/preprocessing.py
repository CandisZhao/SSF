import networkx as nx
import collections
import numpy as np

class GraphProcess:
    protein = {}
    rvt_protein = {}


    def protein_index(self, fin):
        cnt = 0
        tmp = set()
        # l = []
        for line in fin:
            elements = line.split()
            # if elements[0].strip() not in l:
            #     l.append(elements[0].strip())
            # if elements[1].strip() not in l:
            #     l.append(elements[1].strip())
            tmp.add(elements[0].strip())
            tmp.add(elements[1].strip())
        tmp_sort = sorted(list(tmp))

        for item in tmp_sort:
            self.protein[item] = cnt
            cnt += 1

        for k, v in self.protein.items():
            self.rvt_protein[v] = k

    def protein_index_edgelist(self, edgelist):
        cnt = 0
        tmp = set()
        for item in edgelist:
            tmp.add(item[0])
            tmp.add(item[1])
        tmp_sort = sorted(list(tmp))

        # print tmp

        for item in tmp_sort:
            self.protein[item] = cnt
            cnt += 1
        # print self.protein

        for k, v in self.protein.items():
            self.rvt_protein[v] = k

    def protein_hash(self, fin):
        protein_pair = []
        for line in fin:
            elements = line.split()
            protein_pair.append((self.protein[elements[0].strip()],self.protein[elements[1].strip()]))
        return protein_pair

    def protein_edgelist_hash(self, edgelist):
        protein_pair = []
        for elements in edgelist:
            protein_pair.append((self.protein[elements[0]],self.protein[elements[1]]))
        return protein_pair

    def protein_weighted_hash(self, fin):
        protein_weighted_pair = []
        for line in fin:
            elements = line.split()
            protein_weighted_pair.append((self.protein[elements[0].strip()], \
                self.protein[elements[1].strip()], float(elements[2].strip())))
        return protein_weighted_pair

    def construct_weighted_graph(self, fin):
        self.protein_index(fin)
        protein_weighted_pair = self.protein_weighted_hash(fin)

        graph = nx.Graph()
        graph.add_weighted_edges_from(protein_weighted_pair)
        return graph

    def construct_graph(self, fin):
        self.protein_index(fin)
        protein_pair = self.protein_hash(fin)

        graph = nx.from_edgelist(protein_pair)
        return graph

    def construct_edgelist_graph(self, edgelist):
        self.protein_index_edgelist(edgelist)
        protein_pair = self.protein_edgelist_hash(edgelist)

        graph = nx.Graph()
        graph.add_edges_from(protein_pair)
        return graph

    def hash_protein(self, result):
        #print rvt_protein
        protein_result = []
        for cluster in result:
            complex = []
            for item in cluster:
                complex.append(self.rvt_protein[item])
            protein_result.append(complex)
        return protein_result

    def flatten(self, l):
        for el in l:
            if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
                for sub in self.flatten(el):
                    yield sub
            else:
                yield el
    def setcover(self, result):
        _proteins = set(self.flatten(result))
        proteins = dict(zip(_proteins, range(len(_proteins))))
        flags = np.zeros(len(proteins), dtype=bool)
        records = np.array([np.array([proteins[label] for label in labels]) for labels in result])
        ids = []
        while not np.all(flags):
            max_val = 0
            max_id = -1
            for i, record in enumerate(records):
                tmp = np.sum(~flags[record])
                if tmp > max_val:
                    max_val = tmp
                    max_id = i
            flags[records[max_id]] = True
            ids.append(max_id)
        result_setcover = []
        for i in ids:
            result_setcover.append(result[i])
        return result_setcover

    def layout(self, result, fout):
        result_tmp = self.hash_protein(result)    # for this class construct graph
        # result_tmp = result   # for networkx construct
        for cluster in result_tmp:
            #print " ".join(cluster),"cluster"
            #print sorted([int(n) for n in cluster])
            # fout.write(" ".join([str(n) for n in cluster]) + "\n")
            fout.write(" ".join(map(str, cluster)) + "\n")
        return result_tmp










# def protein_index(fin):
#     protein = {}
#     cnt = 0
#     for line in fin:
#         elements = line.split()
#         if not protein.has_key(elements[0].strip()):
#             protein[elements[0].strip()] = cnt
#             cnt += 1
#         if not protein.has_key(elements[1].strip()):
#             protein[elements[1].strip()] = cnt
#             cnt += 1
#     return protein

# def protein_hash(fin):
#     protein = protein_index(fin)
#     #lines = fin.readlines()
#     #print len(fin), "lines"
#     protein_hashpair = []
#     for line in fin:
#         elements = line.split()
#         #print elements,"elements"
#         protein_hashpair.append((protein.get(elements[0].strip()),protein.get(elements[1].strip())))
#     return protein_hashpair

# def construct_graph(fin):
#     #print "hello"
#     protein_pair = protein_hash(fin)
#     #print len(protein_pair)
#     graph = nx.from_edgelist(protein_pair)
#     return graph

# def hash_protein(fin, result):
#     protein = protein_index(fin)
#     #print protein
#     rvt_protein = {}
#     for k, v in protein.items():
#         rvt_protein[v] = k
#     #print rvt_protein
#     protein_result = []
#     for cluster in result:
#         complex = []
#         for item in cluster:
#             complex.append(rvt_protein.get(item))
#         protein_result.append(complex)
#     return protein_result

# def flatten(l):
#     for el in l:
#         if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
#             for sub in flatten(el):
#                 yield sub
#         else:
#             yield el

# def setcover(result):
#     _proteins = set(flatten(result))
#     proteins = dict(zip(_proteins, range(len(_proteins))))
#     flags = np.zeros(len(proteins), dtype=bool)
#     records = np.array([np.array([proteins[label] for label in labels]) for labels in result])
#     ids = []
#     while not np.all(flags):
#         max_val = 0
#         max_id = -1
#         for i, record in enumerate(records):
#             tmp = np.sum(~flags[record])
#             if tmp > max_val:
#                 max_val = tmp
#                 max_id = i
#         flags[records[max_id]] = True
#         ids.append(max_id)
#     result_setcover = []
#     for i in ids:
#         result_setcover.append(result[i])
#     return result_setcover

# def layout(result, fout):
#     for cluster in result:
#         #print " ".join(cluster),"cluster"
#         fout.write(" ".join(cluster) + "\n")
