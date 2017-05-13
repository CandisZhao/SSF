from __future__ import division
from matplotlib import pyplot as plt
import networkx as nx
import numpy as np
from itertools import combinations
import math
import random
import operator
from preprocessing import *
from decimal import *
from operator import itemgetter, attrgetter
from scipy import stats
import pandas as pd
from copy import deepcopy
import os
import sys
import datetime
from compiler.ast import flatten
from scipy.special import gammaln as gamln

factorial_cache = dict()

def factorial(n = 1000):
    for i in range(n+1):
        if i == 0:
            factorial_cache[i] = 0
        else:
            factorial_cache[i] = factorial_cache[i-1] + np.log(i) 


def p_pmf(k, M, n, N):
    return np.exp(factorial_cache[n] - factorial_cache[n-k] - factorial_cache[k] \
        + factorial_cache[M-n] - factorial_cache[M-n-N+k] - factorial_cache[N-k] \
        + factorial_cache[M-N] + factorial_cache[N] - factorial_cache[M])

def expand_community(graph, subgraph, alpha):
    node_sig = np.zeros(len(degrees), dtype=float) + 1
    sub_nodes = subgraph.nodes()
    pre_nodes = set()
    for id in sub_nodes:
        pre_nodes = pre_nodes.union(adj_list[id].keys())

    pre_nodes = list(set(sub_nodes).union(pre_nodes))
    pre_duBs = []
    for id in list(pre_nodes):
        pre_duBs.append(len(set(adj_list[id].keys()) & set(sub_nodes)))
    pre_sig = []

    sum_glo = sum(degrees)
    sum_loc = sum(np.array(degrees)[sub_nodes]) 

    i = 0
    for id in list(pre_nodes):
        rho = pre_duBs[i] * sum_glo /(degrees[id] * sum_loc)
        pvalue_up = 0.0
        
        if id not in sub_nodes:
            table = [[pre_duBs[i], (sum_loc - pre_duBs[i])], [(degrees[id] - pre_duBs[i]), (sum_glo - sum_loc - 2 * degrees[id] + pre_duBs[i])]]
        else:
            table = [[pre_duBs[i], (sum_loc - degrees[id] - pre_duBs[i])], [(degrees[id] - pre_duBs[i]), (sum_glo - sum_loc - degrees[id] + pre_duBs[i])]]

        c = np.asarray(table, dtype=np.int64)
        
        n1 = c[0,0] + c[0,1]
        n2 = c[1,0] + c[1,1]
        n = c[0,0] + c[1,0]
        
        if rho >= 2.0:
            
            p0 = p_pmf(c[0,0], n1 + n2, n1, n)
            if(c[0,1] == 0 or c[1,0] == 0):
                pvalue_up = p0
            else:
                oddsritio = c[0,0] * c[1,1] / (c[0,1]*c[1,0])
                pvalue_up = oddsritio * p0 / (oddsritio - 1)
        else:
            q1 = c[0,1] * c[1,0] / (c[0,0]*c[1,1])
            if q1 < 1.0: 
                p0 = p_pmf(c[0,0], n1 + n2, n1, n)
                pvalue_up = min([p0 * (1-q1**(c[1,0]+1))/(1-q1), 1.0])
            else:    
                pvalue_up = 1.0
        pre_sig.append(pvalue_up)
        i += 1

    ## controlling fwer
    node_sig[pre_nodes] = np.array(pre_sig)
    nodes = np.array(graph.nodes())
    fit_subnodes = nodes[node_sig <= (alpha / len(degrees))]
    xx = graph.subgraph(fit_subnodes.tolist())

    ## controlling fdr
    # node_sig[pre_nodes] = pre_sig

    # temp = node_sig.argsort()
    # ranks = np.empty(len(node_sig), int)
    # ranks[temp] = np.arange(len(node_sig)) + 1

    # node_sig_bh = node_sig * len(degrees) / ranks
    # if sum(node_sig_bh <= alpha) == 0:
    #     return nx.Graph()

    # threshold = max(node_sig[node_sig_bh<=alpha])

    # nodes = np.array(graph.nodes())
    # new_subnodes = nodes[node_sig<=threshold]
    # xx = graph.subgraph(new_subnodes.tolist())
    return xx


def search_one_subgraph(graph, seed, alpha):
    sub_g = graph.subgraph(list(seed))

    pre_sub = deepcopy(sub_g)

    cnt = 0
    while True:
        exp_sub = expand_community(graph, pre_sub, alpha)
        
        if cnt > 30 or len(exp_sub.nodes()) == 0 or set(pre_sub.nodes()) == set(exp_sub.nodes()):
            break
        else:
            pre_sub = deepcopy(exp_sub)
        cnt += 1

    return exp_sub


def search_subgraph(graph, alpha ):

    factorial(graph.number_of_edges()*2)
    pre_adj_list = nx.to_dict_of_dicts(graph)

    global adj_list
    adj_list = [1] * len(pre_adj_list)
    for id, val in pre_adj_list.iteritems():
        adj_list[id] = val

    global degrees
    degrees = [len(items) for items in adj_list]
    print np.mean(degrees), "average degrees"
    print np.var(degrees), "variance degrees"

    

    visited_flags = np.zeros(len(degrees), dtype=bool)
    seed_node = degrees.index(max(degrees))
    
    final_subgraphs = []
    pre_seednode = seed_node
    cnt = 0
    while ~visited_flags[seed_node]:

        visited_flags[seed_node] = True
        seed_neighbors = graph.neighbors(seed_node)
        seed_neighbors.append(seed_node)

        temp = search_one_subgraph(graph, set(seed_neighbors), alpha)
        if temp.number_of_nodes() > 2:

            final_subgraphs.append(temp.nodes())
            visited_flags[temp.nodes()] = True

        temp_degree = list(np.array(degrees) * ~visited_flags)
        seed_node = temp_degree.index(max(temp_degree))
    return final_subgraphs



def main():
    global graph
    if len(sys.argv) < 2:
        return

    fin = open(sys.argv[1], "r")
    fout = open(sys.argv[2], "w")
    input_file = fin.readlines()
    global GP
    GP = GraphProcess()
    graph = GP.construct_graph(input_file)
    alpha = float(sys.argv[3])
    result_f = search_subgraph(graph, alpha)
    result1_f = GP.setcover(result_f)
    GP.layout(result1_f, fout)

    fin.close()
    fout.close()

    print sys.argv[1], "gold_standard/CYC2008_complex"
    os.system("python scripts/match_standalone.py -q -n" + sys.argv[1] + " gold_standard/CYC2008_complex.txt " + sys.argv[2])
    print sys.argv[1], "gold_standard/mips_3_100"
    os.system("python scripts/match_standalone.py -q -n" + sys.argv[1] + " gold_standard/mips_3_100.txt " + sys.argv[2])
    print sys.argv[1], "gold_standard/sgd"
    os.system("python scripts/match_standalone.py -q -n" + sys.argv[1] + " gold_standard/sgd.txt " + sys.argv[2])


if __name__ == '__main__':
    starttime = datetime.datetime.now()
    main()
    endtime = datetime.datetime.now()
    interval = (endtime - starttime).seconds / 3600
    print interval, "hours"