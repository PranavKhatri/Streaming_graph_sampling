# -*- coding: utf-8 -*-
"""
Created on Fri May 19 11:32:14 2017

@author: jzhang4
"""
import networkx as nx
G=nx.erdos_renyi_graph(100, 0.01)
dendrogram = nx.generate_dendrogram(G)
for level in range(len(dendrogram) - 1) :
    print("partition at level", level, "is", nx.partition_at_level(dendrogram, level)) 