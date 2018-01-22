'''
This script is used for:
1.generate 'complete' corresponding .comms file for every snapshot under ....
2.generate 'complete' corresponding .tcomms file for every snapshot with the same sample size
'''
from __future__ import division

import argparse
import os
import shutil
import errno
import subprocess
import networkx as nx
import random
import time
import math
from sampling.sampling_algorithms import *
import analytics
import networkx as nx
from collections import Counter
from math import *
import matplotlib.pyplot as plt
import numpy as np
import DivergenceMetrics2 as dm
import community
import igraph



__all__=['add_path_to_graph','add_node_to_graph','add_edge_to_graph','generate_edge']

def drop_zeros(a_list):
    return [i for i in a_list if i>0]

def log_binning(counter_dict,bin_count=35):

    max_x = log10(max(counter_dict.keys()))
    max_y = log10(max(counter_dict.values()))
    max_base = max([max_x,max_y])

    min_x = log10(min(drop_zeros(counter_dict.keys())))

    bins = np.logspace(min_x,max_base,num=bin_count)

    # Based off of: http://stackoverflow.com/questions/6163334/binning-data-in-python-with-scipy-numpy
    bin_means_y = (np.histogram(counter_dict.keys(),bins,weights=counter_dict.values())[0] / np.histogram(counter_dict.keys(),bins)[0])
    bin_means_x = (np.histogram(counter_dict.keys(),bins,weights=counter_dict.keys())[0] / np.histogram(counter_dict.keys(),bins)[0])

    return bin_means_x,bin_means_y
    
def handleArgs():
    """Handle command-line input arguments."""

    parser = argparse.ArgumentParser(description="Sample graphs.")
    parser.add_argument("-n", "--nodes", type=int, required=True, help="the number of nodes", dest="N")
    parser.add_argument("-s", "--start", default=1, type=int, help="the file number at which to start, inclusive", dest="start")
    parser.add_argument("-e", "--end", default=10, type=int, help="the file number at which to end, inclusive", dest="end")
    parser.add_argument("-o", "--output", default="generated_benches/", help="the output path, defaults to 'generated_benches/'", dest="out_directory_stem")
    parser.add_argument("-percentages", "--sample_percentage", nargs="+", default=[0.1, 0.3, 0.5, 0.7], help="Sample percentage of the whole graph", dest="percentages")
    parser.add_argument("-sac","--sampling_condition", default="", help="choose the sampling algorithm", dest="sampling_condition")

    global args
    args = parser.parse_args()


def generate_edge(G, with_replacement):
    if with_replacement:
        while True:
            yield random.choice(G.edges())
    else:
        edge_list=random.sample(G.edges(),G.number_of_edges())#this will shuffle edges
        for e in edge_list:
            yield e
         
def add_path_to_graph(G,path):
    if len(path)==1:
        add_node_to_graph(G,path[0])
    else:
        G.add_path(path)
        for n in path:
            G.node[n]['times_selected']=G.node[n].get('times_selected',0)+1
        u=path[0]
        for v in path[1:]:
            G.edge[u][v]['times_selected']=G.edge[u][v].get('times_selected',0)+1
            u=v
        G.graph['number_of_nodes_repeated']=G.graph.get('number_of_nodes_repeated',0)+len(path)
        G.graph['number_of_edges_repeated']=G.graph.get('number_of_edges_repeated',0)+len(path)-1
        
def add_node_to_graph(G,n):
    G.add_node(n)
    G.node[n]['times_selected']=G.node[n].get('times_selected',0)+1
    G.graph['number_of_nodes_repeated']=G.graph.get('number_of_nodes_repeated',0)+1

def add_edge_to_graph(G,e, add_nodes=True):
    G.add_edge(*e)
    G.edge[e[0]][e[1]]['times_selected']=G.edge[e[0]][e[1]].get('times_selected',0)+1
    if add_nodes:
        add_node_to_graph(G, e[0])
        add_node_to_graph(G, e[1])
    G.graph['number_of_edges_repeated']=G.graph.get('number_of_edges_repeated',0)+1

def createPathIfNeeded(path):
    """Credits to user 'Heikki Toivonen' on SO: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary"""
    try:
        os.makedirs(path)
    except OSError as error:
        if error.errno != errno.EEXIST:
            raise    

def maximum_connected_components(H):
    """creat maximum_connected_components"""
        
    while True:
        connectedComponents = nx.number_connected_components(H)
         
        print "Num components: "+str(connectedComponents)
        if connectedComponents > 1:
            #print "Taking the highest sub graph"
            nbef = len(H.nodes())
            #print "Nodes before - "+str(len(H.nodes()))
            highestCompNodes = 0
            for comp in nx.connected_component_subgraphs(H):
                compNodes = len(comp.nodes())
                if compNodes > highestCompNodes:
                    highestCompNodes = compNodes
                    H = comp
            #print "Nodes after - "+str(len(H.nodes()))
            naft = len(H.nodes())
            if naft > int(nbef/connectedComponents):
                return H
                break

            else:
                print "try again"
                #G_ = nx.convert_node_labels_to_integers(G, first_label=start)
                continue
        else:
            return H
            break

###add by zhukaijie#####
def expand_partition_with_other_components(part, graph, main_components):
    max_cluster_id = -1
    for i in part:
        if (part[i] > max_cluster_id):
            max_cluster_id = part[i]
            
    connectedComponents = nx.number_connected_components(graph)
    #print "connectedComponents = %d" %connectedComponents
    
    current_cluster_id = max_cluster_id + 1
    if connectedComponents > 1:  
        for comp in nx.connected_component_subgraphs(graph):
            compNodesList = comp.nodes()
            
            maincompNodesList = main_components.nodes()
            if( compNodesList[0] in  maincompNodesList): #indicating this comp is the main comp
                continue;
                #print "this comp is the main one!"
            else:
                #print compNodesList
                for node in compNodesList:
                    part[node] = current_cluster_id
                    
            current_cluster_id = current_cluster_id + 1
                        
###add by zhukaijie#####
def add_other_components_partition(comm_list, graph, main_components):
    
    max_cluster_id = -1
    for i in comm_list:
        if (i > max_cluster_id):
            max_cluster_id = i
    #print "max_cluster_id = %d" %max_cluster_id
        
    connectedComponents = nx.number_connected_components(graph)
    #print "connectedComponents = %d" %connectedComponents
    
    current_cluster_id = max_cluster_id + 1
    if connectedComponents > 1:  
        for comp in nx.connected_component_subgraphs(graph):
            compNodesList = comp.nodes()
            print "len = %d" %(len(compNodesList))
            maincompNodesList = main_components.nodes()
            if( compNodesList[0] in  maincompNodesList): #indicating this comp is the main comp
                continue;
                #print "this comp is the main one!"
            else:
                for node in compNodesList:
                  if (node in maincompNodesList):
                      print "That's strange!"
                  comm_list[node] = current_cluster_id
            current_cluster_id = current_cluster_id + 1
            
    #print "current_cluster_id = %d" %current_cluster_id
    
def removeDuplicateEdges(filename, separator, assume_one_max = False):
    """Given a network file separated by separator, removes edges such that the final network file_name
    contains no two edges that connect the same pair of nodes.
    Assumes node ids and cluster ids are integers.
    If assume_one_max, the function will assume that there are at most two
    edges in the original file connecting the same pair of nodes."""

    read_file = filename
    os.remove("temporary_function_execution.dat")
    write_file = "temporary_function_execution.dat"
#    assert not os.path.isfile(write_file)

    with open(read_file, 'r') as read_f:
        with open(write_file, 'wb') as write_f:
            redundant_edges = {}
            empty_set = set() 
            for line in read_f:
                source, destination = line.split(separator)
                source = int(source)
                destination = int(destination.rstrip()) # remove newline character and trailing spaces
                if not destination in redundant_edges.get(source, empty_set):
                    write_f.write(str(source) + separator + str(destination) + '\n')
                    redundant_edges[destination] = redundant_edges.get(destination, empty_set)
                    redundant_edges[destination].add(source)
                    empty_set = set() # reverse mutation due to previous line
                elif assume_one_max:
                    redundant_edges[source].remove(destination)
    
    #shutil.move(read_file, read_file + 'a')
    shutil.move(write_file, read_file)
    
def sampleCommunities(sample, clustering_file, write_file, separator):
    """Given a network file separated by separator, removes edges such that the final network file_name
    contains no two edges that connect the same pair of nodes.
    Assumes node ids and cluster ids are integers.
    If assume_one_max, the function will assume that there are at most two
    edges in the original file connecting the same pair of nodes."""
    sample_nodes=set(sample.nodes())
    read_file = clustering_file
    #write_file ="sample_"+read_file
    #assert not os.path.isfile(write_file)

    with open(read_file, 'r') as read_f:
        with open(write_file, 'wb') as write_f:
            for line in read_f:
                node_id, cluster_id = line[:-1].split(separator)
                node_id = int(node_id)
                if  node_id in sample_nodes:
                    write_f.write(str(node_id) + separator + str(cluster_id) + '\n')                 
    
    #shutil.move(write_file, read_file)  
    
def buildCommunitieslist(part, total_nodes):
    
    pos = 0
    comm_list = []
    
    i = 0
    while (i < total_nodes):
        comm_list.append(-1)
        i = i + 1
        
    for node in part:                    
        comm_list[node] = part[node]
    
    return comm_list
    
def buildFixedlist(comm_list):
    
    fixed_list = []
    for i in comm_list:
        if (i == -1):
            fixed_list.append(0)
        else:
            fixed_list.append(1) 
       
    return fixed_list


def getMinEdgelistId(edgelist_file, separator):
    """"""

    with open(edgelist_file, 'r') as f:
        source_id, destination_id = f.readline().split(separator)
        destination_id = destination_id[:-1] #remove newline character
        min_id = min(int(source_id), int(destination_id))
        for line in f:
            source_id, destination_id = line[:-1].split(separator) #line[:-1] removes newline character from destination_id 
            min_id = min(int(source_id), int(destination_id), min_id)

    return min_id

def rewriteEdgelistFromZero(graph_file, separator):
    """"""

    temporary_file = 'temporary_program_file.dat'
#    assert not os.path.isfile(temporary_file)

    min_id = getMinEdgelistId(graph_file, separator)
    source = open(graph_file, 'r')
    destination = open(temporary_file, 'wb')

    for line in source:
        source_id, destination_id = line[:-1].split(separator) #line[:-1] removes newline character from destination_id
        source_id = str(int(source_id) - min_id)
        destination_id = str(int(destination_id) - min_id)
        destination.write(source_id + separator + destination_id + '\n')

    source.close()
    destination.close()

    shutil.move(temporary_file, graph_file)

def getMinClusteringId(clustering_file, separator):
    """"""
    
    with open(clustering_file, 'r') as f:
        min_id = int(f.readline().split(separator)[0])
        for line in f:
            node_id = int(line.split(separator)[0])
            min_id = min(node_id, min_id)

    return min_id

def rewriteClusteringFromZero(clustering_file, separator):
    """"""

    temporary_file = 'temporary_program_file.dat'
    assert not os.path.isfile(temporary_file)

    min_id = getMinClusteringId(clustering_file, separator)
    source = open(clustering_file, 'r')
    destination = open(temporary_file, 'wb')

    for line in source:
        node_id, cluster_id = line[:-1].split(separator) #line[:-1] removes newline character from cluster_id
        node_id = str(int(node_id) - min_id)
        destination.write(node_id + separator + cluster_id + '\n')

    source.close()
    destination.close()

    shutil.move(temporary_file, clustering_file)
                    
def verifyCommunitiesliat(comm_list):
    
    commfile = './' + snapshot_directory[j] + '/' + sample_size[k] + '/output-prefix.t0000'+ str(i) + '.comms'
    fp_commfile = open(commfile, 'r')   
    print commfile
    
    for line in fp_commfile:
        line = line.strip('\n')
        node, community = line.split(' ')
        
        if(comm_list[int(node)] != int(community)):
            print "comm_list is not well genrated"
            print line
            print node + ' ' + str(comm_list[int(node)])
            time.sleep(10)
        
    fp_commfile.close() 
    
def kCorePreProcess(comm_list, original_snapshot_graph, subgraph):

    unlabeled_nodes = 0
    unlabeled_list = []
    
    subgraph_nodes_list = subgraph.nodes()
    
    for i in range(1,len(comm_list)):

        if (comm_list[i] == -1):
            
            unlabeled_nodes = unlabeled_nodes + 1
            
            #calculate neighbor proportion
            neighbors_in_original_graph = original_snapshot_graph.neighbors(i)
            neighbors_in_subgraph = []
            for node in neighbors_in_original_graph:
                if node in subgraph_nodes_list:
                    neighbors_in_subgraph.append(node)
            proportion = len(neighbors_in_subgraph) / len(neighbors_in_original_graph)
        
            #generate triples
            element = (i, proportion,  neighbors_in_original_graph)
            unlabeled_list.append(element)
          
    sorted_unlabeled_list =  sorted(unlabeled_list, key=lambda item:item[1], reverse = True )       

    for element in sorted_unlabeled_list:    
#    for element in unlabeled_list:
        node = element[0]
        proportion = element[1]
        neighbors_in_original_graph = element[2]
                
        #generate_neighbors_list
        neighbors_in_subgraph = []
        for neighbor in neighbors_in_original_graph:
            if (comm_list[neighbor] != -1):
                neighbors_in_subgraph.append(neighbor)
        
        #count label of neighbors
        labels_dic = {}
        for neighbor in neighbors_in_subgraph:
            label = comm_list[i]
            if label in labels_dic:
                labels_dic[label] = labels_dic[label] + 1
            else:
                labels_dic[label] = 1
        if len(labels_dic) !=0:
            #select the most popular one as label of unlabeled node 
            labels_list = sorted(labels_dic.iteritems(), key=lambda item:item[1], reverse = True)
            selected_element = labels_list[0]
            selected_label = selected_element[0]
            comm_list[node] = selected_label
        
        
        
def file_save_t_communities(filename,complete_tcomms):
#    fp_complete_tcomms = open(outputfile_prefix + '/snapshot_entire.tcomms', 'w')
    fp_complete_tcomms = open(filename, 'w')
    list_map= [complete_tcomms[mykey] for mykey in sorted(complete_tcomms.keys())]      

    for edge_connects in list_map:
#    for filename in complete_tcomms.keys():
        for line in edge_connects:
            fp_complete_tcomms.write(line)
    fp_complete_tcomms.close()
    
    
    
if __name__ == "__main__":

    snapshot_directory = ['new_pies'];
    dataset=['mergesplit_1000']
    sample_size = ['0.200000']
    
    total_directories = len(snapshot_directory)

    total_samplesizes = len(sample_size)
    total_snapshots = 5
    total_nodes = 1000
    i = 0
    complete_tcomms = {}
    org_complete_tcomms = {}        
    while ( i < total_snapshots ):
        
        originalfile_prefix = './powlaw_degree_small_snapshot_graph_for_streaming_sampling/'+'new_'+ dataset[0] +'_u_20/'
#        rewriteEdgelistFromZero(originalfile_prefix+ 'output-prefix.t0000' + str(i) + '.graph', ' ')
#        rewriteClusteringFromZero(originalfile_prefix+ 'output-prefix.t0000'+ str(i) + '.comms', ' ')        
        
        originalfile = originalfile_prefix + 'output-prefix.t0000' + str(i) + '.graph'
        print "original_snapshot" + originalfile
        fp_originalfile = open(originalfile, 'r')        
        original_snapshot_graph =nx.read_edgelist(fp_originalfile, nodetype=int)
        fp_originalfile.close() 
        igr = igraph.Graph(n = original_snapshot_graph.number_of_nodes(), edges = nx.convert_node_labels_to_integers(original_snapshot_graph).edges())        
#        igr = igraph.Graph(n = original_snapshot_graph.number_of_nodes(), edges = nx.convert_node_labels_to_integers(original_snapshot_graph, first_label=0).edges())


        org_complete_tcommsfile =  originalfile 
        if (org_complete_tcommsfile not in org_complete_tcomms.keys()) :
            org_complete_tcomms[org_complete_tcommsfile] = []
            
        ground_truthfile= originalfile_prefix + 'output-prefix.t0000' + str(i) + '.comms'            
        fp_complete_commsfile = open(ground_truthfile, 'r')
        for line in fp_complete_commsfile.readlines():
            org_complete_tcomms[org_complete_tcommsfile].append(str(i) + ' ' + line)
        fp_complete_commsfile.close()
        
        
        j = 0
        
        while ( j < total_directories ):           
            k = 0           
            while ( k < total_samplesizes ):
                
                outputfile_prefix = './powlaw_degree_benchmark_results/'+ snapshot_directory[j] +'_'+ dataset[0] +'_result/' + sample_size[k]              
                edgefile_prefix = outputfile_prefix + '/snapshot_t'+ str(i)                 
#                rewriteEdgelistFromZero(edgefile_prefix + '.graph', '\t')
#                rewriteClusteringFromZero(edgefile_prefix + '.comms', ' ')                  
                
                edgefile = edgefile_prefix + '.graph'

                print "pies_snapshotfile" + edgefile
                
                fp_edgefile = open(edgefile, 'r')
                sample_=nx.read_edgelist(fp_edgefile, nodetype=int)
                fp_edgefile.close() 
                print edgefile
                
                #partition sample
                sample=maximum_connected_components(sample_)
                part = community.best_partition(nx.Graph(sample))
                
                #expand partition
                expand_partition_with_other_components(part, sample_, sample)
                    
                #generate comm_list for label propagation algorithm
                comm_list = buildCommunitieslist(part, total_nodes)
                
                #complete comm_list
                #add_other_components_partition(comm_list, sample_, sample)

                fixed_list = buildFixedlist(comm_list)
                
                #pre-process here
                kCorePreProcess(comm_list, original_snapshot_graph, sample_)
                #################
                
                #generate complete comm_lsit
                complete_part = igr.community_label_propagation(initial = comm_list, fixed = fixed_list).membership
                ####add more strategies                  
             
                
                complete_tcommsfile =  edgefile_prefix 
                if (complete_tcommsfile not in complete_tcomms.keys()) :
                    complete_tcomms[complete_tcommsfile] = []
             
                complete_commsfile =  edgefile_prefix+ '_entire.comms'
                fp_complete_commsfile = open(complete_commsfile, 'w')
                
                l = 0
                while ( l < total_nodes):
                    fp_complete_commsfile.write( str(l) + ' ' + str(complete_part[l]) + '\n')
                    complete_tcomms[complete_tcommsfile].append( str(i) + ' ' + str(l) + ' ' + str(complete_part[l]) + '\n')
                    l = l + 1
                    
                fp_complete_commsfile.close()
                
                k = k + 1
            
            j = j + 1
            
        i = i + 1
        
    file_save_t_communities(originalfile_prefix + '/snapshot_entire.tcomms',org_complete_tcomms)
    file_save_t_communities(outputfile_prefix + '/snapshot_entire.tcomms',complete_tcomms)

#============================================================================== 
#     file_save_t_communities(filename)
#     fp_complete_tcomms = open(outputfile_prefix + '/snapshot_entire.tcomms', 'w')
#     list_map= [complete_tcomms[mykey] for mykey in sorted(complete_tcomms.keys())]      
# 
#     for edge_connects in list_map:
# #    for filename in complete_tcomms.keys():
#         for line in edge_connects:
#             fp_complete_tcomms.write(line)
#     fp_complete_tcomms.close()
#==============================================================================
 