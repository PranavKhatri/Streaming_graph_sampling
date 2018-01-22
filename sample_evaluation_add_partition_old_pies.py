'''
Created on Aug 30, 2013

@author: Emrah Cem{emrah.cem@utdallas.edu}
'''
import copy
import string
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
import supervised_evaluation as se

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
            print "Taking the highest sub graph"
            nbef = len(H.nodes())
            print "Nodes before - "+str(len(H.nodes()))
            highestCompNodes = 0
            for comp in nx.connected_component_subgraphs(H):
                compNodes = len(comp.nodes())
                if compNodes > highestCompNodes:
                    highestCompNodes = compNodes
                    H = comp
            print "Nodes after - "+str(len(H.nodes()))
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
    

    
if __name__ == "__main__":

#    handleArgs()
#    createPathIfNeeded(args.out_directory_stem)

    import sampling as smp
    sampling_algorithm=['old_pies','new_pies','new_isolated_pies','streamES','streamNS','original_pies']
    sa=0
    ground_truth = 1#no ground-truth
    snapshot= 9
    percentages = [0.15,0.30,0.45,0.60]    
#    percentages = [0.60] 
    deltas=[1.0, 0.9, 0.8, 0.7, 0.6, 0.5]
    

#    percentages = [0.15]
    for per in percentages:        
        results=[]
        listEdges = []
        listEccentricity = []
        listRadius = []
        listModularity = []
        
        y_value_NMI_micro=[]
        y_value_ARS_micro=[]
        y_value_precision_micro=[]
        y_value_recall_micro=[]
        y_value_9_precision_micro=[]
        y_value_9_recall_micro=[]
        y_value_8_precision_micro=[]
        y_value_8_recall_micro=[]
        y_value_7_precision_micro=[]
        y_value_7_recall_micro=[]
        y_value_6_precision_micro=[]
        y_value_6_recall_micro=[]
        y_value_5_precision_micro=[]
        y_value_5_recall_micro=[] 
        #stastics information        
        y_value_sample_cluster_num_micro=[]
        y_value_original_cluster_num_micro=[]
        y_value_mediate_size_micro=[]
        y_value_ANC=[]
        y_value_NLS=[]
        
        results.append('Snapshot'+','+'JS_degree'+ ','+'KD_degree'+','+'JS_CC'+ ','+'KD_CC'+','+'JS_path'+ ','+'KD_path'+','+'num_edges' + ','+'eccentricity'+ ','+'radius'+','+'modularity'+
        ','+'1.0-precision'+ ','+'1.0-recall'  + ','+'0.9-precision'+ ','+'0.9-recall'  + ','+'0.8-precision'+ ','+'0.8-recall'+ ','+'0.7-precision'+ ','+'0.7-recall'+ ','+
'0.6-precision'+ ','+'0.6-recall'+ ','+'0.5-precision'+ ','+'0.5-recall'+ ','+'mediate_size_micro'+ ','+               
'sample_cluster_num_micro'+ ','+'original_cluster_num_micro'+ ','+'ANC'+ ','+'NLS' +','+ '\n')
        
        
        for i in xrange(0, snapshot):
    #    for i in xrange(args.start, args.end + 1):     
            ###zjp add sampling stragegy      
    #        re=open(args.out_directory_stem + 'network_v' + str(i) + '.dat', 'rb')
        
            t=time.time()    

#            removeDuplicateEdges('original_snapshot/' + 'output-prefix.t0000'+ str(i) +'.graph', ' ', assume_one_max = True)
            ge=open('original_snapshot_20w/' + 'output-prefix.t0000'+ str(i) +'.graph', 'rb')
#            ge=open('new_pies_snapshot/0.15/' + 'output-prefix.t0000'+ str(i) +'.graph', 'rb')
            original_=nx.read_edgelist(ge, nodetype=int, create_using=nx.Graph())
#            num_nodes = nx.number_of_nodes(original_)           
            ge.close()
            
           #partition original             
           #G = nx.read_edgelist(path=readfile, delimiter=",", nodetype=int,  create_using=nx.Graph())
           #start = 0
           #original_ = nx.convert_node_labels_to_integers(original, first_label=start)
            original=maximum_connected_components(original_)
#            rewriteEdgelistFromZero(args.out_directory_stem + 'network_v' + str(i) + '.dat', '\t')
#            # Rewrite clustering file such that node ids start from zero to maintain consistency with edgelist file node ids
#            rewriteClusteringFromZero(args.out_directory_stem + 'community_v' + str(i) + '.dat', '\t')

#            numNodes = len(nx.nodes(original))
#            print "For Original Network"
#            print "Num Edges - "+str(nx.number_of_edges(original))
            #print "Center of the graph "+str(nx.center(G_))
#            print "Eccentricity - "+str(nx.diameter(original))
#            print "Radius - "+str(nx.radius(original))
            if ground_truth==0:
                partition_ground_truth = community.best_partition(nx.Graph(original))            
            else:
                partition_set = se.ReadInData('original_snapshot_20w/' + 'output-prefix.t0000'+ str(i) +'.comms') 
                partition_ground_truth = se.partitionFromFile('original_snapshot_20w/' + 'output-prefix.t0000'+ str(i) +'.comms', ' ')
                
            print "Modularity - "+str(community.modularity(partition_ground_truth, original))
            
            #draw graph
#==============================================================================
#             count = 0
#             pos = nx.spring_layout(original)
#             colors = ['#660066' ,'#eeb111' ,'#4bec13' ,'#d1d1d1' ,'#a3a3a3' ,'#c39797' ,'#a35f0c' ,'#5f0ca3' ,'#140ca3' ,'#a30c50' ,'#a30c50' ,'#0ca35f' ,'#bad8eb' ,'#ffe5a9' ,'#f5821f' ,'#00c060' ,'#00c0c0' ,'#b0e0e6' ,'#999999' ,'#ffb6c1' ,'#6897bb']
#             for com in set(partition.values()) :
#                 count = count + 1
#                 list_nodes = [nodes for nodes in partition.keys() if partition[nodes] == com]
#                 nx.draw_networkx_nodes(original, pos, list_nodes, node_color = colors[count-1])
#             
#             nx.draw_networkx_edges(original,pos, alpha=0.5)
#             nx.draw(original, pos, alpha=0.5)
#             plt.show()
#==============================================================================
           
            #original degree disribution
            ba_g = nx.degree_centrality(original)
            ba_g2 = dict(Counter(ba_g.values()))
            ba_gx,ba_gy = log_binning(ba_g2,50)
#==============================================================================            
#             plt.figure(1)
#             plt.xscale('log')
#             plt.yscale('log')
#             plt.scatter(ba_gx,ba_gy,c='r',marker='s',s=50)
#             plt.scatter(ba_g2.keys(),ba_g2.values(),c='b',marker='x')
#             plt.xlabel('Connections (normalized)')
#             plt.ylabel('Frequency')
#             plt.xlim((1e-4,1e-1))
#             plt.ylim((.9,1e4))
#             plt.show()
#==============================================================================          
            
            original_degree=smp.SimpleGraphDegree()
            original_degree_dis =original_degree.compute_frontend_distribution(original)
            print('original_degree_dis',original_degree_dis)
            print('\n')
        
            original_clusteringcoefficient=smp.SimpleGraphClusteringCoefficient()
            original_clusteringcoefficient_dis = original_clusteringcoefficient.compute_frontend_distribution(original)
            print('original_clusteringcoefficient_dis',original_clusteringcoefficient_dis)
            print('\n')
            
            original_pathlength=smp.SimpleGraphPathLength()
            original_pathlength_dis = original_pathlength.compute_frontend_distribution(original)
            print('original_pathlength_dis',original_pathlength_dis)
            print('\n')
            
            
            #sample degree disribution   
#            float('%.2f'%per)  round(per,2)
#            print round(per,2)
            per_format = ("%.6f" % per)
            full_partition_prefix = sampling_algorithm[sa]+'_snapshot/'+str(per_format)+'/output-prefix.t0000'+ str(i) 
            re=open(full_partition_prefix +'.graph', 'rb')
            sample_=nx.read_edgelist(re, nodetype=int)
            re.close() 
            print 'number of unique nodes:',sample_.number_of_nodes()
            print 'number of unique edges:',sample_.number_of_edges()
            #print 'nodes',original.nodes()
            
            sample=maximum_connected_components(sample_)
            #partition sample
         
            edges = nx.number_of_edges(sample)
            listEdges.append(edges)
            eccentricity = nx.diameter(sample)
            listEccentricity.append(eccentricity)
            radius = nx.radius(sample)
            listRadius.append(radius)
            
            
            part = community.best_partition(nx.Graph(sample))
            
            #collect the partion of snapshot t
            comms = []
            for node in part.keys():
                comms.append(str(node) + ' ' + str(part[node]) + '\n')
            
            #write the partition of snapshot t into separate file
            commsfile = full_partition_prefix + '.comms'
            fp_commsfile = open(commsfile, 'w')
            for line in comms:
                fp_commsfile.write(line)
            fp_commsfile.close()  
            print commsfile + " completed"

            partition_sample = se.ReadInData(commsfile)
            first_result_lines = se.Compare3(copy.deepcopy(partition_set), partition_sample, deltas, per, commsfile)
            
            precision_micro=first_result_lines[0][2]  
            recall_micro=first_result_lines[1][2]
            p9_precision_micro = first_result_lines[2][2]
            r9_recall_micro=first_result_lines[3][2]
            p8_precision_micro=first_result_lines[4][2]
            r8_recall_micro=first_result_lines[5][2]
            p7_precision_micro=first_result_lines[6][2]
            r7_recall_micro=first_result_lines[7][2]
            p6_precision_micro=first_result_lines[8][2]
            r6_recall_micro=first_result_lines[9][2]
            p5_precision_micro=first_result_lines[10][2]
            r5_recall_micro=first_result_lines[11][2]
            mediate_size_micro=first_result_lines[12][2]                   
            sample_cluster_num_micro=first_result_lines[13][2]
            original_cluster_num_micro=first_result_lines[14][2]
            ANC=first_result_lines[15][2]
            NLS=first_result_lines[16][2] 
            
            print "precision_micro - "+str(precision_micro)
            print "recall_micro - "+str(recall_micro)
            print "p8_precision_micro - "+str(p8_precision_micro)
            print "r8_recall_micro - "+str(r8_recall_micro)
            print "p5_precision_micro - "+str(p5_precision_micro)
            print "r5_recall_micro - "+str(r5_recall_micro)
            print "mediate_size_micro - "+str(mediate_size_micro)
            print "sample_cluster_num_micro - "+str(sample_cluster_num_micro)
            print "ANC - "+str(ANC)
            print "NLS - "+str(NLS)
            
                                        
            y_value_precision_micro.append(string.atof(first_result_lines[0][2]))  
            y_value_recall_micro.append(string.atof(first_result_lines[1][2]))
            y_value_9_precision_micro.append(string.atof(first_result_lines[2][2]))
            y_value_9_recall_micro.append(string.atof(first_result_lines[3][2]))
            y_value_8_precision_micro.append(string.atof(first_result_lines[4][2]))
            y_value_8_recall_micro.append(string.atof(first_result_lines[5][2]))
            y_value_7_precision_micro.append(string.atof(first_result_lines[6][2]))
            y_value_7_recall_micro.append(string.atof(first_result_lines[7][2]))
            y_value_6_precision_micro.append(string.atof(first_result_lines[8][2]))
            y_value_6_recall_micro.append(string.atof(first_result_lines[9][2]))
            y_value_5_precision_micro.append(string.atof(first_result_lines[10][2]))
            y_value_5_recall_micro.append(string.atof(first_result_lines[11][2]))
            y_value_mediate_size_micro.append(string.atof(first_result_lines[12][2]))                    
            y_value_sample_cluster_num_micro.append(string.atof(first_result_lines[13][2]))
            y_value_original_cluster_num_micro.append(string.atof(first_result_lines[14][2]))
            y_value_ANC.append(string.atof(first_result_lines[15][2]))
            y_value_NLS.append(string.atof(first_result_lines[16][2]))  
            
#            partition = partitionFromFile(full_name, '\t')
#            result_lines_reference = runComparisonGraphMetrics(partition, gold_standard_first_level, False)             
            result_lines_reference = se.runComparisonGraphMetrics(part, partition_ground_truth, False) 
            y_value_NMI_micro.append(string.atof(result_lines_reference[0][2])) 
            y_value_ARS_micro.append(string.atof(result_lines_reference[1][2]))  
                      
            modularity = community.modularity(part, sample)
            listModularity.append(modularity)
            
            print "Num Edges - "+str(edges)
            #print "Center of the graph "+str(nx.center(G_))
            print "Eccentricity - "+str(eccentricity)
            print "Radius - "+str(radius)
            print "Modularity - "+str(modularity)
            
            
            
#==============================================================================
#             plt.figure(2)
#             ba_c = nx.degree_centrality(sample)
#             ba_c2 = dict(Counter(ba_c.values()))
#             ba_x,ba_y = log_binning(ba_c2,50)
#     
#             plt.xscale('log')
#             plt.yscale('log')
#             plt.scatter(ba_x,ba_y,c='r',marker='s',s=50)
#             plt.scatter(ba_c2.keys(),ba_c2.values(),c='b',marker='x')
#             plt.xlim((1e-4,1e-1))
#             plt.ylim((.9,1e4))
#             plt.xlabel('Connections (normalized)')
#             plt.ylabel('Frequency')
#             plt.show()
#==============================================================================
            
            print time.time()-t  
            print 'number of unique nodes:',sample.number_of_nodes()
            print 'number of unique edges:',sample.number_of_edges()
            print 'nodes',sample.nodes()
                            
            sample_degree=smp.SimpleGraphDegree()
            sample_degree_dis =sample_degree.compute_frontend_distribution(sample)
            print('sample_degree_dis',sample_degree_dis)
            print('\n')
    
            
            sample_clusteringcoefficient=smp.SimpleGraphClusteringCoefficient()
            sample_clusteringcoefficient_dis = sample_clusteringcoefficient.compute_frontend_distribution(sample)
            print('sample_clusteringcoefficient_dis',sample_clusteringcoefficient_dis)
            print('\n')
            
            sample_pathlength=smp.SimpleGraphPathLength()
            sample_pathlength_dis = sample_pathlength.compute_frontend_distribution(sample)
            print('sample_pathlength_dis',sample_clusteringcoefficient_dis)
            print('\n')
            
            print('original_degree_dis',original_degree_dis)
            print('sample_degree_dis',sample_degree_dis)
    
    
            js_result_pies =dm.JensenShannonDivergence.compute(original_degree_dis, sample_degree_dis)
            ks_result_pies =dm.KolmogorovSmirnovDistance.compute(original_degree_dis, sample_degree_dis)
            print js_result_pies
            print ks_result_pies            
            cc_js_result_pies =dm.JensenShannonDivergence.compute(original_clusteringcoefficient_dis, sample_clusteringcoefficient_dis)
            cc_ks_result_pies =dm.KolmogorovSmirnovDistance.compute(original_clusteringcoefficient_dis, sample_clusteringcoefficient_dis)
            print cc_js_result_pies
            print cc_ks_result_pies
            
            path_js_result_pies =dm.JensenShannonDivergence.compute(original_pathlength_dis, sample_pathlength_dis)
            path_ks_result_pies =dm.KolmogorovSmirnovDistance.compute(original_pathlength_dis, sample_pathlength_dis)
            print path_js_result_pies
            print path_ks_result_pies
                                            
            results.append(str(i)+ ','+str(js_result_pies)+ ','+str(ks_result_pies)+ ','+str(cc_js_result_pies)+ ','+str(cc_ks_result_pies)+ ','+str(path_ks_result_pies)+ ','+str(path_js_result_pies)+ ','+str(edges) + ','+str(eccentricity)+ ','+str(radius)+ ','+str( modularity)+ ','+str(precision_micro)+ ','+ str(recall_micro)+ ','+ str(p9_precision_micro)+ ','+str(r9_recall_micro)+ ','+str(p8_precision_micro)+ ','+str(r8_recall_micro)+ ','+str(p7_precision_micro)+ ','+str(r7_recall_micro)+ ','+
str(p6_precision_micro)+ ','+str(r6_recall_micro)+ ','+str(p5_precision_micro)+ ','+str(r5_recall_micro)+ ','+str(mediate_size_micro)+ ','+               
str(sample_cluster_num_micro)+ ','+str(original_cluster_num_micro)+ ','+str(ANC)+ ','+str(NLS) +','+ '\n')
              
        we=open(sampling_algorithm[sa]+'_summary_Divergence_sample_p'+str(int(100*float(per))) + '.csv','wb')
        for line in results:
            we.write(line)
        we.close()
        
    print "writing completed"

        
        

        

 
   
