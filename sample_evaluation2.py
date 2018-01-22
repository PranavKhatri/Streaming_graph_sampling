'''
Created on Aug 30, 2013

@author: Emrah Cem{emrah.cem@utdallas.edu}
'''
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
    
#if __name__ == "__main__":
#    handleArgs()
#    createPathIfNeeded(args.out_directory_stem)
  
#    for i in xrange(args.start, args.end + 1):
        # Does seed file need to be handled here?
        ###zjp add sampling stragegy
#        sample=induced_random_vertex_sampler(G, sample_size, with_replacement=False)
#        fh=open(args.out_directory_stem + 'network_sample_v' + str(i) + '.dat','wb')
#        nx.write_edgelist(sample, fh)
#        G=nx.read_edgelist(args.out_directory_stem + 'network_v' + str(i) + '.dat', nodetype=int)
 
#    shutil.move(args.bench_directory_stem + flag_file_name, args.out_directory_stem + "flags.dat")
    
if __name__ == "__main__":

#    handleArgs()
#    createPathIfNeeded(args.out_directory_stem)

    import sampling as smp
    sampling_algorithm=['old_pies','new_pies']
    snapshot= 9
    percentages = [0.15,0.30,0.45,0.60]
    for per in percentages:        
        results=[]
        results.append('Snapshot'+ ' '+'JS_degree'+ ' '+'KD_degree'+ ' '+'JS_CC'+ ' '+'KD_CC'+ ' '+'JS_path'+ ' '+'KD_path'+ '\n')  
        for i in xrange(0, snapshot):
    #    for i in xrange(args.start, args.end + 1):     
            ###zjp add sampling stragegy      
    #        re=open(args.out_directory_stem + 'network_v' + str(i) + '.dat', 'rb')
        
            t=time.time()    
            ge=open('original_snapshot/' + 'output-prefix.t0000'+ str(i) +'.graph', 'rb')
#            ge=open('new_pies_snapshot/0.15/' + 'output-prefix.t0000'+ str(i) +'.graph', 'rb')
            original=nx.read_edgelist(ge, nodetype=int)
            ge.close()
           
            #original degree disribution
            ba_g = nx.degree_centrality(original)
            ba_g2 = dict(Counter(ba_g.values()))
            ba_gx,ba_gy = log_binning(ba_g2,50)
            
            plt.figure(1)
            plt.xscale('log')
            plt.yscale('log')
            plt.scatter(ba_gx,ba_gy,c='r',marker='s',s=50)
            plt.scatter(ba_g2.keys(),ba_g2.values(),c='b',marker='x')
    #        plt.xlim((1e-4,1e-1))
    #        plt.ylim((.9,1e4))
            plt.xlabel('Connections (normalized)')
            plt.ylabel('Frequency')
    #        plt.show()
            
            
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
            print('original_pathlength_dis',original_clusteringcoefficient_dis)
            print('\n')
            
            
            #sample degree disribution   
#            float('%.2f'%per)  round(per,2)
#            print round(per,2)
            per_format = ("%.2f" % per)
            re=open(sampling_algorithm[1]+'_snapshot/'+str(per_format)+'/output-prefix.t0000'+ str(i) +'.graph', 'rb')
            sample=nx.read_edgelist(re, nodetype=int)
            re.close() 
            print 'number of unique nodes:',sample.number_of_nodes()
            print 'number of unique edges:',sample.number_of_edges()
    #        print 'nodes',original.nodes()
            
            plt.figure(2)
            ba_c = nx.degree_centrality(sample)
            ba_c2 = dict(Counter(ba_c.values()))
            ba_x,ba_y = log_binning(ba_c2,50)
    
            plt.xscale('log')
            plt.yscale('log')
            plt.scatter(ba_x,ba_y,c='r',marker='s',s=50)
            plt.scatter(ba_c2.keys(),ba_c2.values(),c='b',marker='x')
            plt.xlim((1e-4,1e-1))
            plt.ylim((.9,1e4))
            plt.xlabel('Connections (normalized)')
            plt.ylabel('Frequency')
            #plt.show()
            
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
            results.append(str(i)+ ' '+str(js_result_pies)+ ' '+str(ks_result_pies)+ ' '+str(cc_js_result_pies)+ ' '+str(cc_ks_result_pies)+ ' '+str(path_ks_result_pies)+ ' '+str(path_js_result_pies)+ '\n')
              
        we=open(sampling_algorithm[1]+'_summary_Divergence_sample_p'+str(int(100*float(per))) + '.dat','wb')
        for line in results:
            we.write(line)
        we.close()
        
    print "writing completed"

        
        
#        sampleCommunities(sample,args.out_directory_stem + 'community_v' + str(i) + '.dat', args.out_directory_stem + 'community_sample_p'+str(int(100*float(val)))+'_v' + str(i) + '.dat','\t')
#        we=open(args.out_directory_stem + 'network_sample_p'+str(int(100*float(val)))+'_v' + str(i) + '.dat','wb')
#        nx.write_edgelist(sample, we ,data=False)
#        we.close()
 

               
        
#==============================================================================
#         #G = nx.read_edgelist(path=readfile, delimiter=",", nodetype=int,  create_using=nx.Graph())
#         start = 0
#         G_ = nx.convert_node_labels_to_integers(G, first_label=start)
#         numNodes = len(nx.nodes(G_))
#         
# #        percentages = [0.1, 0.3, 0.5, 0.7]
#         sampling_conditions = ['induced_random_edge','induced_random_vertex','induced_weighted_random_vertex','kk_path','km_path','random_path','random_vertex','random_edge','random_walk','metropolis_subgraph','metropolized_random_walk','weighted_vertex']
#         
#         #sampling strategy
#         #sampling_condition = str(sampling_conditions[10])
#         sampling_condition = args.sampling_condition
#         if not sampling_condition in sampling_conditions:
#             raise ValueError('Invalid stopping criteria, please choose one from ['+'"UNIQUE_NODES", "UNIQUE_EDGES", "NODES", "EDGES"'+']')
#         for val in args.percentages:       
#             sample_size = int(math.ceil(float(numNodes)*float(val)))        
#             t=time.time()
#             if sampling_condition=='induced_random_vertex':
#                 sample = induced_random_vertex_sampler(G, sample_size, with_replacement=False)
#             elif sampling_condition=='induced_random_edge':
#                 sample = induced_random_edge_sampler(G, sample_size, stopping_condition='UNIQUE_NODES', with_replacement=True)
#             elif sampling_condition=='induced_weighted_random_vertex':#invalid
#                 sample = sampling.sampling_algorithms.induced_weighted_random_vertex_sampler(G, sample_size, weights=None, with_replacement=True)
#             elif sampling_condition=='kk_path':#invalid                
#                 sample = kk_path_sampler(G, sample_size, K=None, vantage_points=None, stopping_condition='UNIQUE_NODES', fuzzy_select=True, include_last_path_when_exceeds=True)
#             elif sampling_condition=='km_path':#invalid
#                 sample = km_path_sampler(G, sample_size, K=None, M=None, source_nodes=None, destination_nodes=None, source_destination_nodes_can_overlap=False, stopping_condition='UNIQUE_NODES', fuzzy_select=True, include_last_path_when_exceeds=True)
#             elif sampling_condition=='random_path':
#                 sample = random_path_sampler(G, sample_size, stopping_condition='UNIQUE_NODES', include_last_path_when_exceeds=True)
#             elif sampling_condition=='random_vertex': #just nodes 
#                 sample = random_vertex_sampler(G, sample_size, with_replacement=False)
#             elif sampling_condition=='random_edge':
#                 sample = random_edge_sampler(G, sample_size, stopping_condition='UNIQUE_NODES', with_replacement=True, include_last_edge_when_exceeds=True)
#             elif sampling_condition=='random_walk':
#                 sample = random_walk_sampler(G, sample_size, initial_node=None, stopping_condition='UNIQUE_NODES', metropolized=False, excluded_initial_steps=0)
#             elif sampling_condition=='metropolis_subgraph':
#                 p=10*G.number_of_edges()*log10(G.number_of_nodes())/G.number_of_nodes() 
#                 best, div=metropolis_subgraph_sampler(G, sample_size, analytics.DivergenceMetrics.JensenShannonDivergence, smp.SimpleGraphDegree(), 1000, p, 10, 2)
#                 sample = best
#             elif sampling_condition=='metropolized_random_walk':
#                 sample = metropolized_random_walk_sampler(G, sample_size, stopping_condition='UNIQUE_NODES', excluded_initial_steps=0)
#             elif sampling_condition=='weighted_vertex':
#                 sample = weighted_vertex_sampler(G, sample_size, weights, with_replacement=True)
#         
#==============================================================================

        

 
   
