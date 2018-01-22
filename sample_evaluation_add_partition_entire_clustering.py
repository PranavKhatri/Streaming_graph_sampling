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
import vbmod as vd

__all__=['add_path_to_graph','add_node_to_graph','add_edge_to_graph','generate_edge']

def handleArgs():
    """Handle command-line input arguments."""

    parser = argparse.ArgumentParser(description="Sample graphs.")
#    parser.add_argument("-n", "--nodes", type=int, required=True, help="the number of nodes", dest="N")
#    parser.add_argument("-s", "--start", default=1, type=int, help="the file number at which to start, inclusive", dest="start")
#    parser.add_argument("-e", "--end", default=10, type=int, help="the file number at which to end, inclusive", dest="end")
    parser.add_argument("-i", "--input", default="original_snapshot_merge/", help="the input path, defaults to 'generated_benches/'", dest="input_directory_stem")
    parser.add_argument("-o", "--output", default="powlaw_degree_benchmark_results/", help="the output path, defaults to powlaw_degree_benchmark_results", dest="output_directory_stem")
    parser.add_argument("-percentages", "--sample_percentage", nargs="+", default=[0.1, 0.3, 0.5, 0.7], help="Sample percentage of the whole graph", dest="percentages")
    parser.add_argument("-sa","--sampling_algorithm", default="", help="choose the sampling algorithm", dest="sampling_algorithm")
    parser.add_argument("-ds","--data_set",  default="", help="choose the data set", dest="ds")
    parser.add_argument("-ss", "--snapshot", default=1, type=int, help="the number of snopshots", dest="snapshot")
    global args
    args = parser.parse_args()

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

def snapshotCommunities(full_partition_prefix, part):
    """write the partition of snapshot t into separate file."""
    comms = []
    for node in part.keys():
        comms.append(str(node) + ' ' + str(part[node]) + '\n')


    commsfile = full_partition_prefix + '.comms'
    fp_commsfile = open(commsfile, 'w')
    for line in comms:
        fp_commsfile.write(line)
    fp_commsfile.close() 

    print commsfile + " completed"



#./bigclam -i:'./as20graph.txt' -c:200  -o:'jianpeng'
# def runBigclam(input_file, number_of_clusters, output_prefix):
#     """Normalized Mutual Information as defined here: https://sites.google.com/site/andrealancichinetti/mutual"""

#     result_lines = []
#     bigclam_directory="./bigclam"
#     process = subprocess.Popen(['./bigclam', '-i:', input_file, '-c:', str(number_of_clusters), '-o', output_prefix], cwd=bigclam_directory, stdout=subprocess.PIPE)
# #    output = process.communicate()[0].split()
# #    value = float(output[1])

#     print input_file + " completed"

def parseBigclamResults(output_prefix, out_file_separator):
    output_result=output_prefix+'cmtyvv.txt'
    commsfile = output_prefix + '.comms'
    # clustering file in (f_path + 'results_1/tp')

    read_file = open(output_result, 'r')
    write_file = open(commsfile, 'wb')

    cluster_number_string = '1'
    for line in read_file:
        if line[0] == '#':
            continue
        nodes = line.split()
        for node in nodes:
            write_file.write(node + out_file_separator + cluster_number_string + '\n')
        cluster_number_string = str(int(cluster_number_string) + 1)

    print('\nSuccessfully ran clustering and wrote results to file ' + output_result + 'to' + commsfile + '\n')

    read_file.close()
    write_file.close()

def runLancichNormMutInf(gold_standard_vector, partition_vector):
    """Normalized Mutual Information as defined here: https://sites.google.com/site/andrealancichinetti/mutual"""

    result_lines = []

    gold_standard_cluster_name = "file1.dat"
    partition_cluster_name = "file2.dat"

    writeLineDefinedClusterFile(gold_standard_vector, args.lnmi_directory + gold_standard_cluster_name)
    writeLineDefinedClusterFile(partition_vector, args.lnmi_directory + partition_cluster_name)

    process = subprocess.Popen(['./mutual', gold_standard_cluster_name, partition_cluster_name], cwd=args.lnmi_directory, stdout=subprocess.PIPE)
    output = process.communicate()[0].split()
    assert output[0] == 'mutual3:'
    value = float(output[1])

    os.remove(args.lnmi_directory + gold_standard_cluster_name)
    os.remove(args.lnmi_directory + partition_cluster_name)

    result_lines.append(['Lancichinetti NMI', 'Entire Graph', value])

    return result_lines

def runGMapAnalysis(relative_dotfile_path):
    """"""

    result_lines = []

    modularity = None
    conductance = None
    coverage = None

    gmap_metric_directory = args.gmap_directory + 'external/eba'
    absolute_dotfile_path = os.getcwd() + '/' + relative_dotfile_path
    subprocess.call(['./kmeans', '-action=metrics', '-o=metric_results.txt', absolute_dotfile_path], cwd = gmap_metric_directory)
    order_kmeans = ' '.join(['./kmeans', '-action=metrics', '-o=metric_results.txt', absolute_dotfile_path]) 
    gmap_output_path = gmap_metric_directory + '/metric_results.txt'
    f = open(gmap_output_path, 'r')
    for line in f:
        pieces = line.split()
        metric = pieces[0]
        value = pieces[1]
        if metric[:10] == 'Modularity':
            result_lines.append(['Modularity', 'Entire Graph', value])
        elif metric[:11] == 'Conductance':
            try:
                float(value)
                result_lines.append(['Conductance', 'Entire Graph', value])
            except ValueError:
                assert value == 'undefined'
                logfile_lines.append('Undefined conductance for ' + relative_dotfile_path)
                result_lines.append(['Conductance', 'Entire Graph', 0.0])
        elif metric[:8] == 'Coverage':
            result_lines.append(['Coverage', 'Entire Graph', value])

    f.close()
    os.remove(gmap_output_path)

    return result_lines

#./bigclam -i:'./as20graph.txt' -c:200  -o:'jianpeng'
def runBigclam(input_file, number_of_clusters, outupt_prefix):
    """Normalized Mutual Information as defined here: https://sites.google.com/site/andrealancichinetti/mutual"""

    result_lines = []
#    bigclam_directory="./bigclam/"
    bigclam_directory="."
    aaa= '-i:'+input_file
    bbb= '-c:'+str(number_of_clusters)
    ccc= '-o:'+outupt_prefix
#    order_list= './bigclam -i'+':'+input_file+' -c'+':'+str(number_of_clusters), '-o'+':'+outupt_prefix
    order_list=' '.join(['./bigclam2', aaa, bbb, ccc])
    print order_list
#    process = subprocess.Popen(['./bigclam2', '-i:', input_file, '-c:', str(number_of_clusters), '-o', outupt_prefix], cwd=bigclam_directory, stdout=subprocess.PIPE)
    process = subprocess.Popen(['./bigclam2', aaa, bbb, ccc], cwd=bigclam_directory, stdout=subprocess.PIPE)
    
#    retcode=subprocess.call(['./bigclam2', '-i:', input_file, '-c:', str(number_of_clusters), '-o', outupt_prefix])
#    print retcode    
    #process = subprocess.Popen([order_list], cwd=bigclam_directory, stdout=subprocess.PIPE)
    # output = process.communicate()[0].split()
    # value = output[1]
    output = process.communicate()
    print output
    print '\n\n\n'

    print input_file + " completed"

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

def parseLancichinettiResults(input_file, out_file):

    # clustering file in (f_path + 'results_1/tp')

    read_file = open(input_file, 'r')
    write_file = open(out_file, 'wb')

    cluster_number_string = '1'
    for line in read_file:
        if line[0] == '#':
            continue
        nodes = line.split()
        for node in nodes:
            write_file.write(node + ' ' + cluster_number_string + '\n')
        cluster_number_string = str(int(cluster_number_string) + 1)

    #print('\nSuccessfully ran transfering and wrote results to file ' + out_file + '\n')

    read_file.close()
    write_file.close() 

def writePartitionToFile(part, file_name):
    fileObject = open(file_name, 'w')  
    for node_affiliation in part:  
        fileObject.write(str(node_affiliation)+' '+ str(part[node_affiliation]))  
        fileObject.write('\n')  
    fileObject.close() 
    # for i in part: 
    #     print "dict[%s]=" % i,part[i] 
    #     fileObject.write(i+' '+ part[i])  
    #     fileObject.write('\n') 

    # print "###########items#####################" 
    # for (k,v) in  part.items:  .items(): 
    #     print "dict[%s]=" % k,v 
    #     fileObject.write(k+' '+ v)  
    #     fileObject.write('\n') 
     
    # print "###########iteritems#################" 
    # for k,v in part.iteritems(): 
    #     print "dict[%s]=" % k,v 
    #     fileObject.write(k +' '+ v)  
    #     fileObject.write('\n')              
    # print "###########iterkeys,itervalues#######" 
    # for k,v in zip(part.iterkeys(),part.itervalues()): 
    #     print "dict[%s]=" % k,v 
    #     fileObject.write(k+' '+ v)  
    #     fileObject.write('\n') 

    
if __name__ == "__main__":
    import sampling as smp
    handleArgs()
#    createPathIfNeeded(args.out_directory_stem)


    # sampling_algorithm=['old_pies','new_pies','new_isolated_pies','streamES','streamNS','original_pies']
    # data_set=['merge','grow','mixed']
    # snapshot= 9
    # percentages = [0.15,0.30,0.45,0.60]   


    sampling_algorithm= args.sampling_algorithm
    percentages=args.percentages
    snapshot=args.snapshot
    data_set=args.ds
    clustering_algorithms = ['vbmod','blondel','bigclam','bigclam_isolated_nodes']
    clustering_algorithm = clustering_algorithms[3]
    ground_truth = 1#no ground-truth
    deltas=[1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]

    # prefix_list=["birthdeath","expand","hide","mergesplit","switch"]
    # prefix= prefix_list[3]
    output_directory_stem= args.output_directory_stem #"./powlaw_degree_benchmark_results/"
    #rootdir = "../test/"


    suffix=".graph"
    

#    percentages = [0.15]
    for per in percentages:        
        results=[]
        listEdges = []
        listEccentricity = []
        listRadius = []
        listModularity = []
        
#==============================================================================
#         y_value_VNMI_micro=[]        
#         y_value_NMI_micro=[]
#         y_value_ARS_micro=[]
#         y_value_precision_micro=[]
#         y_value_recall_micro=[]
#         y_value_9_precision_micro=[]
#         y_value_9_recall_micro=[]
#         y_value_8_precision_micro=[]
#         y_value_8_recall_micro=[]
#         y_value_7_precision_micro=[]
#         y_value_7_recall_micro=[]
#         y_value_6_precision_micro=[]
#         y_value_6_recall_micro=[]
#         y_value_5_precision_micro=[]
#         y_value_5_recall_micro=[] 
#        
#         y_value_4_precision_micro=[]
#         y_value_4_recall_micro=[]
#         y_value_3_precision_micro=[]
#         y_value_3_recall_micro=[]
#         y_value_2_precision_micro=[]
#         y_value_2_recall_micro=[]
#         y_value_1_precision_micro=[]
#         y_value_1_recall_micro=[]
#         y_value_0_precision_micro=[]
#         y_value_0_recall_micro=[] 
#         
#         #stastics information        
#         y_value_sample_cluster_num_micro=[]
#         y_value_original_cluster_num_micro=[]
#         y_value_mediate_size_micro=[]
#         y_value_ANC=[]
#         y_value_NLS=[]
#==============================================================================
        
        results.append('Snapshot'+','+'JS_degree'+ ','+'KD_degree'+','+'JS_CC'+ ','+'KD_CC'+','+'JS_path'+ ','+'KD_path'+','+'num_edges' + ','+'eccentricity'+ ','+'radius'+','+'modularity'+
        ','+'1.0-precision'+ ','+'1.0-recall'  + ','+'0.9-precision'+ ','+'0.9-recall'  + ','+'0.8-precision'+ ','+'0.8-recall'+ ','+'0.7-precision'+ ','+'0.7-recall'+ ','+
'0.6-precision'+ ','+'0.6-recall'+ ','+'0.5-precision'+ ','+'0.5-recall'+ ','+'0.4-precision'+ ','+'0.4-recall'+ ','+'0.3-precision'+ ','+'0.3-recall'+ ','+'0.2-precision'+ ','+'0.2-recall'+ ','+'0.1-precision'+ ','+'0.1-recall' + ','+'0-precision'+ ','+'0-recall'+','+'mediate_size_micro'+ ','+               
'sample_cluster_num_micro'+ ','+'original_cluster_num_micro'+ ','+'ANC'+ ','+'NLS' + ','+ 'NMI'+ ','+'ARS' + ','+'VNMI'+','+ '\n')
        
        
        for i in xrange(0, snapshot): 

            t=time.time()                    
            ge=open(args.input_directory_stem + 'output-prefix.t0000' + str(i) +'.graph', 'rb')
            original_=nx.read_edgelist(ge, nodetype=int, create_using=nx.Graph())
#            num_nodes = nx.number_of_nodes(original_)           
            ge.close()
            
           # partition original             
           # G = nx.read_edgelist(path=readfile, delimiter=",", nodetype=int,  create_using=nx.Graph())
           # start = 0
           # original_ = nx.convert_node_labels_to_integers(original, first_label=start)
            original=maximum_connected_components(original_)
#           rewriteEdgelistFromZero(args.out_directory_stem + 'network_v' + str(i) + '.dat', '\t')
#           Rewrite clustering file such that node ids start from zero to maintain consistency with edgelist file node ids
#           rewriteClusteringFromZero(args.out_directory_stem + 'community_v' + str(i) + '.dat', '\t')

#            numNodes = len(nx.nodes(original))
#            print "For Original Network"
#            print "Num Edges - "+str(nx.number_of_edges(original))
#            print "Center of the graph "+str(nx.center(G_))
#            print "Eccentricity - "+str(nx.diameter(original))
#            print "Radius - "+str(nx.radius(original))

            if ground_truth==0:
                partition_ground_truth = community.best_partition(nx.Graph(original))
                snapshotCommunities(args.input_directory_stem + 'output-prefix.t0000'+ str(i), partition_ground_truth)            
            
            #tranverse the form of clustering
            #parseLancichinettiResults(args.input_directory_stem + prefix+'.t0'+ str(i) + '.comm', args.input_directory_stem + 'output-prefix.t0000'+ str(i) + '.comm')

            partition_set = se.ReadInData(args.input_directory_stem  + 'output-prefix.t0000'+ str(i) +'.comms')
            partition_ground_truth = se.partitionFromFile(args.input_directory_stem + 'output-prefix.t0000'+ str(i) +'.comms', ' ')


            number_of_gt_partition = len(partition_set) 
            print "number_of_gt_partition:" + str(number_of_gt_partition)
#            time.sleep(5)                
            print "Modularity - "+str(community.modularity(partition_ground_truth, original))
            
            #draw graph
            
#            count = 0
##            pos = nx.spring_layout(original)
#            pos = nx.spectral_layout(original)
#            colors = ['#660066' ,'#eeb111' ,'#4bec13' ,'#d1d1d1' ,'#a3a3a3' ,'#c39797' ,'#a35f0c' ,'#5f0ca3' ,'#140ca3' ,'#a30c50' ,'#a30c50' ,'#0ca35f' ,'#bad8eb' ,'#ffe5a9' ,'#f5821f' ,'#00c060' ,'#00c0c0' ,'#b0e0e6' ,'#999999' ,'#ffb6c1' ,'#6897bb']
#            for com in set(partition_ground_truth.values()) :
#                count = count + 1
#                list_nodes = [nodes for nodes in partition_ground_truth.keys() if partition_ground_truth[nodes] == com]
#                nx.draw_networkx_nodes(original, pos, list_nodes, node_color = colors[count-1])
#            
#            nx.draw_networkx_edges(original,pos, alpha=0.5)
#            nx.draw(original, pos, alpha=0.5)
#            plt.show()
                      
            #original degree disribution
#            ba_g = nx.degree_centrality(original)
#            ba_g2 = dict(Counter(ba_g.values()))
#            ba_gx,ba_gy = log_binning(ba_g2,50)

#             plt.figure(1)
#             plt.xscale('log')
#             plt.yscale('log')
# #            plt.scatter(ba_gx,ba_gy,c='r',marker='s',s=50)
#             plt.scatter(ba_g2.keys(),ba_g2.values(),c='b',marker='x')
#             plt.xlabel('Connections (normalized)')
#             plt.ylabel('Frequency')
# #            plt.xlim((1e-4,1e-1))
# #            plt.ylim((.9,1e4))
#             plt.show()
            
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
            per_format = ("%.6f" % float(per))
            full_partition_prefix = output_directory_stem + sampling_algorithm+'_'+ data_set+'_result/'+str(per_format)+'/output-prefix.t0000'+ str(i) 
            
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
            
            if clustering_algorithm == 'vbmod':            
                #run vbmod
    
                Sample_test=nx.read_edgelist(full_partition_prefix +"_test.graph")
                
                # convert networkx graph object to sparse matrix
                node_list=Sample_test.nodes()
                A=nx.to_scipy_sparse_matrix(Sample_test, nodelist=node_list)
                print(A.todense())
                
                N=A.shape[0]        # number of nodes
                Kvec=range(number_of_gt_partition,number_of_gt_partition+1) # range of K values over which to search
                
                # hyperparameters for priors
                net0={}
                net0['ap0']=N*1.
                net0['bp0']=1.
                net0['am0']=1.
                net0['bm0']=N*1.
                
                # options
                opts={}
                opts['NUM_RESTARTS']=1
                
                # run vb
                (net,net_K)=vd.learn_restart(A.tocsr(),Kvec,net0,opts)
                
                q = net['Q']
                filename =full_partition_prefix +'.comms';
                fp = open(filename, "w+");
                u = 0
                r = 0;
                max_probability = 0.0;
                max_index = 0;
                K=number_of_gt_partition
                
                for u in range(len(q)):
                    fp.write(str(node_list[u]))
                #initialization
                    max_probability = 0.0;
                    max_index = 0;
                    for r in range(K): 
                        if q[u,r] > max_probability:
                            max_probability = q[u,r]
                            max_index = r 
                #printf(" %.6f",q[u][r])
                    fp.write(" "+ str(max_index))
                    fp.write("\n")
                fp.close()
                part = se.partitionFromFile(full_partition_prefix +'.comms', ' ') 
                
            elif clustering_algorithm == 'blondel':
            #run blondel
                part = community.best_partition(nx.Graph(sample))
                snapshotCommunities(full_partition_prefix, part)
                
            elif clustering_algorithm == 'bigclam':           
            #run bigclam
            #case 1
                runBigclam(full_partition_prefix +'.graph', number_of_gt_partition, full_partition_prefix)
                parseBigclamResults(full_partition_prefix, ' ')   
                part = se.partitionFromFile(full_partition_prefix +'.comms', ' ')
            elif clustering_algorithm == 'bigclam_isolated_nodes':
            #case 2
                nx.write_edgelist(sample, full_partition_prefix +"_test.graph",delimiter='\t',data=False)
                runBigclam(full_partition_prefix +"_test.graph", number_of_gt_partition, full_partition_prefix)
                parseBigclamResults(full_partition_prefix, ' ')
                part = se.partitionFromFile(full_partition_prefix +'.comms', ' ')            


            expand_partition_with_other_components(part, sample_, sample)
            writePartitionToFile(part, full_partition_prefix +'_plus_isolated.comms')
            
            commsfile = full_partition_prefix + '_plus_isolated.comms'                       
            partition_sample = se.ReadInData(commsfile)
            first_result_lines = se.Compare3(copy.deepcopy(partition_set), partition_sample, deltas, float(per), commsfile)
            result_lines_reference = se.runComparisonGraphMetrics(part, partition_ground_truth, False) 
            
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

            p4_precision_micro = first_result_lines[12][2]
            r4_recall_micro=first_result_lines[13][2]
            p3_precision_micro=first_result_lines[14][2]
            r3_recall_micro=first_result_lines[15][2]
            p2_precision_micro=first_result_lines[16][2]
            r2_recall_micro=first_result_lines[17][2]
            p1_precision_micro=first_result_lines[18][2]
            r1_recall_micro=first_result_lines[19][2]
            p0_precision_micro=first_result_lines[20][2]
            r0_recall_micro=first_result_lines[21][2]
            

            
            mediate_size_micro=first_result_lines[22][2]                   
            sample_cluster_num_micro=first_result_lines[23][2]
            original_cluster_num_micro=first_result_lines[24][2]
            ANC=first_result_lines[25][2]
            NLS=first_result_lines[26][2] 
            
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
            
            NMI = result_lines_reference[0][2]
            ARS = result_lines_reference[1][2]
            VNMI = result_lines_reference[2][2]            
                                        
#==============================================================================
#             y_value_precision_micro.append(string.atof(first_result_lines[0][2]))  
#             y_value_recall_micro.append(string.atof(first_result_lines[1][2]))
#             y_value_9_precision_micro.append(string.atof(first_result_lines[2][2]))
#             y_value_9_recall_micro.append(string.atof(first_result_lines[3][2]))
#             y_value_8_precision_micro.append(string.atof(first_result_lines[4][2]))
#             y_value_8_recall_micro.append(string.atof(first_result_lines[5][2]))
#             y_value_7_precision_micro.append(string.atof(first_result_lines[6][2]))
#             y_value_7_recall_micro.append(string.atof(first_result_lines[7][2]))
#             y_value_6_precision_micro.append(string.atof(first_result_lines[8][2]))
#             y_value_6_recall_micro.append(string.atof(first_result_lines[9][2]))
#             y_value_5_precision_micro.append(string.atof(first_result_lines[10][2]))
#             y_value_5_recall_micro.append(string.atof(first_result_lines[11][2]))
#             y_value_mediate_size_micro.append(string.atof(first_result_lines[12][2]))                    
#             y_value_sample_cluster_num_micro.append(string.atof(first_result_lines[13][2]))
#             y_value_original_cluster_num_micro.append(string.atof(first_result_lines[14][2]))
#             y_value_ANC.append(string.atof(first_result_lines[15][2]))
#             y_value_NLS.append(string.atof(first_result_lines[16][2]))  
# 
#             y_value_NMI_micro.append(string.atof(result_lines_reference[0][2])) 
#             y_value_ARS_micro.append(string.atof(result_lines_reference[1][2])) 
#             y_value_VNMI_micro.append(string.atof(result_lines_reference[2][2])) 
#==============================================================================
            
            
            print "VNMI: - "+str(VNMI) 
#            time.sleep(5)
            print "part: - "+str(list(part))  
#            time.sleep(5)
            

#            modularity = community.modularity(part, sample)
            modularity = 0 
            listModularity.append(modularity)
            print "Modularity - "+str(modularity)
            
            print "Num Edges - "+str(edges)
            #print "Center of the graph "+str(nx.center(G_))
            print "Eccentricity - "+str(eccentricity)
            print "Radius - "+str(radius)

            
            
            
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
                                            
            results.append(str(i)+ ','+str(js_result_pies)+ ','+str(ks_result_pies)+ ','+str(cc_js_result_pies)+ ','+str(cc_ks_result_pies)+ ','
                +str(path_ks_result_pies)+ ','+str(path_js_result_pies)+ ','+str(edges) + ','+str(eccentricity)+ ','+str(radius)+ ','
                +str( modularity)+ ','+str(precision_micro)+ ','+ str(recall_micro)+ ','+ str(p9_precision_micro)+ ','+str(r9_recall_micro)
                + ','+str(p8_precision_micro)+ ','+str(r8_recall_micro)+ ','+str(p7_precision_micro)+ ','+str(r7_recall_micro)+ ','
                + str(p6_precision_micro)+ ','+str(r6_recall_micro)+ ','+str(p5_precision_micro)+ ','+str(r5_recall_micro)+ ','
                +str(p4_precision_micro)+ ','+str(r4_recall_micro)+ ','
                +str(p3_precision_micro)+ ','+str(r3_recall_micro)+ ','
                +str(p2_precision_micro)+ ','+str(r2_recall_micro)+ ','
                +str(p1_precision_micro)+ ','+str(r1_recall_micro)+ ','
                +str(p0_precision_micro)+ ','+str(r0_recall_micro)+ ','
                +str(mediate_size_micro)+ ','+ str(sample_cluster_num_micro)+ ','+str(original_cluster_num_micro)+ ','+str(ANC)+ ','+str(NLS)+ ','+ str(NMI)+ ','+str(ARS)+ ','+str(VNMI)+','+ '\n')
              
        we=open('fixed_'+data_set+'_'+clustering_algorithm+'_'+sampling_algorithm+'_summary_Divergence_sample_p_'+str(int(100*float(per))) + '.csv','wb')
        for line in results:
            we.write(line)
        we.close()
        
    print "writing completed"

        
        

        

 
   
