## Author: Scott Emmons (scott@scottemmons.com)
## Purpose: A script to generate multiple LFR benchmark graphs based on the given parameters
## Date: January 2, 2014

import argparse
import os
import shutil
import errno
import subprocess
import networkx as nx
import numpy as np

####################
# Global Variables #
####################

flag_file_name = "myflags.dat"

##################
# Main Functions #
##################

def handleArgs():
    """Handle command-line input arguments."""

    parser = argparse.ArgumentParser(description="Generate dynamic LFR benchmark graphs.")
    parser.add_argument("-n", "--nodes", type=int, required=True, help="the number of nodes", dest="N")
#    parser.add_argument("-k", "--avgdegree", default=25, type=int, help="the average degree of the nodes, defaults to 25", dest="k")
#    parser.add_argument("--maxk", "--maxdegree", type=int, required=True, help="the maximum degree of the nodes", dest="maxk")
    parser.add_argument("--mu", type=float, required=True, help="the mixing parameter", dest="mu")
#    parser.add_argument("--minc", default=50, type=int, help="the minimum community size, defaults to 50", dest="minc")
#    parser.add_argument("--maxc", type=int, required=True, help="the maximum community size", dest="maxc")
    parser.add_argument("-ds", "--dataset", default="", help="the dataset name", dest="dataset")
    parser.add_argument("-sp", "--snapshots", default=5, type=int, help="the total numbers of snapshots, inclusive", dest="snapshots")
#    parser.add_argument("-s", "--start", default=1, type=int, help="the file number at which to start, inclusive", dest="start")
#    parser.add_argument("-e", "--end", default=10, type=int, help="the file number at which to end, inclusive", dest="end")
    parser.add_argument("-b", "--benchmark", default="binary_networks/", help="the path to the installed LFR generation software", dest="bench_directory_stem")
    parser.add_argument("-o", "--output", default="generated_benches/", help="the output path, defaults to 'generated_benches/'", dest="out_directory_stem")

    global args
    args = parser.parse_args()

def deletePathIfNeeded(path):
    try:
        shutil.rmtree(path)
    except OSError as error:
        if error.errno != errno.ENOENT:
            raise
            
def createPathIfNeeded(path):
    """Credits to user 'Heikki Toivonen' on SO: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary"""
    try:
        os.makedirs(path)
    except OSError as error:
        if error.errno != errno.EEXIST:
            raise

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

    temporary_file = 'rewriteEdgelistFromZero_temporary_program_file.dat'
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

def snapshot_rewriteEdgelistFromZero(graph_file, min_id, separator):
    """"""

    temporary_file = 'rewriteEdgelistFromZero_temporary_program_file.dat'
#    assert not os.path.isfile(temporary_file)

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

    temporary_file = 'rewriteClusteringFromZero_temporary_program_file.dat'
#    assert not os.path.isfile(temporary_file)

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

def snapshot_rewriteClusteringFromZero(clustering_file, min_id, separator):
    """"""
    temporary_file = 'rewriteClusteringFromZero_temporary_program_file.dat'
#    assert not os.path.isfile(temporary_file)

    # min_id = getMinClusteringId(clustering_file, separator)
    source = open(clustering_file, 'r')
    destination = open(temporary_file, 'wb')

    for line in source:
        node_id, cluster_id = line[:-1].split(separator) #line[:-1] removes newline character from cluster_id
        node_id = str(int(node_id) - min_id)
        destination.write(node_id + separator + cluster_id + '\n')

    source.close()
    destination.close()

    shutil.move(temporary_file, clustering_file)

def generateFlagFile(file_name, out_directory_stem, N, k, maxk, mu, minc, maxc):
    """file_name: String
    out_directory_stem: String
    N: int
    mu: float"""

    to_write = ""

    to_write += "-N " + str(N) + "\n"
    to_write += "-k " + str(k) + "\n"
    to_write += "-maxk " + str(maxk) + "\n"
    to_write += "-mu " + str(mu) + "\n"
    to_write += "-t1 2\n"
    to_write += "-t2 1\n"
    to_write += "-minc " + str(minc) + "\n"
    to_write += "-maxc " + str(maxc) + "\n"
    to_write += "-on 0\n"
    to_write += "-om 0\n"

    f = open(out_directory_stem + file_name, 'w')
    f.write(to_write)

def removeDuplicateEdges(filename, separator, assume_one_max = False):
    """Given a network file separated by separator, removes edges such that the final network file_name
    contains no two edges that connect the same pair of nodes.
    Assumes node ids and cluster ids are integers.
    If assume_one_max, the function will assume that there are at most two
    edges in the original file connecting the same pair of nodes."""

    read_file = filename
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
 
#python generate_real_graph.py -n 128 --maxk 7295 --mu 0.2 --minc 50 --maxc 7295 -s 1 -e 1 -b gen_dynamic_20151102/ -o /powlaw_degree_small_snapshot_graph_for_streaming_sampling/
#python uniform_graph.py -n 128 --mu 0.45 -ds mergesplit -sp 5 -b gen_dynamic_20151102/ -o powlaw_degree_small_snapshot_graph_for_streaming_sampling/

if __name__ == "__main__":

    handleArgs()
    deletePathIfNeeded(args.out_directory_stem +'new_'+args.dataset+'_'+str(args.N)+'_u_'+ str(int(100*args.mu))+'/')  
    createPathIfNeeded(args.out_directory_stem +'new_'+args.dataset+'_'+str(args.N)+'_u_'+ str(int(100*args.mu))+'/')
    #move the snapshots and caculate the min_id of each snapshots.
    min_id_snapshot=[]
    for snapshot in range(1, args.snapshots+1):
        org_data_suffix= 'new_'+args.dataset +'_'+str(args.N)+'_u_'+ str(int(100*args.mu))+'/output-prefix.t0000'+ str(snapshot-1)+'.graph'
        org_community_suffix= 'new_'+args.dataset+'_'+str(args.N)+'_u_'+ str(int(100*args.mu))+'/output-prefix.t0000'+ str(snapshot-1)+'.dat'
        normalized_community_suffix= 'new_'+args.dataset+'_'+str(args.N)+'_u_'+ str(int(100*args.mu))+'/output-prefix.t0000'+ str(snapshot-1)+'.comms'

        shutil.copy(args.bench_directory_stem + args.dataset+'.t0'+ str(snapshot)+'.edges', args.out_directory_stem + org_data_suffix)
        shutil.copy(args.bench_directory_stem + args.dataset+'.t0'+ str(snapshot)+'.comm', args.out_directory_stem + org_community_suffix)
        
        min_id_each = getMinEdgelistId(args.out_directory_stem + org_data_suffix,' ')
        min_id_snapshot.append(min_id_each)

    min_id = np.min(np.array(min_id_snapshot))
    print("min_id:"+str(min_id_snapshot))

    #normalize the snapshots 
    for snapshot in range(1, args.snapshots+1):
        #normalize the edges 
        org_data_suffix= 'new_'+args.dataset +'_'+str(args.N)+'_u_'+ str(int(100*args.mu))+'/output-prefix.t0000'+ str(snapshot-1)+'.graph'
        normalized_community_suffix= 'new_'+args.dataset+'_'+str(args.N)+'_u_'+ str(int(100*args.mu))+'/output-prefix.t0000'+ str(snapshot-1)+'.comms'
       
        snapshot_rewriteEdgelistFromZero(args.out_directory_stem + org_data_suffix, min_id,' ')
        removeDuplicateEdges(args.out_directory_stem + org_data_suffix, ' ', assume_one_max = False)

        #rewrite the clustering results 
        parseLancichinettiResults(args.out_directory_stem + org_community_suffix, args.out_directory_stem + normalized_community_suffix)
        snapshot_rewriteClusteringFromZero( args.out_directory_stem + normalized_community_suffix, min_id, ' ')   
        
        # elif args.dataset == "expand":
        #     shutil.copy(args.bench_directory_stem + args.dataset+'.t0'+ snapshot+'.edges', args.out_directory_stem + data_suffix)
        #     shutil.copy(args.bench_directory_stem + args.dataset+'.t0'+ snapshot+'.comm', args.out_directory_stem + org_community_suffix)

        #     rewriteEdgelistFromZero(args.out_directory_stem + data_suffix, ' ')
            
        #     parseLancichinettiResults(args.out_directory_stem + org_community_suffix, args.out_directory_stem + normalized_community_suffix)
        #     rewriteClusteringFromZero(originalfile_prefix+ 'output-prefix.t0000'+ str(i) + '.comm', ' ') 

        # elif args.dataset == "birthdeath":
        #     shutil.copy(args.bench_directory_stem + args.dataset+'.t0'+ snapshot+'.edges', args.out_directory_stem + data_suffix)
        #     shutil.copy(args.bench_directory_stem + args.dataset+'.t0'+ snapshot+'.comm', args.out_directory_stem + org_community_suffix)

        #     rewriteEdgelistFromZero(args.out_directory_stem + data_suffix, ' ')
            
        #     parseLancichinettiResults(args.out_directory_stem + org_community_suffix, args.out_directory_stem + normalized_community_suffix)
        #     rewriteClusteringFromZero(originalfile_prefix+ 'output-prefix.t0000'+ str(i) + '.comm', ' ') 

        # elif args.dataset == "switch":
        #     shutil.copy(args.bench_directory_stem + args.dataset+'.t0'+ snapshot+'.edges', args.out_directory_stem + data_suffix)
        #     shutil.copy(args.bench_directory_stem + args.dataset+'.t0'+ snapshot+'.comm', args.out_directory_stem + org_community_suffix)

        #     rewriteEdgelistFromZero(args.out_directory_stem + data_suffix, ' ')
            
        #     parseLancichinettiResults(args.out_directory_stem + org_community_suffix, args.out_directory_stem + normalized_community_suffix)
        #     rewriteClusteringFromZero(originalfile_prefix+ 'output-prefix.t0000'+ str(i) + '.comm', ' ') 

        # elif args.dataset == "hide":
        #     shutil.copy(args.bench_directory_stem + args.dataset+'.t0'+ snapshot+'.edges', args.out_directory_stem + data_suffix)
        #     shutil.copy(args.bench_directory_stem + args.dataset+'.t0'+ snapshot+'.comm', args.out_directory_stem + org_community_suffix)

        #     rewriteEdgelistFromZero(args.out_directory_stem + data_suffix, ' ')
            
        #     parseLancichinettiResults(args.out_directory_stem + org_community_suffix, args.out_directory_stem + normalized_community_suffix)
        #     rewriteClusteringFromZero(originalfile_prefix+ 'output-prefix.t0000'+ str(i) + '.comm', ' ') 
        # else:
        #     print("no dataset need to copy")



    #        if args.mu >= 0.1 and args.mu <= 0.6:
        if 0:
    #            parseLancichinettiResults(args.out_directory_stem + 'intermediate_community_v' + str(i) + '.dat', args.out_directory_stem + 'community_v' + str(i) + '.dat')
             
            re=open(args.out_directory_stem + 'network_v' + str(i) + '.dat', 'rb')
            G=nx.read_edgelist(re, nodetype=int)
            re.close()             
            #G = nx.read_edgelist(path=readfile, delimiter=",", nodetype=int,  create_using=nx.Graph())
            start = 0
            # G_ = nx.convert_node_labels_to_integers(G, first_label=start)
            G_ = nx.convert_node_labels_to_integers(G, first_label=0, ordering='default', label_attribute='old_label')
            #numNodes = len(nx.nodes(G_))
            print 'number of unique nodes:',G_.number_of_nodes()
            print 'number of unique edges:',G_.number_of_edges()
            #print 'number of nodes:',sample.graph.get('number_of_nodes_repeated',0)
            #print 'number of edges:',sample.graph.get('number_of_edges_repeated',0)
            #print 'nodes',sample.nodes()
            
            
            
       # Remove duplicate edges from edgelist file and rewrite edgelist file such that node ids start from zero for compatibility with clustering program input formats
       # removeDuplicateEdges(args.out_directory_stem + 'network_v' + str(i) + '.dat', '\t', assume_one_max = False)
       # rewriteEdgelistFromZero(args.out_directory_stem + 'network_v' + str(i) + '.dat', '\t')
       # Rewrite clustering file such that node ids start from zero to maintain consistency with edgelist file node ids
       #rewriteClusteringFromZero(args.out_directory_stem + 'community_v' + str(i) + '.dat', '\t')
        
        
        
        
        
        
        
    
    
    
    
    
    
    
    
    
    