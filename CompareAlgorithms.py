__author__ = 'Chulaka Gunasekara'
import networkx as nx
import math
import SamplingAlgorithms
import community
import matplotlib.pyplot as plt
import os
import shutil

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
            if naft > int(0.95*nbef):
                break
            else:
                print "try again"
                #G_ = nx.convert_node_labels_to_integers(G, first_label=start)
                continue
        else:
            break
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
def average(s):
    return sum(s) * 1.0 / len(s)

def standard_deviation(s):
    avg = average(s)
    variance = map(lambda x: (x - avg)**2, s)
    standard_deviation = math.sqrt(average(variance))
    return standard_deviation

i=0
readfile ='new_pies_snapshot/0.30/' + 'output-prefix.t0000'+ str(i) +'.graph'
#readfile = 'facebook.csv'
print "for "+readfile

removeDuplicateEdges('new_pies_snapshot/0.30/' + 'output-prefix.t0000'+ str(i) +'.graph', ' ', assume_one_max = True)

G = nx.read_edgelist(path=readfile, delimiter=",", nodetype=int,  create_using=nx.Graph())
start = 0
G_ = nx.convert_node_labels_to_integers(G, first_label=start)
numNodes = len(nx.nodes(G_))
G_=maximum_connected_components(G_)
partition = community.best_partition(G_)
print "For Original Network"
#
print "Num Edges - "+str(nx.number_of_edges(G_))
#print "Center of the graph "+str(nx.center(G_))
print "Eccentricity - "+str(nx.diameter(G_))
print "Radius - "+str(nx.radius(G_))
print "Modularity - "+str(community.modularity(partition, G_))
percentages = [0.1, 0.3, 0.5, 0.7]

#print size
count = 0
pos = nx.spring_layout(G_)
colors = ['#660066' ,'#eeb111' ,'#4bec13' ,'#d1d1d1' ,'#a3a3a3' ,'#c39797' ,'#a35f0c' ,'#5f0ca3' ,'#140ca3' ,'#a30c50' ,'#a30c50' ,'#0ca35f' ,'#bad8eb' ,'#ffe5a9' ,'#f5821f' ,'#00c060' ,'#00c0c0' ,'#b0e0e6' ,'#999999' ,'#ffb6c1' ,'#6897bb']
for com in set(partition.values()) :
    count = count + 1
    list_nodes = [nodes for nodes in partition.keys() if partition[nodes] == com]
    nx.draw_networkx_nodes(G_, pos, list_nodes, node_color = colors[count-1])

nx.draw_networkx_edges(G_,pos, alpha=0.5)
nx.draw(G_, pos, alpha=0.5)
plt.show()


listEdges = []
listEccentricity = []
listRadius = []
listModularity = []

for val in percentages:
    listEdges = []
    listEccentricity = []
    listRadius = []
    listModularity = []
    for j in range(0, 10):
        keptNodes = int(math.ceil(float(numNodes)*val))
        print "For Sample "+str(j+1)
        while True:
            nbunch = SamplingAlgorithms.multiRandomWalk(G_, keptNodes)
            #nbunch = SamplingAlgorithms.topDegreeSubset("All_facebook.csv.csv", keptNodes)
            H = G_.subgraph(nbunch)
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
                if naft > int(0.95*nbef):
                    break
                else:
                    print "try again"
                    #G_ = nx.convert_node_labels_to_integers(G, first_label=start)
                    continue
            else:
                break
        part = community.best_partition(H)

        edges = nx.number_of_edges(H)
        listEdges.append(edges)
        eccentricity = nx.diameter(H)
        listEccentricity.append(eccentricity)
        radius = nx.radius(H)
        listRadius.append(radius)
        modularity = community.modularity(part, H)
        listModularity.append(modularity)


        print "Num Edges - "+str(edges)
        #print "Center of the graph "+str(nx.center(G_))
        print "Eccentricity - "+str(eccentricity)
        print "Radius - "+str(radius)
        print "Modularity - "+str(modularity)
        print


    print "-------------------------------------------------------"
    print "Summary For Random Walk Sampling at percentage "+str(val)
    print "-------------------------------------------------------"

    print "Num Edges - "+str(average(listEdges))+ " - " + str(max(listEdges)) + " - " +str(min(listEdges)) + "-" + str(standard_deviation(listEdges))
    print "Eccentricity - "+str(average(listEccentricity))+ " - " + str(max(listEccentricity)) + " - " +str(min(listEccentricity)) + "-" + str(standard_deviation(listEccentricity))
    print "Radius - "+str(average(listRadius))+ " - " + str(max(listRadius)) + " - " +str(min(listRadius)) + "-" + str(standard_deviation(listRadius))
    print "Modularity - "+str(average(listModularity))+ " - " + str(max(listModularity)) + " - " +str(min(listModularity)) + "-" + str(standard_deviation(listModularity))
    print "------------------------------------------------------"













