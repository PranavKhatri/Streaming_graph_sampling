from __future__ import division
'''
Created on 25 Nov 2015

@author: Administrator
'''
'''
Created on 24 Nov 2015

@author: Administrator
'''

import numpy as np
import argparse
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import adjusted_rand_score
import matplotlib.pyplot as plt
import subprocess

import pylab as pl
import time
import string
import copy
import os


def Compare4(s1, s2, deltas, sample):
    result_lines = []
    precision = 0.0
    recall = 0.0

    sample_block_num= len(s2)
    original_block_num = len(s1)

    for delta in deltas:
        res1 = set()       
        res2 = set()
        res2.clear()
        res1.clear()
        c1={}
        c2={}
        c1=copy.deepcopy(s1)
        c2=copy.deepcopy(s2)
#        print('c1 length is %d\n' %len(c1))
#        print('c2 length is %d\n' %len(c2))
#        time.sleep(5)
        
        for k2, v2 in c2.items():
            for k1, v1 in c1.items():
                if len(v2)>3: 
                    if Coverage(v1, v2):
                        res1.add(k1)
                    if DeltaCoverage(v1, v2, float(delta), sample):
                        res2.add(k1)
                        c1.pop(k1)
#                        c2.pop(k2)

        delta_precison = len(res2) * 1.0 / sample_block_num
        delta_recall = len(res2) * 1.0 / original_block_num
        delta_composite = 2.0*delta_precison*delta_recall* 1.0 / (delta_precison+delta_recall)
        print (len(res1), len(res2))
        print ('---------------------------------')
        print ('Delta Precision with delta (' + str(delta) + '): ' + str(delta_precison))
        print ('Delta Recall with delta (' + str(delta) + '): ' + str(delta_recall))
        print ('Delta F-measure with delta (' + str(delta) + '): ' + str(delta_composite))      
        result_lines.append([str(delta)+'_precison', 'Entire Graph', delta_precison])
        result_lines.append([str(delta)+'_recall', 'Entire Graph', delta_recall])
        result_lines.append([str(delta)+'_composite', 'Entire Graph', delta_composite])
    
    precision = len(res1) * 1.0 / sample_block_num
    recall = len(res1) * 1.0 / original_block_num
    print ('Precision: ' + str(precision))
    print ('Recall: ' + str(recall))
    result_lines.append(['Precison', 'Entire Graph', precision])
    result_lines.append(['Recall', 'Entire Graph', recall])
    return result_lines
    
def handleArgs():
    """Handles the command-line input arguments, placing them in the global Namespace variable 'args'."""

    parser = argparse.ArgumentParser(description="Automates the command line calls necessary to execute the full clustering analysis script workflow")
    parser.add_argument("--gpre", required=True, help="the stem for the path and filename of the graph files, before the file number", dest="graph_file_prefix")       
    parser.add_argument("--spre", default="", help="the stem for the path and filename of the 'gold_standard' clutsering files, before the file number", dest="gold_standard_file_prefix")
    parser.add_argument("--cnames", nargs="+", default=[], type=str.lower, required=True, help="the names of the clustering methods that will be evaluated, to be used in the naming of output files", dest="clustering_file_names")
    parser.add_argument("--cpre", nargs="+", default=[], help="the stem for the path and filename of the to-be-evaluated clustering files, before the file number", dest="clustering_file_prefixes")
    parser.add_argument("--csuf", nargs="+", default=[], help="the ending to the filename of the to-be-evaluated clustering files, after the file number, including the file extension; must be either a list matching the length of --cpre, or one value that is universal to all in --cpre", dest="clustering_file_suffixes") 
    parser.add_argument("--cnum", type=int, default=1, help="the number of clusterings that exist for each graph", dest="clusterings_per_graph")
    parser.add_argument("-o", "--out", help="the directory to which to write the program output files, defaults to 'metric_results/'", dest="output_path")   
    parser.add_argument("-delta","--deltavalue", nargs="+", default="0.5 0.6 0.7 0.8 0.9", help="the delta parameter(the default value is 0.9)", dest="delta")
    parser.add_argument("-sr", "--samplerate", default=0.7, type=float, help="the sample percentage of the entire graph", dest="samplerate")
   
    parser.add_argument("-b", "--benchmark", default=os.getcwd()+"/generated_benches_u1_10_u2_20_p_50_condition_metropolis_subgraph/n_2000/", help="the path to the process folder. Defaults to the current working directory + '/binary_networks/'", dest="bench_directory_stem")
    parser.add_argument("-gt", "--groundtruth", nargs="+", default=[], help="ground-truth file", dest="ground_truth_file")
#    parser.add_argument("--cpre", nargs="+", default=[], help="the stem for the path and filename of the to-be-evaluated clustering files, before the file number", dest="clustering_file_prefixes")
    global args
    args = parser.parse_args()

def appendLines(to_add, add_to):
    """Append lines to_add, a list, to add_to, a list."""
    
    for line in to_add:
        add_to.append(line)

def partitionFromFile(partition_file, partition_file_separator):
    """Create a partition object from the given file.
    Return dictionary object assigning nodes to clusters based on
    partition_file, a string.
    Returns dictionary object assigning nodes to clusters."""

    partition = {}

    f = open(partition_file, 'r')
    for line in f:
        node, cluster = line.split(partition_file_separator)
        partition[int(node)] = int(cluster.rstrip())
    f.close()

    return partition

def runScikitNormMutInf(gold_standard_vector, partition_vector):
    """Normalized mutual information as defined here: http://scikit-learn.org/stable/modules/generated/sklearn.metrics.normalized_mutual_info_score.html#sklearn.metrics.normalized_mutual_info_score""" 

    result_lines = []

    value = normalized_mutual_info_score(gold_standard_vector, partition_vector)

    result_lines.append(['Scikit-learn NMI', 'Entire Graph', value])

    return result_lines

def runAdjRandScr(gold_standard_vector, partition_vector):
    """Adjusted Rand Score as defined here: http://scikit-learn.org/stable/modules/generated/sklearn.metrics.adjusted_rand_score.html#sklearn.metrics.adjusted_rand_score"""

    result_lines = []

    value = adjusted_rand_score(gold_standard_vector, partition_vector)

    result_lines.append(['Adjusted Rand Score', 'Entire Graph', value])

    return result_lines
    
def generateCorrespondingVectors(partition_1, partition_2):
    """"""

    vector_1 = []
    vector_2 = []

    for key_node in partition_1.keys():
        if partition_2.has_key(key_node):
            vector_1.append(partition_1[key_node])
            vector_2.append(partition_2[key_node])

    return vector_1, vector_2
def writeLineDefinedClusterFile(vector, to_write_path):
    """From vector of cluster definitions which is of the form
    node i assigned to value of vector[i], write the clutser assignments
    to a file assigning the nodes on a line to the same cluster."""

    cluster_to_node = {}

    for i in xrange(len(vector)):
        try:
            cluster_to_node[vector[i]].append(i + 1)
        except KeyError:
            cluster_to_node[vector[i]] = [i + 1]

    with open(to_write_path, 'wb') as f:
        for cluster in cluster_to_node.keys():
            node_list = cluster_to_node[cluster]
            first_node = True
            for node in node_list:
                if first_node:
                    f.write(str(node))
                    first_node = False
                else:
                    f.write(' ' + str(node))
            f.write('\n')

   
def runLancichNormMutInf(gold_standard_vector, partition_vector):
    """Normalized Mutual Information as defined here: https://sites.google.com/site/andrealancichinetti/mutual"""

    result_lines = []

    gold_standard_cluster_name = "file1.dat"
    partition_cluster_name = "file2.dat"
    lnmi_directory="mutual3/"

    writeLineDefinedClusterFile(gold_standard_vector, lnmi_directory + gold_standard_cluster_name)
    writeLineDefinedClusterFile(partition_vector, lnmi_directory + partition_cluster_name)

    process = subprocess.Popen(['./mutual', gold_standard_cluster_name, partition_cluster_name], cwd=lnmi_directory, stdout=subprocess.PIPE)
    output = process.communicate()[0].split()
    assert output[0] == 'mutual3:'
    value = float(output[1])

    os.remove(lnmi_directory + gold_standard_cluster_name)
    os.remove(lnmi_directory + partition_cluster_name)

    result_lines.append(['Lancichinetti NMI', 'Entire Graph', value])

    return result_lines


def runComparisonGraphMetrics(partition, gold_standard, is_directed):
    """"""

    result_lines = []


    gold_standard_vector, partition_vector = generateCorrespondingVectors(gold_standard, partition)

    # Metrics from scikit-learn

    scikit_norm_mut_inf_lines = runScikitNormMutInf(gold_standard_vector, partition_vector)
    appendLines(scikit_norm_mut_inf_lines, result_lines)

    adj_rand_scr_lines = runAdjRandScr(gold_standard_vector, partition_vector)
    appendLines(adj_rand_scr_lines, result_lines)

    # Lancichinetti's NMI measure

    lancich_norm_mut_inf_lines = runLancichNormMutInf(gold_standard_vector, partition_vector)
    appendLines(lancich_norm_mut_inf_lines, result_lines)

    return result_lines

def NMI_score(ground_file, result_file):
    """
    Calculate the NMI score for the given cluster assignments and ground truth classes
    
    :param clusters: the cluster assignments array
    :type clusters: numpy.array
    
    :param classes: the ground truth classes
    :type classes: numpy.array
    
    :returns: the purity score
    :rtype: float
    """
    result_lines = []
    NMI =[]
    purity =[]
    comm = np.loadtxt(ground_file,dtype='int32') 
    myresult = np.loadtxt(result_file,dtype='int32')
    myresult_sorted= np.array( sorted(myresult, key=lambda myresult : myresult[0]))

    myresult_label= np.array(myresult_sorted[:,1])
    
    comm_label= np.array(comm[:,1])
    """
    print (len(myresult))
    comm_sorted =np.array( sorted(comm, key=lambda myresult_sorted : myresult_sorted[0]))

    comm_selected= np.array(comm_sorted[myresult[:,0]-1,:])

    truelabel= np.array(comm_selected[:,1])
    print (len(truelabel))
    print (purity_score(myresult_label, truelabel))
    #    print (normalized_mutual_info_score([0, 0, 0, 0], [0, 1, 2, 3]))
    print ("resut is ",normalized_mutual_info_score(comm_label.tolist(), myresult_label.tolist()))
   
    """
    NMI = normalized_mutual_info_score(comm_label.tolist(), myresult_label.tolist())
    ARS = adjusted_rand_score(comm_label.tolist(), myresult_label.tolist())
    purity = purity_score(myresult_label, comm_label)

    result_lines.append(['NMI', 'Entire Graph', NMI])
    result_lines.append(['Adjusted Rand Score', 'Entire Graph', ARS])
    result_lines.append(['purity', 'Entire Graph', purity])
    return result_lines

def purity_score(clusters, classes):
    """
    Calculate the purity score for the given cluster assignments and ground truth classes
    
    :param clusters: the cluster assignments array
    :type clusters: numpy.array
    
    :param classes: the ground truth classes
    :type classes: numpy.array
    
    :returns: the purity score
    :rtype: float
    """
    
    A = np.c_[(clusters,classes)]
    
    n_accurate = 0.

    for j in np.unique(A[:,0]):
        z = A[A[:,0] == j, 1]
        x = np.argmax(np.bincount(z))
        n_accurate += len(z[z == x])

    return n_accurate / A.shape[0]

def ReadInData(filename):
    clusters = {}
    fin = open(filename, 'r')
    for line in fin.readlines():
        res = line.strip().split()
        if clusters.has_key(int(res[1])):
            ids = clusters.get(int(res[1]))
            ids.add(int(res[0]))
            clusters[int(res[1])] = ids
        else:
            ids = set()
            ids.add(int(res[0]))
            clusters[int(res[1])] = ids
    fin.close()
    return clusters

def Coverage(s1, s2):
    if s1.issuperset(s2) and len(s2)>3:
        return True
    else:
        return False
'''
def DeltaCoverage(s1, s2, threshold, sample):
#    delta = len(s1 & s2) * 1.0 / len(s2 | s1)
#    print("intersection of two set is %d" % (len(s1 & s2)))
#    print("union of two set is %d" % (len(s1 | s2)))

    delta = len(s1 & s2) * 1.0 / len(s2)
    delta2 = len(s1 & s2) * 1.0 / (len(s1) * sample)
#    delta3 = len(s1 & s2)
#    if delta>= threshold and delta3>= 5:#threshold:       
    if delta>= threshold and delta2>= 0.5:#threshold:    
        print("\n\n\nnew original set:delta=%0.2f, delta2=%0.2f\n"%(delta, delta2))
#        print(list(s1))
#        print("sample set \n")
#        print(list(s2))   
#        time.sleep(5)
        return True
        
#    if delta>= threshold or delta2>= 0.5:#threshold:
#        return True 
    else:
        return False



def Compare3(s1, s2, deltas, sample):
    result_lines = []
    precision = 0.0
    recall = 0.0

    sample_block_num= 0
    original_block_num = 0
    block_size=[]
    block_size_median=0
    
    for kk1, vv1 in s1.items():
        if len(vv1)>2:
            original_block_num=original_block_num+1
            
    
    for kk2, vv2 in s2.items():
        if len(vv2)>2:
            sample_block_num=sample_block_num+1
            block_size.append(len(vv2))
    block_size_median = np.median(block_size)

    for delta in deltas:
        res1 = set()       
        res2 = set()
        res3 = set()
        res2.clear()
        res1.clear()
        c1={}
        c2={}
        c1=copy.deepcopy(s1)
        c2=copy.deepcopy(s2)
#        print('c1 length is %d\n' %len(c1))
#        print('c2 length is %d\n' %len(c2))
#        time.sleep(5)
       
        for k2, v2 in c2.items():
            for k1, v1 in c1.items():
                if len(v2)>2 and len(v1)>2: 
                    if Coverage(v1, v2):
                        res1.add(k1)
                    if DeltaCoverage(v1, v2, float(delta), sample):
                        res2.add(k1)
                        res3.add(k2)
#                        c1.pop(k1)
#                        c2.pop(k2)

        delta_precison = len(res3) * 1.0 / sample_block_num
        delta_recall = len(res2) * 1.0 / original_block_num
        print (len(res1), len(res2))
        print ('---------------------------------')
        print ('Delta Precision with delta (' + str(delta) + '): ' + str(delta_precison))
        print ('Delta Recall with delta (' + str(delta) + '): ' + str(delta_recall))    
        result_lines.append([str(delta)+'_precison', 'Entire Graph', delta_precison])
        result_lines.append([str(delta)+'_recall', 'Entire Graph', delta_recall])
    
#    precision = len(res1) * 1.0 / sample_block_num
#    recall = len(res1) * 1.0 / original_block_num
#    print ('Precision: ' + str(precision))
#    print ('Recall: ' + str(recall))
#    result_lines.append(['Precison', 'Entire Graph', precision])
#    result_lines.append(['Recall', 'Entire Graph', recall])

    print ('block_size_median: ' + str(block_size_median))
    print ('sample_block_num: ' + str(sample_block_num))
    print ('original_block_num: ' + str(original_block_num))
    result_lines.append(['Block_size_median', 'Entire Graph', block_size_median])
    result_lines.append(['Sample_block_num', 'Entire Graph', sample_block_num])
    result_lines.append(['Original_block_num', 'Entire Graph', original_block_num])
    return result_lines
'''
def DeltaCoverage_weighted(s1, s2, threshold, sample):
#    delta = len(s1 & s2) * 1.0 / len(s2 | s1)
#    print("intersection of two set is %d" % (len(s1 & s2)))
#    print("union of two set is %d" % (len(s1 | s2)))

    delta = len(s1 & s2) * 1.0 / min(len(s1),len(s2))
    delta2 = len(s1 & s2) * 1.0 / (len(s1) * sample)
    delta3 = len(s1 & s2)
    if delta>= threshold and delta3>= 2:#threshold:       
#    if delta>= threshold and delta2>= 0.5:#threshold:    
#        print("\n\n\nnew original set:delta=%0.2f, delta2=%0.2f\n"%(delta, delta2))
#        print(list(s1))
#        print("sample set \n")
#        print(list(s2))   
#        time.sleep(5)
        return True
        
#    if delta>= threshold or delta2>= 0.5:#threshold:
#        return True 
    else:
        return False

def delta_rank_match_plot(block_match_ratio,clustering_current):
    
    block_match_ratio_sequence=sorted(block_match_ratio,reverse=True) # degree sequence
    #print "Degree sequence", degree_sequence
    dmax=max(block_match_ratio_sequence)
    
    num_of_block=len(block_match_ratio_sequence)
    
#    fig, ax= plt.figure() #size of the whole graph
    fig = plt.figure()
    ax = fig.add_subplot(111)
#    fig, ax= plt.figure(figsize=(16, 20)) #size of the whole graph
#    fig.subplots_adjust(wspace=0.5, hspace=0.5, top=0.90, bottom=0.05)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(18)

#    plt.loglog(degree_sequence,'b-',marker='o')
    x=[(i+1)/float(num_of_block) for i in range(num_of_block)]
    plt.plot(x, block_match_ratio_sequence,'b-',marker='o', markersize=8)
#    plt.bar(block_match_ratio_sequence)
    plt.title("Block match ratio plot", fontsize=18)
    plt.ylabel("value")
    plt.xlabel("rank")
    plt.ylim([0,1.1])
    plt.xlim([1/float(num_of_block) ,1.0])
    
    output_file = clustering_current.split('.')[0]+"_rank_match_plot.pdf"
    #plt.savefig('Imdr_baselines.pdf', format='pdf', dpi=1200)
    plt.savefig(output_file, format='pdf', dpi=1200)
#    plt.show()
    
def match_ratio(c1,c2):
    res_match=[]
    for k1, v1 in c1.items():
        max_match_sample= 0#for each sample block, initiate the max_mathch=0
        match_sample_block=set()
        for k2, v2 in c2.items():
            if len(v2)>2 and len(v1)>2:                  
                intection_set_len = len(v1 & v2)                        
                if intection_set_len > max_match_sample:
                    max_match_sample = intection_set_len
                    match_sample_block=v2
                else:
#                   print('no lager match')
                    continue                                                                              
        res_match.append((len(match_sample_block & v1) * 1.0)/len(v1))#original
    return res_match
    
def Compare3_backup(s1, s2, deltas, sample, clustering_current):
    result_lines = []
    precision = 0.0
    recall = 0.0

    sample_block_num= 0
    original_block_num = 0
    block_size=[]
    block_size_median=0
    
    for kk1, vv1 in s1.items():
        if len(vv1)>2:
            original_block_num=original_block_num+1
            
    
    for kk2, vv2 in s2.items():
        if len(vv2)>2:
            sample_block_num=sample_block_num+1
            block_size.append(len(vv2))
    block_size_median = np.median(block_size)

    for delta in deltas:
        res1 = set()       
        res2 = set()
        res3 = set()
        res_orig=list()
        res_samp=list()
        res2.clear()
        res1.clear()
        c1={}
        c2={}
        c1=copy.deepcopy(s1)
        c2=copy.deepcopy(s2)
#        print('c1 length is %d\n' %len(c1))
#        print('c2 length is %d\n' %len(c2))
#        time.sleep(5)

#        delta-coverage(pi(S)) for precision       
        for k2, v2 in c2.items():
            max_match_orig= 0#for each sample block, initiate the max_mathch=0
            match_orig_block=set()
            for k1, v1 in c1.items():
                if len(v2)>2 and len(v1)>2: 
#                    if Coverage(v1, v2):
#                        res1.add(k1)
#                    if DeltaCoverage(v1, v2, float(delta), sample):
#                        res2.add(k1)
#                        res3.add(k2)
                        #c1.pop(k1)
#                    if DeltaCoverage_sep(v1, v2, float(delta), sample)==1:
#                        res2.add(k1)# original, because of sampling, mapbe it is too low.
#                        res3.add(k2)# sampling 
#                    elif DeltaCoverage_sep(v1, v2, float(delta), sample)==2:
#                        res3.add(k2)
#                    else:
#                        continue                    
                    if DeltaCoverage_weighted(v1, v2, float(delta), sample):
                        intection_set_len = len(v1 & v2)                        
                        if intection_set_len > max_match_orig:
                            max_match_orig = intection_set_len
                            match_orig_block=v1
                        else:
                        
#                            print('no lager match')
#                            print(intection_set_len)
                            continue
                else:
#                    print('less than 3 clusters.')
#                    print("sample set")
#                    print(list(v2))
                    continue                                                                                
            res_samp.append((len(match_orig_block & v2) * 1.0)/len(v2))#sample
            
#        delta-coverage(pi(G)) for recall.  
        for k1, v1 in c1.items():
            max_match_sample= 0#for each sample block, initiate the max_mathch=0
            match_sample_block=set()
            for k2, v2 in c2.items():
                if len(v2)>2 and len(v1)>2:                  
                    if DeltaCoverage_weighted(v1, v2, float(delta), sample):
                        intection_set_len = len(v1 & v2)                        
                        if intection_set_len > max_match_sample:
                            max_match_sample = intection_set_len
                            match_sample_block=v2
                        else:
#                            print('no lager match')
                            continue
                else:
#                    print('less than 3 clusters.')
                    continue                                                                                
            res_orig.append((len(match_sample_block & v1) * 1.0)/len(v1))#original


#        if delta == deltas[-1]:
        match_block_ratio=match_ratio(c1,c2)
        sample_list = [1]*len(match_block_ratio)
        block_match_ave=np.array(match_block_ratio)-np.asarray(sample_list)*sample
            
        block_match_sum=np.square(block_match_ave)
        ase= block_match_sum.sum()/len(match_block_ratio)
        print("ASE:",ase)
            
        block_match_sqrt=np.sqrt(block_match_sum.sum())/sample
        NLS=1-block_match_sqrt/len(match_block_ratio)
        print("NLS:",NLS)
        
        if delta == deltas[0]:
            delta_rank_match_plot(match_block_ratio,clustering_current)
        
        sum_weight_coverage_S = sum(res_samp)
        sum_weight_coverage_G = sum(res_orig)
        delta_precison = sum_weight_coverage_S * 1.0 / sample_block_num
        delta_recall = sum_weight_coverage_G * 1.0 / original_block_num        
#        delta_precison = len(res3) * 1.0 / sample_block_num
#        delta_recall = len(res2) * 1.0 / original_block_num
        
        if delta_precison+delta_recall !=0:
            delta_composite = 2.0*delta_precison*delta_recall* 1.0 / (delta_precison+delta_recall)
        else:
            delta_composite = 0
        
        print (len(res2), len(res3))
        print ('---------------------------------')
  
        print ('Delta Precision with delta (' + str(delta) + '): ' + str(delta_precison))
        print ('Delta Recall with delta (' + str(delta) + '): ' + str(delta_recall)) 
        print ('Delta F-measure with delta (' + str(delta) + '): ' + str(delta_composite)) 
        result_lines.append([str(delta)+'_precison', 'Entire Graph', delta_precison])
        result_lines.append([str(delta)+'_recall', 'Entire Graph', delta_recall])
#        result_lines.append([str(delta)+'_composite', 'Entire Graph', delta_composite])
    
#    precision = len(res1) * 1.0 / sample_block_num
#    recall = len(res1) * 1.0 / original_block_num
#    print ('Precision: ' + str(precision))
#    print ('Recall: ' + str(recall))
#    result_lines.append(['Precison', 'Entire Graph', precision])
#    result_lines.append(['Recall', 'Entire Graph', recall])

    print ('block_size_median: ' + str(block_size_median))
    result_lines.append(['Block_size_median', 'Entire Graph', block_size_median])
    result_lines.append(['Sample_block_num', 'Entire Graph', sample_block_num])
    result_lines.append(['Original_block_num', 'Entire Graph', original_block_num])
    
#    print('difference is:'+ str(abs(original_block_num-sample_block_num)))
#    print('max is:'+ str(max(original_block_num-sample_block_num)))
    print ('sample_block_num: ' + str(sample_block_num))
    print ('original_block_num: ' + str(original_block_num))
    ANC=1.0-(abs(original_block_num-sample_block_num)*1.0/max(original_block_num,sample_block_num))
    print ('ANC: ' + str(ANC))
    result_lines.append(['ANC', 'Entire Graph', ANC])   
    result_lines.append(['NLS', 'Entire Graph', NLS])   
    return result_lines

def Compare3(s1, s2, deltas, sample, clustering_current):
    result_lines = []
    precision = 0.0
    recall = 0.0

    sample_block_num= 0
    original_block_num = 0
    block_size=[]
    block_size_median=0
    
    for kk1, vv1 in s1.items():
        if len(vv1)>=1:
            original_block_num=original_block_num+1
            
    
    for kk2, vv2 in s2.items():
        if len(vv2)>=1:
            sample_block_num=sample_block_num+1
            block_size.append(len(vv2))
    block_size_median = np.median(np.array(block_size))

    for delta in deltas:
        res1 = set()       
        res2 = set()
        res3 = set()
        res_orig=list()
        res_samp=list()
        res2.clear()
        res1.clear()
        c1={}
        c2={}
        c1=copy.deepcopy(s1)
        c2=copy.deepcopy(s2)
#        print('c1 length is %d\n' %len(c1))
#        print('c2 length is %d\n' %len(c2))
#        time.sleep(5)

#        delta-coverage(pi(S)) for precision       
        for k2, v2 in c2.items():
            max_match_orig= 0#for each sample block, initiate the max_mathch=0
            match_orig_block=set()
            for k1, v1 in c1.items():
                if len(v2)>=2 and len(v1)>=2: 
#                    if Coverage(v1, v2):
#                        res1.add(k1)
#                    if DeltaCoverage(v1, v2, float(delta), sample):
#                        res2.add(k1)
#                        res3.add(k2)
                        #c1.pop(k1)
#                    if DeltaCoverage_sep(v1, v2, float(delta), sample)==1:
#                        res2.add(k1)# original, because of sampling, mapbe it is too low.
#                        res3.add(k2)# sampling 
#                    elif DeltaCoverage_sep(v1, v2, float(delta), sample)==2:
#                        res3.add(k2)
#                    else:
#                        continue                    
                    if DeltaCoverage_weighted(v1, v2, float(delta), sample):
                        intection_set_len = len(v1 & v2)                        
                        if intection_set_len > max_match_orig:
                            max_match_orig = intection_set_len
                            match_orig_block=v1
                        else:
                        
#                            print('no lager match')
#                            print(intection_set_len)
                            continue
                else:
#                    print('less than 3 clusters.')
#                    print("sample set")
#                    print(list(v2))
                    continue                                                                                
            res_samp.append((len(match_orig_block & v2) * 1.0)/len(v2))#sample
            
#        delta-coverage(pi(G)) for recall.  
        for k1, v1 in c1.items():
            max_match_sample= 0#for each sample block, initiate the max_mathch=0
            match_sample_block=set()
            for k2, v2 in c2.items():
                if len(v2)>=2 and len(v1)>=2:                  
                    if DeltaCoverage_weighted(v1, v2, float(delta), sample):
                        intection_set_len = len(v1 & v2)                        
                        if intection_set_len > max_match_sample:
                            max_match_sample = intection_set_len
                            match_sample_block=v2
                        else:
#                            print('no lager match')
                            continue
                else:
#                    print('less than 3 clusters.')
                    continue                                                                                
            res_orig.append((len(match_sample_block & v1) * 1.0)/len(v1))#original
        
        sum_weight_coverage_S = sum(res_samp)
        sum_weight_coverage_G = sum(res_orig)
        
        if sample_block_num!=0:
            delta_precison = sum_weight_coverage_S * 1.0 / sample_block_num
        else:
            delta_precison = 0

        if original_block_num !=0:           
            delta_recall = sum_weight_coverage_G * 1.0 / original_block_num 
        else:
            delta_recall=0

#        delta_precison = len(res3) * 1.0 / sample_block_num
#        delta_recall = len(res2) * 1.0 / original_block_num
        
        if delta_precison+delta_recall !=0:
            delta_composite = 2.0*delta_precison*delta_recall* 1.0 / (delta_precison+delta_recall)
        else:
            delta_composite = 0
        
        print (len(res2), len(res3))
        print ('---------------------------------')
  
        print ('Delta Precision with delta (' + str(delta) + '): ' + str(delta_precison))
        print ('Delta Recall with delta (' + str(delta) + '): ' + str(delta_recall)) 
        print ('Delta F-measure with delta (' + str(delta) + '): ' + str(delta_composite)) 
        result_lines.append([str(delta)+'_precison', 'Entire Graph', delta_precison])
        result_lines.append([str(delta)+'_recall', 'Entire Graph', delta_recall])
#        result_lines.append([str(delta)+'_composite', 'Entire Graph', delta_composite])
    
    if 1:
        match_block_ratio=match_ratio(s1,s2)
        sample_list = [1]*len(match_block_ratio)
        block_match_ave=np.array(match_block_ratio)-np.asarray(sample_list)*sample          
        abs_difference = np.array(map(abs, list(block_match_ave)))
        sum_all=0
        for i in range(len(block_match_ave)-1):
            if block_match_ave[i]<0:
                sum_all = sum_all+abs_difference[i]/sample;
            else:
                sum_all = sum_all+abs_difference[i]/match_block_ratio[i];
        NLS=1-sum_all/len(match_block_ratio) 
        print("NLS:",NLS)
        if 0:
            delta_rank_match_plot(match_block_ratio,clustering_current) 

#calculate NLS divied by p                  
#        aver_difference=abs_difference/sample
#        NLS=1-aver_difference.sum()/len(match_block_ratio)
#        print("NLS:",NLS)
            
#calculate ASE            
#        block_match_sum=np.square(block_match_ave)
#        ase= block_match_sum.sum()/len(match_block_ratio)
#        print("ASE:",ase)
            
#calculate real precision recall               
#    precision = len(res1) * 1.0 / sample_block_num
#    recall = len(res1) * 1.0 / original_block_num
#    print ('Precision: ' + str(precision))
#    print ('Recall: ' + str(recall))
#    result_lines.append(['Precison', 'Entire Graph', precision])
#    result_lines.append(['Recall', 'Entire Graph', recall])          

    print ('block_size_median: ' + str(block_size_median))
    result_lines.append(['Block_size_median', 'Entire Graph', block_size_median])
    result_lines.append(['Sample_block_num', 'Entire Graph', sample_block_num])
    result_lines.append(['Original_block_num', 'Entire Graph', original_block_num])
    
#    print('difference is:'+ str(abs(original_block_num-sample_block_num)))
#    print('max is:'+ str(max(original_block_num-sample_block_num)))
    print ('sample_block_num: ' + str(sample_block_num))
    print ('original_block_num: ' + str(original_block_num))
    ANC=1.0-(abs(original_block_num-sample_block_num)*1.0/max(original_block_num,sample_block_num))
    print ('ANC: ' + str(ANC))
    result_lines.append(['ANC', 'Entire Graph', ANC])  
    result_lines.append(['NLS', 'Entire Graph', NLS])  
 
    return result_lines
    
def sort(a):
    for k in range(len(a)):
        (a[k][0],a[k][1]) = (a[k][1],a[k][0])
    a.sort()
    for k in range(len(a)):
        (a[k][0],a[k][1]) = (a[k][1],a[k][0]) 
        
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

def Get_x_value_from_filename(filename):
    if 'v1_' in filename:
        v_pos = filename.find('v1_');
        v_len = 0;
        x_value = string.atof(filename[( v_pos+3) : (len(filename) - len('.coms'))])  
        print x_value
        return x_value

def ReadInData_coms(coms_filename):
    coms_file = open(coms_filename, 'r');
    coms_dir = {};
    i = 1;

    for line in coms_file:

        cluster = [];
        nodes = line.split(' ')
       
        for node in nodes:
            cluster.append(string.atoi(node))
            
        set_cluster = set(cluster)
        coms_dir[i] = set_cluster
        i = i + 1;

    return coms_dir
    
def collect_result_data (handle_collect_file, data): #to collect data from a single file and output into a summarized file


    for i in range(len(data)):
        for j in range(len(data[i])):
#            print(data[i])
            handle_collect_file.write(str(data[i][j])+ ', ')
        handle_collect_file.write('\n')


if __name__ == '__main__':

    #ground_first_level = ReadInData('binary_networks/community_first_level_sample_p100_v1.dat')
    #ground_second_level= ReadInData('binary_networks/community_second_level_sample_p100_v1.dat')

    handleArgs()
    deltas=[1.0, 0.9, 0.8, 0.7, 0.6, 0.5]
    sample = args.samplerate
   
    x = [];
    y_precision_first_level = [];
    y_precision_second_level = [];
    y_NMI_first_level = [];
    y_NMI_second_level = [];
    y_recall_first_level = [];
    y_recall_second_level = [];
    y_cluster = [];
    
#    p = 50    
#    rootdir = "generated_benches_u1_10_u2_20_p_50_condition_metropolis_subgraph/n_2000";

    
#    algorithm_list=['blondel', 'infomap','label_propagation','oslom','mod_opt']
    algorithm_list = args.clustering_file_names
#    metric_list=['NMI', 'ARS','1.0-precision','1.0-recall','0.9-precision','0.9-recall','0.8-precision','0.8-recall','0.7-precision','0.7-recall','0.6-precision','0.6-recall','0.5-precision','0.5-recall']
   

    rootdir = args.bench_directory_stem
    for ground_truth_file_name in args.ground_truth_file:

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
        
        average_NMI = []
        average_ARS = []
        average_precision= []
        average_recall= []
        average_9_precision= []
        average_9_recall= []
        average_8_precision= []
        average_8_recall= []
        average_7_precision= []
        average_7_recall= []
        average_6_precision= []
        average_6_recall= []
        average_5_precision= []
        average_5_recall= []
        #stastics information          
        average_sample_cluster_num=[]
        average_original_cluster_num=[]
        average_mediate_size=[]
        average_ANC=[] 
        average_NLS=[]         


        
#        ground_truth_file_name = args.ground_truth_file
        ground_first_level = ReadInData(rootdir+ground_truth_file_name)  
        gold_standard_first_level = partitionFromFile(rootdir+ground_truth_file_name, '\t')
#    gold_standard_second_level = partitionFromFile(rootdir+'/community_second_level_sample_p' + str(p) +'_v1.dat', '\t')
        for i in range(len(algorithm_list)):
            for parent, dirnames, filenames in os.walk(rootdir):
                for filename in filenames:
                    
                    split_filename = filename.split('.');
#                    parent_path = os.path.dirname(split_filename)
#                    current_file= os.path.join(parent_path,filename)
                    
                    
                    if split_filename[0].find(algorithm_list[i]) !=-1 and split_filename[0].find('clustering_v') !=-1 and split_filename[len(split_filename) - 1] == 'dat':
                        print("match %s algorithm is the file %s" %(algorithm_list[i],filename))    
                        full_name = os.path.join(parent, filename);
                        partition = partitionFromFile(full_name, '\t')
                        result_lines_reference = runComparisonGraphMetrics(partition, gold_standard_first_level, False) 
                        y_value_NMI_micro.append(string.atof(result_lines_reference[0][2])) 
                        y_value_ARS_micro.append(string.atof(result_lines_reference[1][2]))                    
                                          
                        partition2 = ReadInData(full_name)
                        first_result_lines = Compare3(copy.deepcopy(ground_first_level), partition2, deltas, sample, full_name)
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

                      
                        
#                    print("add one more file")
            print("processing %s:"%algorithm_list[i])             
                      
#            average_NMI.append( str(np.array(y_value_NMI_micro).mean())+'+'+ str(np.array(y_value_NMI_micro).std()))          
#            average_ARS.append(str(np.array(y_value_ARS_micro).mean())+'+'+str(np.array(y_value_ARS_micro).std()))           
#            average_precision.append(str(np.array(y_value_precision_micro).mean())+'+'+str(np.array(y_value_precision_micro).std()))           
#            average_recall.append(str(np.array(y_value_recall_micro).mean())+'+'+ str(np.array(y_value_recall_micro).std()))
#            average_9_precision.append(str(np.array(y_value_9_precision_micro).mean())+'+'+str(np.array(y_value_precision_micro).std()))
#            average_9_recall.append(str(np.array(y_value_9_recall_micro).mean())+'+'+ str(np.array(y_value_recall_micro).std()))
#            average_8_precision.append(str(np.array(y_value_8_precision_micro).mean())+'+'+str(np.array(y_value_precision_micro).std()))
#            average_8_recall.append(str(np.array(y_value_8_recall_micro).mean())+'+'+ str(np.array(y_value_recall_micro).std()))
#            average_7_precision.append(str(np.array(y_value_7_precision_micro).mean())+'+'+str(np.array(y_value_precision_micro).std()))
#            average_7_recall.append(str(np.array(y_value_7_recall_micro).mean())+'+'+ str(np.array(y_value_recall_micro).std()))
#            average_6_precision.append(str(np.array(y_value_6_precision_micro).mean())+'+'+str(np.array(y_value_precision_micro).std()))
#            average_6_recall.append(str(np.array(y_value_6_recall_micro).mean())+'+'+ str(np.array(y_value_recall_micro).std()))
#            average_5_precision.append(str(np.array(y_value_5_precision_micro).mean())+'+'+str(np.array(y_value_precision_micro).std()))
#            average_5_recall.append(str(np.array(y_value_5_recall_micro).mean())+'+'+ str(np.array(y_value_recall_micro).std()))
#            average_sample_cluster_num.append(str(np.array(y_value_sample_cluster_num_micro).mean())+'+'+ str(np.array(y_value_sample_cluster_num_micro).std()))
#            average_original_cluster_num.append(str(np.array(y_value_original_cluster_num_micro).mean())+'+'+ str(np.array(y_value_original_cluster_num_micro).std()))
#            average_mediate_size.append(str(np.array(y_value_mediate_size_micro).mean())+'+'+ str(np.array(y_value_mediate_size_micro).std()))
#            average_ANC.append(str(np.array(y_value_ANC).mean())+'+'+ str(np.array(y_value_ANC).std()))
#            average_NLS.append(str(np.array(y_value_ANC).mean())+'+'+ str(np.array(y_value_ANC).std()))

            
            average_NMI.append( str(np.array(y_value_NMI_micro).mean()))        
            average_ARS.append(str(np.array(y_value_ARS_micro).mean()))         
            average_precision.append(str(np.array(y_value_precision_micro).mean()))         
            average_recall.append(str(np.array(y_value_recall_micro).mean()))
            average_9_precision.append(str(np.array(y_value_9_precision_micro).mean()))
            average_9_recall.append(str(np.array(y_value_9_recall_micro).mean()))
            average_8_precision.append(str(np.array(y_value_8_precision_micro).mean()))
            average_8_recall.append(str(np.array(y_value_8_recall_micro).mean()))
            average_7_precision.append(str(np.array(y_value_7_precision_micro).mean()))
            average_7_recall.append(str(np.array(y_value_7_recall_micro).mean()))
            average_6_precision.append(str(np.array(y_value_6_precision_micro).mean()))
            average_6_recall.append(str(np.array(y_value_6_recall_micro).mean()))
            average_5_precision.append(str(np.array(y_value_5_precision_micro).mean()))
            average_5_recall.append(str(np.array(y_value_5_recall_micro).mean()))
            average_sample_cluster_num.append(str(np.array(y_value_sample_cluster_num_micro).mean()))
            average_original_cluster_num.append(str(np.array(y_value_original_cluster_num_micro).mean()))
            average_mediate_size.append(str(np.array(y_value_mediate_size_micro).mean()))
            average_ANC.append(str(np.array(y_value_ANC).mean()))
            average_NLS.append(str(np.array(y_value_NLS).mean()))
            
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
            y_value_sample_cluster_num_micro=[]
            y_value_original_cluster_num_micro=[]
            y_value_mediate_size_micro=[]
            y_value_ANC=[]
            y_value_NLS=[]
            
#        total_result=[average_NMI,average_ARS,average_precision,average_recall,average_9_precision,average_9_recall,average_8_precision, average_8_recall,average_7_precision,average_7_recall,average_6_precision,average_6_recall,average_5_precision,average_5_recall]
        
        total_result=[algorithm_list, average_original_cluster_num,average_sample_cluster_num,average_mediate_size, average_NMI,average_ARS,average_precision,average_recall,average_9_precision,average_9_recall,average_8_precision, average_8_recall,average_7_precision,average_7_recall,average_6_precision,average_6_recall,average_5_precision,average_5_recall,average_ANC, average_NLS]
    #    output_file='zjpcommunity.txt'
        output_file = ground_truth_file_name.split('.')[0]+'_results.csv'
        
        handle_collect_file = open(rootdir+output_file, 'w');
        collect_result_data(handle_collect_file, total_result)
               
                  
'''                    
            if split_filename[len(split_filename) - 1] == 'coms':
                full_name = os.path.join(parent, filename);

                x_value = Get_x_value_from_filename(filename);              
                pos = 0;
                for pos in range(0, len(x)):
		    if x[pos] > x_value:                
                        break;
                    pos = pos + 1;
                x.insert(pos, x_value);
                
                partition = partitionFromFile(full_name + '.dat', ' ')
                result_lines2 = runComparisonGraphMetrics(partition, gold_standard_first_level, False)
                 
         
#                ground_first_level = ReadInData(rootdir+'/community_first_level_sample_p' + str(p) +'_v1.dat')
#                ground_second_level= ReadInData(rootdir+'/community_second_level_sample_p' + str(p) +'_v1.dat')

                parseLancichinettiResults(full_name, full_name + '.dat')
#                partition = ReadInData(full_name + '.dat')

                partition = ReadInData_coms(full_name)
                y_cluster_num = len(partition)
                y_cluster.insert(pos, y_cluster_num)
                #print('The number of clustres in first level is %s and in second level is %s' %(len(ground_first_level), len(ground_second_level)))
                
                
                #the micro results. 
                #NMI 
                NMI_first_result_lines = NMI_score(rootdir+'/community_first_level_sample_p' + str(p) +'_v1.dat', full_name + '.dat',)                                          
                first_result_lines = Compare3(copy.deepcopy(ground_first_level), partition, deltas, sample)
                
                y_value_NMI_micro = string.atof(NMI_first_result_lines[0][2])  
                y_NMI_first_level.insert(pos, y_value_NMI_micro)
                y_value_precision_micro = string.atof(first_result_lines[4][2])  
                #y_precision_first_level.append(y_value)
                y_precision_first_level.insert(pos, y_value_precision_micro)
                y_value_recall_micro = string.atof(first_result_lines[5][2])  
                #y_precision_first_level.append(y_value)
                y_recall_first_level.insert(pos, y_value_recall_micro)
                
                #the macro results.
                #NMI
                NMI_second_result_lines = NMI_score(rootdir+'/community_second_level_sample_p' + str(p) +'_v1.dat', full_name + '.dat',)                
                second_result_lines = Compare3(copy.deepcopy(ground_second_level), partition, deltas, sample)
                
                
                y_value_NMI_macro = string.atof(NMI_second_result_lines[0][2])  
                y_NMI_second_level.insert(pos, y_value_NMI_macro)
                y_value_precision_macro = string.atof(second_result_lines[4][2])  
                y_precision_second_level.insert(pos, y_value_precision_macro)
                y_value_recall_macro = string.atof(second_result_lines[5][2])  
                y_recall_second_level.insert(pos, y_value_recall_macro)

                #print '%s is completed' %filename
                #print 'the full name of the file is' + (os.path.join(parent, filename))

    #print "Finish"
    #print "x:[]"
    #print x
    #print "y_first_level:[]"
    #print y_precision_first_level
    #print "y_second_level:[]"
    #print y_precision_second_level

    #there are two ground-truth clustering in the graph.
    #ground_first_level = ReadInData('binary_networks/community_first_level_sample_p100_v1.dat')
    #ground_second_level= ReadInData('binary_networks/community_second_level_sample_p100_v1.dat')

    #we have several clustering results, and we want to matcha all the results with regard to the two ground truth
    #parseLancichinettiResults('binary_networks/network_v1_0.coms','binary_networks/network_v1_0.dat')
    #parseLancichinettiResults('binary_networks/network_v1_0.9.coms','binary_networks/network_v1_0.9.dat')
    #partition1 = ReadInData('binary_networks/network_v1_0.dat')
    #partition2 = ReadInData('binary_networks/network_v1_0.9.dat')
    #partition = ReadInData('network500/network_v1_0.coms')
    #print('The number of clustres in first level is %s and in second level is %s' %(len(ground_first_level), len(ground_second_level)))
    
    #deltas=[0.8,0.9, 1.0]
    #sample =1.0


    #Compare3(ground_first_level, partition1, deltas, sample)   
    #Compare3(ground_first_level, partition2, deltas, sample)
    #Compare3(ground_second_level, partition1, deltas, sample)
    #Compare3(ground_second_level, partition2, deltas, sample)
  
    #x = [4,1,3,2 ];
    #y_precision_first_level = [1,2,3,4];

#    pl.title('Compare precision between ground first and second level')# give plot a title
#    pl.xlabel('Degree of Precision')# make axis labels
#    pl.ylabel('Precision')
    #pl.xlim(0.8, 1.0)# set axis limits
#    pl.ylim(0.0, 1.1)
#    precision_first_level_line = pl.plot(x, y_precision_first_level, label = 'ground_first_level', color = 'red')# use pylab to plot x and y
#    precision_second_level_line = pl.plot(x, y_precision_second_level, label = 'ground_second_level', color = 'blue')
#    pl.legend(loc = 'best')
    #pl.legend([precision_first_level_line, precision_second_level_line], ('red line', 'green circles'), 'best', numpoints=1)
#    pl.show()# show the plot on the screen
  
#    fig.suptitle('test title', fontsize=20)
#    plt.xlabel('xlabel', fontsize=18)
#    plt.ylabel('ylabel', fontsize=16)

    fig1=plt.figure(1)
    ax =fig1.add_subplot(111, xlabel='x', ylabel='y', title='title')
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(20)
    plt.title('1.0-precision comparison w.r.t. mirco and macro ground-truths', fontsize=15)# give plot a title
    plt.xlabel('Scale of parameters', fontsize=18)# make axis labels
    plt.ylabel('1.0-precision', fontsize=18)
    #pl.xlim(0.8, 1.0)# set axis limits
    plt.ylim(0.0, 1.1)
    precision_first_level_line = plt.plot(x, y_precision_first_level, label = '1.0-precision w.r.t. micro ground-truth', color = 'red')# use pylab to plot x and y
    precision_second_level_line = plt.plot(x, y_precision_second_level, label = '1.0-precision w.r.t. macro ground-truth', color = 'blue')
#    plt.legend(loc = 'upper right')
    plt.legend(loc = 'best', prop={'size':15})
    #pl.legend([precision_first_level_line, precision_second_level_line], ('red line', 'green circles'), 'best', numpoints=1)
    plt.savefig(rootdir+'_precision_result.pdf', format='pdf', dpi=1200)
    plt.show()# show the plot on the screen

    fig2=plt.figure(2)
    ax =fig2.add_subplot(111, xlabel='x', ylabel='y', title='title')
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(20)
    plt.title('1.0-recall comparison w.r.t. mirco and macro ground-truths', fontsize=15)# give plot a title
    plt.xlabel('Scale of parameters', fontsize=18)# make axis labels
    plt.ylabel('1.0-recall', fontsize=18)
    plt.ylim(0.0, 1.1)
    recall_first_level_line = plt.plot(x, y_recall_first_level, label = '1.0-recall w.r.t. micro ground-truth', color = 'red')# use pylab to plot x and y
    recall_second_level_line = plt.plot(x, y_recall_second_level, label = '1.0-recall w.r.t. macro ground-truth', color = 'blue')  
#    plt.legend(loc = 'upper right')
    plt.legend(loc = 'best', prop={'size':15})
    #pl.legend([recall_first_level_line, recall_second_level_line], ('red line', 'green circles'), 'best', numpoints=1)
    plt.savefig(rootdir+'_recall_result.pdf', format='pdf', dpi=1200)
    plt.show()# show the plot on the screen

    fig3=plt.figure(3)
    ax =fig3.add_subplot(111, xlabel='x', ylabel='y', title='title')
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(20)
    plt.title('NMI comparison w.r.t. mirco and macro ground-truths', fontsize=15)# give plot a title
    plt.xlabel('Scale of parameters', fontsize=18)# make axis labels
    plt.ylabel('NMI', fontsize=18)
    plt.ylim(0.0, 1.1)
    NMI_first_level_line = plt.plot(x, y_NMI_first_level, label = 'NMI w.r.t. micro ground-truth', color = 'red')# use pylab to plot x and y
    NMI_second_level_line = plt.plot(x, y_NMI_second_level, label = 'NMI w.r.t. macro ground-truth', color = 'blue')  
    plt.legend(loc = 'best', prop={'size':15})
    #pl.legend([NMI_first_level_line, NMI_second_level_line], ('red line', 'green circles'), 'best', numpoints=1)
    plt.savefig(rootdir+'_NMI_result.pdf', format='pdf', dpi=1200)
    plt.show()# show the plot on the screen

    fig4=plt.figure(4)
    ax =fig4.add_subplot(111, xlabel='x', ylabel='y', title='title')
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(20)
    plt.title('The number of clusters w.r.t. mirco and macro ground-truths', fontsize=15)# give plot a title
    plt.xlabel('Scale of parameters', fontsize=18)# make axis labels
    plt.ylabel('Number of clusters', fontsize=18)
    plt.ylim(0.0, 600)
    cluster_first_level = len(ground_first_level)*np.ones((len(x),1))
    cluster_second_level = len(ground_second_level)*np.ones((len(x),1))
    cluster_num_line = plt.plot(x, y_cluster, label = 'cluster numbers of clustering result', color = 'green') 
    cluster_second_level_line = plt.plot(x, cluster_second_level, label = 'cluster numbers of macro ground-truth', color = 'red')
    cluster_first_level_line = plt.plot(x, cluster_first_level, label = 'cluster numbers of micro ground-truth', color = 'blue')    
    plt.legend(loc = 'best', prop={'size':15})
    #pl.legend([recall_first_level_line, recall_second_level_line], ('red line', 'green circles'), 'best', numpoints=1)
    plt.savefig(rootdir+'_clustering_number_result.pdf', format='pdf', dpi=1200)
    plt.show()# show the plot on the screen

    fig5=plt.figure(5)
    ax =fig5.add_subplot(111, xlabel='x', ylabel='y', title='title')
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(20)

    plt.title('1.0-precision comparison w.r.t. mirco and macro ground-truths', fontsize=15)# give plot a title
    plt.xlabel('Scale of parameters', fontsize=18)# make axis labels
    plt.ylabel('Precision', fontsize=18)
    #pl.xlim(0.8, 1.0)# set axis limits
    plt.ylim(0.0, 1.1)
    precision_first_level_line = plt.plot(x, y_precision_first_level, label = '1.0-precision w.r.t. micro ground-truth', color = 'red')# use pylab to plot x and y (color="red",linewidth=2)
    precision_second_level_line = plt.plot(x, y_precision_second_level, label = '1.0-precision w.r.t. macro ground-truth', color = 'blue')
#    plt.legend(loc = 'upper right')
    plt.legend(loc = 'best', prop={'size':15})
    #pl.legend([precision_first_level_line, precision_second_level_line], ('red line', 'green circles'), 'best', numpoints=1)
    plt.savefig(rootdir+'_precision_result2.pdf', format='pdf', dpi=1200)
    plt.show()# show the plot on the screen
'''



    


    



"""    
#     a = [[1,2,4],[6,5,6],[2,5,9]]
#     sort(a)
#     print(a)
    
#     x = np.array([[1,2,4],[6,5,6],[2,5,9]])
#    x[np.argsort(x)]
#    print x
    comm = np.loadtxt('community.dat',dtype='int32') 
    myresult = np.loadtxt('network_metisResult',dtype='int32')
    myresult_sorted= np.array( sorted(myresult, key=lambda myresult : myresult[0]))
    myresult_sampled= np.array(myresult_sorted[:,0])
    myresult_label= np.array(myresult_sorted[:,1])
    print (len(myresult))
#comm_sorted =np.array( sorted(comm, key=lambda comm : comm[0]))
    
    comm_sorted =np.array( sorted(comm, key=lambda myresult_sorted : myresult_sorted[0]))
    a= myresult[:,0]
    b=comm_sorted[myresult[:,0],:]
    
    comm_selected= np.array(comm_sorted[myresult[:,0],:])
    
    truelabel_whole= np.array(comm_sorted[:,0])

    truelabel= np.array(comm_selected[:,1])
    print (len(truelabel))
    print (purity_score(myresult_label, truelabel))
    
    print (normalized_mutual_info_score([0, 0, 0, 0], [0, 1, 2, 3]))
    

    print ("resut is ",normalized_mutual_info_score(truelabel.tolist(), myresult_label.tolist()))
# 
#     x = [1, 2, 3, 4, 5]# Make an array of x values

y = [1, 4, 9, 16, 25]# Make an array of y values for each x value

pl.plot(x, y)# use pylab to plot x and y

pl.show()# show the plot on the screen
#     clus = np.array([1, 4, 4, 4, 4, 4, 3, 3, 2, 2, 3, 1, 1])
#     clas = np.array([5, 1, 2, 2, 2, 3, 3, 3, 1, 1, 1, 5, 2])
#     print purity_score(clus, clas)
"""

    
    
    
