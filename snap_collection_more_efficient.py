import os
import time
import copy

def update_dynamic_edge_stream(dynamic_edge_stream, last_snapshot, current_snapshot):
    #firstly, we collect all the edges from last_snapshot which have already been deleted in current_snapshot
    #these edges should been added into dynamic_edge_stream as the edge-deleting request in the time window of current_snapshot
    
    j = 0; ##
    
    lower_bound = 0
    
   
    time.sleep(1)
    for edge in last_snapshot:
        while (edge[0] > current_snapshot[lower_bound][0]): #aimed range of edges has not been reached, we should increase lower bound
            lower_bound = lower_bound + 1

        if (edge[0] == current_snapshot[lower_bound][0]): #aimed range has been reached , we should check whether the edge in last snapshot is still in current snapshot
            i = lower_bound
            while (edge[0] == current_snapshot[i][0]): 
                if (edge[1] == current_snapshot[i][1]): #edge in last snapshot is still in current snapshot
                    break; 
                i = i + 1           
            if (edge[0] < current_snapshot[i][0]): #edge in last snapshot is not in current snapshot
                del_edge_stream = "! " + str(edge[0]) + ' ' +  str(edge[1]) + '\n' 
                dynamic_edge_stream.append(del_edge_stream)
            
        elif (edge[0] < current_snapshot[lower_bound][0]):   #aimed range has been passed, edge in last snapshot is not in current snapshot
            del_edge_stream = "! " + str(edge[0]) + ' ' +  str(edge[1]) + '\n'
            dynamic_edge_stream.append(del_edge_stream)

        if (j % 10000 == 0): ##
            print "j = %d" %j
        j = j + 1

    #then, we collect all the edges from current_snapshot which does not exist in last_snapshot
    #these edges should been added into dynamic_edge_stream as the edge-adding request in the time window of current_snapshot
    
    k = 0; ##

    lower_bound = 0
        
    for edge in current_snapshot:

        while (edge[0] > last_snapshot[lower_bound][0]):
            lower_bound = lower_bound + 1

        if (edge[0] == last_snapshot[lower_bound][0]): #aimed range has been reached , we should check whether the edge in current snapshot is still in last snapshot
            i = lower_bound
            while (edge[0] == last_snapshot[i][0]): 
                if (edge[1] == last_snapshot[i][1]): #edge in last snapshot is still in last snapshot
                    break; 
                i = i + 1           
            if (edge[0] < last_snapshot[i][0]): #edge in last snapshot is not in last snapshot, so it a new edge
                edge_stream = str(edge[0]) + ' ' +  str(edge[1]) + '\n'
                dynamic_edge_stream.append(edge_stream)
            
        elif (edge[0] < last_snapshot[lower_bound][0]):   #aimed range has been passed, edge in last snapshot is not in last snapshot, so it is a new edge
            edge_stream = str(edge[0]) + ' ' +  str(edge[1]) + '\n'
            dynamic_edge_stream.append(edge_stream)
	   
        if (k % 10000 == 0):
            print "k = %d" %k ##   
        k = k + 1 ##


    return dynamic_edge_stream

#######################################################
#function: build a edge snapshot based on snapshot file
#######################################################
def build_snapshot(fp_snapshot):
    
    i = 0 ##
	
    snapshot = [];
    for line in fp_snapshot:
        line = line.strip('\n')
        from_id, to_id = line.split(' ')

        edge = []
        edge.append(int(from_id))
        edge.append(int(to_id))
        snapshot.append(tuple(edge))

        if (i % 10000 == 0): ##
            print "i = %d" %i#
        i = i + 1 ##        
        
    return snapshot

def generate_zero_series(index):
    num_of_zero = 4 - ( index / 10 )
    i = 0;
    zero_series = ''
    while ( i < num_of_zero ):
        zero_series = zero_series + '0'
        i = i + 1
    return zero_series


if __name__ == '__main__':

    dynamic_edge_stream = []

    last_snapshot = []
    current_snapshot = []
#    rootdir = "./original_snapshot_20w/"
    
#    rootdir = "./original_snapshot_merge/"
#    rootdir = "./original_snapshot_grow/"
    rootdir = "./original_snapshot_mixed/"
    #rootdir = "../test/"
    

    prefix = "output-prefix.t"
    suffix = ".graph"
     
    i = 0;
 
    print "start collection process\n"
 
    while ( i < 10):
        
        zero_series = generate_zero_series(i)
        filename = prefix + zero_series + str(i) + suffix
        fp_snapshot = open(rootdir + filename, 'r')
 
        print "start collectiong file " + rootdir + filename


        if (i == 0):
            current_snapshot = build_snapshot(fp_snapshot)
            for edge in current_snapshot:
                dynamic_edge_stream.append(str(edge[0]) + ' ' + str(edge[1]) + '\n')
        else:
            last_snapshot = copy.deepcopy(current_snapshot)
            current_snapshot = build_snapshot(fp_snapshot)
            dynamic_edge_stream = update_dynamic_edge_stream(dynamic_edge_stream, last_snapshot, current_snapshot)
            
        fp_snapshot.close()
        dynamic_edge_stream.append(filename + '\n')

        i = i + 1
        
        print filename + " collection completed"

    print "star writing process"
    collection_file = rootdir + "collection.graph"
    fp_collection = open(collection_file, 'w')
    for line in dynamic_edge_stream:
	      fp_collection.write(line)
    fp_collection.close()
    print "writing completed"