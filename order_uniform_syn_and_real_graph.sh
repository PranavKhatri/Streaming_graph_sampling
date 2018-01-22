###real_network u is fixed to 0.45, enron, facebook, reality, slashdot
#python uniform_graph_realgraph.py -n 128 --mu 0.45 -ds enron -sp 9 -b datasets/Dataset_snopshots/ -o powlaw_degree_small_snapshot_graph_for_streaming_sampling/
#python uniform_graph_realgraph.py -n 128 --mu 0.45 -ds facebook -sp 12 -b datasets/Dataset_snopshots/ -o powlaw_degree_small_snapshot_graph_for_streaming_sampling/
#python uniform_graph_realgraph.py -n 128 --mu 0.45 -ds reality -sp 10 -b datasets/Dataset_snopshots/ -o powlaw_degree_small_snapshot_graph_for_streaming_sampling/
#python uniform_graph_realgraph.py -n 128 --mu 0.45 -ds slashdot -sp 12 -b datasets/Dataset_snopshots/ -o powlaw_degree_small_snapshot_graph_for_streaming_sampling/
#
###synthetic graph (n=128)（注意check gen_dynamic的数据使用后面多少数目的顶点指令生成的）
python uniform_graph_syn.py -n 128 --mu 0.45 -ds mergesplit -sp 5 -b gen_dynamic_20151102/ -o powlaw_degree_small_snapshot_graph_for_streaming_sampling/
python uniform_graph_syn.py -n 128 --mu 0.45 -ds expand -sp 5 -b gen_dynamic_20151102/ -o powlaw_degree_small_snapshot_graph_for_streaming_sampling/
#python uniform_graph_syn.py -n 128 --mu 0.45 -ds birthdeath -sp 5 -b gen_dynamic_20151102/ -o powlaw_degree_small_snapshot_graph_for_streaming_sampling/
#python uniform_graph_syn.py -n 128 --mu 0.45 -ds switch -sp 5 -b gen_dynamic_20151102/ -o powlaw_degree_small_snapshot_graph_for_streaming_sampling/
#python uniform_graph_syn.py -n 128 --mu 0.45 -ds hide -sp 5 -b gen_dynamic_20151102/ -o powlaw_degree_small_snapshot_graph_for_streaming_sampling/
#
###synthetic graph (n=1000)（注意check gen_dynamic的数据使用后面多少数目的顶点指令生成的）
python uniform_graph_syn.py -n 1000 --mu 0.20 -ds mergesplit -sp 5 -b gen_dynamic_20151102/ -o powlaw_degree_small_snapshot_graph_for_streaming_sampling/
python uniform_graph_syn.py -n 1000 --mu 0.20 -ds expand -sp 5 -b gen_dynamic_20151102/ -o powlaw_degree_small_snapshot_graph_for_streaming_sampling/
#
################################################################################
###############################################################################
########our benchmark graph of my streaming clustering paper.##################
#./bench_expand -N 1000 -muw 0.2 -k 9 -maxk 15 -s 5 -minc 40 -maxc 60 -expand 2 -contract 2 -r 0.25
#./bench_mergesplit -N 1000 -muw 0.2 -k 9 -maxk 15 -s 5 -minc 40 -maxc 60 -merge 2 -split 2
#./bench_birthdeath -N 1000 -muw 0.2 -k 9 -maxk 15 -s 5 -minc 40 -maxc 60 -birth 2 -death 2
#./bench_switch -N 1000 -muw 0.2 -k 9 -maxk 15 -s 5 -minc 40 -maxc 60 -p 0.1
#./bench_hide -N 1000 -muw 0.2 -k 9 -maxk 15 -s 5 -minc 40 -maxc 60 -hide 2
########new benchmark graph.###################
#./bench_expand -N 128 -muw 0.45 -k 16 -maxk 31 -s 5 -minc 32 -maxc 32 -expand 2 -contract 2 -r 0.25
#./bench_mergesplit -N 128 -muw 0.45 -k 16 -maxk 31 -s 5 -minc 32 -maxc 32 -merge 2 -split 2
#./bench_birthdeath -N 128 -muw 0.45 -k 16 -maxk 31 -s 5 -minc 32 -maxc 32 -birth 2 -death 2
#./bench_switch -N 128 -muw 0.45 -k 16 -maxk 31 -s 5 -minc 32 -maxc 32 -p 0.1
#./bench_hide -N 128 -muw 0.45 -k 16 -maxk 31 -s 5 -minc 32 -maxc 32 -hide 0.3
