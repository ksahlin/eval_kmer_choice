
DATASETS='reads1 reads2'.split()
TOOLS='optimal_k kmergenie'.split()
OUTBASE='/Users/ksahlin/_tmp/Optimal_k/OUT/'
#####################################
# standard python functions

import re
def get_kmer_genie_params():
    pass
    return 4,2

def get_optimal_k_params():
    pass
    return 22,3

def get_memory_and_runtime():
    pass


#####################################


rule all:
    input:
        # 'out/memory_table.tex',
        # 'out/k_choice.tex',
        # 'out/eval_table.tex',    
        expand("/Users/ksahlin/_tmp/Optimal_k/OUT/{tool}_{dataset}.dat", tool=TOOLS, dataset=DATASETS)    


rule optimal_k:
    input: "/Users/ksahlin/_tmp/Optimal_k/testdata_optimal_k/{dataset}.fa"
    output: csv=OUTBASE+"optimal_k_{dataset}.dat", stderr=OUTBASE+"optimal_k_{dataset}.stderr"
    run:
        for i in range(10):
            pass
        shell("/Users/ksahlin/_tmp/Optimal_k/./test_prgrm1.sh 1> {output.csv} 2> {output.stderr}")

rule kmergenie:
    input: "/Users/ksahlin/_tmp/Optimal_k/testdata_optimal_k/{dataset}.fa"
    output: csv=OUTBASE+"kmergenie_{dataset}.dat", stderr=OUTBASE+"kmergenie_{dataset}.stderr"
    run:
        for i in range(10):
            pass
        shell("/Users/ksahlin/_tmp/Optimal_k/./test_prgrm2.sh 1> {output.csv} 2> {output.stderr}")

rule unitiger:
    input: csv=OUTBASE+"{tool}_{dataset}.dat", reads="/Users/ksahlin/_tmp/Optimal_k/testdata_optimal_k/{dataset}.fa"
    output: stats=OUTBASE+"{tool}_{dataset}.stdout", contigs=OUTBASE+"{tool}_{dataset}.contigs.fa"
    run:
        if {tool} == "optimal_k":
            k,a = get_optimal_k_params()
        elif {tool} == "kmergenie":
            k,a = get_kmer_genie_params()

        shell('unitiger {input.reads} {0} {1}'.format(k,a))    