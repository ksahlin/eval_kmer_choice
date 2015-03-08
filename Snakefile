
DATASETS='reads1 reads2'.split()
TOOLS='optimal_k kmergenie'.split()
OUTBASE='/Users/ksahlin/_tmp/Optimal_k/OUT/'

PYTHON2="/usr/bin/python2.7"
GNUTIME="/usr/bin/time -lp" # "/usr/bin/time -v" if Linux
#####################################
# standard python functions

import re
import os
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
        #expand("/Users/ksahlin/_tmp/Optimal_k/OUT/{tool}_{dataset}.dat", tool=TOOLS, dataset=DATASETS)    
        expand(OUTBASE+"{tool}_{dataset}.stdout", tool=TOOLS, dataset=DATASETS),
        expand(OUTBASE+"{tool}_{dataset}.contigs.fa", tool=TOOLS, dataset=DATASETS)

# rule dataset:
#     input: "/Users/ksahlin/_tmp/Optimal_k/testdata_optimal_k/{dataset}.fa"
#     output: "/Users/ksahlin/_tmp/Optimal_k/testdata_optimal_k/{dataset}.fa", temp("/Users/ksahlin/_tmp/Optimal_k/testdata_optimal_k/{dataset}.check")
#     run:
#         if os.path.isfile(input) :
#             print(input,'found')
#         else:
#             print(input,'not found')


rule optimal_k:
    input: "/Users/ksahlin/_tmp/Optimal_k/testdata_optimal_k/{dataset}.fa"
    output: csv=OUTBASE+"optimal_k_{dataset}.dat", stderr=OUTBASE+"optimal_k_{dataset}.stderr"
    run:
        for i in range(10):
            pass
        shell(" {GNUTIME} /Users/ksahlin/_tmp/Optimal_k/./test_prgrm1.sh 1> {output.csv} 2> {output.stderr}")

rule kmergenie:
    input: "/Users/ksahlin/_tmp/Optimal_k/testdata_optimal_k/{dataset}.fa"
    output: csv=OUTBASE+"kmergenie_{dataset}.dat", stderr=OUTBASE+"kmergenie_{dataset}.stderr"
    run:
        for i in range(10):
            pass
        shell(" {GNUTIME} /Users/ksahlin/_tmp/Optimal_k/./test_prgrm2.sh 1> {output.csv} 2> {output.stderr}")

rule unitiger:
    input: csv=OUTBASE+"{tool}_{dataset}.dat", reads="/Users/ksahlin/_tmp/Optimal_k/testdata_optimal_k/{dataset}.fa"
    output: stats=OUTBASE+"{tool}_{dataset}.stdout", contigs=OUTBASE+"{tool}_{dataset}.contigs.fa"
    run:
        print(re.search("optimal_k", input.csv))
        print(re.search("kmergenie", input.csv))
        print(input.csv)
        if re.search("optimal_k", input.csv):
            k,a = get_optimal_k_params()
        elif re.search("kmergenie", input.csv):
            k,a = get_kmer_genie_params()
        else:
            pass #print(shell("{tool}"))

        shell("echo {{input.reads}} {0} {1}".format(k,a)) 
        shell("unitiger {{input.reads}} {0} {1}".format(k,a))    
