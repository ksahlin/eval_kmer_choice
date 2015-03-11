DATASETS='staph rhodo plasm hs14 spruce'.split()
TOOLS='optimal_k kmergenie'.split()
INBASE='/home/kris/Work/optimal_k/config/'
OUTBASE='/proj/b2013169/private/data/optimal_k/'

PYTHON2="/usr/bin/python2.7"
GNUTIME="/usr/bin/time -lp" # "/usr/bin/time -v" if Linux

configfile: "config.json"

#####################################
# standard python functions

import re
import os
def get_kmer_genie_params():
    pass
    return 31,2

def get_optimal_k_params():
    pass
    return 22,3

def get_memory_and_runtime():
    pass


#####################################


rule all:
    input:
        OUTBASE+"performance_table.tex",        
        OUTBASE+"quality_table.tex"

rule optimal_k_index:
    input: reads=INBASE+"{dataset}.cfg"
    output: index=OUTBASE+"optimal_k/{dataset}/index", 
            stderr=OUTBASE+"optimal_k/{dataset}/index.stderr", 
            stdout=OUTBASE+"optimal_k/{dataset}/index.stdout",
            temp_csv=temp("/tmp/optimal_k/{dataset}.csv")
    run:
        shell(" {GNUTIME} optimal-k -r {input.reads}  --buildindex {output.index} -k 2 -K 1 -o {output.temp_csv} 1> {output.stdout} 2> {output.stderr}")

rule optimal_k_sampling:
    input: reads=INBASE+"{dataset}.cfg", index=OUTBASE+"optimal_k/{dataset}/index"
    output: stderr=OUTBASE+"optimal_k/{dataset}/sampling.stderr", 
            stdout=OUTBASE+"optimal_k/{dataset}/sampling.stdout",
            results=OUTBASE+"optimal_k/{dataset}/sampling.csv",
            best_params=OUTBASE+"optimal_k/{dataset}/best_params.txt"
    run:
        shell(" {GNUTIME} optimal-k -r {input.reads}  --readindex {input.index} -a 1 -A 5 -o {output.results} 1> {output.stdout} 2> {output.stderr}")
        for result_file in OUTBASE+"optimal_k/{dataset}":
            k,a = parse_file_here()

        print("{{0}}\t{{1}}".format(k,a),output.best_params)

rule kmergenie:
    input: reads=INBASE+"{dataset}.cfg"
    output: csv=OUTBASE+"kmergenie/{dataset}/kmergenie.dat",
            stderr=OUTBASE+"kmergenie/{dataset}/kmergenie.stderr", 
            stdout=OUTBASE+"kmergenie/{dataset}/kmergenie.stdout",
            best_params=OUTBASE+"kmergenie/{dataset}/best_params.txt"
    run:
        shell(" {GNUTIME} kmergenie -o {OUTBASE}kmergenie/{dataset} 1> {output.stdout} 2> {output.stderr}")
        for result_file in OUTBASE+"kmergenie/{dataset}":
            k,a = parse_file_here()

        print("{{0}}\t{{1}}".format(k,a),output.best_params)

rule unitiger:
    input:  csv=OUTBASE+"{tool}/{dataset}/best_params.txt", 
            reads=INBASE+"{dataset}.cfg"
    output: stdout=OUTBASE+"{tool}/{dataset}.unitiger.stdout",
            stderr=OUTBASE+"{tool}/{dataset}.unitiger.stdout",
            contigs=OUTBASE+"{tool}/{dataset}.contigs.fa"
    run:

        print(re.search("optimal_k", input.csv))
        print(re.search("kmergenie", input.csv))
        print(input.csv)
        if re.search("optimal_k", input.csv):
            k,a = get_optimal_k_params()
        elif re.search("kmergenie", input.csv):
            k,a = get_kmer_genie_params()
        else:
            print(shell("{tool}"))

        shell("echo {{input.reads}} {0} {1}".format(k,a)) 
        shell("{config[unitiger_rules][load_env]}")
        shell("unitiger {{input.reads}} {0} {1} 1> {{output.stats}}".format(k,a))   


rule QUAST:
    input: contigs=OUTBASE+"{tool}/{dataset}.contigs.fa"
    output: results=OUTBASE+"{tool}/{dataset}/report.txt"
    run:
        shell("mkdir -p {output.folder}") 
        shell("/Users/ksahlin/_tmp/Optimal_k/./test_prgrm2.sh < {input.stats} > {output.folder}/results.txt")   


rule time_and_mem:
    input: stderr=OUTBASE+"{tool}_{dataset}.stderr"
    output: folder=OUTBASE+"{tool}_{dataset}/time_and_mem.txt"
    run:
        shell("mkdir -p {output.folder}") 
        shell("touch {output.folder}")   

rule performace_latex_table:
    input: expand(OUTBASE+"{tool}_{dataset}/time_and_mem.txt", tool=TOOLS, dataset=DATASETS)
    output: table=OUTBASE+"performance_table.tex"
    run:
        shell("touch {output.table}") 


rule quality_latex_table:
    input: expand(OUTBASE+"{tool}_{dataset}/report.txt", tool=TOOLS, dataset=DATASETS)
    output: table=OUTBASE+"quality_table.tex"
    run:
        shell("touch {output.table}") 
