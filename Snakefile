DATASETS='staph rhodo plasm hs14 spruce'.split()
TOOLS='optimal_k kmergenie'.split()
METHODS='sampling index kmergenie'.split()
INBASE='/Users/ksahlin/_tmp/Optimal_k/' # '/home/kris/Work/optimal_k/config/' # local testing: /Users/ksahlin/_tmp/Optimal_k/
OUTBASE='/Users/ksahlin/_tmp/Optimal_k/OUT/' # '/proj/b2013169/private/data/optimal_k/' # local testing: /Users/ksahlin/_tmp/Optimal_k/OUT/

PYTHON2="/usr/bin/python2.7"
GNUTIME= "/usr/bin/time -lp" #"/usr/bin/time -v" #if Mac 

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
    output: index=OUTBASE+"{dataset}/optimal_k/index", 
            stderr=OUTBASE+"{dataset}/optimal_k/index.stderr", 
            stdout=OUTBASE+"{dataset}/optimal_k/index.stdout",
            temp_csv=temp("/tmp/{dataset}/index.csv")
    run:
        print('hello')
        #shell(" {GNUTIME} optimal-k -r {input.reads}  --buildindex {output.index} -k 2 -K 1 -o {output.temp_csv} 1> {output.stdout} 2> {output.stderr}")

        ###########
        # for testing on mac:
        shell("/usr/bin/time -lp touch {output.index} 1> {output.stdout} 2> {output.stderr} " )
        shell("touch {output.temp_csv}")
        ###########

rule optimal_k_sampling:
    input: reads=INBASE+"{dataset}.cfg", index=OUTBASE+"{dataset}/optimal_k/index"
    output: stderr=OUTBASE+"{dataset}/optimal_k/sampling.stderr", 
            stdout=OUTBASE+"{dataset}/optimal_k/sampling.stdout",
            results=OUTBASE+"{dataset}/optimal_k/sampling.csv",
            best_params=OUTBASE+"{dataset}/optimal_k/best_params.txt"
    run:
        print("hello")
        # shell(" {GNUTIME} optimal-k -r {input.reads}  --readindex {input.index} -a 1 -A 5 -o {output.results} 1> {output.stdout} 2> {output.stderr}")
        # for result_file in OUTBASE+"{dataset}/sampling":
        #     k,a = parse_file_here()

        # print("{{0}}\t{{1}}".format(k,a),output.best_params)

        ###########
        # for testing on mac:
        shell("/usr/bin/time -lp touch {output.results} 1> {output.stdout} 2> {output.stderr} " )
        shell("touch {output.best_params}")
        ###########

rule kmergenie:
    input: reads=INBASE+"{dataset}.cfg"
    output: csv=OUTBASE+"{dataset}/kmergenie/kmergenie.dat",
            stderr=OUTBASE+"{dataset}/kmergenie/kmergenie.stderr", 
            stdout=OUTBASE+"{dataset}/kmergenie/kmergenie.stdout",
            best_params=OUTBASE+"{dataset}/kmergenie/best_params.txt"
    run:
        print("hello")
        # shell(" {GNUTIME} {PYTHON2} kmergenie -o {OUTBASE}{wildcards.dataset}/kmergenie {input.reads} 1> {output.stdout} 2> {output.stderr}")
        # for result_file in OUTBASE+"{dataset}/kmergenie":
        #     k,a = parse_file_here()

        # print("{{0}}\t{{1}}".format(k,a),output.best_params)

        ###########
        # for testing on mac:
        shell("/usr/bin/time -lp touch {output.csv} 1> {output.stdout} 2> {output.stderr} " )
        shell("touch {output.best_params}")
        ###########


rule unitiger:
    input:  csv=OUTBASE+"{dataset}/{tool}/best_params.txt", 
            reads=INBASE+"{dataset}.cfg"
    output: stdout=OUTBASE+"{dataset}/{tool}.unitiger.stdout",
            stderr=OUTBASE+"{dataset}/{tool}.unitiger.stdout",
            contigs=OUTBASE+"{dataset}/{tool}.contigs.fa"
    run:
        print("hello")
        # print(re.search("optimal_k", input.csv))
        # print(re.search("kmergenie", input.csv))
        # print(input.csv)
        # if re.search("optimal_k", input.csv):
        #     k,a = get_optimal_k_params()
        # elif re.search("kmergenie", input.csv):
        #     k,a = get_kmer_genie_params()
        # else:
        #     print(shell("{tool}"))

        # shell("echo {{input.reads}} {0} {1}".format(k,a)) 
        # shell("{config[unitiger_rules][load_env]}")
        # shell("unitiger {{input.reads}} {0} {1} 1> {{output.stats}}".format(k,a))   

        ###########
        # for testing on mac:
        shell("/usr/bin/time -lp touch {output.contigs} 1> {output.stdout} 2> {output.stderr} " )
        ###########

rule QUAST:
    input: contigs=OUTBASE+"{dataset}/{tool}.contigs.fa"
    output: results=OUTBASE+"{dataset}/{tool}/QUAST/report.txt"
    run:
        shell("mkdir -p {OUTBASE}/{wildcards.dataset}/{wildcards.tool}/QUAST/") 
        shell("/Users/ksahlin/_tmp/Optimal_k/./test_prgrm2.sh < {input.contigs} > {output.results}") 
  
        ###########
        # for testing on mac:
        shell(" touch {output.results} " )
        ###########

rule time_and_mem:
    input: stderr=OUTBASE+"{dataset}/{tool}/{method}.stderr"
    output: outfile=OUTBASE+"{dataset}/{tool}/{method}_time_and_mem.txt"
    run:
        shell("touch {OUTBASE}{wildcards.dataset}/{wildcards.tool}/{wildcards.method}_time_and_mem.txt ") 
        shell("echo {output.outfile}")   

        ###########
        # for testing on mac:
        shell(" touch {output.outfile} " )
        ###########

rule performace_latex_table:
    input: files="/Users/ksahlin/_tmp/Optimal_k/OUT/hs14/kmergenie/kmergenie_time_and_mem.txt" #expand(OUTBASE+"{dataset}/{tool}/{method}_time_and_mem.txt ", dataset=DATASETS, tool=TOOLS,method=METHODS)
    output: table=OUTBASE+"performance_table.tex"
    run:
        #for file in input.files:
        #    print("{0}{1}".format(shell("echo {input.file}"),'lol'))
        shell("touch {output.table}") 

        ###########
        # for testing on mac:
        shell(" touch {output.table}  " )
        ###########

rule quality_latex_table:
    input: expand(OUTBASE+"{dataset}/{tool}/QUAST/report.txt",  dataset=DATASETS, tool=TOOLS)
    output: table=OUTBASE+"quality_table.tex"
    run:
        shell("touch {output.table}") 

        ###########
        # for testing on mac:
        shell(" touch {output.table}  " )
        ###########
