DATASETS='staph rhodo plasm hs14 spruce'.split()
TOOLS='optimal_k kmergenie'.split()
METHODS='sampling index kmergenie'.split()
INBASE='/Users/ksahlin/_tmp/Optimal_k/' # '/home/kris/Work/optimal_k/config/' # local testing: /Users/ksahlin/_tmp/Optimal_k/
OUTBASE='/Users/ksahlin/_tmp/Optimal_k/OUT/' # '/proj/b2013169/private/data/optimal_k/' # local testing: /Users/ksahlin/_tmp/Optimal_k/OUT/

PYTHON2="/usr/bin/python2.7"
GNUTIME= "/usr/bin/time -lp" #"/usr/bin/time -v" #if Mac 

configfile: "config.json"

QUASTSTRING= """
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                     contigs broken  contigs   
# contigs (>= 0 bp)          9318            9318      
# contigs (>= 1000 bp)       5587            5587      
Total length (>= 0 bp)       17827917        17827917  
Total length (>= 1000 bp)    15832094        15832094  
# contigs                    7430            7430      
Largest contig               33090           33090     
Total length                 17194249        17194249  
Reference length             23328019        23328019  
GC (%)                       20.51           20.51     
Reference GC (%)             19.34           19.34     
N50                          3109            3109      
NG50                         2157            2157      
N75                          1822            1822      
L50                          1651            1651      
LG50                         2842            2842      
L75                          3462            3462      
# misassemblies              0               0         
# misassembled contigs       0               0         
Misassembled contigs length  0               0         
# local misassemblies        0               0         
# unaligned contigs          0 + 0 part      0 + 0 part
Unaligned length             0               0         
Genome fraction (%)          73.390          73.390    
Duplication ratio            1.004           1.004     
# N's per 100 kbp            0.00            0.00      
# mismatches per 100 kbp     0.05            0.05      
# indels per 100 kbp         0.28            0.28      
Largest alignment            33090           33090     
NA50                         3109            3109      
NGA50                        2157            2157      
NA75                         1822            1822      
LA50                         1651            1651      
LGA50                        2842            2842      
LA75                         3462            3462    

"""

#####################################
# standard python functions

import re
import os
def get_kmer_genie_params(csv_file_path):
    max_genomic_kmers = 0
    best_k = 0
    best_a = 0
    for line in open(csv_file_path, 'r'):
        k,nr_genomic_kmers, a = map(lambda x: int(x), line.strip().split())
        if nr_genomic_kmers > max_genomic_kmers:
            max_genomic_kmers = nr_genomic_kmers
            best_k = k
            best_a = a
    return best_k, best_a

def get_optimal_k_params(csv_file_path):
    max_genomic_kmers = 0
    best_k = 0
    best_a = 0
    for line in open(csv_file_path, 'r'):
        vals = line.strip().split()
        k,nr_genomic_kmers, a = int(vals[0]), int(vals[1]), int(vals[2]) 
        if nr_genomic_kmers > max_genomic_kmers:
            max_genomic_kmers = nr_genomic_kmers
            best_k = k
            best_a = a
    return best_k, best_a

def get_k_and_a_for_assembler(csv_file_path):
    line = open(csv_file_path, 'r').readlines()[0]
    k, a = map(lambda x: int(x), line.strip().split())
    return k, a

def get_esize(fastafile):
    ctgs = []
    contig = ''
    for line in open(fastafile, 'r'):
        if line[0] == '>':
            if contig:
                ctgs.append(contig)
            contig = ''
        else:
            contig += line.strip()
    ctgs.append(contig)
    lengths = map(lambda x: len(x), ctgs)
    square_lengths = map(lambda x: len(x)**2, ctgs)
    e_size = sum(square_lengths) / float(sum(lengths))

    return(e_size)

def parse_quast(quast_report):
    lines = open(quast_report, 'r').readlines()
    for line in lines:
        if re.match(r"Total length \(\>\= 0 bp\)",line):
            genome_size = int(line.strip().split()[5])
        if re.match(r"N50",line):
            N50 = int(line.strip().split()[1])
        if re.match(r"NA50",line):
            NA50 = int(line.strip().split()[1]) 
        if re.match(r"# misassemblies",line):
            large_misassm = int(line.strip().split()[2]) 
        if re.match(r"# local misassemblies",line):
            local_misassm = int(line.strip().split()[3]) 

    misassm = large_misassm + local_misassm
    return(misassm, N50, NA50, genome_size)

def get_memory_and_runtime(stderr_file):
    # lines = open(stderr_file, 'r').readlines()
    # for line in lines:
    #     if re.match(r"Total length \(\>\= 0 bp\)",line):
    #         genome_size = int(line.strip().split()[5])
    #     if re.match(r"N50",line):
    #         N50 = int(line.strip().split()[1])
    #     if re.match(r"NA50",line):
    #         NA50 = int(line.strip().split()[1]) 
    #     if re.match(r"# misassemblies",line):
    #         large_misassm = int(line.strip().split()[2]) 
    #     if re.match(r"# local misassemblies",line):
    #         local_misassm = int(line.strip().split()[3]) 

    # misassm = large_misassm + local_misassm
    # return(max_gb_used,user_time,sys_time, real_time)
    return 0,0,0,0

def myfunc(wildcards):
    input_list_to_performace_latex_table = []
    for dataset in DATASETS:
        input_list_to_performace_latex_table.append(OUTBASE+"{0}/kmergenie/kmergenie_time_and_mem.txt".format(dataset) )
        input_list_to_performace_latex_table.append(OUTBASE+"{0}/optimal_k/index_time_and_mem.txt".format(dataset) )
        input_list_to_performace_latex_table.append(OUTBASE+"{0}/optimal_k/sampling_time_and_mem.txt".format(dataset) )
    return input_list_to_performace_latex_table

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
        pass #print('hello')
        #shell(" {GNUTIME} optimal-k -r {input.reads}  --buildindex {output.index} -o {output.temp_csv} 1> {output.stdout} 2> {output.stderr}")

        ###########
        # for testing on mac:
        shell("/usr/bin/time -lp touch {output.index} 1> {output.stdout} 2> {output.stderr} " )
        shell("touch {output.temp_csv}")
        ###########

rule optimal_k_sampling:
    input: reads=INBASE+"{dataset}.cfg", index=OUTBASE+"{dataset}/optimal_k/index"
    output: stderr=OUTBASE+"{dataset}/optimal_k/sampling.stderr", 
            stdout=OUTBASE+"{dataset}/optimal_k/sampling.stdout",
            csv=OUTBASE+"{dataset}/optimal_k/sampling.csv",
            best_params=OUTBASE+"{dataset}/optimal_k/best_params.txt"
    run:
        pass #print("hello")
        # shell(" {GNUTIME} optimal-k -r {input.reads}  --loadindex {input.index} -a 1 -A 5 -o {output.csv} 1> {output.stdout} 2> {output.stderr}")
        # for result_file in OUTBASE+"{dataset}/sampling":
        #     k,a = parse_file_here()

        # print("{{0}}\t{{1}}".format(k,a),output.best_params)
        ###########
        # for testing on mac:
        shell("/usr/bin/time -lp touch {output.csv} 1> {output.stdout} 2> {output.stderr} " )
        shell("echo 2 34 22 99 > {output.csv}")
        k,a = get_optimal_k_params(output.csv)
        shell("echo {0} {1} > {{output.best_params}} ".format(k,a))
        ###########

rule kmergenie:
    input: reads=INBASE+"{dataset}.cfg"
    output: csv=OUTBASE+"{dataset}/kmergenie/kmergenie.dat",
            stderr=OUTBASE+"{dataset}/kmergenie/kmergenie.stderr", 
            stdout=OUTBASE+"{dataset}/kmergenie/kmergenie.stdout",
            best_params=OUTBASE+"{dataset}/kmergenie/best_params.txt"
    run:
        pass #print("hello")
        # shell(" {GNUTIME} {PYTHON2} kmergenie -o {OUTBASE}{wildcards.dataset}/kmergenie {input.reads} 1> {output.stdout} 2> {output.stderr}")
        # for result_file in OUTBASE+"{dataset}/kmergenie":
        #     k,a = parse_file_here()

        # print("{{0}}\t{{1}}".format(k,a),output.best_params)

        ###########
        # for testing on mac:
        shell("/usr/bin/time -lp touch {output.csv} 1> {output.stdout} 2> {output.stderr} " )
        shell("echo 1 3 19  > {output.csv}")
        ###########
        k, a = get_kmer_genie_params(output.csv)
        shell("echo {0} {1} > {{output.best_params}} ".format(k,a))


rule unitiger:
    input:  reads=INBASE+"{dataset}.cfg", 
            params=OUTBASE+"{dataset}/{tool}/best_params.txt" # #rules.kmergenie.output.best_params, rules.optimal_k_sampling.output.best_params,
    output: stdout=OUTBASE+"{dataset}/{tool}.unitiger.stdout",
            stderr=OUTBASE+"{dataset}/{tool}.unitiger.stderr",
            contigs=OUTBASE+"{dataset}/{tool}.unitiger.fa"
    run:
        k,a = get_k_and_a_for_assembler(input.params)

        # shell("unitiger {{input.reads}} {0} {1} 1> {{output.stats}}".format(k,a))   

        ###########
        # for testing on mac:
        shell("/usr/bin/time -lp touch {output.contigs} 1> {output.stdout} 2> {output.stderr} " )
        print("{0}\n{1}\n{2}\n{3}\n{4}\n{5}".format('>ctg1','ACGT','>ctg2','AA','>ctg3','C'), file=open(output.contigs, 'w') ) 
        ###########

rule QUAST:
    input: contigs=OUTBASE+"{dataset}/{tool}.unitiger.fa"
    output: results=OUTBASE+"{dataset}/{tool}/QUAST/report.txt",
            nice_format=OUTBASE+"{dataset}/{tool}/result_metrics.csv"
    run:
        shell("mkdir -p {OUTBASE}/{wildcards.dataset}/{wildcards.tool}/QUAST/") 
        # shell("/Users/ksahlin/_tmp/Optimal_k/./test_prgrm2.sh < {input.contigs} > {output.results}") 

        print("{0}".format(QUASTSTRING), file=open(output.results, 'w') )
        misassm, N50, NA50, tot_length = parse_quast(output.results)
        e_size = get_esize(input.contigs)
        print("{0}\n{1}\n{2}\n{3}\n{4}".format(misassm, N50, NA50, tot_length, e_size), file=open(output.nice_format, 'w'))    
        ###########
        # for testing on mac:
        #shell(" tou {output.results} " )
        ###########




rule time_and_mem:
    input:  stderr=OUTBASE+"{dataset}/{tool}/{method}.stderr" #rules.optimal_k_index.output.stderr, rules.optimal_k_sampling.output.stderr,rules.kmergenie.output.stderr #,
    output: outfile=OUTBASE+"{dataset}/{tool}/{method}_time_and_mem.txt"
    run:
        max_gb_used,user_time,sys_time, real_time = get_memory_and_runtime(input.stderr)
        #shell("touch {OUTBASE}{wildcards.dataset}/{wildcards.tool}/{wildcards.method}_time_and_mem.txt ") 
        
        ###########
        # for testing on mac:
        shell(" touch {output.outfile} " )
        ###########

rule performace_latex_table:
    input: files= myfunc  #rules.kmergenie.output.stderr, rules.optimal_k_index.output.stderr, rules.optimal_k_sampling.output.stderr
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
    input: expand(OUTBASE+"{dataset}/{tool}/result_metrics.csv",  dataset=DATASETS, tool=TOOLS)
    output: table=OUTBASE+"quality_table.tex"
    run:
        shell("touch {output.table}") 

        ###########
        # for testing on mac:
        shell(" touch {output.table}  " )
        ###########