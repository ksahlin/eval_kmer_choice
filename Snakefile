 
configfile: "config_uppmax.json"

STDERRSTRING="""
    Command being timed: "optimal-k -r /home/kris/Work/optimal_k/kmergenie/input/spruce_subset.cfg -b /proj/b2013169/nobackup/optimal_k/index -o /tmp/optimal_k_subset.csv"
    User time (seconds): 180706.80
    System time (seconds): 1885.94
    Percent of CPU this job got: 732%
    Elapsed (wall clock) time (h:mm:ss or m:ss): 6:55:33
    Average shared text size (kbytes): 0
    Average unshared data size (kbytes): 0
    Average stack size (kbytes): 0
    Average total size (kbytes): 0
    Maximum resident set size (kbytes): 423691152
    Average resident set size (kbytes): 0
    Major (requiring I/O) page faults: 6
    Minor (reclaiming a frame) page faults: 583402290
    Voluntary context switches: 1850827
    Involuntary context switches: 1591152
    Swaps: 0
    File system inputs: 5608
    File system outputs: 53109072
    Socket messages sent: 0
    Socket messages received: 0
    Signals delivered: 0
    Page size (bytes): 4096
    Exit status: 0
"""
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

####################################################
########## standard python functions ###############
####################################################

import re
import os
def get_kmer_genie_params(csv_file_path):
    max_genomic_kmers = 0
    best_k = 0
    best_a = 0
    lines = open(csv_file_path, 'r').readlines()
    for line in lines[1:]: # removed header
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
    lines = open(csv_file_path, 'r').readlines()
    for line in lines[1:]: # removed header
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

def  parse_gnu_time(stderr_file):
    lines = open(stderr_file, 'r').readlines()

    for l in lines:
        usertime_match =  re.search('User time \(seconds\): [\d.]+', l)
        wct_match = re.search('Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): [\d.:]+', l) 
        mem_match = re.search('Maximum resident set size \(kbytes\): [\d.:]+', l) 
        if usertime_match:
            usertime = float(usertime_match.group().split(':')[1].strip())
        if wct_match:
            wallclocktime = wct_match.group().split()[7]
        if mem_match:
            mem_tmp = int(mem_match.group().split()[5])
            memory_gb = mem_tmp / 4000000.0 

    h,m,s = map(lambda x: int(x), wallclocktime.split(":") )
    tot_wallclock_secs = h*3600 + m*60 + s
    return usertime, tot_wallclock_secs, memory_gb

def myfunc(wildcards):
    input_list_to_performace_latex_table = []
    for dataset in config["DATASETS"]:
        input_list_to_performace_latex_table.append(config["OUTBASE"]+"{0}/kmergenie/default_time_and_mem.txt".format(dataset) )
        input_list_to_performace_latex_table.append(config["OUTBASE"]+"{0}/optimal_k/index_time_and_mem.txt".format(dataset) )
        input_list_to_performace_latex_table.append(config["OUTBASE"]+"{0}/optimal_k/sampling_time_and_mem.txt".format(dataset) )
    return input_list_to_performace_latex_table

###########################################################
###########################################################

rule all:
    input:
        config["OUTBASE"]+"performance_table.tex",        
        config["OUTBASE"]+"quality_table.tex"

rule optimal_k_index:
    input: reads=config["INBASE"]+"{dataset}.cfg"
    output: index=config["OUTBASE"]+"{dataset}/optimal_k/index", 
            stderr=config["OUTBASE"]+"{dataset}/optimal_k/index.stderr", 
            stdout=config["OUTBASE"]+"{dataset}/optimal_k/index.stdout",
            temp_csv=temp("/tmp/{dataset}/index.csv")
    run:
        time=config["GNUTIME"]
        shell(" {time} optimal-k -r {input.reads}  --buildindex {output.index} -o {output.temp_csv} 1> {output.stdout} 2> {output.stderr}")

        ###########
        # for testing on mac:
        #shell("/usr/bin/time -lp touch {output.index} 1> {output.stdout} 2> {output.stderr} " )
        shell("touch {output.index}")
        shell("touch {output.temp_csv}")
        print("{0}".format(STDERRSTRING), file=open(output.stderr, 'w') ) 
        ###########

rule optimal_k_sampling:
    input: reads=config["INBASE"]+"{dataset}.cfg", index=config["OUTBASE"]+"{dataset}/optimal_k/index"
    output: stderr=config["OUTBASE"]+"{dataset}/optimal_k/sampling.stderr", 
            stdout=config["OUTBASE"]+"{dataset}/optimal_k/sampling.stdout",
            csv=config["OUTBASE"]+"{dataset}/optimal_k/sampling.csv",
            best_params=config["OUTBASE"]+"{dataset}/optimal_k/best_params.txt"
    run:
        time=config["GNUTIME"]
        shell(" {time} optimal-k -r {input.reads}  --loadindex {input.index} -a 1 -A 5 -o {output.csv} 1> {output.stdout} 2> {output.stderr}")


        # print("{{0}}\t{{1}}".format(k,a),output.best_params)
        ###########
        # for testing on mac:
        #shell("/usr/bin/time -lp touch {output.csv} 1> {output.stdout} 2> {output.stderr} " )
        shell("echo 2 34 22 99 > {output.csv}")
        print("{0}".format(STDERRSTRING), file=open(output.stderr, 'w') ) 
        ###########

        k,a = get_optimal_k_params(output.csv)
        shell("echo {0} {1} > {{output.best_params}} ".format(k,a))

rule kmergenie:
    input: reads=config["INBASE"]+"{dataset}.cfg"
    output: csv=config["OUTBASE"]+"{dataset}/kmergenie/default.dat",
            stderr=config["OUTBASE"]+"{dataset}/kmergenie/default.stderr", 
            stdout=config["OUTBASE"]+"{dataset}/kmergenie/default.stdout",
            best_params=config["OUTBASE"]+"{dataset}/kmergenie/best_params.txt"
    run:
        env = config["LOAD_PYTHON_ENV"]
        shell("{env}")
        time = config["GNUTIME"]
        out = config["OUTBASE"]
        python = config["PYTHON2"]
        path=config["kmergenie_rules"]["path"]
        shell(" {time} {python} {path}kmergenie -o {out}{wildcards.dataset}/kmergenie/default {input.reads} 1> {output.stdout} 2> {output.stderr}")
        k, a = get_kmer_genie_params(output.csv)
        shell("echo {0} {1} > {{output.best_params}} ".format(k,a))

        # print("{{0}}\t{{1}}".format(k,a),output.best_params)

        ###########
        # for testing on mac:
        # shell("/usr/bin/time -lp touch {output.csv} 1> {output.stdout} 2> {output.stderr} " )
        # shell("echo 1 3 19  > {output.csv}")
        #print("{0}".format(STDERRSTRING), file=open(output.stderr, 'w') ) 
        ###########




rule unitiger:
    input:  reads=config["INBASE"]+"{dataset}.cfg", 
            params=config["OUTBASE"]+"{dataset}/{tool}/best_params.txt" # #rules.kmergenie.output.best_params, rules.optimal_k_sampling.output.best_params,
    output: stdout=config["OUTBASE"]+"{dataset}/{tool}/unitiger.stdout",
            stderr=config["OUTBASE"]+"{dataset}/{tool}/unitiger.stderr",
            contigs=config["OUTBASE"]+"{dataset}/{tool}/unitiger.fa"
    run:
        time=config["GNUTIME"]
        k,a = get_k_and_a_for_assembler(input.params)
        shell("{time} unitiger -r {{input.reads}} -o {{output.contigs}} -k {0} -a {1} 1> {{output.stdout}} 2> {{output.stderr}}".format(k,a))   

        ###########
        # for testing on mac:
        #shell("/usr/bin/time -lp touch {output.contigs} 1> {output.stdout} 2> {output.stderr} " )
        print("{0}\n{1}\n{2}\n{3}\n{4}\n{5}".format('>ctg1','ACGT','>ctg2','AA','>ctg3','C'), file=open(output.contigs, 'w') ) 
        print("{0}".format(STDERRSTRING), file=open(output.stderr, 'w') ) 

        ###########

rule QUAST:
    input: contigs=config["OUTBASE"]+"{dataset}/{tool}/unitiger.fa"
    output: results=config["OUTBASE"]+"{dataset}/{tool}/QUAST/report.txt",
            nice_format=config["OUTBASE"]+"{dataset}/{tool}/result_metrics.csv"
    run:
        out=config['OUTBASE']
        shell("mkdir -p {out}/{wildcards.dataset}/{wildcards.tool}/QUAST/") 
        # shell("/Users/ksahlin/_tmp/Optimal_k/./test_prgrm2.sh < {input.contigs} > {output.results}") 

        print("{0}".format(QUASTSTRING), file=open(output.results, 'w') )
        misassm, N50, NA50, tot_length = parse_quast(output.results)
        e_size = get_esize(input.contigs)
        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(wildcards.dataset, wildcards.tool, e_size, tot_length, N50, misassm,  NA50), file=open(output.nice_format, 'w'))    
        ###########
        # for testing on mac:
        #shell(" tou {output.results} " )
        ###########




rule time_and_mem:
    input:  stderr=config["OUTBASE"]+"{dataset}/{tool}/{method}.stderr" #rules.optimal_k_index.output.stderr, rules.optimal_k_sampling.output.stderr,rules.kmergenie.output.stderr #,
    output: outfile=config["OUTBASE"]+"{dataset}/{tool}/{method}_time_and_mem.txt"
    run:
        usertime, wallclocktime, memory_gb =  parse_gnu_time(input.stderr)
        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(wildcards.dataset, wildcards.tool, wildcards.method, usertime, wallclocktime, memory_gb), file=open(output.outfile, 'w') )
        #shell("touch {config["OUTBASE"]}{wildcards.dataset}/{wildcards.tool}/{wildcards.method}_time_and_mem.txt ") 
        

rule performace_latex_table:
    input: files= myfunc  #rules.kmergenie.output.stderr, rules.optimal_k_index.output.stderr, rules.optimal_k_sampling.output.stderr
    output: table=config["OUTBASE"]+"performance_table.tex"
    run:
        table_file = open(output.table, 'w')
        print("{0} & {1} & {2} & {3} & {4} & {5} \\\ \hline".format('organism', 'tool','method', 'wall clock time', 'user time', 'peak memory'), file=table_file)
        for file_ in input.files:
            line=open(file_,'r').readlines()[0]
            #print("{0} & {1} & {2} & {3} & {4} & {5} \\\ \hline".format(*line.strip().split()))
            print("{0} & {1} & {2} & {3} & {4} & {5} \\\ \hline".format(*line.strip().split()), file=table_file)

rule quality_latex_table:
    input: expand(config["OUTBASE"]+"{dataset}/{tool}/result_metrics.csv",  dataset=config["DATASETS"], tool=config["TOOLS"])
    output: table=config["OUTBASE"]+"quality_table.tex"
    run:
        table_file = open(output.table, 'w')
        print("{0} & {1} & {2} & {3} & {4} & {5} & {6} \\\ \hline".format('organism', 'tool', 'E-size', 'genome size', 'N50', 'misassmblies', 'NA50'), file=table_file)
        for file_ in input:
            line=open(file_,'r').readlines()[0]
            #print("{0} & {1} & {2} & {3} & {4} & {5} & {6} \\\ \hline".format(*line.strip().split()))
            print("{0} & {1} & {2} & {3} & {4} & {5} & {6} \\\ \hline".format(*line.strip().split()), file=table_file)




# rule clean:
#     input:
#     output:
#     run:

# rule test:
#     input:
#     output:
#     run: