"""
Submit this job on uppmax as:
    snakemake --debug --keep-going -j 999 --cluster "sbatch -A {params.account} -p {params.partition} -n {params.n}  -t {params.runtime} -C {params.memsize} -J {params.jobname} --mail-type={params.mail_type} --mail-user={params.mail}"
"""
configfile: "config_uppmax.json"
KMERGENIE_VERSION = str(os.path.getmtime(config["kmergenie_rules"]["path"]))
OPTIMAL_K_VERSION = str(os.path.getmtime(config["optimal_k_rules"]["path"]))

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

OPTIMAL_K_CSV =  """k,a,nr_nodes,nr_edges,avg_internal_nodes,avg_length_unitigs,est_sample_size,nr_unitigs,e_size
15,3,437352359,.,0.00358669,16.0036,464,1.60754e+09,16
16,3,1216502641,.,0.0109691,17.011,100313,3.75641e+09,17.0002
17,3,2443481585,.,0.0488373,18.0488,35368,5.82534e+09,18.0464
18,3,3768844239,.,0.182801,19.1828,17191,5.99015e+09,19.2521
19,3,4704694437,.,0.546444,20.5464,10865,4.36806e+09,20.6209
20,3,5211940469,.,1.30607,22.3061,9197,2.77122e+09,22.5209
21,3,5665059339,.,2.29156,24.2916,14520,1.95253e+09,27.0185
22,3,6122024152,.,3.3553,26.3553,20489,1.5399e+09,31.1681
23,3,6313373920,.,4.21291,28.2129,25727,1.30167e+09,30.9663
24,3,6567472278,.,4.71688,29.7169,30086,1.22502e+09,34.98
25,3,6785297598,.,4.98792,30.9879,31546,1.20132e+09,34.9378
26,3,7024418450,.,5.40689,32.4069,31720,1.15555e+09,49.4987
27,3,7221929603,.,5.86643,33.8664,32518,1.10085e+09,39.7129
28,3,7447682542,.,6.30231,35.3023,33656,1.06454e+09,40.1304
29,3,7539847126,.,6.91963,36.9196,34303,9.8899e+08,40.6148
30,3,7690333419,.,7.19395,38.1939,36400,9.7111e+08,43.1674
31,3,7787962384,.,7.47512,39.4751,36514,9.48968e+08,44.2661
32,3,7864607027,.,7.86582,40.8658,36798,9.12671e+08,52.7146
33,3,7986105483,.,8.36688,42.3669,37677,8.72602e+08,51.3415
34,3,8169930991,.,8.93188,43.9319,38798,8.40919e+08,52.0367
35,3,8219392645,.,9.37237,45.3724,39623,8.07662e+08,51.1146
36,3,8231085100,.,9.91595,46.9159,40593,7.67035e+08,53.4829
37,3,8268502254,.,10.1374,48.1374,42050,7.53975e+08,55.4124
38,3,8289915001,.,10.7897,49.7897,42059,7.10986e+08,66.9164
39,3,8385132393,.,11.0446,51.0446,43856,7.0225e+08,58.4082
40,3,8258890510,.,11.0619,52.0619,43627,6.90701e+08,60.505
41,3,8254125823,.,11.7644,53.7644,43571,6.50215e+08,61.9607
42,3,8224185878,.,11.9208,54.9208,45466,6.38854e+08,58.88
43,3,8202666475,.,12.5658,56.5658,45425,6.05024e+08,70.5345
44,3,8145177434,.,13.0641,58.0641,47082,5.78001e+08,64.6231
45,3,8088032675,.,13.1284,59.1284,48355,5.71473e+08,74.3066
46,3,8037344669,.,13.4943,60.4943,47954,5.52152e+08,74.7161
47,3,7938712600,.,13.6505,61.6505,48654,5.39269e+08,66.2646
48,3,7821017762,.,13.6594,62.6594,48811,5.29512e+08,73.9954
49,3,7713731781,.,13.6663,63.6663,48684,5.21706e+08,75.5299
50,3,7595217688,.,14.0146,65.0146,48370,5.00709e+08,74.6965
51,3,7503177547,.,14.4534,66.4534,49321,4.7952e+08,76.5618
52,3,7347297344,.,14.2832,67.2832,50376,4.74276e+08,77.4555
53,3,7194314832,.,14.1638,68.1638,49784,4.66927e+08,78.3002
54,3,7033237869,.,13.8036,68.8036,49402,4.66858e+08,79.1174
55,3,6857146979,.,14.4255,70.4255,48237,4.36629e+08,75.7662
56,3,6675288775,.,14.68,71.68,50350,4.17033e+08,85.6819
57,3,6542615072,.,14.0498,72.0498,51422,4.24712e+08,82.6553
58,3,6329028753,.,14.3668,73.3668,49197,4.01944e+08,79.7925
59,3,6157790792,.,14.1499,74.1499,50644,3.95916e+08,79.683
60,3,5975523577,.,13.7698,74.7698,50039,3.93178e+08,79.587
61,3,5794708388,.,13.9073,75.9073,48998,3.77565e+08,79.5925
62,3,5579471311,.,13.748,76.748,49591,3.6622e+08,80.8042
63,3,5378535529,.,13.4555,77.4555,49645,3.59955e+08,83.5624
64,3,5175042627,.,12.9112,77.9112,48995,3.59008e+08,82.7664
65,3,4973971800,.,12.7409,78.7409,47601,3.48534e+08,82.9562
66,3,4750674167,.,12.8302,79.8302,47473,3.2996e+08,84.3142
67,3,4516125031,.,12.3665,80.3665,48509,3.23811e+08,82.5627
68,3,4313067214,.,11.9531,80.9531,47748,3.18811e+08,85.0652
69,3,4098827816,.,11.5729,81.5729,46786,3.11285e+08,84.5252
70,3,3898831540,.,11.1252,82.1252,46169,3.06269e+08,86.0131
71,3,3671677869,.,11.2897,83.2897,45144,2.83891e+08,86.1643
72,3,3453354418,.,10.8683,83.8683,46806,2.756e+08,85.271
73,3,3242878142,.,10.1872,84.1872,46240,2.737e+08,86.6188
74,3,3035289562,.,9.85378,84.8538,44564,2.63183e+08,86.3265
75,3,2824665810,.,9.39888,85.3989,44281,2.54992e+08,86.4971
76,3,2612059392,.,9.28966,86.2897,43570,2.37719e+08,87.5998
77,3,2419349767,.,8.84682,86.8468,44462,2.29105e+08,88.1281
78,3,2219116990,.,8.31198,87.312,43760,2.21631e+08,89.1473
79,3,2027339366,.,7.81151,87.8115,42780,2.12786e+08,89.207
80,3,1842703978,.,7.38636,88.3864,42002,2.03048e+08,89.2883
81,3,1668872886,.,7.30262,89.3026,41340,1.84382e+08,90.6264
82,3,1491715339,.,6.81901,89.819,42596,1.74221e+08,91.3292
83,3,1326779450,.,6.28336,90.2834,41962,1.6536e+08,90.7131
84,3,1169431071,.,5.84098,90.841,40922,1.54643e+08,91.7433
85,3,1025375006,.,5.46153,91.4615,40247,1.42995e+08,91.9885"""

####################################################
########## standard python functions ###############
####################################################

import re
import os
import os.path
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
    max_objective = 0
    best_k = 0
    lines = open(csv_file_path, 'r').readlines()
    for line in lines[1:]: # removed header
        vals = line.strip().split(',')
        k,nr_genomic_kmers, e_size = int(vals[0]), int(vals[2]), float(vals[8]) 
        objective = e_size
        if objective > max_objective:
            max_objective = objective
            best_k = k
    return best_k , max_objective

def get_k_and_a_for_assembler(csv_file_path):
    line = open(csv_file_path, 'r').readlines()[0]
    k, a = map(lambda x: int(x), line.strip().split())
    return k, a

# def get_esize(fastafile):
#     ctgs = []
#     contig = ''
#     for line in open(fastafile, 'r'):
#         if line[0] == '>':
#             if contig:
#                 ctgs.append(contig)
#             contig = ''
#         else:
#             contig += line.strip()
#     ctgs.append(contig)
#     lengths = map(lambda x: len(x), ctgs)
#     square_lengths = map(lambda x: len(x)**2, ctgs)
#     e_size = sum(square_lengths) / float(sum(lengths))

#     return(e_size)

def parse_quast(quast_report):
    lines = open(quast_report, 'r').readlines()
    ESIZE_GENOME = '-'
    CORR_ESIZE_GENOME  = '-'

    has_ref = False
    for line in lines:
        if re.match(r"Total length \(\>\= 0 bp\)",line):
            genome_size = int(line.strip().split()[5])
        if re.match(r"NG50",line):
            NG50 = int(line.strip().split()[1])
            has_ref = True
        if re.match(r"NGA50",line):
            NGA50 = int(line.strip().split()[1]) 
            has_ref = True
        if re.match(r"# misassemblies",line):
            large_misassm = int(line.strip().split()[2]) 
            has_ref = True
        if re.match(r"# local misassemblies",line):
            local_misassm = int(line.strip().split()[3]) 
            has_ref = True
        
        # for spruce
        if re.match(r"N50",line):
            N50 = int(line.strip().split()[1])

        if re.match(r"ESIZE_GENOME",line):
            ESIZE_GENOME = float(line.strip().split()[1]) 
        if re.match(r"CORR_ESIZE_GENOME",line):
            CORR_ESIZE_GENOME = float(line.strip().split()[1]) 


    if has_ref:
        misassm = large_misassm + local_misassm
        try:
            NGA50
        except UnboundLocalError:
            NGA50 = "."
        try:
            NG50
        except UnboundLocalError:
            NG50 = "."
        return(misassm, NG50, NGA50, genome_size, ESIZE_GENOME, CORR_ESIZE_GENOME)
    else:
        return('.', N50, '.', genome_size, '-', '-')

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

    vals = list(map(lambda x: float(x), wallclocktime.split(":") ))
    if len(vals) == 3:
        h,m,s = vals
        tot_wallclock_secs = h*3600.0 + m*60.0 + s
    elif len(vals) == 2:
        m,s = vals
        tot_wallclock_secs = m*60.0 + s

    return usertime, tot_wallclock_secs, memory_gb

def myfunc(wildcards):
    input_list_to_performace_latex_table = []
    for dataset in config["DATASETS"]:
        input_list_to_performace_latex_table.append(config["OUTBASE"]+"{0}/kmergenie/default_time_and_mem.txt".format(dataset) )
        input_list_to_performace_latex_table.append(config["OUTBASE"]+"{0}/optimal_k/index_time_and_mem.txt".format(dataset) )
        input_list_to_performace_latex_table.append(config["OUTBASE"]+"{0}/optimal_k/sampling_contigs_time_and_mem.txt".format(dataset) )
        input_list_to_performace_latex_table.append(config["OUTBASE"]+"{0}/optimal_k/sampling_unitigs_time_and_mem.txt".format(dataset) )
        input_list_to_performace_latex_table.append(config["OUTBASE"]+"{0}/preqc/preqc_time_and_mem.txt".format(dataset) )
    return input_list_to_performace_latex_table

def assembly_targets(wildcards):
  input_= []

  for dataset in config["DATASETS"]:
    for tool in config["TOOLS"]:
      for assembler in config["ASSEMBLERS"]:
        if tool == 'preqc' and assembler == 'minia_utg':
            continue
        if tool == 'preqc' and assembler == 'unitiger':
            continue
        elif tool == 'optimal_k' and assembler == 'minia_utg':
          input_.append(config["OUTBASE"]+"{0}/{1}/unitigs_{2}.fasta".format(dataset, tool, assembler) )
        elif tool == 'optimal_k' and assembler == 'unitiger':
          input_.append(config["OUTBASE"]+"{0}/{1}/unitigs_{2}.fasta".format(dataset, tool, assembler) )
        elif tool == 'preqc': 
          input_.append(config["OUTBASE"]+"{0}/{1}/default_{2}.fasta".format(dataset, tool, assembler) )
        elif tool == 'kmergenie':
          input_.append(config["OUTBASE"]+"{0}/{1}/default_{2}.fasta".format(dataset, tool, assembler) )
        elif tool == 'optimal_k':
          input_.append(config["OUTBASE"]+"{0}/{1}/contigs_{2}.fasta".format(dataset, tool, assembler) )

  return input_

def prediction_targets(wildcards):
  input_= []

  for dataset in config["DATASETS"]:
    for tool in config["TOOLS"]:
      if tool == 'optimal_k':
        input_.append(config["OUTBASE"]+"{0}/{1}/best_params_unitigs.txt".format(dataset, tool))
        input_.append(config["OUTBASE"]+"{0}/{1}/best_params_contigs.txt".format(dataset, tool))
      else:
        input_.append(config["OUTBASE"]+"{0}/{1}/best_params_default.txt".format(dataset, tool) )

  return input_

def evaluation_targets(wildcards):
  input_= []

  for dataset in config["DATASETS"]:
    for tool in config["TOOLS"]:
      for assembler in config["ASSEMBLERS"]:
        if tool == 'preqc' and assembler == 'minia_utg':
            continue
        if tool == 'preqc' and assembler == 'unitiger':
            continue
        elif tool == 'optimal_k' and assembler == 'minia_utg':
          input_.append(config["OUTBASE"]+"{0}/{1}/result_metrics_unitigs_{2}.csv".format(dataset, tool, assembler) )
        elif tool == 'optimal_k' and assembler == 'unitiger':
          input_.append(config["OUTBASE"]+"{0}/{1}/result_metrics_unitigs_{2}.csv".format(dataset, tool, assembler) )
        elif tool == 'preqc':
          input_.append(config["OUTBASE"]+"{0}/{1}/result_metrics_default_{2}.csv".format(dataset, tool, assembler) )
        elif tool == 'kmergenie':
          input_.append(config["OUTBASE"]+"{0}/{1}/result_metrics_default_{2}.csv".format(dataset, tool, assembler) )
        elif tool == 'optimal_k':
          input_.append(config["OUTBASE"]+"{0}/{1}/result_metrics_contigs_{2}.csv".format(dataset, tool, assembler) )


        #input_.append(config["OUTBASE"]+"{0}/{1}/result_metrics_{2}.csv".format(dataset, tool, assembler) )

  return input_

def sampling_targets(wildcards):
  input_= []

  for dataset in config["DATASETS"]:
    input_.append(config["OUTBASE"]+"{0}/optimal_k/best_params_unitigs.txt".format(dataset) )
    input_.append(config["OUTBASE"]+"{0}/optimal_k/best_params_contigs.txt".format(dataset) )

  return input_

###########################################################
###########################################################

rule all:
    input:
        config["OUTBASE"]+"performance_table.tex",
        expand(config["OUTBASE"]+"quality_table_{assembler}.tex", assembler=config["ASSEMBLERS"])
    params:
        runtime="15:00",
        memsize = "'mem128GB|mem256GB|mem512GB'",
        partition = "core",
        n = "1",
        jobname="all",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"]

# sub build targets

rule sampling:
    input: sampling_targets   
    params: 
        runtime="15:00",
        memsize = "'mem128GB|mem256GB|mem512GB'",
        partition = "core",
        n = "1",
        jobname="sampling",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"]

rule predict:
    input: prediction_targets
    params: 
        runtime="15:00",
        memsize = "mem128GB",
        partition = "core",
        n = "1",
        jobname="all",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"]

rule assemble:
    input: assembly_targets
    params: 
        runtime="15:00",
        memsize = "mem128GB",
        partition = "core",
        n = "1",
        jobname="all",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"]

rule evaluate:
    input: evaluation_targets
    params: 
        runtime="15:00",
        memsize = "mem128GB",
        partition = "core",
        n = "1",
        jobname="all",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"]

#rule all:
#    input: rules.predict.input, rules.assemble.input, rules.evaluate.input,
#        config["OUTBASE"]+"performance_table.tex",
#        expand(config["OUTBASE"]+"quality_table_{assembler}.tex", assembler=config["ASSEMBLERS"])
#    params:
#        runtime="15:00",
#        memsize = "'mem128GB|mem256GB|mem512GB'",
#        partition = "core",
#        n = "1",
#        jobname="all",
#        account=config["SBATCH"]["ACCOUNT"],
#        mail=config["SBATCH"]["MAIL"],
#        mail_type=config["SBATCH"]["MAIL_TYPE"]


rule optimal_k_index:
    input: reads=config["INBASE"]+"{dataset}.cfg"
    output: index=protected(config["OUTBASE"]+"{dataset}/optimal_k/index.rlcsa.array"), 
            stderr=config["OUTBASE"]+"{dataset}/optimal_k/index.stderr", 
            stdout=config["OUTBASE"]+"{dataset}/optimal_k/index.stdout"
    version: OPTIMAL_K_VERSION
    params: 
        runtime = lambda wildcards: config["SBATCH"][wildcards.dataset]["optimalk_index_time"],
        memsize = lambda wildcards: config["SBATCH"][wildcards.dataset]["memsize"],
        partition = lambda wildcards: config["SBATCH"][wildcards.dataset]["partition"],
        n = lambda wildcards: config["SBATCH"][wildcards.dataset]["n"],
        jobname = "{dataset}"+"_optimalk_indexing",
        account = config["SBATCH"]["ACCOUNT"],
        mail = config["SBATCH"]["MAIL"],
        mail_type = config["SBATCH"]["MAIL_TYPE"]
    run:
        time=config["GNUTIME"]
        runtime=config["SBATCH"][wildcards.dataset]["optimalk_index_time"]
        index_path=config["OUTBASE"]+"{0}/optimal_k/index".format(wildcards.dataset)
        stderr=config["OUTBASE"]+"{0}/optimal_k/index_tmp.stderr".format(wildcards.dataset) 
        stdout=config["OUTBASE"]+"{0}/optimal_k/index_tmp.stdout".format(wildcards.dataset)
        shell(" {time} optimal-k -r {input.reads}  --buildindex {index_path} 1> {stdout} 2> {stderr}")
        shell("cp {stdout} {output.stdout} ")
        shell("cp {stderr} {output.stderr} ")

        ###########
        # for testing on mac:
        # shell("touch {output.index}")
        # shell("touch {output.temp_csv}")
        # print("{0}".format(STDERRSTRING), file=open(output.stderr, 'w') ) 
        ###########

rule optimal_k_sampling_contigs:
    input: reads=config["INBASE"]+"{dataset}.cfg", index=config["OUTBASE"]+"{dataset}/optimal_k/index.rlcsa.array"
    output: stderr=config["OUTBASE"]+"{dataset}/optimal_k/sampling_contigs.stderr", 
            stdout=config["OUTBASE"]+"{dataset}/optimal_k/sampling_contigs.stdout",
            complete=config["OUTBASE"]+"{dataset}/optimal_k/sampling_contigs_ok.txt",
            #best_params=config["OUTBASE"]+"{dataset}/optimal_k/best_params.txt"
    version: OPTIMAL_K_VERSION
    params: 
        runtime= lambda wildcards: config["SBATCH"][wildcards.dataset]["optimalk_sample_time"],
        memsize = lambda wildcards: config["SBATCH"][wildcards.dataset]["memsize"],
        partition = lambda wildcards: config["SBATCH"][wildcards.dataset]["partition"],
        n = lambda wildcards: config["SBATCH"][wildcards.dataset]["n"],
        jobname="{dataset}"+"_optimalk_sampling_contigs",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"]
    run:
        time=config["GNUTIME"]
        prefix=config["OUTBASE"]+"{0}/optimal_k/sampling_contigs".format(wildcards.dataset)
        #min_a = config["optimal_k_rules"]["min_abundance"]
        #max_a = config["optimal_k_rules"]["max_abundance"]
        index_path=config["OUTBASE"]+"{0}/optimal_k/index".format(wildcards.dataset)

        shell(" {time} optimal-k -r {input.reads}  --loadindex {index_path} -a 2 -A 12 -o {prefix} 1> {output.stdout} 2> {output.stderr}")
        shell("touch {output.complete}")

rule optimal_k_sampling_unitigs:
    input: reads=config["INBASE"]+"{dataset}.cfg", index=config["OUTBASE"]+"{dataset}/optimal_k/index.rlcsa.array"
    output: stderr=config["OUTBASE"]+"{dataset}/optimal_k/sampling_unitigs.stderr", 
            stdout=config["OUTBASE"]+"{dataset}/optimal_k/sampling_unitigs.stdout",
            complete=config["OUTBASE"]+"{dataset}/optimal_k/sampling_unitigs_ok.txt",
            #best_params=config["OUTBASE"]+"{dataset}/optimal_k/best_params.txt"
    version: OPTIMAL_K_VERSION
    params: 
        runtime= lambda wildcards: config["SBATCH"][wildcards.dataset]["optimalk_sample_time"],
        memsize = lambda wildcards: config["SBATCH"][wildcards.dataset]["memsize"],
        partition = lambda wildcards: config["SBATCH"][wildcards.dataset]["partition"],
        n = lambda wildcards: config["SBATCH"][wildcards.dataset]["n"],
        jobname="{dataset}"+"_optimalk_sampling_unitigs",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"]
    run:
        time=config["GNUTIME"]
        prefix=config["OUTBASE"]+"{0}/optimal_k/sampling_unitigs".format(wildcards.dataset)
        #min_a = config["optimal_k_rules"]["min_abundance"]
        #max_a = config["optimal_k_rules"]["max_abundance"]
        index_path=config["OUTBASE"]+"{0}/optimal_k/index".format(wildcards.dataset)

        shell(" {time} optimal-k -r {input.reads}  --loadindex {index_path} -a 2 -A 12 -u -o {prefix} 1> {output.stdout} 2> {output.stderr}")
        shell("touch {output.complete}")

rule get_best_params:
    input: config["OUTBASE"]+"{dataset}/optimal_k/sampling_{type}_ok.txt"
    output: best_params= config["OUTBASE"]+"{dataset}/optimal_k/best_params_{type}.txt"  #config["OUTBASE"]+"{dataset}/optimal_k/best_params.txt"
    params: 
        runtime="10:00",
        memsize = "'mem128GB|mem256GB|mem512GB'",
        partition = "core",
        n = "1",
        jobname="get_best_params_{type}",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"]
    run:
        prefix=config["OUTBASE"]+"{0}/optimal_k/sampling_{1}".format(wildcards.dataset, wildcards.type)
        min_a = config["optimal_k_rules"]["{0}".format(wildcards.type)]["min_abundance"]
        max_a = config["optimal_k_rules"]["{0}".format(wildcards.type)]["max_abundance"]
        max_objective = 0
        best_k = 0
        best_a = 0
        find_k = config["optimal_k_rules"]["script_path"]+"fit_curve.py"
        python = config["PYTHON2"]
        for abundance in range(int(min_a),int(max_a)+1):
            csv_file = prefix+".a{0}.csv".format(abundance)
            k, objective = list(shell("{python} {find_k} {csv_file}", iterable=True))[0].split()
            k = int(float(k))
            objective = float(objective)

            if objective > max_objective:
                best_k = k
                best_a = abundance
                max_objective = objective

        shell("echo {0} {1} > {{output.best_params}} ".format(best_k,best_a))


rule kmergenie:
    input: reads=config["INBASE"]+"{dataset}.cfg"
    output: csv=config["OUTBASE"]+"{dataset}/kmergenie/default.dat",
            stderr=config["OUTBASE"]+"{dataset}/kmergenie/default.stderr", 
            stdout=config["OUTBASE"]+"{dataset}/kmergenie/default.stdout",
            best_params=config["OUTBASE"]+"{dataset}/kmergenie/best_params_default.txt"
    version: KMERGENIE_VERSION

    params: 
        runtime=lambda wildcards: config["SBATCH"][wildcards.dataset]["kmergenie_time"],
        memsize = lambda wildcards: config["SBATCH"][wildcards.dataset]["small_memsize"],
        partition = lambda wildcards: config["SBATCH"][wildcards.dataset]["small_partition"],
        n = lambda wildcards: config["SBATCH"][wildcards.dataset]["small_n"],
        jobname="{dataset}"+"_kmergenie",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"]
    run:
        env = config["LOAD_PYTHON_ENV"]
        shell("{env}")
        time = config["GNUTIME"]
        out = config["OUTBASE"]
        python = config["PYTHON2"]
        path=config["kmergenie_rules"]["path"]
        if (wildcards.dataset == "spruce" or wildcards.dataset =="hs14"):
            model = "--diploid" 
        else:
            model = ""
        shell(" {time} {python} {path}kmergenie {model} -o {out}{wildcards.dataset}/kmergenie/default {input.reads} 1> {output.stdout} 2> {output.stderr}")
        k, a = get_kmer_genie_params(output.csv)
        shell("echo {0} {1} > {{output.best_params}} ".format(k,a))


rule preqc:
    input: reads=config["INBASE"]+"{dataset}.cfg"
    output: stderr=config["OUTBASE"]+"{dataset}/preqc/preqc.stderr", 
            best_params=config["OUTBASE"]+"{dataset}/preqc/best_params_default.txt"
    params: 
        runtime=lambda wildcards:  config["SBATCH"][wildcards.dataset]["preqc_time"],
        memsize = lambda wildcards: config["SBATCH"][wildcards.dataset]["memsize"],
        partition = lambda wildcards: config["SBATCH"][wildcards.dataset]["partition"],
        n = lambda wildcards: config["SBATCH"][wildcards.dataset]["n"],
        jobname="{dataset}_"+"preQC",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"]
    run:
        pass
        time = config["GNUTIME"]
        stderr=config["OUTBASE"]+"{0}/preqc/sga_pipeline.stderr".format(wildcards.dataset) 
        file1 = list(shell("head -n 1 {input.reads}", iterable=True))[0]
        file2 = list(shell("head -n 2 {input.reads}", iterable=True))[1]
        shell("mkdir -p /tmp/{wildcards.dataset}/preqc/")
        tmp_prefix = "/tmp/{0}/preqc/preqc_genome".format(wildcards.dataset) 
        
        shell("{time} SGA_optk_pipeline {tmp_prefix} {file1} {file2} 2> {stderr}")
        shell("cp {stderr} {output.stderr}")
        best_k = list(shell("SGA_wrapper get_k {tmp_prefix}.preqc",iterable=True))[0]
        shell("echo {0} {1} > {{output.best_params}} ".format(best_k,3))



rule unitiger:
    input:  reads=config["INBASE"]+"{dataset}.cfg", 
            params=config["OUTBASE"]+"{dataset}/{tool}/best_params_{type}.txt" # #rules.kmergenie.output.best_params, rules.optimal_k_sampling.output.best_params,
    output: unitigs=config["OUTBASE"]+"{dataset}/{tool}/{type}_unitiger.fasta"
    params: 
        runtime=lambda wildcards: config["SBATCH"][wildcards.dataset]["unitiger_time"],
        memsize = lambda wildcards: config["SBATCH"][wildcards.dataset]["memsize"],
        partition = lambda wildcards: config["SBATCH"][wildcards.dataset]["partition"],
        n = lambda wildcards: config["SBATCH"][wildcards.dataset]["n"],
        jobname="{dataset}_{tool}"+"_unitiger",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"]
    run:
        time = config["GNUTIME"]
        prefix=config["OUTBASE"]+"{0}/{1}/unitiger".format(wildcards.dataset, wildcards.tool)
        k,a = get_k_and_a_for_assembler(input.params)

        env = config["LOAD_PYTHON_ENV"]
        shell("{env}")
        python = config["PYTHON2"]
        path=config["unitiger_rules"]["path"]
        stderr=config["OUTBASE"]+"{0}/{1}/unitiger.stderr".format(wildcards.dataset, wildcards.tool)
        #shell("{time} {python} {path}Unitiger_wrapper.py -r {input.reads} -o {prefix} -k {k} -K {k} -a {a} -A {a} 2>&1 | tee -a {stderr}")   
        shell("{time} {path}Unitiger -r {input.reads} -o {prefix} -k {k} -a {a} -t 0 2>&1 | tee -a {stderr}")   
        shell("mv {0} {1}".format(prefix+".k{0}.a{1}.unitigs".format(k,a), output.unitigs))

        ###########
        # for testing on mac:

        # print("{0}\n{1}\n{2}\n{3}\n{4}\n{5}".format('>ctg1','CTAGCTCTACGTCACTCACGCCCCGCTTTCTATTGATGGAAGTCGTCTAATTCACTATAACAGCGAATCGGGGCCCCTCAGCCCATATGCTGAGCCCTCCTGTACGTGATCTATACTGGCTTTTAATACAGAAGGCCACCACTA',\
        #     '>ctg2','CCCTAGCACCGTCACTCTATTTTGTACCCTTGAACTTCTCGACATTCTATTTCGGCCAGGCGTACAAACCTGCGGTGATGGGCCTGCTAAACACCACC',\
        #     '>ctg3','CGGCGAAGCTTAGGCGCTTCAAAAGCCAAAACATCAACGATGGTTCGCCGGGGGTGACGCCTTACCTATCTAGCGTGCGTCGCGTGATGCACGCTGGCATTGAGGTGAATTCGGCCTAGGATGCTTAATCAGAGCATGTTCCATCGTTAGCGGCTGCCAGGAAGGTGGTATATCACCTCCGGGGTGGTCAAAAATGGGGTCGC'),\
        #      file=open(prefix+".fasta", 'w') ) 
        # print("{0}".format(STDERRSTRING), file=open(output.stderr, 'w') ) 
        # print("{0}".format(STDERRSTRING), file=open(output.stdout, 'w') ) 
        ###########

rule minia:
    input:  reads=config["INBASE"]+"{dataset}.cfg", 
            params=config["OUTBASE"]+"{dataset}/{tool}/best_params_{type}.txt" # #rules.kmergenie.output.best_params, rules.optimal_k_sampling.output.best_params,
    output: contigs=config["OUTBASE"]+"{dataset}/{tool}/{type}_minia.fasta"
    params: 
        runtime=lambda wildcards:  config["SBATCH"][wildcards.dataset]["minia_time"],
        memsize = lambda wildcards: config["SBATCH"][wildcards.dataset]["small_memsize"],
        partition = lambda wildcards: config["SBATCH"][wildcards.dataset]["partition"],
        n = lambda wildcards: config["SBATCH"][wildcards.dataset]["n"],
        jobname="{dataset}_{tool}"+"_minia",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"]
    run:
        time = config["GNUTIME"]
        prefix=config["OUTBASE"]+"{0}/{1}/{2}_minia".format(wildcards.dataset, wildcards.tool, wildcards.type)
        k,a = get_k_and_a_for_assembler(input.params)
        # stdout=config["OUTBASE"]+"{0}/{1}/minia.stdout".format(wildcards.dataset, wildcards.tool)
        stderr=config["OUTBASE"]+"{0}/{1}/minia.stderr".format(wildcards.dataset, wildcards.tool) 
        shell("{time} minia -in {input.reads} -kmer-size {k} -abundance-min 3 -no-length-cutoff -out {prefix} 2>&1 | tee -a {stderr}")   
        shell("mv {0} {1}".format(prefix+".contigs.fa", output.contigs))

        # clean
        folder = config["OUTBASE"]+"{0}/{1}/".format(wildcards.dataset, wildcards.tool)
        shell("rm -f {folder}*.parts.* ")
        shell("rm -f {folder}*.h5")


rule minia_utg:
    input:  reads=config["INBASE"]+"{dataset}.cfg", 
            params=config["OUTBASE"]+"{dataset}/{tool}/best_params_{type}.txt" 
    output: contigs=config["OUTBASE"]+"{dataset}/{tool}/{type}_minia_utg.fasta"
    params: 
        runtime=lambda wildcards:  config["SBATCH"][wildcards.dataset]["minia_time"],
        memsize = lambda wildcards: config["SBATCH"][wildcards.dataset]["small_memsize"],
        partition = lambda wildcards: config["SBATCH"][wildcards.dataset]["small_partition"],
        n = lambda wildcards: config["SBATCH"][wildcards.dataset]["small_n"],
        jobname="{dataset}_{tool}_"+"_minia_utg",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"]
    run:
        time = config["GNUTIME"]
        prefix=config["OUTBASE"]+"{0}/{1}/{2}_minia_utg".format(wildcards.dataset, wildcards.tool, wildcards.type)
        k,a = get_k_and_a_for_assembler(input.params)
        # stdout=config["OUTBASE"]+"{0}/{1}/minia_utg.stdout".format(wildcards.dataset, wildcards.tool)
        stderr=config["OUTBASE"]+"{0}/{1}/minia_utg.stderr".format(wildcards.dataset, wildcards.tool) 
        shell("{time} minia -in {input.reads} -traversal unitig -starter simple -no-length-cutoff -kmer-size {k} -abundance-min {a} -out {prefix} 2>&1 | tee -a {stderr}")   
        shell("mv {0} {1}".format(prefix+".contigs.fa", output.contigs))

        # clean
        folder = config["OUTBASE"]+"{0}/{1}/".format(wildcards.dataset, wildcards.tool)
        shell("rm -f {folder}*.parts.* ")
        shell("rm -f {folder}*.h5")

rule abyss:
    input:  reads=config["INBASE"]+"{dataset}.cfg", 
            params=config["OUTBASE"]+"{dataset}/{tool}/best_params_{type}.txt" # #rules.kmergenie.output.best_params, rules.optimal_k_sampling.output.best_params,
    output: contigs=config["OUTBASE"]+"{dataset}/{tool}/{type}_abyss.fasta"
    params: 
        runtime=lambda wildcards:  config["SBATCH"][wildcards.dataset]["abyss_time"],
        memsize = lambda wildcards: config["SBATCH"][wildcards.dataset]["memsize"],
        partition = lambda wildcards: config["SBATCH"][wildcards.dataset]["partition"],
        n = lambda wildcards: config["SBATCH"][wildcards.dataset]["n"],
        jobname="{dataset}_{tool}_"+"abyss",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"]
    run:
        time = config["GNUTIME"]
        prefix=config["OUTBASE"]+"{0}/{1}/{2}_abyss".format(wildcards.dataset, wildcards.tool, wildcards.type)
        k,a = get_k_and_a_for_assembler(input.params)
        # stdout=config["OUTBASE"]+"{0}/{1}/abyss.stdout".format(wildcards.dataset, wildcards.tool)
        stderr=config["OUTBASE"]+"{0}/{1}/abyss.stderr".format(wildcards.dataset, wildcards.tool) 
        file1 = list(shell("head -n 1 {input.reads}", iterable=True))[0]
        file2 = list(shell("head -n 2 {input.reads}", iterable=True))[1]
        shell("abyss-pe {prefix} {k} {file1} {file2} --always-make 2>&1 | tee -a {stderr}")  
        # Apperently, abyss-contigs.fa is a symlink to abyss-6.fa, so we 
        # mv abyss-6.fa file instead
        shell("mv {0} {1}".format(prefix+"-6.fa", prefix+'.fasta'))
        # we update the timestamp on the output of rule abyss so the snakemake does not return
        # "Output files ... are older than input files.". Unsure why this happens.
        shell("touch {0} ".format(prefix+'.fasta'))
        abyss_output= prefix+'-*'
        shell("rm {abyss_output}")

rule velvet:
    input:  reads=config["INBASE"]+"{dataset}.cfg", 
            params=config["OUTBASE"]+"{dataset}/{tool}/best_params_{type}.txt" # #rules.kmergenie.output.best_params, rules.optimal_k_sampling.output.best_params,
    output: contigs=config["OUTBASE"]+"{dataset}/{tool}/{type}_velvet.fasta"
    params: 
        runtime=lambda wildcards:  config["SBATCH"][wildcards.dataset]["velvet_time"],
        memsize = lambda wildcards: config["SBATCH"][wildcards.dataset]["memsize"],
        partition = lambda wildcards: config["SBATCH"][wildcards.dataset]["partition"],
        n = lambda wildcards: config["SBATCH"][wildcards.dataset]["n"],
        jobname="{dataset}_{tool}_"+"velvet",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"],
        insertsize= lambda wildcards: config["inserts"][wildcards.dataset]
    run:
        time = config["GNUTIME"]
        prefix=config["OUTBASE"]+"{0}/{1}/velvet_tmp".format(wildcards.dataset, wildcards.tool)
        file1 = list(shell("head -n 1 {input.reads}", iterable=True))[0]
        file2 = list(shell("head -n 2 {input.reads}", iterable=True))[1]
        k,a = get_k_and_a_for_assembler(input.params)

        stdout=config["OUTBASE"]+"{0}/{1}/velvet.stdout".format(wildcards.dataset, wildcards.tool)
        stderr=config["OUTBASE"]+"{0}/{1}/velvet.stderr".format(wildcards.dataset, wildcards.tool) 
        shell("velvet-pe {prefix} {k} {file1} {file2} {wildcards.dataset} {params.insertsize} 1> {stdout} 2> {stderr}")  
        #velvet_out=config["OUTBASE"]+"{0}/{1}/velvet.fasta".format(wildcards.dataset, wildcards.tool)
        shell("mv {0} {1}".format(prefix+"/contigs.fa", output.contigs))
        #shell("rm -r {prefix}")

rule QUAST:
    input: contigs=config["OUTBASE"]+"{dataset}/{tool}/{type}_{assembler}.fasta",
            param=config["OUTBASE"]+"{dataset}/{tool}/best_params_{type}.txt"
    output: #results=config["OUTBASE"]+"{dataset}/{tool}/QUAST/report.txt",
            nice_format=config["OUTBASE"]+"{dataset}/{tool}/result_metrics_{type}_{assembler}.csv"
    log: config["OUTBASE"]+"{dataset}/{tool}/{assembler}/quast.output"
    params: 
        runtime=lambda wildcards: config["SBATCH"][wildcards.dataset]["quast_time"],
        memsize = lambda wildcards: config["SBATCH"][wildcards.dataset]["small_memsize"],
        partition = lambda wildcards: config["SBATCH"][wildcards.dataset]["small_partition"],
        n = lambda wildcards: config["SBATCH"][wildcards.dataset]["small_n"],
        jobname="quast_{dataset}_{tool}_{assembler}",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"]
    run:
        env = config["LOAD_PYTHON_ENV"]
        shell("{env}")
        python = config["PYTHON2"]
        out=config['OUTBASE']
        path=config["quast_rules"]["path"]
        min_contig =  config["quast_rules"]["min_contig"]
        outpath="{0}/{1}/{2}/{3}/QUAST/".format(out, wildcards.dataset,wildcards.tool, wildcards.assembler)
        reference = config["REFERENCES"][wildcards.dataset]
        fail_flag = "/tmp/{0}_{1}_{2}_quastfail.err".format( wildcards.dataset,wildcards.tool, wildcards.assembler)
        if wildcards.dataset == "spruce":
            shell("if {python} {path}quast.py -o {outpath} --min-contig 200 --no-plots {input.contigs} &> {log}; then : nothing; else touch {fail_flag} ")
        elif wildcards.dataset == "hs14":
            shell(" {python} {path}quast.py -o {outpath} --min-contig 100 --no-plots {input.contigs} &> {log} ") 

        else:
            shell(" {python} {path}quast.py -R  {reference} -o {outpath} --min-contig 30 --no-plots {input.contigs} &> {log} ") 

        if os.path.isfile(fail_flag):
            os.remove(fail_flag)
            k,a = get_k_and_a_for_assembler(input.param)
            print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}".format(wildcards.dataset, wildcards.tool, k, a, '-', '-', '-', '-', '-', '-'), file=open(output.nice_format, 'w'))    

        else:
            misassm, N50, NA50, tot_length, ESIZE_GENOME, CORR_ESIZE_GENOME = parse_quast(outpath+"report.txt")
            k,a = get_k_and_a_for_assembler(input.param)
            print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}".format(wildcards.dataset, wildcards.tool, k, a, ESIZE_GENOME, CORR_ESIZE_GENOME, tot_length, N50, misassm,  NA50), file=open(output.nice_format, 'w'))    



rule time_and_mem:
    input:  stderr=config["OUTBASE"]+"{dataset}/{tool}/{method}_{type}.stderr" #rules.optimal_k_index.output.stderr, rules.optimal_k_sampling.output.stderr,rules.kmergenie.output.stderr #,
    output: outfile=config["OUTBASE"]+"{dataset}/{tool}/{method}_{type}_time_and_mem.txt"
    params: 
        runtime="15:00",
        memsize = "'mem128GB|mem256GB|mem512GB'",
        partition = "core",
        n = "1",
        jobname="{dataset}"+"_time_and_mem",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"]
    run:
        usertime, wallclocktime, memory_gb =  parse_gnu_time(input.stderr)
        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(wildcards.dataset, wildcards.tool, wildcards.method, usertime, wallclocktime, memory_gb), file=open(output.outfile, 'w') )
        

rule performace_latex_table:
    input: files= myfunc  #rules.kmergenie.output.stderr, rules.optimal_k_index.output.stderr, rules.optimal_k_sampling.output.stderr
    output: table=config["OUTBASE"]+"performance_table.tex"
    params: 
        runtime="15:00",
        memsize = "'mem128GB|mem256GB|mem512GB'",
        partition = "core",
        n = "1",
        jobname="performace_latex_table",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"]
    run:
        table_file = open(output.table, 'w')
        print("{0} & {1} & {2} & {3} & {4} & {5} \\\ \hline".format('organism', 'tool','method', 'user time', 'wall clock time', 'peak memory'), file=table_file)
        for file_ in input.files:
            line=open(file_,'r').readlines()[0]
            print("{0} & {1} & {2} & {3} & {4} & {5} \\\ \hline".format(*line.strip().split()), file=table_file)

rule quality_latex_table:
    input: evaluation_targets #map(lambda x: x+"{assembler}.csv", expand(config["OUTBASE"]+"{dataset}/{tool}/result_metrics_",  dataset=config["DATASETS"], tool=config["TOOLS"]) ) 
    output: table=config["OUTBASE"]+"quality_table_{assembler}.tex"
    params: 
        runtime="15:00",
        memsize = "'mem128GB|mem256GB|mem512GB'",
        partition = "core",
        n = "1",
        jobname="quality_latex_table",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"]
    run:
        table_file = open(output.table, 'w')
        print("{0} & {1} & {2} & {3} & {4} & {5} & {6} & {7} & {8} \\\ \hline".format('organism', 'tool', 'k', 'a', 'E-size', 'Total contigs size', 'NG50', 'misassmblies', 'NGA50'), file=table_file)
        for file_ in input:
            line=open(file_,'r').readlines()[0]
            print("{0} & {1} & {2} & {3} & {4} & {5} & {6} & {7} & {8} \\\ \hline".format(*line.strip().split()), file=table_file)




# rule clean:
#     input:
#     output:
#     run:

# rule test:
#     input:
#     output:
#     run:
