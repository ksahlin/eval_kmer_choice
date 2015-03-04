# The following test will take a two read files estimate best k with kmer_genie

#import os
configfile: "config.json"
REF="/Users/ksahlin/Documents/workspace/data/data/real_genomes/staph/data/reference.fasta"

DATASETS="reads1 reads2".split()
#paths to infiles
#INFILES=expand("/Users/ksahlin/_tmp/testdata_optimal_k/{dataset}.fa", dataset=DATASETS)

# a pseudo-rule that collects the target files 
rule all:
   input:  expand("/tmp/kmergenie_{dataset}.dat", dataset=DATASETS)

# a general rule using wildcards that does the work


rule kmergenie:
	input: "/Users/ksahlin/_tmp/testdata_optimal_k/{dataset}.fa"
	output: "/tmp/kmergenie_{dataset}.dat" #k,a,realtime,usertime,memory
	shell:
		"""
			{config[kmergenie_rules][load_env]}
			kmergenie {input} -o {output}
		"""

