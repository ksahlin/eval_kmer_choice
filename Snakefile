# path to track and reference
TRACK   = 'hg19.gtf'
REF     = 'hg19.fa'


# sample names and classes
CLASS1  = '101 102'.split()
CLASS2  = '103 104'.split()
SAMPLES = CLASS1 + CLASS2


# path to bam files
CLASS1_BAM = expand('mapped/{sample}.bam', sample=CLASS1)
CLASS2_BAM = expand('mapped/{sample}.bam', sample=CLASS2)


rule all:
    input:
        'out/memory_table.tex',
        'out/k_choice.tex'
        'out/eval_table.tex'        


rule optimal_k:
    input: "/Users/ksahlin/_tmp/testdata_optimal_k/{dataset}.fa"
    output: "/Users/ksahlin/_tmp/kmergenie_{dataset}.dat", "/tmp/kmergenie_{dataset}.stderr"
        csv="{opt_k_folder}/{tool}_{dataset}_{stepsize}.dat"
    output:
        'assembly/{sample}/transcripts.gtf',
        dir='assembly/{sample}'
    threads: 4
    shell:
        'cufflinks --num-threads {threads} -o {output.dir} '
        '--frag-bias-correct {REF} {input}'

rule kmergenie:
    input:
        'mapped/{sample}.bam'
    output:
        'assembly/{sample}/transcripts.gtf',
        dir='assembly/{sample}'
    threads: 4
    shell:
        'cufflinks --num-threads {threads} -o {output.dir} '
        '--frag-bias-correct {REF} {input}'

rule compose_merge:
    input:
        expand('assembly/{sample}/transcripts.gtf', sample=SAMPLES)
    output:
        txt='assembly/assemblies.txt'
    run:
        with open(output.txt, 'w') as out:
            print(*input, sep="\n", file=out)


rule merge_assemblies:
    input:
        'assembly/assemblies.txt'
    output:
        'assembly/merged/merged.gtf', dir='assembly/merged'
    shell:
        'cuffmerge -o {output.dir} -s {REF} {input}'


rule compare_assemblies:
    input:
        'assembly/merged/merged.gtf'
    output:
        'assembly/comparison/all.stats',
        dir='assembly/comparison'
    shell:
        'cuffcompare -o {output.dir}all -s {REF} -r {TRACK} {input}'


rule diffexp:
    input:
        class1=CLASS1_BAM,
        class2=CLASS2_BAM,
        gtf='assembly/merged/merged.gtf'
    output:
        'diffexp/gene_exp.diff', 'diffexp/isoform_exp.diff'
    params:
        class1=",".join(CLASS1_BAM),
        class2=",".join(CLASS2_BAM)
    threads: 8
    shell:
        'cuffdiff --num-threads {threads} {gtf} {params.class1} {params.class2}'