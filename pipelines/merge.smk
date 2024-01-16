configfile: "pipelines/config.json"

samples_hc = list(pd.read_table(config['samples_hc'], header = None, names = ['Code'])['Code'].values)
samples_hc_split = {}
for i in samples_hc:
    path = "data/file_lsts/hc_fastq_split/" + i + "_split.tsv"
    if os.path.exists(path):
        samples_hc_split[i] = list(pd.read_table(path, header = None, names = ['Code'])['Code'].values)

def sample_hc_to_sample_hc_lst(wildcards):
    return expand("data/bams/{id_ary}.bam", id_ary = samples_hc_split[wildcards.id])

# Merging bams
rule merge_bam:
    input:
        bams = sample_hc_to_sample_hc_lst
    output:
        bam = "data/merge_bams/{id}.bam"
    threads: 8
    resources:
        mem = '50G'
    shell: """
        samtools merge -O BAM {output.bam} {input.bams}
    """
rule index_merge_bam:
    input:
        bam = rules.merge_bam.output.bam
    output:
        bai = "data/merge_bams/{id}.bam.bai"
    threads: 8
    resources:
        mem = '50G'
    shell: """
        samtools index {input.bam}
    """