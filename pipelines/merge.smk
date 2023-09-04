configfile: "pipelines/config.json"

# Borrowed from Annie Froster
rule merge_bams:
    input:
        bam1 = "data/bams/{id}.bam",
        bam2 = "../GAMCC_0.1x/data/bams/{id}.bam"
    output:
        merged_bams = "data/merged_bams/{id}.bam",
        merged_bais = "data/merged_bams/{id}.bam.bai"
    threads: 1
    resources:
        mem = '30G'
    shell: """
        samtools merge -o {output.merged_bams} {input.bam1} {input.bam2}
        samtools index {output.merged_bams}
    """
