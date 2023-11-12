configfile: "pipelines/config.json"

rule rmdup:
    input:
        bam = "data/bams/{id}.bam"
    output:
        dedup_bam = "data/dedup_bams/{id}.bam"
    shell: """
        samtools rmdup {input.bam} {output.dedup_bam}
    """

rule rmdup_index:
    input:
        dedup_bam = rules.rmdup.output.dedup_bam
    output:
        dedup_bai = "data/dedup_bams/{id}.bam.bai"
    shell: """
        tabix {input.dedup_bam}
    """
