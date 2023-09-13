configfile: "pipelines/config.json"

rule rmdup:
    input:
        bam = "data/bams/{id}.bam"
    output:
        dedup_bam = "data/dedup_bams/{id}.bam"
    shell: """
        samtools rmdup {input.bam} {output.dedup_bam}
    """
