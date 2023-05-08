configfile: "pipelines/config.json"

rule index_reference:
    input:
        reference = config["ref38"]
    output:
        amb = "data/reference/GRCh38.fa.amb",
        ann = "data/reference/GRCh38.fa.ann",
        bwt = "data/reference/GRCh38.fa.bwt",
        pac = "data/reference/GRCh38.fa.pac",
        sa = "data/reference/GRCh38.fa.sa"
    resources: mem_mb = 5000
    shell: """
        bwa index {input.reference}
    """
