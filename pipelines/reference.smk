configfile: "pipelines/config.json"

# Borrowed from Annie Froster
rule concatenate_refs:
    input:
        pf = "data/references/Pf3D7_v3.fasta",
        human = "data/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
        phi = "data/references/NC_001422.fna"
    output:
        human_unzip = temp("data/tmp/human.fa"),
        ref = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta"
    threads: 1
    resources:
        mem = '3G'
    shell: """
        pigz -p {threads} -dc {input.human} > {output.human_unzip}
        cat {output.human_unzip} {input.pf} {input.phi} > {output.ref}
    """

rule index_reference:
    input:
        reference = rules.concatenate_refs.output.ref
        # reference = config["ref38"]
    output:
        # amb = "data/reference/GRCh38.fa.amb",
        # ann = "data/reference/GRCh38.fa.ann",
        # bwt = "data/reference/GRCh38.fa.bwt",
        # pac = "data/reference/GRCh38.fa.pac",
        # sa = "data/reference/GRCh38.fa.sa"
        amb = "data/references/GRCh38_no_alt_Pf3D7_v3_phiX.fa.amb",
        ann = "data/references/GRCh38_no_alt_Pf3D7_v3_phiX.fa.ann",
        bwt = "data/references/GRCh38_no_alt_Pf3D7_v3_phiX.fa.bwt",
        pac = "data/references/GRCh38_no_alt_Pf3D7_v3_phiX.fa.pac",
        sa = "data/references/GRCh38_no_alt_Pf3D7_v3_phiX.fa.sa"
    resources: 
        mem = '3G'
    shell: """
        bwa index {input.reference}
    """
