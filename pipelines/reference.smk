configfile: "pipelines/config.json"

concatenate = config['concatenate']

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
        reference = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta" if concatenate else config["ref38"]
    output:
        # amb = "data/reference/GRCh38.fa.amb",
        # ann = "data/reference/GRCh38.fa.ann",
        # bwt = "data/reference/GRCh38.fa.bwt",
        # pac = "data/reference/GRCh38.fa.pac",
        # sa = "data/reference/GRCh38.fa.sa"
        amb = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.amb" if concatenate else "data/references/GRCh38.fa.amb",
        ann = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.ann" if concatenate else "data/references/GRCh38.fa.ann",
        bwt = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.bwt" if concatenate else "data/references/GRCh38.fa.bwt",
        pac = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.pac" if concatenate else "data/references/GRCh38.fa.pac",
        sa = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.sa" if concatenate else "data/references/GRCh38.fa.sa"
    resources:
        mem = '30G'
    shell: """
        bwa index {input.reference}
    """
