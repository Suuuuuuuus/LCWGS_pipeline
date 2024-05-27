configfile: "pipelines/config.json"
include: "software.smk"
include: "auxiliary.smk"

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

# ref_indir = ["malariaGen_v1_b37"]
ref_outdir = ["malariaGen_v1_b38"]

def get_indir_vcf(wildcards):
    d = wildcards.ref_outdir.replace("38", "37")
    return "data/ref_panel/" + d + "/" + d + ".chr" + wildcards.chr + ".vcf.gz"

rule lift_over_malariaGen:
    input:
        vcf = get_indir_vcf,
        reference = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta",
        chain = "data/ref_panel/hg19ToHg38.over.chain.gz"
    output:
        tmp_vcf = temp("data/ref_panel/{ref_outdir}/{ref_outdir}.chr{chr}.tmp.vcf"),
        lifted = "data/ref_panel/{ref_outdir}/{ref_outdir}.chr{chr}.vcf.gz",
        rejected = "data/ref_panel/{ref_outdir}/{ref_outdir}.chr{chr}.rejected.vcf.gz"
    resources: mem = '50G'
    threads: 4
    params:
        picard = tools["picard_plus"]
    shell: """
        gunzip -c {input.vcf} | sed 's/Type=String,Number=1/Number=1,Type=String/g' > {output.tmp_vcf}
        bgzip {output.tmp_vcf}
        touch {output.tmp_vcf}
        
        {params.picard} LiftoverVcf \
        -I {output.vcf}.gz \
        -O {output.lifted} \
        -CHAIN {input.chain} \
        -REJECT {output.rejected} \
        -R {input.reference}
    """
