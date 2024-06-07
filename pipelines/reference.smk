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
ref_outdirs = ["malariaGen_v1_b38"]

rule create_ref_dict:
    input:
        reference = "data/references/GRCh38_with_alt.fa"
    output:
        dictionary = "data/references/GRCh38_with_alt.dict"
    resources: mem = '50G'
    threads: 4
    params:
        picard = tools["picard_plus"]
    shell: """
        {params.picard} CreateSequenceDictionary \
        -R {input.reference} \
        -O {output.dictionary}
    """ 

def get_indir_vcf(wildcards):
    d = wildcards.ref_outdir.replace("38", "37")
    c = wildcards.chr
    return "data/ref_panel/" + d + "/" + d + ".chr" + c + ".vcf.gz"

rule lift_over_malariaGen_v1:
    input:
        vcf = get_indir_vcf,
        reference = "data/references/GRCh38_with_alt.fa",
        chain = "data/ref_panel/b37ToHg38.over.chain",
        dictionary = "data/references/GRCh38_with_alt.dict"
    output:
        tmp1_vcf = temp("data/ref_panel/{ref_outdir}/{ref_outdir}.chr{chr}.tmp1.vcf.gz"),
        tmp_vcf = temp("data/ref_panel/{ref_outdir}/{ref_outdir}.chr{chr}.tmp.vcf.gz"),
        tmp2_vcf = temp("data/ref_panel/{ref_outdir}/{ref_outdir}.chr{chr}.tmp2.vcf.gz"),
        lifted = "data/ref_panel/{ref_outdir}/{ref_outdir}.chr{chr}.vcf.gz",
        rejected = "data/ref_panel/{ref_outdir}/{ref_outdir}.chr{chr}.rejected.vcf.gz"
    resources: mem = '90G'
    threads: 4
    params:
        picard = tools["picard_ppplus"],
        c = "{chr}",
        rename_chr = "data/ref_panel/{ref_outdir}/rename_chr{chr}.tsv"
    shell: """
        gunzip -c {input.vcf} | sed 's/Type=String,Number=1/Number=1,Type=String/g' | bcftools view -Oz -o {output.tmp1_vcf}
        tabix -f {output.tmp1_vcf}

        bcftools view -Oz -o {output.tmp_vcf} {output.tmp1_vcf}
        tabix -f {output.tmp_vcf}

        if [ {params.c} -lt 10 ]; then
            echo "0{params.c}    {params.c}" > {params.rename_chr}
            bcftools annotate --rename-chrs {params.rename_chr} -Oz -o {output.tmp2_vcf} {output.tmp_vcf}
            rm {params.rename_chr}
        else
            cp {output.tmp1_vcf}.gz {output.tmp2_vcf}
        fi
        tabix {output.tmp2_vcf}
        
        {params.picard} LiftoverVcf \
        -I {output.tmp2_vcf} \
        -O {output.lifted} \
        -CHAIN {input.chain} \
        -REJECT {output.rejected} \
        -WMC true \
        -R {input.reference}
    """

rule convert_shapeit_to_vcf:
    input:
        shapeit = "data/ref_panel/shapeit/MalariaGEN_combined_reference_panel_v3_chr{chr}_biallelic_snps.shapeit.haps.gz",
        samples = "data/ref_panel/shapeit/malariaGen_v3_b37.samples"
    output:
        vcf = temp("data/ref_panel/malariaGen_v3_b37_alone/malariaGen_v3_b37_alone.chr{chr}.tmp.vcf.gz"),
        tmp1 = temp("data/ref_panel/malariaGen_v3_b37_alone/malariaGen_v3_b37_alone.chr{chr}.tmp1.vcf.gz")
    resources: mem = '50G'
    threads: 4
    params:
        qctool = tools["qctool"]
    shell: """
        mkdir -p data/ref_panel/malariaGen_v3_b37_alone/

        {params.qctool} \
        -filetype shapeit_haplotypes \
        -g {input.shapeit} \
        -s {input.samples} \
        -og {output.tmp1}

        gunzip -c {output.tmp1} | sed 's/Type=1,Number=String/Number=1,Type=String/g' | bcftools view -Oz -o {output.vcf}
        tabix -f {output.vcf}
    """

rule lift_over_malariaGen_v3:
    input:
        vcf = "data/ref_panel/malariaGen_v3_b37_alone/malariaGen_v3_b37_alone.chr{chr}.tmp.vcf.gz",
        reference = "data/references/GRCh38_with_alt.fa",
        chain = "data/ref_panel/b37ToHg38.over.chain",
        dictionary = "data/references/GRCh38_with_alt.dict"
    output:
        tmp_vcf = temp("data/ref_panel/malariaGen_v3_b38_alone/malariaGen_v3_b38_alone.chr{chr}.tmp.vcf.gz"),
        lifted = "data/ref_panel/malariaGen_v3_b38_alone/malariaGen_v3_b38_alone.chr{chr}.vcf.gz",
        rejected = "data/ref_panel/malariaGen_v3_b38_alone/malariaGen_v3_b38_alone.chr{chr}.rejected.vcf.gz"
    resources: mem = '80G'
    threads: 4
    params:
        picard = tools["picard_pplus"]
    shell: """
        mkdir -p data/ref_panel/malariaGen_v3_b38_alone/

        bcftools view -Oz -o {output.tmp_vcf} {input.vcf}

        {params.picard} LiftoverVcf \
        -I {output.tmp_vcf} \
        -O {output.lifted} \
        -CHAIN {input.chain} \
        -REJECT {output.rejected} \
        -WMC true \
        -R {input.reference}
    """

# ";NA;" is a common feature for a match in both ref panels
rule merge_malariaGen_v3_with_oneKG:
    input:
        mg = "data/ref_panel/malariaGen_v3_b38_alone/malariaGen_v3_b38_alone.chr{chr}.vcf.gz",
        oneKG = "data/ref_panel/oneKG/oneKG.chr{chr}.vcf.gz"
    output:
        vcf = "data/ref_panel/malariaGen_v3_b38/malariaGen_v3_b38.chr{chr}.vcf.gz",
        tbi = "data/ref_panel/malariaGen_v3_b38/malariaGen_v3_b38.chr{chr}.vcf.gz.tbi"
    resources: mem = '70G'
    threads: 4
    shell: """
        mkdir -p data/ref_panel/malariaGen_v3_b38/
        
        bcftools merge {input.mg} {input.oneKG} -Ou | \
        bcftools view -i 'ID~";NA;"' -Oz -o {output.vcf}

        tabix {output.vcf} 
    """