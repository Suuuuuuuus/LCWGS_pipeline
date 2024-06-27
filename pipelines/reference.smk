configfile: "pipelines/config.json"
include: "software.smk"
include: "auxiliary.smk"

concatenate = config['concatenate']
RECOMB_POP = config["RECOMB_POP"]

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

rule lift_over_malariaGen_v1:
    input:
        vcf = "data/ref_panel/malariaGen_v1_b37/malariaGen_v1_b37.chr{chr}.vcf.gz",
        reference = "data/references/GRCh38_with_alt.fa",
        chain = "data/ref_panel/b37ToHg38.over.chain",
        dictionary = "data/references/GRCh38_with_alt.dict"
    output:
        tmp1_vcf = temp("data/ref_panel/malariaGen_v1_b38/malariaGen_v1_b38.chr{chr}.tmp1.vcf.gz"),
        tmp_vcf = temp("data/ref_panel/malariaGen_v1_b38/malariaGen_v1_b38.chr{chr}.tmp.vcf.gz"),
        tmp2_vcf = temp("data/ref_panel/malariaGen_v1_b38/malariaGen_v1_b38.chr{chr}.tmp2.vcf.gz"),
        lifted = "data/ref_panel/malariaGen_v1_b38/malariaGen_v1_b38.chr{chr}.vcf.gz",
        rejected = "data/ref_panel/malariaGen_v1_b38/malariaGen_v1_b38.chr{chr}.rejected.vcf.gz"
    resources: mem = '110G'
    threads: 8
    params:
        picard = tools["picard_ppplus"],
        c = "{chr}",
        rename_chr = "data/ref_panel/malariaGen_v1_b38/rename_chr{chr}.tsv"
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
            cp {output.tmp_vcf} {output.tmp2_vcf}
        fi
        tabix {output.tmp2_vcf}
        
        {params.picard} LiftoverVcf \
        -I {output.tmp2_vcf} \
        -O {output.tmp1_vcf} \
        -CHAIN {input.chain} \
        -REJECT {output.rejected} \
        -WMC true \
        --MAX_RECORDS_IN_RAM 50000 \
        -R {input.reference}

        tabix -f {output.tmp1_vcf}

        bcftools view -r chr{wildcards.chr} {output.tmp1_vcf} | \
        bcftools sort -Oz -o {output.lifted}
        tabix -f {output.lifted}
    """

rule convert_shapeit_to_vcf:
    input:
        shapeit = "data/ref_panel/shapeit/MalariaGEN_combined_reference_panel_v3_chr{chr}_biallelic_snps.shapeit.haps.gz",
        samples = "data/ref_panel/shapeit/malariaGen_v3_b37.samples"
    output:
        vcf = "data/ref_panel/malariaGen_v3_b37_alone/malariaGen_v3_b37_alone.chr{chr}.vcf.gz",
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
        vcf = "data/ref_panel/malariaGen_v3_b37_alone/malariaGen_v3_b37_alone.chr{chr}.vcf.gz",
        reference = "data/references/GRCh38_with_alt.fa",
        chain = "data/ref_panel/b37ToHg38.over.chain",
        dictionary = "data/references/GRCh38_with_alt.dict"
    output:
        tmp1_vcf = temp("data/ref_panel/malariaGen_v3_b38_alone/malariaGen_v3_b38_alone.chr{chr}.tmp1.vcf.gz"),
        tmp2_vcf = temp("data/ref_panel/malariaGen_v3_b38_alone/malariaGen_v3_b38_alone.chr{chr}.tmp2.vcf.gz"),
        lifted = "data/ref_panel/malariaGen_v3_b38_alone/malariaGen_v3_b38_alone.chr{chr}.vcf.gz",
        rejected = "data/ref_panel/malariaGen_v3_b38_alone/malariaGen_v3_b38_alone.chr{chr}.rejected.vcf.gz"
    resources: mem = '80G'
    threads: 4
    params:
        picard = tools["picard_pplus"]
    shell: """
        mkdir -p data/ref_panel/malariaGen_v3_b38_alone/

        bcftools view -Oz -o {output.tmp1_vcf} {input.vcf}

        {params.picard} LiftoverVcf \
        -I {output.tmp1_vcf} \
        -O {output.tmp2_vcf} \
        -CHAIN {input.chain} \
        -REJECT {output.rejected} \
        -WMC true \
        --MAX_RECORDS_IN_RAM 50000 \
        -R {input.reference}

        bcftools view -r chr{wildcards.chr} {output.tmp2_vcf} | \
        bcftools sort -Oz -o {output.lifted}
        tabix -f {output.lifted}
    """

to_merge = ['malariaGen_v3_b38_alone', 'oneKG']

rule prepare_merge_1KGmGenv3_vcf:
    input:
        vcf = "data/ref_panel/{to_merge}/{to_merge}.chr{chr}.vcf.gz"
    output:
        haps = temp("data/ref_panel/malariaGen_v3_b38/tmp/{to_merge}.chr{chr}.hap"),
        legend = temp("data/ref_panel/malariaGen_v3_b38/tmp/{to_merge}.chr{chr}.legend")
    resources: mem = '30G'
    threads: 4
    params: outdir = "data/ref_panel/malariaGen_v3_b38/tmp/"
    shell: """
        mkdir -p {params.outdir}

        bcftools norm -m+ {input.vcf} | \
        bcftools view -m2 -M2 -v snps | \
        bcftools sort | \
        bcftools convert -h {params.outdir}{wildcards.to_merge}.chr{wildcards.chr}

        gunzip {params.outdir}{wildcards.to_merge}.chr{wildcards.chr}.hap.gz
        gunzip {params.outdir}{wildcards.to_merge}.chr{wildcards.chr}.legend.gz
    """

rule prepare_merge_1KGmGenv3_sample:
    input:
        mg = "data/ref_panel/malariaGen_v3_b38_alone/malariaGen_v3_b38_alone.chr22.vcf.gz",
        oneKG = "data/ref_panel/oneKG/oneKG.chr22.vcf.gz"
    output:
        tmp_sample = temp("data/ref_panel/malariaGen_v3_b38/tmp/tmp.sample"),
        sample = temp("data/ref_panel/malariaGen_v3_b38/tmp/merged.samples")
    resources: mem = '30G'
    threads: 1
    params: outdir = "data/ref_panel/malariaGen_v3_b38/tmp/"
    shell: """
        mkdir -p {params.outdir}

        echo "sample population group sex" >> {output.sample}

        bcftools query -l {input.oneKG} >> {output.tmp_sample}
        bcftools query -l {input.mg} >> {output.tmp_sample}

        for i in $(cat {output.tmp_sample}); do
            echo $i $i $i 2 >> {output.sample}
        done
    """

# -k_hap = 2018 (N = 1009) was calculated from this https://www.internationalgenome.org/data-portal/population
rule merge_1KGmGenv3_per_chunk:
    input:
        haps = expand("data/ref_panel/malariaGen_v3_b38/tmp/{to_merge}.chr{chr}.hap", to_merge = to_merge, allow_missing = True),
        legends = expand("data/ref_panel/malariaGen_v3_b38/tmp/{to_merge}.chr{chr}.legend", to_merge = to_merge, allow_missing = True),
        gen_map = f"data/imputation_accessories/maps/{RECOMB_POP}-chr{{chr}}-final.b38.txt",
        sample = rules.prepare_merge_1KGmGenv3_sample.output.sample
    output:
        haps = temp("data/ref_panel/malariaGen_v3_b38/regions/chr{chr}.{regionStart}.{regionEnd}.hap"),
        legend = temp("data/ref_panel/malariaGen_v3_b38/regions/chr{chr}.{regionStart}.{regionEnd}.legend"),
        vcf = "data/ref_panel/malariaGen_v3_b38/regions/chr{chr}.{regionStart}.{regionEnd}.vcf.gz",
    resources: mem = '70G'
    threads: 4
    params: 
        impute2 = tools['impute2'],
        outdir = "data/ref_panel/malariaGen_v3_b38/regions/",
        output_prefix = "data/ref_panel/malariaGen_v3_b38/regions/chr{chr}.{regionStart}.{regionEnd}",
        mGen_haps = "data/ref_panel/malariaGen_v3_b38/tmp/malariaGen_v3_b38_alone.chr{chr}.hap",
        mGen_legend = "data/ref_panel/malariaGen_v3_b38/tmp/malariaGen_v3_b38_alone.chr{chr}.legend",
        oneKG_haps = "data/ref_panel/malariaGen_v3_b38/tmp/oneKG.chr{chr}.hap",
        oneKG_legend = "data/ref_panel/malariaGen_v3_b38/tmp/oneKG.chr{chr}.legend",
    shell: """
        mkdir -p {params.outdir}

       {params.impute2} \
        -merge_ref_panels_output_ref {params.output_prefix}.tmp \
        -m {input.gen_map} \
        -h {params.oneKG_haps} \
           {params.mGen_haps} \
        -l {params.oneKG_legend} \
           {params.mGen_legend} \
        -int {wildcards.regionStart} {wildcards.regionEnd} \
        -k_hap 2018 1510 \
        -Ne 20000

        awk -F ' ' 'NR==1 {{print; next}} {{$1 = "chr{wildcards.chr}:"$2"_"$3"_"$4; print $0}}' \
        {params.output_prefix}.tmp.legend > {output.legend}
        mv {params.output_prefix}.tmp.hap {output.haps}

        bgzip -f {output.legend}
        bgzip -f {output.haps}
        touch {output.legend}
        touch {output.haps}

        cp {input.sample} {params.output_prefix}.samples

        bcftools convert -H {params.output_prefix} | bcftools sort -Oz -o {output.vcf}
        tabix -f {output.vcf}
    """

region_file = "data/imputation_accessories/5Mb_chunks.json"
mGen_vcf_prefix = "data/ref_panel/malariaGen_v3_b38/regions/chr"
mGen_chunk_RData, mGen_chunk_vcf_lst, mGen_chunk_vcf_dict = get_vcf_concat_lst(region_file, '', mGen_vcf_prefix)

def get_input_vcfs_as_list(wildcards):
    return(mGen_chunk_vcf_dict[str(wildcards.chr)])

def get_input_vcfs_as_string(wildcards):
    return(" ".join(map(str, mGen_chunk_vcf_dict[str(wildcards.chr)])))

rule merge_1KGmGenv3_chunks:
    input:
        vcfs = get_input_vcfs_as_list
    output:
        vcf = "data/ref_panel/malariaGen_v3_b38/malariaGen_v3_b38.chr{chr}.vcf.gz",
    resources: mem = '30G'
    threads: 1
    params: 
        input_string = get_input_vcfs_as_string,
    shell: """
        bcftools concat --ligate-force -Oz -o {output.vcf} {params.input_string}
        tabix {output.vcf}
    """