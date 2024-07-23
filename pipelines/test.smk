configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"

import json
import pandas as pd
import numpy as np
import sys
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

rule all:
    input:
        lifted = "data/hla_from_Ruth/multiEth_sites.b38.vcf.gz",
        vcf = "data/hla_from_Ruth/chr6.vcf.gz"
        
rule liftover:
    input:
        vcf = "data/hla_from_Ruth/clean_GGVP_hla_GRCh37.vcf.gz",
        reference = "data/references/GRCh38_with_alt.fa",
        chain = "data/ref_panel/b37ToHg38.over.chain",
        dictionary = "data/references/GRCh38_with_alt.dict"
    output:
        tmp_vcf = temp("data/hla_from_Ruth/multiEth_sites.b38.tmp.vcf.gz"),
        lifted = "data/hla_from_Ruth/multiEth_sites.b38.vcf.gz",
        rejected = "data/hla_from_Ruth/multiEth_sites.b38.rejected.vcf.gz"
    resources: mem = '30G'
    threads: 4
    params:
        picard = tools["picard"]
    shell: """
        {params.picard} LiftoverVcf \
        -I {input.vcf} \
        -O {output.tmp_vcf} \
        -CHAIN {input.chain} \
        -REJECT {output.rejected} \
        -WMC true \
        --MAX_RECORDS_IN_RAM 50000 \
        -R {input.reference}

        tabix -f {output.tmp_vcf}

        bcftools view -r chr6 {output.tmp_vcf} | \
        bcftools sort -Oz -o {output.lifted}
        tabix -f {output.lifted}
    """

rule test:
    input:
        lifted = "data/hla_from_Ruth/multiEth_sites.b38.vcf.gz",
        b38_scaffold = "results/hla/reference/multiEth_sites.b38.vcf.gz"
    output:
        vcf = "data/hla_from_Ruth/chr6.vcf.gz",
        scaffold = temp("data/hla_from_Ruth/multiEth_sites.b38.txt")
    shell: """
        bcftools query -f '%CHROM\t%POS\n' {input.b38_scaffold} > {output.scaffold}

        bcftools view -R {output.scaffold} -m2 -M2 -v snps \
        -Oz -o {output.vcf} {input.lifted}
    """