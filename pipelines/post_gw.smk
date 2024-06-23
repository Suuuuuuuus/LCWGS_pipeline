include: "auxiliary.smk"
include: "software.smk"
configfile: "pipelines/config.json"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus
from lcwgsus.variables import *

rule lift_over_server_vcf:
    input:
        vcf = "results/two-stage-imputation/vanilla/oneKG_hrc/vcf/chr{chr}.dose.vcf.gz",
        reference = "data/references/GRCh38_with_alt.fa",
        chain = "data/ref_panel/b37ToHg38.over.chain",
        dictionary = "data/references/GRCh38_with_alt.dict"
    output:
        tmp1_vcf = temp("results/two-stage-imputation/vanilla/oneKG_hrc/lifted_vcf/chr{chr}.tmp1.vcf.gz"),
        lifted = "results/two-stage-imputation/vanilla/oneKG_hrc/lifted_vcf/chr{chr}.dose.vcf.gz",
        rejected = "results/two-stage-imputation/vanilla/oneKG_hrc/lifted_vcf/chr{chr}.rejected.vcf.gz"
    resources: mem = '80G'
    threads: 4
    params:
        picard = tools["picard_pplus"]
    shell: """
        mkdir -p results/two-stage-imputation/vanilla/oneKG_hrc/lifted_vcf/

        tabix -f {input.vcf}

        {params.picard} LiftoverVcf \
        -I {input.vcf} \
        -O {output.tmp1_vcf} \
        -CHAIN {input.chain} \
        -REJECT {output.rejected} \
        -WMC true \
        --MAX_RECORDS_IN_RAM 50000 \
        -R {input.reference}

        bcftools view -r chr{wildcards.chr} {output.tmp1_vcf} | \
        bcftools sort -Oz -o {output.lifted}
        tabix -f {output.lifted}
    """