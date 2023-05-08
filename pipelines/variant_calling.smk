configfile: "pipelines/config.json"

ids_1x = [samples_1x['id'] for samples_1x in config['samples_1x']]
ids_20x = [samples_20x['id'] for samples_20x in config['samples_20x']]
chromosome = [i for i in range(1,23)]
id_all = ids_1x + ids_20x

rule imputation_accuracy_NA12878:
        input:
                script = "scripts/calculate_imputation_accuracy_NA12878_v2.py"
        output:
                graph = "graphs/NA12878_imputation_accuracy.png"
        resources:
                mem_mb = 300000
        threads:
                16
        shell: """
                python {input.script}
        """

rule imputation_accuracy_lcwgs:
        input:
                script = "scripts/calculate_imputation_accuracy_lcwgs_v2.py"
        output:
                graph = "graphs/lcwgs_imputation_accuracy.png"
        resources:
                mem_mb = 300000
        threads:
                16
        shell: """
                python {input.script}
        """

'''
rule extract_gnomAD_MAF:
	input:	
		gnomAD_vcf = "data/gnomAD_vcf/gnomad.genomes.v3.1.2.sites.chr{chr}.vcf.bgz"
	output:
		gnomAD_MAF = "results/variant_calling/gnomAD_MAFs/gnomAD_MAF_chr{chr}.txt"
	threads: 4
	shell: """
		/bin/bash -c "zgrep -v '#' {input.gnomAD_vcf} | \
		awk -v "OFS=\t" '{$8=$8;sub(/^.*AF=/, "", $8); sub(/;.*/, "", $8); print $1,$2,$4,$5,$8}' | \
		awk '$5 !~ /^[[:alpha:]]/' | \
		awk '!(length($3)>1 || length($4)>1)' | \
		awk '!($5==0.00000)' > {output.gnomAD_MAF}"
	"""
'''
