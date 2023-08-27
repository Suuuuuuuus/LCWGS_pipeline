configfile: "pipelines/config.json"

clean_fastq = config['clean_fastq']
reheader = config['reheader']
concatenate = config['concatenate']

rule alignment:
    input:
        fastq1 = "data/fastq_cleaned/{id}_1.fastq.gz" if clean_fastq else "data/fastq/{id}_1.fastq.gz",
        fastq2 = "data/fastq_cleaned/{id}_2.fastq.gz" if clean_fastq else "data/fastq/{id}_2.fastq.gz",
        reference = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta" if concatenate else config["ref38"]
        # fastq1 = "data/fastq/{id}_1.fastq.gz",
        # fastq2 = "data/fastq/{id}_2.fastq.gz", 
        # reference = config["ref38"]
    output: 
        bam = temp("data/bams/tmp/{id}.bam")
    resources: 
        mem_mb = 50000
    threads: 4
    params: 
        reheader = reheader
    shell: """
        bwa mem -t {threads} {input.reference} {input.fastq1} {input.fastq2} | samtools view -b -o {output.bam}
	    if [[ -d "data/bam_headers" && {params.reheader} == "True" ]]
        then
            mkdir -p "data/tmp"
            grep -v PG "data/bam_headers/{wildcards.id}.header.txt" > "data/tmp/{wildcards.id}.reheader.txt"
            samtools reheader "data/tmp/{wildcards.id}.reheader.txt" {output.bam} > "data/tmp/{wildcards.id}.bam"
            rm {output.bam}
            mv "data/tmp/{wildcards.id}.bam" {output.bam}
        fi
    """

rule fixmate:
	input:
		bam = rules.alignment.output.bam
	output:
		fixmate = temp("data/bams/tmp/{id}.fixmate.bam")
	resources: mem_mb = 50000
	shell: """
		samtools fixmate -m {input.bam} {output.fixmate}
	"""

rule sort:
	input:	
		fixmate = rules.fixmate.output.fixmate
	output:
		sorted = temp("data/bams/tmp/{id}.sorted.bam")
	resources: mem_mb = 50000
	shell: """
		samtools sort -o {output.sorted} {input.fixmate}
	"""

rule index:
	input:
		sorted = rules.sort.output.sorted
	output:
		bai = temp("data/bams/tmp/{id}.sorted.bam.bai")
	resources: mem_mb = 50000
	shell: """
		samtools index {input.sorted}
	"""

rule rename_alignment_files:
	input:
		bam = rules.index.input.sorted,
		bai = rules.index.output.bai
	output:
		bam = "data/bams/{id}.bam",
		bai = "data/bams/{id}.bam.bai"
	shell: """
		mv {input.bam} {output.bam}
		mv {input.bai} {output.bai}
	"""

