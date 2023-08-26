configfile: "pipelines/config.json"

rule alignment:
    input:
        fastq1 = "data/fastq/{id}_1.fastq.gz",
        fastq2 = "data/fastq/{id}_2.fastq.gz", 
        reference = config["ref38"]
    output: 
        bam = temp("data/bams/tmp/{id}.bam")
    resources: 
        mem_mb = 50000
    threads: 4
    shell: """
	bwa mem -t {threads} {input.reference} {input.fastq1} {input.fastq2} | samtools view -b -o {output.bam}
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

