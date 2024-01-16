configfile: "pipelines/config.json"

# Removing duplicates
rule fastuniq: # Currently deprecated as we are basically gonna remove in markdup
    input:
        fastq1 = "data/fastq/{id}_1.fastq.gz",
        fastq2 = "data/fastq/{id}_2.fastq.gz"
    output:
        fastq1 = temp("data/tmp/{id}_fast_uniq_1.fastq.gz"),
        fastq2 = temp("data/tmp/{id}_fast_uniq_2.fastq.gz"),
        filelist = temp("data/tmp/{id}.list"),
        fastq1_uncompress = temp("data/tmp/input_{id}_fast_uniq_1.fastq.gz"),
        fastq2_uncompress = temp("data/tmp/input_{id}_fast_uniq_2.fastq.gz"),
        fastq1_unzip = temp("data/tmp/{id}_fast_uniq_1.fastq"),
        fastq2_unzip = temp("data/tmp/{id}_fast_uniq_2.fastq")
    threads: 2
    resources:
        mem = '10G'
    shell: """
        pigz -p {threads} -dc {input.fastq1} > {output.fastq1_uncompress}
        pigz -p {threads} -dc {input.fastq2} > {output.fastq2_uncompress}
        echo {output.fastq1_uncompress} > {output.filelist}
        echo {output.fastq2_uncompress} >> {output.filelist}
        fastuniq -i {output.filelist} -o {output.fastq1_unzip} -p {output.fastq2_unzip}
        gzip -c {output.fastq1_unzip} > {output.fastq1}
        gzip -c {output.fastq2_unzip} > {output.fastq2}
    """

# Adapter trimming
rule trimmomatic:
    input:
        fastq1 = "data/fastq/tmp/{id}/{id}_1.fastq.gz", # For chunked fastq files
        fastq2 = "data/fastq/tmp/{id}/{id}_2.fastq.gz"
        # fastq1 = "data/fastq/{id}_1.fastq.gz", # For non-fastuniq fastq files
        # fastq2 = "data/fastq/{id}_2.fastq.gz"
        # fastq1 = rules.fastuniq.output.fastq1, # For fast-uniq'ed fastq files
        # fastq2 = rules.fastuniq.output.fastq2
    output:
        fwd_pair = "data/fastq_cleaned/{id}_1.fastq.gz",
        rev_pair = "data/fastq_cleaned/{id}_2.fastq.gz",
        fwd_unpair = "data/fastq_cleaned/{id}_unpaired_1.fastq.gz",
        rev_unpair = "data/fastq_cleaned/{id}_unpaired_2.fastq.gz"
    params:
        adapters = config['adapter']
    threads: 2
    resources:
        mem = '10G'
    shell: """
        trimmomatic PE \
        {input.fastq1} {input.fastq2} \
        {output.fwd_pair} {output.fwd_unpair} \
        {output.rev_pair} {output.rev_unpair} \
        ILLUMINACLIP:{params.adapters}:2:30:10:2:true MINLEN:50
    """