# LCWGS_pipeline
Sus' lcwgs pipeline to process low coverage Illumina sequencing data. Current analysis includes:
* Breadth of coverage (see `mosdepth`)
* Depth of coverage
* Jellyfish kmer analysis of error rate
* Fastqc
* Samtools duplication rate

Future analysis to be added:
* Imputation accuracy examination (see `QUILT`)
* (Potentially) Lab prep analysis to examine data quality

To Do:
* Add QUILT imputation and corresponding imputation analysis
* Modify interaction between smk files and python scripts so that they can better accommodate different file names, etc.
* Add a global parameter specification file for users
* Modify breadth of coverage graph to solve memory issues
* Integrate per_bin_kmer_error_rate calculation bash file to the pipeline

Notes:
* `bcftools concat` command has a `--ligate` option to join separatedly imputed regions together. Either should increase the buffer size or use `regionStart = 1e6 - 1e5, regionEnd = 2e6 + 1e5, buffer = 0` rather than `regionStart = 1e6, regionEnd = 2e6, buffer = 1e5` (and then use `--ligate-warn`)
