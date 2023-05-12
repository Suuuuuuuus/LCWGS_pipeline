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
