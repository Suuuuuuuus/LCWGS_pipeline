# LCWGS_pipeline
Sus' lcwgs pipeline to process low coverage Illumina sequencing data. Current analysis includes:
* Breadth of coverage (see `mosdepth`)
* Depth of coverage
* Jellyfish kmer analysis of error rate
* Fastqc
* Samtools duplication rate
* QUILT imputation

To Do:
* Add a global parameter specification file for users (less priority)
* Optimisation: codes, names, structures, separation punctures, etc
* Modify breadth of coverage graph to solve memory issues ("graphs/fig8_prop_genome_at_least_coverage.png" this maybe problematic)

Notes:
* `bcftools concat` command has a `--ligate` option to join separatedly imputed regions together. Either should increase the buffer size or use `regionStart = 1e6 - 1e5, regionEnd = 2e6 + 1e5, buffer = 0` rather than `regionStart = 1e6, regionEnd = 2e6, buffer = 1e5` (and then use `--ligate-warn`)

Inputs:
* A `config.json` file to specify sample names, etc. You should modify this file which is under directory `pipelines`
* A bunch of `fastq` files
* An index `fa` file
* If kmer analysis is to be performed, high coverage jellyfish `jf` files should be provided for classify-kmers to find errors
* If imputation is to be performed, reference panel (as `vcf` files) and allele frequency files (a bash script under `scripts/extract_gnomAD_MAF.sh` is provided) should be provided. This part will be facilitated with (Prof. Robert Davis' QUILT package (https://github.com/rwdavies/QUILT) and a QUILT-wrap code (https://github.com/rwdavies/QUILT-wrap)).

Outputs:
* `alignment.smk`:
    * Aligned and sorted bam files and there indices
* `index_reference.smk`:
    * Index the provided reference file
* `subsample.smk`:
    * By default, these rules subsample the original `fastq` files to 1x and align them to the reference
* `coverage.smk`:
    * Coverage bedgraphs and subsampled coverage bedgraphs
    * At-least-coverage graph
    * By chromosome, average coverage for each window of genomic positions and graphs
    * Average coverage of each sample
    * Average coverage for each chromosome of each sample (and subsampled sample)
    * Percentage of genome (and subsampled genome) uncovered (termed uncoverage rate) and a graph
* `dup_rate.smk`:
    * Samtools duplication rate and a graph
    * Average fragment length for each sample
    * Total bases of overlapping
    * Proportion of fragments larger/smaller than a threshold for each sample (and subsampled sample)
* `fastqc.smk`:
    * Fastqc results
    * Fastqc duplication rates
* `kmer.smk`:
    * Classify-kmers results on subsampled samples
    * Overall and by-fragment-size kmer accuracy rate
    * Graphs that show kmer accuracy rate for each sample accross read length

Run the Pipeline:
* For now, the whole pipeline is separated into different snakemake files that groups a bunch of jobs together. To run a specific file, use, for example, `snakemake -s pipelines/master.smk -c 1 alignment_all` (needless to say, don't forget to dry-run first by `-n`).
* I don't know if there is a single command that you can run through the whole pipeline, but I deliberately leave a result-wrapper under `rule all` (which basically merge all results into a huge dataframe that hopefully eases downstream analysis).

