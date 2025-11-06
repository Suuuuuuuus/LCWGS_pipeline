LCWGS_pipeline
=======================================

LCWGS_pipeline is a low-coverage whole-genome sequencing (lcWGS) pipeline, mainly developed for processing lcWGS paired-end data for genome-wide imputation and HLA imputation. In this repository, pipelines for processing high-coverage data and microarray data have also been included.

Overview
-------------
* Low-coverage (lc):
    * Preprocessing: adapter trimming using `trimmomatic` and duplicates removal using `fastuniq`
    * Fastqc
    * Subsampling using `seqtk`
    * Coverage analysis
    * Sample sequencing QC to assess lcWGS data quality
    * `QUILT` imputation
    * HLA imputation
* High-coverage (hc):
    * Preprocessing
    * Chunking up large fastqs for separate alignment using `seqkit`, followed by merging
    * Variant calling with `GATK`
* DNA microarray (chip):
    * Convertion to vcf file format and basic qc
    * Preparing for imputation on the servers

Since this pipeline is mainly developed to process lcWGS data, description below focuses on such.

Inputs
-------------

* A `pipelines/config.json` file to specify sample names, etc. You should modify this file which is under directory `pipelines`.
* A `data/sample_tsvs/samples_lc.tsv` file that stores sample name information for lc samples. Note that each bit of analysis can be run separately, and three corresponding snakemake master files are in place for separate data processing.
* Sequencing `fastq` files, must be placed and named `data/fastq/{id}_[1|2].fastq.gz` with `{id}` matching the sample identifiers in the specified sample files.
* An index `fa` file.
* If k-mer analysis is to be performed, high coverage jellyfish `jf` files should be provided for `classify-kmers` to find errors.
* If imputation is to be performed, reference panel (as `vcf` files) and allele frequency files (a bash script under `scripts/extract_gnomAD_MAF.sh` is provided) should be provided. This part will be facilitated with (Prof. Robert Davis' QUILT package (https://github.com/rwdavies/QUILT)).


Before running the pipeline
-------------

This section is an example of getting the pipeline ready in use.

* Git clone this repo, optionally rename it to whatever you like, and cd into this directory.

```
git clone https://github.com/Suuuuuuuus/LCWGS_pipeline.git
cd LCWGS_pipeline
```

* Use a conda environment for software management. An example of environment fyi is in `files/conda.yml`. Additional software that is not installable via conda can be placed at `software/` (parallel to `LCWGS_pipeline`), and you can refer to `pipelines/software.smk` as additional examples. 

* Create directories to store data and results.

`mkdir -p results/ data/`

* Transfer fastq files to `data/fastq/`, and create the `samples_lc.tsv` file under `data/sample_tsvs/` which is a text file of no header and single column storing sample IDs per row. If your samples come in batches, you can process each batch of samples separately and merge them post-imputation. This feature is supported in the pipeline, and to achieve this you can simply create a tsv file for each batch and store them in `data/sample_tsvs/batches/samples_<batch>.tsv`. In this case, your `samples_lc.tsv` should contain all samples in this project (which will also be used in the pipeline!), and you can create this file with something like:

```
mkdir -p  data/sample_tsvs/batches/
cat data/sample_tsvs/batches/* > data/sample_tsvs/samples_lc.tsv
```

NB: It is your responsibility to make sure `samples_lc.tsv` has unique sample IDs (in the case of previous sequencing per batch contains duplicate samples).

* You need a reference `.fa` file where the `fastq` files are aligned to. This file should be transfered to `data/references/` and indexed with `bwa index`. You are welcome to do this by hand or modify `pipelines/reference.smk` to meet your own need.

* You need a file that records adapters used in sequencing under `data/miscellaneous/`. An example of IDT adapters is provided in `files/`:
```
mkdir -p data/miscellaneous/
cp files/IDT-Stubby-PE.fa data/miscellaneous/
```

* You need some `.bed` files under `data/bedgraph/` to perform certain QC. These files are provided in `files/`:

```
mkdir -p data/bedgraph/
cp files/*.bed data/bedgraph/
```

Run the pipeline
-------------

* Before running the pipeline, make sure this variable is made global in your environment as well as correctly specified in `pipelines/config.json`, from which the pipeline will learn the home directory.

`export HOME_DIR="<path_to_dir>"`

* For now, the whole pipeline is separated into different snakemake files that groups a bunch of jobs together. There are three master files for either lc, hc and chip analysis are in place. To run a specific rule in a specific file, use, for example, `snakemake -s pipelines/master_lc.smk -c 1 alignment_all` (needless to say, don't forget to dry-run first by `-n`).
* Alternatively, a `submit_snakemake.sh` submission script is provided specifically for cluster users. You should first modify the files in the `slurm/` folder to enable a correct communication between snakemake and your job management system. After that, you can run this file with 0, 3 or 4 parameters:
    * If no parameter is passed in, it generates the DAG for all rules.
    * If multiple parameters are passed in, the first should be the name of the master file to run (lc, hc or chip), the second to be number of cores required, and the third should be the name of the snakemake rule file to be run (e.g.: alignment, chip_qc, etc.). The fourth is for other options that are passed to snakemake, like a dry-run option `-nr`. For example, this will run all rules till the end of everything in `alignment_all`:

    `./submit_snakemake.sh lc 1000 alignment_all -nr`

Explanation of entries in the config file
-------------

* General:
    * `ref38`: Ready-in-use reference file
    * `concatenate`: `True` if several reference files are to be concatenated
    * `sample_linker`: A linker file that tells correspondance of samples if they have multiple names
    * `clean_fastq`: `True` if preprocessing of fastq files is required - drop duplicates, trim adapters, etc.
    * `adapter`: File for performing adapter trimming
    * `reheader`: `True` if user has self-defined headers to replace after alignment, these header files should be put in `data/bam_headers/`
* Lc:
    * `samples_lc`: All low-coverage (lc) sample names
    * `subsample_depth_1x`: `=10666667`, number of reads for a read of length 150 to be 1x
    * `subsample_depth`: Depth to which subsample is performed
    * `num_coverage`: Length of cumsum coverage array (skew plot)
    * `rmdup`: `True` if duplicates are to be removed before lc imputation
    * `rename_samples`: Used to be in rule concat in `imputation.smk`, but currently deprecated
    * `rename_samples_file`: Used to be in rule concat in `imputation.smk`, but currently deprecated
    * `access_bed`: Accessible region from the 1KG accessibility track. This is more generally considered in coverage analysis
    * `hla_bed`: Genomic coordinates of HLA-A, B, C, DQB1, and DRB1
    * `panels`: List of reference panels to be used for imputation
    * `ANALYSIS_DIR`: Address for the analysis folder. Should be `results/imputation/`
    * `RECOMB_POP`: Three-letter 1KG abbreviation for the population to be used
    * `NGEN`: Number of generations
    * `WINDOWSIZE`: Size of imputation chunks 
    * `BUFFER`: Imputation buffer size
    * `bam_batch`: Number of QUILT job to submit for parallelisation purposes

Lots of the other entries in this file are currently for private use and can be neglected.

Dependencies
------------

The pipeline requires a variety of utilities and packages, including but not limited to:
* A Python 3.8+ environment with the following installations [numpy](https://numpy.org/), [pandas](https://pandas.pydata.org/), [matplotlib](https://matplotlib.org/), [seaborn](https://seaborn.pydata.org/), [statsmodels](https://statsmodels.org/), [scipy](https://scipy.org/) and [lcwgsus](https://github.com/Suuuuuuuus/lcwgsus).
* Other packages include [QUILT](https://github.com/rwdavies/QUILT), [QCTOOL_v2](https://www.chg.ox.ac.uk/~gav/qctool_v2/), [Iorek](https://enkre.net/cgi-bin/code/iorek/dir?ci=trunk), [GATK4](https://github.com/broadinstitute/gatk), [bwa], [samtools], [R], [sqlite3], [seqkit], [bedtools], [fastqc], [multiqc], [jellyfish], [fastuniq], [trimmomatic], [seqtk], [bcftools], etc.

Citation
------------

LCWGS_pipeline is not yet a published work, so citing the GitHub repo suffices.

Development
-----------

See the main site of LCWGS_pipeline: https://github.com/Suuuuuuuus/LCWGS_pipeline.

Bugs shall be submitted to the [issue tracker](https://github.com/Suuuuuuuus/LCWGS_pipeline/issues).
