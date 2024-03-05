LCWGS_pipeline
=======================================

LCWGS_pipeline is a low-coverage whole-genome sequencing (lcWGS) pipeline to process low, high coverage paired Illumina short-read sequencing data and DNA microarray data. 

Current analysis
-------------
* Low-coverage (lc):
    * Preprocessing: adapter trimming using `trimmomatic` and duplicates removal using `fastuniq`
    * Fastqc
    * Subsampling using `seqtk`
    * Breadth of coverage (see `mosdepth`)
    * Depth of coverage
    * Jellyfish k-mer analysis of error rate
    * Duplication rate
    * `QUILT` imputation with user-supplied reference panel
* High-coverage (hc):
    * Preprocessing
    * Chunking up large fastqs for separate alignment using `seqkit`, followed by merging
    * Variant calling with `GATK`
* DNA microarray (chip):
    * Convertion to vcf file format and basic qc
    * Preparing for imputation on the servers
* Others:
    * Lots of plotting, filtering and calculation utilities (though yet under-developped)

Inputs
-------------

* A `pipelines/config.json` file to specify sample names, etc. You should modify this file which is under directory `pipelines`.
* A `data/sample_tsvs/samples_lc.tsv`, a `data/sample_tsvs/chip.tsv` and a `data/sample_tsvs/samples_hc.tsv` file that stores sample name information for lc, hc and chip samples, respectively. Note that each bit of analysis can be run separately, and three corresponding snakemake master files are in place for separate data processing.
* A bunch of `fastq` files for lc and hc samples.
* A bunch of chip files for data processing. Our data requires a genotype file, a sample file and an annotation file.
* An index `fa` file, or several `fa` files if concatenate is set `True` to join them up.
* If k-mer analysis is to be performed, high coverage jellyfish `jf` files should be provided for `classify-kmers` to find errors.
* If imputation is to be performed, reference panel (as `vcf` files) and allele frequency files (a bash script under `scripts/extract_gnomAD_MAF.sh` is provided) should be provided. This part will be facilitated with (Prof. Robert Davis' QUILT package (https://github.com/rwdavies/QUILT) and a QUILT-wrap code (https://github.com/rwdavies/QUILT-wrap)).

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
    * `rm_bed_regions`: `True` if regions from the genome are to be removed in coverage analysis. Specific for the plot_sequencing_skew utility (likely to be optimised later)
    * `bed_regions`: Regions to be removed
    * `access_bed`: Accessible region from the 1KG accessibility track. This is more generally considered in coverage analysis
    * `panels`: List of reference panels to be used for imputation
    * `QUILT_HOME`: Address where the `QUILT` software resides
    * `ANALYSIS_DIR`: Address for the analysis folder. Should be `results/imputation/`
    * `RECOMB_POP`: Three-letter 1KG abbreviation for the population to be used
    * `NGEN`: Number of generations
    * `WINDOWSIZE`: Size of imputation chunks 
    * `BUFFER`: Imputation buffer size
    * `BAMLIST`: Address of the `bamlist.txt` file which stores all bam files to be imputed
    * `PANEL_NAME`: Name of the current imputation panel
* Hc:
    * `samples_hc`: All high-coverage (hc) sample names
    * `make_chunk`: `True` if the hc files need to be chunked
    * `fastq_chunk_size`: Chunk size
    * `hc_panel`: Name of the reference panel to be used for GATK HapCaller to call variants. Only sites in these vcfs will be called to ease comparison
    * `bqsr_known_sites`: For GATK HapCaller to call variants
* Chip:
    * `chip_genotypes`: Chip genotype file
    * `chip_samples`: Chip sample file
    * `chip_annotation`: Chip annotation file

Run the pipeline
-------------

* For now, the whole pipeline is separated into different snakemake files that groups a bunch of jobs together. There are three master files for either lc, hc and chip analysis are in place. To run a specific rule in a specific file, use, for example, `snakemake -s pipelines/master_lc.smk -c 1 alignment_all` (needless to say, don't forget to dry-run first by `-n`).
* Alternatively, a `submit_snakemake.sh` submission script is provided specifically for cluster users. You should first modify the files in the `slurm/` folder to enable a correct communication between snakemake and your job management system. After that, you can run this file with 0, 3 or 4 parameters:
    * If no parameter is passed in, it generates the DAG for all rules.
    * If multiple parameters are passed in, the first should be the name of the master file to run (lc, hc or chip), the second to be number of cores required, and the third should be the name of the snakemake rule file to be run (e.g.: alignment, chip_qc, etc.). The fourth is for other options that are passed to snakemake, like a dry-run option `-nr`.
For example, you can submit by `./submit_snakemake.sh lc 8 alignment_all` which will run all rules till the end of everything in `alignment_all`. It will also run preprocessing rules as alignment requires pre-cleaning of the data. However, running multiple snakemake files are disencouraged, as between some rules the user might need to manually input something or run some scripts. It is always suggested to run the pipeline as granular (per snakemake rule file) as possible.

Documentation
-------------

Online documentation is currently UNAVAILABLE, due to its lack of priority :) 

However, feel free to reach out, submit bugs, and/or make suggestions on improvement of this library. All are welcomed!

Dependencies
------------

The pipeline requires a variety of utilities and packages, including but not limited to:
* A Python 3.8+ environment with the following installations [numpy](https://numpy.org/), [pandas](https://pandas.pydata.org/), [matplotlib](https://matplotlib.org/), [seaborn](https://seaborn.pydata.org/), [statsmodels](https://statsmodels.org/), [scipy](https://scipy.org/) and [lcwgsus](https://github.com/Suuuuuuuus/lcwgsus).
* Other packages include [QUILT](https://github.com/rwdavies/QUILT), [QCTOOL_v2](https://www.chg.ox.ac.uk/~gav/qctool_v2/), [Iorek](https://enkre.net/cgi-bin/code/iorek/dir?ci=trunk), [GATK4](https://github.com/broadinstitute/gatk), [bwa], [samtools], [R], [sqlite3], [seqkit], [bedtools], [fastqc], [multiqc], [jellyfish], [fastuniq], [trimmomatic], [seqtk], [bcftools], etc.

Citation
------------

LCWGS_pipeline is not a published work, so citing the GitHub repo suffices.

Development
-----------

See the main site of LCWGS_pipeline: https://github.com/Suuuuuuuus/LCWGS_pipeline.

Bugs shall be submitted to the [issue tracker](https://github.com/Suuuuuuuus/LCWGS_pipeline/issues). Please kindly provide a REPRODUCIBLE example demonstrating the problem, otherwise debugging per se is foolhardy :)
