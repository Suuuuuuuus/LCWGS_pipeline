# $1 = number of cores
# $2 = rules to run

snakemake -s pipelines/master.smk -c "$1" --profile slurm/ "$2"

#snakemake --forceall --rulegraph | dot -Tpdf > dag.pdf
