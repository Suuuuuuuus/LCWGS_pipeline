# $1 = number of cores
# $2 = rules to run
cores=$1
rule="$2"

snakemake -s pipelines/master.smk -c "$cores" --profile slurm/ "$rule"

#snakemake --forceall --rulegraph | dot -Tpdf > dag.pdf
