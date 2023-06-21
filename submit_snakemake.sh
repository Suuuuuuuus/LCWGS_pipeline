snakemake -s pipelines/master.smk -c "$1" --profile slurm/

#snakemake --forceall --rulegraph | dot -Tpdf > dag.pdf
