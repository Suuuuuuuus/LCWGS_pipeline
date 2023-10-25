if [ "$#" -eq 1 ]; then
    snakemake --forceall --rulegraph | dot -Tpdf > dag.pdf
elif [ "$#" -eq 2 ]; then
    cores="$1"
    rule="$2"
    snakemake -s pipelines/master.smk -c "$cores" --profile slurm/ "$rule"
elif [ "$#" -eq 3 ]; then
    cores="$1"
    rule="$2"
    option="$3"
    snakemake -s pipelines/master.smk -c "$cores" "$option" --profile slurm/ "$rule"
else
    echo "Invalid inputs."
    exit 1  # Exit with an error status
fi
