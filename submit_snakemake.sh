if [ "$#" -eq 0 ]; then
    snakemake --forceall --rulegraph | dot -Tpdf > dag.pdf
elif [ "$#" -eq 3 ]; then
    master="$1"
    cores="$2"
    rule="$3"
    snakemake -s pipelines/master_"$master".smk -c "$cores" --profile slurm/ "$rule"
elif [ "$#" -eq 4 ]; then
    master="$1"
    cores="$2"
    rule="$3"
    option="$4"
    snakemake -s pipelines/master_"$master".smk -c "$cores" "$option" --profile slurm/ --scheduler greedy "$rule"
else
    echo "Invalid inputs."
    exit 1  # Exit with an error status
fi
