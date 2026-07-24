START=0     # indice iniziale
END=191     # indice finale

# Numero di core per job
NCORES=1
#CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-72}  # 72 di default se non definito
CPUS_PER_TASK=4
max_parallel=$(( CPUS_PER_TASK / NCORES ))

# Funzione per eseguire la simulazione in una cartella
run_job() {
    DIR=$1
    cd "$DIR"
    echo "==> Starting job in $DIR"
    ./Allrun 
    cd ..
}

dirs=()
for i in $(seq -f "%05g" $START $END); do
    dirs+=(ensemble.$i/)
done

# Lancia i job mancanti
for dir in "${dirs[@]}"; do
    LOGFILE="${dir}output.csv"
    
    # Verifica se run.log esiste e contiene la riga con il tempo
    if [[ -f "$LOGFILE" ]]; then
        if grep -q "Total elapsed real time" "$LOGFILE"; then
            echo "==> Skipping $dir (already completed)"
            continue
        fi
    fi

    # Esegui in parallelo se non completato
    run_job "$dir" &

    # Aspetta se si raggiunge il massimo numero di job attivi
    while [ "$(jobs -r | wc -l)" -ge "$max_parallel" ]; do
        sleep 5
    done
done

wait  # aspetta che tutti finiscano
