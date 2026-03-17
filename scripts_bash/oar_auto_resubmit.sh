#!/bin/bash
# submit_until_done.sh
MAX_RETRIES=10
SCRIPT=$1
NEXTFLOW_DIR=${2:-.}  # répertoire contenant .nextflow.log, par défaut .

for i in $(seq 1 $MAX_RETRIES); do
    echo "=== Attempt $i / $MAX_RETRIES ==="

    # ── 1. Soumettre le job et capturer l'output ──────────────────────────────
    OAR_OUTPUT=$(oarsub -S "$SCRIPT")
    echo "$OAR_OUTPUT"

    # ── 2. Extraire le job ID ─────────────────────────────────────────────────
    JOB_ID=$(echo "$OAR_OUTPUT" | grep "OAR_JOB_ID" | cut -d'=' -f2 | tr -d ' ')

    if [ -z "$JOB_ID" ]; then
        echo "❌ Failed to get job ID, aborting"
        exit 1
    fi
    echo "✅ Job submitted: $JOB_ID"

    # ── 3. Polling : attendre la fin du job ───────────────────────────────────
    echo "⏳ Waiting for job $JOB_ID to finish..."
    while true; do
        sleep 60

        # oarstat retourne rien si le job n'existe plus (= terminé)
        STATUS=$(oarstat -j "$JOB_ID" -s 2>/dev/null)

        if [ -z "$STATUS" ]; then
            echo "Job $JOB_ID no longer in oarstat → finished"
            break
        fi

        echo "  $(date '+%H:%M:%S') Job $JOB_ID status: $STATUS"

        # Sortie anticipée si erreur OAR
        if echo "$STATUS" | grep -qiE "Error|Finishing|Terminated"; then
            echo "Job ended with status: $STATUS"
            break
        fi
    done

    # # ── 4. Vérifier si Nextflow a terminé avec succès ─────────────────────────
    # LOG="$NEXTFLOW_DIR/.nextflow.log"

    # if [ ! -f "$LOG" ]; then
    #     echo "❌ No .nextflow.log found at $LOG"
    #     exit 1
    # fi

    # if grep -q "Workflow complete" "$LOG"; then
    #     echo "🎉 Workflow complete after $i attempt(s)!"
    #     exit 0
    # else
    #     echo "⚠️  Workflow not complete yet, re-submitting..."
    # fi

    # ── 4. Vérifier si le HTML final existe ───────────────────────────────────
    Meta_file="$NEXTFLOW_DIR/08_metaanalysis.html"

    if [ ! -f "$Meta_file" ]; then
        echo "⚠️  $Meta_file not found yet, re-submitting..."
    else
        echo "🎉 Found report: $Meta_file after $i attempt(s)!"
        exit 0
    fi
done


echo "❌ Max retries ($MAX_RETRIES) reached without completion"
exit 1