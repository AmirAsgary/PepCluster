#!/bin/bash
# submit_pepcluster.sh
#
# Refinement has three sub-steps; at the threshold extremes two of them go
# quadratic on 11M peptides and blow past a 24h wall-clock:
#   * low  threshold (0.3) -> a few GIANT clusters -> the O(k^2) medoid explodes
#                           -> fix with --fast-medoid  (needs pepcluster >= 0.1.3)
#   * high threshold (0.8) -> ~millions of clusters  -> the O(clusters) merge
#                           -> fix with --merge-cap 32 (needs pepcluster >= 0.1.3)
# The middle thresholds (0.4-0.7) sit below both cliffs and finish fine as-is,
# so we only add the fast flags for 0.3 and 0.8. Defaults are unchanged, so the
# middle-threshold runs stay byte-for-byte identical to before.
#
# pepcluster is single-threaded, so 2 cores is plenty (1 to run + headroom).

set -euo pipefail

BASE_OUT="/cbscratch/amirasgary2/pepcluster_out"
INPUT="${BASE_OUT}/peptides.fasta"

for t in 0.3 0.4 0.5 0.6 0.7 0.8; do
    dir=$(printf "%.0f" "$(echo "$t * 100" | bc)")
    outdir="${BASE_OUT}/${dir}"
    mkdir -p "$outdir"

    # Only the two extremes get the fast-refinement flags.
    extra=""
    if [[ "$t" == "0.3" ]]; then
        extra="--fast-medoid"          # kills the giant-cluster medoid blow-up
    elif [[ "$t" == "0.8" ]]; then
        extra="--merge-cap 32"         # kills the many-cluster merge blow-up
    fi

    sbatch \
        -p soeding \
        --mem=64G \
        -c 2 \
        -t 01-00:00:00 \
        -J "pepclust_${dir}" \
        -o "${outdir}/slurm_%j.out" \
        -e "${outdir}/slurm_%j.err" \
        --wrap="pepcluster -i ${INPUT} -o ${outdir} --refinement -t ${t} --backend rust ${extra}"
done
