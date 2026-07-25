#!/bin/bash

set -euo pipefail

REPO="/scratch/s.2337862/CBS3D++_SI_DD"
SUITE_ROOT="/scratch/s.2337862/CBS3D_cases/blanket/partition_sets/mesh_3eaed573_20260724"
CAMPAIGN_ROOT="/scratch/s.2337862/CBS3D_benchmarks/blanket_strong_scaling_10000_20260725"
RUNNER="$REPO/jobs/mpi_strong_scaling_steady_10000.slurm"
ANALYZER="$REPO/tools/analyze_mpi_strong_scaling.py"

mkdir -p "$CAMPAIGN_ROOT/logs" "$CAMPAIGN_ROOT/P0040" "$CAMPAIGN_ROOT/P0080"

[[ -f "$RUNNER" ]] || { echo "ERROR: missing runner: $RUNNER" >&2; exit 1; }
[[ -f "$ANALYZER" ]] || { echo "ERROR: missing analyzer: $ANALYZER" >&2; exit 1; }
[[ -d "$SUITE_ROOT/P0040/cbs_partitions_40" ]] || { echo "ERROR: missing P40 partition set" >&2; exit 1; }
[[ -d "$SUITE_ROOT/P0080/cbs_partitions_80" ]] || { echo "ERROR: missing P80 partition set" >&2; exit 1; }

existing=$(squeue -h -u "$USER" -n cbs_p40_scale,cbs_p80_scale,cbs_scale_plot -o '%A' || true)
if [[ -n "$existing" ]]; then
    echo "ERROR: scaling campaign jobs already exist: $existing" >&2
    exit 1
fi

cat > "$CAMPAIGN_ROOT/campaign.txt" <<META
case=blanket
simulation_mode=steady_state
iterations=10000
partition=compute
rank_levels=40,80
nodes=1,2
ranks_per_node=40
omp_threads=1
repetitions=1_warmup_plus_3_measured
vtu_output=off
residual_log_every=1000
console_log_every=1000
mesh_tetrahedra=2383581
mesh_nodes=407458
mesh_sha256=3eaed57362bc68144d9fc0ff9bc1dd898e112e3df116d2bb18fad125e39e0136
created=$(date --iso-8601=seconds)
META

P40_JOB=$(sbatch --parsable \
    --partition=compute \
    --nodes=1 \
    --ntasks=40 \
    --ntasks-per-node=40 \
    --cpus-per-task=1 \
    --time=01:15:00 \
    --exclusive \
    --job-name=cbs_p40_scale \
    --output="$CAMPAIGN_ROOT/logs/p40-%j.out" \
    --error="$CAMPAIGN_ROOT/logs/p40-%j.err" \
    --export=ALL,CBS3D_RANKS=40,CBS3D_NODES=1,CBS3D_PARTITION_ROOT="$SUITE_ROOT/P0040/cbs_partitions_40",CBS3D_RUN_ROOT="$CAMPAIGN_ROOT/P0040",CBS3D_REPEATS=4 \
    "$RUNNER")
P40_JOB=${P40_JOB%%;*}

P80_JOB=$(sbatch --parsable \
    --dependency="afterok:$P40_JOB" \
    --partition=compute \
    --nodes=2 \
    --ntasks=80 \
    --ntasks-per-node=40 \
    --cpus-per-task=1 \
    --time=01:00:00 \
    --exclusive \
    --job-name=cbs_p80_scale \
    --output="$CAMPAIGN_ROOT/logs/p80-%j.out" \
    --error="$CAMPAIGN_ROOT/logs/p80-%j.err" \
    --export=ALL,CBS3D_RANKS=80,CBS3D_NODES=2,CBS3D_PARTITION_ROOT="$SUITE_ROOT/P0080/cbs_partitions_80",CBS3D_RUN_ROOT="$CAMPAIGN_ROOT/P0080",CBS3D_REPEATS=4 \
    "$RUNNER")
P80_JOB=${P80_JOB%%;*}

ANALYSIS_JOB=$(sbatch --parsable \
    --dependency="afterok:$P80_JOB" \
    --partition=development \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --time=00:10:00 \
    --job-name=cbs_scale_plot \
    --output="$CAMPAIGN_ROOT/logs/analysis-%j.out" \
    --error="$CAMPAIGN_ROOT/logs/analysis-%j.err" \
    --wrap="python3 '$ANALYZER' --campaign-root '$CAMPAIGN_ROOT'")
ANALYSIS_JOB=${ANALYSIS_JOB%%;*}

cat > "$CAMPAIGN_ROOT/job_chain.txt" <<JOBS
p40_job=$P40_JOB
p80_job=$P80_JOB
analysis_job=$ANALYSIS_JOB
p80_dependency=afterok:$P40_JOB
analysis_dependency=afterok:$P80_JOB
JOBS

echo "P40 job      : $P40_JOB"
echo "P80 job      : $P80_JOB (after P40 passes)"
echo "Analysis job : $ANALYSIS_JOB (after P80 passes)"
echo "Campaign     : $CAMPAIGN_ROOT"
echo
squeue -j "$P40_JOB,$P80_JOB,$ANALYSIS_JOB" \
    -o '%.18i %.20j %.10P %.2t %.19S %.6D %R'
