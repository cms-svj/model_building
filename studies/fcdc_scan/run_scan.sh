#!/bin/bash
# Generate the Lambda / probVector scan samples with run_model (needs a machine where the repo environment, init.sh, works):
#   8 mpi/Lambda points at probVector 0.333, plus probVector 0.1, 0.2, 0.5, 0.75 at mpi/Lambda = 0.6 and 1.7,
#   each for the complete (fcdc) and simplified (fcdcSimp) model: 32 samples.
# Settings (environment variables):
#   OUT     output directory for the samples (default models/fcdc)
#   STEPS   run_model steps (default "all"; e.g. "cards delphes hist" to rerun Delphes on existing Pythia output, "hist" for histograms only)
#   EVENTS  events per sample (default 10000)
#   NPROC   samples processed in parallel (default 12)
#   ONLY    regex to select samples by name, e.g. "r1p7" or "pv0.75" (default: all)
#   DRY     if set, only print the commands
set -e
cd "$(dirname "$0")/../.."
source init.sh >/dev/null 2>&1

OUT=${OUT:-models/fcdc}
STEPS=${STEPS:-all}
EVENTS=${EVENTS:-10000}
NPROC=${NPROC:-12}
CONFIG=configs/model_fcdc_lambdascan.py
mkdir -p $OUT/logs
CMDS=$OUT/logs/commands_$(date +%Y%m%d_%H%M%S).txt

# sample list: config object (from $CONFIG) and an optional probVector override
{
for r in r0p3 r0p6 r1 r1p4 r1p6 r1p7 r2 r2p5; do for m in full simp; do echo "${r}_${m}"; done; done
for pv in 0.1 0.2 0.5 0.75; do for r in r0p6 r1p7; do for m in full simp; do echo "${r}_${m} $pv"; done; done; done
} | while read obj pv; do
    name=${obj}${pv:+_pv$pv}
    echo "./run_model helper -C $CONFIG -O config.$obj ${pv:+--pvector $pv} --dir $OUT --steps $STEPS --events $EVENTS > $OUT/logs/$name.log 2>&1; echo $name exit=\$?"
done | grep -E -- "${ONLY:-.}" > $CMDS || true

echo "$(wc -l < $CMDS) samples, commands in $CMDS"
if [ -n "$DRY" ]; then
    cat $CMDS
    exit 0
fi
RESULTS=${CMDS%.txt}_results.txt
time xargs -P $NPROC -I{} bash -c '{}' < $CMDS | tee $RESULTS
# run_model exits non-zero if a step fails; list those samples (details in $OUT/logs/<sample>.log)
FAILED=$(grep -v "exit=0$" $RESULTS || true)
if [ -n "$FAILED" ]; then
    echo "failed samples:"; echo "$FAILED"
    exit 1
fi
echo "all samples done"
