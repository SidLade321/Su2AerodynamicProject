#!/bin/bash
# Runs the NACA 0012 Euler case at M = 0.3 for several angles of attack.
# Each run's convergence history is saved as history_aoa<angle>.csv
set -e
for AOA in 0 2 4 6 8; do
    echo "=== AoA = $AOA deg ==="
    sed "s/^AOA=.*/AOA= $AOA/" liftcurve_base.cfg > run.cfg
    SU2_CFD run.cfg > log_aoa$AOA.txt
    cp history.csv history_aoa$AOA.csv
    tail -1 history_aoa$AOA.csv
done
echo "Sweep finished."
