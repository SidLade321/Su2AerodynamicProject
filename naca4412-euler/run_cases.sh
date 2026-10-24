#!/bin/bash
# NACA 4412 Euler runs.
#   1. Mesh-convergence study at AoA = 5 deg on the coarse, medium and fine meshes
#   2. Lift-curve sweep on the medium mesh (AoA = -4 ... 8 deg)
# Each case runs in its own folder under runs/, so cases can run side by side
# (set JOBS to the number of physical CPU cores). The histories and surface data that
# analyse.py needs are copied into results/.
# Usage:  python make_mesh.py && ./run_cases.sh      (SU2_CFD must be on the PATH)
set -e
JOBS=${JOBS:-1}
mkdir -p runs results

run() {   # run <mesh level> <AoA>
    local tag="$1_a$2" dir="runs/$1_a$2"
    mkdir -p "$dir"
    sed -e "s/^MESH_FILENAME=.*/MESH_FILENAME= ..\/..\/mesh_naca4412_$1.su2/" \
        -e "s/^AOA=.*/AOA= $2/" naca4412.cfg > "$dir/run.cfg"
    (cd "$dir" && SU2_CFD run.cfg > log.txt)
    cp "$dir/history.csv"      "results/history_$tag.csv"
    cp "$dir/surface_flow.csv" "results/surface_$tag.csv"
    echo "$tag: $(tail -1 "$dir/history.csv" | tr -s ' ')"
}
export -f run

# Longest case first so the parallel jobs finish together
printf "%s\n" "fine 5" "medium 5" "coarse 5" \
       "medium -4" "medium -2" "medium 0" "medium 2" "medium 4" "medium 6" "medium 8" |
    xargs -P "$JOBS" -L 1 bash -c 'run $0 $1'
echo "All runs finished."
