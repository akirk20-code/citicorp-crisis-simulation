#!/usr/bin/env bash
# ============================================================
# mesh_all_cases.sh
# Generate STLs and mesh for all 6 geometry variants.
# Run from WSL: bash /mnt/c/.../mesh_all_cases.sh
# ============================================================

set -e  # Exit on any error

PYTHON=/mnt/c/Users/kirka/AppData/Local/Programs/Python/Python313/python.exe
WIN_CFD="/mnt/c/Users/kirka/OneDrive/Documents/GMU PhD/OR 750 Reliability, Safety, and Risk/Citi Sim/citicorp_cfd"
WSL_HOME=/home/kirka

# Load OpenFOAM Foundation v13
source /opt/openfoam13/etc/bashrc

# Cases: (script_name  case_dir)
declare -A CASES=(
    [osm]="citicorp_cfd_osm"
    [nyc3d]="citicorp_cfd_nyc3d"
    [lod2]="citicorp_cfd_lod2"
    [osm_1978]="citicorp_cfd_osm_1978"
    [nyc3d_1978]="citicorp_cfd_nyc3d_1978"
    [lod2_1978]="citicorp_cfd_lod2_1978"
)

NPROCS=8  # Meshing cores (conservative for snappyHexMesh memory)

for key in "${!CASES[@]}"; do
    case_dir="$WSL_HOME/${CASES[$key]}"
    stl_src="$WIN_CFD/constant/triSurface"
    stl_dst="$case_dir/constant/triSurface"

    echo ""
    echo "========================================================"
    echo "Case: $key  →  $case_dir"
    echo "========================================================"

    # Step 1: Generate STLs (run Python script from Windows-accessible location)
    echo "[1/5] Generating STLs with generate_stl_${key}.py..."
    cd "$WIN_CFD"
    $PYTHON "generate_stl_${key}.py" 2>&1 | tail -20

    # Step 2: Copy STLs to case
    echo "[2/5] Copying STLs to $stl_dst..."
    mkdir -p "$stl_dst"
    cp "$stl_src/citicorp_tower.stl"  "$stl_dst/"
    cp "$stl_src/citicorp_stilts.stl" "$stl_dst/"
    # Surroundings STL may have different name depending on script
    if ls "$stl_src"/surroundings*.stl 2>/dev/null | grep -q .; then
        cp "$stl_src"/surroundings*.stl "$stl_dst/surroundings.stl"
    fi
    ls -lh "$stl_dst/"

    # Step 3: Extract surface features
    echo "[3/5] Running surfaceFeatureExtract..."
    cd "$case_dir"
    surfaceFeatureExtract > log.surfaceFeatureExtract 2>&1 || \
        echo "  Warning: surfaceFeatureExtract failed (non-fatal)"

    # Step 4: Block mesh
    echo "[4/5] Running blockMesh..."
    blockMesh > log.blockMesh 2>&1
    echo "  blockMesh complete"

    # Step 5: snappyHexMesh (parallel)
    echo "[5/5] Running snappyHexMesh ($NPROCS cores)..."
    cp system/decomposeParDict system/decomposeParDict.bak 2>/dev/null || true
    sed -i "s/numberOfSubdomains.*/numberOfSubdomains $NPROCS;/" system/decomposeParDict
    decomposePar > log.decomposePar 2>&1
    mpirun -np $NPROCS snappyHexMesh -parallel -overwrite > log.snappyHexMesh 2>&1
    reconstructParMesh -constant > log.reconstructParMesh 2>&1
    rm -rf processor*

    echo ""
    echo "  Mesh statistics:"
    grep -E "cells|Mesh OK|FAILED" log.snappyHexMesh | tail -5 || true

    echo "  Done: $key"
done

echo ""
echo "========================================================"
echo "All 6 cases meshed!"
echo "========================================================"
echo ""
echo "To view in ParaView:"
for key in "${!CASES[@]}"; do
    foam_file="$WSL_HOME/${CASES[$key]}/${CASES[$key]}.foam"
    echo "  $key:  $foam_file"
done
echo ""
echo "Or open the Windows-side .foam files:"
for key in "${!CASES[@]}"; do
    echo "  \\\\wsl\$\\Ubuntu\\home\\kirka\\${CASES[$key]}\\${CASES[$key]}.foam"
done
