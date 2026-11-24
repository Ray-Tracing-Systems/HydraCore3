#!/bin/bash
# Generates ISPC implementation of 'Integrator' with kernel_slicer (megakernel and wavefront variants, configs
# kmake_mega_ispc.json and kmake_wave_ispc.json) and builds 'bin-release/hydra_ispc_mega' and 'bin-release/hydra_ispc_wave'.
# Run them with '--ispc'. Needs 'ispc' in PATH.
#
# usage: bash ispc_build.sh <kernel_slicer folder> [mega|wave|all] [additional kslicer options, e.g. -ispc_wave_chunk 1024]

set -e
KSLICER_DIR=$(realpath "$1")
VARIANTS=${2:-all}
EXTRA_ARGS="${@:3}"
HYDRA_DIR=$(cd "$(dirname "$0")" && pwd)
if [ "$VARIANTS" == "all" ]; then VARIANTS="mega wave"; fi

for VARIANT in $VARIANTS; do
  echo "== ISPC, $VARIANT: kernel_slicer"
  (cd "$KSLICER_DIR" && ./cmake-build-release/kslicer "$HYDRA_DIR/kmake_${VARIANT}_ispc.json" $EXTRA_ARGS > "$HYDRA_DIR/kslicer_ispc_${VARIANT}.log" 2>&1) ||
    { echo "kslicer failed, see kslicer_ispc_${VARIANT}.log"; exit 1; }
  cd "$HYDRA_DIR"
  echo "== ISPC, $VARIANT: ispc"
  bash z_build_ispc.sh
  mkdir -p "ispc_$VARIANT"
  mv integrator_pt_ispc.cpp integrator_pt_kernels.ispc integrator_pt_kernels.h integrator_pt_kernels.o z_build_ispc.sh "ispc_$VARIANT/"
  echo "== ISPC, $VARIANT: cmake"
  cmake -S . -B "build_ispc_$VARIANT" -DCMAKE_BUILD_TYPE=Release -DUSE_ISPC=ON -DISPC_DIR="$HYDRA_DIR/ispc_$VARIANT" -DHYDRA_EXE_NAME="hydra_ispc_$VARIANT" > /dev/null
  cmake --build "build_ispc_$VARIANT" -j 16
done
