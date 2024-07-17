#!/bin/bash
# slicer_execute absolute/path/to/slicer/folder absolute/path/to/slicer/executable
#e.g. bash slicer_execute.sh ~/kernel_slicer/ ~/kernel_slicer/kslicer
start_dir=$PWD
cd $1

$2 $start_dir/integrator_pt.cpp \
   $start_dir/integrator_pt_lgt.cpp \
   $start_dir/integrator_pt_mat.cpp \
   $start_dir/integrator_rt.cpp \
   $start_dir/integrator_spectrum.cpp \
   $start_dir/external/LiteRT/BVH/BVH2Common.cpp \
-mainClass Integrator \
-composInterface ISceneObject \
-composImplementation BVHRT \
-stdlibfolder $PWD/TINYSTL \
-pattern rtv \
-I$PWD/TINYSTL                                  ignore \
-I$start_dir/external/LiteMath                  ignore \
-I$start_dir/external/LiteScene                 ignore \
-I$start_dir/cam_plugin                        process \
-I$start_dir/external                          process \
-I$start_dir/external/LiteRT                   process \
-I$start_dir/external/LiteRT/BVH               process \
-I$start_dir/external/LiteRT/sdfScene          process \
-I$start_dir/external/CrossRT                  process \
-shaderCC glsl \
-enable_motion_blur 0 \
-gen_gpu_api 0 \
-megakernel 1 \
-DPUGIXML_NO_EXCEPTIONS -DKERNEL_SLICER -DUSE_LITERT -v \
-DDISABLE_SDF_GRID \
-DDISABLE_SDF_FRAME_OCTREE \
-DDISABLE_SDF_SBS \
-DDISABLE_SDF_SVS \
-DDISABLE_SDF_HP \
-DDISABLE_RF_GRID \
-DDISABLE_GS_PRIMITIVE
#-DDISABLE_SDF_GRID \