#!/bin/bash

cd external/LiteRT
make -j15
./render_app -benchmark $3 build 
cd ../..

#cd external/LiteRT
#bash scripts/slicer_preprocess_mesh_only.sh $1 $2
#bash scripts/slicer_build_shaders.sh $1 $2
#make -j15
#./render_app -benchmark $3 mesh_lod
#cd ../..

#bash scripts/slicer_preprocess_mesh_only.sh $1 $2
#bash scripts/slicer_build_shaders.sh $1 $2
#make -j15
#./bin-release/hydra -benchmark $3 mesh_lod

cd external/LiteRT
bash scripts/slicer_preprocess_mesh_only.sh $1 $2
bash scripts/slicer_build_shaders.sh $1 $2
make -j15
./render_app -benchmark $3 mesh 
cd ../..

bash scripts/slicer_preprocess_mesh_only.sh $1 $2
bash scripts/slicer_build_shaders.sh $1 $2
make -j15
./bin-release/hydra -benchmark $3 mesh

cd external/LiteRT
bash scripts/slicer_preprocess_grid_only.sh $1 $2
bash scripts/slicer_build_shaders.sh $1 $2
make -j15
./render_app -benchmark $3 sdf_grid 
cd ../..

bash scripts/slicer_preprocess_grid_only.sh $1 $2
bash scripts/slicer_build_shaders.sh $1 $2
make -j15
./bin-release/hydra -benchmark $3 sdf_grid

cd external/LiteRT
bash scripts/slicer_preprocess_octree_only.sh $1 $2
bash scripts/slicer_build_shaders.sh $1 $2
make -j15
./render_app -benchmark $3 sdf_octree 
cd ../..

bash scripts/slicer_preprocess_octree_only.sh $1 $2
bash scripts/slicer_build_shaders.sh $1 $2
make -j15
./bin-release/hydra -benchmark $3 sdf_octree

cd external/LiteRT
bash scripts/slicer_preprocess_framed_only.sh $1 $2
bash scripts/slicer_build_shaders.sh $1 $2
make -j15
./render_app -benchmark $3 sdf_frame_octree 
cd ../..

bash scripts/slicer_preprocess_framed_only.sh $1 $2
bash scripts/slicer_build_shaders.sh $1 $2
make -j15
./bin-release/hydra -benchmark $3 sdf_frame_octree

cd external/LiteRT
bash scripts/slicer_preprocess_sbs_only.sh $1 $2
bash scripts/slicer_build_shaders.sh $1 $2
make -j15
./render_app -benchmark $3 sdf_SBS-2-1 
./render_app -benchmark $3 sdf_SBS-2-2 
cd ../..

bash scripts/slicer_preprocess_sbs_only.sh $1 $2
bash scripts/slicer_build_shaders.sh $1 $2
make -j15
./bin-release/hydra -benchmark $3 sdf_SBS-2-1
./bin-release/hydra -benchmark $3 sdf_SBS-2-2

cd external/LiteRT
bash scripts/slicer_preprocess_svs_only.sh $1 $2
bash scripts/slicer_build_shaders.sh $1 $2
make -j15
./render_app -benchmark $3 sdf_SVS 
cd ../..

bash scripts/slicer_preprocess_svs_only.sh $1 $2
bash scripts/slicer_build_shaders.sh $1 $2
make -j15
./bin-release/hydra -benchmark $3 sdf_SVS