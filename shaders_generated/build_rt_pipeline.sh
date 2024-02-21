#!/bin/sh
glslangValidator -V --target-env vulkan1.2 -S rgen CastSingleRayRGEN.glsl -o CastSingleRayRGEN.glsl.spv -DGLSL -I.. -I/home/frol/PROG/HydraRepos/HydraCore3/external/LiteScene -I/home/frol/PROG/HydraRepos/HydraCore3/external/LiteMath 
