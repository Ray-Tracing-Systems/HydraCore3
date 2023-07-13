# HydraCore3
Modern rendering core: spec, vulkan (by kernel_slicer) and other

# Build
1) clone https://github.com/Ray-Tracing-Systems/kernel_slicer
2) clone this repo inside "apps" folder of kernel_slicer
3) if build only CPU version: disable Vulkan by '-DUSE_VULKAN=OFF'
4) if build both CPU and GPU versions, run kernel_slicer with 'HydraCore3/GLSL' config and then build solution normally with Cmake

# Development pipeline
1) Select/Find/Make a reference image to you feature
2) Implement it in renderer on CPU
3) run kernel_slicer to get GPU version and be sure that code succesefully transformed to shaders
  3.1) you may work with CPU build only, but this is long ...  
4) Add test to python script 'testing/run_tests.py':
  4.1) You have to clone https://github.com/Ray-Tracing-Systems/HydraAPI-tests
  4.2) You have to run tests from 'HydraAPI-tests' repo to generate scene files
  4.3) If you don't have access to some closed test repo which is used in test, please comment out such tests. Otherwise clone these repos.
  4.4) Currently you have to set two variables in script: 'PATH_TO_HYDRA2_TESTS' and 'PATH_TO_HYDRA3_SCENS' (the last one is currently closed)
5) Add you feature to specitication in specification: 'HydraAPI-tests/doc/doc_hydra_standart/hydra_spec.tex' 

