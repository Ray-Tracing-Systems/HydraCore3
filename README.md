# HydraCore3
Modern rendering core: spec, vulkan (by kernel_slicer) and other

# Build
1) clone https://github.com/Ray-Tracing-Systems/kernel_slicer
2) clone this repo inside "apps" folder of kernel_slicer
3) if build only CPU version: comment all Vulkan files in Cmake and then build with cmake
4) if build both CPU and GPU versions, run 'kslicer.sh/kslicer.bat' and then build solution normally with Cmake
