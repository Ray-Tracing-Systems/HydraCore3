# Hydra3 Camera Plugin example

cam_pligin folder is an example of interaction between optics sim and hydra integrator.
In order to simplify interaction with the rendering system, the camera is implemented outside the renderer in separate class.
In addition, to keep things simply stupid, for this demo you have 2 different 'main' which is separated from normal hydra main.cpp. They are:
* main_with_cam.cpp     -- for CPU rendering with camera plugin
* main_with_cam_gpu.cpp -- for CPU rendering with camera plugin

So, it is supposed that you can take hydra Integrator source code for both CPU and GPU and integrate it to any custom application you like.

## Build (CPU and GPU) 
1) git clone https://github.com/Ray-Tracing-Systems/HydraCore3.git // clone this repo
2) git checkout gpu  // select gpu branch
3) clone submodules
   * git submodule init
   * git submodule init
4) build volk with Cmake inside folder 'external/volk':
   * cd external/volk
   * cmake .
   * make      // or in any other way 
6) clone 'comparisonrender' repo in any desired folder
  * in this example is it "../comparisonrender" with reletion to HydraCore3 directory
7) select some build directory 
  * mkdir build && cd build      
8) build solution normally with CMake for CPU:
  * cmake -DCMAKE_BUILD_TYPE=Release -DUSE_VULKAN=OFF -DCAM_PLUGIN=ON ..
  * make -j 8  // or in any other way 
  * please see VS code config, "CamTest(Tests/Blend/0002)" for example, you need "Tests/Blend/0002" scene from 'comparisonrender' repo
9) build solution normally with CMake for GPU:
  * cmake -DCMAKE_BUILD_TYPE=Release -DUSE_VULKAN=ON -DCAM_PLUGIN=ON ..
  * make -j 8  // or in any other way 
  * please see VS code config, "CamTest(Tests/Blend/0002)" for example, you need "Tests/Blend/0002" scene from 'comparisonrender' repo
  
  