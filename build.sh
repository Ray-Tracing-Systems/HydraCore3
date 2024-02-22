#!/bin/bash


export CC=/usr/bin/clang
export CXX=/usr/bin/clang++
cmake -G "Ninja" -DCMAKE_BUILD_TYPE=RelWithDbgInfo -B ./bin-release 
cd ./bin-release
ninja 