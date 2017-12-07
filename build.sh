#!/bin/bash

# rm -rf build 
mkdir -p build
cd build
cmake .. -G "Ninja"
cmake --build . --target install
cd ..

