---
documentclass: physycomen
title:  "tools-piccante"
author: "Fedeli, Sgattoni, Sinigardi"
---

[![Build Status Master](https://travis-ci.org/ALaDyn/tools-piccante.svg?branch=master)](https://travis-ci.org/ALaDyn/tools-piccante "master")
[![Build status](https://ci.appveyor.com/api/projects/status/vipcofq9725tll9d?svg=true)](https://ci.appveyor.com/project/cenit/tools-piccante)


### Purpose
These tools have been written in order to analyse input and output of the `piccante` code.

### Installation
**CMake** and a **C++11** compatible compiler are required. To build the executable, clone the repository and then type  
```
mkdir build ; cd build ; cmake .. ; cmake --build . --target install
```
With CMake you can also deploy projects for the most common IDEs.  
Uses MPI, please install a compatible library on your system.
There is also a makefile, which is now deprecated.

### Synopsis
As an example, you can run one of the tools in this way
```
./newReader input_file [-cutx $value -integratex]
```
where *input_file* must be an existing and valid binary file produced by `piccante` and the parameters inside the `[]` are not mandatory but useful to get a specific insight.
