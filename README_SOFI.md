# SOFI Package

We are currently in the process of cleaning up this repository 
and adding a proper user manual and testing - stay tuned for updates

TODO: update/adapt all links and descriptions etc

TODO: insert zenodo/doi

SOFI Package is a toolkit for super-resolution optical fluctuation imaging analysis in 2D and 3D.

More information can be found in the [SOFI documentation].

Please cite as: TODO

## Quick start

SOFI Package is developed in Matlab with the core cumulant
calculation provided for GPU (CUDA) and for CPU (C++ code implemented via Mexx
files). 

The code was tested in Matlab R2016b, under Windows 7, Windows 10 and Mac and requires the XXX toolboxes. No additional software installation is required.

The following commands will clone the repository

```sh
$ git clone https://github.com/acts-project/acts <source-dir>
```

To launch the calculation:

1.  Please make sure that test data are in the same folder together with the code. 
2.  Run `multicolor_sofi_main.m`
    
    All tests should pass, results and figures will be saved in the
    automatically created results folder. The analysis of the test data should
    have a runtime of TODO on a standard desktop computer.

If you find a bug, have a feature request, or want to contribute to the SOFI Analysis Software, TODO.

## Using GPU algorithms

-   the cumulant calculation can be accelerated by GPU. It requires that the
    computer has a CUDA-enabled NVIDIA GPU
    (http://ch.mathworks.com/discovery/matlab-gpu.html).

    If you have the Matlab Parallel Processing Toolbox and a NVIDIA graphics card, write the following command on the "Command Window":

    ```MATLAB
    gpuDevice
    ```

It should display:

```MATLAB
ans = 
	CUDADevice with properties:
```

with a list of properties. Note the 'ComputeCapibility' of your graphics card.


 - Compilation for your GPU card: 
Step 1: note the 'ComputeCapibility' of your graphics card which is for example '2.0' in the case of NVIDIA Geforce GTX 480. 
The compute capability of a graphics card can also be found here: https://developer.nvidia.com/cuda-gpus
Step 2: Set the Matlab current folder to "multicolor_sofi\funcs\private". The tab "Current folder" should display files including 'gpu.cu', 'gpu.ptx', 'nvcc.m' and 'nvccbat.bat'.
Step 3: Execute the following command on the "Command Window": 

```sh
>> nvcc -arch=sm_20 -ptx gpu.cu
```

This command launches the NVIDIA compiler to recompile gpu.cu. Make sure that the -arch option is set to the compute capability of your CUDA-capable graphics card. 
In the example displayed above, it is 2.0 (the compute capibility of NVIDIA Geforce GTX 480).

## Using CPU algorithms

We recommend using the GPU version of the SOFI algorithm (see above). The Mexx
files included in the current software package are compiled for Win32bit
systems and we provide the original files that the user can compile for their
computer system.

## Repository organization

The repository contains all code of the SOFI Analysis Software project.

TODO: adapt to the final structure and describe

## Authors and license

Contributors to the SOFI Analysis Software project are: Marcel Leutenegger,
Stefan Geissbühler, Tomas Lukes, Kristin Grußmayer, Adrien Descloux, Vytautas
Navikas

The SOFI Package project is published under the terms of the GNU general public licence v3. A
copy of the license can be found in the [LICENSE](LICENSE) file or at https://www.gnu.org/licenses/gpl-3.0.de.html.

The SOFI Package project contains copies of the following external
packages:

TODO: insert
