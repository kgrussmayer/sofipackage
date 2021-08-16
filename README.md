# SOFI Package

SOFI Package is a toolkit for super-resolution optical fluctuation imaging analysis in 2D and 3D.

More information can be found in the [SOFI documentation].

Please cite as: TODO

User manual and documentation is currently under development - stay tuned for updates

SOFI Package is developed in Matlab with the core calculation of cumulants optimized for 
GPU (implemented in CUDA)and for CPU (C++ code compiled as mex files). 

The code was tested in Matlab R2016b, under Windows 7, Windows 10 and Mac and requires Image Processing Toolbox. 
No additional software installation is required. If you want to use fast CUDA implementation, please make sure that 
you have a CUDA compatible graphics card and CUDA related software installed on your computer. For more details 
about using CUDA, please see notes below.

## Quick start

1. Clone the repository

```sh
$  git clone https://github.com/kgrussmayer/sofipackage.git
```

2. Download test data https://1drv.ms/f/s!Ak84leQ_7uk50BVykhJYSrw5Utgq and move 
them into ./data folder into the sofipackage repository.
There should be: 
* `./sofi/sofi2d`
* `./sofi/sofi3d`
* `./sofi/sofi3mc`
* `./sofi/sofibiplane`
3. Open sofipackage directory in MATLAB. Please make sure that it is set as 
"current folder" and it appears in the MATLAB path.  

4. Launch tests `./tests/test_expected_results.m`
All tests should pass, results and figures will be saved in
automatically created results folder. The analysis of all the test data should
take approx. 3min on a standard desktop computer (with at least 16GB RAM).

5. Run your own experiments. Example experiments for all 4 SOFI modalities are 
in `./experiments` this includes SOFI2D, SOFI3D, SOFI mutlicolor and SOFI biplane. 
For a new experiment with SOFI3D, create a new copy of `experiment_sofi3d.m` into a new folder 
for example `./experiment/sofi3d/20210630`. Edit configuration in the new experiment 
(change input files, output path etc.). Run the experiment while always keeping sofipackage 
as the "current folder" in MATLAB.

If you find a bug, have a feature request, or want to contribute to the SOFI Analysis Software, get in touch.

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
