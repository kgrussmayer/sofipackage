Software for multicolor SOFI imaging. More details are described in the paper: Spectral Cross-Cumulants for Multicolor Super-resolved SOFI Imaging

Developed and tested in Matlab R2016b, under Windows 7 and Windows 10. No additional software installation is required.

To launch the calculation:

1) Please make sure that test data are in the same folder together with the code. 
   There should be the following folder structure: 

funcs    	... 		mutlicolor sofi interface and auxiliary functions
test_data 	... 		simulated data for testing
test_files 	... 		expected results saved as mat files to test if the processing get results as expected.
config.m  	... 		configuration file where user can specify names of files to be processed and processing settings
multicolor_sofi_main.m		main file

2) Run multicolor_sofi_main.m


All tests should pass, results and figures will be saved in the automatically created results folder. The analysis of the test data should have a runtime of less than 2min on a standard desktop computer.


/////////////////////////////////////////////////////////////////
Using GPU algorithms

 - SOFI algorithm can be accelerated by GPU. It requires that the computer has a CUDA-enabled NVIDIA GPU (http://ch.mathworks.com/discovery/matlab-gpu.html).
If you have the Matlab Parallel Processing Toolbox and a NVIDIA graphics card, write the following command on the "Command Window":

>> gpuDevice

It should display:

ans = 
	CUDADevice with properties:

with a list of properties. Note the 'ComputeCapibility' of your graphics card.


 - Compilation for your GPU card: 
Step 1: note the 'ComputeCapibility' of your graphics card which is for example '2.0' in the case of NVIDIA Geforce GTX 480. 
The compute capability of a graphics card can also be found here: https://developer.nvidia.com/cuda-gpus
Step 2: Set the Matlab current folder to "multicolor_sofi\funcs\private". The tab "Current folder" should display files including 'gpu.cu', 'gpu.ptx', 'nvcc.m' and 'nvccbat.bat'.
Step 3: Execute the following command on the "Command Window": 
>> nvcc -arch=sm_20 -ptx gpu.cu
This command launches the NVIDIA compiler to recompile gpu.cu. Make sure that the -arch option is set to the compute capability of your CUDA-capable graphics card. 
In the example displayed above, it is 2.0 (the compute capibility of NVIDIA Geforce GTX 480).

/////////////////////////////////////////////////////////////////
Using CPU algorithms

- we recommend using the GPU version of the SOFI algorithm (see above). The mex files included in the current software package are compiled for Win32bit systems 
and we provide the original files that the user can compile for their computer system.
