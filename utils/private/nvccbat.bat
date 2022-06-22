echo off
set PATH=C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Auxiliary\Build\;%PATH%
set PATH=C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.5\bin;%PATH%
call vcvars64.bat
nvcc %1 %2 %3 %4 %5 %6 %7 %8 %9 
