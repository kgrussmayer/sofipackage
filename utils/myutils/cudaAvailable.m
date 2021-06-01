function gpu=cudaAvailable
try
   gpu=parallel.gpu.GPUDevice.current();
   gpu=isa(gpu,'parallel.gpu.CUDADevice') && gpu.DeviceSupported;
catch msg
   gpu=false;
end