function gpu=cudaAvailable
try
    gpu=parallel.gpu.GPUDevice.current();
    gpu=isa(gpu,'parallel.gpu.CUDADevice') && gpu.DeviceSupported;
catch
    gpu=false;
end
% eof