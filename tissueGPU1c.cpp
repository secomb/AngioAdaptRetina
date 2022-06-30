/**************************************************************************
tissueGPU1.cpp
program to call tissueGPU1.cu on GPU
TWS, December 2011
Cuda 10.1 Version, August 2019
**************************************************************************/
//#include <shrUtils.h>
//#include <cutil_inline.h>
#include <cuda_runtime.h>

extern "C" void tissueGPU1(int* tisspoints, float* dtt000, float* pt000, float* qt000, int nnt, int blockstep, int stepGPU, int iblockstep);

void tissueGPU1c(int stepGPU, int blockstep, int iblockstep)
{
	extern int nnt, * d_tisspoints;
	extern float* pt000, * qt000, * d_qt000, * d_pt000, * d_dtt000;
	cudaMemcpy(d_qt000, qt000, nnt * sizeof(float), cudaMemcpyHostToDevice);
	tissueGPU1(d_tisspoints, d_dtt000, d_pt000, d_qt000, nnt, stepGPU, blockstep, iblockstep);
	cudaMemcpy(pt000, d_pt000, nnt * sizeof(float), cudaMemcpyDeviceToHost);
}
