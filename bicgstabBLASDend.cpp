/**************************************************************************
bicgstabBLASDend - double precision
end bicgstabBLASD
TWS, March 2011
Cuda 10.1 Version, August 2019
**************************************************************************/
//#include <shrUtils.h>
//#include <cutil_inline.h>
#include <cuda_runtime.h>
#include "nrutil.h"

void bicgstabBLASDend(int n)
{
	extern double* h_x, * h_b, * h_a, * h_rs;
	extern double* d_a, * d_res, * d_x, * d_b, * d_r, * d_rs, * d_v, * d_s, * d_t, * d_p, * d_er;

	cudaFree(d_res);
	cudaFree(d_b);
	cudaFree(d_x);
	cudaFree(d_er);
	cudaFree(d_p);
	cudaFree(d_t);
	cudaFree(d_s);
	cudaFree(d_v);
	cudaFree(d_rs);
	cudaFree(d_r);
	cudaFree(d_a);

	free_dvector(h_rs, 0, n - 1);
	free_dvector(h_b, 0, n - 1);
	free_dvector(h_x, 0, n - 1);
	free_dvector(h_a, 0, n * n - 1);
}
