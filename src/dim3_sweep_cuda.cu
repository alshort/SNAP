#include <cuda.h>
#include <cuda_runtime.h>

__global__
void add(int n, float *x, float *y)
{
  for (int i = 0; i < n; i++)
      y[i] = x[i] + y[i];
}