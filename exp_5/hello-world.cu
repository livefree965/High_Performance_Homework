#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>
#include <numeric>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <cuda_runtime.h>
#include <time.h>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
using namespace std;

__global__ void hello()
{
	printf("hello from block %d %d thread %d  %d \n", blockIdx.x,blockIdx.y,threadIdx.x,threadIdx.y);
}

int main()
{
	dim3 dimBlock(8, 16);
	dim3 dimGrid(2, 4);
	hello << <dimGrid, dimBlock >> > ();
	return 0;
}