
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
const int N = 16;
const int grid_x_size = 5000;
const int grid_y_size = 5000;
const int block_x_size = 1;
const int block_y_size = 1;
const int mat_size = 5000;
using namespace std;
__global__ void gpuMatMultKernel(const float *a, const float *b, float *result)
{
	int threadId = (blockIdx.y * blockDim.y + threadIdx.y) * gridDim.x * blockDim.x
		+ blockIdx.x * blockDim.x + threadIdx.x;
	if (threadId < mat_size*mat_size)
	{
		int row = threadId / mat_size;
		int column = threadId % mat_size;

		result[threadId] = 0;
		for (int i = 0; i < N; i++)
		{
			result[threadId] += a[row * mat_size + i] * b[i * mat_size + column];
		}
	}
}

__global__ void hello(char *a, int *b)
{
	printf("hello %d %d %d \n", threadIdx.x, b[threadIdx.x], a[0]);
	a[threadIdx.x] += b[threadIdx.x];
	printf("hello %d %d %d \n", threadIdx.x, b[threadIdx.x], a[0]);
}

void print_mat(float* src) {
	for (int i = 0; i < mat_size; i++)
	{
		for (int j = 0; j < mat_size; j++)
		{
			printf("%f ", src[i*mat_size + j]);
		}
		printf("\n");
	}
}
void save_mat(float* src) {
	ofstream outfile("res.txt");
	for (int i = 0; i < mat_size; i++)
	{
		for (int j = 0; j < mat_size; j++)
		{
			outfile << i << " " << j << " " << src[i*mat_size + j] << endl;
		}
	}
	outfile.close();
}
int main()
{
	float* mat_a = (float*)malloc(sizeof(float)*mat_size*mat_size);
	float* mat_b = (float*)malloc(sizeof(float)*mat_size*mat_size);
	float* mat_c = (float*)malloc(sizeof(float)*mat_size*mat_size);
	for (int i = 0; i < mat_size; i++)
	{
		for (int j = 0; j < mat_size; j++)
		{
			mat_a[i*mat_size + j] = i - 0.1*j + 1;
			mat_b[i*mat_size + j] = 0.2*j - 0.1*i;
		}
	}
	//print_mat(mat_a);
	//print_mat(mat_b);
	float* mat_a_cuda;
	float* mat_b_cuda;
	float* mat_c_cuda;
	cudaMalloc((void**)&mat_a_cuda, sizeof(float)*mat_size*mat_size);
	cudaMalloc((void**)&mat_b_cuda, sizeof(float)*mat_size*mat_size);
	cudaMalloc((void**)&mat_c_cuda, sizeof(float)*mat_size*mat_size);
	cudaMemcpy(mat_a_cuda, mat_a, sizeof(float)*mat_size*mat_size, cudaMemcpyHostToDevice);
	cudaMemcpy(mat_b_cuda, mat_b, sizeof(float)*mat_size*mat_size, cudaMemcpyHostToDevice);
	dim3 dimBlock(block_x_size, block_y_size);
	dim3 dimGrid(grid_x_size, grid_y_size);
	gpuMatMultKernel <<<dimGrid, dimBlock >> > (mat_a_cuda, mat_b_cuda,mat_c_cuda);
	printf("finish\n");
	cudaMemcpy(mat_c, mat_c_cuda, sizeof(float)*mat_size*mat_size, cudaMemcpyDeviceToHost);
	//print_mat(mat_c);
	cudaFree(mat_a_cuda);
	cudaFree(mat_b_cuda);
	save_mat(mat_c);
	return 0;
}
