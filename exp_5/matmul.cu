
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
using namespace std;
template <int BLOCK_SIZE> __global__ void
MatMul_CUDA(float *A, float *B, float *C, int mat_a_col, int mat_b_col)
{
	
	int block_x = blockIdx.x, block_y = blockIdx.y; // Block index
	int thread_x = threadIdx.x, thread_y = threadIdx.y; // Thread index

	
	int mat_a_begin = mat_a_col * BLOCK_SIZE * block_y;// Index of the first sub-matrix of A processed by the block
	int mat_a_end = mat_a_begin + mat_a_col - 1;// Index of the last sub-matrix of A processed by the block
	int mat_a_step = BLOCK_SIZE;// Step size used to iterate through the sub-matrices of A

	
	int mat_b_begin = BLOCK_SIZE * block_x; // Index of the first sub-matrix of B processed by the block
	int mat_b_step = BLOCK_SIZE * mat_b_col; // Step size used to iterate through the sub-matrices of B

	float c_sub_res = 0; // c_sub_res is used to store the element of the block sub-matrix

	// Loop  all the sub-matrices of A and B
	for (int a = mat_a_begin, b = mat_b_begin; a <= mat_a_end; a += mat_a_step, b += mat_b_step)
	{

		// A_sub,B_sub used to store the sub-matrix of A and B
		__shared__ float A_sub[BLOCK_SIZE][BLOCK_SIZE],B_sub[BLOCK_SIZE][BLOCK_SIZE];
		A_sub[thread_y][thread_x] = A[a + mat_a_col * thread_y + thread_x]; // Load A,B sub-matrices from device memory to share memory
		B_sub[thread_y][thread_x] = B[b + mat_b_col * thread_y + thread_x]; // Each thread loadsone element of each matrix

		__syncthreads(); // Make sure the matrix are loaded
		
		// Multiply the two matrix
		for (int k = 0; k < BLOCK_SIZE; ++k)
			c_sub_res += A_sub[thread_y][k] * B_sub[k][thread_x];

		// Synchronize to make sure that the preceding computation is finish
		__syncthreads();
	}

	// Write the block sub-matrix to device memory, each thread writes one element
	int c = mat_b_col * BLOCK_SIZE * block_y + BLOCK_SIZE * block_x;
	C[c + mat_b_col * thread_y + thread_x] = c_sub_res;
}

void show_mat(float**src, int col, int row) {
	/*
		src:address of the begin of matrix
		col:col size
		row:row size
	*/
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			printf("%f ", (*src)[i*col + j]);
		}
		printf("\n");
	}
}

void save_mat(float**src, int col, int row) {
	/*
		src:address of the begin of matrix
		col:col size
		row:row size
	*/
	ofstream outfile("mat_c.txt");
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			outfile << i << " " << j << " " << (*src)[i*col + j] << endl;
		}
	}
	outfile.close();
	printf("save finish");
}
void InitMatA(float **A, int col, int row) {
	/*
		A:address of begin of matrix
		col:col size
		row:row size
	*/
	*A = (float*)malloc(col*row * sizeof(float));
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			(*A)[i*col + j] = i - 0.1*j + 1;
		}
	}
}
void InitMatB(float **B, int col, int row) {
	/*
		B:address of begin of matrix
		col:col size
		row:row size
	*/
	*B = (float*)malloc(col*row * sizeof(float));
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			(*B)[i*col + j] = 0.2*j - 0.1*i;
		}
	}
}

void MatMul_CPU(float *A, float *B, float **C, int left_col, int left_row, int right_col ) {
	/*
		Mat a multi b serial version
		A: left matrix
		B: right matrix
		C: result matrix
		left_col: left matrix column size
		right_col: like above
		left row: like above
	*/
	clock_t beg, end;
	printf("Computing in CPU...\n");
	beg = clock();
	*C = (float*)malloc(left_row*right_col * sizeof(float));
	float *tmp = *C;
	float sum;
	int mat_a_index_base;
	for (int i = 0; i < left_row; i++) {
		for (int j = 0; j < right_col; j++) {
			sum = 0.0f;
			mat_a_index_base = i * left_col;
			for (int k = 0; k < left_col; k++) {
				sum += A[mat_a_index_base + k] * B[k*right_col + j];
			}
			tmp[i*right_col + j] = sum;
		}
	}
	end = clock();
	printf("Time= %.3f msec\n", (double)(end - beg)*1000 / CLOCKS_PER_SEC);
}

void MatMul_GPU(float *A, float *B, float **C, int block_size, dim3 &dims_A, dim3 &dims_B) {
	// Calculate malloc size for matrix A, B, C
	unsigned int size_A = dims_A.x * dims_A.y; //Size of matrix A
	unsigned int mem_size_A = sizeof(float) * size_A; //Memory size of matrix A
	unsigned int size_B = dims_B.x * dims_B.y;
	unsigned int mem_size_B = sizeof(float) * size_B;
	dim3 dimsC(dims_B.x, dims_A.y, 1);
	unsigned int mem_size_C = dimsC.x * dimsC.y * sizeof(float);
	*C = (float*)malloc(mem_size_C);

	//Above allocate memory for host,now malloc for device
	float *CUDA_A, *CUDA_B, *CUDA_C;

	cudaError_t error_A, error_B, error_C;
	error_A = cudaMalloc((void **)&CUDA_A, mem_size_A);
	error_B = cudaMalloc((void **)&CUDA_B, mem_size_B);
	error_C = cudaMalloc((void **)&CUDA_C, mem_size_C);

	//Allocate memory to store matrix in GPU

	if (error_A != cudaSuccess || error_B != cudaSuccess || error_C != cudaSuccess)
	{
		printf("cudaMalloc failed");
		exit(EXIT_FAILURE);
	}

	//Copy matrix from A to CUDA_A
	error_A = cudaMemcpy(CUDA_A, A, mem_size_A, cudaMemcpyHostToDevice);
	error_B = cudaMemcpy(CUDA_B, B, mem_size_B, cudaMemcpyHostToDevice);
	
	if (error_A != cudaSuccess || error_B != cudaSuccess)
	{
		printf("cudaMemcpy failed");
		exit(EXIT_FAILURE);
	}

	// To one thread,solve one ceil in matrix c
	//All thread num is matrix c size
	dim3 blocks(block_size, block_size);
	dim3 grid(dims_B.x / blocks.x, dims_A.y / blocks.y);

	// Begin compute
	printf("Computing in GPU...\n");

	//Warm up
	MatMul_CUDA<16> << < grid, blocks >> > (CUDA_A, CUDA_B, CUDA_C, dims_A.x, dims_B.x);
	printf("Warmup Finish\n");

	//Synchronize the gpu
	cudaDeviceSynchronize();

	//Record the time
	cudaEvent_t begin, end;
	cudaEventCreate(&begin);
	cudaEventRecord(begin, NULL);

	MatMul_CUDA<16> << < grid, blocks >> > (CUDA_A, CUDA_B, CUDA_C, dims_A.x, dims_B.x);

	printf("Finish\n");
	cudaEventCreate(&end);
	cudaEventRecord(end, NULL);

	// Wait for the stop event to complete
	cudaEventSynchronize(end);

	float time_total = 0.0f;
	cudaEventElapsedTime(&time_total, begin, end);

	// Compute and print the performance
	printf(
		"Time= %.3f msec \n",
		time_total);

	// Copy result from gpu
	error_C = cudaMemcpy(*C, CUDA_C, mem_size_C, cudaMemcpyDeviceToHost);

	if (error_C != cudaSuccess)
	{
		printf("cudaMemcpy (h_C,d_C) returned error %s (code %d), line(%d)\n",
			cudaGetErrorString(error_C), error_C, __LINE__);
		exit(EXIT_FAILURE);
	}

	// Free memory
	cudaFree(CUDA_A);
	cudaFree(CUDA_B);
	cudaFree(CUDA_C);

	// Use driver to clean up all state
	cudaDeviceReset();
}

bool Check_res(float *serial_result, float *cuda_result, int length) {
	bool correct = true;
	// test relative error by (<x, y>_cpu - <x,y>_gpu)  < eps
	double eps = 1.e-3; // machine zero
	double percent = 0.0f;
	for (int i = 0; i < length; i++)
	{
		double abs_err = fabs(cuda_result[i] - serial_result[i]);
		percent = abs_err / serial_result[i];
		if (percent > eps)
		{
			printf("Wrong at [%d] cpu=%.8f, gpu=%.8f percent=%.2f \n",
				i, serial_result[i], cuda_result[i], percent * 100);
			correct = false;
		}
	}
	printf("%s\n", correct ? "Result PASS" : "Result FAIL");
	return correct;
}
/**
 * Program main
 */
int main(int argc, char **argv)
{
	// Use Geforce 960M
	int devID = 0;

	// input five parameters
	int block_size = strtol(argv[1], NULL, 10);
	int mat_a_col = strtol(argv[2], NULL, 10);
	int mat_a_row = strtol(argv[3], NULL, 10);
	int mat_b_col = strtol(argv[4], NULL, 10);
	int mat_b_row = strtol(argv[5], NULL, 10);

	//see if we calculate right
	bool res;


	float *A, *B, *C1, *C2;
	//C1: calculate by gpu
	//C2: calculate by cpu

	cudaError_t error;
	error = cudaSetDevice(devID);

	if (error != cudaSuccess)
	{
		printf("cudaSetDevice wrong");
	}

	dim3 dims_A(mat_a_col, mat_a_row, 1);
	dim3 dims_B(mat_b_col, mat_b_row, 1);

	//to partition,we require that
	if (mat_a_col % block_size != 0 || mat_a_row % block_size != 0 ||
		mat_b_col % block_size != 0 || mat_b_row % block_size != 0) {
		printf("Dimension size must be dividable by block size!\n");
		exit(1);
	}

	//the condition to satisfiy a*b
	if (dims_A.x != dims_B.y)
	{
		printf("Outer matrix dimensions must be equal");
		exit(1);
	}

	// Init matrix
	InitMatA(&A, mat_a_col, mat_a_row);
	InitMatB(&B, mat_b_col, mat_b_row);
	
	
	//show_mat(&A, mat_a_w, mat_a_h);
	//printf("\n");
	//show_mat(&B, mat_b_w, mat_b_h);
	
	//GPU calculate
	MatMul_GPU(A, B, &C1, block_size, dims_A, dims_B);
	
	//CPU calculate
	printf("\n");
	MatMul_CPU(A, B, &C2, mat_a_col, mat_a_row, mat_b_col);
	
	//compare if the res is the same
	res = Check_res(C2, C1, mat_a_row * mat_b_col);
	
	save_mat(&C2, mat_a_col, mat_b_row);
	//free memory
	free(A); 
	free(B);
	free(C1);
	free(C2);
	exit(res);
}

