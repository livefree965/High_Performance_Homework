#include <mpi.h>
#include <stdio.h>
#include <time.h>
#include <fstream>
#include <stdlib.h>
#include <string>
using namespace std;
int mat_info[3];
double get_element(char* src, int* res) {
	const char *d = "\t";
	char *p;
	p = strtok(src, d);
	for (int i = 0; i < 2; i++)
	{
		res[i] = atoi(p);
		p = strtok(NULL, d);
	}
	return atof(p);
}
void update_info(string filename, int* res) {
	ifstream infile;
	char buffer[256];
	infile.open(filename);
	infile.getline(buffer, 256, '\n');
	get_element(buffer, mat_info);
	res[0] = mat_info[0];
	res[1] = mat_info[1];
}
double** init_mat(int row, int col) {
	double** temp = (double**)malloc(sizeof(double *)*row);
	for (int i = 0; i < row; i++)
	{
		temp[i] = (double*)malloc(sizeof(double)*col);
	}
	return temp;
}
double** get_mat(string filename) {
	ifstream infile;
	char buffer[256];
	infile.open(filename);
	infile.getline(buffer, 256, '\n');
	get_element(buffer, mat_info);
	double** mat = init_mat(mat_info[0], mat_info[1]);
	while (!infile.eof())
	{
		infile.getline(buffer, 256, '\n');
		mat[mat_info[0] - 1][mat_info[1] - 1] = get_element(buffer, mat_info);
	}
	infile.close();
	return mat;
}
int main(int argc, char *argv[])
{
	string vec_file = "实验3的向量数据.mtx";
	double** vec_data = get_mat(vec_file);
	string mat_file = "实验3的矩阵数据.mtx";
	double** mat_data = get_mat(mat_file);
	int mat_size[2], vec_size[2];
	update_info(vec_file, vec_size);
	update_info(mat_file, mat_size);
	double** res = init_mat(mat_size[0], vec_size[1]);
	for (int row = 0; row < mat_size[0]; row++) {
		for (int col = 0; col < vec_size[1]; col++) {
			res[row][col]=0;
		}
	}
	for (int row = 0; row < mat_size[0]; row++) {
		for (int col = 0; col < vec_size[1]; col++) {  
			for (int inner = 0; inner < mat_size[1]; inner++) {
				res[row][col] += mat_data[row][inner] * vec_data[inner][col];
			}
			printf("%lf ", res[row][col]);
		}
		printf("\n");
	}
	printf("%d %d\n", mat_size[0], mat_size[1]);
	printf("%d %d\n", vec_size[0], vec_size[1]);
	/*
	for (int row = 0; row < 2; row++) {
		for (int col = 0; col < 1; col++) {
			// Multiply the row of A by the column of B to get the row, column of product.  
			for (int inner = 0; inner < 2; inner++) {
				product[row][col] += aMatrix[row][inner] * bMatrix[inner][col];
			}
			std::cout << product[row][col] << "  ";
		}
		std::cout << "\n";
	}*/
	/*
	double beg = MPI_Wtime();
	double a, b, h;
	int n;
	a = 0;
	b = 1000;
	n = 128;
	h = (b - a) / n;
	int local_n;
	double local_a, local_b;
	int my_rank, comm_size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	if (my_rank == 0) {
		local_n = n / comm_size;
		for (int recv = 1; recv < comm_size; recv++)
		{
			local_a = a + recv * local_n * h;
			local_b = local_a + local_n * h;
			MPI_Send(&local_n, 1, MPI_INT, recv, 0, MPI_COMM_WORLD);
			MPI_Send(&local_a, 1, MPI_DOUBLE, recv, 0, MPI_COMM_WORLD);
			MPI_Send(&local_b, 1, MPI_DOUBLE, recv, 0, MPI_COMM_WORLD);
		}
		local_a = a + my_rank * local_n * h;
		local_b = local_a + local_n * h;
	}
	else
	{
		MPI_Recv(&local_n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&local_a, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&local_b, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	local_n = n / comm_size;
	local_a = a + my_rank * local_n * h;
	local_b = local_a + local_n * h;
	double local_intergral = Trap(local_a, local_b, local_n, h);
	if (my_rank != 0) {
		MPI_Send(&local_intergral, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}

	else
	{
		double total_int = local_intergral;

		for (int source = 1; source < comm_size; source++)
		{
			MPI_Recv(&local_intergral, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			total_int += local_intergral;
		}
		printf("Total is %lf\n", total_int);
		printf("Time use : %lf\n", MPI_Wtime()-beg);
	}

	MPI_Finalize();
	*/
	return 0;
}
