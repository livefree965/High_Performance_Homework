#include <stdio.h>
#include <iostream>
#include <time.h>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <random>
#include <vector>
#include <string.h>
#include <mpi.h>
#include <map>
using namespace std;

void generate_mat(const string &filename, int row, int col, int seed) {
	ofstream outfile;
	outfile.open(filename);
	outfile << row << "\t" << col << "\t" << row * col << "\n";
	srand(seed);
	for (int i = 0; i < row * col; i++) {
		outfile << i / col + 1 << "\t" << i % col + 1 << "\t" << rand() % 10;
		if (i != row * col - 1)
			outfile << "\n";
	}
	outfile.close();
}

double get_element(char *src, int *res) {
	/*
	* 获取一行矩阵的数据。
	*/
	const char *d = "\t";
	char *p;
	p = strtok(src, d);
	for (int i = 0; i < 2; i++) {
		res[i] = strtol(p, NULL, 10);
		p = strtok(NULL, d);
	}
	return strtod(p, NULL);
}

void update_info(const string &filename, int *res) {
	ifstream infile;
	int mat_info[3];
	char buffer[256];
	infile.open(filename);
	infile.getline(buffer, 256, '\n');
	get_element(buffer, mat_info);
	res[0] = mat_info[0];
	res[1] = mat_info[1];
}

map<int, map<int, double >> get_mat(const string &filename) {
	ifstream infile;
	char buffer[256];
	int element[3];
	infile.open(filename);
	infile.getline(buffer, 256, '\n');
	get_element(buffer, element);
	map<int, map<int, double >> mat_data;
	while (!infile.eof()) {
		infile.getline(buffer, 256, '\n');
		double data = get_element(buffer, element);
		mat_data[element[0] - 1][element[1] - 1] = data;
	}
	infile.close();
	return mat_data;
}

void show_mat(map<int, map<int, double >> mat) {
	for (map<int, map<int, double> >::iterator i = mat.begin(); i != mat.end(); i++) {
		for (map<int, double>::iterator j = (i->second).begin(); j != (i->second).end(); j++) {
			cout << j->second << ends;
		}
		cout << endl;
	}
}

void save_mat(map<int, map<int, double >> mat,int row,int col) {
	ofstream out;
	out.open("res_mat.mtx");
	int all_size = 0;
	for (int i = 0; i < row; i++)
	{
		all_size += mat[i].size();
	}
	out << row << "\t" << col << "\t" << all_size << endl;
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			if (mat[i][j]!=0)
				out << i+1 << "\t" << j+1 << "\t" << mat[i][j] << endl;
		}
	}
	out.close();
}

int main(int argc, char *argv[]) {
	string mat_a = "test_mat.mtx";
	string mat_b = "test_vec.mtx";
	map<int, map<int, double >> vec_data = get_mat(mat_b);
	map<int, map<int, double >> mat_data = get_mat(mat_a);
	int mat_size[2], vec_size[2];
	update_info(mat_b, vec_size);
	update_info(mat_a, mat_size);
	map<int, map<int, double >> res;
	//serial verison
	/*for (int row = 0; row < mat_size[0]; row++) {
		for (int col = 0; col < vec_size[1]; col++) {
			for (int inner = 0; inner < mat_size[1]; inner++) {
				res[row][col] += mat_data[row][inner] * vec_data[inner][col];
			}
		}
	}
	show_mat(res);*/
	int my_rank, comm_size;
	MPI_Init(&argc, &argv);
	double beg = MPI_Wtime();
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	int part = mat_size[0] / comm_size;
	if (my_rank != 0) {
		for (int row = my_rank * part; row < (my_rank + 1)*part; row++) {
			for (int col = 0; col < vec_size[1]; col++) {
				for (int inner = 0; inner < mat_size[1]; inner++) {
					res[row][col] += mat_data[row][inner] * vec_data[inner][col];
				}
			}
		}
		double temp;
		for (int row = my_rank * part; row < (my_rank + 1)*part; row++) {
			for (int col = 0; col < vec_size[1]; col++)
			{
				MPI_Send(&res[row][col], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			}
		}
	}
	else
	{
		for (int row = my_rank * part; row < (my_rank + 1)*part; row++) {
			for (int col = 0; col < vec_size[1]; col++) {
				for (int inner = 0; inner < mat_size[1]; inner++) {
					res[row][col] += mat_data[row][inner] * vec_data[inner][col];
				}
			}
		}
		printf("rows :%d\n", mat_size[0]);
		printf("Core for row :%d\n", part);
		//show_mat(mat_a);
		//show_mat(mat_b);
		double temp;
		for (int i = 1; i < comm_size; i++)
		{
			for (int row = i * part; row < (i + 1)*part; row++) {
				for (int col = 0; col < vec_size[1]; col++)
				{
					MPI_Recv(&temp, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					res[row][col] = temp;
				}
				
			}
		}
	}
	if (my_rank == 0)
	{
		printf("Res: \n");
		for (int row = 0; row < mat_size[0]; row++) {
			for (int col = 0; col < vec_size[1]; col++) {
				printf("%lf ", res[row][col]);
			}
			printf("\n");
		}
		printf("Time use : %lf\n", MPI_Wtime() - beg);
		save_mat(res, mat_size[0], vec_size[1]);
	}
	MPI_Finalize();
	return 0;
}

