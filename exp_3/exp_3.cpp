//
// Created by Mat_size337259 on 11/7/18.
//
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
#include <pthread.h>
#include <semaphore.h>
#include <iomanip>
#include <memory.h>
#include <math.h>

using namespace std;
int Mat_size;
int Row_char_count = 16;

int generate_mat(int row, int col, string filename) {
    ofstream outfile;
    outfile.open(filename);
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            long beg = outfile.tellp();
            outfile << setw(5) << i << setw(5) << j << setw(5) << rand() % 100;
            outfile << endl;
        }
    }
    outfile.close();
}

void show_mat(int **mat, int size, string filename) {
    ofstream outfile;
    outfile.open(filename);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            outfile << setw(5) << mat[i][j] << " ";
        }
        outfile << endl;
    }
    outfile.close();
}

void show_vec(int *vec, int size, string filename) {
    ofstream outfile;
    outfile.open(filename);
    for (int i = 0; i < size; ++i) {
        outfile << setw(5) << vec[i] << " ";
        outfile << endl;
    }
    outfile.close();
}

void get_element(char *src, int *res) {
    /*
    * 获取一行矩阵的数据。
    */
    const char *d = " ";
    char *p;
    p = strtok(src, d);
    for (int i = 0; i < 3; i++) {
        res[i] = strtol(p, NULL, 10);
        p = strtok(NULL, d);
    }
    return;
}

int **read_mat(int beg_row, int beg_col, int part_size, string filename) {
    // beg_row：开始读取的行数
    //beg_col：开始读取的列数
    //part_size: 读取的矩阵小方块的大小(边长)
    //filename： 读取的文件名
    ifstream infile;
    infile.open(filename);
    int tmp[3]; //用来存储解析到的三元组数据
    char data[20]; //用来存储读取的每一行的数据
    int **res = new int *[Mat_size]; //用来保存矩阵数据，是一个二元数组
    for (int i = 0; i < Mat_size; ++i) {
        res[i] = new int[Mat_size];
    }
    for (int i = 0; i < part_size; ++i) {
        infile.seekg(((beg_row + i) * Mat_size + beg_col) * 16);
        //定位到矩阵元素所在的位置
        for (int j = 0; j < part_size; ++j) {
            //超出边界则设定值为0
            if (beg_row + i >= Mat_size || beg_col + j >= Mat_size) {
                res[i][j] = 0;
                continue;
            }
            infile.getline(data, 20); //读取一行数据
            get_element(data, tmp); //解析数据并保存到tmp中
            res[i][j] = tmp[2]; //将数据放入二元数组
        }
    }
    return res; //返回二元数组
}

int *read_vec(int *res, int beg_row, int part_size, string filename) {
    /*
    res:向量保存的一维数组指针
    beg_row:开始的行数
    part_size:读取的行数大小
    filename:读取的文件名
    */
    ifstream infile;
    infile.open(filename);
    int tmp[3];
    char data[20];
    for (int i = 0; i < part_size; ++i) {
        if (beg_row + i >= Mat_size) {
            //如果读取超出边界则设为0
            res[i] = 0;
            continue;
        }
        infile.seekg((beg_row + i) * Row_char_count);
        infile.getline(data, 20);
        get_element(data, tmp); //解析当行数据
        res[i] = tmp[2]; //存入一维数组
    }
    return res;
}

int main(int argc, char *argv[]) {
    if (argc == 2) {
        Mat_size = int(pow(2, atoi(argv[1])));
    } else
        Mat_size = 64;
//    int **mat1_data = read_mat(0, 0, Mat_size, "mat.data");
//    int *vec1_data = new int[Mat_size];
//    read_vec(vec1_data, 0, Mat_size, "vec.data");
//    int *res1 = new int[Mat_size];
//    for (int l = 0; l < Mat_size; ++l) {
//        res1[l] = 0;
//    }
//    for (int i = 0; i < Mat_size; ++i) {
//        for (int j = 0; j < Mat_size; ++j) {
//            res1[i] += mat1_data[i][j] * vec1_data[j];
//        }
//    }
//    show_vec(res1, Mat_size, "standard.mtx");

//    int *tmp = read_vec(1, 2, "vec.data");
//    int **all = read_mat(0, 0, Mat_size, "mat.data");
//    show_mat(all, 16, "all.data");
//    show_mat(res, 3);
    int my_rank, comm_size;
    int rankX, rankY;
    int ndims = 2; //维数
    int periods[2] = {0, 0}; //关乎进程的位置关系，我们默认第一个和最后一个不连接
    int reorder = 0;
    int remainX[2] = {1, 0}; //用于创建子通信域
    int remainY[2] = {0, 1}; //用于创建子通信域
    MPI_Comm comm2d; //笛卡尔拓扑所有进程的通信域
    MPI_Comm commX, commY;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    if (my_rank == 0) {
        srand(3);
        generate_mat(Mat_size, Mat_size, "mat.data");
        generate_mat(Mat_size, 1, "vec.data");
    }
    int dims[2] = {int(sqrt(comm_size)), int(sqrt(comm_size))}; //我们拓扑成正方形的
    int process_dim = dims[0]; //拓扑的一行/列进程数
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &comm2d); //创建拓扑
    MPI_Cart_sub(comm2d, remainX, &commX); //创建列的子通信域
    MPI_Cart_sub(comm2d, remainY, &commY); //创建行的子通信域
    MPI_Comm_rank(commX, &rankX); //获取进程在当列的下标
    MPI_Comm_rank(commY, &rankY); //获取进程在当行的下标
    double begin = MPI_Wtime();
    int process_mat_size = int(ceil(Mat_size / sqrt(comm_size)));
//    Mat_size = process_mat_size * int(sqrt(comm_size));
//    int process_mat_size = int(Mat_size / sqrt(comm_size));
    int **mat_data = read_mat(rankX * process_mat_size, rankY * process_mat_size, process_mat_size, "mat.data");
    int *vec_data = new int[process_mat_size];
    if (rankX == rankY) {
        read_vec(vec_data, rankX * process_mat_size, process_mat_size, "vec.data");
    }
    MPI_Bcast(vec_data, process_mat_size, MPI_INT, int(my_rank % process_dim), commX);
    int *res = new int[process_mat_size];
    int *res_sum = new int[process_mat_size];
    int *all_sum = new int[Mat_size];
    for (int k = 0; k < process_mat_size; ++k) {
        res[k] = 0;
    }
    for (int i = 0; i < process_mat_size; ++i) {
        for (int j = 0; j < process_mat_size; ++j) {
            res[i] += mat_data[i][j] * vec_data[j];
        }
    }
    MPI_Reduce(res, res_sum, process_mat_size, MPI_INT, MPI_SUM, 0, commY);
    MPI_Gather(res_sum, process_mat_size, MPI_INT, all_sum, process_mat_size, MPI_INT, 0, commX);
    if (my_rank == 0) {
        double end = MPI_Wtime() - begin;
        cout << "Mat size: " << Mat_size << " Cores:" << comm_size << " Use time: " << end << endl;
//        show_vec(all_sum, Mat_size, "sum.mtx");
    }
//    printf("rank = %d;   X = %d;  Y = %d\n", my_rank, rankX, rankY);
//    MPI_Comm top_comm;
//    int Mat_size = int(sqrt(comm_size));
//    int dims[2] = {Mat_size, Mat_size};
//    int periods[2] = {0, 0};
//    for (int i = 0; i < 2; ++i) {
//        dims[i] = Mat_size;
//    }
//    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &top_comm);
//    int rankx, ranky;
//    if (my_rank == 1) {
//    }
    MPI_Finalize();
    return 0;
}