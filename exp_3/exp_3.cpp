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
short Mat_size;
short Row_char_count = 16;

short generate_mat(short row, short col, string filename) {
    ofstream outfile;
    outfile.open(filename);
    for (short i = 0; i < row; ++i) {
        for (short j = 0; j < col; ++j) {
            long beg = outfile.tellp();
            outfile << setw(5) << i << setw(5) << j << setw(5) << rand() % 100;
            outfile << endl;
        }
    }
    outfile.close();
}

void show_mat(short **mat, short size, string filename) {
    ofstream outfile;
    outfile.open(filename);
    for (short i = 0; i < size; ++i) {
        for (short j = 0; j < size; ++j) {
            outfile << setw(5) << mat[i][j] << " ";
        }
        outfile << endl;
    }
    outfile.close();
}

void show_vec(short *vec, short size, string filename) {
    ofstream outfile;
    outfile.open(filename);
    for (short i = 0; i < size; ++i) {
        outfile << setw(5) << vec[i] << " ";
        outfile << endl;
    }
    outfile.close();
}

void get_element(char *src, short *res) {
    /*
    * 获取一行矩阵的数据。
    */
    const char *d = " ";
    char *p;
    p = strtok(src, d);
    for (short i = 0; i < 3; i++) {
        res[i] = strtol(p, NULL, 10);
        p = strtok(NULL, d);
    }
    return;
}

short **read_mat(short beg_row, short beg_col, short part_size, string filename) {
    // beg_row：开始读取的行数
    //beg_col：开始读取的列数
    //part_size: 读取的矩阵小方块的大小(边长)
    //filename： 读取的文件名
    ifstream infile;
    infile.open(filename);
    short tmp[3]; //用来存储解析到的三元组数据
    char data[20]; //用来存储读取的每一行的数据
    short **res = new short *[Mat_size]; //用来保存矩阵数据，是一个二元数组
    for (short i = 0; i < Mat_size; ++i) {
        res[i] = new short[Mat_size];
    }
    for (short i = 0; i < part_size; ++i) {
        infile.seekg(((beg_row + i) * Mat_size + beg_col) * 16);
        //定位到矩阵元素所在的位置
        for (short j = 0; j < part_size; ++j) {
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

short *read_vec(short *res, short beg_row, short part_size, string filename) {
    /*
    res:向量保存的一维数组指针
    beg_row:开始的行数
    part_size:读取的行数大小
    filename:读取的文件名
    */
    ifstream infile;
    infile.open(filename);
    short tmp[3];
    char data[20];
    for (short i = 0; i < part_size; ++i) {
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
        Mat_size = short(pow(2, atoi(argv[1])));
    } else
        Mat_size = 64;
//    short **mat1_data = read_mat(0, 0, Mat_size, "mat.data");
//    short *vec1_data = new short[Mat_size];
//    read_vec(vec1_data, 0, Mat_size, "vec.data");
//    short *res1 = new short[Mat_size];
//    for (short l = 0; l < Mat_size; ++l) {
//        res1[l] = 0;
//    }
//    for (short i = 0; i < Mat_size; ++i) {
//        for (short j = 0; j < Mat_size; ++j) {
//            res1[i] += mat1_data[i][j] * vec1_data[j];
//        }
//    }
//    show_vec(res1, Mat_size, "standard.mtx");

//    short *tmp = read_vec(1, 2, "vec.data");
//    short **all = read_mat(0, 0, Mat_size, "mat.data");
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
    int              dims[2] = {short(sqrt(comm_size)), short(sqrt(comm_size))}; //我们拓扑成正方形的
    short process_dim = dims[0]; //拓扑的一行/列进程数
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &comm2d); //创建拓扑
    MPI_Cart_sub(comm2d, remainX, &commX); //创建列的子通信域
    MPI_Cart_sub(comm2d, remainY, &commY); //创建行的子通信域
    MPI_Comm_rank(commX, &rankX); //获取进程在当列的下标
    MPI_Comm_rank(commY, &rankY); //获取进程在当行的下标
    double begin = MPI_Wtime();
    short process_mat_size = short(ceil(Mat_size / sqrt(comm_size)));
//    Mat_size = process_mat_size * short(sqrt(comm_size));
//    short process_mat_size = short(Mat_size / sqrt(comm_size));
    short **mat_data = read_mat(rankX * process_mat_size, rankY * process_mat_size, process_mat_size, "mat.data");
    short *vec_data = new short[process_mat_size];
    if (rankX == rankY) {
        read_vec(vec_data, rankX * process_mat_size, process_mat_size, "vec.data");
    }
    MPI_Bcast(vec_data, process_mat_size, MPI_SHORT, short(my_rank % process_dim), commX);
    short *res = new short[process_mat_size];
    short *res_sum = new short[process_mat_size];
    short *all_sum = new short[Mat_size];
    for (short k = 0; k < process_mat_size; ++k) {
        res[k] = 0;
    }
    for (short i = 0; i < process_mat_size; ++i) {
        for (short j = 0; j < process_mat_size; ++j) {
            res[i] += mat_data[i][j] * vec_data[j];
        }
    }
    MPI_Reduce(res, res_sum, process_mat_size, MPI_SHORT, MPI_SUM, 0, commY);
    MPI_Gather(res_sum, process_mat_size, MPI_SHORT, all_sum, process_mat_size, MPI_SHORT, 0, commX);
    if (my_rank == 0) {
        double end = MPI_Wtime() - begin;
        cout << "Mat size: " << Mat_size << " Cores:" << comm_size << " Use time: " << end << endl;
//        show_vec(all_sum, Mat_size, "sum.mtx");
    }
//    prshortf("rank = %d;   X = %d;  Y = %d\n", my_rank, rankX, rankY);
//    MPI_Comm top_comm;
//    short Mat_size = short(sqrt(comm_size));
//    short dims[2] = {Mat_size, Mat_size};
//    short periods[2] = {0, 0};
//    for (short i = 0; i < 2; ++i) {
//        dims[i] = Mat_size;
//    }
//    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &top_comm);
//    short rankx, ranky;
//    if (my_rank == 1) {
//    }
    MPI_Finalize();
    return 0;
}