#include <stdio.h>
#include <iostream>
#include <time.h>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <random>
#include <vector>
#include <string.h>
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
        res[i] = atoi(p);
        p = strtok(NULL, d);
    }
    return atof(p);
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
        mat_data[element[0] - 1][element[1] - 1] = get_element(buffer, element);
    }
    infile.close();
    return mat_data;
}

void show_mat(const string &filename) {
    ifstream infile;
    char buffer[256];
    infile.open(filename);
    infile.getline(buffer, 256, '\n');
    int temp[3];
    get_element(buffer, temp);
    int row = temp[0];
    int col = temp[1];
    cout << filename << " ";
    printf("size: %d*%d\n", row, col);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            infile.getline(buffer, 256, '\n');
            printf("%f ", get_element(buffer, temp));
        }
        printf("\n");
    }
    infile.close();
}

int main(int argc, char *argv[]) {
    string mat_a = "mat_data.mtx";
    string mat_b = "vec_data.mtx";
    map<int, map<int, double >> vec_data = get_mat(mat_b);
    map<int, map<int, double >> mat_data = get_mat(mat_a);
    int mat_size[2], vec_size[2];
    update_info(mat_b, vec_size);
    update_info(mat_a, mat_size);
    map<int, map<int, double >>res;
    for (int row = 0; row < mat_size[0]; row++) {
        for (int col = 0; col < vec_size[1]; col++) {
            for (int inner = 0; inner < mat_size[1]; inner++) {
                res[row][col] += mat_data[row][inner] * vec_data[inner][col];
            }
        }
    }
    /*
    printf("Res: \n");
    for (int row = 0; row < mat_size[0]; row++) {
        for (int col = 0; col < vec_size[1]; col++) {
            for (int inner = 0; inner < mat_size[1]; inner++) {
                res[row][col] += mat_data[row][inner] * vec_data[inner][col];
            }
            printf("%lf ", res[row][col]);
        }
        printf("\n");
    }
    */
    return 0;
}

