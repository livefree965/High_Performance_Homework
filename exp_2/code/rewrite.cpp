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

struct Mat_Data {
    vector<int> **idx;
    vector<double> **val;
    int idx_size;
    int val_size;

    Mat_Data(int idx_size, int val_size) {
        this->idx_size = idx_size;
        this->val_size = val_size;
        this->idx = new vector<int> *[idx_size];
        this->val = new vector<double> *[val_size];
        for (int i = 0; i < idx_size; ++i) {
            this->idx[i] = new vector<int>;
        }
        for (int i = 0; i < val_size; ++i) {
            this->val[i] = new vector<double>;
        }
    }

    void show_array() {
        cout << "idx:" << endl;
        for (int i = 0; i < idx_size; ++i) {
            if (i > 10)
                break;
            for (int j = 0; j < idx[i]->size(); ++j) {
                cout << (*idx[i])[j] << " ";
            }
            cout << endl;
        }
        cout << endl;
        cout << "val:" << endl;
        for (int i = 0; i < val_size; ++i) {
            if (i > 10)
                break;
            for (int j = 0; j < val[i]->size(); ++j) {
                cout << (*val[i])[j] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
};

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

Mat_Data get_mat(const string &filename) {
    ifstream infile;
    char buffer[256];
    int element[3];
    infile.open(filename);
    infile.getline(buffer, 256, '\n');
    element[2] = get_element(buffer, element);
    struct Mat_Data mat_data(element[0], element[0]);
    while (!infile.eof()) {
        infile.getline(buffer, 256, '\n');
        double data = get_element(buffer, element);
        mat_data.idx[element[0] - 1]->push_back(element[1] - 1);
        mat_data.val[element[0] - 1]->push_back(data);
    }
    infile.close();
    return mat_data;
}

void save_mat(map<int, map<int, double >> mat, int row, int col) {
    ofstream out;
    out.open("res_mat.mtx");
    int all_size = 0;
    for (int i = 0; i < row; i++) {
        all_size += mat[i].size();
    }
    out << row << "\t" << col << "\t" << all_size << endl;
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            if (mat[i][j] != 0)
                out << i + 1 << "\t" << j + 1 << "\t" << mat[i][j] << endl;
        }
    }
    out.close();
}

int main(int argc, char *argv[]) {
    // string mat_a = "mat_data.mtx";
    // string mat_b = "vec_data.mtx";
    string mat_a = "test_mat.mtx";
    string mat_b = "test_vec.mtx";
    Mat_Data vec_data = get_mat(mat_b);
    Mat_Data mat_data = get_mat(mat_a);
    int mat_size[2], vec_size[2];
    update_info(mat_b, vec_size);
    update_info(mat_a, mat_size);
    double *res = new double[mat_size[0]];
//    serial verison
//    for (int row = 0; row < mat_size[0]; row++) {
//        for (int inner = 0; inner < mat_data.idx[row]->size(); inner++) {
//            if (inner == 0)
//                res[row] = (*mat_data.val[row])[inner] * (*vec_data.val[(*mat_data.idx[row])[inner]])[0];
//            else
//                res[row] += (*mat_data.val[row])[inner] * (*vec_data.val[(*mat_data.idx[row])[inner]])[0];
//        }
//    }
//    for (int j = 0; j < mat_size[0]; ++j) {
//        cout << res[j] << ' ';
//    }
    int my_rank, comm_size;
    MPI_Init(&argc, &argv);
    double beg = MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    int part = mat_size[0] / comm_size;
    if (part == 0)
        part++;
    if (my_rank != 0) {
        int row_end = (my_rank + 2) * part >= mat_size[0] ? mat_size[0] : (my_rank + 1) * part;
        for (int row = my_rank * part; row < row_end; row++) {
            for (int inner = 0; inner < mat_data.idx[row]->size(); inner++) {
                if (inner == 0)
                    res[row] = (*mat_data.val[row])[inner] * (*vec_data.val[(*mat_data.idx[row])[inner]])[0];
                else
                    res[row] += (*mat_data.val[row])[inner] * (*vec_data.val[(*mat_data.idx[row])[inner]])[0];
            }
        }
        for (int row = my_rank * part; row < row_end; row++) {
            MPI_Send(&res[row], mat_size[1], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    } else {
        printf("Rows :%d\n", mat_size[0]);
        printf("Core %d for row :%d-%d\n", my_rank, my_rank * part, (my_rank + 1) * part);
        for (int row = my_rank * part; row < (my_rank + 1) * part; row++) {
            for (int inner = 0; inner < mat_data.idx[row]->size(); inner++) {
                if (inner == 0)
                    res[row] = (*mat_data.val[row])[inner] * (*vec_data.val[(*mat_data.idx[row])[inner]])[0];
                else
                    res[row] += (*mat_data.val[row])[inner] * (*vec_data.val[(*mat_data.idx[row])[inner]])[0];
            }
        }
        double temp;
        for (int i = 1; i < comm_size; i++) {
            int row_end = (i + 2) * part >= mat_size[0] ? mat_size[0] : (i + 1) * part;
            printf("Core %d for row :%d-%d\n", i, i * part, row_end);
            for (int row = i * part; row < row_end; row++) {
                MPI_Recv(&res[row], mat_size[1], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
    if (my_rank == 0) {
        printf("Res: \n");
        for (int row = 0; row < mat_size[0]; row++) {
            printf("%lf \n", res[row]);
        }
        printf("Time use : %lf\n", MPI_Wtime() - beg);
    }
    MPI_Finalize();
    return 0;
}

