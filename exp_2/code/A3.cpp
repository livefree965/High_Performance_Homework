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
    res[2] = int(get_element(buffer, mat_info));
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

void save_mat(const string &filename, double *mat, int row, int col) {
    ofstream out;
    out.open(filename);
    out << row << "\t" << col << "\t" << row << endl;
    for (int i = 0; i < row; i++) {
        out << i + 1 << "\t" << 1 << "\t" << mat[i] << endl;
    }
    out.close();
}

int main(int argc, char *argv[]) {
    int my_rank, comm_size;
    int mat_size[3], vec_size[3];
    int row, idx;
    double val;
    double *res;
    MPI_Init(&argc, &argv);
    double beg = MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    int own_beg_row, own_end_row;
    if (my_rank != 0) {
        MPI_Recv(&mat_size, 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//        cout << "recv" << mat_size[0] << ' ' << mat_size[1] << ' ' << mat_size[2] << endl;
        MPI_Recv(&vec_size, 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//        cout << "recv" << vec_size[0] << ' ' << vec_size[1] << ' ' << vec_size[2] << endl;
        Mat_Data vec_data(vec_size[0], vec_size[0]);
        Mat_Data mat_data(mat_size[0], mat_size[0]);
        for (int i = 0; i < mat_data.idx_size; ++i) {
            for (int j = 0; j < mat_data.idx[i]->size(); ++j) {
                MPI_Scatter(&i, 1, MPI_INT, &row, 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Scatter(&(*mat_data.idx[i])[j], 1, MPI_INT, &idx, 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Scatter(&(*mat_data.val[i])[j], 1, MPI_INT, &val, 1, MPI_INT, 0, MPI_COMM_WORLD);
                mat_data.idx[row]->push_back(idx);
                mat_data.val[row]->push_back(val);

            }
        }
        for (int i = 0; i < vec_data.idx_size; ++i) {
            for (int j = 0; j < vec_data.idx[i]->size(); ++j) {
                MPI_Scatter(&i, 1, MPI_INT, &row, 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Scatter(&(*vec_data.idx[i])[j], 1, MPI_INT, &idx, 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Scatter(&(*vec_data.val[i])[j], 1, MPI_INT, &val, 1, MPI_INT, 0, MPI_COMM_WORLD);
                vec_data.idx[row]->push_back(idx);
                vec_data.val[row]->push_back(val);
            }
        }
        res = new double[mat_size[0]];
//        mat_data.show_array();
        /*
         * get distribute task
         */
        MPI_Recv(&own_beg_row, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&own_end_row, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        /*
         * below to cal
         */
        for (int row = own_beg_row; row < own_end_row; row++) {
            for (int inner = 0; inner < mat_data.idx[row]->size(); inner++) {
                if (inner == 0)
                    res[row] = (*mat_data.val[row])[inner] * (*vec_data.val[(*mat_data.idx[row])[inner]])[0];
                else
                    res[row] += (*mat_data.val[row])[inner] * (*vec_data.val[(*mat_data.idx[row])[inner]])[0];
            }
        }
        for (int row = own_beg_row; row < own_end_row; row++) {
            MPI_Send(&res[row], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    } else {
        string mat_a = "mat_data.mtx";
        string mat_b = "vec_data.mtx";
        Mat_Data vec_data = get_mat(mat_b);
        Mat_Data mat_data = get_mat(mat_a);
        update_info(mat_b, vec_size);
        update_info(mat_a, mat_size);
        res = new double[mat_size[0]];
        for (int l = 0; l < comm_size; ++l) {
            MPI_Send(&mat_size, 3, MPI_INT, l, 0, MPI_COMM_WORLD);
            MPI_Send(&vec_size, 3, MPI_INT, l, 0, MPI_COMM_WORLD);
        }
        for (int i = 0; i < mat_data.idx_size; ++i) {
            for (int j = 0; j < mat_data.idx[i]->size(); ++j) {
                MPI_Scatter(&i, 1, MPI_INT, &row, 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Scatter(&(*mat_data.idx[i])[j], 1, MPI_INT, &idx, 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Scatter(&(*mat_data.val[i])[j], 1, MPI_INT, &val, 1, MPI_INT, 0, MPI_COMM_WORLD);
            }
        }
        for (int i = 0; i < vec_data.idx_size; ++i) {
            for (int j = 0; j < vec_data.idx[i]->size(); ++j) {
                for (int k = 0; k < comm_size; ++k) {
                    MPI_Scatter(&i, 1, MPI_INT, &row, 1, MPI_INT, 0, MPI_COMM_WORLD);
                    MPI_Scatter(&(*vec_data.idx[i])[j], 1, MPI_INT, &idx, 1, MPI_INT, 0, MPI_COMM_WORLD);
                    MPI_Scatter(&(*vec_data.val[i])[j], 1, MPI_INT, &val, 1, MPI_INT, 0, MPI_COMM_WORLD);
                }
            }
        }
        /*
         * distribute task
         */
        vector<pair<int, int>> dis_row;
        int solv_size = mat_size[2] / comm_size;
        int beg_row = 0;
        int acl = 0;
        for (int j = 0; j < mat_data.idx_size; ++j) {
            acl += mat_data.idx[j]->size();
            if (acl >= solv_size) {
                acl = 0;
                dis_row.push_back(pair<int, int>(beg_row, j + 1));
                beg_row = j + 1;
            }
            if (dis_row.size() == comm_size - 1) {
                dis_row.push_back(pair<int, int>(beg_row, mat_size[0]));
                break;
            }
        }
        for (int k = 1; k < comm_size; ++k) {
            MPI_Send(&dis_row[k].first, 1, MPI_INT, k, 0, MPI_COMM_WORLD);
            MPI_Send(&dis_row[k].second, 1, MPI_INT, k, 0, MPI_COMM_WORLD);
        }
        own_beg_row = dis_row[0].first;
        own_end_row = dis_row[0].second;
        /*
         * below to cal
         */
        printf("Rows :%d\n", mat_size[0]);
        printf("Core %d for row :%d-%d row_sum: %d\n", my_rank, dis_row[my_rank].first, dis_row[my_rank].second,
               dis_row[my_rank].second - dis_row[my_rank].first);
        for (int row = own_beg_row; row < own_end_row; row++) {
            for (int inner = 0; inner < mat_data.idx[row]->size(); inner++) {
                if (inner == 0) {
                    res[row] = (*mat_data.val[row])[inner] * (*vec_data.val[(*mat_data.idx[row])[inner]])[0];
                } else {
                    res[row] += (*mat_data.val[row])[inner] * (*vec_data.val[(*mat_data.idx[row])[inner]])[0];
                }
            }
        }
        for (int i = 1; i < comm_size; i++) {
            printf("Core %d for row :%d-%d row_sum: %d\n", i, dis_row[i].first, dis_row[i].second,
                   dis_row[i].second - dis_row[i].first);
            for (int row = dis_row[i].first; row < dis_row[i].second; row++) {
                MPI_Recv(&res[row], 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
    if (my_rank == 0) {
        printf("Res: \n");
        for (int row = 0; row < mat_size[0]; row++) {
            if (row < 10) {
                printf("%lf \n", res[row]);
            } else
                break;
        }
        printf("Time use : %lf\n", MPI_Wtime() - beg);
        save_mat("A1_mat.mtx", res, mat_size[0], 1);
    }
    MPI_Finalize();
    return 0;
}