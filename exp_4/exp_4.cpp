#include <iostream>
#include <fstream>
#include <algorithm>
#include <mpi.h>
#include <stdio.h>

using namespace std;

void gen_data() {
    FILE *fp;
    fp = fopen("test_data", "wb");
    long size = 18;
    fwrite(&size, 8, 1, fp);
    long arr[18] = {48, 39, 6, 72, 91, 14, 69, 40, 89, 61, 12, 21, 84, 58, 32, 33, 72, 20};
    for (int i = 0; i < 18; ++i) {
        fwrite(&arr[i], 8, 1, fp);
    }
    fclose(fp);
}

void first_sample(long *datas, long data_size, long *output, long step) {
    sort(datas, datas + data_size);
    long output_pos = 0;
    for (long i = 0; i < data_size; i += step) {
        output[output_pos] = datas[i];
        output_pos += 1;
    }
}


void second_sample(long *datas, long data_size, long *output, int process_size) {
    long step = process_size;
    sort(datas, datas + data_size);
    long output_pos = 0;
    for (long i = 0; i < process_size - 1; i++) {
        output[output_pos] = datas[(((i + 1) * process_size) + (process_size / 2)) - 1];
        output_pos += 1;
    }
}

void read_datas(long *datas, long data_size, FILE *fp) {
    for (int i = 0; i < data_size; ++i) {
        fread(datas + i, 8, 1, fp);
    }
}

int main(int argc, char *argv[]) {
//    gen_data();

    int element_num, process_element_num;
    FILE *fp = NULL;
//    fp = fopen("/share/home/shared_dir/psrs_data", "rb");
    fp = fopen("/share/home/16337259/CLionProjects/High_Performance_Homework/cmake-build-debug/test_data", "rb");
    fread(&element_num, sizeof(long), 1, fp);
    //open file,read length


    MPI_Init(&argc, &argv);
    int my_rank, comm_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    process_element_num = element_num / comm_size;
    long element_beg = process_element_num * my_rank;
    long element_end = process_element_num * (my_rank + 1);
    long element_len = element_end - element_beg;
    fseek(fp, element_beg * 8, SEEK_CUR);
    long *datas = new long[element_len];
    read_datas(datas, element_len, fp);
    //read process data


    long *sample_data = new long[comm_size];
    first_sample(datas, element_len, sample_data, element_num / (comm_size * comm_size));
    printf("rank %d sample %d %d %d\n", my_rank, sample_data[0], sample_data[1], sample_data[2]);
    printf("rank %d get %ld-%ld\n", my_rank, element_beg, element_end);
    //sort and get sample



    long *sample_merge = new long[comm_size * comm_size];
    MPI_Gather(sample_data, comm_size, MPI_LONG, sample_merge, comm_size, MPI_LONG, 0, MPI_COMM_WORLD);
    //merge smaple to 0 process

    long *pivots = new long[comm_size - 1];
    if (my_rank == 0) {
        second_sample(sample_merge, comm_size * comm_size, pivots, comm_size);
    }
    //sort and get pivot

    MPI_Bcast(pivots, comm_size - 1, MPI_LONG, 0, MPI_COMM_WORLD);
    if (my_rank == 1) {
        cout << "my pivots is " << pivots[0] << " " << pivots[1] << endl;
    }
    //send pivot

    int *pivots_index = new int[comm_size];
    pivots_index[0] = 0;
    int find_pivot_index = 0;
    for (int i = 0; i < process_element_num; ++i) {
        if (datas[i] <= pivots[find_pivot_index])
            pivots_index[find_pivot_index]++;
        else if (find_pivot_index == comm_size - 1)
            pivots_index[comm_size - 1] = process_element_num - i + 1;
        else {
            find_pivot_index++;
            pivots_index[find_pivot_index] = 1;
        }
    }
    if (my_rank == 1) {
        for (int i = 0; i < comm_size; ++i) {
            cout << pivots_index[i] << endl;
        }
    }

//    int gather_sample_size = comm_size * comm_size;
//    cout << "hello" << endl;
    MPI_Finalize();
    return 0;
}

//
// Created by 16337259 on 11/27/18.
//

