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
    fp = fopen("/root/CLionProjects/High_Performance_Homework/cmake-build-debug/test_data", "rb");
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
//    printf("rank %d sample %ld %ld %ld\n", my_rank, sample_data[0], sample_data[1], sample_data[2]);
//    printf("rank %d get %ld-%ld\n", my_rank, element_beg, element_end);
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
//    if (my_rank == 1) {
//        cout << "my pivots is " << pivots[0] << " " << pivots[1] << endl;
//    }
    //send pivot

    int *pivots_sum = new int[comm_size];
    pivots_sum[0] = 0;
    int find_pivot_index = 0;
    for (int i = 0; i < process_element_num; ++i) {
        if (datas[i] <= pivots[find_pivot_index])
            pivots_sum[find_pivot_index]++;
        else if (find_pivot_index == comm_size - 1)
            pivots_sum[comm_size - 1] = process_element_num - i + 1;
        else {
            find_pivot_index++;
            pivots_sum[find_pivot_index] = 1;
        }
    }

    int *pivots_newsum = new int[comm_size];
    MPI_Alltoall(pivots_sum, 1, MPI_INT, pivots_newsum, 1, MPI_INT, MPI_COMM_WORLD);
    int total_size = 0;
    for (int i = 0; i < comm_size; ++i) {
        total_size += pivots_newsum[i];
    }
    long *new_datas = new long[total_size];
    int *senddisp = new int[comm_size];
    int *recvdisp = new int[comm_size];
    senddisp[0] = 0;
    recvdisp[0] = 0;
    for (int i = 1; i < comm_size; ++i) {
        senddisp[i] = pivots_sum[i - 1] + senddisp[i - 1];
        recvdisp[i] = pivots_newsum[i - 1] + recvdisp[i - 1];
    }
    MPI_Alltoallv(datas, pivots_sum, senddisp, MPI_LONG, new_datas, pivots_newsum, recvdisp, MPI_LONG,
                  MPI_COMM_WORLD);
    //redistribute data
    long *sort_new_datas = new long[total_size];
    int *indexes = new int[comm_size];
    int *partitionEnds = new int[comm_size];
    memcpy(indexes, recvdisp, sizeof(int) * comm_size);
    memcpy(partitionEnds, recvdisp + 1, sizeof(int) * (comm_size - 1));
    partitionEnds[comm_size - 1] = total_size;
//    sort(new_datas, new_datas + total_size);

    for (int i = 0; i < total_size; ++i) {
        long lowest = 9223372036854775807;
        int ind = 0;
        for (int j = 0; j < comm_size; ++j) {
//            if (my_rank == 0)
//                cout << "indexs[j]" << indexes[j] << endl;
            if (indexes[j] < partitionEnds[j] && new_datas[indexes[j]] < lowest) {
                lowest = new_datas[indexes[j]];
                ind = j;
            }
        }
        sort_new_datas[i] = lowest;
        indexes[ind] += 1;
    }
    //merge_sort

    int *subListSize = new int[comm_size];
    MPI_Gather(&total_size, 1, MPI_INT, subListSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    delete[]datas;
    if (my_rank == 0) {
        recvdisp[0] = 0;
        for (int i = 1; i < comm_size; ++i) {
            recvdisp[i] = subListSize[i - 1] + recvdisp[i - 1];
        }
        datas = new long[element_num];
    }
    MPI_Gatherv(sort_new_datas, total_size, MPI_LONG, datas, subListSize, recvdisp, MPI_LONG, 0, MPI_COMM_WORLD);
    if (my_rank == 0) {
//        printf("pivots_sum %d %d %d\n", pivots_sum[0], pivots_sum[1], pivots_sum[2]);
//        printf("send disp %d %d %d\n", senddisp[0], senddisp[1], senddisp[2]);
//        printf("pivots_newsum %d %d %d\n", pivots_newsum[0], pivots_newsum[1], pivots_newsum[2]);
//        printf("recv disp %d %d %d\n", recvdisp[0], recvdisp[1], recvdisp[2]);
//        for (int i = 0; i < total_size; ++i) {
//            cout << partitionEnds[i] << endl;
//        }
//        for (int i = 0; i < comm_size; ++i) {
//            cout << recvdisp[i] << endl;
//        }
//        for (int i = 0; i < comm_size; ++i) {
//            cout << subListSize[i] << endl;
//        }
//        cout << element_num << endl;
        for (int j = 0; j < element_num; ++j) {
            cout << datas[j] << endl;
        }
//            cout << pivots_sum[i] << endl;
    }
//    int gather_sample_size = comm_size * comm_size;
//    cout << "hello" << endl;
    MPI_Finalize();
    return 0;
}

//
// Created by 16337259 on 11/27/18.
//

