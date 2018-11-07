//
// Created by 16337259 on 11/7/18.
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

using namespace std;

int main(int argc, char *argv[]) {
    int my_rank, comm_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    int row_size = int(sqrt(comm_size));
    int *dims = new int[int(sqrt(comm_size))];
//    MPI_Cart_create(MPI_COMM_WORLD, 2,dims,)
    cout << my_rank << endl;
    MPI_Finalize();
    return 0;
}