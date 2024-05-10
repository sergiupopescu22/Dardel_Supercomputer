
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char* argv[]) {

    int rank, size, i, provided;

    // number of cells (global)
    int nxc = 128; // make sure nxc is divisible by size
    double L = 2 * 3.1415; // Length of the domain

    MPI_Request request1, request2, request3, request4;
    MPI_Status status;


    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // number of nodes (local to the process): 0 and nxn_loc-1 are ghost cells 
    int nxn_loc = nxc / size + 3; // number of nodes is number cells + 1; we add also 2 ghost cells
    double L_loc = L / ((double)size);
    double dx = L / ((double)nxc);

    // define out function
    double* f = calloc(nxn_loc, sizeof(double)); // allocate and fill with z
    double* f_cos = calloc(nxn_loc, sizeof(double)); // allocate and fill with z
    double* dfdx = calloc(nxn_loc, sizeof(double)); // allocate and fill with z

    for (i = 1; i < (nxn_loc - 1); i++)
    {
        f[i] = sin(L_loc * rank + (i - 1) * dx);
        f_cos[i] = cos(L_loc * rank + (i - 1) * dx);
    }


    // need to communicate and fill ghost cells f[0] and f[nxn_loc-1]
    // communicate ghost cells

    // we will have to send to rank-1 the value of f[1]
    // we will have to send to rank+1 the value of f[nxn_loc-2]
    // we will use tag 1 for the right value, and tag 0 for the left value

    // HERE WE SEND THE GHOST CELLS
    //-----------------------------------------------------------------------------
    if (rank != 0)
        MPI_Send(f + 1, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
    else
        MPI_Send(f + 1, 1, MPI_DOUBLE, size - 1, 0, MPI_COMM_WORLD);

    if (rank != size - 1)
        MPI_Send(f + nxn_loc - 2, 1, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
    else
        MPI_Send(f + nxn_loc - 2, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    //-----------------------------------------------------------------------------


    // HERE WE RECEIVE THE GHOST CELLS
    //-----------------------------------------------------------------------------
    double right_ghost_cell_recv = 0;
    double left_ghost_cell_recv = 0;
    if (rank != size - 1)
        MPI_Recv(&right_ghost_cell_recv, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &status);
    else
        MPI_Recv(&right_ghost_cell_recv, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    //MPI_Wait(&request3, &status);

    if (rank != 0)
        MPI_Recv(&left_ghost_cell_recv, 1, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &status);
    else
        MPI_Recv(&left_ghost_cell_recv, 1, MPI_DOUBLE, size - 1, 1, MPI_COMM_WORLD, &status);
    //MPI_Wait(&request4, &status);
    //-----------------------------------------------------------------------------

    printf("\nMy rank %d of %d\n", rank, size);
    printf("received left ghost cell = %f\n", left_ghost_cell_recv);
    printf("received right ghost cell = %f\n", right_ghost_cell_recv);

    f[0] = left_ghost_cell_recv;
    f[nxn_loc - 1] = right_ghost_cell_recv;
    // here we finish the calculations

    // calculate first order derivative using central difference
    // here we need to correct value of the ghost cells!
    for (i = 1; i < (nxn_loc - 1); i++)
        dfdx[i] = (f[i + 1] - f[i - 1]) / (2 * dx);

    printf("Here are my values for f including ghost cells\n");

    for (i = 0; i < nxn_loc; i++)
        printf("%f\n", f[i]);
    printf("\n");

    printf("NOW WE WILL COMPARE THE RESULTS:\n");
    printf("%f  %f\n", f_cos[1], dfdx[1]);
    printf("%f  %f\n", f_cos[nxn_loc - 2], dfdx[nxn_loc - 2]);

    printf("\n");

    MPI_Finalize();
}

