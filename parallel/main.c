#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "helpers.h"

int main(int argc, char *argv[]) {
    int rank, nproc, M, N, niter;
    float *pold, *pnew, *plim, *data, *masterdata;

    // initializeaza mediul mpi si aloca resursele necesare
    initialize_mpi(argc, argv, &rank, &nproc, &M, &N, &pold, &pnew, &plim, &data, &niter);

    double start_time, end_time;  // variabile pentru timp

    if (rank == 0) {
        start_time = MPI_Wtime();  // incepe masurarea timpului
    }

    // distribuie datele initiale intre procese
    distribute_data(rank, nproc, M, N, &masterdata, data, plim);

    // proceseaza datele in paralel si aplica algoritmi pe segmentele locale
    process_data(rank, nproc, M/nproc, N, niter, pold, pnew, plim);

    // colecteaza si scrie datele procesate intr-o imagine finala
    finalize_and_write_data(rank, nproc, M, N, M/nproc, N, niter, pnew, data, masterdata);

    if (rank == 0) {
        end_time = MPI_Wtime();  // termina masurarea timpului
        printf("Durata totala de executie: %f secunde.\n", end_time - start_time);
    }

    // elibereaza memoria alocata si finalizeaza mediul mpi
    free(pold);
    free(pnew);
    free(plim);
    free(data);
    MPI_Finalize();
    return 0;
}
