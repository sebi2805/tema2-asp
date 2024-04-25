#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "pgm_IO.h"  // Asigură-te că acest header este configurat pentru a lucra cu MPI
void updateImage(float *pold, float *pnew, float *plim, int rows, int cols) {
    for (int i = 1; i < rows ; i++) {
        for (int j = 1; j < cols ; j++) {
            pnew[i * cols + j] = 0.25 * (pold[(i - 1) * cols + j] + pold[(i + 1) * cols + j] +
                                         pold[i * cols + (j - 1)] + pold[i * cols + (j + 1)] - plim[i * cols + j]);
        }
    }
}
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, nproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    int M, N;  // Dimensiuni totale ale imaginii
    if (rank == 0) {
        pgm_size("image_640x480.pgm", &M, &N);
    }

    // Distribuie dimensiunile la toate procesele
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int MP = M / nproc;  // Segmentul de rânduri alocat fiecărui proces
    int NP = N;          // Toate coloanele sunt procesate de fiecare proces
    float *masterdata = NULL;
    int niter=1;
    if (rank == 0) {
        masterdata = malloc(M * N * sizeof(float));  // Doar procesul master alocă memoria pentru întreaga imagine
    }
    float (*pold)[NP + 2] = malloc((MP + 2) * (NP + 2) * sizeof(float));
    float (*pnew)[NP + 2] = malloc((MP + 2) * (NP + 2) * sizeof(float));
    float (*plim)[NP + 2] = malloc((MP + 2) * (NP + 2) * sizeof(float));
    float (*data)[NP] = malloc(MP * NP * sizeof(float));

    if (!pold || !pnew || !plim || !data) {
        printf("Eroare la alocarea memoriei.\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
        MPI_Finalize();
        return -1;
    }

    // Inițializare matrici cu 255, inclusiv halo
    for (int i = 0; i < MP + 2; i++) {
        for (int j = 0; j < NP + 2; j++) {
            pold[i][j] = 255.0f;
            pnew[i][j] = 255.0f;
            plim[i][j] = 255.0f;
        }
    }

    // Citirea și distribuția datelor din imagine
    if (rank == 0) {
        float *full_data = malloc(M * N * sizeof(float));
        pgm_read("image_640x480.pgm", full_data, M, N);

        // Distribuția segmentelor de date către fiecare proces
        for (int i = 1; i < nproc; i++) {
            MPI_Send(&full_data[i * MP * N], MP * N, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
        }

        memcpy(data, full_data, MP * N * sizeof(float));
        free(full_data);
    } else {
        MPI_Recv(data, MP * N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // Copiem 'data' în 'plim', respectând indexarea corectă și păstrând halo-ul la 255
    for (int i = 1; i <= MP; i++) {
        for (int j = 1; j <= NP; j++) {
            plim[i][j] = data[i-1][j-1];
        }
    }

    int iter;
    MPI_Status status;

    for (iter = 0; iter < niter; iter++) {
        // Halo exchange
        if (rank > 0) {  // Nu primul proces
            MPI_Sendrecv(&pold[1][1], NP, MPI_FLOAT, rank - 1, 0,
                        &pold[0][1], NP, MPI_FLOAT, rank - 1, 1, MPI_COMM_WORLD, &status);
        }
        if (rank < nproc - 1) {  // Nu ultimul proces
            MPI_Sendrecv(&pold[MP][1], NP, MPI_FLOAT, rank + 1, 1,
                        &pold[MP + 1][1], NP, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, &status);
        }

        // Aplicați algoritmul de reconstrucție pe datele locale
        updateImage(pold, pnew, plim, MP + 2, NP + 2);

        // Copiați rezultatele în pold pentru următoarea iterație
        for (int i = 1; i <= MP; i++) {
            memcpy(&pold[i][1], &pnew[i][1], NP * sizeof(float));
        }
    }

    MPI_Gather(data, MP * NP, MPI_FLOAT, masterdata, MP * NP, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // Procesul master scrie imaginea finală
    if (rank == 0) {
        pgm_write("processed_image1.pgm", masterdata, M, N);
        free(masterdata);
    }

    free(data);

    // Restul procesării și MPI_Finalize
    MPI_Finalize();
    return 0;
}
