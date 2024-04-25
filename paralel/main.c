#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "pgm_IO.h"

int main(int argc, char *argv[]) {
    const char *input_filename = "image_640x480.pgm";
    const char *output_filename = "output.pgm";
    int iterations[] = {800, 1000, 1200, 1400, 1600};
    int num_iterations = sizeof(iterations) / sizeof(iterations[0]);

    int M, N; // Dimensiunile imaginii
    float *masterdata, *data, *pold, *pnew, *plim;
    int rank, nproc;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    if (rank == 0) {
        // Procesul master citeste imaginea si initializeaza masterdata
        pgm_size(input_filename, &M, &N);
        masterdata = (float *)malloc(M * N * sizeof(float));
        pgm_read(input_filename, masterdata, M, N);
    }

    // Distribuirea dimensiunilor imaginii catre toate procesele
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Calcularea dimensiunilor blocurilor pentru fiecare proces
    int MP = M / nproc;
    int NP = N;

    // Initializarea matricilor locale
    data = (float *)malloc(MP * NP * sizeof(float));
    pold = (float *)malloc((MP + 2) * (NP + 2) * sizeof(float));
    pnew = (float *)malloc((MP + 2) * (NP + 2) * sizeof(float));
    plim = (float *)malloc((MP + 2) * (NP + 2) * sizeof(float));

    // Distribuirea segmentelor de date catre toate procesele
    MPI_Scatter(masterdata, MP * NP, MPI_FLOAT, data, MP * NP, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // Implementarea algoritmului de reconstructie a imaginii
    // (Aceasta parte va fi completata in continuare)

       // Inițializarea valorilor matricelor
    for (int i = 0; i < MP + 2; i++) {
        for (int j = 0; j < NP + 2; j++) {
            pold[i * (NP + 2) + j] = 255;  // Exemplu de inițializare
            pnew[i * (NP + 2) + j] = 255;
            plim[i * (NP + 2) + j] = (i > 0 && i <= MP && j > 0 && j <= NP) ? data[(i - 1) * NP + (j - 1)] : 255;
        }
    }

    // Procesarea iterativă a reconstrucției imaginii
    for (int n = 0; n < num_iterations; n++) {
        int niter = iterations[n];

        for (int iter = 0; iter < niter; iter++) {
            // Halo exchange la marginile blocului
            if (rank < nproc - 1) { // Trimite și primește de la procesul de sub
                MPI_Sendrecv(&pold[MP * (NP + 2) + 1], NP, MPI_FLOAT, rank + 1, 0,
                             &pold[(MP + 1) * (NP + 2) + 1], NP, MPI_FLOAT, rank + 1, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (rank > 0) { // Trimite și primește de la procesul de deasupra
                MPI_Sendrecv(&pold[1 * (NP + 2) + 1], NP, MPI_FLOAT, rank - 1, 0,
                             &pold[0 * (NP + 2) + 1], NP, MPI_FLOAT, rank - 1, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            // Aplică algoritmul de reconstrucție
            for (int i = 1; i <= MP; i++) {
                for (int j = 1; j <= NP; j++) {
                    int index = i * (NP + 2) + j;
                    pnew[index] = 0.25 * (pold[index - (NP + 2)] + pold[index + (NP + 2)] + pold[index - 1] + pold[index + 1] - plim[index]);
                }
            }

            // Copierea matricei pnew în pold pentru următoarea iterație
            for (int i = 1; i <= MP; i++) {
                for (int j = 1; j <= NP; j++) {
                    pold[i * (NP + 2) + j] = pnew[i * (NP + 2) + j];
                }
            }
        }
        printf("Procesul %d a finalizat iteratia %d.\n", rank, niter);
    }

    // Adunarea tuturor rezultatelor locale în matricea master
    MPI_Gather(data, MP * NP, MPI_FLOAT, masterdata, MP * NP, MPI_FLOAT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // Procesul master scrie imaginea finală
        pgm_write(output_filename, masterdata, M, N);
    }

    // Eliberarea resurselor și finalizarea MPI
    free(data);
    free(pold);
    free(pnew);
    free(plim);


    MPI_Finalize();

    return 0;
}
