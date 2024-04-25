#include "helpers.h"

void updateImage(float *pold, float *pnew, float *plim, int rows, int cols) {
    for (int i = 1; i < rows ; i++) {
        for (int j = 1; j < cols ; j++) {
            pnew[i * cols + j] = 0.25 * (pold[(i - 1) * cols + j] + pold[(i + 1) * cols + j] +
                                         pold[i * cols + (j - 1)] + pold[i * cols + (j + 1)] - plim[i * cols + j]);
        }
    }
}

void copyWithoutHalo(float *src, float *dest, int rows, int cols) {
    for (int i = 1; i < rows ; i++) {
      memcpy(dest + i * cols + 1, src + i * cols + 1, (cols - 2) * sizeof(float));
    }
}

void initialize_mpi(int argc, char *argv[], int *rank, int *nproc, int *M, int *N, float **pold, float **pnew, float **plim, float **data, int *niter) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, rank);
    MPI_Comm_size(MPI_COMM_WORLD, nproc);

    if (*rank == 0) {
        pgm_size("image_640x480.pgm", M, N);
    }
    MPI_Bcast(M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int MP = *M / *nproc;
    int NP = *N;

       if (argc != 2 || (atoi(argv[1]) <= 0)) {
        if (*rank == 0) {
            if (argc != 2) {
                printf("Utilizare: %s numar_de_iteratii\n", argv[0]);
                printf("Numar de iteratii nu a fost specificat, folosim valoarea implicita de 800.\n");
            } else {
                printf("Numar invalid de iteratii specificat: '%s'. Folosim valoarea implicita de 800.\n", argv[1]);
            }
        }
        *niter = 800; // seteaza valoarea implicita daca argumentul este invalid sau lipsa
    } else {
        *niter = atoi(argv[1]); // foloseste valoarea specificata ca argument
    }


    *pold = malloc((MP + 2) * (NP + 2) * sizeof(float));
    *pnew = malloc((MP + 2) * (NP + 2) * sizeof(float));
    *plim = malloc((MP + 2) * (NP + 2) * sizeof(float));
    *data = malloc(MP * NP * sizeof(float));

    if (!*pold || !*pnew || !*plim || !*data) {
        printf("Eroare la alocarea memoriei.\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
        MPI_Finalize();
        exit(-1);
    }

    for (int i = 0; i < MP + 2; i++) {
        for (int j = 0; j < NP + 2; j++) {
            int index = i * (NP + 2) + j;
            (*pold)[index] = 255.0f;
            (*pnew)[index] = 255.0f;
            (*plim)[index] = 255.0f;
        }
    }
}

void distribute_data(int rank, int nproc, int M, int N, float **masterdata, float *data, float *plim) {
    int MP = M / nproc;  // segmentele de randuri sunt distribuite uniform
    int NP = N;          

    if (rank == 0) {
        *masterdata = malloc(M * N * sizeof(float));  // doar masterul aloca masterdata si citeste imaginea
        if (*masterdata == NULL) {
            printf("Eroare la alocarea memoriei pentru masterdata.\n");
            MPI_Abort(MPI_COMM_WORLD, -1);
            MPI_Finalize();
            exit(-1);
        }
        pgm_read("image_640x480.pgm", *masterdata, M, N);
    }

    // utilizarea MPI_Scatter pentru a distribui datele
    MPI_Scatter(*masterdata, MP * NP, MPI_FLOAT, data, MP * NP, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // Copiem 'data' in 'plim', respectand conturul
    for (int i = 1; i <= MP; i++) {
        for (int j = 1; j <= NP; j++) {
            int index_plim = (i * (NP + 2)) + j; // NP + 2 pentru că avem un halo de 1 element pe fiecare latură
            int index_data = ((i - 1) * NP) + (j - 1);
            plim[index_plim] = data[index_data];
        }
    }
}

void process_data(int rank, int nproc, int MP, int NP, int niter, float *pold, float *pnew, float *plim) {
    MPI_Status status;

    for (int iter = 0; iter < niter; iter++) {
        if (rank > 0) {  // Nu primul proces
            MPI_Sendrecv(pold + (1 * (NP + 2) + 1), NP, MPI_FLOAT, rank - 1, 0,
                         pold + (0 * (NP + 2) + 1), NP, MPI_FLOAT, rank - 1, 1, MPI_COMM_WORLD, &status);
        }
        if (rank < nproc - 1) {  // Nu ultimul proces
            MPI_Sendrecv(pold + (MP * (NP + 2) + 1), NP, MPI_FLOAT, rank + 1, 1,
                         pold + ((MP + 1) * (NP + 2) + 1), NP, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, &status);
        }

        // Actualizarea imaginii
        updateImage(pold, pnew, plim, MP + 2, NP + 2);
        copyWithoutHalo(pnew, pold, MP + 2, NP + 2);
    }
}

void finalize_and_write_data(int rank, int nproc, int M, int N, int MP, int NP, int niter, float *pnew, float *data, float *masterdata) {
    // dupa iteratii, copiem datele procesate inapoi in data
    for (int i = 1; i <= MP; i++) {
        for (int j = 1; j <= NP; j++) {
            int index_data = ((i - 1) * NP) + (j - 1);
            int index_pold = (i * (NP + 2)) + j;
            data[index_data] = pnew[index_pold];
        }
    }

    // colectarea datelor procesate de la toate procesele în masterdata
    MPI_Gather(data, MP * NP, MPI_FLOAT, masterdata, MP * NP, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // Procesul master scrie imaginea finala
   if (rank == 0) {
        char filename[256];
        sprintf(filename, "processed_image_%d.pgm", niter);  // am adaugat numarul de iteratii in numele fisierului
        pgm_write(filename, masterdata, M, N);
        free(masterdata);
    }
}
