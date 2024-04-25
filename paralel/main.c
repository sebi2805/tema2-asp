#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "pgm_IO.h"  // Include header-ul pentru funcțiile de lucru cu fișiere PGM

void initialize(float *array, int rows, int cols, float value);
void updateImage(float *pold, float *pnew, float *plim, int rows, int cols);
void copyWithoutHalo(float *src, float *dest, int rows, int cols) ;

int main(int argc, char *argv[]) {
    int M, N;  // Dimensiunile imaginii
    float *data, *pold, *pnew, *plim;
    int niter = 10;  // Numărul default de iterații
    char *filename = "image_640x480.pgm";  // Numele fișierului default
    int rank, nprocs;  // Rank-ul procesului curent și numărul total de procese
    int MP, NP;  // Dimensiunile per proces

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Prelucrarea argumentelor liniei de comandă
    if (argc > 1) {
        niter = atoi(argv[1]);  // Primul argument este numărul de iterații
        if (argc > 2) {
            filename = argv[2];  // Al doilea argument este numele fișierului
        }
    }

    // Only the master process reads the size of the image
    if (rank == 0) {
        // Citirea dimensiunilor imaginii din fișier
        pgm_size(filename, &M, &N);
    }

    // Broadcast the dimensions to all processes
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Calculate local dimensions assuming M is divisible by nprocs
    MP = M / nprocs;
    NP = N;

    // Allocate memory for local buffers including the halo
    int extendedRows = MP + 2;
    int extendedCols = NP + 2;

    pold = (float*)malloc(extendedRows * extendedCols * sizeof(float));
    pnew = (float*)malloc(extendedRows * extendedCols * sizeof(float));
    plim = (float*)malloc(extendedRows * extendedCols * sizeof(float));
    data = (float*)malloc(MP * NP * sizeof(float));  // Only the core part without halo

    if (!pold || !pnew || !plim || !data) {
        fprintf(stderr, "Error in memory allocation\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // Initialize matrices with some default values (e.g., 255)
    initialize(pold, extendedRows, extendedCols, 255.0f);
    initialize(pnew, extendedRows, extendedCols, 255.0f);
    initialize(plim, extendedRows, extendedCols, 255.0f);

    // Rest of the program (scatter, compute, gather, etc.)
    if (rank == 0) {
    // Citirea datelor imaginii în 'data' principal
    pgm_read(filename, data, M, N);  // Se presupune că data este suficient de mare pentru a ține întreaga imagine
}

// Distribuția segmentelor de imagine către fiecare proces
MPI_Scatter(data, MP * NP, MPI_FLOAT, 
            pold + extendedCols + 1, // Offset pentru a lăsa spațiu pentru halo
            MP * NP, MPI_FLOAT, 
            0, MPI_COMM_WORLD);
            // Exemplu de actualizare a halo-urilor laterale
// Trimiterea și primirea coloanelor de margine
float* send_right = (float*)malloc((MP + 2) * sizeof(float)); // Buffer pentru coloana dreaptă
float* recv_left = (float*)malloc((MP + 2) * sizeof(float));  // Buffer pentru coloana stângă

// Copiem coloana dreaptă în buffer-ul de trimitere
for (int i = 0; i < MP + 2; i++) {
    send_right[i] = pold[i * extendedCols + NP];  // ultima coloană validă
}

// Send to right, receive from left
MPI_Sendrecv(send_right, MP + 2, MPI_FLOAT, (rank + 1) % nprocs, 0,
             recv_left, MP + 2, MPI_FLOAT, (rank - 1 + nprocs) % nprocs, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

// Actualizăm coloana stângă de halo cu datele primite
for (int i = 0; i < MP + 2; i++) {
    pold[i * extendedCols] = recv_left[i];
}

// Similar, actualizările pentru celelalte halo-uri

free(send_right);
free(recv_left);
// Procesarea locală a imaginii
// (Implementează logica specifică aici)
   for (int iter = 0; iter < niter; iter++) {
        updateImage(pold, pnew, plim, extendedRows, extendedCols);
        copyWithoutHalo(pnew, pold, extendedRows, extendedCols);
    }
// Colectarea datelor procesate înapoi la procesul principal
MPI_Gather(pnew + extendedCols + 1, MP * NP, MPI_FLOAT,
           data, MP * NP, MPI_FLOAT,
           0, MPI_COMM_WORLD);

if (rank == 0) {
    // Procesul principal poate acum rescrie imaginea procesată înapoi într-un fișier
    pgm_write("processed_image.pgm", data, M, N);
}

    // Finalize MPI
    MPI_Finalize();

    return 0;
}

void initialize(float *array, int rows, int cols, float value) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            array[i * cols + j] = value;
        }
    }
}
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
