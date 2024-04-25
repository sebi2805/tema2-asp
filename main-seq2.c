#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pgm_IO.h"  // Include header-ul pentru funcțiile de lucru cu fișiere PGM

void initialize(float *array, int rows, int cols, float value);
void updateImage(float *pold, float *pnew, float *plim, int rows, int cols);
void copyWithoutHalo(float *src, float *dest, int rows, int cols);

int main(int argc, char *argv[]) {
    int M, N;  // Dimensiunile imaginii
    float *data, *pold, *pnew, *plim;
    int niter = 1;  // Numărul default de iterații
    char *filename = "image_640x480.pgm";  // Numele fișierului default

    // Prelucrarea argumentelor liniei de comandă
    if (argc > 1) {
        niter = atoi(argv[1]);  // Primul argument este numărul de iterații
        if (argc > 2) {
            filename = argv[2];  // Al doilea argument este numele fișierului
        }
    }

    // Citirea dimensiunilor imaginii din fișier
    pgm_size(filename, &M, &N);

    // Alocare memorie pentru matricea inițială și cele intermediare, incluzând halo
    int extendedRows = M + 2;
    int extendedCols = N + 2;

    pold = (float *)malloc(extendedRows * extendedCols * sizeof(float));
    pnew = (float *)malloc(extendedRows * extendedCols * sizeof(float));
    plim = (float *)malloc(extendedRows * extendedCols * sizeof(float));
    data = (float *)malloc(M * N * sizeof(float));

    if (!pold || !pnew || !plim || !data) {
        printf("Eroare la alocarea memoriei.\n");
        return -1;
    }

    // Inițializare matrici cu 255, inclusiv halo
    initialize(pold, extendedRows, extendedCols, 255.0f);
    initialize(pnew, extendedRows, extendedCols, 255.0f);
    initialize(plim, extendedRows, extendedCols, 255.0f);

    // Citirea datelor imaginii în 'data'
    pgm_read(filename, data, M, N);

    // Copiem 'data' în 'plim' respectând indexarea corectă și păstrând halo-ul la 255
    for (int i = 1; i <= M; i++) {
        for (int j = 1; j <= N; j++) {
            plim[i * extendedCols + (j)] = data[(i-1) * N + j-1];
        }
    }
     FILE *out = fopen("plim_matrix.txt", "w");
    if (out == NULL) {
        printf("Eroare la deschiderea fisierului!\n");
        return 1;
    }
    for (int i = 0; i < extendedRows; i++) {
        for (int j = 0; j < extendedCols; j++) {
            fprintf(out, "%f ", plim[i * extendedCols + j]);
        }
        fprintf(out, "\n");
    }
    fclose(out);

    // Iterarea algoritmului
    for (int iter = 0; iter < niter; iter++) {
        updateImage(pold, pnew, plim, extendedRows, extendedCols);
        copyWithoutHalo(pnew, pold, extendedRows, extendedCols);
    }

    // Scrierea imaginii procesate înapoi în fișierul PGM
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            data[i * N + j] = pnew[(i+1) * extendedCols + (j+1)];
        }
    }
    pgm_write("processed_image.pgm", data, M, N);

    // Eliberează memoria alocată
    free(pold);
    free(pnew);
    free(plim);
    free(data);

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
