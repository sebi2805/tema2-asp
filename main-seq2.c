#include <stdlib.h>
#include <stdio.h>
#include "pgm_IO.h"  // Asigură-te că ai acest header și funcțiile implementate

void initialize(float *array, int rows, int cols, float value);
void updateImage(float *pold, float *pnew, float *plim, int rows, int cols);

int main() {
    int M, N;  // Dimensiunile imaginii
    float *data, *pold, *pnew, *plim;

    // Citirea dimensiunilor imaginii din fișier
    pgm_size("image_640x480.pgm", &N, &M);

    // Alocare memorie pentru matricea inițială și cele intermediare, incluzând halo
    int extendedRows = M + 2;
    int extendedCols = N + 2;

    pold = (float *)malloc(extendedRows * extendedCols * sizeof(float));
    pnew = (float *)malloc(extendedRows * extendedCols * sizeof(float));
    plim = (float *)malloc(extendedRows * extendedCols * sizeof(float));
    data = (float *)malloc(M * N * sizeof(float));  // fără halo

    if (!pold || !pnew || !plim || !data) {
        printf("Eroare la alocarea memoriei.\n");
        return -1;
    }

    // Inițializare matrici cu 255, inclusiv halo
    initialize(pold, extendedRows, extendedCols, 255.0f);
    initialize(pnew, extendedRows, extendedCols, 255.0f);
    initialize(plim, extendedRows, extendedCols, 255.0f);

    // Citirea datelor imaginii în 'data'
    pgm_read("image_640x480.pgm", data, N, M);

    // Copiem 'data' în 'plim' respectând indexarea corectă și păstrând halo-ul la 255
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            plim[(i+1) * extendedCols + (j+1)] = data[i * N + j];
        }
    }

    // Aplică algoritmul de actualizare a imaginii
    updateImage(pold, pnew, plim, extendedRows, extendedCols);

    // Scrierea imaginii procesate înapoi în fișierul PGM
    // Copiem 'pnew' înapoi în 'data' pentru a scrie
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            data[i * N + j] = pnew[(i+1) * extendedCols + (j+1)];
        }
    }

    pgm_write("processed_image.pgm", data, N, M);

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
    for (int i = 1; i < rows - 1; i++) {
        for (int j = 1; j < cols - 1; j++) {
            pnew[i * cols + j] = 0.25 * (pold[(i - 1) * cols + j] + pold[(i + 1) * cols + j] +
                                         pold[i * cols + (j - 1)] + pold[i * cols + (j + 1)] - plim[i * cols + j]);
        }
    }
}
