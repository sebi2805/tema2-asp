#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "pgm_IO.h"  // Include the header file for PGM file I/O functions

// Function prototypes
void initialize(float *array, int rows, int cols, float value);
void updateImage(float *pold, float *pnew, float *plim, int rows, int cols);
void copyWithoutHalo(float *src, float *dest, int rows, int cols);

int main(int argc, char *argv[]) {
    int M, N;  // Dimensions of the image
    float *data, *pold, *pnew, *plim; // Pointers for image data buffers
    int niter = 1;  // Default number of iterations
    char *filename = "image_640x480.pgm";  // Default filename
    clock_t start, end;  // Variables to measure execution time
    double cpu_time_used;

    start = clock();

    // Process command-line arguments
    if (argc > 1) {
        niter = atoi(argv[1]);  // First argument is number of iterations
        if (argc > 2) {
            filename = argv[2];  // Second argument is the filename
        }
    }

    // Read the image dimensions from the file
    pgm_size(filename, &M, &N);

    // Allocate memory for the initial and intermediate matrices, including halo
    int extendedRows = M + 2;
    int extendedCols = N + 2;

    pold = (float *)malloc(extendedRows * extendedCols * sizeof(float));
    pnew = (float *)malloc(extendedRows * extendedCols * sizeof(float));
    plim = (float *)malloc(extendedRows * extendedCols * sizeof(float));
    data = (float *)malloc(M * N * sizeof(float));

    if (!pold || !pnew || !plim || !data) {
        printf("Error allocating memory.\n");
        return -1;
    }

    // Initialize matrices to 255, including halo
    initialize(pold, extendedRows, extendedCols, 255.0f);
    initialize(pnew, extendedRows, extendedCols, 255.0f);
    initialize(plim, extendedRows, extendedCols, 255.0f);

    // Read the original gray-scale image
    pgm_read(filename, data, M, N);

    // Map original image to plim buffer without halo
    for (int i = 1; i <= M; i++) {
        for (int j = 1; j <= N; j++) {
            plim[i * extendedCols + j] = data[(i-1) * N + j-1];
        }
    }

    // Image processing iterations
    for (int iter = 0; iter < niter; iter++) {
        updateImage(pold, pnew, plim, extendedRows, extendedCols);
        copyWithoutHalo(pnew, pold, extendedRows, extendedCols);
    }

    // Copy processed data back to original data buffer
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            data[i * N + j] = pnew[(i+1) * extendedCols + (j+1)];
        }
    }
    pgm_write("processed_image.pgm", data, M, N);

    // Free memory
    free(pold);
    free(pnew);
    free(plim);
    free(data);

    end = clock();  
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Execution time is: %f seconds.\n", cpu_time_used);
    return 0;
}

// Initialize the array with the specified value
void initialize(float *array, int rows, int cols, float value) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            array[i * cols + j] = value;
        }
    }
}

// Update image based on old and limit values
void updateImage(float *pold, float *pnew, float *plim, int rows, int cols) {
    for (int i = 1; i < rows ; i++) {
        for (int j = 1; j < cols ; j++) {
            pnew[i * cols + j] = 0.25 * (pold[(i - 1) * cols + j] + pold[(i + 1) * cols + j] +
                                         pold[i * cols + (j - 1)] + pold[i * cols + (j + 1)] - plim[i * cols + j]);
        }
    }
}

// Copy data without including the halo region
void copyWithoutHalo(float *src, float *dest, int rows, int cols) {
    for (int i = 1; i < rows ; i++) {
      memcpy(dest + i * cols + 1, src + i * cols + 1, (cols - 2) * sizeof(float));
    }
}
