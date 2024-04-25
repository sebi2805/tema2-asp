#include "pgm_IO.h"
#include <stdlib.h>
#include <stdio.h>

int main() {
    const char *input_filename = "image_640x480.pgm";
    const char *output_filename = "output.pgm";
    int iterations[] = {800, 1000, 1200, 1400, 1600};  // Array de iterații
    int num_iterations = sizeof(iterations) / sizeof(iterations[0]);

    int M, N;
    float *data, *pold, *pnew, *plim;

    pgm_size(input_filename, &M, &N);
    data = (float *) malloc(M * N * sizeof(float));
    pold = (float *) malloc((M + 2) * (N + 2) * sizeof(float));
    pnew = (float *) malloc((M + 2) * (N + 2) * sizeof(float));
    plim = (float *) malloc((M + 2) * (N + 2) * sizeof(float));

    pgm_read(input_filename, data, M, N);

    for (int i = 0; i < M + 2; i++) {
        for (int j = 0; j < N + 2; j++) {
            pold[i * (N + 2) + j] = 255;
            pnew[i * (N + 2) + j] = 255;
            plim[i * (N + 2) + j] = (i > 0 && i <= M && j > 0 && j <= N) ? data[(i - 1) * N + (j - 1)] : 255;
        }
    }

    for (int n = 0; n < num_iterations; n++) {
        int niter = iterations[n];
        for (int iter = 0; iter < niter; iter++) {
            for (int i = 1; i <= M; i++) {
                for (int j = 1; j <= N; j++) {
                    int index = i * (N + 2) + j;
                    pnew[index] = 0.25 * (pold[index - (N + 2)] + pold[index + (N + 2)] + pold[index - 1] + pold[index + 1] - plim[index]);
                }
            }
            for (int i = 1; i <= M; i++) {
                for (int j = 1; j <= N; j++) {
                    pold[i * (N + 2) + j] = pnew[i * (N + 2) + j];
                }
            }
        }
        printf("Iteratia %d finalizata.\n", niter);
    }

    pgm_write(output_filename, pold + (N + 2) + 1, M, N); // Sări peste halo

    free(data);
    free(pold);
    free(pnew);
    free(plim);

    return 0;
}
