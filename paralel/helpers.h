#ifndef HELPERS_H
#define HELPERS_H

#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "pgm_IO.h"

// initializeaza mpi si aloca resursele necesare
void initialize_mpi(int argc, char *argv[], int *rank, int *nproc, int *M, int *N, float **pold, float **pnew, float **plim, float **data, int *niter);

// distribuie datele de la procesul master catre toate procesele
void distribute_data(int rank, int nproc, int M, int N, float **masterdata, float *data, float *plim);

// proceseaza datele, aplica algoritmul si realizeaza schimbul de halo
void process_data(int rank, int nproc, int MP, int NP, int niter, float *pold, float *pnew, float *plim);

// finalizeaza prelucrarea datelor, aduna toate bucatile procesate si scrie imaginea finala
void finalize_and_write_data(int rank, int nproc, int M, int N, int MP, int NP, int niter, float *pnew, float *data, float *masterdata);

// copiaza datele fara halo
void copyWithoutHalo(float *src, float *dest, int rows, int cols);

// actualizeaza imaginea folosind algoritmul de reconstructie
void updateImage(float *pold, float *pnew, float *plim, int rows, int cols);

#endif // HELPERS_H
