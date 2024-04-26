# Image Reconstruction from Edge Detection

This project focuses on reconstructing images that have been transformed to highlight edges. This transformation allows for more efficient storage using compression techniques such as Huffman coding. The project includes both sequential and parallel implementations, with the parallel version utilizing MPI for distributed processing.

## Project Overview

The application processes edge-detected images to reconstruct the original images. It uses the following formula to calculate the edge values, which helps in distinguishing between edge and non-edge pixels:
plim[i][j] = img[i−1][j] + img[i+1][j] + img[i][j−1] + img[i][j+1] - 4 \* img[i][j]

If a pixel's value is significantly different from its neighbors, it is likely part of an edge. The project provides a sequential version for simplicity and a parallel version using MPI to demonstrate performance improvements in a distributed environment.

## File Structure

- `main.c`: Main program file for both sequential and parallel implementations.
- `pgm_IO.h`/`pgm_IO.c`: Functions for reading and writing PGM files.
- `helpers.h`/`helpers.c`: Helper functions used in the parallel version.
- `Makefile`: Instructions to build both the sequential and parallel versions.
- `image_640x480.pgm`: Sample PGM file used as input.
- `parallel-allotted-matrix.png`: Diagram showing the data partitioning in the MPI implementation.

## Prerequisites

To compile and run this project, you need:

- A C compiler like GCC.
- Make.
- An MPI implementation (e.g., MPICH or OpenMPI) for the parallel version.
- Proper setup of all dependencies for handling PGM files.

## Compilation

To compile the project, use the provided Makefile for each folder.

## Execution

It is recommended to have a Linux distribution for this, you will need the mpi and imagegick.
For the sequential program you just run:

- ./reconstruct [number of iterations]
- convert image_640x480.pgm imagine.png (so you can view the image in a png format)

For the parallel program you just run:

- mpirun -np [number of processes it is recommended 4 | at 26 April the program still does not take in account some cases of matrix divizibility] ./program [number of iterations]
- convert image_640x480.pgm imagine.png (so you can view the image in a png format)

## Output

The output of the program is a PGM file named processed_image.pgm, which contains the reconstructed image from its edge-detected format.

# MPI Parallel Implementation

In the parallel version, MPI is used to distribute the image processing workload across multiple processors. The strategy for data partitioning is illustrated in parallel-allotted-matrix.png, showing how the image matrix is divided among processors to optimize parallel processing and minimize communication overhead.
