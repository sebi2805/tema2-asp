#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void pgm_size (const char *filename, int *nx, int *ny);
void pgm_read (const char *filename, void *vx, int nx, int ny);
void pgm_write(const char *filename, void *vx, int nx, int ny);
