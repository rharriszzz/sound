/*
gcc -g -O2 -std=c99 -Wall -DMAIN -o fft.exe fft.c -lm
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <ctype.h>
#include <time.h>

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

#include "fft.h"

complex DOUBLE *copy_complex_double_array(complex DOUBLE *array, int n)
{
  int size = n * sizeof(complex DOUBLE);
  complex DOUBLE *result = (complex DOUBLE *)malloc(size);
  memcpy(result, array, size);
  return result;
}

void dft(complex DOUBLE *x, complex DOUBLE *X, int N, int inverse, int balanced)
{
  DOUBLE f = balanced ? SQRT(1.0/N) : inverse ? (1.0/N) : 1.0;
  complex DOUBLE W = CEXP((inverse ? 1.0 : -1.0) * I * 2 * M_PI / N);
  for (int k = 0; k < N; k++) {
    complex DOUBLE sum = 0;
    for (int n = 0; n < N; n++) {
      sum += x[n] * CPOW(W, n * k);
    }
    X[k] = f * sum;
  }
}

static void _fft(complex DOUBLE *buf, complex DOUBLE *out, int n, int step, int inverse)
{
  if (step < n) {
    _fft(out,        buf,        n, step * 2, inverse);
    _fft(out + step, buf + step, n, step * 2, inverse);
 
    for (int i = 0; i < n; i += 2 * step) {
      complex DOUBLE t = CEXP((inverse ? 1.0 : -1.0) * I * M_PI * i / n) * out[i + step];
      buf[(i + 0)/2] = out[i] + t;
      buf[(i + n)/2] = out[i] - t;
    }
  }
}

void fft(complex DOUBLE *buf, int n, int inverse, int balanced)
{
  int m = 1;
  while (m < n) {m *= 2;}
  if (m != n) {
    printf("Error: this fft routine requires that n be a power of two.\n");
    return;
  }
  complex DOUBLE *out = (complex DOUBLE *)calloc(n, sizeof(complex DOUBLE));
  for (int i = 0; i < n; i++) out[i] = buf[i];
  _fft(buf, out, n, 1, inverse);
  free(out);
  DOUBLE f = balanced ? (1.0/SQRT(n)) : inverse ? (1.0/n) : 1.0;
  for (int i=0; i < n; i++) {
    buf[i] *= f;
  }
}


void chirp_z(complex DOUBLE *x, int N,
             complex DOUBLE A, complex DOUBLE W,
             complex DOUBLE *X, int M)
{
  int Lmin = N + M - 1;
  int L = 1;
  int Lexp = 0;
  while (L < Lmin) {
    L *= 2;
    Lexp += 1;
  }
  if (0)
    printf("N=%d, M=%d, Lmin=%d, L=%d, A amp=%e, phase=%e, W amp=%e, phase=%e\n",
           N, M, Lmin, L, cabs(A), carg(A)/(2*M_PI), cabs(W), carg(W)/(2*M_PI)); fflush(stdout);
  complex DOUBLE *y = calloc(L, sizeof(complex DOUBLE));
  complex DOUBLE *v = calloc(L, sizeof(complex DOUBLE));
  complex DOUBLE *g = calloc(L, sizeof(complex DOUBLE));
  for (int n = 0; n < L; n++) {
    y[n] = (n < N) ? (CPOW(A, -n) * CPOW(W, n*n/2.0) * x[n]) : 0.0;
  }
  for (int n = 0; n < L; n++) {
    int d = (n < M) ? n : (n > (L - N)) ? (L - n) : 0;
    v[n] = CPOW(W, -d*d/2.0);
  }
  fft(y, L, FALSE, FALSE);
  fft(v, L, FALSE, FALSE);
  for (int n = 0; n < L; n++) {
    g[n] = y[n] * v[n];
  }
  fft(g, L, TRUE, FALSE);
   for (int k = 0; k < M; k++) {
    X[k] = g[k] * CPOW(W, k*k/2.0);
  }
  free(y); 
  free(v); 
  free(g); 
 }

void chirp_z_fft(complex DOUBLE *x, int n)
{
  chirp_z(x, n, (complex DOUBLE)1.0, CEXP(-I * 2 * M_PI / n), x, n);
}

#if MAIN
void show_complex_values(complex DOUBLE *X, int N)
{
  for (int i = 0; i < N; i++) {
    DOUBLE magnitude = CABS(X[i]);
    DOUBLE phase = CARG(X[i]) / (2 * M_PI);
    printf("%4d " F10_6 " " F10_6 "\n", i, magnitude, phase);
  }

}

int main(int argc, char **argv)
{
  int test = (argc > 1) ? strtol(argv[1], NULL, 10) : 0;
  int M = 9;
  int N = lrint(exp2((double)M));
  printf("data\n");
  complex DOUBLE *a = (complex DOUBLE *)calloc(N, sizeof(complex DOUBLE));
  for (int i = 0; i < N; i++) {
    if (test == 0) {
      DOUBLE c = 0;
      for (int j = 2; j < 12; j += 2) {
        c += SIN(2*M_PI*j*i/N) / j;
      }
      a[i] = c;
    } else if (test == 1) {
      DOUBLE c = 0;
      for (int j = 1; j < 12; j += 2) {
        c += SIN(2*M_PI*j*i/N) / (j + 1);
      }
      a[i] = c;
    } else if (test == 2) {
      a[i] = SIN(2*M_PI*16*i/N);
    } else if (test == 3) {
      a[i] = EXP(-i*4.0/N);
    } else if (test == 4) {
      a[i] = SIN(2*M_PI*16*i/N) * EXP(-i*4.0/N);
    }
  }
  show_complex_values(a, N);
  printf("\n" "dft\n");
  complex DOUBLE *X = (complex DOUBLE *)calloc(N, sizeof(complex DOUBLE));
  dft(a, X, N, FALSE, FALSE);
  show_complex_values(X, N);
  printf("\n" "inverse dft\n");
  complex DOUBLE *x = (complex DOUBLE *)calloc(N, sizeof(complex DOUBLE));
  dft(X, x, N, TRUE, FALSE);
  show_complex_values(x, N);
  printf("\n" "fft\n");
  complex DOUBLE *A = copy_complex_double_array(a, N);
  fft(A, N, FALSE, FALSE);
  show_complex_values(A, N);
  printf("\n" "fft2\n");
  A = copy_complex_double_array(a, N);
  chirp_z_fft(A, N);
  show_complex_values(A, N);
}
#endif
