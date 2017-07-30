
#ifndef _FFT_
#define _FFT_ 1

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#if 0 /* defined(_HAVE_LONG_DOUBLE) && !defined(_LDBL_EQ_DBL) */
#define DOUBLE long double
#define F10_6 "%10.6Lf"
#define SQRT sqrtl
#define SIN  sinl
#define COS  cosl
#define EXP  expl
#define CEXP cexpl
#define CPOW cpowl
#define CABS cabsl
#define CARG cargl
#else
#define DOUBLE double
#define F10_6 "%10.6f"
#define SQRT sqrt
#define SIN  sin
#define COS  cos
#define EXP  exp
#define CEXP cexp
#define CPOW cpow
#define CABS cabs
#define CARG carg
#endif

complex DOUBLE *copy_complex_double_array(complex DOUBLE *array, int n);
void dft(complex DOUBLE *x, complex DOUBLE *X, int N, int inverse, int balanced);
void fft(complex DOUBLE buf[], int n, int inverse, int balanced);
void chirp_z(complex DOUBLE *x, int N,
             complex DOUBLE A, complex DOUBLE W,
             complex DOUBLE *X, int M);
void chirp_z_fft(complex DOUBLE *x, int n);
#endif
