#ifdef DOUBLE_PREC
#define REAL double precision
#define MPI_REAL MPI_DOUBLE_PRECISION
#define COMPLEX double complex
#define MPI_COMPLEX MPI_DOUBLE_COMPLEX
#define float dfloat
#define acos dacos
#define cos dcos
#define cmplx dcmplx
#define cabs cdabs
#define sqrt dsqrt
#define conjg dconjg
#define real dreal
#define imag dimag
#define exp dexp
#define cexp cdexp
#define csqrt cdsqrt
#define clog cdlog
#define cgemm zgemm
#define cgetrf zgetrf
#define cgetri zgetri
#define cgetrs zgetrs
#define cgehrd zgehrd
#define chseqr zhseqr
#define chegv zhegv
#define cgesv zgesv
#define cgeev zgeev
#define sfftw_plan_dft_1d dfftw_plan_dft_1d
#define sfftw_execute dfftw_execute
#define sfftw_destroy_plan dfftw_destroy_plan
#endif /* DOUBLE PRECISION */
