#include "../convert.F90"

subroutine init_environ(rank, size, starttime)

  USE CONSTANTS
  IMPLICIT NONE
       
#ifdef USE_MPI
  include 'mpif.h'
#endif /* USE_MPI */

  Real starttime
  integer rank, size
  integer ierr

  call cpu_time(starttime)

#ifdef USE_MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
#else
  rank = 0
  size = 1
#endif /* USE_MPI */

  if (rank .eq. 0) then
#if defined (FLEX)
     print *, 'Code compiled for the fluctuation exchange approximation.'
#elif defined (THIRD_ORDER)
     print *, 'Code compiled for third order perturbation theory.'
#elif defined (SECOND_ORDER)
     print *, 'Code compiled for second order perturbation theory.'
#else
     print *, 'Code compiled for first order perturbation theory.'
#endif /* defined (FEA) */

     print *
     print *, 'Compiled for', np, ' processes.'
     print *, 'Running on', size, ' processes.'

  endif

  if (np .eq. size)  then

     if (rank == 0) then 
        print *, 'Compiled and run-time size are the same.'
        print *, 'Very good.'
        print *
     endif

  else 

     if (rank .eq. 0) then
        print *, 'np not equal size.  Stop.' 
     endif

#ifdef USE_MPI
     call MPI_Finalize(ierr)
#endif /* USE_MPI */

     stop

  endif

  if (rank .eq. 0) then
#ifdef SINGLE_PREC
     write(6,*) 'Compiled for SINGLE PRECISION.'
#else
     write(6,*) 'Compiled for DOUBLE PRECISION.'
#endif /* SINGLE_PREC */

  endif

  return
end subroutine init_environ
	 
