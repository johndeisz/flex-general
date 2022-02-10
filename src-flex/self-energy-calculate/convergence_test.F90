#include "../convert.F90"

subroutine convergence_test(convergence, rank, iteration, &
     sigma, sigma_old, convg_crit, last_it_time)

  USE CONSTANTS
  IMPLICIT NONE

#ifdef USE_MPI
  include "mpif.h"
#endif /* USE_MPI */

  LOGICAL convergence
  INTEGER rank, iteration 
  COMPLEX, dimension (0:4*nb-1,0:4*nb-1,0:mp1,0:nc1) :: sigma, sigma_old
  REAL convg_crit, last_it_time

  INTEGER i, j, l, k

  INTEGER tmp, converged_count, points
  REAL current_time

#ifdef USE_MPI
  INTEGER stat(MPI_STATUS_SIZE)
  INTEGER ierr
#endif /* USE_MPI */

  points = 4 * 4 * nb * nb * nc * m
  converged_count = 0

  do i = 0, 4*nb-1
     do j = 0, 4*nb-1

        do l = 0, mp1
           do k = 0, nc1

              if ( ( cabs( sigma(i,j,l,k) - sigma_old(i,j,l,k) ) .lt. &
                   convg_crit * cabs( sigma_old(i,j,l,k) ) ) .or. &
                   ( cabs( sigma(i,j,l,k) ) .lt. 1.0e-07 ) ) then
                
                 converged_count = converged_count + 1
                
              endif

           enddo
        enddo
     enddo
  enddo

#ifdef USE_MPI
  if (rank .eq. 0) then

     do i = 1, np - 1
        tmp = 0
        call MPI_Send(tmp, 1, MPI_INTEGER, i, 0, MPI_COMM_WORLD, ierr)
        call MPI_Recv(tmp, 1, MPI_INTEGER, i, 1, MPI_COMM_WORLD, stat, ierr)
        converged_count = converged_count + tmp
     enddo

  else 

     call MPI_Recv(tmp, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, stat, ierr)
     call MPI_Send(converged_count, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, ierr)

  endif

#endif /* USE_MPI */

  if (rank .eq. 0) then

     write(6,fmt=600,advance='no') 100.0d0*dfloat(converged_count)/dfloat(points)

     if (converged_count .eq. points) then 
        convergence = .true.
     else 
        convergence = .false.
     endif

  endif

#ifdef USE_MPI
  call MPI_Bcast(convergence, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif /* USE_MPI */
  
600 format(f10.6)

  return
end subroutine convergence_test
