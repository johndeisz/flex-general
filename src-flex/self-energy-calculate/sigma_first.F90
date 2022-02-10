#include "../convert.F90"

subroutine sigma_first(sigma1, sigma1_old, delta_sigma1, &
     delta_sigma1_old, alpha, g_tau0_local, gamma0_ph, sigma_tol, &
     sigma_converged, alpha_scheme, iteration)

  USE CONSTANTS
  IMPLICIT NONE

#ifdef USE_MPI
  include 'mpif.h'
#endif /* USE_MPI */


  COMPLEX, dimension (0:4*nb-1, 0:4*nb-1) :: sigma1, &
       sigma1_old, delta_sigma1, delta_sigma1_old, g_tau0_local

  REAL alpha
  COMPLEX gamma0_ph(0:16*nb*nb-1, 0:16*nb*nb-1)
  REAL sigma_tol
  LOGICAL sigma_converged
  INTEGER alpha_scheme
  INTEGER iteration

  INTEGER rank
  INTEGER na, nap, na0, na1

  COMPLEX gamma0

  INTEGER points

#ifdef USE_MPI
  INTEGER ierr
#endif
      
#ifdef USE_MPI
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
#else
  rank = 0
#endif /* USE_MPI */

  if (rank .eq. 0) then

     sigma1_old = sigma1
     if (iteration .gt. 1) then
        delta_sigma1_old = delta_sigma1
     endif
      
     sigma1 = cmplx(0.0d0, 0.0d0)

     do na = 0, 4*nb-1
        do nap = 0, 4*nb-1

           do na0 = 0, 4*nb-1
              do na1 = 0, 4*nb-1

                 gamma0 = gamma0_ph(4*nb*na0+na1,4*nb*nap+na)

                 sigma1(na,nap) = sigma1(na,nap) + &
                      0.5d0 * gamma0 * g_tau0_local(na1,na0)

                 if (na0 .eq. na1) then
                  
                    if (mod(na0,4) .le. 1) then
                     
                       sigma1(na,nap) = sigma1(na,nap) + &
                            0.5d0 * gamma0 
                    
                    endif
                  
                 endif
                 
              enddo
           enddo

        enddo
     enddo

     points = 0

     do na = 0, 4*nb-1
        do nap = 0, 4*nb-1

           if (cabs(sigma1(na,nap) - sigma1_old(na,nap)) .lt. sigma_tol) then

              points = points + 1

           endif

        enddo
     enddo

     if (points .lt. 16*nb*nb) then
        sigma_converged = .false.
     endif

     write(6,fmt=400,advance='NO') 100.0d0*dfloat(points)/dfloat(16*nb*nb)
        
     delta_sigma1 = sigma1 - sigma1_old
     sigma1 = alpha*sigma1 + (1.0d0 - alpha)*sigma1_old

  endif

#ifdef USE_MPI
  call MPI_Bcast(alpha, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(sigma1, 16*nb*nb, MPI_COMPLEX, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(sigma_converged, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif /* USE_MPI */

400 format('       ', f8.4)

  return
end subroutine sigma_first
