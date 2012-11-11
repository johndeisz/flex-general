#include "../convert.F90"

subroutine generate_tau_eps_omega(t, tau, epsilon, omega)

  USE CONSTANTS
  IMPLICIT NONE

#ifdef USE_MPI
  include 'mpif.h'
#endif /* USE_MPI */

  REAL t
  REAL tau(0:mp1), epsilon(0:mp1), omega(0:mp1)

  INTEGER rank, ierr

  REAL delta_tau
  REAL pi

  INTEGER l

#ifdef USE_MPI
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
#else
  rank = 0
#endif /* USE_MPI */

  delta_tau = (1.d0 / t ) / float(m)
  pi = acos(-1.0d0)

  do l = 0, mp1
     tau(l) = float(l + rank*mp) * delta_tau  
  enddo

  if (np .eq. 1) then

     do l = 0, m / 2 - 1
        epsilon(l) = float(2 * l + 1) * pi * t
        omega(l) = float(2 * l) * pi * t
     enddo
     do l = m/2, m1
        epsilon(l) = float(2 * (l-m) + 1) * pi * t
        omega(l) = float(2 * (l-m)) * pi * t
     enddo

  else 
     
     if ( rank .lt.  np/2 ) then

        do l = 0, mp1
           epsilon(l) = float(2 * (l+rank*mp) + 1) * pi * t
           omega(l) = float(2 * (l+rank*mp) ) * pi * t
        enddo
        
     else

        do l = 0, mp1
           epsilon(l) = float(2 * (l+rank*mp - m) + 1) * pi * t
           omega(l) = float(2 * (l+rank*mp - m) ) * pi * t           
        enddo
        
     endif
  endif

  return
end subroutine generate_tau_eps_omega
