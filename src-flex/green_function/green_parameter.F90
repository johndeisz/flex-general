#include "../convert.F90"

subroutine green_parameter(rank, g, t, x, c_r, delta_g_r, delta_gp_r)

  USE CONSTANTS
  IMPLICIT NONE

#ifdef USE_MPI
  include 'mpif.h'
#endif /* USE_MPI */

  INTEGER rank
  COMPLEX g(0:4*nb-1,0:4*nb-1,0:mp1,0:nc1)
  REAL t, x(0:1,0:1)
  COMPLEX c_r(0:1,0:1,0:4*nb-1,0:4*nb-1,0:nc1)
  COMPLEX, dimension (0:4*nb-1, 0:4*nb-1, 0:nc1) :: delta_g_r, delta_gp_r

  REAL delta_tau

  INTEGER k, na1, na2

  COMPLEX gt(0:2)

  INTEGER ierr

  delta_tau = (1.0d0 / t) / float(m)

  do na1 = 0, 4*nb-1
     do na2 = 0, 4*nb-1

        do k = 0, nc1
           
           if (rank .eq. 0) then

              gt(0) = g(na1,na2,0,k)
              gt(1) = g(na1,na2,1,k)
              gt(2) = g(na1,na2,2,k)
              
           endif

#ifdef USE_MPI
           call MPI_Bcast( gt, 3, MPI_COMPLEX, 0, MPI_COMM_WORLD, ierr)
#endif /* USE_MPI */

           c_r(1,1,na1,na2,k) = &
                ( gt(0) - delta_g_r(na1,na2,k) / 2.0d0 + & 
                delta_gp_r(na1,na2,k) * & 
                tanh(0.5d0*x(1,0) / t) / (2.0d0 * x(1,0) ) ) / & 
                ( tanh(0.5d0 * x(1,0) / t) /  (2.0d0 * x(1,0) ) - & 
                tanh(0.5d0 * x(1,1) / t) /  (2.0d0 * x(1,1) ) )

           c_r(1,0,na1,na2,k) = delta_gp_r(na1,na2,k) - c_r(1,1,na1,na2,k)
            
           c_r(0,1,na1,na2,k) = & 
                -( ( 2.0d0*gt(1) - 0.5d0*gt(2)- 1.5d0*gt(0) ) / & 
                delta_tau  -  0.5d0 * delta_gp_r(na1,na2,k) + &
                0.5d0 * delta_g_r(na1,na2,k) * & 
                x(0,0) * tanh( 0.5d0*x(0,0) / t) ) / &
                ( 0.5d0 * x(0,0) * tanh(0.5d0 * x(0,0) / t)  - & 
                0.5d0 * x(0,1) * tanh(0.5 * x(0,1) / t)  )
                     
           c_r(0,0,na1,na2,k) = -delta_g_r(na1,na2,k) - c_r(0,1,na1,na2,k)

                
        enddo
     enddo
  enddo

  return
end subroutine green_parameter
