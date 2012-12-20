#include "../convert.F90"

subroutine tmat_param(rank, method, t_mat, t, y, d, delta_t, delta_t_prime)

  USE CONSTANTS

#ifdef USE_MPI
  include 'mpif.h'
#endif /* USE_MPI */

  INTEGER rank
  INTEGER method
  COMPLEX t_mat(0:16*nb*nb-1,0:16*nb*nb-1,0:mp1,0:nc1) 
  REAL t, y(0:1, 0:1) 
  COMPLEX d(0:1,0:1,0:16*nb*nb-1,0:16*nb*nb-1,0:nc1)
  COMPLEX, dimension (0:16*nb*nb-1,0:16*nb*nb-1,0:nc1) :: &
       delta_t, delta_t_prime

  REAL delta_tau
  COMPLEX temp_d(0:1,0:1)
  
  INTEGER k
  INTEGER i, j, i1, j1

  INTEGER ierr

   if (rank .eq. 0) then
     open(unit=9,file='my_error_file',status='old', access='append')
     write(9,*) '    starting tmat_param'
     close(unit=9)
  endif


  delta_tau = (1.0d0 / t) / float(m)

  if (rank .eq. 0) then
  do k = 0, nc1

     do i = 0, 16*nb*nb-1
        do j = 0, 16*nb*nb-1

!           if (rank .eq. 0) then

              if (method .eq. 1) then

                 temp_d(1,1) = ( t_mat(i,j,0,k) - delta_t(i,j,k) / 2.0d0 + & 
                      delta_t_prime(i,j,k) / & 
                      ( tanh(0.5d0*y(1,0)/t) * 2.0d0 * y(1,0) ) ) / &
                      ( 1.0d0 / ( tanh(0.5d0 * y(1,0) / t) * & 
                      2.0d0 * y(1,0) ) -  1.0d0 / & 
                      ( tanh(0.5d0 * y(1,1) / t) * 2.0d0 * y(1,1) ) ) 

              else 

                 temp_d(1,1) = cmplx(0.0d0,0.0d0)

              endif
		  
              temp_d(1,0) = delta_t_prime(i,j,k) - temp_d(1,1)
		
              if (method .eq. 1) then

                 temp_d(0,1) = -( ( 2.0d0 * t_mat(i,j,1,k) - & 
                      0.5d0 * t_mat(i,j,2,k) - & 
                      1.5d0 * t_mat(i,j,0,k) ) / delta_tau - &
                      0.5d0 * delta_t_prime(i,j,k) + &
                      0.5d0 * delta_t(i,j,k) * y(0,0) / &
                      tanh( 0.5d0 * y(0,0) / t) ) / &
                      ( 0.5d0 * y(0,0) / tanh(0.5d0 * y(0,0) / t) - & 
                      0.5d0 * y(0,1) / tanh(0.5d0 * y(0,1) / t)  )
	  
              else 

                 temp_d(0,1) = cmplx(0.0d0,0.0d0)

              endif

              temp_d(0,0) = -delta_t(i,j,k) - temp_d(0,1)

!           endif

! #ifdef USE_MPI
!            call MPI_Bcast( temp_d, 4, MPI_COMPLEX, 0, MPI_COMM_WORLD, ierr)
! #endif	/* USE_MPI */
           
           do i1 = 0, 1
              do j1 = 0, 1

                 d(i1,j1,i,j,k) = temp_d(i1,j1) 

              enddo
           enddo

        enddo
     enddo

  enddo
 endif

#ifdef USE_MPI
  call MPI_Bcast( d, 4*16*nb*nb*16*nb*nb*nc, MPI_COMPLEX, 0, &
    MPI_COMM_WORLD, ierr)
#endif

   if (rank .eq. 0) then
     open(unit=9,file='my_error_file',status='old', access='append')
     write(9,*) '    leaving tmat_param'
     close(unit=9)
  endif

  return
end subroutine tmat_param
