#include "../convert.F90"

subroutine calc_new_alpha(alpha, delta_sigma_e0_old, sigma, sigma_old)

  USE CONSTANTS
  IMPLICIT NONE

  REAL alpha
  COMPLEX delta_sigma_e0_old(0:4*nb-1,0:4*nb-1,0:nc1)
  COMPLEX, dimension (0:4*nb-1,0:4*nb-1,0:mp1,0:nc1) :: sigma, sigma_old

  REAL sum
  integer i, j, k

  sum = 0.0d0

  do i = 0, 4*nb-1
     do j = 0, 4*nb-1
        do k = 0, nc1

            sum = sum + real( sigma(i,j,0,k) - sigma_old(i,j,0,k) ) * & 
                 real( delta_sigma_e0_old(i,j,k) ) + &
                 imag( sigma(i,j,0,k) - sigma_old(i,j,0,k) ) * & 
                 imag( delta_sigma_e0_old(i,j,k) )

            delta_sigma_e0_old(i,j,k) = sigma(i,j,0,k) - sigma_old(i,j,0,k)

         enddo
      enddo
   enddo
   
   !     If the change in sigma at zero frequency is "mostly parallel"
   !     with the same from the previous iteration, then increase 
   !     the mixing.  If the change is "anti-parallel" in comparison
   !     to the previous iteration (i.e. oscillating), then lower the
   !     mixing.

   if (sum .gt. 0.0d0) then
      if (alpha .le. 1.3d0) then
         alpha = 1.1d0 * alpha
      endif
   else
      alpha = alpha / 1.5d0
   endif

   return
 end subroutine calc_new_alpha
