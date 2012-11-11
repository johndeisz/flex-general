#include "../convert.F90"

subroutine init_sigma1( uu, density, sigma1)

  USE CONSTANTS
  IMPLICIT NONE

  REAL uu, density
  COMPLEX sigma1(0:4*nb-1,0:4*nb-1)
  INTEGER i

  sigma1 = 0.0d0

  do i = 0, nb-1

     sigma1( i*4+0, i*4+0) =  uu * density / float(2*nb)
     sigma1( i*4+1, i*4+1) =  uu * density / float(2*nb)
     sigma1( i*4+2, i*4+2) = -uu * density / float(2*nb)
     sigma1( i*4+3, i*4+3) = -uu * density / float(2*nb)

  enddo

  return
end subroutine init_sigma1
