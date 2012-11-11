#include "../convert.F90"

subroutine psi_transform(psi, isignv)

  USE CONSTANTS
  IMPLICIT NONE

  COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)
  INTEGER isignv(0:3)
      
  COMPLEX temp_3D(0:nl-1)

  INTEGER i, j

  do i = 0, 2*nb-1
     do j = 0, 2*nb-1

        temp_3D = psi(i,j,:)

        call fft_3D_lattice(temp_3D, isignv)  
        
        psi(i,j,:) = temp_3D
          
     enddo
  enddo
  
  return
end subroutine psi_transform
