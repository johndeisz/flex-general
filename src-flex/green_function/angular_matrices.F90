#include "../convert.F90"

subroutine angular_matrices(Lvec)

  USE CONSTANTS
  IMPLICIT NONE
  COMPLEX Lvec(1:3,0:nb-1,0:nb-1)

  COMPLEX Lvec_ml(1:3,-2:2,-2:2)
  COMPLEX sum(1:3)
  REAL L_plus(-2:2,-2:2)
  REAL L_minus(-2:2,-2:2)

  COMPLEX psi(-2:2, 0:10) ! Fix this for multiband systems

  integer ml, nu1, nu2, ml1, ml2
  
  psi = cmplx(0.0d0, 0.0d0)

  ! Define spatial wave functions
  psi(-2,0) = cmplx(0.0d0, 0.0d0)
  psi(-1,0) = 1.0d0 / sqrt(2.0d0) 
  psi(0,0) = 0.0d0
  psi(1,0) = -1.0d0 / sqrt(2.0d0)
  psi(2,0) = 0.0d0

  psi(-2,1) = 0.0d0 
  psi(-1,1) = cmplx(0.0d0,1.0d0 / sqrt(2.0d0))
  psi(0,1) = 0.0d0
  psi(1,1) = cmplx(0.0d0, 1.0d0 / sqrt(2.0d0))
  psi(2,1) = 0.0d0

  psi(-2,2) = cmplx(0.0d0, 1.0d0 / sqrt(2.0d0))
  psi(-1,2) = 0.0d0
  psi(0,2) = 0.d0
  psi(1,2) = 0.d0
  psi(2,2) = cmplx( 0.0d0, -1.0d0 / sqrt(2.0d0))

  Lvec_ml = cmplx(0.0d0, 0.d0)
  L_plus = 0.0d0
  L_minus = 0.0d0

  Lvec = cmplx(0.0d0,0.0d0)

  if (nb .eq. 3) then

     do ml = -2, 2
        Lvec_ml(3,ml,ml) = float(ml)
     enddo

     do ml = -2,1
        L_plus(ml+1,ml) = sqrt(6.0d0 - float(ml*(ml+1)))
     enddo

     do ml = -1, 2
        L_minus(ml-1,ml) = sqrt(6.0d0 - float(ml*(ml-1)))
     enddo

     Lvec_ml(1,:,:) = 0.5d0*cmplx(L_plus+L_minus,0.0d0)
     Lvec_ml(2,:,:) = 0.5d0*cmplx(0.0d0,L_minus-L_plus)

     ! Compute the spin-orbit matrix elements
	
     do nu1=0,nb-1
        do nu2 = 0,nb-1

           sum = cmplx(0.0d0, 0.0d0)

           do ml1 = -2,2
              do ml2 = -2,2

                 sum = sum + conjg(psi(ml1,nu1))* &
                      Lvec_ml(:,ml1,ml2)*psi(ml2,nu2)
              
              enddo
           enddo
            
           Lvec(:,nu1,nu2) = sum

        enddo
     enddo
 
  endif
     
  return
end subroutine angular_matrices
