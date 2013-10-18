#include "../convert.F90"

subroutine gamma0_define(rank, gamma0_ph, uu, up, uj)

  USE CONSTANTS
  IMPLICIT NONE

  INTEGER rank
  COMPLEX gamma0_ph(0:16*nb*nb-1, 0:16*nb*nb-1)
  REAL uu, up, uj

  INTEGER na1, na2, na3, na4
  INTEGER nu1, nu2, nu3, nu4
  INTEGER is1, is2, is3, is4
  INTEGER ind1, ind2, ind3, ind4
  INTEGER ind1p, ind2p, ind3p, ind4p

  COMPLEX gamma0(0:4*nb-1, 0:4*nb-1, 0:4*nb-1, 0:4*nb-1)

  COMPLEX, dimension (1:3, 0:1, 0:1) :: pauli
  REAL, dimension (0:1, 0:1) :: dirac
  REAL, dimension (0:nb-1) :: uum
  REAL, dimension (0:nb-1, 0:nb-1) :: upm, ujm

  REAL s0, c0

  !     First, initialize the gamma array
  gamma0 = cmplx(0.0d0, 0.0d0)

  dirac = 0.0d0
  do is1 = 0, 1
     dirac(is1,is1) = 1.0d0
  enddo
  
  pauli(1,0,0) = 0.0d0
  pauli(1,1,0) = 1.0d0
  pauli(1,0,1) = 1.0d0
  pauli(1,1,1) = 0.0d0

  pauli(2,0,0) = 0.0d0
  pauli(2,1,0) = cmplx(0.0d0,1.0d0)
  pauli(2,0,1) = cmplx(0.0d0,-1.0d0)
  pauli(2,1,1) = 0.0d0

  pauli(3,0,0) = 1.0d0
  pauli(3,1,0) = 0.0d0
  pauli(3,0,1) = 0.0d0
  pauli(3,1,1) = -1.d0

!  if (rank .eq. 0) then
!    write(6,*) 'Cubic symmetrized vertex used'
!  endif
!  do nu1 = 0, nb-1
!    uum(nu1) = uu
!  enddo

!  do nu1 = 0, nb-1
!   do nu2= 0, nb-1
!     upm(nu1,nu2) = up
!     ujm(nu1,nu2) = uj
!   enddo
!  enddo

   if (rank .eq. 0) then
     write(6,*) 'SrRuO Band-dependent vertex used'
   endif
   if (nb .ne. 3) then
   if (rank .eq. 0) then
     write(6,*) 'Nb does not equal 3 - stopping'
   endif
  endif
  uum(0) = 0.969d0*uu
  uum(1) = 0.969d0*uu
  uum(2) = 1.063d0*uu

  upm(0,1) = 0.974d0*up
  upm(1,0) = upm(0,1)
  upm(0,2) = 1.015d0*up
  upm(2,0) = upm(0,2)
  upm(1,2) = upm(0,2)
  upm(2,1) = upm(0,2)

  ujm(0,1) = 0.923d0*uj
  ujm(1,0) = ujm(0,1)
  ujm(0,2) = 1.000d0*uj
  ujm(2,0) = ujm(0,2)
  ujm(1,2) = ujm(0,2)
  ujm(2,1) = ujm(0,2)

  do nu1 = 0, nb-1
     do nu2 = 0, nb-1
        do nu3 = 0, nb-1
           do nu4 = 0, nb-1

              s0 = 0.0d0
              c0 = 0.0d0

              if (nu1 .eq. nu2) then
                 if (nu3 .eq. nu2) then
                    if (nu4 .eq. nu3) then

                       s0 = uum(nu1)
                       c0 = uum(nu1)

                    endif
                 else
                    if (nu4 .eq. nu3) then

                       s0 = ujm(nu1,nu3)
                       c0 = ujm(nu1,nu3)

                    endif
                 endif

              else

                 if (nu1 .eq. nu3) then
                    if (nu2 .eq. nu4) then

                       s0 = ujm(nu1,nu2)
                       c0 = -ujm(nu1,nu2) + 2.0d0*upm(nu1,nu2)
                    
                    endif
                 else
                    if (nu1 .eq. nu4) then
                       if (nu2 .eq. nu3) then

                          s0 = upm(nu1,nu2)
                          c0 = 2.0*ujm(nu1,nu2) - upm(nu1,nu2)

                       endif
                    endif
                 endif
              endif
                    
              do is1 = 0,1
                 ind1 = 4*nu1 + is1

                 do is2 = 0,1
                    ind2 = 4*nu2 + is2

                    do is3 = 0,1
                       ind3 = 4*nu3 + is3

                       do is4 = 0,1
                          ind4 = 4*nu4 + is4

                          gamma0(ind1,ind2,ind3,ind4) = -0.5d0 * s0 * &
                               (pauli(1,is1,is3)*pauli(1,is2,is4) + &
                               pauli(2,is1,is3)*pauli(2,is2,is4) + &
                               pauli(3,is1,is3)*pauli(3,is2,is4)) + &
                               0.5d0 * c0 * dirac(is1,is3)*dirac(is2,is4)

                      
                       enddo
                    enddo
                 enddo
              enddo

           enddo
        enddo
     enddo
  enddo

  do nu1=0,nb-1
     do nu2=0, nb-1
        do nu3=0, nb-1
           do nu4=0,nb-1
              
              do is1 = 0,1
                 ind1 = 4*nu1 + is1
                 ind1p = ind1 + 2

                 do is2 = 0,1
                    ind2 = 4*nu2 + is2
                    ind2p = ind2 + 2

                    do is3 = 0,1
                       ind3 = 4*nu3 + is3
                       ind3p = ind3 + 2

                       do is4 = 0,1
                          ind4 = 4*nu4 + is4
                          ind4p = ind4 + 2

                          gamma0(ind1p,ind2,ind3p,ind4) = &
                               -gamma0(ind2,ind3,ind4,ind1)

                          gamma0(ind1,ind2p,ind3,ind4p) = &
                               -gamma0(ind1,ind4,ind3,ind2)

                          gamma0(ind1,ind2p,ind3p,ind4) = &
                               gamma0(ind1,ind3,ind4,ind2)

                          gamma0(ind1p,ind2,ind3,ind4p) = &
                               gamma0(ind2,ind4,ind3,ind1)

                          gamma0(ind1p,ind2p,ind3p,ind4p) = &
                               gamma0(ind4,ind3,ind2,ind1)

                       enddo
                    enddo
                 enddo
              enddo

           enddo
        enddo
     enddo
  enddo

  do na1 = 0, 4*nb-1
     do na2 = 0, 4*nb-1
        do na3 = 0, 4*nb-1
           do na4 = 0, 4*nb-1

              gamma0_ph( 4*nb*na1+na3, 4*nb*na4+na2) = gamma0(na1,na2,na3,na4)

           enddo
        enddo
     enddo
  enddo

  return
end subroutine gamma0_define
