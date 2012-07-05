#include "../convert.F"

      subroutine gamma0_define(gamma0_ph, uu, up, uj)

#include "../constants.F"

c /***************************************************************
c * gamma_chi -- program to calculate the gamma and chi arrays  *
c *                                                             *
c * Author: Tom Slife                                           *
c ***************************************************************/

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

      REAL s0, c0

c     First, initialize the gamma array
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

      do nu1 = 0, nb-1
        do nu2 = 0, nb-1
          do nu3 = 0, nb-1
            do nu4 = 0, nb-1

              s0 = 0.0d0
              c0 = 0.0d0

              if (nu1 .eq. nu2) then
                if (nu3 .eq. nu2) then
                  if (nu4 .eq. nu3) then

                    s0 = uu
                    c0 = uu

                  endif
                else
                  if (nu4 .eq. nu3) then

                    s0 = uj
                    c0 = uj

                  endif
                endif

              else

                if (nu1 .eq. nu3) then
                  if (nu2 .eq. nu4) then

                    s0 = uj
                    c0 = -uj + 2.0d0*up
                    
                  endif
                else
                  if (nu1 .eq. nu4) then
                    if (nu2 .eq. nu3) then

                      s0 = up
                      c0 = 2.0*uj - up

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

                      gamma0(ind1,ind2,ind3,ind4) = 
     $                   -0.5d0 * s0 * 
     $                   (pauli(1,is1,is3)*pauli(1,is2,is4) +
     $                   pauli(2,is1,is3)*pauli(2,is2,is4) +
     $                   pauli(3,is1,is3)*pauli(3,is2,is4)) +
     $                   0.5d0 * c0 * dirac(is1,is3)*dirac(is2,is4)

                      
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

                      gamma0(ind1p,ind2,ind3p,ind4) = 
     $                   -gamma0(ind2,ind3,ind4,ind1)

                      gamma0(ind1,ind2p,ind3,ind4p) = 
     $                   -gamma0(ind1,ind4,ind3,ind2)

                      gamma0(ind1,ind2p,ind3p,ind4) =
     $                   gamma0(ind1,ind3,ind4,ind2)

                      gamma0(ind1p,ind2,ind3,ind4p) =
     $                   gamma0(ind2,ind4,ind3,ind1)

                      gamma0(ind1p,ind2p,ind3p,ind4p) = 
     $                   gamma0(ind4,ind3,ind2,ind1)

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

              gamma0_ph( 4*nb*na1+na3, 4*nb*na4+na2) = 
     $           gamma0(na1,na2,na3,na4)

            enddo
          enddo
        enddo
      enddo

      return
      end
