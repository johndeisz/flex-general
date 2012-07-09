      subroutine h0_generate(h0, ek, psi, h, prfld, sigma1, h_so)

#include "constants.F"

      DOUBLE COMPLEX h0(0:4*nb-1, 0:4*nb-1, 0:nl-1)
      DOUBLE COMPLEX ek(0:nb-1,0:nb-1,0:nl-1)
      DOUBLE COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)
      DOUBLE PRECISION, dimension (0:nb-1,1:3) :: h
      DOUBLE PRECISION prfld
      DOUBLE COMPLEX sigma1(0:4*nb-1,0:4*nb-1)
      DOUBLE COMPLEX h_so(0:2*nb-1, 0:2*nb-1)

      DOUBLE PRECISION, dimension (0:nb-1, 0:nb-1) :: delta
      DOUBLE COMPLEX Lvec(1:3,0:nb-1,0:nb-1)

      INTEGER nu1, nu2, k, k_minus, k1, k2, k3
      INTEGER is1, is2
      INTEGER ind1, ind2

      delta = 0.0d0

      do nu1 = 0, nb-1
        delta(nu1,nu1) = 1.0d0
      enddo

      do k3 = 0, llz1
        do k2 = 0, lly1
          do k1 = 0, llx1

            k = k3*(llx*lly) + k2*llx + k1
            k_minus = mod(llz-k3,llz)*(llx*lly) +
     $         mod(lly-k2,lly)*llx + mod(llx-k1,llx)

            do nu1 = 0, nb-1
              do nu2 = 0, nb-1

                ind1 = 4 * nu1
                ind2 = 4 * nu2

                h0(ind1 + 0, ind2 + 0, k) = ek(nu1,nu2,k)  +
     $             cmplx((-h(nu1,3))*
     $             delta(nu1,nu2), 0.0d0)

                h0(ind1 + 0, ind2 + 1, k) = 
     $             cmplx(-h(nu1,1), h(nu1,2))*
     $              delta(nu1,nu2)

                h0(ind1 + 0, ind2 + 2, k) = prfld *
     $             conjg( psi( 2*nu1, 2*nu2, k_minus) )

                h0(ind1 + 0, ind2 + 3, k) = prfld *
     $             conjg( psi( 2*nu1, 2*nu2+1, k_minus) )
	
                h0(ind1 + 1, ind2 + 0, k) = 
     $             cmplx(-h(nu1,1),-h(nu1,2))*
     $                delta(nu1,nu2)

                h0(ind1 + 1, ind2 + 1, k) =  ek(nu1,nu2,k) +
     $             cmplx( ( h(nu1,3))*
     $             delta(nu1,nu2), 0.0d0)

                h0(ind1 + 1, ind2 + 2, k) = prfld * 
     $             conjg( psi( 2*nu1+1, 2*nu2, k_minus) )

                h0(ind1 + 1, ind2 + 3, k) = prfld *
     $             conjg( psi( 2*nu1+1, 2*nu2+1, k_minus) )

                h0(ind1 + 2, ind2 + 0, k) = prfld *
     $             psi( 2*nu2, 2*nu1, k_minus)

                h0(ind1 + 2, ind2 + 1, k) = prfld *
     $             psi( 2*nu2+1, 2*nu1, k_minus)

                h0(ind1 + 2, ind2 + 2, k) = -ek(nu2,nu1,k_minus) 
     $             + cmplx((h(nu1,3))*
     $             delta(nu1,nu2), 0.0d0 )

                h0(ind1 + 2, ind2 + 3, k) = 
     $             cmplx(h(nu1,1), h(nu1,2))*
     $              delta(nu1,nu2)

                h0(ind1 + 3, ind2 + 0, k) = prfld *
     $             psi( 2*nu2, 2*nu1+1, k_minus)

                h0(ind1 + 3, ind2 + 1, k) = prfld *
     $             psi( 2*nu2+1, 2*nu1+1, k_minus)

                h0(ind1 + 3, ind2 + 2, k) = 
     $             cmplx(h(nu1,1),-h(nu1,2))*
     $              delta(nu1,nu2)

                h0(ind1 + 3, ind2 + 3, k) = 
     $             -ek(nu2,nu1,k_minus) +  
     $             cmplx((- h(nu1,3))*
     $              delta(nu1,nu2), 0.0d0 )

              enddo
              
            enddo

            h0(:,:,k) = h0(:,:,k) + sigma1

          enddo
        enddo
      enddo

c    Remove spin g-factor correction from magnetic field:
      h = h / (gs/2.0d0)

      call angular_matrices(Lvec)

      do nu1 = 0, nb-1
        do nu2 = 0, nb-1
          do is1 = 0, 1

            h0(4*nu1+is1, 4*nu2+is1,:) =
     $         h0(4*nu1+is1, 4*nu2+is1,:) -
     $         h(nu1,1)*Lvec(1,nu1,nu2) -
     $         h(nu1,2)*Lvec(2,nu1,nu2) -
     $         h(nu1,3)*Lvec(3,nu1,nu2)


            h0(4*nu1+is1+2, 4*nu2+is1+2,:) =
     $         h0(4*nu1+is1+2, 4*nu2+is1+2,:) +
     $         h(nu1,1)*conjg(Lvec(1,nu1,nu2)) +
     $         h(nu1,2)*conjg(Lvec(2,nu1,nu2)) +
     $         h(nu1,3)*conjg(Lvec(3,nu1,nu2))
          enddo
        enddo
      enddo

c     Add spin-orbit interaction
      do nu1 = 0, nb-1
        do is1 = 0, 1
          do nu2 = 0, nb-1
            do is2 = 0, 1

              h0(4*nu1+is1, 4*nu2+is2,:) =
     $           h0(4*nu1+is1, 4*nu2+is2,:) +
     $           h_so(2*nu1+is1,2*nu2+is2)

              h0(4*nu1+is1+2, 4*nu2+is2+2,:) =
     $           h0(4*nu1+is1+2, 4*nu2+is2+2,:)  -
     $           conjg(h_so(2*nu1+is1,2*nu2+is2))

            enddo
          enddo
        enddo
      enddo

      return
      end
