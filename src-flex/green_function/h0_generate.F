#include "../convert.F"

      subroutine h0_generate(h0, ek, v_pert_eff, psi, h_eff,  	
     $   prfld_eff, mu, sigma1, h_so, iteration)

c /**************************************************************
c * H0 -- program to calculate the H0 array                    *
c *							      *
c * Author: Tom Slife	                 		      *
c  **************************************************************/

#include "../constants.F"

      COMPLEX h0(0:4*nb-1, 0:4*nb-1, 0:nl-1)
      COMPLEX ek(0:nb-1,0:nb-1,0:nl-1)
      COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)
      REAL mu
      COMPLEX sigma1(0:4*nb-1,0:4*nb-1)
      COMPLEX h_so(0:2*nb-1, 0:2*nb-1)
      COMPLEX Lvec(1:3,0:nb-1,0:nb-1)
      INTEGER iteration

      REAL, dimension (0:nb-1,1:3) :: h_eff
      REAL prfld_eff
      REAL v_pert_eff(0:nb-1)

      REAL, dimension (0:nb-1, 0:nb-1) :: delta

      INTEGER nu1, nu2, k, k_minus, k1, k2, k3
      INTEGER is1, is2
      INTEGER ind1, ind2

      delta = 0.0d0

      do nu1 = 0, nb-1
        delta(nu1,nu1) = 1.0d0
      enddo

c Include anomolous part of the g-factor for spin coupling
      h_eff = h_eff * (gs/2.0d0)

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
     $             cmplx((v_pert_eff(nu1)-mu -h_eff(nu1,3))*
     $             delta(nu1,nu2), 0.0d0)

                h0(ind1 + 0, ind2 + 1, k) = 
     $             cmplx(-h_eff(nu1,1), h_eff(nu1,2))*
     $              delta(nu1,nu2)

                h0(ind1 + 0, ind2 + 2, k) = prfld_eff *
     $             conjg( psi( 2*nu1, 2*nu2, k_minus) )

                h0(ind1 + 0, ind2 + 3, k) = prfld_eff *
     $             conjg( psi( 2*nu1, 2*nu2+1, k_minus) )
	
                h0(ind1 + 1, ind2 + 0, k) = 
     $             cmplx(-h_eff(nu1,1),-h_eff(nu1,2))*
     $                delta(nu1,nu2)

                h0(ind1 + 1, ind2 + 1, k) =  ek(nu1,nu2,k) +
     $             cmplx( (v_pert_eff(nu1)-mu + h_eff(nu1,3))*
     $             delta(nu1,nu2), 0.0d0)

                h0(ind1 + 1, ind2 + 2, k) = prfld_eff * 
     $             conjg( psi( 2*nu1+1, 2*nu2, k_minus) )

                h0(ind1 + 1, ind2 + 3, k) = prfld_eff *
     $             conjg( psi( 2*nu1+1, 2*nu2+1, k_minus) )

                h0(ind1 + 2, ind2 + 0, k) = prfld_eff *
     $             psi( 2*nu2, 2*nu1, k_minus)

                h0(ind1 + 2, ind2 + 1, k) = prfld_eff *
     $             psi( 2*nu2+1, 2*nu1, k_minus)

                h0(ind1 + 2, ind2 + 2, k) = -ek(nu2,nu1,k_minus) 
     $             + cmplx((-v_pert_eff(nu1)+mu + h_eff(nu1,3))*
     $             delta(nu1,nu2), 0.0d0 )

                h0(ind1 + 2, ind2 + 3, k) = 
     $             cmplx(h_eff(nu1,1), h_eff(nu1,2))*
     $              delta(nu1,nu2)

                h0(ind1 + 3, ind2 + 0, k) = prfld_eff *
     $             psi( 2*nu2, 2*nu1+1, k_minus)

                h0(ind1 + 3, ind2 + 1, k) = prfld_eff *
     $             psi( 2*nu2+1, 2*nu1+1, k_minus)

                h0(ind1 + 3, ind2 + 2, k) = 
     $             cmplx(h_eff(nu1,1),-h_eff(nu1,2))*
     $              delta(nu1,nu2)

                h0(ind1 + 3, ind2 + 3, k) = 
     $             -ek(nu2,nu1,k_minus) +  
     $             cmplx((-v_pert_eff(nu1)+mu - h_eff(nu1,3))*
     $              delta(nu1,nu2), 0.0d0 )

              enddo
              
            enddo

            h0(:,:,k) = h0(:,:,k) + sigma1

          enddo
        enddo
      enddo

c    Remove spin g-factor correction from magnetic field:
      h_eff = h_eff / (gs/2.0d0)

      call angular_matrices(Lvec)

      do nu1 = 0, nb-1
        do nu2 = 0, nb-1
          do is1 = 0, 1

            h0(4*nu1+is1, 4*nu2+is1,:) = 
     $         h0(4*nu1+is1, 4*nu2+is1,:) -
     $         h_eff(nu1,1)*Lvec(1,nu1,nu2) -
     $         h_eff(nu1,2)*Lvec(2,nu1,nu2) -
     $         h_eff(nu1,3)*Lvec(3,nu1,nu2)


            h0(4*nu1+is1+2, 4*nu2+is1+2,:) = 
     $         h0(4*nu1+is1+2, 4*nu2+is1+2,:) +
     $         h_eff(nu1,1)*conjg(Lvec(1,nu1,nu2)) +
     $         h_eff(nu1,2)*conjg(Lvec(2,nu1,nu2)) +
     $         h_eff(nu1,3)*conjg(Lvec(3,nu1,nu2))
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
