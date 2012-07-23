module h_zero

#include "../convert.F90"

contains

  function h0(k, tij, ed, v_pert_eff, psi, h_eff, prfld_eff, mu, sigma1, h_so)

    USE CONSTANTS
    USE bare_dispersion

    IMPLICIT NONE

    COMPLEX, dimension (0:4*nb-1, 0:4*nb-1) :: h0
    INTEGER k
    REAL ed(0:nb-1)
    COMPLEX tij(0:nb-1,0:nb-1,-2:2,-2:2,-2:2)
    COMPLEX ek(0:nb-1,0:nb-1), ek_minus(0:nb-1, 0:nb-1)
    COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)
    REAL mu
    COMPLEX sigma1(0:4*nb-1,0:4*nb-1)
    COMPLEX h_so(0:2*nb-1, 0:2*nb-1)
    COMPLEX Lvec(1:3,0:nb-1,0:nb-1)

    REAL, dimension (0:nb-1,1:3) :: h_eff
    REAL prfld_eff
    REAL v_pert_eff(0:nb-1)

    REAL, dimension (0:nb-1, 0:nb-1) :: delta

    INTEGER k_minus, k1, k2, k3
    INTEGER nu1, nu2
    INTEGER is1, is2
    INTEGER ind1, ind2

    delta = 0.0d0

    do nu1 = 0, nb-1
       delta(nu1,nu1) = 1.0d0
    enddo

    ! Include anomolous part of the g-factor for spin coupling
    h_eff = h_eff * (gs/2.0d0)

    k3 = int(k/(llx*lly))
    k2 = int( (k-k3*llx*lly) / llx )
    k1 = k - k3*llx*lly - k2*llx

    k_minus = mod(llz-k3,llz)*(llx*lly) + mod(lly-k2,lly)*llx + mod(llx-k1,llx)
    ek = ekl(k, tij, ed)
    ek_minus = ekl(k_minus, tij, ed)

    do nu1 = 0, nb-1
       do nu2 = 0, nb-1

          ind1 = 4 * nu1
          ind2 = 4 * nu2

          h0(ind1 + 0, ind2 + 0) = ek(nu1,nu2)  + &
               cmplx((v_pert_eff(nu1)-mu -h_eff(nu1,3))* &
               delta(nu1,nu2), 0.0d0)

          h0(ind1 + 0, ind2 + 1) = cmplx(-h_eff(nu1,1), h_eff(nu1,2))* &
               delta(nu1,nu2)

          h0(ind1 + 0, ind2 + 2) = prfld_eff * &
               conjg( psi( 2*nu1, 2*nu2, k_minus) )
          
          h0(ind1 + 0, ind2 + 3) = prfld_eff * &
               conjg( psi( 2*nu1, 2*nu2+1, k_minus) )
	
          h0(ind1 + 1, ind2 + 0) = & 
               cmplx(-h_eff(nu1,1),-h_eff(nu1,2))* delta(nu1,nu2)

          h0(ind1 + 1, ind2 + 1) =  ek(nu1,nu2) + &
               cmplx( (v_pert_eff(nu1)-mu + h_eff(nu1,3))* &
               delta(nu1,nu2), 0.0d0)

          h0(ind1 + 1, ind2 + 2) = prfld_eff * & 
               conjg( psi( 2*nu1+1, 2*nu2, k_minus) )

          h0(ind1 + 1, ind2 + 3) = prfld_eff * &
               conjg( psi( 2*nu1+1, 2*nu2+1, k_minus) )

          h0(ind1 + 2, ind2 + 0) = prfld_eff * &
               psi( 2*nu2, 2*nu1, k_minus)

          h0(ind1 + 2, ind2 + 1) = prfld_eff * &
               psi( 2*nu2+1, 2*nu1, k_minus)

          h0(ind1 + 2, ind2 + 2) = -ek_minus(nu2,nu1) & 
               + cmplx((-v_pert_eff(nu1)+mu + h_eff(nu1,3))* &
               delta(nu1,nu2), 0.0d0 )

          h0(ind1 + 2, ind2 + 3) = cmplx(h_eff(nu1,1), h_eff(nu1,2))* &
               delta(nu1,nu2)

          h0(ind1 + 3, ind2 + 0) = prfld_eff * &
               psi( 2*nu2, 2*nu1+1, k_minus)

          h0(ind1 + 3, ind2 + 1) = prfld_eff * &
               psi( 2*nu2+1, 2*nu1+1, k_minus)

          h0(ind1 + 3, ind2 + 2) = & 
               cmplx(h_eff(nu1,1),-h_eff(nu1,2))* delta(nu1,nu2)

          h0(ind1 + 3, ind2 + 3) = &
               -ek_minus(nu2,nu1) + &  
               cmplx((-v_pert_eff(nu1)+mu - h_eff(nu1,3))* &
               delta(nu1,nu2), 0.0d0 )

       enddo
    enddo
              
    h0 = h0 + sigma1
        
    !    Remove spin g-factor correction from magnetic field:
    h_eff = h_eff / (gs/2.0d0)

    call angular_matrices(Lvec)

    do nu1 = 0, nb-1
       do nu2 = 0, nb-1
          do is1 = 0, 1

             h0(4*nu1+is1, 4*nu2+is1) =  h0(4*nu1+is1, 4*nu2+is1) - &
                  h_eff(nu1,1)*Lvec(1,nu1,nu2) - &
                  h_eff(nu1,2)*Lvec(2,nu1,nu2) - &
                  h_eff(nu1,3)*Lvec(3,nu1,nu2)

             h0(4*nu1+is1+2, 4*nu2+is1+2) = h0(4*nu1+is1+2, 4*nu2+is1+2) + &
                  h_eff(nu1,1)*conjg(Lvec(1,nu1,nu2)) + &
                  h_eff(nu1,2)*conjg(Lvec(2,nu1,nu2)) + &
                  h_eff(nu1,3)*conjg(Lvec(3,nu1,nu2))
          enddo
       enddo
    enddo

    !     Add spin-orbit interaction
    do nu1 = 0, nb-1
       do is1 = 0, 1
          do nu2 = 0, nb-1
             do is2 = 0, 1

                h0(4*nu1+is1, 4*nu2+is2) = h0(4*nu1+is1, 4*nu2+is2) + & 
                     h_so(2*nu1+is1,2*nu2+is2) 

                h0(4*nu1+is1+2, 4*nu2+is2+2) = h0(4*nu1+is1+2, 4*nu2+is2+2)  - &
                     conjg(h_so(2*nu1+is1,2*nu2+is2))

             enddo
          enddo
       enddo
    enddo
              
  end function h0
end module h_zero
