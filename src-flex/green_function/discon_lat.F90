#include "../convert.F90"

subroutine discon_lat(tij, ed, v_pert_eff, psi, h_eff, prfld_eff, mu, &
     sigma1, h_so, delta_gl_k, delta_glp_k)

  USE CONSTANTS
  USE h_zero

  REAL ed(0:nb-1)
  COMPLEX tij(0:nb-1,0:nb-1,-2:2,-2:2,-2:2)
  REAL v_pert_eff(0:nb-1)
  COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)
  REAL h_eff(0:nb-1,1:3)
  REAL prfld_eff
  REAL mu
  COMPLEX, dimension (0:4*nb-1,0:4*nb-1) :: sigma1
  COMPLEX h_so(0:2*nb-1, 0:2*nb-1)

  COMPLEX, dimension (0:4*nb-1,0:4*nb-1,0:nl-1) :: delta_gl_k, delta_glp_k

  INTEGER na1, na2
  INTEGER k

  !     Express the discontinuities in k-space
      
  do k = 0, nl - 1
     delta_glp_k(:,:,k) = h0(k, tij, ed, v_pert_eff, psi, h_eff, & 	
          prfld_eff, mu, sigma1, h_so)
  end do

  do na1 = 0, 4*nb-1
     do na2 = 0, 4*nb-1

        if (na1 .eq. na2) then

           delta_gl_k(na1,na2,:) = cmplx(-1.0d0,0.0d0)

        else

           delta_gl_k(na1,na2,:) = cmplx(0.0d0,0.0d0)

        endif

     enddo
  enddo

  return
end subroutine discon_lat
