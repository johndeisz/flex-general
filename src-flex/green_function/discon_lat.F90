#include "../convert.F90"

      subroutine discon_lat(ek, v_pert_eff, psi, h_eff, 	
     $   prfld_eff, mu, sigma1, h_so, 
     $   delta_gl_k, delta_glp_k)

#include "../constants.F90"

      EXTERNAL h0
      COMPLEX h0(0:4*nb-1,0:4*nb-1)

      COMPLEX ek(0:nb-1,0:nb-1,0:nl-1)
      REAL v_pert_eff(0:nb-1)
      COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)
      REAL h_eff(0:nb-1,1:3)
      REAL prfld_eff
      REAL mu
      COMPLEX, dimension (0:4*nb-1,0:4*nb-1) :: sigma1
      COMPLEX h_so(0:2*nb-1, 0:2*nb-1)

      COMPLEX, dimension (0:4*nb-1,0:4*nb-1,0:nl-1) ::
     $   delta_gl_k, delta_glp_k

      INTEGER na1, na2
      INTEGER k

c     Express the discontinuities in k-space
      
      do k = 0, nl - 1
        delta_glp_k(:,:,k) = h0(k, ek, v_pert_eff, psi, h_eff,  	
     $   prfld_eff, mu, sigma1, h_so)
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
      end
