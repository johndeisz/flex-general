#include "../convert.F"

      subroutine discon_lat(h0, delta_gl_k, delta_glp_k)

#include "../constants.F"

      COMPLEX, dimension (0:4*nb-1,0:4*nb-1,0:nl-1) ::
     $   h0, delta_gl_k, delta_glp_k

      INTEGER na1, na2

c     Express the discontinuities in k-space
      
      delta_glp_k = h0

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
