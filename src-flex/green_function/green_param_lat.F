#include "../convert.F"

      subroutine green_param_lat(cl_k, delta_gl_k, delta_glp_k)  

#include "../constants.F"

      COMPLEX cl_k(0:1,0:1,0:4*nb-1,0:4*nb-1,0:nl-1)
      COMPLEX, dimension (0:4*nb-1,0:4*nb-1, 0:nl-1) ::
     $   delta_gl_k, delta_glp_k

c     Initialize c_k_lat so that the asymptotic 
c     behavior of the lattice green's function is accounted for.

      cl_k(0,0,:,:,:) = -delta_gl_k
      cl_k(0,1,:,:,:) = cmplx(0.0d0, 0.0d0)
      cl_k(1,0,:,:,:) = delta_glp_k
      cl_k(1,1,:,:,:) = cmplx(0.0d0, 0.0d0)

      return
      end



