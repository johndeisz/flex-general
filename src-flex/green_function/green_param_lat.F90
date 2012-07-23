#include "../convert.F90"

MODULE green_param_lat

CONTAINS

  FUNCTION c_lattice_k(kl, tij, ed, v_pert_eff, psi, h_eff, prfld_eff, mu, &
       sigma1, h_so)  

    USE CONSTANTS
    USE h_zero

    IMPLICIT NONE

    COMPLEX, DIMENSION (0:1,0:1,0:4*nb-1,0:4*nb-1) :: c_lattice_k
    INTEGER kl
    REAL ed(0:nb-1)
    COMPLEX tij(0:nb-1,0:nb-1,-2:2,-2:2,-2:2)
    REAL v_pert_eff(0:nb-1)
    COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)
    REAL h_eff(0:nb-1,1:3)
    REAL prfld_eff
    REAL mu
    COMPLEX, dimension (0:4*nb-1,0:4*nb-1) :: sigma1
    COMPLEX h_so(0:2*nb-1, 0:2*nb-1)
    COMPLEX, dimension (0:4*nb-1,0:4*nb-1) :: delta_gl_k, delta_glp_k
    INTEGER na1
    INTEGER ib, jb, is
    
    !     Express the discontinuities in k-space
      
    delta_glp_k(:,:) = h0(kl, tij, ed, v_pert_eff, psi, h_eff, & 	
         prfld_eff, mu, sigma1, h_so)

    delta_gl_k(:,:) = cmplx(0.0d0, 0.0d0)

    do na1 = 0, 4*nb-1
       delta_gl_k(na1,na1) = cmplx(-1.0d0,0.0d0)
    enddo

    c_lattice_k(0,0,:,:) = -delta_gl_k
    c_lattice_k(0,1,:,:) = cmplx(0.0d0, 0.0d0)
    c_lattice_k(1,0,:,:) = delta_glp_k
    c_lattice_k(1,1,:,:) = cmplx(0.0d0, 0.0d0)

  END FUNCTION c_lattice_k

END MODULE green_param_lat



