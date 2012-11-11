#include "../convert.F90"

subroutine discontinuities(tij, ed, v_pert_eff, psi, h_eff, prfld_eff, &
     mu, sigma1, h_so, delta_g_r, delta_g_k, delta_gp_r, delta_gp_k)

  USE CONSTANTS
  USE h_zero
  USE bare_dispersion
  IMPLICIT NONE
  
  REAL ed(0:nb-1)
  COMPLEX tij(0:nb-1,0:nb-1,-2:2,-2:2,-2:2)
  REAL v_pert_eff(0:nb-1)
  COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)
  REAL h_eff(0:nb-1,1:3)
  REAL prfld_eff
  REAL mu
  COMPLEX, dimension (0:4*nb-1,0:4*nb-1) :: sigma1
  COMPLEX h_so(0:2*nb-1, 0:2*nb-1)

  COMPLEX, dimension (0:4*nb-1,0:4*nb-1,0:nc-1) :: &
       delta_g_r, delta_g_k, delta_gp_r, delta_gp_k

  COMPLEX temp_a(0:nc-1)
  COMPLEX temp_b(0:nc-1)

  INTEGER ifftv(0:3)
  INTEGER na1, na2, k1, k2, k3, k

  INTEGER x_stretch, y_stretch, z_stretch
  INTEGER course_grain_points
  INTEGER ix, iy, iz, kl1, kl2, kl3
  INTEGER ix_max, iy_max, iz_max

  REAL wx(-llx/lcx: llx/lcx)
  REAL wy(-lly/lcy: lly/lcy)
  REAL wz(-llz/lcz: llz/lcz)
  
  INTEGER klx, kly, klz, kl
  COMPLEX temp_gp_k(0:4*nb-1,0:4*nb-1)

  !     Express the discontinuities in k-space

  !     Compute weights for k-points within the course graining
  !     cell.  The following assigns a weight of unity to all
  !     points within the cell, a weight of 1/2 on a cube face,
  !     a weight of 1/4 on a cube edge, and 1/8 for a cube corner.

  x_stretch = llx / lcx
  y_stretch = lly / lcy
  z_stretch = llz / lcz

  course_grain_points = x_stretch * y_stretch * z_stretch

  if (x_stretch .gt. 1) then
     ix_max = x_stretch/2
     do ix = 0, ix_max - 1
        wx(ix) = 1.0d0 
        wx(-ix) = 1.0d0 
     enddo
     wx(ix_max) = 0.5d0 
     wx(-ix_max) = 0.5d0 
  else
     ix_max = 0
     wx(0) = 1.0d0
  endif

  if (y_stretch .gt. 1) then
     iy_max = y_stretch/2
     do iy = 0, iy_max - 1
        wy(iy) = 1.0d0 
        wy(-iy) = 1.0d0 
     enddo
     wy(iy_max) = 0.5d0 
     wy(-iy_max) = 0.50d0 
  else
     iy_max = 0
     wy(0) = 1.0d0
  endif

  if (z_stretch .gt. 1) then
     iz_max = z_stretch/2
     do iz = 0, iz_max - 1
        wz(iz) = 1.0d0 
        wz(-iz) = 1.0d0 
     enddo
     wz(iz_max) = 0.5d0 
     wz(-iz_max) = 0.5d0 
  else
     iz_max = 0
     wz(0) = 1.0d0
  endif

  do k1 = 0, lcx-1
     do k2 = 0, lcy-1
        do k3 = 0, lcz-1

           k = k1 + k2*lcx + k3*lcx*lcy
           
           do na1 = 0, 4*nb-1
              do na2 = 0, 4*nb-1
 
                 if (na1 .eq. na2) then
                    delta_g_k(na1,na2,k) = -1.0d0
                 else
                    delta_g_k(na1,na2,k) = -0.0d0
                 endif
                
              enddo
           enddo

           !     For each cluster k-point, sum over all points in the
           !     course-graining cell.

           temp_gp_k = cmplx(0.0d0,0.0d0)

           do ix = -ix_max, ix_max
              do iy = -iy_max, iy_max
                 do iz = -iz_max, iz_max
                   
                    klx = mod(k1*x_stretch + ix + llx, llx)   
                    kly = mod(k2*y_stretch + iy + lly, lly)
                    klz = mod(k3*z_stretch + iz + llz, llz)

                    kl = klx + kly*llx + klz*llx*lly
 
                    temp_gp_k = temp_gp_k +  wx(ix) * wy(iy) * wz(iz) * &
                         h0(kl,  tij, ed, v_pert_eff, psi, h_eff, &  	
                         prfld_eff, mu, sigma1, h_so)

                 enddo
              enddo
           enddo

           delta_gp_k(:,:,k) = temp_gp_k / float(course_grain_points)
            
        enddo
     enddo
  enddo

  !     Fourier transform to obtain the discontinuities in r-space

  ifftv(1) =  1
  ifftv(2) =  1
  ifftv(3) =  1

  do na1 = 0, 4*nb-1
     do na2 = 0, 4*nb-1

        temp_a = delta_g_k(na1,na2,:)
        temp_b = delta_gp_k(na1,na2,:)

        call fft_3D ( temp_a, ifftv)
        call fft_3D ( temp_b, ifftv)

        delta_g_r(na1,na2,:) = temp_a / float(nc)
        delta_gp_r(na1,na2,:) = temp_b / float(nc)

     enddo
  enddo

  return
end subroutine discontinuities
