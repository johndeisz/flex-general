#include "../convert.F90"

subroutine dyson(myrank, g, q_tau, q_epsilon, tij, ed, v_pert_eff, &
     psi, h_eff, prfld_eff, mu, sigma1, h_so, sigma, epsilon, t)

  USE CONSTANTS
  USE h_zero

  IMPLICIT NONE 

#ifdef USE_MPI
  include 'mpif.h'
#endif /* USE_MPI */

  INTEGER myrank
  COMPLEX, dimension (0:4*nb-1,0:4*nb-1,0:mp1,0:nc1) :: g,sigma

  REAL ed(0:nb-1)
  COMPLEX tij(0:nb-1,0:nb-1,-2:2,-2:2,-2:2)
  REAL v_pert_eff(0:nb-1)
  COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)
  REAL h_eff(0:nb-1,1:3)
  REAL prfld_eff
  REAL mu
  COMPLEX, dimension (0:4*nb-1,0:4*nb-1) :: sigma1
  COMPLEX h_so(0:2*nb-1, 0:2*nb-1)

  REAL q_tau(0:1,0:1,0:mp1)
  COMPLEX q_epsilon(0:1,0:1,0:mp1)
  REAL epsilon(0:mp1)
  REAL t

  INTEGER x_stretch, y_stretch, z_stretch
  INTEGER coarse_grain_points
  REAL wx(-llx/lcx: llx/lcx)
  REAL wy(-lly/lcy: lly/lcy)
  REAL wz(-llz/lcz: llz/lcz)
  INTEGER ix_max, iy_max, iz_max
  COMPLEX e_n(0:4*nb-1, 0:4*nb-1)

  !     Counters
  INTEGER ix, iy, iz
  INTEGER i, j
  INTEGER klx, kly, klz, kl
  INTEGER kx, ky, kz, k, l
  INTEGER k1, k2, k3
  INTEGER ia, ib

  !     Arrays
  COMPLEX, dimension (0:4*nb-1,0:4*nb-1) :: identity, temp_gc, temp_gl2, &
       temp_gl

  !   Lapack related
  INTEGER ipiv(4*nb)
  INTEGER info

#ifdef USE_MPI
  INTEGER dummy, i_proc, ierr
  INTEGER stat(MPI_STATUS_SIZE)
  COMPLEX g_temp0(0:4*nb-1,0:4*nb-1,0:nl-1)
#endif /* USE_MPI */

  x_stretch = llx / lcx
  y_stretch = lly / lcy
  z_stretch = llz / lcz

  coarse_grain_points = x_stretch * y_stretch * z_stretch

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

  do i = 0, 4*nb-1
     do j = 0, 4*nb-1
        if (i .eq. j) then
           identity(i,j) = cmplx(1.0d0, 0.0d0)
        else
           identity(i,j) = cmplx(0.0d0, 0.0d0)
        endif
     enddo
  enddo

  !     Loop over cluster k-points.
  do kx = 0, lcx-1
     do ky = 0, lcy-1
        do kz = 0, lcz-1

           k = kx + ky*lcx + kz*lcx*lcy
           
           do l = 0, mp1

              do i = 0, 4*nb-1
                 do j = 0, 4*nb-1
                    if (i .eq. j) then
                       e_n(i,j) = epsilon(l)*cmplx(0.0d0, 1.0d0)
                    else
                       e_n(i,j) = epsilon(l)*cmplx(0.0d0, 0.0d0)
                    endif
                 enddo
              enddo

              temp_gc = cmplx(0.0d0, 0.0d0)

              !     For each cluster k-point, sum over all points in the
              !     coarse-graining cell.

              do ix = -ix_max, ix_max
                 do iy = -iy_max, iy_max
                    do iz = -iz_max, iz_max

                       klx = mod(kx*x_stretch + ix + llx, llx)   
                       kly = mod(ky*y_stretch + iy + lly, lly)
                       klz = mod(kz*z_stretch + iz + llz, llz)
                       
                       kl = klx + kly*llx + klz*llx*lly

                       temp_gl2 = e_n - & 
                            h0(kl,tij, ed, v_pert_eff, psi, h_eff, &  	
                            prfld_eff, mu, sigma1, h_so)  -  sigma(:,:,l,k)
                       temp_gl = identity

                       !     Invert temp_gl
                       call cgesv(4*nb, 4*nb, temp_gl2, 4*nb, & 
                            ipiv, temp_gl, 4*nb, info)

                       if (info .ne. 0) then
                          write(6,*) 'info not equal to zero'
                       endif

                       temp_gc = temp_gc + wx(ix) * wy(iy) * wz(iz) * temp_gl

                    enddo
                 enddo
              enddo

              g(:,:,l,k) = temp_gc / float(coarse_grain_points)

           enddo
        enddo
     enddo
  enddo

  return
end subroutine dyson
