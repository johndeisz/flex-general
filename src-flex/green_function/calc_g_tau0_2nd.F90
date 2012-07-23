#include "../convert.F90"

subroutine calc_g_tau0_2nd(myrank, g_ktau0, q_tau, q_epsilon, tij, ed, &
     v_pert_eff, psi, h_eff, prfld_eff, mu, sigma1, h_so, sigma, epsilon, t)

  USE CONSTANTS
  USE h_zero
  USE green_param_lat

  IMPLICIT NONE

#ifdef USE_MPI
  include 'mpif.h'
#endif /* USE_MPI */

  INTEGER myrank
  COMPLEX, dimension (0:4*nb-1,0:4*nb-1,0:mp1,0:nc1) :: sigma
  COMPLEX, dimension (0:4*nb-1,0:4*nb-1,0:nl-1) :: g_ktau0

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
  COMPLEX cl_k(0:1,0:1,0:4*nb-1,0:4*nb-1)

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
  INTEGER ia, ib, is, jb

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

  g_ktau0 = cmplx(0.0d0, 0.0d0)

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

                       cl_k = c_lattice_k(kl,tij, ed, v_pert_eff, psi, &
                            h_eff, prfld_eff, mu, sigma1, h_so)
                       
                       do ia = 0, 1
                          do ib = 0, 1
                             temp_gl = temp_gl - cl_k(ia,ib,:,:) * &
                                  q_epsilon(ia,ib,l)
                          enddo
                       enddo

                       !     Obtain g_tau0 by summing over all frequencies.  
                       !     Include
                       !     weights wx,wy,wz as the lattice 
                       !     greens function may appear
                       !     in several course graining cells.

                       g_ktau0(:,:,kl) =  g_ktau0(:,:,kl) + & 
                            t * wx(ix) * wy(iy) * wz(iz) * temp_gl

                    enddo
                 enddo
              enddo

           enddo
        enddo

     enddo
  enddo

  !     If NP > 1 collect, let myrank = 0 collect all the contributions
  !     to sum_e g(e,k)

  !     Add contributions to g_tau0 from all other processes
#ifdef USE_MPI
  if (myrank .eq. 0) then

     !     Copy g_temp0 contribution to cl_k(0,0)

     dummy = 0

     do i_proc = 1, np - 1

        call MPI_Send(dummy, 1, MPI_INTEGER, i_proc, 0, MPI_COMM_WORLD, ierr)
        
        call MPI_Recv(g_temp0, 16*nb*nb*nl, MPI_COMPLEX, & 
             i_proc, 1, MPI_COMM_WORLD, stat, ierr)

        g_ktau0 = g_ktau0 + g_temp0

     enddo

  else 

     call MPI_Recv(dummy, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, stat, ierr)
     call MPI_Send(g_ktau0, 16*nb*nb*nl, MPI_COMPLEX, & 
          0, 1, MPI_COMM_WORLD, ierr)

  endif
#endif /* USE_MPI */

  !     Add back analytic contributions

  if (myrank .eq. 0) then

     do kl = 0, nl-1
        cl_k = c_lattice_k(kl,tij, ed, v_pert_eff, psi, &
             h_eff, prfld_eff, mu, sigma1, h_so)
  
        do ia = 0, 1
           do ib = 0, 1
              g_ktau0(:,:,kl) =  g_ktau0(:,:,kl) + &
                   cl_k(ia,ib,:,:) * q_tau(ia,ib,0)
           enddo
        enddo
    
     enddo

  endif

  return
end subroutine calc_g_tau0_2nd
