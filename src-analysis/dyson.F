      subroutine dyson(h0k)

#include "../constants.F"

      DOUBLE COMPLEX, dimension (0:4*nb-1,0:4*nb-1) :: h0k
      DOUBLE COMPLEX, dimension (0:4*nb-1,0:4*nb-1, 0:nc-1, 
     $   0:n_pade_max) :: p, q

      DOUBLE COMPLEX, dimension (0:4*nb-1,0:4*nb-1) :: sigma
 
      REAL epsilon(0:mp1)
      REAL t

      INTEGER x_stretch, y_stretch, z_stretch
      REAL wx(-llx/lcx: llx/lcx)
      REAL wy(-lly/lcy: lly/lcy)
      REAL wz(-llz/lcz: llz/lcz)
      INTEGER ix_max, iy_max, iz_max
      COMPLEX e_n(0:4*nb-1, 0:4*nb-1)

c     Counters
      INTEGER ix, iy, iz
      INTEGER i, j
      INTEGER klx, kly, klz, kl
      INTEGER kx, ky, kz, k, l
      INTEGER k1, k2, k3
      INTEGER ia, ib

c     Arrays
      COMPLEX, dimension (0:4*nb-1,0:4*nb-1) ::
     $   identity, temp_gc, temp_gl2, temp_gl

c     Lapack related
      INTEGER ipiv(4*nb)
      INTEGER info

      p = dcmplx(0.0d0, 0.0d0)
      q = dcmplx(0.0d0, 0.0d0)

c     Read pade coefficients
      do nu1 = 0, nb-1
        do ia1 = 0, 3
          na1 = 4*nu1 + ia1
          do nu2 = 0, nb-1
            do ia2 = 0, 3
              na2 = 4*nu2 + ia2
              do k = 0, nc-1
                read(50,*) id1, id2, id3, id4, id5, n_pade(k)
                if (n_pade(k) .gt. 0) then
                  do l = 0, n_pade(k)
                    read(50,*) p(na1, na2, k, l)
                  enddo
                endif
              enddo
            enddo
          enddo
        enddo
      enddo

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

c$$$      do i = 0, 4*nb-1
c$$$        do j = 0, 4*nb-1
c$$$          if (i .eq. j) then
c$$$            identity(i,j) = cmplx(1.0d0, 0.0d0)
c$$$          else
c$$$            identity(i,j) = cmplx(0.0d0, 0.0d0)
c$$$          endif
c$$$        enddo
c$$$      enddo
c$$$
c$$$c     Loop over cluster k-points.
c$$$      do kx = 0, lcx-1
c$$$        do ky = 0, lcy-1
c$$$          do kz = 0, lcz-1
c$$$
c$$$            k = kx + ky*lcx + kz*lcx*lcy
c$$$
c$$$            do ix = -ix_max, ix_max
c$$$              do iy = -iy_max, iy_max
c$$$                do iz = -iz_max, iz_max
c$$$
c$$$                  klx = mod(kx*x_stretch + ix + llx, llx)   
c$$$                  kly = mod(ky*y_stretch + iy + lly, lly)
c$$$                  klz = mod(kz*z_stretch + iz + llz, llz)
c$$$                  
c$$$                  kl = klx + kly*llx + klz*llx*lly
c$$$
c$$$            do l = 0, mp1
c$$$
c$$$              do i = 0, 4*nb-1
c$$$                do j = 0, 4*nb-1
c$$$                  if (i .eq. j) then
c$$$                    e_n(i,j) = epsilon(l)*cmplx(0.0d0, 1.0d0)
c$$$                  else
c$$$                    e_n(i,j) = epsilon(l)*cmplx(0.0d0, 0.0d0)
c$$$                  endif
c$$$                enddo
c$$$              enddo
c$$$
c$$$              temp_gc = cmplx(0.0d0, 0.0d0)
c$$$
c$$$c     For each cluster k-point, sum over all points in the
c$$$c     coarse-graining cell.
c$$$
c$$$              do ix = -ix_max, ix_max
c$$$                do iy = -iy_max, iy_max
c$$$                  do iz = -iz_max, iz_max
c$$$
c$$$                    klx = mod(kx*x_stretch + ix + llx, llx)   
c$$$                    kly = mod(ky*y_stretch + iy + lly, lly)
c$$$                    klz = mod(kz*z_stretch + iz + llz, llz)
c$$$
c$$$                    kl = klx + kly*llx + klz*llx*lly
c$$$
c$$$                    temp_gl2 = e_n - h0(:,:,kl) -
c$$$     $                 sigma(:,:,l,k)
c$$$                    temp_gl = identity
c$$$
c$$$c     Invert temp_gl
c$$$                    call cgesv(4*nb, 4*nb, temp_gl2, 4*nb, 
c$$$     $                 ipiv, temp_gl, 4*nb, info)
c$$$
c$$$                    if (info .ne. 0) then
c$$$                      write(6,*) 'info not equal to zero'
c$$$                    endif
c$$$
c$$$                    temp_gc = temp_gc +
c$$$     $                 wx(ix) * wy(iy) * wz(iz) * temp_gl
c$$$
c$$$c     Remove the analytic term so as to perform a trace
c$$$c     to get the tau=0+ lattice Green's function.
c$$$
c$$$                    do ia = 0, 1
c$$$                      do ib = 0, 1
c$$$                        temp_gl = temp_gl -
c$$$     $                     cl_k(ia,ib,:,:,kl) *
c$$$     $                     q_epsilon(ia,ib,l)
c$$$                      enddo
c$$$                    enddo
c$$$
c$$$c     Obtain g_tau0 by summing over all frequencies.  Include
c$$$c     weights wx,wy,wz as the lattice greens function may appear
c$$$c     in several course graining cells.
c$$$
c$$$                    g_ktau0(:,:,kl) =  g_ktau0(:,:,kl) + 
c$$$     $                 t * wx(ix) * wy(iy) * wz(iz) * 
c$$$     $                 temp_gl
c$$$
c$$$                  enddo
c$$$                enddo
c$$$              enddo
c$$$
c$$$              g(:,:,l,k) = temp_gc / float(coarse_grain_points)
c$$$
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$      enddo

      return
      end
          







