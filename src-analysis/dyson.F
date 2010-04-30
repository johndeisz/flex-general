      subroutine dyson(h0k, klin, w_min, dw, n_w, eta)

#include "constants.F"

      DOUBLE COMPLEX, dimension (0:4*nb-1,0:4*nb-1) :: h0k
      INTEGER klin
      DOUBLE PRECISION w_min, dw, eta
      INTEGER n_w

      DOUBLE COMPLEX, dimension (0:4*nb-1,0:4*nb-1, 0:nc-1, 
     $   0:n_pade_max) :: p, q

      INTEGER sigma_flag
      INTEGER n_pade(0:4*nb-1,0:4*nb-1,0:nc-1)
      INTEGER nu1, ia1, nu2, ia2, na1, na2

      INTEGER id1, id2, id3, id4, id5

      DOUBLE COMPLEX, dimension (0:4*nb-1,0:4*nb-1) :: sigma
      DOUBLE COMPLEX pade_evaluate
      EXTERNAL pade_evaluate
 
      INTEGER x_stretch, y_stretch, z_stretch
      DOUBLE PRECISION wx(-llx/lcx: llx/lcx)
      DOUBLE PRECISION wy(-lly/lcy: lly/lcy)
      DOUBLE PRECISION wz(-llz/lcz: llz/lcz)
      INTEGER ix_max, iy_max, iz_max
      DOUBLE COMPLEX e_n(0:4*nb-1, 0:4*nb-1)
      DOUBLE PRECISION w
      DOUBLE COMPLEX z

c     Counters
      INTEGER ix, iy, iz
      INTEGER i, j
      INTEGER klx, kly, klz, kl
      INTEGER kx, ky, kz, k, l
      INTEGER k1, k2, k3
      INTEGER ia, ib

c     Arrays
      DOUBLE COMPLEX, dimension (0:4*nb-1,0:4*nb-1) ::
     $   a, identity, temp_gl2, temp_gl

c     Lapack related
      INTEGER ipiv(4*nb)
      INTEGER info

      DOUBLE PRECISION sum

      p = dcmplx(0.0d0, 0.0d0)
      q = dcmplx(0.0d0, 0.0d0)

      write(6,*) "Dynamical self-energy? Yes-1, No-0."
      read(5,*) sigma_flag
      
      if (sigma_flag .eq. 1) then
        read(50,*) 
        read(50,*)
c     Read pade coefficients
        do nu1 = 0, nb-1
          do ia1 = 0, 3
            na1 = 4*nu1 + ia1
            do nu2 = 0, nb-1
              do ia2 = 0, 3
                na2 = 4*nu2 + ia2
                do k = 0, nc-1
                  read(50,*) id1, id2, id3, id4, id5, 
     $               n_pade(na1,na2,k)
                  if (n_pade(na1,na2,k) .gt. 0) then
                    do l = 0, n_pade(na1,na2,k)
                      read(50,*) p(na1, na2, k, l),
     $                   q(na1,na2,k,l)
                    enddo
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      endif

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

      do l = 1, n_w

        w = w_min + dw*dfloat(l-1)
        z = dcmplx(w, eta)
        a = dcmplx(0.0d0,0.0d0)

c     Loop over cluster k-points.
      
        do kx = 0, lcx-1
          do ky = 0, lcy-1
            do kz = 0, lcz-1

              k = kx + ky*lcx + kz*lcx*lcy

              do ix = -ix_max, ix_max
                do iy = -iy_max, iy_max
                  do iz = -iz_max, iz_max

                    klx = mod(kx*x_stretch + ix + llx, llx)   
                    kly = mod(ky*y_stretch + iy + lly, lly)
                    klz = mod(kz*z_stretch + iz + llz, llz)
                  
                    kl = klx + kly*llx + klz*llx*lly
                    
                    if (kl .eq. klin) then
                      
                      do i = 0, 4*nb-1
                        do j = 0, 4*nb-1

                          sigma(i,j) = pade_evaluate(z,n_pade(i,j,k),
     $                       p(i,j,k,:), q(i,j,k,:))

                          if (i .eq. j) then
                            e_n(i,j) = z
                          else
                            e_n(i,j) = dcmplx(0.0d0, 0.0d0)
                          endif

                        enddo
                      enddo

                      temp_gl2 = e_n - h0k - sigma
                      temp_gl = identity

c     Invert temp_gl2

                      call cgesv(4*nb, 4*nb, temp_gl2, 4*nb, 
     $                   ipiv, temp_gl, 4*nb, info)

                      if (info .ne. 0) then
                        write(6,*) 'info not equal to zero'
                      endif

                      a = a + wx(ix)*wy(iy)*wz(iz)*temp_gl

                    endif

                  enddo
                enddo
              enddo

            enddo
          enddo
        enddo

        sum = 0.0d0

        do nu1 = 0, nb-1
          do ia1 = 0, 1

            na1 = 4*nu1 + ia1
            sum = sum + dimag(a(na1,na1))

          enddo
        enddo

        write(6,*) w, -sum / (2.0d0*pi)

      enddo

      return
      end
          







