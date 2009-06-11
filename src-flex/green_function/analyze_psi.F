#include "../convert.F"

      subroutine analyze_psi(psi, tij)

#include "../constants.F"

      COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)
      COMPLEX tij(0:nb-1,0:nb-1,0:nl-1)

      INTEGER isignv(0:3)
      INTEGER i, j, ix, iy, iz, jx, jy, jz
      COMPLEX total
      REAL pi

      COMPLEX dot
      REAL dot_real, dot_imag
      REAL mag

      pi = 4.0*atan(1.0d0)


c     First fft the wave function to r-space
      isignv(0) = 0
      isignv(1) = 1
      isignv(2) = 1
      isignv(3) = 1

      call psi_transform(psi, isignv)

      psi = psi / float(nl)

      if ( (lly .eq. 1) .and. (llz .eq. 1) ) then  ! 1D
        call analyze_psi_1D(psi)
      endif


      if ( (lly .gt. 1) .and. (llz .eq. 1) ) then ! 2D

        if ( (lly .eq. llx) .and. (lcx .eq. lcy) .and.
     $     ( cabs(tij(0,0,1) - tij(0,0,llx)) .lt. 0.00001d0) ) then

          call analyze_psi_2D_tet(psi)
        
        else

          call analyze_psi_2D_ortho(psi)

        endif

      endif
c$$$
c$$$
c$$$      if ( (nly .gt. 1) .and. (nlz .gt. 1) ) then ! 3D
c$$$ 
c$$$        if ( (nly .eq. nlx) .and. (nx .eq. ny) .and.
c$$$     $     (abs(tx - ty) .lt. 0.00001d0) .and.
c$$$     $     (abs(txx - tyy) .lt. 0.00001d0) ) then !  x = y
c$$$
c$$$          if ( (nlz .eq. nlx) .and. (nz .eq. nx) .and.
c$$$     $       (abs(tx - tz) .lt. 0.00001d0) .and.
c$$$     $       (abs(txx - tzz) .lt. 0.00001d0) ) then ! cubic
c$$$
c$$$            call analyze_psi_3D_cubic(psi)
c$$$
c$$$          else                  ! tetragonal
c$$$
c$$$            call analyze_psi_3D_tet(psi)
c$$$
c$$$          endif
c$$$
c$$$        endif
c$$$
c$$$      endif
c$$$
c$$$
c$$$      write(6,*) 
c$$$      write(6,*) "Pair wave function "
c$$$      write(6,*) "i  j   ix   iy    iz    psi(i,j, ix, iy, iz)"
c$$$
c$$$      do i = 0, 1
c$$$        do j = 0, 1
c$$$          do ix = 0, nlx1
c$$$
c$$$            if (ix .gt. nlx/2) then
c$$$              jx = ix - nlx
c$$$            else
c$$$              jx = ix
c$$$            endif
c$$$
c$$$            do iy = 0, nly1
c$$$
c$$$              if (iy .gt. nly/2) then
c$$$                jy = iy - nly
c$$$              else
c$$$                jy = iy
c$$$              endif 
c$$$
c$$$              do iz = 0, nlz1
c$$$
c$$$                if (iz .gt. nlz/2) then
c$$$                  jz = iz - nlz
c$$$                else
c$$$                  jz = iz
c$$$                endif
c$$$
c$$$                if ( (ABS(jz) .le. 4) .and. (ABS(jy) .le. 4)
c$$$     $             .and. (ABS(jx) .le. 4) ) then
c$$$
c$$$                  write(6,*) i, " ", j, "  ", jx, " ", jy, " ",
c$$$     $              jz, "  ", psi(i,j,jx,jy,jz)
c$$$
c$$$                 endif
c$$$
c$$$              enddo
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$      enddo

      return
      end
