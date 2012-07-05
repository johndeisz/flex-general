#include "../convert.F"

      subroutine bare_dispersion(tij, ed, ek, ek_min)

#include "../constants.F"

      COMPLEX, dimension (0:nb-1,0:nb-1,0:nl-1) :: tij, ek
      REAL, dimension (0:nb-1) :: ed
      REAL ek_min, pi

      INTEGER ll(1:3)

      INTEGER ib, ibp
      INTEGER il
      INTEGER isignv(0:3)

      COMPLEX, dimension (0:nl-1) :: tij_tmp

      pi = acos(-1.0d0)

      ll = (/ llx, lly, llz /)

! Perform a 3D Fourier transform to put the hopping matrix into k-space
! Include the on-site potential energy, ed

      isignv(0) = 0
      isignv(1) = -1
      isignv(2) = -1
      isignv(3) = -1

      do ib = 0, nb-1
        do ibp = 0, nb-1
          
          tij_tmp = tij(ib,ibp,:)

          if (ib .eq. ibp) then
            tij_tmp(0) = tij_tmp(0) + ed(ib)
          endif

          call fft_3D_lattice (tij_tmp, isignv)      
          ek(ib,ibp,:) = tij_tmp

        enddo
      enddo

! Determine the smallest value along the band-diagonal

      ek_min = 1.0d8

      do ib = 0, nb-1
        do il = 0, nl-1
          if ( real( ek(ib,ib,il) ) .lt. ek_min) then
            ek_min = real( ek(ib,ib,il) )
c$$$            write(6,*) ib, il, ek_min
          endif
        enddo
      enddo

      return
      end







