      subroutine bare_dispersion(tij, ed, ek)

#include "constants.F"

      DOUBLE COMPLEX, dimension (0:nb-1,0:nb-1,0:nl-1) :: tij, ek
      DOUBLE PRECISION, dimension (0:nb-1) :: ed
      DOUBLE PRECISION pi

      INTEGER ll(1:3)

      INTEGER ib, ibp
      INTEGER il
      INTEGER isignv(0:3)

      DOUBLE COMPLEX, dimension (0:nl-1) :: tij_tmp

      pi = dacos(-1.0d0)

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

      return
      end






