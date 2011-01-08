#include "../convert.F"

      subroutine analyze_psi_exchange(psi, overall)

#include "../constants.F"

      COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)
      REAL overall

      COMPLEX overlap
      INTEGER ix, iy, iz, ind1, ind2
      INTEGER ib, is, ibp, isp
      REAL norm

      write(6,*)
      write(6,*) 'Pair wave function weight by exchange symmetry'

      overlap = cmplx(0.0d0, 0.0d0)
      norm = 0.0d0

      do ix = 0, llx1
        do iy = 0, lly1
          do iz = 0, llz1

            ind1 = ix + iy*llx + iz*llx*lly
            ind2 = mod(llx-ix,llx) +
     $         mod(lly-iy,lly)*llx +
     $         mod(llz-iz,llz)*llx*lly

            do ib = 0, nb-1
              do is = 0, 1
                do ibp = 0, nb-1
                  do isp = 0, 1

                    overlap = (1.0d0/2.0d0) *
     $                 (psi(2*ib+is, 2*ibp+isp, ind1) -
     $                 psi(2*ibp+isp, 2*ib+is, ind2) )

                    norm = norm + (cabs(overlap))**2.0d0
 
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
      enddo

      write(6,*) 'Antisymmetric weight = ', norm/overall

      norm = 0.0d0

      do ix = 0, llx1
        do iy = 0, lly1
          do iz = 0, llz1

            ind1 = ix + iy*llx + iz*llx*lly
            ind2 = mod(llx-ix,llx) +
     $         mod(lly-iy,lly)*llx +
     $         mod(llz-iz,llz)*llx*lly

            do ib = 0, nb-1
              do is = 0, 1
                do ibp = 0, nb-1
                  do isp = 0, 1
                    
                   overlap = (1.0d0/2.0d0) *
     $                 (psi(2*ib+is, 2*ibp+isp, ind1) +
     $                 psi(2*ibp+isp, 2*ib+is, ind2) )

                    norm = norm + (cabs(overlap))**2.0d0

                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
      enddo

      write(6,*) 'Symmetric weight = ', norm/overall

      return
      end
