#include "../convert.F90"

subroutine analyze_psi_exchange(psi, overall)

  USE CONSTANTS
  IMPLICIT NONE

  COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)
  REAL overall

  COMPLEX overlap
  INTEGER ix, iy, iz, ind1, ind2
  INTEGER ib, is, ibp, isp
  REAL norm

  write(6,*)
  write(6,*) 'Pair wave function weight by exchange symmetry'

  norm = 0.0d0

  do ix = 0, llx1
     do iy = 0, lly1
        do iz = 0, llz1

           ind1 = ix + iy*llx + iz*llx*lly
           ind2 = mod(llx-ix,llx) + mod(lly-iy,lly)*llx + &
                mod(llz-iz,llz)*llx*lly

           do ib = 0, nb-1
              do is = 0, 1
                 do ibp = 0, nb-1
                    do isp = 0, 1

                       overlap = (1.0d0/2.0d0) * &
                            (psi(2*ib+is, 2*ibp+isp, ind1) - &
                            psi(2*ibp+isp, 2*ib+is, ind2) )

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
           ind2 = mod(llx-ix,llx) + mod(lly-iy,lly)*llx + &
                mod(llz-iz,llz)*llx*lly

           do ib = 0, nb-1
              do is = 0, 1
                 do ibp = 0, nb-1
                    do isp = 0, 1
                    
                       overlap = (1.0d0/2.0d0) * &
                            (psi(2*ib+is, 2*ibp+isp, ind1) + &
                            psi(2*ibp+isp, 2*ib+is, ind2) ) 

                            norm = norm + (cabs(overlap))**2.0d0

                    enddo
                 enddo
              enddo
           enddo
           
        enddo
     enddo
  enddo

  write(6,*) 'Symmetric weight = ', norm/overall

!!$c     write(6,*) 'ib is ibp isp ind1 ind2'
!!$c     write(6,*) 'psi'
!!$c     write(6,*) 'P psi'
!!$
!!$c     do ix = 0, llx1
!!$c       do iy = 0, lly1
!!$c         do iz = 0, llz1
!!$
!!$c           ind1 = ix + iy*llx + iz*llx*lly
!!$c           ind2 = mod(llx-ix,llx) +
!!$c    $         mod(lly-iy,lly)*llx +
!!$c    $         mod(llz-iz,llz)*llx*lly
!!$
!!$c           do ib = 0, nb-1
!!$c             do is = 0, 1
!!$c               do ibp = 0, nb-1
!!$c                 do isp = 0, 1
!!$
!!$
!!$c                  write(6,58) ib, is, ibp, isp, ind1, ind2
!!$c                   write(6,*)  psi(2*ib+is, 2*ibp+isp, ind1) 
!!$c                   write(6,*)  psi(2*ibp+isp, 2*ib+is, ind2)
!!$
!!$c                 enddo
!!$c                enddo
!!$c               enddo
!!$c              enddo
!!$
!!$c             enddo
!!$c            enddo
!!$c           enddo

58 format(i1,2x,i1,2x,i1,2x,i1,2x,i4,2x,i4)
      

  return
end subroutine analyze_psi_exchange
