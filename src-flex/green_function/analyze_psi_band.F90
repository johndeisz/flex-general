#include "../convert.F90"

subroutine analyze_psi_band(psi, overall)

  USE CONSTANTS

  COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)
  REAL overall

  REAL norm
  COMPLEX overlap

  INTEGER is, isp, ib, ibp, ir

  write(6,*)
  write(6,*) 'Pair wave function weight by band indices'

  do ib = 0, nb-1
     do ibp = 0, nb-1

        norm = 0.0d0

        do is = 0, 1
           do isp = 0, 1
              do ir = 0, nl - 1

                 overlap = psi(2*ib+is, 2*ibp+isp, ir)
                 norm = norm + (cabs(overlap))**2.0d0 
                    
              enddo
           enddo
        enddo

        write(6,159) ib, ibp, norm/overall

     enddo
  enddo

159 format('|ib=',i2,', ibp =',i2,'> = ',f10.8)

  return
end subroutine analyze_psi_band
