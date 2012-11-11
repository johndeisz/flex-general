#include "../convert.F90"

subroutine analyze_psi_spin(psi, overall)

  USE CONSTANTS
  IMPLICIT NONE

  COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)
  REAL overall

  REAL norm
  COMPLEX overlap

  REAL chi_s(0:1, -1:1, 0:1, 0:1)

  INTEGER s, ms, is, isp, ib, ibp, ir

  write(6,*)
  write(6,*) 'Pair wave function weight by pair spin'

  chi_s = 0.0d0
  
  chi_s(0,0,0,0) = 0.0d0
  chi_s(0,0,0,1) = 1.0d0 / sqrt(2.0d0) 
  chi_s(0,0,1,0) = -1.0d0 / sqrt(2.d0)
  chi_s(0,0,1,1) = 0.0d0

  chi_s(1,1,0,0) = 1.0d0
  chi_s(1,1,0,1) = 0.0d0
  chi_s(1,1,1,0) = 0.0d0
  chi_s(1,1,1,1) = 0.0d0

  chi_s(1,0,0,0) = 0.0d0
  chi_s(1,0,0,1) = 1.0d0 / sqrt(2.0d0)
  chi_s(1,0,1,0) = 1.0d0 / sqrt(2.0d0)
  chi_s(1,0,1,1) = 0.0d0

  chi_s(1,-1, 0, 0) = 0.0d0
  chi_s(1,-1, 0, 1) = 0.0d0
  chi_s(1,-1, 1, 0) = 0.0d0
  chi_s(1,-1, 1, 1) = 1.0d0

  do s = 0, 1
     do ms = -s, s
        
        norm = 0.0d0

        do ib = 0, nb-1
           do ibp = 0, nb-1
              do ir = 0, nl - 1

                 overlap = 0.0d0

                 do is = 0, 1
                    do isp = 0, 1
                       
                       overlap = overlap +  &
                            chi_s(s, ms, is, isp)* &
                            psi(2*ib+is, 2*ibp+isp, ir)
                    
                    enddo
                 enddo

                 norm = norm + (cabs(overlap))**2.0d0 
                    
              enddo
           enddo
        enddo

        write(6,159) s,ms, norm/overall

     enddo
  enddo

159 format('|S=',i1,', Ms =',i2,'> = ',f13.10)
  
  return
end subroutine analyze_psi_spin
