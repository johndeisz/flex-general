#include "../convert.F90"

subroutine analyze_psi_1D(psi)

  USE CONSTANTS
  IMPLICIT NONE
  
  COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)

  REAL pi

  COMPLEX dot
  REAL dot_real, dot_imag
  REAL mag
  
  INTEGER ib

  pi = 4.0*atan(1.0d0)

  write(6,*)
  write(6,*) "One dimensional pair wave function amplitudes"

  !    local singlet
  write(6,*) "**** local s-wave ****"
  write(6,*) "ib  amp      phase/pi"
  do ib = 0, nb-1
     dot = sqrt(0.5d0)*(psi(2*ib+0,2*ib+1,0) - psi(2*ib+1,2*ib+0,0))       
     dot_real = real(dot)
     dot_imag = imag(dot)
     mag = sqrt( dot_real**2 + dot_imag**2)
     write(6,200) ib,  mag, atan2(dot_real, dot_imag) / pi
  enddo

  !     extended singlet 
  write(6,*) "**** extended s-wave *****"
  write(6,*) "ib  amp      phase/pi"
  do ib = 0, nb-1
     dot = 0.5d0 * ( psi(2*ib+0,2*ib+1,1) - psi(2*ib+1,2*ib+0,1) + &
          psi(2*ib+0,2*ib+1,llx1) - psi(2*ib+1,2*ib+0,llx1) )
     dot_real = real(dot)
     dot_imag = imag(dot)
     mag = sqrt( dot_real**2 + dot_imag**2)
     write(6,200) ib, mag, atan2(dot_real, dot_imag) / pi      
  enddo

  !     px, Sz=1
  write(6,*) "**** px, Sz=1  *****"
  write(6,*) "ib  amp      phase/pi"
  do ib = 0, nb-1
     dot = sqrt(0.5d0)*(psi(2*ib+0,2*ib+0,1) - psi(2*ib+0,2*ib+0,llx1) )
     dot_real = real(dot)
     dot_imag = imag(dot)
     mag = sqrt( dot_real**2 + dot_imag**2)
     write(6,200) ib, mag, atan2(dot_real, dot_imag) / pi
  enddo

  !     px, Sz=0
  write(6,*) "**** px, Sz=0  *****"
  write(6,*) "ib  amp      phase/pi"
  do ib = 0, nb-1
     dot = 0.5d0 * (psi(2*ib+0,2*ib+1,1) +  psi(2*ib+1,2*ib+0,1) - &
          psi(2*ib+0,2*ib+1,llx1) - psi(2*ib+1,2*ib+0,llx1))
     dot_real = real(dot)
     dot_imag = imag(dot)
     mag = sqrt( dot_real**2 + dot_imag**2)
     write(6,200) ib, mag, atan2(dot_real, dot_imag) / pi
  enddo

  !     px, Sz=-1
  write(6,*) "**** px, Sz=-1  *****"
  write(6,*) "ib  amp      phase/pi"
  do ib = 0, nb-1
     dot = sqrt(0.5d0) * ( psi(2*ib+1,2*ib+1,1) - psi(2*ib+1,2*ib+1,llx1) )
     dot_real = real(dot)
     dot_imag = imag(dot)
     mag = sqrt( dot_real**2 + dot_imag**2)
     write(6,200) ib, mag, atan2(dot_real, dot_imag) / pi
  enddo

200 format(i3, 2x, e16.8, 2x, e16.8)

  return
end subroutine analyze_psi_1D
