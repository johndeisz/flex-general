#include "../convert.F90"

subroutine analyze_psi_2D_ortho(psi)

  USE CONSTANTS
  IMPLICIT NONE

  COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)

  REAL pi

  COMPLEX dot
  REAL dot_real, dot_imag
  REAL mag

  INTEGER ib, jb

  pi = 4.0*atan(1.0d0)

  write(6,*)
  write(6,*) "Symmetry analysis of 2D pair wavefunction using D2H basis."
  write(6,*) 

  !     local singlet
  write(6,*) "**** local s-wave ****"
  write(6,*) "ib  jb  amp      phase/pi"
  do ib = 0, nb-1
     do jb = 0, nb-1
        dot = psi(2*ib+0,2*jb+1,0)- psi(2*ib+1,2*jb+0,0)
        if (ib .eq. jb) then
           dot = dot / sqrt(2.0d0)
        endif
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, jb, mag, atan2(dot_real, dot_imag) / pi
     enddo
  enddo

  !     singlet, x^2
  write(6,*) "**** singlet x^2 ****"
  write(6,*) "ib  jb  amp      phase/pi"
  do ib = 0, nb-1
     do jb = 0, nb-1
        dot = psi(2*ib+0,2*jb+1,1) - psi(2*ib+1,2*jb+0,1) 
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, jb, mag, atan2(dot_real, dot_imag) / pi      
     enddo
  enddo

  !     extended s, 2x^2
  write(6,*) "**** extended s, (2x)^2 ****"
  write(6,*) "ib  jb  amp      phase/pi"
  do ib = 0, nb-1
     do jb = 0, nb-1
        dot = psi(2*ib+0,2*jb+1,2) - psi(2*ib+1,2*jb+0,2) 
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, jb, mag, atan2(dot_real, dot_imag) / pi      
     enddo
  enddo

  !     extended s, y^2
  write(6,*) "**** extended s, y^2 ****"
  write(6,*) "ib  jb   amp      phase/pi"
  do ib = 0, nb-1
     do jb = 0, nb-1
        dot = psi(2*ib+0,2*jb+1,llx) - psi(2*ib+1,2*jb+0,llx) 
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, jb, mag, atan2(dot_real, dot_imag) / pi      
     enddo
  enddo

  !     extended s, 2y^2
  write(6,*) "**** extended s, (2y)^2 ****"
  write(6,*) "ib  jb amp      phase/pi"
  do ib = 0, nb-1
     do jb = 0, nb-1
        dot = psi(2*ib+0,2*jb+1,2*llx) - psi(2*ib+1,2*jb+0,2*llx)
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, jb, mag, atan2(dot_real, dot_imag) / pi      
     enddo
  enddo

  !     extended s, x+y (A2)
  write(6,*) "**** extended s, x+y (A2) ****"
  write(6,*) "ib  jb amp      phase/pi"
  do ib = 0, nb-1
     do jb = 0, nb-1
        dot = psi(2*ib+0,2*jb+1,llx+1) - psi(2*ib+1,2*jb+0,llx+1) 
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, jb, mag, atan2(dot_real, dot_imag) / pi      
     enddo
  enddo

  !     extended s, x-y (A2)
  write(6,*) "**** extended s, x-y (A2) ****"
  write(6,*) "ib  jb amp      phase/pi"
  do ib = 0, nb-1
     do jb = 0, nb-1
        dot = psi(2*ib+0,2*jb+1,llx+llx1) - psi(2*ib+1,2*jb+0,llx+llx1) 
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, jb, mag, atan2(dot_real, dot_imag) / pi      
     enddo
  enddo

  !     extended s, 2x+2y (A2)
  write(6,*) "**** extended s, 2x+2y (A2) ****"
  write(6,*) "ib  jb amp      phase/pi"
  do ib = 0, nb-1
     do jb = 0, nb-1
        dot = psi(2*ib+0,2*jb+1,2*llx+2) - psi(2*ib+1,2*jb+0,2*llx+2) 
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, jb, mag, atan2(dot_real, dot_imag) / pi      
     enddo
  enddo

  !     extended s, 2x-2y (A2)
  write(6,*) "**** extended s, 2x-2y (A2) ****"
  write(6,*) "ib  jb amp      phase/pi"
  do ib = 0, nb-1
     do jb = 0, nb-1
        dot = psi(2*ib+0,2*jb+1,2*llx+llx1-1) - &
             psi(2*ib+1,2*jb+0,2*llx+llx1-1) 
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, jb, mag, atan2(dot_real, dot_imag) / pi      
     enddo
  enddo
  
  write(6,*) "**** local p, Sz=1 ****"
  write(6,*) "ib  jb amp      phase/pi"
  do ib = 0, nb-1
     do jb = 0, nb-1
        dot = sqrt(2.0d0)*psi(2*ib+0,2*jb+0,0) 
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, jb, mag, atan2(dot_real, dot_imag) / pi
     enddo
  enddo

  write(6,*) "**** local p, Sz=0 ****"
  write(6,*) "ib  jb amp      phase/pi"
  do ib = 0, nb-1
     do jb = 0, nb-1
        dot = psi(2*ib+0,2*jb+1,0)+ psi(2*ib+1,2*jb+0,0)
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, jb, mag, atan2(dot_real, dot_imag) / pi
     enddo
  enddo

  write(6,*) "**** local p, Sz=-1 ****"
  write(6,*) "ib  jb amp      phase/pi"
  do ib = 0, nb-1
     do jb = 0, nb-1
        dot = sqrt(2.0d0)*psi(2*ib+1,2*jb+1,1) 
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, jb, mag, atan2(dot_real, dot_imag) / pi
     enddo
  enddo

  !     px
  write(6,*) "**** px, Sz=1 ****"
  write(6,*) "ib  jb amp      phase/pi"
  do ib = 0, nb-1
     do jb = 0, nb-1
        dot = sqrt(2.0d0)*psi(2*ib+0,2*jb+0,1) 
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, jb, mag, atan2(dot_real, dot_imag) / pi
     enddo
  enddo

  write(6,*) "**** px, Sz=0 ****"
  write(6,*) "ib jb  amp      phase/pi"
  do ib = 0, nb-1
     do jb = 0, nb-1
        dot = psi(2*ib+0,2*jb+1,1) + psi(2*ib+1,2*jb+0,1) 
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, jb, mag, atan2(dot_real, dot_imag) / pi
     enddo
  enddo

  write(6,*) "**** px, Sz=-1 ****"
  write(6,*) "ib jb  amp      phase/pi"
  do ib = 0, nb-1
     do jb = 0, nb-1
        dot = sqrt(2.0d0) * psi(2*ib+1,2*jb+1,1)  
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, jb, mag, atan2(dot_real, dot_imag) / pi
     enddo
  enddo
      
  write(6,*) "**** py, Sz=1 ****"
  write(6,*) "ib  jb amp      phase/pi"
  do ib = 0, nb-1
     do jb = 0, nb-1
        dot = sqrt(2.0d0)*psi(2*ib+0,2*jb+0,llx) 
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, jb, mag, atan2(dot_real, dot_imag) / pi
     enddo
  enddo

  write(6,*) "**** py, Sz=0 ****"
  write(6,*) "ib jb  amp      phase/pi"
  do ib = 0, nb-1
     do jb = 0, nb-1
        dot = psi(2*ib+0,2*jb+1,llx) + psi(2*ib+1,2*jb+0,llx) 
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, jb, mag, atan2(dot_real, dot_imag) / pi
     enddo
  enddo

  write(6,*) "**** py, Sz=-1 ****"
  write(6,*) "ib jb  amp      phase/pi"
  do ib = 0, nb-1
     do jb = 0, nb-1
        dot = sqrt(2.0d0) * psi(2*ib+1,2*jb+1,llx)  
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, jb, mag, atan2(dot_real, dot_imag) / pi
     enddo
  enddo

  write(6,*) "**** px+y, Sz=1 ****"
  write(6,*) "ib  jb amp      phase/pi"
  do ib = 0, nb-1
     do jb = 0, nb-1
        dot = sqrt(2.0d0)*psi(2*ib+0,2*jb+0,llx+1) 
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, jb, mag, atan2(dot_real, dot_imag) / pi
     enddo
  enddo

  write(6,*) "**** px+y, Sz=0 ****"
  write(6,*) "ib jb  amp      phase/pi"
  do ib = 0, nb-1
     do jb = 0, nb-1
        dot = psi(2*ib+0,2*jb+1,llx+1) + psi(2*ib+1,2*jb+0,llx+1) 
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, jb, mag, atan2(dot_real, dot_imag) / pi
     enddo
  enddo

  write(6,*) "**** px+y, Sz=-1 ****"
  write(6,*) "ib jb  amp      phase/pi"
  do ib = 0, nb-1
     do jb = 0, nb-1
        dot = sqrt(2.0d0) * psi(2*ib+1,2*jb+1,llx+1)  
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, jb, mag, atan2(dot_real, dot_imag) / pi
     enddo
  enddo

  write(6,*) "**** px-y, Sz=1 ****"
  write(6,*) "ib  jb amp      phase/pi"
  do ib = 0, nb-1
     do jb = 0, nb-1
        dot = sqrt(2.0d0)*psi(2*ib+0,2*jb+0,llx+llx1) 
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, jb, mag, atan2(dot_real, dot_imag) / pi
     enddo
  enddo
  
  write(6,*) "**** px-y, Sz=0 ****"
  write(6,*) "ib jb  amp      phase/pi"
  do ib = 0, nb-1
     do jb = 0, nb-1
        dot = psi(2*ib+0,2*jb+1,llx+llx1) +  psi(2*ib+1,2*jb+0,llx+llx1) 
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, jb, mag, atan2(dot_real, dot_imag) / pi
     enddo
  enddo
  
  write(6,*) "**** px-y, Sz=-1 ****"
  write(6,*) "ib jb  amp      phase/pi"
  do ib = 0, nb-1
     do jb = 0, nb-1
        dot = sqrt(2.0d0) * psi(2*ib+1,2*jb+1,llx+llx1)  
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, jb, mag, atan2(dot_real, dot_imag) / pi
     enddo
  enddo

200 format(i3, 2x, i3, 2x, e16.8, 2x, e16.8)
  
  return
end subroutine analyze_psi_2D_ortho
