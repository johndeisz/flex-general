#include "../convert.F"

      subroutine analyze_psi_2D_tet(psi)

#include "../constants.F"

      COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)

      REAL pi

      COMPLEX dot
      REAL dot_real, dot_imag
      REAL mag

      INTEGER ib

      pi = 4.0*atan(1.0d0)

      write(6,*)
      write(6,*) "Symmetry analysis of 2D pair wavefunction ",
     $   "using D4H basis."

c     local singlet
      write(6,*) "**** local s-wave ****"
      write(6,*) "ib  amp      phase/pi"
      do ib = 0, nb-1
        dot = sqrt(0.5d0)*(psi(2*ib+0,2*ib+1,0) - 
     $     psi(2*ib+1,2*ib+0,0))       
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, mag, atan2(dot_real, dot_imag) / pi
      enddo

c     extended s-wave (nn)
      write(6,*) "**** extended s-wave (nn) ****"
      write(6,*) "ib  amp      phase/pi"
      do ib = 0, nb-1     
        dot = 0.5d0 * sqrt(0.5d0) *
     $     ( psi(2*ib+0,2*ib+1,1) - 
     $     psi(2*ib+1,2*ib+0,1) +
     $     psi(2*ib+0,2*ib+1,llx1) - 
     $     psi(2*ib+1,2*ib+0,llx1) +
     $     psi(2*ib+0,2*ib+1,llx) - 
     $     psi(2*ib+1,2*ib+0,llx) +
     $     psi(2*ib+0,2*ib+1,lly1*llx) - 
     $     psi(2*ib+1,2*ib+0,lly1*llx) )
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, mag, atan2(dot_real, dot_imag) / pi      
      enddo

c     d_{x2-y2}
      write(6,*) "**** d_{x2-y2} ****"
      write(6,*) "ib  amp      phase/pi"
      do ib = 0, nb-1 
        dot = 0.5d0 * sqrt(0.5d0) *
     $     ( psi(2*ib+0,2*ib+1,1) - 
     $     psi(2*ib+1,2*ib+0,1) +
     $     psi(2*ib+0,2*ib+1,llx1) - 
     $     psi(2*ib+1,2*ib+0,llx1) -
     $     psi(2*ib+0,2*ib+1,llx) + 
     $     psi(2*ib+1,2*ib+0,llx) -
     $     psi(2*ib+0,2*ib+1,lly1*llx) + 
     $     psi(2*ib+1,2*ib+0,lly1*llx) )
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, mag, atan2(dot_real, dot_imag) / pi      
      enddo
      
c     d_xy
      write(6,*) "**** d_xy ****"
      write(6,*) "ib  amp      phase/pi"
      do ib = 0, nb-1 
        dot = 0.5d0 * sqrt(0.5d0) * 
     $     ( psi(2*ib+0,2*ib+1,llx+1) - 
     $     psi(2*ib+1,2*ib+0,llx+1) -
     $     psi(2*ib+0,2*ib+1,lly1*llx+1) + 
     $     psi(2*ib+1,2*ib+0,lly1*llx+1) -
     $     psi(2*ib+0,2*ib+1,llx1+llx) + 
     $     psi(2*ib+1,2*ib+0,llx1+llx) +
     $     psi(2*ib+0,2*ib+1,llx1+lly1*llx) - 
     $     psi(2*ib+1,2*ib+0,llx1+lly1*llx) )
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, mag, atan2(dot_real, dot_imag) / pi      
      enddo

c     px
      write(6,*) "**** px, Sz=1 ****"
      write(6,*) "ib  amp      phase/pi"
      do ib = 0, nb-1
        dot = sqrt(0.5d0)*(psi(2*ib+0,2*ib+0,1) - 
     $     psi(2*ib+0,2*ib+0,llx1))
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, mag, atan2(dot_real, dot_imag) / pi
      enddo

      write(6,*) "**** px, Sz=0 ****"
      write(6,*) "ib  amp      phase/pi"
      do ib = 0, nb-1
        dot = 0.5d0 * (psi(2*ib+0,2*ib+1,1) + 
     $     psi(2*ib+1,2*ib+0,1) -
     $     psi(2*ib+0,2*ib+1,llx1) - 
     $     psi(2*ib+1,2*ib+0,llx1))
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, mag, atan2(dot_real, dot_imag) / pi
      enddo

      write(6,*) "**** px, Sz=-1 ****"
      write(6,*) "ib  amp      phase/pi"
      do ib = 0, nb-1
        dot = sqrt(0.5d0) * ( psi(2*ib+1,2*ib+1,1) - 
     $     psi(2*ib+1,2*ib+1,llx1) )
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, mag, atan2(dot_real, dot_imag) / pi
      enddo

c     py
      write(6,*) "**** py, Sz=1 ****"
      write(6,*) "ib  amp      phase/pi"
      do ib = 0, nb-1
        dot = sqrt(0.5d0)*(psi(2*ib+0,2*ib+0,llx) - 
     $     psi(2*ib+0,2*ib+0,lly1*llx))
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, mag, atan2(dot_real, dot_imag) / pi
      enddo

      write(6,*) "**** py, Sz=0 ****"
      write(6,*) "ib  amp      phase/pi"
      do ib = 0, nb-1
        dot = 0.5d0 * (psi(2*ib+0,2*ib+1,llx) + 
     $     psi(2*ib+1,2*ib+0,llx) -
     $     psi(2*ib+0,2*ib+1,lly1*llx) - 
     $     psi(2*ib+1,2*ib+0,lly1*llx))
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, mag, atan2(dot_real, dot_imag) / pi
      enddo

      write(6,*) "**** py, Sz=-1 ****"
      write(6,*) "ib  amp      phase/pi"
      do ib = 0, nb-1
        dot = sqrt(0.5d0) * ( psi(2*ib+1,2*ib+1,llx) - 
     $     psi(2*ib+1,2*ib+1,lly1*llx) )
        dot_real = real(dot)
        dot_imag = imag(dot)
        mag = sqrt( dot_real**2 + dot_imag**2)
        write(6,200) ib, mag, atan2(dot_real, dot_imag) / pi
      enddo

 200  format(i3, 2x, e16.8, 2x, e16.8)

      return
      end
