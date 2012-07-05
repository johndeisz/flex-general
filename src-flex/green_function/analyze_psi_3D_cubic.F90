#include "../convert.F90"

      subroutine analyze_psi_3D_cubic(psi)

#include "../constants.F90"

      COMPLEX psi(0:1, 0:1, 0:nlx1, 0:nly1, 0:nlz1)

      REAL pi

      COMPLEX dot
      REAL dot_real, dot_imag
      REAL mag

      pi = 4.0*atan(1.0d0)

      write(6,*)
      write(6,*) "Symmetry analysis of 3D pair wavefunction ",
     $   "using OH basis."

c     local singlet
      dot = sqrt(0.5d0)*(psi(0,1,0,0,0) - psi(1,0,0,0,0))       
      dot_real = real(dot)
      dot_imag = imag(dot)
      mag = sqrt( dot_real**2 + dot_imag**2)
      write(6,*) "s-wave, local: amp = ", mag,"  phase/pi = ",
     $     atan2(dot_real, dot_imag) / pi

c     extended s-wave (nn)
      dot = (1.0d0 / sqrt(12.0d0)) *
     $   ( psi(0,1,1,0,0) - psi(1,0,1,0,0) +
     $   psi(0,1,nlx1,0,0) - psi(1,0,nlx1,0,0) +
     $   psi(0,1,0,1,0) - psi(1,0,0,1,0) +
     $   psi(0,1,0,nly1,0) - psi(1,0,0,nly1,0) +
     $   psi(0,1,0,0,1) - psi(1,0,0,0,1) +
     $   psi(0,1,0,0,nlz1) - psi(1,0,0,0,nlz1) )
      dot_real = real(dot)
      dot_imag = imag(dot)
      mag = sqrt( dot_real**2 + dot_imag**2)
      write(6,*) "s-wave, x^2 + y^2 + z^2: amp = ", mag,
     $   "  phase/pi = ",
     $   atan2(dot_real, dot_imag) / pi      

c     d_{x2-y2}
      dot = 0.5d0 * sqrt(0.5d0) *
     $   ( psi(0,1,1,0,0) - psi(1,0,1,0,0) +
     $   psi(0,1,nlx1,0,0) - psi(1,0,nlx1,0,0) -
     $   psi(0,1,0,1,0) + psi(1,0,0,1,0) -
     $   psi(0,1,0,nly1,0) + psi(1,0,0,nly1,0) )
      dot_real = real(dot)
      dot_imag = imag(dot)
      mag = sqrt( dot_real**2 + dot_imag**2)
      write(6,*) "d-wave, x^2-y^2: amp = ", mag,
     $   "  phase/pi = ",
     $   atan2(dot_real, dot_imag) / pi      

c     d_{3z2-r2)
      dot = (1.0d0 / sqrt(24.0d0) ) *
     $   (2.0d0 * psi(0,1,0,0,1) - 2.0d0 * psi(1,0,0,0,1) +
     $   2.0d0 * psi(0,1,0,0,nlz1) - 2.0d0 * psi(1,0,0,0,nlz1) -
     $   psi(0,1,1,0,0) + psi(1,0,1,0,0) -
     $   psi(0,1,nlx1,0,0) + psi(1,0,nlx1,0,0) -
     $   psi(0,1,0,1,0) + psi(1,0,0,1,0) -
     $   psi(0,1,0,nly1,0) + psi(1,0,0,nly1,0) )
      dot_real = real(dot)
      dot_imag = imag(dot)
      mag = sqrt( dot_real**2 + dot_imag**2)
      write(6,*) "d-wave, 3z^2-r^2: amp = ", mag,
     $   "  phase/pi = ",
     $   atan2(dot_real, dot_imag) / pi      

c     d_xy
      dot = (1.0d0 / sqrt(8.0d0)) * 
     $   ( psi(0,1,1,1,0) - psi(1,0,1,1,0) -
     $   psi(0,1,1,nly1,0) + psi(1,0,1,nly1,0) -
     $   psi(0,1,nlx1,1,0) + psi(1,0,nlx1,1,0) +
     $   psi(0,1,nlx1,nly1,0) - psi(1,0,nlx1,nly1,0) )
      dot_real = real(dot)
      dot_imag = imag(dot)
      mag = sqrt( dot_real**2 + dot_imag**2)
      write(6,*) "d_{xy}: amp = ", mag,
     $   "  phase/pi = ",
     $   atan2(dot_real, dot_imag) / pi      

c     d_xz
      dot = (1.0d0 / sqrt(8.0d0)) * 
     $   ( psi(0,1,1,0,1) - psi(1,0,1,0,1) -
     $   psi(0,1,1,0,nlz1) + psi(1,0,1,0,nlz1) -
     $   psi(0,1,nlx1,0,1) + psi(1,0,nlx1,0,1) +
     $   psi(0,1,nlx1,0,nlz1) - psi(1,0,nlx1,0,nlz1) )
      dot_real = real(dot)
      dot_imag = imag(dot)
      mag = sqrt( dot_real**2 + dot_imag**2)
      write(6,*) "d_{xz}: amp = ", mag,
     $   "  phase/pi = ",
     $   atan2(dot_real, dot_imag) / pi       

c     d_yz
      dot = (1.0d0 / sqrt(8.0d0)) * 
     $   ( psi(0,1,0,1,1) - psi(1,0,0,1,1) -
     $   psi(0,1,0,1,nlz1) + psi(1,0,0,1,nlz1) -
     $   psi(0,1,0,nly1,1) + psi(1,0,0,nly1,1) +
     $   psi(0,1,0,nly1,nlz1) - psi(1,0,0,nly1,nlz1) )
      dot_real = real(dot)
      dot_imag = imag(dot)
      mag = sqrt( dot_real**2 + dot_imag**2)
      write(6,*) "d_{yz}: amp = ", mag,
     $   "  phase/pi = ",
     $   atan2(dot_real, dot_imag) / pi       

c     px
      dot = sqrt(0.5d0)*(psi(0,0,1,0,0) - psi(0,0,nlx1,0,0))
      dot_real = real(dot)
      dot_imag = imag(dot)
      mag = sqrt( dot_real**2 + dot_imag**2)
      write(6,*) "px, Sz=1: amp = ", mag,"  phase/pi = ",
     $   atan2(dot_real, dot_imag) / pi

      dot = 0.5d0 * (psi(0,1,1,0,0) + psi(1,0,1,0,0) -
     $   psi(0,1,nlx1,0,0) - psi(1,0,nlx1,0,0))
      dot_real = real(dot)
      dot_imag = imag(dot)
      mag = sqrt( dot_real**2 + dot_imag**2)
      write(6,*) "px, Sz=0: amp = ", mag,"  phase/pi = ",
     $   atan2(dot_real, dot_imag) / pi

      dot = sqrt(0.5d0) * ( psi(1,1,1,0,0) - 
     $   psi(1,1,nlx1,0,0) )
      dot_real = real(dot)
      dot_imag = imag(dot)
      mag = sqrt( dot_real**2 + dot_imag**2)
      write(6,*) "px, Sz=-1: amp = ", mag,"  phase/pi = ",
     $   atan2(dot_real, dot_imag) / pi

c     py
      dot = sqrt(0.5d0)*(psi(0,0,0,1,0) - psi(0,0,0,nly1,0))
      dot_real = real(dot)
      dot_imag = imag(dot)
      mag = sqrt( dot_real**2 + dot_imag**2)
      write(6,*) "py, Sz=1: amp = ", mag,"  phase/pi = ",
     $   atan2(dot_real, dot_imag) / pi
          
      dot = 0.5d0 * (psi(0,1,0,1,0) + psi(1,0,0,1,0) -
     $   psi(0,1,0,nly1,0) - psi(1,0,0,nly1,0))
      dot_real = real(dot)
      dot_imag = imag(dot)
      mag = sqrt( dot_real**2 + dot_imag**2)
      write(6,*) "py, Sz=0: amp = ", mag,"  phase/pi = ",
     $   atan2(dot_real, dot_imag) / pi

      dot = sqrt(0.5d0) * ( psi(1,1,0,1,0) - 
     $   psi(1,1,0,nly1,0) )
      dot_real = real(dot)
      dot_imag = imag(dot)
      mag = sqrt( dot_real**2 + dot_imag**2)
      write(6,*) "py, Sz=-1: amp = ", mag,"  phase/pi = ",
     $   atan2(dot_real, dot_imag) / pi

c     pz
      dot = sqrt(0.5d0) * ( psi(0,0,0,0,1) - psi(0,0,0,0,nlz1) )
      dot_real = real(dot)
      dot_imag = imag(dot)
      mag = sqrt( dot_real**2 + dot_imag**2)
      write(6,*) "pz, Sz=1: amp = ", mag,"  phase/pi = ",
     $   atan2(dot_real, dot_imag) / pi
          
      dot = 0.5d0 * (psi(0,1,0,0,1) + psi(1,0,0,0,1) -
     $   psi(0,1,0,0,nlz1) - psi(1,0,0,0,nlz1))
      dot_real = real(dot)
      dot_imag = imag(dot)
      mag = sqrt( dot_real**2 + dot_imag**2)
      write(6,*) "pz, Sz=0: amp = ", mag,"  phase/pi = ",
     $   atan2(dot_real, dot_imag) / pi

      dot = sqrt(0.5d0) * ( psi(1,1,0,0,1) - 
     $   psi(1,1,0,0,nlz1) )
      dot_real = real(dot)
      dot_imag = imag(dot)
      mag = sqrt( dot_real**2 + dot_imag**2)
      write(6,*) "pz, Sz=-1: amp = ", mag,"  phase/pi = ",
     $   atan2(dot_real, dot_imag) / pi

      return
      end
