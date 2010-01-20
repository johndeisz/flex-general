#include "../convert.F"

      subroutine kinetic_energy(g_tau0, ek, kinetic)

#include "../constants.F"

      COMPLEX g_tau0(0:4*nb-1,0:4*nb-1,0:nl-1)
      COMPLEX ek(0:nb-1,0:nb-1,0:nl-1)
      COMPLEX delta(0:nb-1, 0:nb-1)
      REAL kinetic

      INTEGER k, nu1, nu2, is
      COMPLEX kin_temp

      delta = cmplx(0.0d0, 0.0d0)
      do nu1 = 0, nb-1
        delta(nu1, nu1) = cmplx(1.0d0, 0.0d0)
      enddo

      kin_temp = cmplx(0.0d0, 0.0d0)

      do nu1 = 0, nb-1
        do nu2 = 0, nb-1
          do is = 0, 1

            do k = 0, nl-1

              kin_temp = kin_temp + ek(nu1, nu2, k) * 
     $           ( g_tau0(4*nu2+is, 4*nu1+is, k) + delta(nu1,nu2) )

            enddo
            
          enddo
        enddo
      enddo

      write(6,*) 'kin_temp = ', kin_temp
            
      kinetic = real(kin_temp) / float(nl)

      return
      end
