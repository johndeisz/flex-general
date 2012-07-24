#include "../convert.F90"

subroutine kinetic_energy(g_tau0, tij, ed, kinetic)

  USE CONSTANTS
  USE bare_dispersion

  COMPLEX g_tau0(0:4*nb-1,0:4*nb-1,0:nl-1)
  REAL ed(0:nb-1)
  COMPLEX tij(0:nb-1,0:nb-1,-2:2,-2:2,-2:2)
  COMPLEX ek(0:nb-1,0:nb-1)
  COMPLEX delta(0:nb-1, 0:nb-1)
  REAL kinetic

  INTEGER k, nu1, nu2, is
  INTEGER klx, kly, klz, k_minus
  COMPLEX kin_temp

  delta = cmplx(0.0d0, 0.0d0)
  do nu1 = 0, nb-1
     delta(nu1, nu1) = cmplx(1.0d0, 0.0d0)
  enddo

  kin_temp = cmplx(0.0d0, 0.0d0)
  
  do k = 0, nl-1

     ek = ekl(k, tij, ed)

     do nu1 = 0, nb-1
        ek(nu1,nu1) = ek(nu1,nu1) - ed(nu1)  ! Remove local potential
     enddo

     do nu1 = 0, nb-1
        do nu2 = 0, nb-1

           do is = 0, 1

              kin_temp = kin_temp + ek(nu1, nu2) * &
                   ( g_tau0(4*nu2+is, 4*nu1+is, k) + delta(nu1,nu2) )

           enddo
           
        enddo
     enddo
  enddo

  write(6,*) 'Trace of hk in particle basis = ', kin_temp
            
  kinetic = real(kin_temp) / float(nl)

  ! Debug

  kin_temp = 0.0d0

  do klx = 0, llx-1
     do kly = 0, lly - 1
        do klz = 0, llz - 1

           k = klz*llx*lly + kly*llx + klx

           k_minus = mod(llz-klz,llz)*llx*lly + &
                mod(lly-kly,lly)*llx + mod(llx-klx,llx)

           ek = ekl(k, tij, ed)

           do nu1 = 0, nb-1
              ek(nu1,nu1) = ek(nu1,nu1) - ed(nu1)  ! Remove local potential
           enddo

           do nu1 = 0, nb-1
              do nu2 = 0, nb-1

                 do is = 0, 1

                    kin_temp = kin_temp -  ek(nu1, nu2) * &
                         g_tau0(4*nu1+2+is, 4*nu2+2+is, k_minus)

                 enddo
              enddo
           enddo

        enddo
     enddo
  enddo

  write(6,*) 'Trace of hk in hole basis = ', kin_temp

  return
end subroutine kinetic_energy
