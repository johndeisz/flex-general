#include "../convert.F90"

subroutine pair_wave_test(psi, g_tau0, alpha, m_psi)

  USE CONSTANTS
  IMPLICIT NONE

  COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)
  COMPLEX g_tau0(0:4*nb-1, 0:4*nb-1, 0:nl-1)
  REAL m_band, m_psi, alpha

  COMPLEX psi_old(0:2*nb-1,0:2*nb-1, 0:nl-1)
  
  INTEGER k, ib, ibp, is, isp

  REAL phase_old, phase_new

  INTEGER klx, kly, klz, k_minus

  psi_old = psi
  phase_old = atan2(imag(psi(0,1,llx+2)),real(psi(0,1,llx+2)))

  m_psi = 0.0d0
  
  do ib = 0, nb-1

     m_band = 0.0d0         

     do ibp = 0, nb-1

        do is = 0, 1
           do isp = 0, 1

              do k = 0, nl-1

                 m_band = m_band + ( cabs( g_tau0( 4*ibp+2+isp, &
                      4*ib+0+is,k)) )**2

              enddo
                
           enddo
        enddo

     enddo

     if (m_band .gt. m_psi) then
        m_psi = m_band
     endif

  enddo

  m_psi = m_psi / float(nl)
  m_psi = sqrt(m_psi)

  if (m_psi .gt. 1.0d-10) then
    
     do ib = 0, nb-1
        do ibp = 0, nb-1

           do is = 0, 1
              do isp = 0, 1

                 do klx = 0, llx -1
                    do kly = 0, lly-1
                       do klz = 0, llz - 1

                          k = klx + kly*llx + klz*llx*lly

                          k_minus = mod(llz-klz,llz)*llx*lly + &
                               mod(lly-kly,lly)*llx + mod(llx-klx,llx)

                          psi(2*ib+is,2*ibp+isp,k_minus) = &
                               -g_tau0(4*ibp+2+isp, 4*ib+0+is,k) / m_psi

                       enddo
                    enddo
                 enddo
                
              enddo
           enddo
            
        enddo
     enddo

  else

     psi = cmplx(0.0d0,0.0d0)
 
  endif
        
  !    psi = (1.0d0 - alpha)*psi_old + alpha*psi

  !     phase_new = atan2(imag(psi(0,1,llx+2)),real(psi(0,1,llx+2)))
  !   psi = psi*cexp(cmplx(0.0d0,phase_old-phase_new))

  return
end subroutine pair_wave_test



