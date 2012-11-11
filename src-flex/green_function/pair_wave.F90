#include "../convert.F90"

subroutine pair_wave(psi, g_tau0, alpha, m_psi)

  USE CONSTANTS
  IMPLICIT NONE

#ifdef USE_MPI
  include 'mpif.h'
#endif

  COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)
  COMPLEX g_tau0(0:4*nb-1, 0:4*nb-1, 0:nl-1)
  REAL m_band, m_psi, alpha

  COMPLEX psi_old(0:2*nb-1,0:2*nb-1, 0:nl-1)
  COMPLEX tmp_psi

  INTEGER k, ib, ibp, is, isp, jb, js

  INTEGER rank
  REAL phase_old, phase_new

#ifdef USE_MPI
  INTEGER ierr
#endif /* USE_MPI */

#ifdef USE_MPI
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
#else
  rank = 0
#endif /* USE_MPI */

  psi_old = psi
  phase_old = atan2(imag(psi(0,1,llx+2)),real(psi(0,1,llx+2)))

  if (rank .eq. 0) then

     m_psi = 0.0d0

     do ib = 0, nb-1

        m_band = 0.0d0         

        do ibp = 0, nb-1

           do is = 0, 1
              do isp = 0, 1

                 do k = 0, nl-1

                    m_band = m_band + & 
                         ( cabs( g_tau0( 4*ib+2+is, 4*ibp+0+isp,k)) )**2

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

                    do k = 0, nl-1

                       psi(2*ib+is,2*ibp+isp,k) =  &
                            g_tau0( 4*ib+2+is, 4*ibp+0+isp,k) / m_psi

                    enddo
                
                 enddo
              enddo
            
           enddo
        enddo

     else

        psi = cmplx(0.0d0,0.0d0)
 
     endif
        
  endif

  psi = (1.0d0 - alpha)*psi_old + alpha*psi

  !    phase_new = atan2(imag(psi(0,1,llx+2)),real(psi(0,1,llx+2)))
  !     psi = psi*cexp(cmplx(0.0d0,phase_old-phase_new))

#ifdef USE_MPI
  call MPI_Bcast(psi, 4*nb*nb*nl, MPI_COMPLEX, 0, MPI_COMM_WORLD, ierr)
#endif /* USE_MPI */
      
  return
end subroutine pair_wave



