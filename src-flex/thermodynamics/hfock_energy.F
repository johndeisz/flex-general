#include "../convert.F"

      subroutine hfock_energy(sigma1, g_tau0_local, e_hf)

#include "../constants.F"

      COMPLEX, dimension (0:4*nb-1, 0:4*nb-1) :: sigma1
      COMPLEX g_tau0_local(0:4*nb-1,0:4*nb-1)
      REAL e_hf

      INTEGER nua, nuap
      COMPLEX g0, e_hf_temp

      e_hf_temp = cmplx(0.0d0, 0.0d0)

      do nua = 0, 4*nb-1
        do nuap = 0, 4*nb-1

          g0 = g_tau0_local(nuap, nua)

          if ( (nua .eq. nuap) .and. (mod(nua,4) .le. 1) ) then
            g0 = g0 + 1.0d0
          endif

          e_hf_temp = e_hf_temp + sigma1(nua,nuap)*g0

        enddo
      enddo

      e_hf = real(e_hf_temp)/4.0d0
      write(6,*) 'Debug: trace sig1*g = ', e_hf_temp

      return
      end
