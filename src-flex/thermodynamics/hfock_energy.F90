#include "../convert.F90"

subroutine hfock_energy(sigma1, g_tau0_local, e_hf)

  USE CONSTANTS

  COMPLEX, dimension (0:4*nb-1, 0:4*nb-1) :: sigma1
  COMPLEX g_tau0_local(0:4*nb-1,0:4*nb-1)
  REAL e_hf

  INTEGER nua, nuap
  COMPLEX g0, e_hf_temp, e_hf_temp2

  e_hf_temp = cmplx(0.0d0, 0.0d0)
  e_hf_temp2 = cmplx(0.0d0, 0.0d0)

  do nua = 0, 4*nb-1

     do nuap = 0, 4*nb-1

        g0 = g_tau0_local(nuap, nua)

        if ( (nua .eq. nuap) .and. (mod(nua,4) .le. 1) ) then
           g0 = g0 + 1.0d0
        endif

        if ( mod(nua,4) .le. 1) then
           e_hf_temp = e_hf_temp + sigma1(nua,nuap)*g0
        else
           e_hf_temp2 = e_hf_temp2 + sigma1(nua,nuap)*g0 
        endif

     enddo

  enddo

  e_hf = real(e_hf_temp)/2.0d0
  write(6,*) 'Debug1: trace sig1*g = ', e_hf_temp
  write(6,*) 'Debug1: trace sig1*g = ', e_hf_temp2

  return
end subroutine hfock_energy
