#include "../convert.F90"

subroutine pair_field_energy(g_tau0, prfld, psi, pair_energy)

  USE CONSTANTS

  COMPLEX g_tau0(0:4*nb-1,0:4*nb-1,0:nl-1)
  REAL prfld
  COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)
  REAL pair_energy

  COMPLEX tr_phi_psi
  INTEGER nu, is, nup, isp
  INTEGER k1, k2, k3, k, km
  INTEGER ind, indp, nua, nuap

  tr_phi_psi = cmplx(0.0d0, 0.0d0)

  do nu = 0, nb-1
     do is = 0, 1

        ind = 2*nu + is
        nua = 4*nu+is

        do nup = 0, nb-1
           do isp = 0, 1

              indp = 2*nup + isp
              nuap = 4*nup + isp + 2 

              do k1 = 0, llx1
                 do k2 = 0, lly1
                    do k3 = 0, llz1
                       
                       k = k1 + k2*llx + k3*llx*lly
                       km = mod(-k1+llx,llx) + mod(-k2+lly,lly)*llx + &
                            mod(-k3+llz,llz)*llx*lly

                       tr_phi_psi = tr_phi_psi - & 
                            psi(ind,indp,km) * g_tau0(nua,nuap,k) &
                            + conjg(psi(indp,ind,k)) * g_tau0(nuap,nua,k)

                    enddo
                 enddo
              enddo
              
           enddo
        enddo
     enddo
  enddo

  write(6,*) 'Debug: tr_phi_psi ', tr_phi_psi
  pair_energy = -0.5d0*prfld*real(tr_phi_psi) / float(nl)

  return
end subroutine pair_field_energy
