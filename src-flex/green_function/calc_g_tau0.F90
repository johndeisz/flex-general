#include "../convert.F90"

subroutine calc_g_tau0( tij, ed, v_pert_eff, psi, h_eff, prfld_eff, mu, &
     sigma1, h_so, t, g_tau0, g_tau0_local)

  USE CONSTANTS
  USE h_zero
  IMPLICIT NONE

  REAL ed(0:nb-1)
  COMPLEX tij(0:nb-1,0:nb-1,-2:2,-2:2,-2:2)
  COMPLEX ek(0:nb-1,0:nb-1,0:nl-1)
  REAL v_pert_eff(0:nb-1)
  COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)
  REAL h_eff(0:nb-1,1:3)
  REAL prfld_eff
  REAL mu
  COMPLEX, dimension (0:4*nb-1,0:4*nb-1) :: sigma1
  COMPLEX h_so(0:2*nb-1, 0:2*nb-1)
  
  COMPLEX, dimension (0:4*nb-1,0:4*nb-1,0:nl-1) :: g_tau0
  COMPLEX g_tau0_local(0:4*nb-1,0:4*nb-1)
  REAL t

  COMPLEX, dimension (0:4*nb-1, 0:4*nb-1) :: identity, id, h_temp
  REAL eigen_value(0:4*nb-1)
  COMPLEX work(20*nb*nb)
  INTEGER lwork
  REAL rwork(12*nb-2)
  INTEGER info

  INTEGER ipiv(4*nb)

  COMPLEX, dimension (0:4*nb-1,0:4*nb-1) :: a, a_inv, gk_temp1, gk_temp2
  COMPLEX gk_temp3

  INTEGER i, j, k, k1, k2, k3

  lwork = 20*nb*nb

  !     Define the identity matrix.
  identity = cmplx(0.0d0, 0.0d0)
  do i = 0, 4*nb-1
     identity(i,i) = cmplx(1.0d0, 0.0d0)
  enddo

  !     Loop over all lattice k-points.
  do k1 = 0, llx1
     do k2 = 0, lly1
        do k3 = 0, llz1

           k = k3 * llx * lly + k2 * llx + k1

           !     Diagonalize h0 for each k-point.

           h_temp = h0(k,tij, ed, v_pert_eff, psi, h_eff, prfld_eff, mu, &
                sigma1, h_so)
           id = identity

           call chegv(1,'V','U', 4*nb, h_temp, 4*nb, id, 4*nb, & 
                eigen_value, work, lwork, rwork, info) 

           if (info .ne. 0) then
              write(6,*) 'zhegv failed in calc_g_tau0'
           endif

           a = h_temp
           id = identity
           
           call cgesv(4*nb, 4*nb, h_temp, 4*nb, ipiv, id, 4*nb, info)

           if (info .ne. 0) then
              write(6,*) 'zgesv failed in calc_g_tau0'
           endif

           a_inv = id
            
           !     Create the diagonal form of gk_tau0
           gk_temp1 = cmplx(0.0d0, 0.0d0)

           do i = 0, 4*nb-1

              if (eigen_value(i) .gt. 0.0d0) then

                 gk_temp1(i,i) = -1.0d0 / (exp(-eigen_value(i)/t) + 1.0d0)

              else

                 gk_temp1(i,i) = -exp(eigen_value(i)/t) / &
                      (exp(eigen_value(i)/t) + 1.0d0)

              endif

           enddo

           ! Multiple the diagonal force by diagonalization matrix, A_inv
           call cgemm('N','N',4*nb,4*nb,4*nb,cmplx(1.0d0,0.0d0), &
                gk_temp1,4*nb,a_inv,4*nb,cmplx(0.0d0,0.0d0), gk_temp2,4*nb)

           !     Multiply the previous produce by A.
           call cgemm('N','N',4*nb,4*nb,4*nb,cmplx(1.0d0,0.0d0), &
                a,4*nb,gk_temp2,4*nb,cmplx(0.0d0,0.0d0), g_tau0(:,:,k),4*nb)

        enddo
     enddo
  enddo
      
  g_tau0_local = cmplx(0.0d0, 0.0d0)

  do k = 0, nl-1
     g_tau0_local =  g_tau0_local + g_tau0(:,:,k)
  enddo
          
  g_tau0_local = g_tau0_local / float(nl)

  return
end subroutine calc_g_tau0
