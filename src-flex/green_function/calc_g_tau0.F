#include "../convert.F"

      subroutine calc_g_tau0( h0, t, g_tau0, g_tau0_local)

#include "../constants.F"

      COMPLEX, dimension (0:4*nb-1,0:4*nb-1,0:nl-1) :: h0,
     $   g_tau0
      COMPLEX g_tau0_local(0:4*nb-1,0:4*nb-1)
      REAL t

      COMPLEX, dimension (0:4*nb-1, 0:4*nb-1) :: 
     $   identity, id, h_temp
      REAL eigen_value(0:4*nb-1)
      COMPLEX work(20*nb*nb)
      INTEGER lwork
      REAL rwork(12*nb-2)
      INTEGER info

      INTEGER ipiv(4*nb)

      COMPLEX, dimension (0:4*nb-1,0:4*nb-1) ::
     $   a, a_inv, gk_temp1, gk_temp2
      COMPLEX gk_temp3

      INTEGER i, j, k, k1, k2, k3

      lwork = 20*nb*nb

c     Define the identity matrix.
      identity = cmplx(0.0d0, 0.0d0)
      do i = 0, 4*nb-1
        identity(i,i) = cmplx(1.0d0, 0.0d0)
      enddo

c     Loop over all lattice k-points.
      do k1 = 0, llx1
        do k2 = 0, lly1
          do k3 = 0, llz1

            k = k3 * llx * lly + k2 * llx + k1

c     Diagonalize h0 for each k-point.

            h_temp = h0(:,:,k)
            id = identity

            call chegv(1,'V','U', 4*nb, h_temp, 4*nb, id, 4*nb, 
     $         eigen_value, work, lwork, rwork, info) 

            if (info .ne. 0) then
              write(6,*) 'zhegv failed in calc_g_tau0'
            endif

            a = h_temp
            id = identity
           
            call cgesv(4*nb, 4*nb, h_temp, 4*nb, 
     $         ipiv, id, 4*nb, info)

            if (info .ne. 0) then
              write(6,*) 'zgesv failed in calc_g_tau0'
            endif

            a_inv = id
            
c     Create the diagonal form of gk_tau0
            gk_temp1 = cmplx(0.0d0, 0.0d0)

            do i = 0, 4*nb-1

              if (eigen_value(i) .gt. 0.0d0) then

                gk_temp1(i,i) = -1.0d0 / 
     $             (exp(-eigen_value(i)/t) + 1.0d0)

              else

                gk_temp1(i,i) = -exp(eigen_value(i)/t) / 
     $             (exp(eigen_value(i)/t) + 1.0d0)

              endif

            enddo

c     Multiple the diagonal force by diagonalization matrix, A_inv
            call cgemm('N','N',4*nb,4*nb,4*nb,cmplx(1.0d0,0.0d0),
     $         gk_temp1,4*nb,a_inv,4*nb,cmplx(0.0d0,0.0d0),
     $         gk_temp2,4*nb)

c     Multiply the previous produce by A.
            call cgemm('N','N',4*nb,4*nb,4*nb,cmplx(1.0d0,0.0d0),
     $         a,4*nb,gk_temp2,4*nb,cmplx(0.0d0,0.0d0),
     $         g_tau0(:,:,k),4*nb)

          enddo
        enddo
      enddo
      
      g_tau0_local = cmplx(0.0d0, 0.0d0)

      do k = 0, nl-1
        g_tau0_local =  g_tau0_local + g_tau0(:,:,k)
      enddo
          
      g_tau0_local = g_tau0_local / float(nl)

      return
      end
