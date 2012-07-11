#include "../convert.F90"

subroutine t_transform(rank, t, t_mat, delta_t, delta_t_prime, y, &
     r_tau, r_omega)

  USE CONSTANTS

  INTEGER rank
  REAL t
  COMPLEX t_mat(0:16*nb*nb-1,0:16*nb*nb-1,0:mp1,0:nc1)
  COMPLEX, dimension (0:16*nb*nb-1,0:16*nb*nb-1,0:nc1) :: &
       delta_t, delta_t_prime
  REAL y(0:1,0:1)
  REAL r_tau(0:1,0:1,0:mp1)
  COMPLEX r_omega(0:1,0:1,0:mp1)

  COMPLEX d(0:1,0:1,0:16*nb*nb-1,0:16*nb*nb-1,0:nc1)
  COMPLEX sum_dr
  COMPLEX temp_3D(0:nc1)
  COMPLEX temp_4D(0:mp1, 0:nc1)      

  INTEGER i, j, l, k
  INTEGER ia, ib

  INTEGER ifftv(0:3)
  REAL delta_tau
  COMPLEX t_prime

  delta_tau = (1.0d0 / t) / float(m)

  !     Generate the new D arrays
  d(0,1,:,:,:) = cmplx(0.0d0, 0.0d0)
  d(0,0,:,:,:) = -delta_t 
  d(1,1,:,:,:) = cmplx(0.0d0, 0.0d0)
  d(1,0,:,:,:) = delta_t_prime 

  do i = 0, 16*nb*nb-1
     do j = 0, 16*nb*nb-1

        do l = 0, mp1
           do k = 0, nc1

              sum_dr = cmplx(0.0d0, 0.0d0)
              do ia = 0, 1
                 do ib = 0, 1
                    sum_dr = sum_dr + R_omega(ia,ib,l) * d(ia,ib,i,j,k)
                 enddo
              enddo

              temp_4D(l,k) = t_mat(i,j,l,k) - sum_dr

           enddo
        enddo

        !     Transform on the time-index
        ifftv(0) = -1
        ifftv(1) = 0
        ifftv(2) = 0
        ifftv(3) = 0

        call fft_4D(rank, temp_4D, ifftv)

        !     Add the sum back to the transformed T arrays
        do l = 0, mp1
           do k = 0, nc1
 
              sum_dr = cmplx(0.0d0, 0.0d0)
              do ia = 0, 1
                 do ib = 0, 1
                    sum_dr = sum_dr + r_tau(ia,ib,l) * d(ia,ib,i,j,k)
                 enddo
              enddo

              temp_4D(l,k) = t * temp_4D(l,k) + sum_dr

           enddo
        enddo
  
        !     Transform on the space-indices
        ifftv(0) = 0
        ifftv(1) = 1
        ifftv(2) = 1
        ifftv(3) = 1

        call fft_4D(rank, temp_4D, ifftv)

        t_mat(i,j,:,:) = temp_4D / float(nc)

     enddo
  enddo

  !     Transform delta_t and delta_t_prime
  ifftv(0) = 0
  ifftv(1) = 1
  ifftv(2) = 1
  ifftv(3) = 1

  do i = 0, 16*nb*nb-1
     do j = 0, 16*nb*nb-1
          
        temp_3D = delta_t(i,j,:)
        call fft_3D( temp_3D, ifftv)
        delta_t(i,j,:) = temp_3D / float(nc)

        temp_3D = delta_t_prime(i,j,:)
        call fft_3D( temp_3D, ifftv)
        delta_t_prime(i,j,:) = temp_3D / float(nc)

     enddo
  enddo

  return
end subroutine t_transform
