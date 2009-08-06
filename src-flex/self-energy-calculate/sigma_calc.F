#include "../convert.F"

      subroutine sigma_calc(rank, t, sigma, chi, g, g_mtau,
     $   c_r, tau, epsilon, q_tau, q_epsilon, x, y, 
     $   r_tau, r_omega, a_int, gamma_ph,
     $   overall_eigenvalue_max,  dominant_chi_eigenvector,
     $   dominant_chi_index)

#include "../constants.F"

      INTEGER rank
      REAL t
      COMPLEX, dimension (0:4*nb-1,0:4*nb-1,0:mp1,0:nc1) ::
     $   sigma, g, g_mtau
      COMPLEX chi(0:16*nb*nb-1,0:16*nb*nb-1,0:mp1,0:nc1)
      COMPLEX c_r(0:1,0:1,0:4*nb-1,0:4*nb-1,0:nc1)
      REAL tau(0:mp1), epsilon(0:mp1)
      REAL q_tau(0:1,0:1,0:mp1)
      COMPLEX q_epsilon(0:1,0:1,0:mp1)
      REAL x(0:1,0:1), y(0:1,0:1)
      REAL r_tau(0:1,0:1,0:mp1)
      COMPLEX r_omega(0:1,0:1,0:mp1)
      COMPLEX a_int(0:1,0:1,0:1,0:1,0:mp1)
      COMPLEX gamma_ph(0:16*nb*nb-1, 0:16*nb*nb-1)

      COMPLEX, dimension (0:16*nb*nb-1,0:16*nb*nb-1,0:nc1) ::
     $   delta_chi, delta_chi_prime

      COMPLEX d(0:1,0:1,0:16*nb*nb-1,0:16*nb*nb-1,0:nc1)
      COMPLEX sum
      COMPLEX temp_sum

      INTEGER h, i, j, k
      INTEGER l, ix, iy, iz, ir
      INTEGER kx, ky, kz
      INTEGER ifftv(0:3)

      INTEGER i0, i1, j0, j1

      COMPLEX temp_4D(0:mp1, 0:nc1)

      REAL delta_tau 
      REAL pi

      COMPLEX dominant_chi_eigenvector(0:16*nb*nb-1)
      REAL overall_eigenvalue_max
      INTEGER dominant_chi_index(0:1)

      pi = 4.0d0 * atan(1.0d0)
      delta_tau = (1.0d0 / t) / float(m)

      call chi_calc(rank, t, chi, g, g_mtau, delta_chi,
     $   delta_chi_prime, gamma_ph)

#ifdef THIRD_ORDER      
      call transform_chi(rank, t, chi, delta_chi, 
     $   delta_chi_prime, r_tau, r_omega,
     $   overall_eigenvalue_max, dominant_chi_eigenvector, 
     $   dominant_chi_index)
#endif /* THIRD_ORDER */

      call t_generate(chi, gamma_ph, delta_chi,
     $   delta_chi_prime)

#ifdef THIRD_ORDER
      call t_transform(rank, t, chi, delta_chi, 
     $   delta_chi_prime, y, r_tau, r_omega)
#endif /* THIRD_ORDER */

      call tmat_param(rank, 1, chi, t, y, 
     $  d, delta_chi, delta_chi_prime)

c     Form sigma(tau, r)
      do l = 0, mp1
        do ir = 0, nc1

          do h = 0, 4*nb-1
            do i = 0, 4*nb-1

              temp_sum = cmplx( 0.0d0, 0.0d0)

              do j = 0, 4*nb-1
                do k = 0, 4*nb-1
                  temp_sum = temp_sum -
     $               g(k,j,l,ir) * 
     $               chi(4*nb*h+k,4*nb*i+j,l,ir)
                enddo
              enddo

              sigma(h,i,l,ir) = temp_sum
              
            enddo
          enddo

        enddo
      enddo

c     Subtract analytic portion from Sigma(tau, r)
      do h = 0, 4*nb-1
        do i = 0, 4*nb-1

          do l = 0, mp1
            do ir = 0, nc1

              temp_sum = cmplx(0.0d0, 0.0d0)
                  
              do j = 0, 4*nb-1
                do k = 0, 4*nb-1

                  do i0 = 0, 1
                    do j0 = 0, 1
                      do i1 = 0, 1
                        do j1 = 0, 1
                          temp_sum = temp_sum -
     $                       c_r(i0,j0,k,j,ir) *
     $                       d(i1,j1,4*nb*h+k,4*nb*i+j,ir) *
     $                       q_tau(i0,j0,l) * r_tau(i1,j1,l)
                        enddo
                      enddo
                    enddo
                  enddo
                              
                enddo
              enddo

              temp_4D(l,ir) = 
     $           ( sigma(h,i,l,ir) - temp_sum ) *
     $           delta_tau * 
     $           cexp(cmplx(0.0d0, pi*tau(l)*t))

            enddo
          enddo

          ifftv(0) =  1
          ifftv(1:3) =  0

          call fft_4D(rank, temp_4D, ifftv)

          do l = 0, mp1
            do ir = 0, nc1

              temp_sum = cmplx(0.0d0, 0.0d0)

              do j = 0, 4*nb-1
                do k = 0, 4*nb-1
                  do i0 = 0, 1
                    do j0 = 0, 1
                      do i1 = 0, 1
                        do j1 = 0, 1

                          temp_sum = temp_sum -
     $                       c_r(i0,j0,k,j,ir) *
     $                       d(i1,j1,4*nb*h+k,4*nb*i+j,ir) *
     $                       a_int(i0,j0,i1,j1,l)

                        enddo
                      enddo
                    enddo
                  enddo
                enddo
              enddo

              temp_4D(l,ir) = 
     $           temp_4D(l,ir) + temp_sum

            enddo
          enddo

          ifftv(0) = 0
          ifftv(1:3) = -1
          call fft_4D(rank, temp_4D, ifftv)

          sigma(h,i,:,:) = temp_4D

        enddo
      enddo

      return
      end
