#include "../convert.F90"

subroutine sigma_calc(rank, t, sigma, chi, g, g_mtau, c_r, &
     tau, epsilon, q_tau, q_epsilon, x, y, r_tau, r_omega, &
     a_int, gamma_ph, overall_eigenvalue_max, dominant_chi_eigenvector, &
     dominant_chi_index, ft_sigma, d)

  USE CONSTANTS
  IMPLICIT NONE

#ifdef USE_MPI
  include 'mpif.h'
#endif /* USE_MPI */


  INTEGER rank
  INTEGER method
  REAL t
  COMPLEX, dimension (0:4*nb-1,0:4*nb-1,0:mp1,0:nc1) :: sigma, g, g_mtau
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
  LOGICAL ft_sigma

  COMPLEX, dimension (0:16*nb*nb-1,0:16*nb*nb-1,0:nc1) :: delta_chi, &
       delta_chi_prime

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

  ! Test of phase 
  INTEGER ib, ibp, is, isp, ind, indp
  REAL phase, phase_fix, d_phase
  COMPLEX fac

#ifdef USE_MPI
  INTEGER ierr
#endif

  pi = 4.0d0 * atan(1.0d0)
  delta_tau = (1.0d0 / t) / float(m)

 if (rank .eq. 0) then
     open(unit=9,file='my_error_file',status='old', access='append')
     write(9,*) '  start sigma_calc'
     close(unit=9)
  endif


  call chi_calc(rank, t, chi, g, g_mtau, delta_chi, delta_chi_prime, gamma_ph)

 if (rank .eq. 0) then
     open(unit=9,file='my_error_file',status='old', access='append')
     write(9,*) '  chi_calc'
     close(unit=9)
  endif


#ifdef THIRD_ORDER      
  call transform_chi(rank, t, chi, delta_chi, delta_chi_prime, &
       r_tau, r_omega, overall_eigenvalue_max, dominant_chi_eigenvector, & 
       dominant_chi_index)

 if (rank .eq. 0) then
     open(unit=9,file='my_error_file',status='old', access='append')
     write(9,*) '  transform_chi'
     close(unit=9)
  endif


#endif /* THIRD_ORDER */

  call t_generate(chi, gamma_ph, delta_chi, delta_chi_prime)

 if (rank .eq. 0) then
     open(unit=9,file='my_error_file',status='old', access='append')
     write(9,*) '  t_generate'
     close(unit=9)
  endif


#ifdef THIRD_ORDER
  call t_transform(rank, t, chi, delta_chi, delta_chi_prime, y, r_tau, r_omega)

 if (rank .eq. 0) then
     open(unit=9,file='my_error_file',status='old', access='append')
     write(9,*) '  t_transform'
     close(unit=9)
  endif

#endif /* THIRD_ORDER */

  method = 1
  call tmat_param(rank, method, chi, t, y, d, delta_chi, delta_chi_prime)

 if (rank .eq. 0) then
     open(unit=9,file='my_error_file',status='old', access='append')
     write(9,*) '  tmat_param'
     close(unit=9)
  endif


  !     Form sigma(tau, r)
  do l = 0, mp1
     do ir = 0, nc1

        do h = 0, 4*nb-1
           do i = 0, 4*nb-1

              temp_sum = cmplx( 0.0d0, 0.0d0)

              do j = 0, 4*nb-1
                 do k = 0, 4*nb-1
                    temp_sum = temp_sum - g(k,j,l,ir) * & 
                         chi(4*nb*h+k,4*nb*i+j,l,ir)
                 enddo
              enddo
              
              sigma(h,i,l,ir) = temp_sum
              
           enddo
        enddo

     enddo
  enddo

 if (rank .eq. 0) then
     open(unit=9,file='my_error_file',status='old', access='append')
     write(9,*) '  form sigma_tau'
     close(unit=9)
  endif



!!$c     if (rank. eq. 0) then
!!$
!!$c        phase = atan2(imag(sigma(8,11,0,1)), 
!!$c    $   real(sigma(8,11,0,1)))  
!!$
!!$     
!!$c      write(6,*)
!!$c      write(6,*) 'phase = ', phase
!!$c      phase_fix = 0.32d0
!!$c      d_phase = 0.1d0*(phase_fix - phase)
!!$c      write(6,*) 'd_phase = ', d_phase
!!$c      fac = cexp(cmplx(0.0d0,d_phase))
!!$c      write(6,*) 'fac = ', fac
!!$c     endif
!!$
!!$#ifdef USE_MPI
!!$c     call MPI_Bcast(fac, 1, MPI_COMPLEX, 0,
!!$c    $   MPI_COMM_WORLD, ierr)
!!$#endif /* USE_MPI */
!!$
!!$
!!$c     do ib = 0, nb-1
!!$c      do is = 0, 1
!!$c        do ibp = 0, nb-1
!!$c          do isp = 0,1
!!$
!!$c            ind = 4*ib+is
!!$c            indp = 4*ibp + isp + 2
!!$
!!$c            sigma(ind,indp,:,:) = sigma(ind,indp,:,:)*fac
!!$
!!$c            ind = 4*ib+is + 2
!!$c            indp = 4*ibp + isp
!!$
!!$c            sigma(ind,indp,:,:) = sigma(ind,indp,:,:)*conjg(fac)
!!$ 
!!$c          enddo
!!$c        enddo
!!$c     enddo
!!$c     enddo
             
  if (ft_sigma) then 
         !     Subtract analytic portion from Sigma(tau, r)
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
                                   temp_sum = temp_sum - &
                                        c_r(i0,j0,k,j,ir) * &
                                        d(i1,j1,4*nb*h+k,4*nb*i+j,ir) * &
                                        q_tau(i0,j0,l) * r_tau(i1,j1,l)
                                enddo
                             enddo
                          enddo
                       enddo
                              
                    enddo
                 enddo

                 temp_4D(l,ir) = ( sigma(h,i,l,ir) - temp_sum ) * &
                      delta_tau * cexp(cmplx(0.0d0, pi*tau(l)*t))

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

                                   temp_sum = temp_sum - &
                                        c_r(i0,j0,k,j,ir) * &
                                        d(i1,j1,4*nb*h+k,4*nb*i+j,ir) * &
                                        a_int(i0,j0,i1,j1,l)

                                enddo
                             enddo
                          enddo
                       enddo
                    enddo
                 enddo

                 temp_4D(l,ir) = temp_4D(l,ir) + temp_sum

              enddo
           enddo

           ifftv(0) = 0
           ifftv(1:3) = -1
           call fft_4D(rank, temp_4D, ifftv)

           sigma(h,i,:,:) = temp_4D

        enddo
     enddo

  endif

 if (rank .eq. 0) then
     open(unit=9,file='my_error_file',status='old', access='append')
     write(9,*) '  finish sigma_calc'
     close(unit=9)
  endif


  return
end subroutine sigma_calc
