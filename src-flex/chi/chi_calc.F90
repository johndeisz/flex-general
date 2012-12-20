#include "../convert.F90"

subroutine chi_calc(rank, t, chi, g, g_mtau, delta_chi, &
     delta_chi_prime, gamma_ph)

  USE CONSTANTS
  IMPLICIT NONE

#ifdef USE_MPI
  include 'mpif.h'
#endif

  INTEGER rank
  REAL t
  COMPLEX chi(0:16*nb*nb-1,0:16*nb*nb-1,0:mp1,0:nc1)
  COMPLEX, dimension (0:4*nb-1,0:4*nb-1,0:mp1,0:nc1) :: g, g_mtau
  COMPLEX, dimension (0:16*nb*nb-1,0:16*nb*nb-1,0:nc1) :: delta_chi, &
       delta_chi_prime
  COMPLEX gamma_ph(0:16*nb*nb-1, 0:16*nb*nb-1)

  REAL delta_tau
  COMPLEX product
  COMPLEX matrix_product(0:16*nb*nb-1, 0:16*nb*nb-1)
  COMPLEX chi_mtau(0:2)
  COMPLEX chi_prime_plus, chi_prime_minus

  INTEGER h, i, j, k
  INTEGER l, ix, iy, iz, ir, irm

  INTEGER ierr

  delta_tau = (1.0d0 / t) / float(m)

  !     Calculate chi as a function of r, tau

  do h = 0, 4*nb-1
     do i = 0, 4*nb-1
        do j = 0, 4*nb-1
           do k = 0, 4*nb-1
              chi(4*nb*h+j, 4*nb*i+k, :, :) = -g(h,i,:,:)*g_mtau(k,j,:,:)
           enddo
        enddo
     enddo
  enddo
  

  !     Have the processor with the smallest tau values find 
  !     the change close to 0
  
  if (rank .eq. 0) then

     do ix = 0, lcx1
        do iy = 0, lcy1
           do iz = 0, lcz1

              ir = ix + iy*lcx + iz*lcx*lcy
              irm = mod(-ix+lcx,lcx) + mod(-iy+lcy,lcy)*lcx + &
                   mod(-iz+lcz,lcz)*lcx*lcy

              do h = 0, 4*nb-1
                 do i = 0, 4*nb-1
                    do j = 0, 4*nb-1
                       do k = 0, 4*nb-1
                          
                          do l = 0, 2
                             chi_mtau(l) = -g_mtau(h,i,l,irm)*g(k,j,l,irm)
                          enddo
		  
                          delta_chi(4*nb*h+j,4*nb*i+k,ir) = &
                               chi(4*nb*h+j,4*nb*i+k,0,ir) - chi_mtau(0)
		  
                          chi_prime_plus = &
                               2.0d0 * chi(4*nb*h+j,4*nb*i+k,1,ir) - &
                               0.5d0 * chi(4*nb*h+j,4*nb*i+k,2,ir) - &
                               1.5d0 * chi(4*nb*h+j,4*nb*i+k,0,ir)
                  
                          chi_prime_plus = chi_prime_plus / delta_tau
		  
                          chi_prime_minus =  -2.0d0 * chi_mtau(1) + & 
                               0.5d0 * chi_mtau(2) + 1.5d0 * chi_mtau(0)
                  
                          chi_prime_minus = chi_prime_minus / delta_tau

                          delta_chi_prime(4*nb*h+j,4*nb*i+k,ir) = &
                               chi_prime_plus - chi_prime_minus
                    
                       enddo
                    enddo
                 enddo
              enddo

              ! Multiply the results by 
              !gamma to obtain deltachi and deltachiprime

              call cgemm('N','N', 16*nb*nb, 16*nb*nb, 16*nb*nb, &
                   cmplx(1.0d0,0.0d0), gamma_ph, 16*nb*nb,  &
                   delta_chi(:,:,ir), 16*nb*nb, cmplx(0.0d0,0.0d0), &
                   matrix_product, 16*nb*nb)

              delta_chi(:,:,ir) = 0.5d0 * matrix_product

              call cgemm('N','N', 16*nb*nb, 16*nb*nb, 16*nb*nb, &
                   cmplx(1.0d0,0.0d0), gamma_ph, 16*nb*nb, & 
                   delta_chi_prime(:,:,ir), 16*nb*nb, cmplx(0.0d0,0.0d0), &
                   matrix_product, 16*nb*nb)
              
              delta_chi_prime(:,:,ir) = 0.5d0 * matrix_product

           enddo
        enddo
     enddo

  endif
      
#ifdef USE_MPI
  do i=0, 16*nb*nb-1 
    do j=0, 16*nb*nb-1
!      call MPI_Bcast(delta_chi, 16*16*nb*nb*nb*nb*nc, &
!         MPI_COMPLEX, 0, MPI_COMM_WORLD, ierr)
!      call MPI_Bcast(delta_chi_prime, 16*16*nb*nb*nb*nb*nc, &
!         MPI_COMPLEX, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(delta_chi(i,j,:), nc, &
         MPI_COMPLEX, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(delta_chi_prime(i,j,:), nc, &
         MPI_COMPLEX, 0, MPI_COMM_WORLD, ierr)
    enddo
   enddo
#endif

  !     Multiply chi by gamma

  do l = 0, mp1
     do ir = 0, nc1

        call cgemm('N','N', 16*nb*nb, 16*nb*nb, 16*nb*nb, &
             cmplx(1.0d0,0.0d0), gamma_ph, 16*nb*nb, & 
             chi(:,:,l,ir), 16*nb*nb, cmplx(0.0d0,0.0d0), &
             matrix_product, 16*nb*nb)         

        chi(:,:,l,ir) = 0.5d0 * matrix_product

     enddo
  enddo

  !     Multiplye results by (-1) to account for the number of
  !     appearances of the vertex function

  chi = -chi
  delta_chi = -delta_chi
  delta_chi_prime = -delta_chi_prime

  return
end subroutine chi_calc
