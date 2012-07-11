#include "../convert.F90"

subroutine t_generate(t_mat, gamma_ph, delta_t, delta_t_prime)

  USE CONSTANTS

  COMPLEX t_mat(0:16*nb*nb-1, 0:16*nb*nb-1,0:mp1,0:nc1)
  COMPLEX gamma_ph(0:16*nb*nb-1, 0:16*nb*nb-1)
  COMPLEX, dimension (0:16*nb*nb-1,0:16*nb*nb-1,0:nc1) :: delta_t, delta_t_prime

  INTEGER ipiv(16*nb*nb), info

  COMPLEX, dimension (0:16*nb*nb-1, 0:16*nb*nb-1) :: delta, chi_squared, &
       temp_product, temp_sum, temp_chi, temp_array, temp_inverse

  INTEGER i, j, k, l

  delta = cmplx(0.0d0, 0.0d0)
  do i = 0, 16*nb*nb-1
     delta(i,i) = cmplx(1.0d0, 0.0d0)
  enddo

  do k = 0, nc1
     do l = 0, mp1
          
        temp_chi = t_mat(:,:,l,k)

#ifdef FLEX
        temp_array = delta - temp_chi
        temp_inverse = delta

        call cgesv(16*nb*nb, 16*nb*nb, temp_array, & 
             16*nb*nb, ipiv, temp_inverse, 16*nb*nb, info)
#endif /* FLEX */

#ifdef THIRD_ORDER
        !     Calculate chi_squared
        call cgemm('N','N', 16*nb*nb, 16*nb*nb, 16*nb*nb, &
             cmplx(1.0d0,0.0d0), temp_chi, 16*nb*nb, &
             temp_chi, 16*nb*nb, cmplx(0.0d0,0.0d0), &
             chi_squared, 16*nb*nb)

        !     Calculate the product of chi_squared with minus chi
#ifdef FLEX
        call cgemm('N','N', 16*nb*nb, 16*nb*nb, 16*nb*nb, &
             cmplx(1.0d0,0.0d0), chi_squared, 16*nb*nb, &
             temp_inverse, 16*nb*nb, cmplx(0.0d0,0.0d0), &
             temp_product, 16*nb*nb)
#else
        temp_product = chi_squared
#endif /* FLEX */
        !     Calculate the temporary sum of 1/2*chi with temp_product
        temp_sum = temp_chi / 3.0d0 + temp_product
#else
        temp_sum = temp_chi / 3.0d0
#endif /* THIRD_ORDER */              

        call cgemm('N','N', 16*nb*nb, 16*nb*nb, 16*nb*nb, &
             cmplx(1.0d0,0.0d0), temp_sum, 16*nb*nb, &
             gamma_ph, 16*nb*nb, cmplx(0.0d0,0.0d0), temp_product, 16*nb*nb)

        t_mat(:,:,l,k) = temp_product

     enddo
  enddo

  !     Calculate the discontinuities in T and T_prime
  do k = 0, nc1

     call cgemm('N','N', 16*nb*nb, 16*nb*nb, 16*nb*nb, &
          cmplx(1.0d0,0.0d0), delta_t(:,:,k), 16*nb*nb, &
          gamma_ph, 16*nb*nb, cmplx(0.0d0,0.0d0), &
          temp_product, 16*nb*nb)

     temp_product = temp_product / 3.0d0
     delta_t(:,:,k) = temp_product

  enddo

  do k = 0, nc1

#ifdef THIRD_ORDER
     call cgemm('N','N', 16*nb*nb, 16*nb*nb, 16*nb*nb, &
          cmplx(1.0d0,0.0d0), delta_t(:,:,k), 16*nb*nb, &
          delta_t(:,:,k), 16*nb*nb, cmplx(0.0d0,0.0d0), chi_squared, 16*nb*nb)
    
     temp_sum = delta_t_prime(:,:,k) / 3.0d0 +  chi_squared
#else
     temp_sum = delta_t_prime(:,:,k) / 3.0d0   
#endif /* THIRD_ORDER */

     call cgemm('N','N', 16*nb*nb, 16*nb*nb, 16*nb*nb, &
          cmplx(1.0d0,0.0d0), temp_sum, 16*nb*nb, &
          gamma_ph, 16*nb*nb, cmplx(0.0d0,0.0d0), &
          delta_t_prime(:,:,k), 16*nb*nb)
     
  enddo

  return
end subroutine t_generate

  









