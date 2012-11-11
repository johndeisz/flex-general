#include "../convert.F90"

subroutine g_minus_tau(rank, g, t, c_r, q_epsilon, q_tau, tau)

  USE CONSTANTS
  IMPLICIT NONE

  INTEGER rank
  COMPLEX g(0:4*nb-1, 0:4*nb-1, 0:mp1,0:nc1)
  REAL t 
  COMPLEX c_r(0:1,0:1,0:4*nb-1,0:4*nb-1,0:nc1)
  COMPLEX q_epsilon(0:1,0:1,0:mp1) 
  REAL q_tau(0:1,0:1,0:mp1) 
  REAL tau(0:mp1)

  COMPLEX c_k(0:1,0:1,0:4*nb-1,0:4*nb-1,0:nc1)
  COMPLEX temp_g(0:mp1, 0:nc1)

  INTEGER ifftv(0:3), jfftv(0:3)

  COMPLEX c_r_temp(0:1,0:1,0:4*nb-1,0:4*nb-1,0:nc1)

  COMPLEX temp_ck(0:nc1)
  REAL conjg_fac

  INTEGER ia, ib, i, j, l, k1, k2, k3, kx, ky, kz
  INTEGER na1, na2
  REAL pi

  COMPLEX g11

  !     constants
  pi =  acos(-1.0d0)

  !     First conjugate g
  g = conjg(g)

  !     Form c_k

  jfftv(1) =  -1
  jfftv(2) =  -1
  jfftv(3) =  -1

  do ia = 0, 1
     do ib = 0, 1
        do na1 = 0, 4*nb-1
           do na2 = 0, 4*nb-1

              temp_ck = c_r(ia,ib,na1,na2,:)
              call fft_3D ( temp_ck, jfftv )
              c_k(ia,ib,na1,na2,:) = temp_ck

           enddo
        enddo
     enddo
  enddo

  !     Subtractions.
  !     The order here needs to be rearranged.  See the disorder code.

  do ia = 0, 1
     do ib = 0, 1
        do l = 0, mp1
           g(:,:,l,:) = g(:,:,l,:) - &
                conjg(c_k(ia,ib,:,:,:))*conjg(q_epsilon(ia,ib,l))
        enddo
     enddo
  enddo
  
  !     setup for fft

  ifftv(0) = -1 
  ifftv(1) =  1
  ifftv(2) =  1
  ifftv(3) =  1

  do na1 = 0, 4*nb-1
     do na2 = 0, 4*nb-1

        temp_g = g(na1,na2,:,:)
        call fft_4D (rank, temp_g, ifftv) ! Performs the fft.
        do l = 0, mp1
           g(na1,na2,l,:) = t * &  
                cexp(-cmplx( 0.0d0, pi * t * tau(l) )) * &  
                temp_g(l,:) / float(nc)
        enddo

     enddo
  enddo
      
  !     Form the FFT of the conjugate of c_k.

  jfftv(1) =  1
  jfftv(2) =  1
  jfftv(3) =  1

  do ia = 0, 1
     do ib = 0, 1
        do na1 = 0, 4*nb-1
           do na2 = 0, 4*nb-1
              temp_ck = conjg(c_k(ia,ib,na1,na2,:))
              call fft_3D ( temp_ck, jfftv )
              c_r_temp(ia,ib,na1,na2,:) = temp_ck / float(nc)
           enddo
        enddo
     enddo
  enddo
  
  !     Add back analytic terms. Note that the FFT of Q0* is the negative
  !     of the FFT of Q0, but the FFT of Q1* is the same as the FFT of Q1

  do ia = 0, 1

     if (ia .eq. 0) then 
        conjg_fac = -1.0d0
     else
        conjg_fac = 1.0d0
     endif

     do ib = 0, 1

        do l = 0, mp1

           g(:,:,l,:) = g(:,:,l,:) + c_r_temp(ia,ib,:,:,:) * &
                conjg_fac * q_tau(ia,ib,l)

        enddo

     enddo
  enddo

  !     Reconjugate to obtain the final result.
  g = conjg(g)

  return
end subroutine g_minus_tau
