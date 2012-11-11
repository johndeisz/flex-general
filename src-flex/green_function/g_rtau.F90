#include "../convert.F90"

subroutine g_rtau(rank, g, t, c_r, q_epsilon, q_tau, tau)

  USE CONSTANTS
  IMPLICIT NONE

  INTEGER rank
  COMPLEX g(0:4*nb-1,0:4*nb-1,0:mp1,0:nc1) 
  REAL t
  COMPLEX c_r(0:1,0:1,0:4*nb-1,0:4*nb-1,0:nc1) 
  COMPLEX q_epsilon(0:1,0:1,0:mp1)
  REAL q_tau(0:1,0:1,0:mp1), tau(0:mp1)

  COMPLEX c_k(0:1,0:1,0:4*nb-1,0:4*nb-1,0:nc1) 

  INTEGER ifftv(0:3), jfftv(0:3) 

  COMPLEX temp_ck(0:nc1)
  
  INTEGER ia, ib, l
  REAL pi
  INTEGER na1, na2

  COMPLEX temp_g(0:mp1, 0:nc1)

  pi =  acos(-1.0d0)

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
  
  !     The order here needs to be rearranged.  See the disorder code.

  do ia = 0, 1
     do ib = 0, 1

        do l = 0, mp1

           g(:,:,l,:) = g(:,:,l,:) - c_k(ia,ib,:,:,:) * q_epsilon(ia,ib,l) 

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

           g(na1,na2,l,:) =  t * cexp(-cmplx( 0.0d0, pi * t * tau(l) )) * &  
                temp_g(l,:) / float(nc)

        enddo

     enddo
  enddo

  !     Restore analytic terms.

  do ia = 0, 1
     do ib = 0, 1

        do l = 0, mp1

           g(:,:,l,:) = g(:,:,l,:) + c_r(ia,ib,:,:,:) * q_tau(ia,ib,l)

        enddo
        
     enddo
  enddo
  
  return
end subroutine g_rtau






























