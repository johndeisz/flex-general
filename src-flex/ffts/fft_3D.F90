#include "../convert.F90"

subroutine fft_3D (mat, isignv)

  USE CONSTANTS
  IMPLICIT NONE

  COMPLEX mat(0:nc-1)
  INTEGER isignv(0:3)

#include "../fftw3.f"

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !     isignv controls whether or not an fft is done on a 
  !     given axis and the direction of the fft.   
  !     nu-k -> tau-r  is (1, -1, -1, -1)
  !     no scaling is done -- this is done by hand in the flux code

  
  COMPLEX matx_in(0:lcx-1), matx_out(0:lcx-1)
  COMPLEX maty_in(0:lcy-1), maty_out(0:lcy-1)
  COMPLEX matz_in(0:lcz-1), matz_out(0:lcz-1)

  INTEGER j1, j2, j3, ind
  INTEGER fftw_plan
  INTEGER*8 plan
  
  !     Do x-axis Fourier transform.
        
  if ( (lcx .gt. 1) .and. (isignv(1) .ne. 0) ) then

     call sfftw_plan_dft_1d(plan,lcx,matx_in,matx_out, isignv(1),FFTW_ESTIMATE)

     do j2 = 0, lcy-1
        do j3 = 0, lcz-1

           do j1 = 0, lcx-1
              ind = j1 + j2*lcx + j3*lcx*lcy
              matx_in(j1) = mat(ind)
           enddo
            
           call sfftw_execute(plan,matx_in, matx_out)

           do j1 = 0, lcx-1
              ind = j1 + j2*lcx + j3*lcx*lcy            
              mat(ind) = matx_out(j1)
           enddo

        enddo
     enddo

     call sfftw_destroy_plan(plan)

  endif

  if ( (lcy .gt. 1) .and. (isignv(2) .ne. 0) ) then

     call sfftw_plan_dft_1d(plan,lcy,maty_in,maty_out, isignv(2),FFTW_ESTIMATE)

     do j1 = 0, lcx-1
        do j3 = 0, lcz-1

           do j2 = 0, lcy-1
              ind = j1 + j2*lcx + j3*lcx*lcy
              maty_in(j2) = mat(ind)
           enddo

           call sfftw_execute(plan,maty_in, maty_out)
           
           do j2 = 0, lcy-1
              ind = j1 + j2*lcx + j3*lcx*lcy           
              mat(ind) = maty_out(j2)
           enddo
            
        enddo
     enddo

     call sfftw_destroy_plan(plan)
        
  endif ! (ny > 1 .and. isignv(2) != 0)

  if ( (lcz .gt. 1) .and. (isignv(3) .ne. 0) ) then

     call sfftw_plan_dft_1d(plan,lcz,matz_in,matz_out, isignv(3),FFTW_ESTIMATE)

     do j1 = 0, lcx-1
        do j2 = 0, lcy-1

           do j3 = 0, lcz-1
              ind = j1 + j2*lcx + j3*lcx*lcy              
              matz_in(j3) = mat(ind)
           enddo

           call sfftw_execute(plan,matz_in, matz_out)

           do j3 = 0, lcz-1
              ind = j1 + j2*lcx + j3*lcx*lcy
              mat(ind) = matz_out(j3)
           enddo

        enddo
     enddo

     call sfftw_destroy_plan(plan)
     
  endif

  return
end subroutine fft_3D
     
