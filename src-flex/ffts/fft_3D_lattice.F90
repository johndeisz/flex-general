#include "../convert.F90"

subroutine fft_3D_lattice (mat, isignv)

  USE CONSTANTS
  IMPLICIT NONE

  COMPLEX mat(0:nl-1)
  INTEGER isignv(0:3)

#include "../fftw3.f"

  !     isignv controls whether or not an fft is done on a 
  !     given axis and the direction of the fft.   
  !     nu-k -> tau-r  is (1, -1, -1, -1)
  !     no scaling is done -- this is done by hand in the flux code

      
  COMPLEX matx_in(0:llx-1), matx_out(0:llx-1)
  COMPLEX maty_in(0:lly-1), maty_out(0:lly-1)
  COMPLEX matz_in(0:llz-1), matz_out(0:llz-1)
  
  INTEGER j1, j2, j3
  INTEGER fftw_plan
  INTEGER*8 plan  

  !     Do x-axis Fourier transform.
        
  if ( (llx .gt. 1) .and. (isignv(1) .ne. 0) ) then
     
     call sfftw_plan_dft_1d(plan, llx, matx_in, matx_out, &
          isignv(1), FFTW_ESTIMATE)
     
     do j2 = 0, lly-1
        do j3 = 0, llz-1
           
           do j1 = 0, llx-1
              matx_in(j1) = mat(j1 + j2*llx + j3*llx*lly)
           enddo

           call sfftw_execute(plan,matx_in, matx_out)

           do j1 = 0, llx-1 
              mat(j1 + j2*llx + j3*llx*lly) = matx_out(j1)
           enddo

        enddo
     enddo

     call sfftw_destroy_plan(plan)

  endif

  if ( (lly .gt. 1) .and. (isignv(2) .ne. 0) ) then
     
     call sfftw_plan_dft_1d(plan, lly,maty_in,maty_out,isignv(2), &
          FFTW_ESTIMATE)

     do j1 = 0, llx - 1
        do j3 = 0, llz - 1
           
           do j2 = 0, lly - 1
              maty_in(j2) = mat(j1 + j2*llx + j3*llx*lly)
           enddo

           call sfftw_execute(plan,maty_in,maty_out)
           
           do j2 = 0, lly - 1
              mat(j1 + j2*llx + j3*llx*lly) = maty_out(j2)
           enddo
            
        enddo
     enddo

     call sfftw_destroy_plan(plan)
     
  endif ! (lly > 1 .and. isignv(2) != 0)
  
  if ( (llz .gt. 1) .and. (isignv(3) .ne. 0) ) then

     call sfftw_plan_dft_1d(plan, llz, matz_in,matz_out, isignv(3), &
          FFTW_ESTIMATE)

     do j1 = 0, llx - 1
        do j2 = 0, lly - 1

           do j3 = 0, llz - 1
              matz_in(j3) = mat(j1 + j2*llx + j3*llx*lly)
           enddo

           call sfftw_execute(plan,matz_in,matz_out)

           do j3 = 0, llz - 1
              mat(j1 + j2*llx + j3*llx*lly) = matz_out(j3)
           enddo

        enddo
     enddo

     call sfftw_destroy_plan(plan)            

  endif

  return
end subroutine fft_3D_lattice
      
