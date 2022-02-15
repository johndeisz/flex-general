#include "../convert.F90"

subroutine fft_4D (myid, a, isignv)

  USE CONSTANTS
  IMPLICIT NONE

  INTEGER myid
  COMPLEX a(0:mp1,0:nc1)
  INTEGER isignv(0:3)
      
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !  constants needed for the fortran wrapper functions!
#include "../fftw3.f"

  COMPLEX temp_row(0:mp1,0:nc1)
  COMPLEX temp_col(0:m1,0:ncp1)

  COMPLEX temp_t_in(0:m1), temp_t_out(0:m1)
  COMPLEX temp_x_in(0:lcx1), temp_x_out(0:lcx1)
  COMPLEX temp_y_in(0:lcy1), temp_y_out(0:lcy1)
  COMPLEX temp_z_in(0:lcz1), temp_z_out(0:lcz1)

  INTEGER l, j, j1, j2, j3

  INTEGER fftw_plan
  INTEGER*8 plan

  !     FFT with respect to time or frequency.

  if (isignv(0) .ne. 0) then

     call sfftw_plan_dft_1d(plan,m,temp_t_in,temp_t_out, isignv(0), &
          FFTW_ESTIMATE)

     temp_row = a
     call row_dist_to_col_dist(myid, temp_row, temp_col)

     do j = 0, ncp1
        temp_t_in = temp_col(:,j)
        call sfftw_execute(plan,temp_t_in,temp_t_out)
        temp_col(:,j) = temp_t_out
     enddo

     call col_dist_to_row_dist(myid, temp_col, temp_row)
     a = temp_row

     call sfftw_destroy_plan(plan)

  endif

  !     do x-axis dft
        
  if ( (lcx .gt. 1) .and. (isignv(1) .ne. 0) ) then

     call sfftw_plan_dft_1d(plan,lcx,temp_x_in,temp_x_out, &
          isignv(1),FFTW_ESTIMATE)

     do l = 0, mp1
        do j2 = 0, lcy1
           do j3 = 0, lcz1

              do j1 = 0, lcx1
                 temp_x_in(j1) = a(l, j1+j2*lcx+j3*lcx*lcy)
              enddo

              call sfftw_execute(plan,temp_x_in, temp_x_out)

              do j1 = 0, lcx1
                 a(l, j1+j2*lcx+j3*lcx*lcy) = temp_x_out(j1)
              enddo

           enddo
        enddo
     enddo
     
     call sfftw_destroy_plan(plan)             
     
  endif

  !     do y dft

  if ( (lcy .gt. 1) .and. (isignv(2) .ne. 0) ) then

     call sfftw_plan_dft_1d(plan,lcy,temp_y_in,temp_y_out, &
          isignv(2),FFTW_ESTIMATE)
       
     do l = 0, mp1
        do j1 = 0, lcx1
           do j3 = 0, lcz1

              do j2 = 0, lcy1
                 temp_y_in(j2) = a(l,j1+j2*lcx+j3*lcx*lcy)
              enddo
              call sfftw_execute(plan,temp_y_in,temp_y_out)

              do j2 = 0, lcy1
                 a(l,j1+j2*lcx+j3*lcx*lcy) = temp_y_out(j2)
              enddo
            
           enddo
        enddo
     enddo

     call sfftw_destroy_plan(plan)
     
  endif                     ! (ny > 1 .and. isignv(2) != 0)

  !     do z dft

  if ( (lcz .gt. 1) .and. (isignv(3) .ne. 0) ) then

     call sfftw_plan_dft_1d(plan,lcz,temp_z_in,temp_z_out, &
          isignv(3),FFTW_ESTIMATE)

     do l = 0, mp1
        do j1 = 0, lcx1
           do j2 = 0, lcy1

              do j3 = 0, lcz1
                 temp_z_in(j3) = a(l,j1+j2*lcx+j3*lcx*lcy)
              enddo
              call sfftw_execute(plan,temp_z_in,temp_z_out)

              do j3 = 0, lcz1
                 a(l,j1+j2*lcx+j3*lcx*lcy) = temp_z_out(j3)
              enddo

           enddo
        enddo
     enddo

     call sfftw_destroy_plan(plan)

  endif

  return
end subroutine fft_4D
