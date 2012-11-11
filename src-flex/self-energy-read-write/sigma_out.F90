#include "../convert.F90"

subroutine sigma_out(rank, pade_max, sigma_output_file, t, sigma)

  USE CONSTANTS
  USE pade_eval
  IMPLICIT NONE

#ifdef USE_MPI
  include 'mpif.h'
#endif /* USE_MPI */

  INTEGER rank
  INTEGER pade_max
  CHARACTER*128 sigma_output_file
  REAL t
  COMPLEX sigma(0:4*nb-1,0:4*nb-1,0:mp1,0:nc1)
  
  REAL pi
  COMPLEX a(0:m/2-1), z(0:m/2-1)
  COMPLEX p_temp(0:m/2-1), q_temp(0:m/2-1)
  COMPLEX p(0:ncp1,0:m/2-1), q(0:ncp1,0:m/2-1)

  INTEGER k_count, n_pade(0:ncp1)

  INTEGER i, j, l, k
  INTEGER nu1, ia1, nu2, ia2
  INTEGER ind1, ind2

  COMPLEX freq
  INTEGER i_pade
  LOGICAL causal, previous_causal, causal_help

  COMPLEX temprow(0:mp1,0:nc1)
  COMPLEX tempcol(0:m1,0:ncp1)

  COMPLEX temprow_00(0:mp1,0:nc1), temprow_11(0:mp1,0:nc1)

  COMPLEX tempcol_00(0:m1,0:ncp1), tempcol_11(0:m1,0:ncp1)

  INTEGER dummy, i_tag
  INTEGER i_proc, ierr

#ifdef USE_MPI
  INTEGER stat(MPI_STATUS_SIZE)
#endif /* USE_MPI */

  if (rank .eq. 0) then
     !       open(unit=88,file=sigma_output_file,status='old',
#ifdef f90_flag
     !    $     access='sequential',position='append')
#else
     !    $     access='append')
#endif
     write(40,*) '-------Sigma (n >= 2)------------ '
     write(40,*) 'i  j  site_index   pade_order : followed by coefficients'
  endif

  pi = 4.0d0 * atan(1.0d0)

  do l = 0, m/2 - 1
     z(l) = cmplx( 0.0d0, pi * t * float(2*l + 1) )
  enddo

  i_pade = (m/4) - 1        ! Maximum possible Pade order

  if ( i_pade > pade_max + 1) then
     i_pade = pade_max +1
  endif
      
  previous_causal = .false.
  causal = .false.

  do while ( (i_pade .ge. 0) .and. &
       ( (.not. previous_causal) .or. (.not. causal) ) )

     previous_causal = causal
     causal = .true.
        
     temprow_00 = sigma(0,0,:,:)
     temprow_11 = sigma(1,1,:,:)

     call row_dist_to_col_dist(rank, temprow_00, tempcol_00)
     call row_dist_to_col_dist(rank, temprow_11, tempcol_11) 

     do k = 0, ncp1

        do l = 0, m/2 -1 
           a(l) = tempcol_00(l,k)
        enddo

        call pade(i_pade, a, z, p_temp, q_temp)

        freq = cmplx(-20.0d0, 0.1d0)

        do while ( real(freq) < 20.0d0 ) 
           if ( imag( pade_evaluate(freq, i_pade, p_temp, & 
                q_temp) ) .gt. 0.0d0 ) then
              causal = .false.
           endif
           freq = freq + 0.1d0
        enddo

        do l = 0, m/2 - 1
           a(l) = tempcol_11(l,k)
        enddo
              
        call pade(i_pade, a, z, p_temp, q_temp)
        freq = cmplx(-20.0d0, 0.1d0)

        do while ( real(freq) < 20.0 )
              
           if ( imag(pade_evaluate(freq, i_pade, p_temp, &
                q_temp) ) .gt. 0.0d0 ) then
              causal = .false.
           endif
           freq = freq + 0.1d0
            
        enddo
          
     enddo

#ifdef USE_MPI
     if (rank .eq. 0) then
        do i_proc = 1, np - 1
           call MPI_Send(dummy, 1, MPI_INTEGER, i_proc, 500+i_pade, &
                MPI_COMM_WORLD, ierr)
            call MPI_Recv(causal_help, 1, MPI_LOGICAL, i_proc, & 
                 500+i_pade, MPI_COMM_WORLD, stat, ierr)
            causal = causal_help .and. causal
         enddo
      else 
         call MPI_Recv(dummy, 1, MPI_INTEGER, 0, 500+i_pade, &
              MPI_COMM_WORLD, stat, ierr)
         call MPI_Send(causal, 1, MPI_LOGICAL, 0, 500+i_pade, & 
              MPI_COMM_WORLD, ierr)
      endif

      call MPI_Bcast(causal, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr) 
#endif /* USE_MPI */
 
      i_pade = i_pade - 1

   enddo

   i_pade = i_pade + 1

   if ( i_pade .lt. pade_max ) then
      if (i_pade .le. 4) then
         pade_max = 5   
      else 
         pade_max = i_pade
      endif
   endif

   if (rank .eq. 0) then
      write(6,*)
      write(6,*)
      write(6,*) 'Optimal number of pade coefficients = ', i_pade
      write(6,*) 'Number of pade coefficients used = ', pade_max
      write(6,*)
   endif

   !     Obtain and print the pade coefficients.
      
   do nu1 = 0, nb-1
      do ia1 = 0, 3
         ind1 = 4*nu1 + ia1

         do nu2 = 0, nb-1
            do ia2 = 0, 3
               ind2 = 4*nu2 + ia2
               
               temprow = sigma(ind1,ind2,:,:)
              
               call row_dist_to_col_dist(rank, temprow, tempcol)

               do k = 0, ncp1
                
                  k_count = 0

                  !     Determine how many coefficients to use for off-diagonal
                  !     parts of sigma.  (diagonal determined by causality)

                  if (ind1 .ne. ind2) then

                     do l = 0, m/2 - 1
                        if (cabs(tempcol(l,k)) .gt. 1.0d-10) then
                           k_count = k_count + 1
                        endif
                     enddo

                     if (k_count .lt. 2*pade_max + 2) then
                        i_pade = (k_count/2) - 1
                     else 
                        i_pade = pade_max
                     endif
                  
                  else

                     i_pade = pade_max

                  endif

                  if ( i_pade .gt. 0) then
                     
                     do l = 0, 2*i_pade + 1
                        a(l) = tempcol(l,k)
                     enddo

                     call pade(i_pade, a, z, p_temp, q_temp)
                  
                     do l = 0, i_pade
                        q(k,l) = q_temp(l)
                        p(k,l) = p_temp(l)
                     enddo

                  endif
                  
                  n_pade(k) = i_pade

               enddo

               if (rank .eq. 0) then
	
                  do k = 0, ncp1

                     write(40,300)  nu1, ia1, nu2,ia2, k, n_pade(k)
                  
                     if (n_pade(k) > 0) then
                        do l = 0, n_pade(k)
                           write(40,*)  p(k,l), ' ', q(k,l)
                        enddo
                     endif
                     
                  enddo

               endif

#ifdef USE_MPI
               if (np .gt. 1) then

                  if (rank .eq. 0) then

                     do i_proc = 1, np - 1

                        i_tag = 1
                        dummy = 0
                      
                        call MPI_Send(dummy, 1, MPI_INTEGER, & 
                             i_proc, i_tag, MPI_COMM_WORLD, ierr)

                        call MPI_Recv(n_pade, ncp, MPI_INTEGER, & 
                             i_proc, i_tag, MPI_COMM_WORLD, stat, ierr) 
              
                        i_tag = 2
                        call MPI_Recv(p, ncp*m/2, MPI_COMPLEX, & 
                             i_proc, i_tag, MPI_COMM_WORLD, stat, ierr)
                      
                        i_tag = 3
                        call MPI_Recv(q, ncp*m/2, MPI_COMPLEX, & 
                             i_proc, i_tag, MPI_COMM_WORLD, stat, ierr)

                        do k = 0, ncp1
                        
                           write(40,300) nu1,ia1,nu2,ia2,&
                                k + i_proc*ncp, n_pade(k)
                        
                           if (n_pade(k) .gt. 0) then

                              do l = 0, n_pade(k)
                                 write(40,*) p(k,l), ' ', q(k,l)
                              enddo
                    
                           endif

                        enddo

                     enddo

                  else

                     i_tag = 1
                  
                     call MPI_Recv(dummy, 1, MPI_INTEGER, 0, i_tag, &
                          MPI_COMM_WORLD, stat, ierr)
                  
                     call MPI_Send(n_pade, ncp, MPI_INTEGER, 0, i_tag, &
                          MPI_COMM_WORLD, ierr) 

                     i_tag = 2
                     call MPI_Send(p, ncp*m/2, MPI_COMPLEX, 0, i_tag, & 
                          MPI_COMM_WORLD, ierr)
             
                     i_tag = 3
                     call MPI_Send(q, ncp*m/2, MPI_COMPLEX, 0, i_tag, &
                          MPI_COMM_WORLD, ierr)

                  endif

               endif
#endif /* USE_MPI */

            enddo
         enddo
      enddo
   enddo

   if (rank .eq. 0) then
      !       close(unit=88)
   endif
   
300 format(i3,2x,i2,2x,i3,2x,i2, 2x, i6, 2x,i4)

   return
 end subroutine sigma_out
