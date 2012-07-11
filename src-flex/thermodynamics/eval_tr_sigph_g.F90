#include "../convert.F90"

subroutine eval_tr_sigph_g(rank, tr_sig_g, sigma, &
     g_mtau, t, c_r, d_r, q_tau, q_mtau, r_tau, x, y) 

  USE CONSTANTS

#ifdef USE_MPI
  include 'mpif.h'
#endif /* USE_MPI */
      
  INTEGER rank
  COMPLEX tr_sig_g
  COMPLEX, dimension (0:4*nb-1,0:4*nb-1,0:mp-1,0:nc1) :: sigma, g_mtau
  REAL t

  COMPLEX c_r(0:1,0:1,0:4*nb-1,0:4*nb-1,0:nc1)
  COMPLEX d_r(0:1,0:1,0:16*nb*nb-1,0:16*nb*nb-1,0:nc1)
  REAL q_tau(0:1,0:1,0:mp1), q_mtau(0:1,0:1,0:mp1)
  REAL r_tau(0:1,0:1,0:mp1)

  REAL x(0:1,0:1), y(0:1,0:1)

  COMPLEX temp

  INTEGER l, ix, iy, iz, ir, irm
  INTEGER nua, nuap, a0, b0
  REAL weight
  COMPLEX my_sum, sum
  INTEGER dummy

  INTEGER i1, j1, i2, j2, i3, j3

#ifdef USE_MPI
  INTEGER stat(MPI_STATUS_SIZE), ierr
#endif /* USE_MPI */

  REAL l_integral(0:1,0:1,0:1,0:1,0:1,0:1)
  INTEGER i

  !     Form the integrand which will be traced.

  my_sum = 0.0d0

  do nua = 0, 4*nb-1
     do nuap = 0, 4*nb-1

        do ix = 0, lcx1
           do iy = 0, lcy1
              do iz = 0, lcz1

                 ir = ix + iy*lcx + iz*lcx*lcy
                 irm = mod(-ix+lcx,lcx) + mod(-iy+lcy,lcy)*lcx + &
                      mod(-iz+lcz,lcz)*lcx*lcy

                 do l = 0, mp1

                    temp = sigma(nua,nuap,l,ir)*g_mtau(nuap,nua,l,ir) 
                  
                    !     Now subtract the analytic part of the integrand

                    do a0 = 0, 4*nb-1
                       do b0 = 0, 4*nb-1

                          do i1 = 0, 1
                             do j1 = 0, 1
                                do i2 = 0, 1
                                   do j2 = 0, 1
                                      do i3 = 0, 1
                                         do j3 = 0, 1
                              
                                            
                                            temp = temp + &
                                                 c_r(i1,j1,b0,a0,ir) * &
                                                 d_r(i2,j2,4*nb*nua+b0, &
                                                 4*nb*nuap+a0,ir) * &
                                                 c_r(i3,j3,nuap,nua,irm) * &
                                                 q_tau(i1,j1,l) * &
                                                 r_tau(i2,j2,l) * &
                                                 q_mtau(i3,j3,l)

                                         enddo
                                      enddo
                                   enddo
                                enddo
                             enddo
                          enddo

                       enddo
                    enddo

                    weight = (1.0d0 / t) / float(m)

                    if ( mod(l,2) .eq. 0 ) then
                       weight = 2.0d0 * weight / 3.0d0
                    else  
                       weight = 4.0d0 * weight / 3.0d0
                    endif

                    my_sum = my_sum + temp * weight

                 enddo

              enddo
           enddo
        enddo
     enddo
  enddo

  if (rank .eq. 0) then

     sum = my_sum
     dummy = 0

#ifdef USE_MPI
     do i = 1, np - 1

        call MPI_Send(dummy, 1, MPI_INTEGER, i, 0,  MPI_COMM_WORLD, ierr)
        call MPI_Recv(my_sum, 1, MPI_COMPLEX, i, 1, &
             MPI_COMM_WORLD, stat, ierr)

        sum = sum + my_sum

     enddo

  else 

     call MPI_Recv(dummy, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, stat, ierr)
     call MPI_Send(my_sum, 1, MPI_COMPLEX, 0, 1, MPI_COMM_WORLD, ierr)

#endif /* USE_MPI */
  endif

  if (rank .eq. 0) then

     call l_ints(l_integral, t, x, y)

     my_sum = 0.0

     do ix = 0, lcx1
        do iy = 0, lcy1
           do iz = 0, lcz1

              ir = ix + iy*lcx + iz*lcx*lcy
              irm = mod(-ix+lcx,lcx) + mod(-iy+lcy,lcy)*lcx + &
                   mod(-iz+lcz,lcz)*lcx*lcy
              
              do nua = 0, 4*nb-1
                 do nuap = 0, 4*nb-1
                    do a0 = 0, 4*nb-1
                       do b0 = 0, 4*nb-1
                      
                          do i1 = 0, 1
                             do j1 = 0, 1
                                do i2 = 0, 1
                                   do j2 = 0, 1
                                      do i3 = 0, 1
                                         do j3 = 0, 1
                                  
                                            my_sum = my_sum + &
                                                 c_r(i1,j1,b0,a0,ir) * &
                                                 d_r(i2,j2,4*nb*nua+b0, &
                                                 4*nb*nuap+a0,ir) * &
                                                 c_r(i3,j3,nuap,nua,irm) * &
                                                 l_integral(i1,j1,i2,j2,i3,j3)

                                         enddo
                                      enddo
                                   enddo
                                enddo
                             enddo
                          enddo

                       enddo
                    enddo
                 enddo
              enddo

           enddo
        enddo
     enddo

     tr_sig_g = sum - my_sum

  endif
      
  return
end subroutine eval_tr_sigph_g
