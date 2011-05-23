#include "../convert.F"

      subroutine eval_tr_sigph_g(rank, tr_sigph_g, sigma, 
     $   g_mtau, t, c_r, d_r, q_tau, q_mtau, r_tau, x, y) 

#include "../constants.F"

#ifdef USE_MPI
      include 'mpif.h'
#endif /* USE_MPI */
      
      INTEGER rank
      COMPLEX tr_sigph_g
      COMPLEX sigma(0:1, 0:1, 0:mp1, 0:nx1, 0:ny1, 0:nz1)
      COMPLEX g_mtau(0:1, 0:1, 0:mp1, 0:nx1, 0:ny1, 0:nz1)
      REAL t

      COMPLEX c_r(0:1,0:1,0:1,0:1,0:nx1,0:ny1,0:nz1)  
      COMPLEX d_r(0:1,0:1,0:3,0:3,0:nx1,0:ny1,0:nz1)
      REAL q_tau(0:1,0:1,0:mp1), q_mtau(0:1,0:1,0:mp1)
      REAL r_tau(0:1,0:1,0:mp1)

      REAL x(0:1,0:1), y(0:1,0:1)

      COMPLEX temp(0:mp1, 0:nx1, 0:ny1, 0:nz1)

      INTEGER l, k1, k2, k3
      REAL weight
      COMPLEX my_sum, sum
      INTEGER dummy

      INTEGER a, ap, a0, b0
      INTEGER i1, j1, i2, j2, i3, j3

#ifdef USE_MPI
      INTEGER stat(MPI_STATUS_SIZE), ierr
#endif /* USE_MPI */

      REAL l_integral(0:1,0:1,0:1,0:1,0:1,0:1)
      INTEGER i

c     Form the integrand which will be traced.

      do l = 0, mp1
        do k1 = 0, nx1
          do k2 = 0, ny1
            do k3 = 0, nz1

              temp(l,k1,k2,k3) = 
     $           sigma(0,0,l,k1,k2,k3) * g_mtau(0,0,l,k1,k2,k3) + 
     $           sigma(0,1,l,k1,k2,k3) * g_mtau(1,0,l,k1,k2,k3) +
     $           sigma(1,0,l,k1,k2,k3) * g_mtau(0,1,l,k1,k2,k3) +
     $           sigma(1,1,l,k1,k2,k3) * g_mtau(1,1,l,k1,k2,k3)

            enddo
          enddo
        enddo
      enddo
     
c     Now subtract the analytic part of the integrand

      do l = 0, mp1
        do k1 = 0, nx1
          do k2 = 0, ny1
            do k3 = 0, nz1

              do a = 0, 1
                do ap = 0, 1
                  do a0 = 0, 1
                    do b0 = 0, 1

                      do i1 = 0, 1
                        do j1 = 0, 1
                          do i2 = 0, 1
                            do j2 = 0, 1
                              do i3 = 0, 1
                                do j3 = 0, 1
                      

                                  temp(l,k1,k2,k3) =
     $                               temp(l,k1,k2,k3) +
     $                               c_r(i1,j1,b0,a0,k1,k2,k3) *
     $                               d_r(i2,j2,2*a+b0,2*ap+a0,
     $                               k1,k2,k3) *
     $                               c_r(i3,j3,ap,a,mod(nx-k1,nx),
     $                               mod(ny-k2,ny), mod(nz-k3,nz)) *
     $                               q_tau(i1,j1,l) * r_tau(i2,j2,l)
     $                               * q_mtau(i3,j3,l)

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
      enddo

      my_sum = 0.0d0

      do l = 0, mp1

        weight = (1.0d0 / t) / float(m)

        if ( mod(l,2) .eq. 0 ) then
          weight = 2.0d0 * weight / 3.0d0
        else  
          weight = 4.0d0 * weight / 3.0d0
        endif

        do k1 = 0, nx1
          do k2 = 0, ny1
            do k3 = 0, nz1
		   
              my_sum = my_sum + temp(l,k1,k2,k3) * weight

            enddo
          enddo
        enddo

      enddo

      if (rank .eq. 0) then

        sum = my_sum
        dummy = 0

#ifdef USE_MPI
        do i = 1, n_proc - 1

          call MPI_Send(dummy, 1, MPI_INTEGER, i, 0, 
     $       MPI_COMM_WORLD, ierr)
          call MPI_Recv(my_sum, 1, MPI_COMPLEX, i, 1,
     $       MPI_COMM_WORLD, stat, ierr)

          sum = sum + my_sum

        enddo

      else 

        call MPI_Recv(dummy, 1, MPI_INTEGER, 0, 0, 
     $     MPI_COMM_WORLD, stat, ierr)
        call MPI_Send(my_sum, 1, MPI_COMPLEX, 0, 1,
     $     MPI_COMM_WORLD, ierr)

#endif /* USE_MPI */
      endif

      if (rank .eq. 0) then

        call l_ints(l_integral, t, x, y)

        my_sum = 0.0

        do k1 = 0, nx1
          do k2 = 0, ny1
            do k3 = 0, nz1

              do a = 0, 1
                do ap = 0, 1
                  do a0 = 0, 1
                    do b0 = 0, 1

                      do i1 = 0, 1
                        do j1 = 0, 1
                          do i2 = 0, 1
                            do j2 = 0, 1
                              do i3 = 0, 1
                                do j3 = 0, 1

                                  my_sum = my_sum +
     $                               c_r(i1,j1,b0,a0,k1,k2,k3) *
     $                               d_r(i2,j2,2*a+b0,2*ap+a0,
     $                               k1,k2,k3) * 
     $                               c_r(i3,j3,ap,a,mod(nx-k1,nx),
     $                               mod(ny-k2,ny), mod(nz-k3,nz))
     $                               * l_integral(i1,j1,i2,j2,i3,j3)

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

        tr_sigph_g = sum - my_sum

      endif
      
      return
      end
