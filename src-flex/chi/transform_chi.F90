#include "../convert.F90"	

subroutine transform_chi(rank, t, chi, delta_chi, delta_chi_prime, &
     r_tau, r_omega,  overall_eigenvalue_max, dominant_chi_eigenvector, & 
     dominant_chi_index)

  USE CONSTANTS
  IMPLICIT NONE

#ifdef USE_MPI
  include 'mpif.h'
#endif

  INTEGER rank
  REAL t
  COMPLEX chi(0:16*nb*nb-1,0:16*nb*nb-1,0:mp1,0:nc1)
  COMPLEX, dimension(0:16*nb*nb-1,0:16*nb*nb-1,0:nc1) :: &
       delta_chi, delta_chi_prime
  REAL r_tau(0:1,0:1,0:mp1)
  COMPLEX r_omega(0:1,0:1,0:mp1)

  COMPLEX d(0:1,0:1,0:16*nb*nb-1,0:16*nb*nb-1,0:nc1)
  INTEGER i, j, l, k, ia, ib
  INTEGER ifftv(0:3)
  COMPLEX sum_dr, sum_mdr
  COMPLEX temp_4D(0:mp1, 0:nc1)
  COMPLEX temp_3D(0:nc1)
  REAL delta_tau 

  CHARACTER*1 jobvl, jobvr
  COMPLEX, dimension (0:16*nb*nb-1, 0:16*nb*nb-1) :: temp_matrix, vl, vr
  COMPLEX eigenvalue(0:16*nb*nb-1)
  COMPLEX work(512*nb*nb)
  REAL rwork(32*nb*nb)
  INTEGER info

  COMPLEX dominant_chi_eigenvector(0:16*nb*nb-1)
  REAL overall_eigenvalue_max
  REAL this_eigenvalue_max
  INTEGER dominant_chi_index(0:1)

#ifdef USE_MPI
  INTEGER dummy
  INTEGER tmp_index(0:1)
  REAL tmp_eigenvalue_max
  COMPLEX tmp_eigenvector(0:16*nb*nb-1)
  INTEGER ierr
  INTEGER stat(MPI_STATUS_SIZE)
#endif /* USE_MPI */

  delta_tau = (1.0d0 / t) / float(m)
  
  d(0,1,:,:,:) = cmplx( 0.0d0, 0.0d0)
  d(0,0,:,:,:) = -delta_chi
  d(1,1,:,:,:) = cmplx( 0.0d0, 0.0d0)
  d(1,0,:,:,:) = delta_chi_prime

  do i = 0, 16*nb*nb-1
     do j = 0, 16*nb*nb-1

        ifftv(0) =  1 
        ifftv(1) =  0
        ifftv(2) =  0
        ifftv(3) =  0

        do l = 0, mp1
           do k = 0, nc1
              sum_dr = cmplx(0.0d0, 0.0d0) 
              do ia = 0, 1
                 do ib = 0, 1
                    sum_dr = sum_dr + r_tau(ia,ib,l) * d(ia,ib,i,j,k)
                 enddo
              enddo
              temp_4D(l,k) = (chi(i,j,l,k) - sum_dr) * delta_tau
           enddo
        enddo

        call fft_4D(rank, temp_4D, ifftv)

        do l = 0, mp1
           do k = 0, nc1
                  
              sum_dr = cmplx(0.0d0, 0.0d0)
              do ia = 0, 1
                 do ib = 0, 1
                    sum_dr = sum_dr + d(ia,ib,i,j,k) * r_omega(ia,ib,l)
                 enddo
              enddo
              temp_4D(l,k) = temp_4D(l,k) + sum_dr

           enddo
        enddo

        ifftv(0) = 0 
        ifftv(1) = -1
        ifftv(2) = -1
        ifftv(3) = -1

        call fft_4D(rank, temp_4D, ifftv)
        chi(i,j,:,:) = temp_4D

     enddo
  enddo

  !     Transform discontinuities
  ifftv(0) = 0
  ifftv(1) = -1
  ifftv(2) = -1
  ifftv(3) = -1

  do i = 0, 16*nb*nb-1
     do j = 0, 16*nb*nb-1
          
        temp_3D = delta_chi(i,j,:)
        call fft_3D( temp_3D, ifftv)
        delta_chi(i,j,:) = temp_3D

        temp_3D = delta_chi_prime(i,j,:)
        call fft_3D( temp_3D, ifftv)
        delta_chi_prime(i,j,:) = temp_3D
        
     enddo
  enddo

  !     Determine the dominant eigenvector of chi and
  !     artificially rescale singular chi(w,k) matrices.
  jobvl = 'N'
  jobvr = 'V' 

  overall_eigenvalue_max = -10.0d0

  do l = 0, mp1
     do k = 0, nc1
              
        this_eigenvalue_max = 0.0d0

        temp_matrix = chi(:,:,l,k)

        call cgeev(jobvl, jobvr, 16*nb*nb, temp_matrix, & 
             16*nb*nb, eigenvalue, vl, 16*nb*nb, vr, & 
             16*nb*nb, work, 512*nb*nb, rwork, info)

        do i = 0, 16*nb*nb-1

           if (real(eigenvalue(i)) .gt. this_eigenvalue_max) then
                  
              this_eigenvalue_max = real(eigenvalue(i))

              if (real(eigenvalue(i)) .gt. overall_eigenvalue_max) then

                 overall_eigenvalue_max = real(eigenvalue(i))
                 do j = 0, 16*nb*nb-1
                    dominant_chi_eigenvector(j) = vr(j,i)
                 enddo
                 dominant_chi_index(0) = l
                 dominant_chi_index(1) = k

              endif
              
           endif

        enddo

#ifdef FLEX
        if (this_eigenvalue_max .gt. 0.9999999d0) then

           chi(:,:,l,k) = chi(:,:,l,k) * 0.99995d0 / this_eigenvalue_max

        endif
#endif /* FLEX */

     enddo
  enddo
      
#ifdef USE_MPI
  if (rank .eq. 0) then

     do i = 1, np-1

        dummy = 0
        call MPI_Send(dummy, 1, MPI_INTEGER, i, 0, &
             MPI_COMM_WORLD, ierr)
        call MPI_Recv(tmp_eigenvalue_max, 1, MPI_REAL, & 
             i, 1, MPI_COMM_WORLD, stat, ierr)
        call MPI_Send(dummy, 1, MPI_INTEGER, i, 2, &
             MPI_COMM_WORLD, ierr)  
        call MPI_Recv(tmp_eigenvector, 16*nb*nb, & 
             MPI_COMPLEX, i, 3, MPI_COMM_WORLD, stat, ierr)
        call MPI_Send(dummy, 1, MPI_INTEGER, i, 4, MPI_COMM_WORLD, ierr)  
        call MPI_Recv(tmp_index, 2, MPI_INTEGER, i, 5, MPI_COMM_WORLD, &
             stat, ierr)

        if (tmp_eigenvalue_max .gt. overall_eigenvalue_max) then

           overall_eigenvalue_max = tmp_eigenvalue_max

           do j = 0, 16*nb*nb-1
              dominant_chi_eigenvector(j) = tmp_eigenvector(j)
           enddo

           dominant_chi_index(0) = tmp_index(0) + mp * i
           dominant_chi_index(1) = tmp_index(1)

        endif

     enddo

  else

     call MPI_Recv(dummy, 1, MPI_INTEGER, 0, 0, &
          MPI_COMM_WORLD, stat, ierr)
     call MPI_Send(overall_eigenvalue_max, 1, & 
          MPI_REAL, 0, 1, MPI_COMM_WORLD, ierr)
     call MPI_Recv(dummy, 1, MPI_INTEGER, 0, 2, &
          MPI_COMM_WORLD, stat, ierr)          
     call MPI_Send(dominant_chi_eigenvector, 16*nb*nb, & 
          MPI_COMPLEX, 0, 3, MPI_COMM_WORLD, ierr)
     call MPI_Recv(dummy, 1, MPI_INTEGER, 0, 4, &
          MPI_COMM_WORLD, stat, ierr)           
     call MPI_Send(dominant_chi_index, 2, MPI_INTEGER, & 
          0, 5, MPI_COMM_WORLD, ierr)

  endif

#endif /* USE_MPI */

  return
end subroutine transform_chi




