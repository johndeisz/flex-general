#include "../convert.F"

      subroutine init_pair_wave(psi)

#include "../constants.F"

#ifdef USE_MPI
      include 'mpif.h'
#endif /* USE_MPI */

      COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)

      INTEGER rank, ierr
      INTEGER i, j, ix, iy, iz
      INTEGER ib, is, jb, js, il
      INTEGER ind1, ind2
      
      REAL real_psi, imag_psi

      COMPLEX tmp_psi
      REAL psi_norm

#ifdef USE_MPI
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
#else
      rank = 0
#endif /* USE_MPI */
     
      if (rank .eq. 0) then

        do i = 0, 2*nb-1
          do j = 0, 2*nb-1
            do il = 0, nl-1
              real_psi = rand() - 0.5d0
              imag_psi = rand() - 0.5d0
	    
              psi(i,j,il) = cmplx(real_psi, imag_psi)

            enddo
          enddo
        enddo

c     Antisymmetrize
        do ib = 1, nb
          do is = 0, 1
            do jb = 1, nb
              do js = 0, 1
                do ix = 0, llx1
                  do iy = 0, lly1
                    do iz = 0, llz1

                      ind1 = ix + iy*llx + iz*llx*lly
                      ind2 = mod(llx-ix,llx) +
     $                   mod(lly-iy,lly)*llx +
     $                   mod(llz-iz,llz)*llx*lly

                      tmp_psi = 0.5d0 *  
     $                   ( psi(2*(ib-1)+is,2*(jb-1)+js,ind1) - 
     $                   psi( 2*(jb-1)+js,2*(ib-1)+is,ind2) )

                      psi(2*(ib-1)+is,2*(jb-1)+js, ind1) = 
     $                   tmp_psi

                      psi(2*(jb-1)+js,2*(ib-1)+is, ind2) = 
     $                   -tmp_psi
                      
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo

        psi_norm = 0.0d0

        do i = 0, 2*nb-1
          do j = 0, 2*nb-1
            do il = 0, nl-1

              psi_norm = psi_norm + 
     $           (real(psi(i,j,il)))**2 +
     $           (imag(psi(i,j,il)))**2

            enddo
          enddo
        enddo

        psi = psi / sqrt(psi_norm)

      endif

#ifdef USE_MPI
      call MPI_Bcast(psi, 4*nb*nb*nl,
     $   MPI_COMPLEX, 0, MPI_COMM_WORLD, ierr)
#endif /* USE_MPI */

      return
      end
