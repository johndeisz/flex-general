#include "../convert.F90"

#ifdef SECOND_ORDER
subroutine sigma_input(sigma_input_file, uu, sigma1, psi, sigma, mu, epsilon)
#else
  subroutine sigma_input(sigma_input_file, uu, sigma1, psi) 
#endif

    USE CONSTANTS
#ifdef SECOND_ORDER
    USE pade_eval
#endif
    IMPLICIT NONE

#ifdef USE_MPI
    include 'mpif.h'
#endif

    CHARACTER*128 sigma_input_file
    REAL uu
    COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)
    COMPLEX sigma1(0:4*nb-1, 0:4*nb-1)
#ifdef SECOND_ORDER
    COMPLEX sigma(0:4*nb-1,0:4*nb-1, 0:mp1, 0:nc-1)
    REAL mu
    REAL epsilon(0:mp1)
#endif 

    INTEGER lcx_in, lcy_in, lcz_in
    INTEGER llx_in, lly_in, llz_in
    INTEGER pair_flag, so_flag
    INTEGER id1, id2, id3, id4, id5

    REAL uu_1, up_1, uj_1

    LOGICAL x_stretch, x_compress
    INTEGER x_stretch_factor, x_compress_factor
    LOGICAL y_stretch, y_compress
    INTEGER y_stretch_factor, y_compress_factor
    LOGICAL z_stretch, z_compress
    INTEGER z_stretch_factor, z_compress_factor

    INTEGER ix(0:2), iy(0:2), iz(0:2)
    REAL wx(0:2), wy(0:2), wz(0:2)
    COMPLEX psi_temp, sigma_temp
    REAL re_part, im_part

    INTEGER ii

#ifdef SECOND_ORDER
    COMPLEX p(0:m/2-1), q(0:m/2-1)
    INTEGER n_pade
    COMPLEX z
    COMPLEX temp_pade
#endif /* SECOND_ORDER */

    INTEGER rank

    INTEGER k1, k2, k3, k, kp
    INTEGER i1, i2, i3
    INTEGER nu1, nu2, is1, is2, ia1, ia2, ind1, ind2
    INTEGER l
    INTEGER ib, ibp
    INTEGER iix, iiy, iiz
    INTEGER max_x, max_y, max_z

#ifdef USE_MPI
    INTEGER ierr
#endif /* USE_MPI */


#ifdef USE_MPI
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
#else
    rank = 0
#endif /* USE_MPI */

    if (rank .eq. 0) then
       open(unit=50, file=sigma_input_file, status='old')
       read(50,*)
       read(50,*)
       read(50,*)
       read(50,*)
       read(50,*)
       read(50,*) ! nb
       read(50,*)
       read(50,*) lcx_in, lcy_in, lcz_in, id1, llx_in, lly_in, llz_in
       read(50,*)  ! string
       read(50,*)  ! temp, mu, flux_x, flux_y, flux_z
       read(50,*)  ! string
       read(50,*)  ! prfld
       read(50,*)  ! string
       do ib = 0, nb - 1
          read(50,*) ! ed(ib) 
       enddo
       read(50,*)  ! string
       do ib = 0, nb - 1
          read(50,*)  ! h(ib,:)
       enddo
       read(50,*) ! string
       read(50,*) uu_1, up_1, uj_1
       read(50,*) ! string

       if (llx_in .gt. 2) then
          max_x = 2
       else
          max_x = llx_in - 1
       endif

       if (lly_in .gt. 2) then
          max_y = 2
       else
          max_y = lly_in - 1
       endif

       if (llz_in .gt. 2) then
          max_z = 2
       else
          max_z = llz_in - 1
       endif

       do iix = -max_x, max_x
          do iiy = -max_y, max_y
             do iiz = -max_z, max_z

                read(50,*)  ! ix, iy, iz

                do ib = 0, nb-1
                   do ibp = 0, nb-1
                      read(50,*) ! ib, ibp, tij
                   enddo
                enddo

             enddo
          enddo
       enddo

       read(50,*)
    endif

#ifdef USE_MPI
    call MPI_Bcast(uu_1, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(lcx_in, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(lcy_in, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(lcz_in, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(llx_in, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(lly_in, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(llz_in, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif /* USE_MPI */

    !     Determine how to stretch/compress the input pair wave function
    !     if the lattice size for the input file
    !     and the present calculation do not match.

    if (llx > llx_in) then
       x_stretch = .true.
       x_stretch_factor = llx / llx_in
       x_compress = .false.
       x_compress_factor = 1
    else
       if (llx < llx_in) then
          x_stretch = .false.
          x_stretch_factor = 1
          x_compress = .true.
          x_compress_factor = llx_in / llx
       else 
          x_stretch = .false.
          x_stretch_factor = 1
          x_compress = .false.
          x_compress_factor = 1
       endif
    endif

    if (lly > lly_in) then
       y_stretch = .true.
       y_stretch_factor = lly / lly_in
       y_compress = .false.
       y_compress_factor = 1
    else 
       if (lly < lly_in) then
          y_stretch = .false.
          y_stretch_factor = 1
          y_compress = .true.
          y_compress_factor = lly_in / lly
       else 
          y_stretch = .false.
          y_stretch_factor = 1
          y_compress = .false.
          y_compress_factor = 1
       endif
    endif

    if (llz > llz_in) then
       z_stretch = .true.
       z_stretch_factor = llz / llz_in
       z_compress = .false.
       z_compress_factor = 1
    else 
       if (llz < llz_in) then
          z_stretch = .false.
          z_stretch_factor = 1
          z_compress = .true.
          z_compress_factor = llz_in / llz
       else 
          z_stretch = .false.
          z_stretch_factor = 1
          z_compress = .false.
          z_compress_factor = 1
       endif
    endif

    !     Pair wave function stuff
    if (rank .eq. 0) then
       read(50,*)
       read(50,*)
       read(50,*) pair_flag

       if (pair_flag > 0) then

          read(50,*) 

          do nu1 = 0, nb-1
             do is1 = 0, 1
                do nu2 = 0, nb-1
                   do is2 = 0, 1

                      do k3 = 0, min(llz, llz_in) - 1
                         do k2 = 0, min(lly, lly_in) - 1
                            do k1 = 0, min(llx, llx_in) - 1

                               k = k1*x_stretch_factor + &
                                    k2*y_stretch_factor*llx + &
                                    k3*z_stretch_factor*llx*lly

                               read(50,*)  id1, id2, id3, id4, id5, & 
                                    psi(2*nu1+is1, 2*nu2+is2,k) 

                               do i1 = 1, x_compress_factor - 1
                                  read(50,*)
                               enddo

                            enddo

                            do i2 = 1, y_compress_factor - 1
                               do i1 = 0, llx_in - 1
                                  read(50,*)
                               enddo
                            enddo

                         enddo

                         do i3 = 1, z_compress_factor-1
                            do i2 = 0, lly_in - 1
                               do i1 = 0, llx_in - 1
                                  read(50, *)
                               enddo
                            enddo
                         enddo

                      enddo

                   enddo
                enddo
             enddo
          enddo

          !     Stretch the pair wave function.
          
          do nu1 = 0, nb-1
             do is1 = 0, 1
                do nu2 = 0, nb-1
                   do is2 = 0, 1

                      do k3 = 0, llz-1
           
                         iz(0) = k3 - mod(k3,z_stretch_factor)
                         iz(1) = mod(iz(0) + z_stretch_factor,llz) 
           
                         wz(0) = 1.0d0 - float(k3-iz(0)) / & 
                              float(z_stretch_factor)
                         wz(1) = 1.0d0 - wz(0)

                         do k2 = 0, lly-1

                            iy(0) = k2 - mod(k2, y_stretch_factor)
                            iy(1) = mod(iy(0) + y_stretch_factor, lly)
                
                            wy(0) = 1.0d0 - float( k2 - iy(0) ) / &
                                 float(y_stretch_factor)
                            wy(1) = 1.0d0 - wy(0)

                            do k1 = 0, llx-1

                               ix(0) = k1 - mod(k1, x_stretch_factor)
                               ix(1) = mod(ix(0) + x_stretch_factor, llx)

                               wx(0) = 1.0d0 - float( k1 - ix(0) ) / & 
                                    float(x_stretch_factor)
                               wx(1) = 1.0d0 - wx(0)

                               psi_temp = cmplx(0.0d0,0.0d0)

                               do i1 = 0, 1
                                  do i2 = 0, 1
                                     do i3 = 0, 1
                              
                                        kp = ix(i1) + iy(i2)*llx + &
                                             iz(i3)*llx*lly

                                        psi_temp = psi_temp + &
                                             wx(i1) * wy(i2) * wz(i3) * &
                                             psi(2*nu1+is1,2*nu2+is2,kp)
                          
                                     enddo
                                  enddo
                               enddo

                               k = k1 + k2*llx + k3*llx*lly

                               psi(2*nu1+is1,2*nu2+is2,k) = psi_temp

                            enddo
                         enddo
                      enddo

                   enddo
                enddo
             enddo
          enddo

        endif

        !     Spin-orbit matrix stuff
        read(50,*)
        read(50,*)
        read(50,*)
        read(50,*) so_flag

        if (so_flag > 0) then
           read(50,*)
           do nu1 = 0, nb-1
              do is1 = 0, 1
                 do nu2 = 0, nb-1
                    do is2 = 0, 1
                       read(50,*)
                    enddo
                 enddo
              enddo
           enddo
        endif

        read(50,*)
        read(50,*)
        read(50,*)

        ! Read and scale the Hartree self-energy
        do nu1 = 0, nb-1
           do ia1 = 0, 3
              do nu2 = 0, nb-1
                 do ia2 = 0, 3

                    read(50,*) id1, id2, id3, id4, sigma1(4*nu1+ia1, 4*nu2+ia2)
                
                 enddo
              enddo
           enddo
        enddo

        if ( abs(uu_1) .gt. 1.0d-4 ) then
          sigma1 = (uu / uu_1) * sigma1
        endif 

     endif

#ifdef USE_MPI
     call MPI_Bcast(psi, 4*nb*nb*nl, MPI_COMPLEX, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(sigma1, 16*nb*nb, MPI_COMPLEX, 0, MPI_COMM_WORLD, ierr)
#endif /* USE_MPI */      

     !     Determine how to stretch/compress the input sigma
     !     if the lattice size for the input file
     !     and the present calculation do not match.

#ifdef SECOND_ORDER

     if (lcx > lcx_in) then
        x_stretch = .true.
        x_stretch_factor = lcx / lcx_in
        x_compress = .false.
        x_compress_factor = 1
     else
        if (lcx < lcx_in) then
           x_stretch = .false.
           x_stretch_factor = 1
           x_compress = .true.
           x_compress_factor = lcx_in / lcx
        else
           x_stretch = .false.
           x_stretch_factor = 1
           x_compress = .false.
           x_compress_factor = 1
        endif
     endif

     if (lcy > lcy_in) then
        y_stretch = .true.
        y_stretch_factor = lcy / lcy_in
        y_compress = .false.
        y_compress_factor = 1
     else
        if (lcy < lcy_in) then
           y_stretch = .false.
           y_stretch_factor = 1
           y_compress = .true.
           y_compress_factor = lcy_in / lcy
        else
           y_stretch = .false.
           y_stretch_factor = 1
           y_compress = .false.
           y_compress_factor = 1
        endif
     endif

     if (lcz > lcz_in) then
        z_stretch = .true.
        z_stretch_factor = lcz / lcz_in
        z_compress = .false.
        z_compress_factor = 1
     else
        if (lcz < lcz_in) then
           z_stretch = .false.
           z_stretch_factor = 1
           z_compress = .true.
           z_compress_factor = lcz_in / lcz
        else
           z_stretch = .false.
           z_stretch_factor = 1
           z_compress = .false.
           z_compress_factor = 1
        endif
     endif

     if (rank .eq. 0) then
        read(50,*)
        read(50,*)
        read(50,*)
     endif

     do nu1 = 0, nb-1
        do ia1 = 0, 3

           ind1 = 4*nu1 + ia1

           do nu2 = 0, nb-1
              do ia2 = 0, 3

                 ind2 = 4*nu2 + ia2
              
                 do k3 = 0, min(lcz, lcz_in) - 1
                    do k2 = 0, min(lcy, lcy_in) - 1
                       do k1 = 0, min(lcx, lcx_in) - 1

                          if (rank .eq. 0) then

                             read(50,*) id1, id2, id3, id4, id5, n_pade

                             if (n_pade > 0) then
                                do l = 0, n_pade
                                   read(50,*) p(l), q(l) 
                                enddo
                             endif

                          endif
		  
#ifdef USE_MPI
                          call MPI_Bcast(n_pade, 1, MPI_INTEGER, 0, &
                               MPI_COMM_WORLD, ierr)

                          if (n_pade > 0) then
                             call MPI_Bcast(p, n_pade+1, MPI_COMPLEX, 0, & 
                                  MPI_COMM_WORLD, ierr)
                             call MPI_Bcast(q, n_pade+1, MPI_COMPLEX, 0, & 
                                  MPI_COMM_WORLD, ierr)
                          endif
#endif /* USE_MPI */


                          k = k1*x_stretch_factor + &
                               k2*y_stretch_factor*lcx + &
                               k3*z_stretch_factor*lcx*lcy

                          do l = 0, mp1

                             z = cmplx(0.0d0, abs( epsilon(l) ) )
                  
                             if (epsilon(l) .gt. 0.0d0) then
                                
                                if (n_pade .gt. 0) then

                                   sigma(ind1,ind2,l,k) = &
                                        pade_evaluate(z, n_pade, p, q) * & 
                                        (uu / uu_1)**2

                                else 
                                   sigma(ind1,ind2,l,k) =  cmplx(0.0d0,0.0d0)
                                endif

                             else 

                                if (n_pade .gt. 0) then
                                   sigma(ind1, ind2, l,k) = &
                                        conjg(pade_evaluate(z, n_pade, p, q) ) &
                                        * (uu / uu_1)**2
                                else 
                                   sigma(ind1,ind2, l,k) = cmplx(0.0d0,0.0d0)
                                endif
                      
                             endif     ! (epsilon .gt. 0)
                    
                          enddo       ! l = 0, mp1

                          if (rank .eq. 0) then
                             do i1 = 1, x_compress_factor - 1

                                read(50,*) id1, id2, id3, id4, id5, n_pade
                      
                                if (n_pade > 0) then
                                   do l = 0, n_pade
                                      read(50,*) p(l), q(l) 
                                   enddo
                                endif

                             enddo
                          endif
              
                       enddo         ! k1

                       if (rank .eq. 0) then
                          do i2 = 1, y_compress_factor - 1
                             do i1 = 0, lcx_in - 1

                                read(50,*) id1,id2,id3,id4,id5, n_pade
              
                                if (n_pade > 0) then
                                   do l = 0, n_pade
                                      read(50,*) p(l), q(l) 
                                   enddo
                                endif

                             enddo
                          enddo
          
                       endif

                    enddo           ! k2
      
                    if (rank .eq. 0) then
                       do i3 = 1, z_compress_factor - 1
                          do i2 = 0, lcy_in - 1
                             do i1 = 0, lcx_in - 1

                                read(50,*) id1,id2,id3,id4,id5, n_pade
                      
                                if (n_pade > 0) then
                                   do l = 0, n_pade
                                      read(50,*) p(l), q(l) 
                                   enddo
                                endif

                             enddo
                          enddo
                       enddo

                    endif
                    
                 enddo             ! k3

                 !     Stretch the self-energy.

                 do l = 0, mp1

                    do k3 = 0, lcz-1
                
                       iz(1) = k3 - mod(k3,z_stretch_factor)
                       iz(2) = mod(iz(1) + z_stretch_factor,lcz)
                
                       wz(1) = 1.0d0 - float( k3 - iz(1) ) / & 
                            float(z_stretch_factor)
                
                       wz(2) = 1.0d0 - wz(1)

                       do k2 = 0, lcy-1

                          iy(1) = k2 - mod(k2,y_stretch_factor)
                          iy(2) = mod(iy(1) + y_stretch_factor,lcy)
                
                          wy(1) = 1.0d0 - float( k2 - iy(1) ) / & 
                               float(y_stretch_factor)
                
                          wy(2) = 1.0d0 - wy(1)
              
                          do k1 = 0, lcx-1

                             ix(1) = k1 - mod(k1,x_stretch_factor)
                             ix(2) = mod(ix(1) + x_stretch_factor,lcx)

                             wx(1) = 1.0d0 - float( k1 - ix(1) ) / & 
                                  float(x_stretch_factor)

                             wx(2) = 1.0d0 - wx(1)
                             
                             sigma_temp = cmplx(0.0d0, 0.0d0)

                             do i1 = 1,2
                                do i2 = 1,2
                                   do i3 = 1,2

                                      kp = ix(i1) + iy(i2)*lcx + &
                                           iz(i3)*lcx*lcy
                  
                                      sigma_temp = sigma_temp + &
                                           wx(i1)*wy(i2)*wz(i3) * &
                                           sigma(ind1,ind2,l,kp )
                        
                                   enddo
                                enddo
                             enddo

                             k = k1 + k2*lcx + k3*lcx*lcy
                             sigma(ind1,ind2,l,k) = sigma_temp

                          enddo
                       enddo
                    enddo
          
                 enddo

              enddo
           enddo
        enddo
     enddo

#endif /* SECOND_ORDER */

     if (rank .eq. 0) then
        close (50)
     endif

     return
   end subroutine sigma_input
