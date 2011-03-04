      subroutine readin(sigma_input_file, t, mu, prfld, 
     $   ed, h, tij, psi, sigma1, h_so)

#include "constants.F"

      CHARACTER*128 sigma_input_file
      INTEGER lcxi, lcyi, lczi, llxi, llyi, llzi
      INTEGER nbi
      INTEGER id1, id2, id3, id4, id5
      DOUBLE PRECISION t, mu, prfld
      DOUBLE PRECISION dum1, dum2, dum3
      DOUBLE PRECISION ed(0:nb-1), h(0:nb-1,1:3)
      DOUBLE COMPLEX tij(0:nb-1,0:nb-1,0:nl-1)
      DOUBLE COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)
      DOUBLE COMPLEX sigma1(0:4*nb-1,0:4*nb-1)
      DOUBLE COMPLEX h_so(0:2*nb-1, 0:2*nb-1)
      INTEGER ib, ibp
      INTEGER max_x, max_y, max_z, ix, iy, iz, k
      INTEGER id, idp, pr_flag, so_flag
      integer nu1, is1, nu2, is2, ia1, ia2
      
      open(unit=50, file=sigma_input_file, status='old')
      read(50,*)
      read(50,*)
      read(50,*)
      read(50,*)
      read(50,*)
      read(50,*) nbi
      read(50,*)
      read(50,*) lcxi, lcyi, lczi, id1, llxi, llyi,llzi

      if ( (lcxi .ne. lcx) .or. (lcyi .ne. lcy) .or. 
     $   (lczi .ne. lcz) .or. (llxi .ne. llx) .or. 
     $   (llyi .ne. lly) .or. (llzi .ne. llz) .or. 
     $   (nbi .ne. nb) ) then
        write(6,*) "Compile time flags do not match input ",
     $     "dimensions. Stopping."
        stop
      endif

      read(50,*)  ! string
      read(50,*)   t, mu, dum1, dum2, dum3 
      read(50,*)  ! string
      read(50,*)  prfld
      read(50,*)  ! string
      do ib = 0, nb - 1
        read(50,*) id1, ed(ib) 
      enddo
      read(50,*)  ! string
      do ib = 0, nb - 1
        read(50,*)  id1, h(ib,:)
      enddo
      read(50,*) ! string
      read(50,*) ! uu_1, up_1, uj_1
      read(50,*) ! string

c-----Hopping matrix elements -------------------

      tij = dcmplx(0.0d0, 0.0d0)

      if (llx .gt. 2) then
        max_x = 2
      else
        max_x = llx - 1
      endif

      if (lly .gt. 2) then
        max_y = 2
      else
        max_y = lly - 1
      endif

      if (llz .gt. 2) then
        max_z = 2
      else
        max_z = llz - 1
      endif
          
      do ix = -max_x, max_x
        do iy = -max_y, max_y
          do iz = -max_z, max_z

            read(50,*)
              
            k = mod(ix+llx,llx) + mod(iy+lly,lly)*llx +
     $           mod(iz+llz,llz)*llx*lly

            do ib = 0, nb-1
              do ibp = 0, nb-1
                read(50,*) id, idp, tij(ib,ibp,k)
              enddo
            enddo

          enddo
        enddo
      enddo

c     Pair wave function

      read(50,*)
      read(50,*)
      read(50,*)

      psi = dcmplx(0.0d0, 0.0d0)
      read(50,*) pr_flag
      if (pr_flag .eq. 1) then
        read(50,*)
        do nu1 = 0, nb-1
          do is1 = 0, 1
            do nu2 = 0,nb-1
              do is2 = 0, 1

                do k = 0, nl-1

                  read(50,*) id1, id2, id3, id4, id5,
     $               psi(2*nu1+is1,2*nu2+is2,k)
                  
                enddo
              enddo
            enddo

          enddo
        enddo
      endif

c Spin-orbit interaction: default to zero
      h_so = dcmplx(0.0d0, 0.d0)
      read(50,*)
      read(50,*)
      read(50,*)
      read(50,*) so_flag
      if (so_flag .eq. 1) then
        read(50,*)
        do nu1 = 0, nb-1
          do is1 = 0, 1
            do nu2 = 0,nb-1
              do is2 = 0, 1

                read(50,*) id1, id2, id3, id4, 
     $             h_so(2*nu1+is1,2*nu2+is2)
                  
              enddo
            enddo
          enddo
        enddo

      endif

c Sigma_1

      read(50,*)
      read(50,*)
      read(50,*)

      do nu1 = 0, nb-1
        do ia1 = 0, 3
          do nu2 = 0, nb-1
            do ia2 = 0, 3
              read(50,*) id1, id2, id3, id4, 
     $           sigma1(nu1*4+ia1, nu2*4+ia2)            
            enddo
          enddo
        enddo
      enddo

      read(50,*)

      return
      end
