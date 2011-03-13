      program electronic_states

#include "constants.F"

      CHARACTER*128 sigma_input_file, output_file
      DOUBLE PRECISION t, mu, prfld
      DOUBLE PRECISION ed(0:nb-1), h(0:nb-1,1:3)
      DOUBLE COMPLEX tij(0:nb-1,0:nb-1,0:nl-1)
      DOUBLE COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)
      DOUBLE COMPLEX sigma1(0:4*nb-1,0:4*nb-1)

      INTEGER ib
      DOUBLE COMPLEX ek(0:nb-1,0:nb-1,0:nl-1)
      DOUBLE COMPLEX h_so(0:2*nb-1, 0:2*nb-1)
      DOUBLE COMPLEX h0(0:4*nb-1, 0:4*nb-1, 0:nl-1)

      INTEGER ikx, iky, ikz, ik

      DOUBLE COMPLEX, dimension (0:4*nb-1, 0:4*nb-1) :: h0k, 
     $   tmp_h0k

      DOUBLE PRECISION eigs_1st(1:4*nb)
      INTEGER lwork, lrwork, liwork
      DOUBLE COMPLEX work(1:16*nb*nb)
      DOUBLE PRECISION rwork(1:16*nb*nb)
      INTEGER iwork(1:16*nb*nb)
      INTEGER info

      DOUBLE PRECISION w_min, w_max, dw, eta
      INTEGER n_w

      DOUBLE COMPLEX, dimension (0:4*nb-1,0:4*nb-1, 0:nc-1, 
     $   0:n_pade_max) :: p, q

      INTEGER sigma_flag
      INTEGER n_pade(0:4*nb-1,0:4*nb-1,0:nc-1)
      INTEGER k, nu1, ia1, nu2, ia2, na1, na2
      INTEGER id1, id2, id3, id4, id5

      INTEGER l
      DOUBLE COMPLEX z
      DOUBLE PRECISION w, a

      INTEGER itype, i_sc_type


      DOUBLE COMPLEX, dimension (0:4*nb-1, 0:4*nb-1) :: hk, hk_herm,
     $   sigma
      DOUBLE PRECISION eigs(0:4*nb-1)
      DOUBLE COMPLEX work2(16*nb+32*nb**2)
      INTEGER lwork2
      DOUBLE PRECISION rwork2(2+20*nb+64*nb**2)
      INTEGER lrwork2
      INTEGER iwork2(6+20*nb)
      INTEGER liwork2, info2

      INTEGER i,j

      DOUBLE COMPLEX a_up_spin
      INTEGER eigs_smaller_new, eigs_smaller_old

      INTEGER ic, icx, icy
      INTEGER x_stretch, y_stretch

      DOUBLE COMPLEX pade_evaluate
      EXTERNAL pade_evaluate
      
      INTEGER out_type
      DOUBLE COMPLEX gamma

      INTEGER ii,jj

      lwork2 = 16*nb + 32*nb**2
      lrwork2 = 2+20*nb+64*nb**2
      liwork2= 6+20*nb

      write(6,*) "Self-energy file to be used"
      read(5,*) sigma_input_file

      call readin(sigma_input_file, t, mu, prfld, ed, h, 
     $   tij, psi, sigma1, h_so)

      write(6,*) "t = ", t
      write(6,*) "mu = ", mu
      write(6,*) "prfld = ", prfld
      do ib = 0, nb-1
        write(6,*) "ib, ed ", ib," ", ed(ib)
      enddo
      do ib = 0, nb-1
        write(6,*) "ib, h = ", ib, " ", h(ib,:)
      enddo

      t = kb*t

c     Bare bandstructure and one-electron matrix
      call bare_dispersion(tij, ed, ek, mu)
      call h0_generate(h0, ek, psi, h, prfld, sigma1, h_so)

      write(6,*) "Dynamical self-energy? Yes-1, No-0."
      read(5,*) sigma_flag

      p = dcmplx(0.0d0, 0.0d0)
      q = dcmplx(0.0d0, 0.0d0)
      
      if (sigma_flag .eq. 1) then
        read(50,*) 
        read(50,*)
c     Read pade coefficients
        do nu1 = 0, nb-1
          do ia1 = 0, 3
            na1 = 4*nu1 + ia1
            do nu2 = 0, nb-1
              do ia2 = 0, 3
                na2 = 4*nu2 + ia2
                do k = 0, nc-1
                  read(50,*) id1, id2, id3, id4, id5, 
     $               n_pade(na1,na2,k)
                  if (n_pade(na1,na2,k) .gt. 0) then
                    do l = 0, n_pade(na1,na2,k)
                      read(50,*) p(na1, na2, k, l),
     $                   q(na1,na2,k,l)
                    enddo
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      endif

      close(unit=50)

      write(6,*) "Calculation type:"
      write(6,*) "1 - Spectral function"
      write(6,*) "2 - A(k,w) vs. k,w scans"
      write(6,*) "3 - QP energy and width vs. k"
      read(5,*) itype

      if (itype .eq. 1) then

        write(6,200) -llx/2,llx/2, llx/2 
        read(5,*) ikx
        write(6,201) -lly/2,lly/2,lly/2 
        read(5,*) iky
        write(6,202) -llz/2,llz/2,llz/2 
        read(5,*) ikz

        if ( (ikx .lt. -llx/2) .or. (ikx .gt. llx/2) .or.
     $     (iky .lt. -lly/2) .or. (iky .gt. lly/2) .or.
     $     (ikz .lt. -llz/2) .or. (ikz .gt. llz/2) ) then
          write(6,*)
          write(6,*) "k-index out of range. Stopping."
          write(6,*)
          stop
        endif

        ik = mod(ikx+llx,llx) + mod(iky+lly,lly)*llx + 
     $     mod(ikz+llz,llz)*llx*lly

        h0k = h0(:,:,ik)

        tmp_h0k = h0k

        lwork = 16*nb*nb
        lrwork = 16*nb*nb
        liwork = 16*nb*nb
        call ZHEEVD('N', 'U',  4*nb,  tmp_h0k,  4*nb,  eigs_1st,  
     $     work,  lwork,  rwork, lrwork, iwork, liwork, info )
 
        if (info .ne. 0) then
          write(6,*) "zheevd failed. stopping."
        else
          write(6,300) eigs_1st(1)
          write(6,301) eigs_1st(4*nb)
        endif

        write(6,*)
        write(6,*) "Input minimum energy value "
        read(5,*) w_min
        write(6,*) "Input maximum energy value "
        read(5,*) w_max
        write(6,*) "Num energy points "
        read(5,*) n_w 
        dw = (w_max - w_min) / dfloat(n_w - 1)
        write(6,*) "energy increment = ", dw

        write(6,*) "Resolution parameter - ",
     $     "should be larger than energy increment."
        read(5,*) eta
        if ( eta .lt. dw ) then
          write(6,*) "Resolution paramater too small. Stopping."
          stop
        endif

        write(6,*) 'Name of output file for a(k,w)'
        read(5,*) output_file
        open(unit=40, file=output_file, status="unknown")

        do l = 1, n_w

          w = w_min + dw*dfloat(l-1)
          z = dcmplx(w, eta)
          call dyson(h0k, ik, p, q, n_pade, sigma_flag, z, a)

          write(40,*) w, a

        enddo

        close(unit=40)

      endif

      if (itype .eq. 2) then

        write(6,*) '1 - Gamma-M scan'
        write(6,*) '2 - Gamma-X scan'
        write(6,*) '3 - M-x scan'
        read(5,*) i_sc_type

        write(6,*)
        write(6,*) "Input minimum energy value "
        read(5,*) w_min
        write(6,*) "Input maximum energy value "
        read(5,*) w_max
        write(6,*) "Num energy points "
        read(5,*) n_w 
        dw = (w_max - w_min) / dfloat(n_w - 1)
        write(6,*) "energy increment = ", dw

        write(6,*) "Resolution parameter - ",
     $     "should be larger than energy increment."
        read(5,*) eta

        write(6,*) 'Name of output file'
        read(5,*) output_file
        open(unit=40, file=output_file, status="unknown")

        do ikx = 0, llx/2

          if (i_sc_type .eq. 1) then
            iky = 0
          endif
          if (i_sc_type .eq. 2) then
            iky = ikx
          endif
          if (i_sc_type .eq. 3) then
            iky = llx/2
          endif

          ik = mod(ikx+llx,llx) + mod(iky+lly,lly)*llx + 
     $       mod(ikz+llz,llz)*llx*lly

          h0k = h0(:,:,ik)

          do l = 1, n_w
            
            w = w_min + dw*dfloat(l-1)
            z = dcmplx(w, eta)

            call dyson(h0k, ik, p, q, n_pade, sigma_flag, z, a)

            write(40,*) pi*float(ikx)/float(llx/2), w, a

          enddo

        enddo

        close(unit=40)

      endif

      if (itype .eq. 3) then

        write(6,*) '1 - Gamma-M scan'
        write(6,*) '2 - Gamma-X scan'
        write(6,*) '3 - M-X scan'
        write(6,*) '4 - Fixed ky scan'
        read(5,*) i_sc_type
        if (i_sc_type .eq. 4) then
          write(6,201)
          read(5,*) iky
        endif
        write(6,*)
        write(6,*) '1 - Eqp and gamma vs. k/(pi/a)'
        write(6,*) '2 - Gamma vs Eqp'
        read(5,*) out_type
        write(6,*)
        write(6,*) "Input minimum energy value "
        read(5,*) w_min
        write(6,*) "Input maximum energy value "
        read(5,*) w_max
        write(6,*) "Num energy points "
        read(5,*) n_w 
        dw = (w_max - w_min) / dfloat(n_w - 1)
        write(6,*) "energy increment = ", dw

        write(6,*) "Resolution parameter - ",
     $     "should be larger than energy increment."
        read(5,*) eta

        write(6,*) 'Name of output file'
        read(5,*) output_file
        open(unit=40, file=output_file, status="unknown")

        x_stretch = llx / lcx
        y_stretch = lly / lcy

        do ikx = 0, llx/2

          if ((2*mod(ikx, x_stretch) .ne. x_stretch) .or.
     $       lcx .eq. 1) then

            icx = nint(float(ikx)/float(x_stretch))
            icx = mod(icx,lcx)
           
            if (i_sc_type .eq. 1) then
              iky = 0
            endif
            if (i_sc_type .eq. 2) then
              iky = ikx
            endif
            if (i_sc_type .eq. 3) then
                iky = lly/2
            endif

            if ( (2*mod(iky, y_stretch) .ne. y_stretch) .or.
     $         lcy .eq. 1) then

              icy = nint(float(iky)/float(y_stretch))
              icy = mod(icy,lcy)

              ik = ikx + iky*llx 

              h0k = h0(:,:,ik)

              ic = icx + icy*lcx

              eigs_smaller_old = 0

              do l = 1, n_w
            
                w = w_min + dw*dfloat(l-1)
                z = dcmplx(w, eta)

                do i = 0, 4*nb-1
                  do j = 0, 4*nb-1

                    if ( (sigma_flag .eq. 1) .and.
     $                 (n_pade(i,j,ic) .gt. 0)) then

                      sigma(i,j) = pade_evaluate(z,n_pade(i,j,ic),
     $                   p(i,j,ic,:), q(i,j,ic,:))

                    else
                      sigma(i,j) = dcmplx(0.0d0,0.0d0)
                    endif

                  enddo
                enddo

                hk = h0k + sigma

                do i = 0, 4*nb-1
                  do j= 0, 4*nb-1

                    hk_herm(i,j) = 
     $                 0.5d0*(hk(i,j) + conjg(hk(j,i)))

                  enddo
                enddo

                call zheevd('V','U',4*nb,hk_herm,4*nb,eigs,work2,
     $             lwork2,rwork2,lrwork2,iwork2,liwork2,info2)

                if (info2 .ne. 0) then
                  write(6,*) "zheevd failed. stopping."
                endif

                eigs_smaller_new = 0

                do i = 0, 4*nb-1

                  if (eigs(i) .lt. w) then
                    eigs_smaller_new = eigs_smaller_new + 1
                  endif

                enddo

                if ( (l .gt. 1) .and. (eigs_smaller_new .gt. 
     $             eigs_smaller_old) ) then

                  do j = eigs_smaller_old, eigs_smaller_new - 1

                    a_up_spin = 0.0d0

                    do ib = 0, nb-1

                      a_up_spin = a_up_spin + hk_herm(4*ib,j)*
     $                   conjg(hk_herm(4*ib,j))

                    enddo

                    if (dreal(a_up_spin) .gt. 0.1d0) then

                      gamma = 0

                      do ii = 0, 4*nb-1
                        do jj = 0, 4*nb-1
                          gamma = gamma + conjg(hk_herm(ii,j))*
     $                       sigma(ii,jj)*hk_herm(jj,j)
                        enddo
                      enddo
                      
                      if (out_type .eq. 1) then
                        write(40,*) 2.0d0*float(ikx)/float(llx),
     $                     eigs(j), -dimag(gamma)
                      else
                        write(40,*) eigs(j), -dimag(gamma)
                      endif

                    endif

                  enddo
                endif

                eigs_smaller_old = eigs_smaller_new
                  
              enddo

            endif

          endif

        enddo
      
        close(unit=40)
        
      endif

 200  format("kx index, [",i3,",",i3,"]:: 0 -> kx = 0, ", 
     $   i3, "-> kx = pi/a")
 201  format("ky index, [",i3,",",i3,"]:: 0 -> ky = 0, ", 
     $   i3, "-> ky = pi/b")
 202  format("kz index, [",i3,",",i3,"]:: 0 -> kz = 0, ", 
     $   i3, "-> kz = pi/c")

 300  format("low eigenval of h0k+sig1 = ", e16.8)
 301  format("high eigenval of h0k+sig1 = ", e16.8)

      stop
      end
