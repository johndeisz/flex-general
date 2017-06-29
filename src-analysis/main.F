      program electronic_states

#include "constants.F"

      CHARACTER*128 sigma_input_file, output_file
      DOUBLE PRECISION t, mu, prfld
      DOUBLE PRECISION ed(0:nb-1), h(0:nb-1,1:3)
      DOUBLE COMPLEX tij(0:nb-1,0:nb-1,0:nl-1)
      DOUBLE COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)
      DOUBLE COMPLEX sigma1(0:4*nb-1,0:4*nb-1)

      INTEGER ib, ibp
      DOUBLE COMPLEX ek(0:nb-1,0:nb-1,0:nl-1)
      DOUBLE COMPLEX h_so(0:2*nb-1, 0:2*nb-1)
      DOUBLE COMPLEX h0(0:4*nb-1, 0:4*nb-1, 0:nl-1)

      INTEGER ikx, iky, ikz, ik, ikx_max, is, isp

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
      DOUBLE COMPLEX z, sigmaij
      DOUBLE PRECISION w, a, a00, a11, a22

      INTEGER itype, i_sc_type


      DOUBLE COMPLEX, dimension (0:4*nb-1, 0:4*nb-1) :: hk, 
     $   hk_herm, hkz,
     $   sigma, sigma_p, sigmaz
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

      DOUBLE COMPLEX pade_evaluate, pade_prime
      EXTERNAL pade_evaluate, pade_prime
      
      INTEGER out_type
      DOUBLE COMPLEX gamma, Z_qp, e_qp, eqp_help

      INTEGER ii,jj, iii

      DOUBLE PRECISION up0, up1, up2, down0, down1, down2
      DOUBLE PRECISION wx(-llx/lcx: llx/lcx)
      DOUBLE PRECISION wy(-lly/lcy: lly/lcy)

      COMPLEX, dimension (0:4*nb-1,0:4*nb-1) ::
     $   identity, temp_gc, temp_gl2, temp_gl

      CHARACTER*128 base_out_file
      INTEGER ipiv(4*nb)
      INTEGER ix, ix_max, iy, iy_max
      INTEGER klx, kly, kl, kx, ky

      lwork2 = 16*nb + 32*nb**2
      lrwork2 = 2+20*nb+64*nb**2
      liwork2= 6+20*nb

       do i = 0, 4*nb-1
          do j = 0, 4*nb-1
            if (i .eq. j) then
              identity(i,j) = cmplx(1.0d0, 0.0d0)
            else
              identity(i,j) = cmplx(0.0d0, 0.0d0)
            endif
          enddo
        enddo

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
      write(6,*) "0 - Self-energy"
      write(6,*) "1 - Spectral function"
      write(6,*) "2 - A(k,w) vs. k,w scans"
      write(6,*) "3 - QP energy and width vs. k"
      write(6,*) "4 - DOS"
      write(6,*) "5 - Z vs QP energy"
      write(6,*) "6 - MDC"
      read(5,*) itype

      if (itype .eq. 0) then

        write(6,200) -lcx/2,lcx/2, lcx/2 
        read(5,*) ikx
        write(6,201) -lcy/2,lcy/2,lcy/2 
        read(5,*) iky
        write(6,202) -lcz/2,lcz/2,lcz/2 
        read(5,*) ikz

        if ( (ikx .lt. -lcx/2) .or. (ikx .gt. lcx/2) .or.
     $     (iky .lt. -lcy/2) .or. (iky .gt. lcy/2) .or.
     $     (ikz .lt. -lcz/2) .or. (ikz .gt. lcz/2) ) then
          write(6,*)
          write(6,*) "k-index out of range. Stopping."
          write(6,*)
          stop
        endif

        ik = mod(ikx+lcx,lcx) + mod(iky+lcy,lcy)*lcx + 
     $     mod(ikz+lcz,lcz)*lcx*lcy


       write(6,*) 'Outgoing band index (0,1,2)'
       read(5,*) ib
       write(6,*) 'Outoing spin index (0,1)'
       read(5,*) is
       write(6,*) 'Incoming band index (0,1,2)'
       read(5,*) ibp
       write(6,*) 'Incoming spin index (0,1)'
       read(5,*) isp
       i = 4*ib+is
       j = 4*ibp + isp


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

        write(6,*) 'Name of base output file for sigma'
        read(5,*) base_out_file
        open(unit=91, file='real_'//base_out_file,
     $     status="unknown")
        open(unit=92, file='imag_'//base_out_file,
     $     status="unknown")

        do l = 1, n_w

          w = w_min + dw*dfloat(l-1)
          z = dcmplx(w, eta)

          if ( (sigma_flag .eq. 1) .and.
     $      (n_pade(i,j,ik) .gt. 0)) then
               sigmaij = pade_evaluate(z,n_pade(i,j,ik),
     $                  p(i,j,ik,:), q(i,j,ik,:))
          else
             sigmaij = dcmplx(0.0d0,0.0d0)
          endif
          write(91,*) w, dreal(sigmaij)
          write(92,*) w, dimag(sigmaij)

        enddo

        close(unit=91)
        close(unit=92)

      endif

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

        write(6,*) "Base output file name "
        read(5,*) base_out_file
        open(unit=40, file='all_'//base_out_file, status="unknown")

        open(unit=91, file='band0_'//base_out_file,
     $     status="unknown")
        open(unit=92, file='band1_'//base_out_file,
     $     status="unknown")
        open(unit=93, file='band2_'//base_out_file,
     $     status="unknown")

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

            call dyson2(h0k, ik, p, q, n_pade, sigma_flag, z, a,
     $         a00, a11, a22)

            write(40,*) pi*float(ikx)/float(llx/2), w, a
            write(91,*) pi*float(ikx)/float(llx/2), w, a00
            write(92,*) pi*float(ikx)/float(llx/2), w, a11
            write(93,*) pi*float(ikx)/float(llx/2), w, a22

          enddo

            write(40,*) 
            write(91,*)
            write(92,*) 
            write(93,*)


        enddo

        close(unit=40)
        close(unit=91)
        close(unit=92)
        close(unit=93)

      endif

      if (itype .eq. 6) then

        write(6,*) '1 - Gamma-M scan'
        write(6,*) '2 - Gamma-X scan'
        write(6,*) '3 - M-x scan'
        read(5,*) i_sc_type

        write(6,*)
        write(6,*) "Input energy value "
        read(5,*) w
        write(6,*) "Energy broadening "
        read(5,*) eta

        write(6,*) "Base output file name "
        read(5,*) base_out_file
        open(unit=40, file='all_'//base_out_file, status="unknown")

        open(unit=91, file='band0_'//base_out_file,
     $     status="unknown")
        open(unit=92, file='band1_'//base_out_file,
     $     status="unknown")
        open(unit=93, file='band2_'//base_out_file,
     $     status="unknown")

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

            z = dcmplx(w, eta)

            call dyson2(h0k, ik, p, q, n_pade, sigma_flag, z, a,
     $         a00, a11, a22)

            write(40,*) pi*float(ikx)/float(llx/2), a
            write(91,*) pi*float(ikx)/float(llx/2), a00
            write(92,*) pi*float(ikx)/float(llx/2), a11
            write(93,*) pi*float(ikx)/float(llx/2), a22

        enddo

        close(unit=40)
        close(unit=91)
        close(unit=92)
        close(unit=93)

      endif

      if (itype .eq. 3) then

        write(6,*) '1 - Gamma-Mx scan'
        write(6,*) '2 - Gamma-X scan'
        write(6,*) '3 - M-X scan'
        write(6,*) '4 - ky = 2*kx scan'
        write(6,*) '5 - Fixed ky scan'
        write(6,*) '6 - Gamma-My scan'
        read(5,*) i_sc_type
        if (i_sc_type .eq. 5) then
          write(6,201) -lly/2, lly/2, lly/2
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

        if (i_sc_type .eq. 4) then
          ikx_max = llx/4
        else
          ikx_max = llx/2
        endif

        do ikx = 0, ikx_max

          if ((2*mod(ikx, x_stretch) .ne. x_stretch) .or.
     $       lcx .eq. 1) then

            if (i_sc_type .eq. 6) then
              icx = 0
            else
              icx = nint(float(ikx)/float(x_stretch))
              icx = mod(icx,lcx)
            endif
           
            if (i_sc_type .eq. 1) then
              iky = 0
            endif
            if (i_sc_type .eq. 2) then
              iky = ikx
            endif
            if (i_sc_type .eq. 3) then
              iky = lly/2
            endif
            if (i_sc_type .eq. 4) then
              iky = 2*ikx
            endif
            if (i_sc_type .eq. 6) then
              iky = ikx
            endif

            if ( (2*mod(iky, y_stretch) .ne. y_stretch) .or.
     $         lcy .eq. 1) then

              icy = nint(float(iky)/float(y_stretch))
              icy = mod(icy,lcy)

              if (i_sc_type .eq. 6) then
                 ik = iky*llx
              else
                 ik = ikx + iky*llx 
              endif

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

      if (itype .eq. 4) then

        write(6,*)

        write(6,*) "Base output file name "
        read(5,*) base_out_file
        open(unit=90, file='dos_'//base_out_file, status="unknown")

        open(unit=80, file='dos_up_'//base_out_file,
     $     status="unknown")
        open(unit=81, file='dos_down_'//base_out_file,
     $     status="unknown")

        open(unit=70, file='dos_band0_'//base_out_file,
     $     status="unknown")
        open(unit=71, file='dos_band1_'//base_out_file,
     $     status="unknown")
        open(unit=72, file='dos_band2_'//base_out_file,
     $     status="unknown")

        open(unit=60, file='dos_band0_up_'//base_out_file,
     $     status="unknown")
        open(unit=61, file='dos_band0_down_'//base_out_file,
     $     status="unknown")
        open(unit=62, file='dos_band1_up_'//base_out_file,
     $     status="unknown")
        open(unit=63, file='dos_band1_down_'//base_out_file,
     $     status="unknown")
        open(unit=64, file='dos_band2_up_'//base_out_file,
     $     status="unknown")
        open(unit=65, file='dos_band2_down_'//base_out_file,
     $     status="unknown")

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

        x_stretch = llx / lcx
        y_stretch = lly / lcy

        if (x_stretch .gt. 1) then
          ix_max = x_stretch/2
          do ix = 0, ix_max - 1
            wx(ix) = 1.0d0
            wx(-ix) = 1.0d0
          enddo
          wx(ix_max) = 0.5d0
          wx(-ix_max) = 0.5d0
        else
          ix_max = 0
          wx(0) = 1.0d0
        endif

        if (y_stretch .gt. 1) then
          iy_max = y_stretch/2
          do iy = 0, iy_max - 1
            wy(iy) = 1.0d0
            wy(-iy) = 1.0d0
          enddo
          wy(iy_max) = 0.5d0
          wy(-iy_max) = 0.50d0
        else
          iy_max = 0
          wy(0) = 1.0d0
        endif

 


        do l = 1, n_w
            
          w = w_min + dw*dfloat(l-1)
          z = dcmplx(w, eta)
  
          temp_gc = cmplx(0.0d0, 0.0d0)

          do kx = 0, lcx-1
            do ky = 0, lcy-1

              k = kx + ky*lcx 

              do i = 0, 4*nb-1
                do j = 0, 4*nb-1

                  if ( (sigma_flag .eq. 1) .and.
     $               (n_pade(i,j,k) .gt. 0)) then

                    sigma(i,j) = pade_evaluate(z,n_pade(i,j,k),
     $                 p(i,j,k,:), q(i,j,k,:))
                    
                  else
                    sigma(i,j) = dcmplx(0.0d0,0.0d0)
                  endif

                enddo
              enddo
                
c     For each cluster k-point, sum over all points in the
c     coarse-graining cell.

              do ix = -ix_max, ix_max
                do iy = -iy_max, iy_max

                  klx = mod(kx*x_stretch + ix + llx, llx)
                  kly = mod(ky*y_stretch + iy + lly, lly)

                  kl = klx + kly*llx 

                  temp_gl2 = z*identity - h0(:,:,kl) - sigma
                  temp_gl = identity

                  call cgesv(4*nb, 4*nb, temp_gl2, 4*nb,
     $               ipiv, temp_gl, 4*nb, info)

                  if (info .ne. 0) then
                    write(6,*) 'info not equal to zero'
                  endif

                  temp_gc = temp_gc +
     $               wx(ix) * wy(iy) * temp_gl

                enddo
              enddo

            enddo
          enddo

          up0 = -aimag(temp_gc(0,0))/pi/float(llx*lly)
          down0 = -aimag(temp_gc(1,1))/pi/float(llx*lly)
          up1 = -aimag(temp_gc(4,4))/pi/float(llx*lly)
          down1 = -aimag(temp_gc(5,5))/pi/float(llx*lly)
          up2 = -aimag(temp_gc(8,8))/pi/float(llx*lly)
          down2= -aimag(temp_gc(9,9))/pi/float(llx*lly)

          write(60,*) w, up0
          write(61,*) w, down0
          write(62,*) w, up1
          write(63,*) w, down1
          write(64,*) w, up2
          write(65,*) w, down2
          
          write(70,*) w, 0.5d0*(up0+down0)
          write(71,*) w, 0.5d0*(up1+down1)
          write(72,*) w, 0.5d0*(up2+down2)

          write(80,*) w, (up0+up1+up2)/3.0d0
          write(81,*) w, (down0+down1+down2)/3.0d0

          write(90,*) w, (up0+up1+up2+down0+down1+down2)/6.0d0

        enddo
      endif

      if (itype .eq. 5) then

        write(6,*) '1 - Gamma-M scan'
        write(6,*) '2 - Gamma-X scan'
        write(6,*) '3 - M-X scan'
        write(6,*) '4 - ky = 2*kx scan'
        write(6,*) '5 - Fixed ky scan'
        read(5,*) i_sc_type
        if (i_sc_type .eq. 5) then
          write(6,201) -lly/2, lly/2, lly/2
          read(5,*) iky
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

        write(6,*) 'Name of output file'
        read(5,*) output_file
        open(unit=40, file=output_file, status="unknown")

        x_stretch = llx / lcx
        y_stretch = lly / lcy

        if (i_sc_type .eq. 4) then
          ikx_max = llx/4
        else
          ikx_max = llx/2
        endif

        do ikx = 0, ikx_max

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
            if (i_sc_type .eq. 4) then
              iky = 2*ikx
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

                      z = dcmplx(w, eta)

c$$$		        do iii = 1,100

                      do ii = 0, 4*nb-1
                        do jj = 0, 4*nb-1
			   

                          if ( (sigma_flag .eq. 1) .and.
     $                       (n_pade(ii,jj,ic) .gt. 0)) then

                            sigmaz(ii,jj) = 
     $                         pade_evaluate(z,n_pade(ii,jj,ic),
     $                         p(ii,jj,ic,:), q(ii,jj,ic,:))

                          else    
                            sigmaz(ii,jj) = dcmplx(0.0d0,0.0d0)
                          endif  

                        enddo
                      enddo

                      e_qp = dcmplx(0.0d0, 0.0d0)

                      hkz = h0k + sigmaz

                      do ii = 0, 4*nb-1
                        do jj = 0, 4*nb-1
                          e_qp = e_qp + conjg(hk_herm(ii,j))*
     $                       hkz(ii,jj)* hk_herm(jj,j)
                        enddo
                      enddo

c$$$                      eqp_help = dconjg(e_qp)

c$$$                      write(6,*) iii, cdabs(e_qp - z)
c$$$                      z = 0.9*z + 0.1*e_qp
c$$$
c$$$                    enddo
		

                      do ii = 0, 4*nb-1
                        do jj = 0, 4*nb-1

                          if ( (sigma_flag .eq. 1) .and.
     $                       (n_pade(ii,jj,ic) .gt. 0)) then

                            sigma_p(ii,jj) = 
     $                         pade_prime(z,n_pade(ii,jj,ic),
     $                         p(ii,jj,ic,:), q(ii,jj,ic,:))

                          else    
                            sigma_p(ii,jj) = dcmplx(0.0d0,0.0d0)
                          endif   

                        enddo   
                      enddo   
c     
c     temp_gl2 = identity - sigma
                      temp_gl2 = identity - sigma_p
                      temp_gl = identity

                      call cgesv(4*nb, 4*nb, temp_gl2, 4*nb,
     $                   ipiv, temp_gl, 4*nb, info)

                      Z_qp = dcmplx(0.0d0, 0.0d0)

                      do ii = 0, 4*nb-1
                        do jj = 0, 4*nb-1
                          Z_qp = Z_qp + conjg(hk_herm(ii,j))*
     $                       temp_gl(ii,jj)*hk_herm(jj,j)
                        enddo
                      enddo
                    
c     write(40,*) eigs(j), real(Z_qp)
                      write(40,*) 
     $                   2.0d0*float(ikx)/float(llx), 
     $                   dreal(e_qp), -dimag(e_qp),
     $                   dreal(Z_qp), dimag(Z_qp) 

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

 200  format("kx index, [",i4,",",i4,"]:: 0 -> kx = 0, ", 
     $   i3, "-> kx = pi/a")
 201  format("ky index, [",i4,",",i4,"]:: 0 -> ky = 0, ", 
     $   i3, "-> ky = pi/b")
 202  format("kz index, [",i4,",",i4,"]:: 0 -> kz = 0, ", 
     $   i3, "-> kz = pi/c")

 300  format("low eigenval of h0k+sig1 = ", e16.8)
 301  format("high eigenval of h0k+sig1 = ", e16.8)

      stop
      end
