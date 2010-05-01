      program electronic_states

#include "constants.F"

      CHARACTER*128 sigma_input_file
      DOUBLE PRECISION t, mu, prfld
      DOUBLE PRECISION ed(0:nb-1), h(0:nb-1,1:3)
      DOUBLE COMPLEX tij(0:nb-1,0:nb-1,0:nl-1)
      DOUBLE COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)
      DOUBLE COMPLEX sigma1(0:4*nb-1,0:4*nb-1)

      INTEGER ib
      DOUBLE COMPLEX ek(0:nb-1,0:nb-1,0:nl-1)
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

      write(6,*) "Self-energy file to be used"
      read(5,*) sigma_input_file

      call readin(sigma_input_file, t, mu, prfld, ed, h, 
     $   tij, psi, sigma1)

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
      call h0_generate(h0, ek, psi, h, prfld, sigma1)

      write(6,200) -llx/2,llx/2, llx/2 
      read(5,*) ikx
      write(6,201) -lly/2,lly/2,lly/2 
      read(5,*) iky
      write(6,201) -llz/2,llz/2,llz/2 
      read(5,*) ikz

      if ( (ikx .lt. -llx/2) .or. (ikx .gt. llx/2) .or.
     $   (iky .lt. -lly/2) .or. (iky .gt. lly/2) .or.
     $   (ikz .lt. -llz/2) .or. (ikz .gt. llz/2) ) then
        write(6,*)
        write(6,*) "k-index out of range. Stopping."
        write(6,*)
        stop
      endif

      ik = mod(ikx+llx,llx) + mod(iky+lly,lly)*llx + 
     $   mod(ikz+llz,llz)*llx*lly

      h0k = h0(:,:,ik)

      tmp_h0k = h0k

      lwork = 16*nb*nb
      call ZHEEVD('N', 'U',  4*nb,  tmp_h0k,  4*nb,  eigs_1st,  
     $   WORK,  LWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )

 
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
     $   "should be larger than energy increment."
      read(5,*) eta
      if ( eta .lt. dw ) then
        write(6,*) "Resolution paramater too small. Stopping."
        stop
      endif

      call dyson(h0k,ik, w_min,dw,n_w,eta)

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
