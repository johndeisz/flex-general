      program spinorbit

      implicit none

      integer ml, ms

c ms = 0 -> spin up
c ms = 1 -> spin down

      complex double hspin(0:9,0:9)
      complex double psi(-2:2, 0:3)
      complex double psi_spin(0:9,0:5)

      hspin = dcmplx(0.0d0,0.0d0)

c Define spatial wave functions
      psi(-2,0) = 0.0d0
      psi(-1,0) = -1.0d0 / dsqrt(2.0d0) 
      psi(0,0) = 0.0d0
      psi(1,0) = -1.0d0 / dsqrt(2.0d0)
      psi(2,0) = 0.0d0

      psi(-2,1) = 0.0d0 
      psi(-1,1) = dcmplx(0.0d0,-1.0d0 / dsqrt(2.0d0))
      psi(0,1) = 0.0d0
      psi(1,1) = dcmplx(0.0d0, 1.0d0 / dsqrt(2.0d0))
      psi(2,1) = 0.0d0

      psi(-2,2) = dcmplx(0.0d0, 1.0d0 / dsqrt(2.0d0))
      psi(-1,2) = 0.0d0
      psi(0,2) = 0.d0
      psi(1,2) = 0.d0
      psi(2,2) = dcmpx( 0.0d0, -1.0d0 / dsqrt(2.0d0))

      psi_spin = cmplx(0.0d0, 0.0d0)

      do nu = 0, 2
        do is = 0, 1
          nus = 2*nu + is

          do ml = -2, 2
            ind = 2*(ml+2) + is

            psi_spin(ind, nus) = psi(ml,nu)
          enddo

        enddo
      enddo

      do ml = -2, 2
        do ms = 0, 1

          ind1 = 2*(ml+2)+ms
          
c Diagonal element
          ind2 = ind1
          hspin(ind2,ind1) = dfloat(ml)*(0.5d0-dfloat(ms))

c Orbital angular momentum raised, spin angular momentum lowered
          if ( ml .ne. 2) then
            if (ms .eq. 0) then

              ind2 = 2*(ml+3) + 1 
              hspin(ind2,ind1) = dsqrt(6.0d0 - dfloat(ml*(ml+1)))

            endif
          endif 

c Orbital angular momentum lowered, spin angular momentum raised
          if (ml .ne. -2) then
            if (ms .eq. 1) then

              ind2 = 2*(ml+1) + 0
              hspin(ind2,ind1) =dsqrt(6.0d0 - dfloat(ml*(ml-1)))

            endif
          endif
               
        enddo
      enddo

       


           
