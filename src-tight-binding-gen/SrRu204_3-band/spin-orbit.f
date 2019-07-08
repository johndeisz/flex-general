      program spinorbit

      implicit none

c ms = 0 -> spin up
c ms = 1 -> spin down

      double complex hspin(0:9,0:9)
      double complex psi(-2:2, 0:2)
      double complex psi_spin(0:9,0:5)

      integer nus, nu, is, ml, ind, ms, ind1, ind2
      integer nus1, nus2

      double complex sum

      hspin = dcmplx(0.0d0,0.0d0)

c Define spatial wave functions
      psi(-2,0) = dcmplx(0.0d0, 0.0d0)
      psi(-1,0) = 1.0d0 / dsqrt(2.0d0) 
      psi(0,0) = 0.0d0
      psi(1,0) = -1.0d0 / dsqrt(2.0d0)
      psi(2,0) = 0.0d0

      psi(-2,1) = 0.0d0 
      psi(-1,1) = dcmplx(0.0d0,1.0d0 / dsqrt(2.0d0))
      psi(0,1) = 0.0d0
      psi(1,1) = dcmplx(0.0d0, 1.0d0 / dsqrt(2.0d0))
      psi(2,1) = 0.0d0

      psi(-2,2) = dcmplx(0.0d0, 1.0d0 / dsqrt(2.0d0))
      psi(-1,2) = 0.0d0
      psi(0,2) = 0.d0
      psi(1,2) = 0.d0
      psi(2,2) = dcmplx( 0.0d0, -1.0d0 / dsqrt(2.0d0))

      psi_spin = dcmplx(0.0d0, 0.0d0)

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
              hspin(ind2,ind1) = 0.5d0*dsqrt(6.0d0 - dfloat(ml*(ml+1)))

            endif
          endif 

c Orbital angular momentum lowered, spin angular momentum raised
          if (ml .ne. -2) then
            if (ms .eq. 1) then

              ind2 = 2*(ml+1) + 0
              hspin(ind2,ind1) = 0.5d0*dsqrt(6.0d0 - dfloat(ml*(ml-1)))

            endif
          endif
               
        enddo
      enddo
       
c Compute the spin-orbit matrix elements
	
      do nus1=0,5
        do nus2 = 0,5

          sum = 0.0d0

          do ind1 = 0,9
            do ind2 = 0, 9

              sum = sum + dconjg(psi_spin(ind1,nus1))*
     $           hspin(ind1,ind2)*psi_spin(ind2,nus2)
              
            enddo
          enddo
            
          write(6,*) nus1, nus2, sum

        enddo
      enddo

      stop
      end



           
