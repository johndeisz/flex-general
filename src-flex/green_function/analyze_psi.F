#include "../convert.F"

      subroutine analyze_psi(psi, tij)

#include "../constants.F"

      COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)
      COMPLEX tij(0:nb-1,0:nb-1,0:nl-1)

      INTEGER isignv(0:3)
      INTEGER i, j, ix, iy, iz, jx, jy, jz
      COMPLEX total
      REAL pi

      COMPLEX dot
      REAL dot_real, dot_imag
      REAL mag

      pi = 4.0*atan(1.0d0)


c     First fft the wave function to r-space
      isignv(0) = 0
      isignv(1) = 1
      isignv(2) = 1
      isignv(3) = 1

      call psi_transform(psi, isignv)

      psi = psi / float(nl)

c Determine overall wavefunction normalization

      overall = 0.0d0

      do is = 0, 1
        do isp = 0, 1
          do ib = 0, nb-1
            do ibp = 0, nb - 1
              do ir = 0, nl - 1
              
                overall = overall + 
     $             cabs(psi(2*ib+is,2*ibp+isp,ir))**2.0d0
                
              enddo
            enddo
          enddo
        enddo
      enddo

c Determine pairing by bind and spin sector

c Singlet-singlet

      norm = 0.0d0

      do ib = 0, nb-1
        do ibp = 0, nb -1
          do ir = 0, nl - 1
  
            overlap = dcmplx(0.0d0, 0.0d0)
            
            do is = 0, 1
              
              if (is .eq. 0) then
                phase = 1.0d0
              else
                phase = -1.0d0
              endif

              do isp = 0, 1
                
                if (isp .eq. 0) then
                  phasep = 1.0d0
                else
                  phasep = -1.0d0
                endif
            
                overlap = overlap + 
     $             (1.0d0/2.0d0)*phase*phasep*
     $             psi(2*ib+is,2*ibp+isp,ir)

              enddo
            enddo

            norm = norm + cabs(overlap)**2.0 
            
          enddo
        enddo
      enddo
        
      write(6,*) '<S=0 Sz=0|psi|S=0 Sz=0> = ', norm / overall

c Triplet-triplet

      norm = 0.0d0

      do ib = 0, nb-1
        do ibp = 0, nb -1
          do ir = 0, nl - 1
            
            norm = norm + 
     $         (cabs(psi(2*ib+0,2*ibp+0,ir)))**2.0d0
            
          enddo
        enddo
      enddo

      write(6,*) '<S=1 Sz=1|psi|S=1 Sz=1> = ', norm / overall

      norm = 0.0d0

      do is = 0, 1

        if (is .eq. 0) then
          phase = 1.0d0
        else
          phase = 1.0d0
        endif

        do isp = 0, 1
          
          if (isp .eq. 0) then
            phasep = 1.0d0
          else
            phasep = 1.0d0
          endif

          do ib = 0, nb-1
            do ibp = 0, nb -1
              do ir = 0, nl - 1
              
                norm = norm + 
     $             (1.0d0/4.0d0)*phase*phasep*
     $             (cabs(psi(2*ib+is,2*ibp+isp,ir)))**2.0d0
                
              enddo
            enddo
          enddo

        enddo
      enddo

      write(6,*) '<S=1 Sz=0|psi|S=1 Sz=0> = ', norm / overall

      norm = 0.0d0

      do ib = 0, nb-1
        do ibp = 0, nb -1
          do ir = 0, nl - 1
            
            norm = norm + 
     $         (cabs(psi(2*ib+1,2*ibp+1,ir)))**2.0d0
            
          enddo
        enddo
      enddo

      write(6,*) '<S=1 Sz=-1|psi|S=1 Sz=-1> = ', norm / overall
 
c Singlet-triplet

      norm = 0.0d0

      do is = 0, 1

        if (is .eq. 0) then
          phase = 1.0d0
        else
          phase = -1.0d0
        endif

        do ib = 0, nb-1
          do ibp = 0, nb -1
            do ir = 0, nl - 1
              
              norm = norm + 
     $           (1.0d0/2.0d0)*phase*
     $           (cabs(psi(2*ib+is,2*ibp+0,ir)))**2.0d0
                
            enddo
          enddo
        enddo

      enddo

      write(6,*) '<S=0 Sz=0|psi|S=1 Sz=1> = ', norm / overall    

      norm = 0.0d0

      do is = 0, 1

        if (is .eq. 0) then
          phase = 1.0d0
        else
          phase = -1.0d0
        endif

        do isp = 0, 1
          
          if (isp .eq. 0) then
            phasep = 1.0d0
          else
            phasep = 1.0d0
          endif

          do ib = 0, nb-1
            do ibp = 0, nb -1
              do ir = 0, nl - 1
              
                norm = norm + 
     $             (1.0d0/4.0d0)*phase*phasep*
     $             (cabs(psi(2*ib+is,2*ibp+isp,ir)))**2.0d0
                
              enddo
            enddo
          enddo

        enddo
      enddo

      write(6,*) '<S=0 Sz=0|psi|S=1 Sz=0> = ', norm / overall

      norm = 0.0d0

      do is = 0, 1

        if (is .eq. 0) then
          phase = 1.0d0
        else
          phase = -1.0d0
        endif

        do ib = 0, nb-1
          do ibp = 0, nb -1
            do ir = 0, nl - 1
              
              norm = norm + 
     $           (1.0d0/2.0d0)*phase*
     $           (cabs(psi(2*ib+is,2*ibp+1,ir)))**2.0d0
                
            enddo
          enddo
        enddo

      enddo

      write(6,*) '<S=0 Sz=0|psi|S=1 Sz=-1> = ', norm / overall    

c Triplet-singlet

      norm = 0.0d0

      do isp = 0, 1

        if (isp.eq. 0) then
          phasep = 1.0d0
        else
          phasep = -1.0d0
        endif

        do ib = 0, nb-1
          do ibp = 0, nb -1
            do ir = 0, nl - 1
              
              norm = norm + 
     $           (1.0d0/2.0d0)*phasep*
     $           (cabs(psi(2*ib+0,2*ibp+isp,ir)))**2.0d0
                
            enddo
          enddo
        enddo

      enddo

      write(6,*) '<S=1 Sz=1|psi|S=0 Sz=0> = ', norm / overall    

      norm = 0.0d0

      do is = 0, 1

        if (is .eq. 0) then
          phase = 1.0d0
        else
          phase = 1.0d0
        endif

        do isp = 0, 1
          
          if (isp .eq. 0) then
            phasep = 1.0d0
          else
            phasep = -1.0d0
          endif

          do ib = 0, nb-1
            do ibp = 0, nb -1
              do ir = 0, nl - 1
              
                norm = norm + 
     $             (1.0d0/4.0d0)*phase*phasep*
     $             (cabs(psi(2*ib+is,2*ibp+isp,ir)))**2.0d0
                
              enddo
            enddo
          enddo

        enddo
      enddo

      write(6,*) '<S=1 Sz=0|psi|S=0 Sz=0> = ', norm / overall

      norm = 0.0d0

      do isp = 0, 1

        if (isp .eq. 0) then
          phasep = 1.0d0
        else
          phasep = -1.0d0
        endif

        do ib = 0, nb-1
          do ibp = 0, nb -1
            do ir = 0, nl - 1
              
              norm = norm + 
     $           (1.0d0/sqrt(2.0d0))*phasep*
     $           (cabs(psi(2*ib+1,2*ibp+isp,ir)))**2.0d0
                
            enddo
          enddo
        enddo

      enddo

      write(6,*) '<S=1 Sz=-1|psi|S=0 Sz=0> = ', norm / overall    
c

c Triplet-triplet (off-diagonal)

      norm = 0.0d0

      do isp = 0, 1
          
        if (isp .eq. 0) then
          phasep = 1.0d0
        else
          phasep = 1.0d0
        endif

        do ib = 0, nb-1
          do ibp = 0, nb -1
            do ir = 0, nl - 1
            
              norm = norm + 
     $           (cabs(psi(2*ib+0,2*ibp+0,ir)))**2.0d0
            
            enddo
          enddo
        enddo

      write(6,*) '<S=1 Sz=1|psi|S=1 Sz=1> = ', norm / overall


      do ib = 0, nb-1
        do ibp = 0, nb -1
          do ir = 0, nl - 1
            
            norm = norm + 
     $         (cabs(psi(2*ib+0,2*ibp+0,ir)))**2.0d0
            
          enddo
        enddo
      enddo

      write(6,*) '<S=1 Sz=1|psi|S=1 Sz=1> = ', norm / overall
      norm = 0.0d0

      do is = 0, 1

        if (is .eq. 0) then
          phase = 1.0d0
        else
          phase = 1.0d0
        endif

        do isp = 0, 1
          
          if (isp .eq. 0) then
            phasep = 1.0d0
          else
            phasep = 1.0d0
          endif

          do ib = 0, nb-1
            do ibp = 0, nb -1
              do ir = 0, nl - 1
              
                norm = norm + 
     $             (1.0d0/2.0d0)*phase*phasep*
     $             (cabs(psi(2*ib+is,2*ibp+isp,ir)))**2.0d0
                
              enddo
            enddo
          enddo

        enddo
      enddo

      write(6,*) '<S=0 Sz=0|psi|S=0 Sz=0> = ', norm / overall

      norm = 0.0d0

      do ib = 0, nb-1
        do ibp = 0, nb -1
          do ir = 0, nl - 1
            
            norm = norm + 
     $         (cabs(psi(2*ib+1,2*ibp+1,ir)))**2.0d0
            
          enddo
        enddo
      enddo

      write(6,*) '<S=1 Sz=-1|psi|S=1 Sz=-1> = ', norm / overall


      if ( (lly .eq. 1) .and. (llz .eq. 1) ) then  ! 1D
        call analyze_psi_1D(psi)
      endif


      if ( (lly .gt. 1) .and. (llz .eq. 1) ) then ! 2D

        if ( (lly .eq. llx) .and. (lcx .eq. lcy) .and.
     $     ( cabs(tij(0,0,1) - tij(0,0,llx)) .lt. 0.00001d0) ) then

          call analyze_psi_2D_tet(psi)
        
        else

          call analyze_psi_2D_ortho(psi)

        endif

      endif
c$$$
c$$$
c$$$      if ( (nly .gt. 1) .and. (nlz .gt. 1) ) then ! 3D
c$$$ 
c$$$        if ( (nly .eq. nlx) .and. (nx .eq. ny) .and.
c$$$     $     (abs(tx - ty) .lt. 0.00001d0) .and.
c$$$     $     (abs(txx - tyy) .lt. 0.00001d0) ) then !  x = y
c$$$
c$$$          if ( (nlz .eq. nlx) .and. (nz .eq. nx) .and.
c$$$     $       (abs(tx - tz) .lt. 0.00001d0) .and.
c$$$     $       (abs(txx - tzz) .lt. 0.00001d0) ) then ! cubic
c$$$
c$$$            call analyze_psi_3D_cubic(psi)
c$$$
c$$$          else                  ! tetragonal
c$$$
c$$$            call analyze_psi_3D_tet(psi)
c$$$
c$$$          endif
c$$$
c$$$        endif
c$$$
c$$$      endif
c$$$
c$$$
c$$$      write(6,*) 
c$$$      write(6,*) "Pair wave function "
c$$$      write(6,*) "i  j   ix   iy    iz    psi(i,j, ix, iy, iz)"
c$$$
c$$$      do i = 0, 1
c$$$        do j = 0, 1
c$$$          do ix = 0, nlx1
c$$$
c$$$            if (ix .gt. nlx/2) then
c$$$              jx = ix - nlx
c$$$            else
c$$$              jx = ix
c$$$            endif
c$$$
c$$$            do iy = 0, nly1
c$$$
c$$$              if (iy .gt. nly/2) then
c$$$                jy = iy - nly
c$$$              else
c$$$                jy = iy
c$$$              endif 
c$$$
c$$$              do iz = 0, nlz1
c$$$
c$$$                if (iz .gt. nlz/2) then
c$$$                  jz = iz - nlz
c$$$                else
c$$$                  jz = iz
c$$$                endif
c$$$
c$$$                if ( (ABS(jz) .le. 4) .and. (ABS(jy) .le. 4)
c$$$     $             .and. (ABS(jx) .le. 4) ) then
c$$$
c$$$                  write(6,*) i, " ", j, "  ", jx, " ", jy, " ",
c$$$     $              jz, "  ", psi(i,j,jx,jy,jz)
c$$$
c$$$                 endif
c$$$
c$$$              enddo
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$      enddo

      return
      end
