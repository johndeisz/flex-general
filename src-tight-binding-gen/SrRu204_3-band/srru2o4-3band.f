      program threeband

c Parameterization taken from Pavarini and Mazin, PRB 74 035115 (2006)
  
      implicit none

      DOUBLE PRECISION t(0:2,0:2,-2:2,-2:2), ed(0:2), flux(1:3)
      DOUBLE PRECISION hop(0:5,0:2,0:2)

      integer ix, iy, i1, i2
      DOUBLE PRECISION t1, t2, t3, t4
      
      t  = 0.0d0
      hop = 0.0d0

      read(5,*) flux(1), flux(2), flux(3)
      read(5,*) ed(0), ed(1), ed(2)
      read(5,*) hop(:,1,0)
      read(5,*) hop(:,0,1)
      read(5,*) hop(:,1,1)
      read(5,*) hop(:,2,0)
      read(5,*) hop(:,0,2)

c Convert to eV

      ed = ed * 13.605698066d0 / 1000.0d0
      hop = hop * 13.605698066d0 / 1000.0d0

      do ix = -2, 2
        do iy = -2, 2

           t(0,0,ix,iy) = hop(0,abs(ix),abs(iy))
           t(1,1,ix,iy) = hop(1,abs(ix),abs(iy))
           t(2,2,ix,iy) = hop(2,abs(ix),abs(iy))
           t(0,1,ix,iy) = hop(3,abs(ix),abs(iy))
           t(1,0,ix,iy) = t(0,1,ix,iy)
           t(0,2,ix,iy) = hop(4,abs(ix),abs(iy))
           t(2,0,ix,iy) = t(0,2,ix,iy)
           t(1,2,ix,iy) = hop(5,abs(ix),abs(iy))
           t(2,1,ix,iy) = t(1,2,ix,iy)

         enddo
        enddo

        t(0,1,-1,1) = -t(0,1,-1,1)
        t(1,0,-1,1) = -t(1,0,-1,1)

        t(0,1,1,-1) = -t(0,1,1,-1)
        t(1,0,1,-1) = -t(1,0,1,-1)

      write(6,*) '--------------flux values --------------'
      write(6,*) flux(1), flux(2), flux(3)
      write(6,*) 
      write(6,*) '--------------orbital level-------------'
      write(6,*) 0, ed(0)
      write(6,*) 1, ed(1)
      write(6,*) 2, ed(2)

      do ix = -2,2
        do iy = -2,2
          write(6,*) 
          write(6,200) ix,iy
          do i1 = 0, 2
           do i2 = 0, 2 
          write(6,300) i1,i2,t(i1,i2,ix,iy),0.0d0
          enddo
          enddo
        enddo
      enddo

 200  format('------------------- [',i3,',',i3,'] hopping',
     $   '------------------------')
 300  format(i3,',',i3,'  ','(',D16.9,',',D16.9,')')


      stop
      end
