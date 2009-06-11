      program oneband
  
      implicit none

      DOUBLE PRECISION t(-2:2,-2:2), ed0, flux(1:3)

      integer ix, iy
      DOUBLE PRECISION hop
      
      t  = 0.0d0

      read(5,*) flux(1), flux(2), flux(3)
      read(5,*) ed0
      read(5,*) 
      read(5,*) t(1,0)
      read(5,*) 
      read(5,*) t(1,1)
      read(5,*)
      read(5,*) t(2,0)
      read(5,*)
      read(5,*) t(2,1)
      read(5,*)
      read(5,*) t(2,2)

      write(6,*) '--------------flux values --------------'
      write(6,*) flux(1), flux(2), flux(3)
      write(6,*) 
      write(6,*) '--------------orbital level-------------'
      write(6,*) 0, ed0

      do iy = 1, 2
       do ix = 0, iy-1
        t(ix,iy) = t(iy,ix)
       enddo
      enddo
      
      do ix = -2,2
        do iy = -2,2
          hop = t(abs(ix),abs(iy))

          write(6,*) 
          write(6,200) ix,iy
          write(6,300) 0,0,hop,0.0d0

        enddo
      enddo

 200  format('------------------- [',i3,',',i3,'] hopping',
     $   '------------------------')
 300  format(i3,',',i3,'  ','(',D16.9,',',D16.9,')')


      stop
      end
