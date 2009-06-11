      program oneband
  
      implicit none

      DOUBLE PRECISION t(0:1,0:1,-2:2,-2:2), ed(0:1), flux(1:3)

      integer ix, iy, i1, i2
      DOUBLE PRECISION t1, t2, t3, t4
      
      t  = 0.0d0

      read(5,*) flux(1), flux(2), flux(3)
      read(5,*) ed(0), ed(1)
      read(5,*) t1, t2, t3, t4

      t(0,0,1,0) = -t1
      t(0,0,-1,0) = -t1
      t(0,0,0,1) = -t2
      t(0,0,0,-1) = -t2
      t(0,0,1,1) = -t3
      t(0,0,1,-1) = -t3
      t(0,0,-1,1) = -t3
      t(0,0,-1,-1) = -t3

      t(1,1,1,0) = -t2
      t(1,1,-1,0) = -t2
      t(1,1,0,1) = -t1
      t(1,1,0,-1) = -t1
      t(1,1,1,1) = -t3
      t(1,1,1,-1) = -t3
      t(1,1,-1,1) = -t3
      t(1,1,-1,-1) = -t3

      t(0,1,1,1) = t4
      t(0,1,1,-1) = -t4
      t(0,1,-1,1) = -t4
      t(0,1,-1,-1) = t4

      t(1,0,1,1) = t4
      t(1,0,1,-1) = -t4
      t(1,0,-1,1) = -t4
      t(1,0,-1,-1) = t4

      write(6,*) '--------------flux values --------------'
      write(6,*) flux(1), flux(2), flux(3)
      write(6,*) 
      write(6,*) '--------------orbital level-------------'
      write(6,*) 0, ed(0)
      write(6,*) 1, ed(1)

      do ix = -2,2
        do iy = -2,2
          write(6,*) 
          write(6,200) ix,iy
          do i1 = 0, 1
           do i2 = 0, 1 
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
