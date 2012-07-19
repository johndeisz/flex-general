      program BiS2_2band
  
      implicit none

      DOUBLE PRECISION t(0:1,0:1,-2:2,-2:2), ed(0:1), flux(1:3)

      integer ix, iy, i1, i2
      DOUBLE PRECISION t1, t2, t3, t4
      
      t  = 0.0d0
      flux = 0.0d0

      ed(0) = 2.811d0
      ed(1) = 2.811d0

      t(0,0,1,0) = -0.167d0
      t(0,0,1,-1) = 0.880d0
      t(0,0,1,1) = 0.094d0
      t(0,0,2,1) = 0.014d0
      t(0,0,2,-1) = 0.069d0

      t(0,1,1,0) = 0.107d0
      t(0,1,2,0) = -0.028d0
      t(0,1,2,1) = 0.020d0
      t(0,1,2,-1) = 0.020d0

      t(1,0,1,0) = 0.107d0
      t(1,0,2,0) = -0.028d0
      t(1,0,2,1) = 0.020d0
      t(1,0,2,-1) = 0.020d0

      t(1,1,1,0) = -0.167d0
      t(1,1,1,-1) = 0.094d0
      t(1,1,1,1) = 0.880d0
      t(1,1,2,1) = 0.069d0
      t(1,1,2,-1) = 0.014d0    

      do ix = 1,2
        do iy = -1,1

          t(0,0,-ix,-iy) = t(0,0,ix,iy)
          t(1,1,-ix,-iy) = t(1,1,ix,iy)
          t(0,1,-ix,-iy) = t(0,1,ix,iy)
          t(1,0,-ix,-iy) = t(1,0,ix,iy)

        enddo
      enddo

      do ix = -1,1
        do iy = 1,2
          t(0,0,ix,iy) = t(0,0,iy,ix)
          t(1,1,ix,iy) = t(1,1,iy,ix)
          t(0,1,ix,iy) = -t(0,1,iy,ix)
          t(1,0,ix,iy) = -t(1,0,iy,ix)
        enddo
      enddo

      do ix = -1,1
        do iy = -2,-1
          t(0,0,ix,iy) = t(0,0,iy,ix)
          t(1,1,ix,iy) = t(1,1,iy,ix)
          t(0,1,ix,iy) = -t(0,1,iy,ix)
          t(1,0,ix,iy) = -t(1,0,iy,ix)
        enddo
      enddo


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
