      program oneband

c By assumption the unit cell contains n atoms along the x-direction
c and 1 atom along the y-direction
  
      implicit none

      DOUBLE PRECISION t(0:1,0:1,-2:2,-2:2), ed(0:1), flux(1:3)

      integer ix, iy, i1, i2
      DOUBLE PRECISION tx, ty, txy, txx, tyy
      
      t  = 0.0d0

      read(5,*) flux(1), flux(2), flux(3)
      read(5,*) ed(0), ed(1)
      read(5,*) tx
      read(5,*) ty
      read(5,*) txy
      read(5,*) txx
      read(5,*) tyy

      t(0,0,1,0) = -txx
      t(0,0,-1,0) = -txx
      t(0,0,0,1) = -ty
      t(0,0,0,-1) = -ty
      t(0,0,0,2) = -tyy
      t(0,0,0,-2) = -tyy

      t(1,1,1,0) = -txx
      t(1,1,-1,0) = -txx
      t(1,1,0,1) = -ty
      t(1,1,0,-1) = -ty
      t(1,1,0,2) = -tyy
      t(1,1,0,-2) = -tyy

      t(1,0,-1,0) = -tx
      t(0,1,1,0) = -tx
      t(1,0,0,0) = -tx
      t(0,1,0,0) = -tx
  
      t(1,0,0,1) = -txy
      t(1,0,0,-1) = -txy
      t(1,0,1,1) = -txy
      t(1,0,1,1) = -txy

      t(0,1,0,1) = -txy
      t(0,1,0,-1) = -txy
      t(0,1,-1,1) = -txy
      t(0,1,-1,-1) = -txy

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
