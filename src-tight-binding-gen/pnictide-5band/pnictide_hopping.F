      program pnictide

      implicit none

      DOUBLE COMPLEX, dimension (1:5,1:5,-2:2,-2:2) :: tij
      DOUBLE PRECISION, dimension (1:5, 1:5, -2:2, -2:2) :: tin
      DOUBLE PRECISION, dimension (1:5) :: ed
      DOUBLE PRECISION, dimension (1:3) :: flux

      INTEGER, dimension (1:5,1:5) :: lookup_l, lookup_r, sym_factor

      INTEGER ib, jb, ix, iy, il, ind_old, ind_new
      INTEGER nb

      nb=5

      read(5,*) flux(1), flux(2), flux(3)
      read(5,*)
      read(5,*) ed(1), ed(2), ed(3), ed(4), ed(5)

      tij  = cmplx(0.0d0, 0.0d0)
      tin = 0.0d0

      do ix = 1, 2
        do iy = 0, ix

          read(5,*)
          do ib = 1, 5
            read(5,*) tin(ib,ib:5,ix,iy)
          enddo

        enddo
      enddo

      tij = 0.1d0*dcmplx(tin, 0.0d0)


!Use symmetries to fill out the hopping matrix
!------------ [Dx,Dy]->[Dx,-Dy] phase factors/transpositions ---------
      lookup_l(1,1) = 1
      lookup_r(1,1) = 1
      sym_factor(1,1) = 1

      lookup_l(1,2) = 1
      lookup_r(1,2) = 3
      sym_factor(1,2) = -1
      
      lookup_l(1,3) = 1
      lookup_r(1,3) = 2
      sym_factor(1,3) = -1

      lookup_l(1,4) = 1
      lookup_r(1,4) = 4
      sym_factor(1,4) = -1

      lookup_l(1,5) = 1
      lookup_r(1,5) = 5
      sym_factor(1,5) = 1

      lookup_l(2,2) = 3
      lookup_r(2,2) = 3
      sym_factor(2,2) = 1
      
      lookup_l(2,3) = 2
      lookup_r(2,3) = 3
      sym_factor(2,3) = 1

      lookup_l(2,4) = 3
      lookup_r(2,4) = 4
      sym_factor(2,4) = 1

      lookup_l(2,5) = 3
      lookup_r(2,5) = 5
      sym_factor(2,5) = -1

      lookup_l(3,3) = 2
      lookup_r(3,3) = 2
      sym_factor(3,3) = 1

      lookup_l(3,4) = 2
      lookup_r(3,4) = 4
      sym_factor(3,4)  = 1

      lookup_l(3,5) = 2
      lookup_r(3,5) = 5
      sym_factor(3,5) = -1

      lookup_l(4,4) = 4
      lookup_r(4,4) = 4
      sym_factor(4,4) = 1

      lookup_l(4,5) = 4
      lookup_r(4,5) = 5
      sym_factor(4,5) = -1

      lookup_l(5,5) = 5
      lookup_r(5,5) = 5
      sym_factor(5,5) = 1

      do ix = 1, 2
        do iy = 0, ix
          do ib = 1, nb
            do jb = ib, nb
              tij(ib, jb, ix, -iy) =
     $           dfloat( sym_factor(ib,jb) ) *
     $           tij(lookup_l(ib,jb), lookup_r(ib,jb), ix, iy)
            enddo
          enddo
        enddo
      enddo

!------------ [Dx,Dy] -> [-Dx,-Dy] phase factors -----------------

      sym_factor(1,1) = 1
      sym_factor(1,2) = -1
      sym_factor(1,3) = -1
      sym_factor(1,4) = 1
      sym_factor(1,5) = 1

      sym_factor(2,2) = 1
      sym_factor(2,3) = 1
      sym_factor(2,4) = -1
      sym_factor(2,5) = -1

      sym_factor(3,3) = 1
      sym_factor(3,4)  = -1
      sym_factor(3,5) = -1

      sym_factor(4,4) = 1
      sym_factor(4,5) = 1

      sym_factor(5,5) = 1

      do ix = 1, 2
        do iy = -ix, ix
          do ib = 1, nb
            do jb = ib, nb
              tij(ib,jb,-ix,-iy) =
     $           float( sym_factor(ib,jb) ) * tij(ib,jb,ix,iy)
            enddo
          enddo
        enddo
      enddo

!------------ [Dx,Dy] -> [Dy,Dx] phase factors -----------------

      sym_factor(1,1) = 1
      sym_factor(1,2) = -1
      sym_factor(1,3) = 1
      sym_factor(1,4) = 1
      sym_factor(1,5) = -1

      sym_factor(2,2) = 1
      sym_factor(2,3) = -1
      sym_factor(2,4) = -1
      sym_factor(2,5) = 1

      sym_factor(3,3) = 1
      sym_factor(3,4)  = 1
      sym_factor(3,5) = -1

      sym_factor(4,4) = 1
      sym_factor(4,5) = -1

      sym_factor(5,5) = 1

      do ix = -2, 2
        do iy = -abs(ix), abs(ix)
          do ib = 1, nb
            do jb = ib, nb
              tij(ib,jb,iy,ix) =
     $           float( sym_factor(ib,jb) ) * tij(ib,jb,ix,iy)
            enddo
          enddo
        enddo
      enddo

! Complete the lower diagonal of the hopping matrix

      do ix = -2, 2
        do iy = -2, 2
          do ib = 1, nb
            do jb = ib, nb
              tij(jb,ib,-ix,-iy) = tij(ib,jb,ix,iy)
            enddo
          enddo
        enddo
      enddo          

! Write output

      write(6,*) '--------------flux values --------------'
      write(6,*) flux(1), flux(2), flux(3)
      write(6,*)
      write(6,*) '--------------orbital level-------------'
      do ib=1,5
        write(6,*) ib-1, ed(ib)
      enddo

      do ix = -2, 2
        do iy = -2, 2

          write(6,*)
          write(6,200) ix,iy
                  
          do ib = 1,5
            do jb = 1,5

              write(6,300) ib, jb, dreal(tij(ib,jb,ix,iy)),
     $           dimag(tij(ib,jb,ix,iy))  

            enddo
          enddo

        enddo
      enddo

 200  format('------------------- [',i3,',',i3,'] hopping',
     $   '------------------------')
 300  format(i3,',',i3,'  ','(',D16.9,',',D16.9,')')


      return
      end







