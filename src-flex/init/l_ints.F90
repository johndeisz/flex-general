#include "../convert.F90"

subroutine l_ints(l_integral, t, x, y)

  USE CONSTANTS
  IMPLICIT NONE

  REAL t
  REAL x(0:1, 0:1), y(0:1, 0:1)
  REAL l_integral(0:1,0:1,0:1,0:1,0:1,0:1), l_help

  REAL ap(0:1,0:1,0:1), am(0:1,0:1,0:1), b(0:1,0:1,0:1)
  REAL ex1, ex3, ey2, fx1, fx3, fy2
  INTEGER i1, j1, i2, j2, i3, j3, k1, k2, k3

  REAL sign1, sign2, sign3
  REAL left1, left2, left3
  REAL right1, right2, right3

  INTEGER i

  do i = 0, 1

     ap(0,0,i) = -0.5d0
     ap(1,0,i) = -0.5d0

     ap(0,1,i) = 0.5d0 / x(1,i)
     ap(1,1,i) = -0.5d0 / x(1,i)

     am(0,0,i) = 0.5d0
     am(1,0,i) = 0.5d0

     am(0,1,i) = 0.5d0 / x(1,i)
     am(1,1,i) = -0.5d0 / x(1,i)

     b(0,0,i) = -0.5d0
     b(1,0,i) = -0.5d0

     b(0,1,i) = -0.5d0 / y(1,i)
     b(1,1,i) = 0.5d0 / y(1,i)

  enddo

  do i1 = 0, 1
     do j1 = 0, 1
 
        ex1 = exp(-x(i1,j1) / t)
        fx1 = 1.0d0 / (1.0d0 + ex1)

        do i2 = 0, 1
           do j2 = 0, 1

              ey2 = exp(-y(i2,j2) / t)
              fy2 = 1.0d0 / (1.0d0 - ey2)

              do i3 = 0, 1
                 do j3 = 0, 1

                    ex3 = exp(-x(i3,j3) / t)
                    fx3 = 1.0d0 / (1.0d0 + ex3)

                    l_help = 0.0d0

                    do k1 = 0, 1

                       if (k1 .eq. 0) then
                          left1 = 1.0d0
                          right1 = ex1
                          sign1 = 1.0d0
                       else
                          left1 = ex1
                          right1 = 1.0d0
                          sign1 = -1.0d0
                       endif

                       do k2 = 0, 1
           
                          if (k2 .eq. 0) then
                             left2 = ey2
                             right2 = 1.0d0
                             sign2 = -1.0d0 
                          else 
                             left2 = 1.0d0
                             right2 = ey2
                             sign2 = 1.0d0
                          endif

                          do k3 = 0, 1
           
                             if (k3 .eq. 0) then
                                left3 = 1.0
                                right3 = ex3
                                sign3 = 1.0
                             else 
                                left3 = ex3
                                right3 = 1.0
                                sign3 = -1.0
                             endif

                             l_help = l_help +  ap(k1,i1,j1) * b(k2,i2,j2) * & 
                                  am(k3,i3,j3) * ( -sign2 ) * & 
                                  ( left1*left2*left3 - &
                                  right1*right2*right3 ) / &
                                  ( sign1*x(i1,j1) + sign2*y(i2,j2) + & 
                                  sign3*x(i3,j3) )

                          enddo
                       enddo
                    enddo

                    l_integral(i1,j1,i2,j2,i3,j3) =  l_help  * fx1 * fy2 * fx3

                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo

  return
end subroutine l_ints
