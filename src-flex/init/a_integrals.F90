#include "../convert.F90"

subroutine a_integrals( t, x, y, epsilon, a_int)

  USE CONSTANTS
  IMPLICIT NONE

  REAL t, x(0:1,0:1), y(0:1,0:1), epsilon(0:mp1)
  COMPLEX a_int(0:1,0:1,0:1,0:1,0:mp1)

  REAL a(0:1,0:1,0:1), b(0:1,0:1,0:1)

  REAL fx, ex, nby, ey
  INTEGER ia, ib, ja, jb

  COMPLEX part1(0:mp1), part2(0:mp1)
  COMPLEX part3(0:mp1), part4(0:mp1)

  INTEGER i, l

  !     Define constants appearing in q and r functions.

  do i = 0,1
     
     a(0,0,i) = -0.5d0
     a(1,0,i) = -0.5d0
     a(0,1,i) = 0.5d0 / x(1,i)
     a(1,1,i) = -0.5d0 / x(1,i)
        
     b(0,0,i) = -0.5d0
     b(1,0,i) = -0.5d0
     b(0,1,i) = -0.5d0 / y(1,i)
     b(1,1,i) = 0.5d0 / y(1,i)

  enddo


  !     Evaluate the integrals for the Fourier transform of  Q(-tau) R(tau)

  do ia = 0, 1
     do ib = 0, 1

        fx = 1.0d0 / ( 1.0d0 + exp( -x(ia,ib) / t) )
        ex = exp( -x(ia,ib) / t)

        do ja = 0, 1
           do jb = 0, 1
 
              nby = 1.0d0 / ( exp( -y(ja,jb) / t ) - 1.0d0 )
              ey = exp( -y(ja,jb) / t )
              
              do l = 0, mp1

                 part1(l) = (ex + ey) * a(0,ia,ib) * b(0,ja,jb) / & 
                      ( cmplx( 0.0d0, epsilon(l) ) +  ( x(ia,ib) - y(ja,jb) ) )
	      
                 part2(l) = (ex*ey + 1.0d0) * &
                      a(1,ia,ib) * b(0,ja,jb) / & 
                      ( cmplx( 0.0d0, epsilon(l) ) + & 
                      ( -x(ia,ib) - y(ja,jb) ) )

                 part3(l) = -( ex*ey + 1.0d0 ) * a(0,ia,ib) * b(1,ja,jb) / & 
                      ( cmplx( 0.0d0,epsilon(l) ) + ( x(ia,ib) + y(ja,jb) ) )
 
                 part4(l) = -(ex + ey) * a(1,ia,ib) * b(1,ja,jb) / & 
                      ( cmplx( 0.0d0, epsilon(l) ) + ( -x(ia,ib) + y(ja,jb) ) )

                 a_int(ia,ib,ja,jb,l) = fx * nby * &
                      ( part1(l) +  part2(l) + part3(l) + part4(l) )

              enddo
              
           enddo
        enddo

     enddo
  enddo


  return
end subroutine a_integrals
