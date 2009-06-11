#include "../convert.F"

      subroutine parameter_init(x, c_r, delta_g_r, delta_gp_r, y)  

#include "../constants.F"

      REAL x(0:1,0:1) 
      COMPLEX c_r(0:1,0:1,0:4*nb-1,0:4*nb-1,0:nc-1)
      COMPLEX delta_g_r(0:4*nb-1,0:4*nb-1,0:nc-1)
      COMPLEX delta_gp_r(0:4*nb-1,0:4*nb-1,0:nc-1)
      REAL y(0:1,0:1)

c     Set x[i][j] and y[i][j].  
c     For simplicity in subsequent manipulations,
c     they should assume unique values.  These will remain constant
c     throughout the calculation.

      x(0,0) = 1.5d0
      x(0,1) = 2.5d0
      x(1,0) = 2.0d0
      x(1,1) = 3.0d0

      y(0,0) = 1.25d0
      y(0,1) = 2.25d0
      y(1,0) = 1.75d0
      y(1,1) = 2.75d0   

c     Initialize c_r so that the asymptotic behvior of g is accounted for.

      c_r(0,0,:,:,:) = -delta_g_r
      c_r(0,1,:,:,:) = cmplx(0.0d0, 0.0d0)
      c_r(1,0,:,:,:) = delta_gp_r
      c_r(1,1,:,:,:) = cmplx(0.0d0, 0.0d0)

      return
      end



