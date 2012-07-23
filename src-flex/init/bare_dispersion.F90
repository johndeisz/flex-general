#include "../convert.F90"

MODULE bare_dispersion

CONTAINS

  FUNCTION ekl(kl, tij, ed)

    USE CONSTANTS
    IMPLICIT NONE
    INTEGER kl
    COMPLEX, dimension (0:nb-1,0:nb-1) :: ekl
    COMPLEX, dimension (0:nb-1,0:nb-1,-2:2,-2:2,-2:2) :: tij   
    REAL, dimension (0:nb-1) :: ed

    REAL pi
    INTEGER k1, k2, k3
    REAL kx, ky, kz
    INTEGER ib, ibp

    INTEGER max_x, max_y, max_z

    COMPLEX ffx, ffy, ffz
    INTEGER ix, iy, iz

    pi = acos(-1.0d0)

    k3 = int(kl/(llx*lly))
    k2 = int( (kl-k3*llx*lly) / llx )
    k1 = kl - k3*llx*lly - k2*llx
    
    kx = float(k1) * (2.0d0*pi/float(llx))
    ky = float(k2) * (2.0d0*pi/float(lly))
    kz = float(k3) * (2.0d0*pi/float(llz))

    if (llx .gt. 2) then
       max_x = 2
    else
       max_x = llx - 1
    endif

    if (lly .gt. 2) then
       max_y = 2
    else
       max_y = lly - 1
    endif

    if (llz .gt. 2) then
       max_z = 2
    else
       max_z = llz - 1
    endif

    ekl = cmplx(0.0d0, 0.0d0)

    do ix = -max_x, max_x
       ffx = cexp(cmplx(0.0d0,-1.0d0)*kx*float(ix))

       do iy = -max_y, max_y
          ffy = cexp(cmplx(0.0d0,-1.0d0)*ky*float(iy))        

          do iz = -max_z, max_z
             ffz = cexp(cmplx(0.0d0,-1.0d0)*kz*float(iz))   

             ekl = ekl + ffx*ffy*ffz*tij(:,:,ix,iy,iz)


          enddo
       enddo
    enddo

    do ib = 0, nb-1
       ekl(ib,ib) = ekl(ib,ib) + ed(ib)
    enddo

    return

  END FUNCTION ekl

  REAL FUNCTION ek_minimum(tij, ed)
 
    USE CONSTANTS
    IMPLICIT NONE
    COMPLEX, dimension (0:nb-1,0:nb-1,-2:2,-2:2,-2:2) :: tij
    REAL, dimension (0:nb-1) :: ed
    COMPLEX, dimension (0:nb-1, 0:nb-1) :: ek
    INTEGER il, ib

    ! Determine the smallest value along the band-diagonal

    ek_minimum = 1.0d8
 
    do il = 0, nl-1
       ek = ekl(il, tij, ed)
       do ib = 0, nb-1 
          if ( real( ek(ib,ib) ) .lt. ek_minimum) then
             ek_minimum = real( ek(ib,ib) )
          endif
       enddo
    enddo
    return
  END FUNCTION ek_minimum

END MODULE bare_dispersion
