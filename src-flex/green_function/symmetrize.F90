#include "../convert.F90"

subroutine symmetrize(rank,g, g_mtau)

  USE CONSTANTS
  IMPLICIT NONE

  INTEGER rank
  COMPLEX g(0:4*nb-1,0:4*nb-1,0:mp1,0:nc1) 
  COMPLEX g_mtau(0:4*nb-1,0:4*nb-1,0:mp1,0:nc1)

  INTEGER l, ir, irm, nua, nuap, nuam, nuapm
  INTEGER nua2, nua2p
  INTEGER ib,is, ibp, isp, ip, ipp
  INTEGER ix, iy, iz

  COMPLEX A, S, B
  INTEGER i,j,k
  COMPLEX sum_test

  do l = 0, mp1
     do ir = 0, nc1
        do ib = 0, nb-1
           do is = 0, 1
              do ibp = 0, nb-1
                 do isp = 0, 1

                    nua = 4*ib+is
                    nuap = 4*ibp+isp+2

                    nuam = 4*ib+is+2
                    nuapm = 4*ibp+isp

                    A = 0.5d0 * (g(nua, nuap, l, ir) - &
                         g_mtau(nuapm,nuam, l, ir) )

                    g(nua, nuap, l, ir) = A
                    g_mtau(nuapm, nuam, l, ir) = -A

                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo

  do l = 0, mp1
     do ir = 0, nc1
        do ib = 0, nb-1
           do is = 0, 1
              do ibp = 0, nb-1
                 do isp = 0, 1

                    nua = 4*ib+is+2
                    nuap = 4*ibp+isp

                    nuam = 4*ib+is
                    nuapm = 4*ibp+isp+2

                    A = 0.5d0 * (g(nua, nuap, l, ir) - &
                         g_mtau(nuapm,nuam, l, ir) )

                    g(nua, nuap, l, ir) = A
                    g_mtau(nuapm, nuam, l, ir) = -A

                    nua = 4*ib+is
                    nuap = 4*ibp+isp

                    nuam = 4*ib+is+2
                    nuapm = 4*ibp+isp+2

                    A = 0.5d0 * (g(nua, nuap, l, ir) - &
                         g_mtau(nuapm,nuam, l, ir) )

                    g(nua, nuap, l, ir) = A
                    g_mtau(nuapm, nuam, l, ir) = -A

                    nua = 4*ib+is+2
                    nuap = 4*ibp+isp+2

                    nuam = 4*ib+is
                    nuapm = 4*ibp+isp

                    A = 0.5d0 * (g(nua, nuap, l, ir) - &
                         g_mtau(nuapm,nuam, l, ir) )

                    g(nua, nuap, l, ir) = A
                    g_mtau(nuapm, nuam, l, ir) = -A

                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo

  do ix = 0, lcx1
   do iy = 0, lcy1
     do iz = 0, lcz1

          ir = ix + iy*lcx + iz*lcy*lcx
          irm = mod(lcx-ix,lcx) + mod(lcy-iy,lcy)*lcx + &
                   mod(lcz-iz,lcz)*lcx*lcy
 
         do l = 0, mp1
           do nua = 0, 4*nb-1
             do nuap = 0, 4*nb-1

               A = 0.5d0 * ( g(nua,nuap, l, ir) + &
                        conjg( g(nuap, nua, l, irm) ))

               g(nua, nuap, l, ir) = A
               g(nuap, nua, l, irm) = conjg(A)

               A = 0.5d0 * ( g_mtau(nua,nuap, l, ir) + &
                        conjg( g_mtau(nuap, nua, l, irm) ))

               g_mtau(nua, nuap, l, ir) = A
               g_mtau(nuap, nua, l, irm) = conjg(A)

            enddo
          enddo
        enddo

     enddo
    enddo
  enddo 

  if (rank .eq. 0) then
     do ix = 0, lcx1
        do iy = 0, lcy1
           do iz = 0, lcz1

              ir = ix + iy*lcx + iz*lcy*lcx
              irm = mod(lcx-ix,lcx) + mod(lcy-iy,lcy)*lcx + &
                   mod(lcz-iz,lcz)*lcx*lcy
 
              do ib = 0, nb-1
                 do is = 0, 1
                   do ip = 0, 1
                    do ibp = 0, nb-1
                       do isp = 0, 1
                         do ipp = 0, 1

                          nua = 4*ib+2*ip+is
                          nuap = 4*ibp+2*ipp+isp

                          nuam = 4*ib+2*mod(ip+1,2)+is
                          nuapm = 4*ibp+2*mod(ipp+1,2)+isp

                          A = 0.5d0 * (g(nua, nuap, 0, ir) - &
                               g(nuapm,nuam, 0, irm) )

                          if (  (nua .eq. nuap) .and. (ir .eq. 0) ) then
                            B = 0.5
                          else
                            B = 0.0
                          endif

                          g(nua, nuap, 0, ir) = A - B
                          g(nuapm, nuam, 0, irm) = -A - B

                          A = 0.5d0 * (g_mtau(nua, nuap, 0, ir) - &
                               g_mtau(nuapm,nuam, 0, irm) )

                          g_mtau(nua, nuap, 0, ir) = A + B
                          g_mtau(nuapm, nuam, 0, irm) = -A + B

                       enddo
                    enddo
                 enddo
              enddo
              enddo
              enddo

           enddo
        enddo
     enddo

  endif


!!$c Project out triplet symmetry
!!$
!!$c$$$      if (rank .eq. 0) then
!!$c$$$       do ix = 0, lcx1
!!$c$$$        do iy = 0, lcy1
!!$c$$$         do iz = 0, lcz1
!!$c$$$           ir = ix + iy*lcx + iz*lcy*lcx
!!$c$$$
!!$c$$$
!!$c$$$           do ib = 0, nb-1
!!$c$$$            do is = 0, 1
!!$c$$$             do ibp = 0, nb-1
!!$c$$$               do isp = 0, 1
!!$c$$$
!!$c$$$               nua = 4*ib+is
!!$c$$$               nuap = 4*ibp+isp+2
!!$c$$$
!!$c$$$               nua2 = 4*ib+isp
!!$c$$$               nua2p = 4*ibp + is+2
!!$c$$$
!!$c$$$               S = 0.5d0 * (g(nua, nuap, 0, ir) +
!!$c$$$     $          g(nua2,nua2p, 0, ir) )
!!$c$$$
!!$c$$$               g(nua, nuap, 0, ir) = S
!!$c$$$               g(nua2, nua2p, 0, ir) = S
!!$c$$$
!!$c$$$               S = 0.5d0 * (g_mtau(nua, nuap, 0, ir) +
!!$c$$$     $          g_mtau(nua2,nua2p, 0, ir) )
!!$c$$$
!!$c$$$               g_mtau(nua, nuap, 0, ir) = S
!!$c$$$               g_mtau(nua2, nua2p, 0, ir) = S
!!$c$$$
!!$c$$$               nua = 4*ib+is+2
!!$c$$$               nuap = 4*ibp+isp
!!$c$$$
!!$c$$$               nua2 = 4*ib+isp+2
!!$c$$$               nua2p = 4*ibp + is
!!$c$$$
!!$c$$$               S = 0.5d0 * (g(nua, nuap, 0, ir) +
!!$c$$$     $          g(nua2,nua2p, 0, ir) )
!!$c$$$
!!$c$$$
!!$c$$$               g(nua, nuap, 0, ir) = S
!!$c$$$               g(nua2, nua2p, 0, ir) = S
!!$c$$$
!!$c$$$
!!$c$$$               S = 0.5d0 * (g_mtau(nua, nuap, 0, ir) +
!!$c$$$     $          g_mtau(nua2,nua2p, 0, ir) )
!!$c$$$
!!$c$$$
!!$c$$$               g_mtau(nua, nuap, 0, ir) = S
!!$c$$$               g_mtau(nua2, nua2p, 0, ir) = S
!!$c$$$
!!$c$$$            enddo
!!$c$$$            enddo
!!$c$$$           enddo
!!$c$$$           enddo
!!$c$$$           enddo
!!$c$$$          enddo
!!$c$$$ 
!!$c$$$         enddo
!!$c$$$
!!$c$$$       endif

  return
end subroutine symmetrize
