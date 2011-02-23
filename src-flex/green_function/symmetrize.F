#include "../convert.F"

      subroutine symmetrize(rank,g, g_mtau)

#include "../constants.F"

      INTEGER rank
      COMPLEX g(0:4*nb-1,0:4*nb-1,0:mp1,0:nc1) 
      COMPLEX g_mtau(0:4*nb-1,0:4*nb-1,0:mp1,0:nc1)

      INTEGER l, ir, irm, nua, nuap, nuam, nuapm
      INTEGER nua2, nua2p
      INTEGER ib,is, ibp, isp
      INTEGER ix, iy, iz

      COMPLEX A, S

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

             A = 0.5d0 * (g(nua, nuap, l, ir) -
     $          g_mtau(nuapm,nuam, l, ir) )

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

             A = 0.5d0 * (g(nua, nuap, l, ir) -
     $          g_mtau(nuapm,nuam, l, ir) )

             g(nua, nuap, l, ir) = A
             g_mtau(nuapm, nuam, l, ir) = -A

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
           irm = mod(lcx-ix,lcx) + mod(lcy-iy,lcy)*lcx +
     $           mod(lcz-iz,lcz)*lcx*lcy
 
           do ib = 0, nb-1
            do is = 0, 1
             do ibp = 0, nb-1
               do isp = 0, 1

               nua = 4*ib+is
               nuap = 4*ibp+isp+2

               nuam = 4*ib+is+2
               nuapm = 4*ibp+isp

             A = 0.5d0 * (g(nua, nuap, 0, ir) -
     $          g(nuapm,nuam, 0, irm) )

             g(nua, nuap, 0, ir) = A
             g(nuapm, nuam, 0, irm) = -A

             A = 0.5d0 * (g_mtau(nua, nuap, 0, ir) -
     $          g_mtau(nuapm,nuam, 0, irm) )

             g_mtau(nua, nuap, 0, ir) = A
             g_mtau(nuapm, nuam, 0, irm) = -A


         enddo
         enddo
         enddo
        enddo 
           do ib = 0, nb-1
            do is = 0, 1
             do ibp = 0, nb-1
               do isp = 0, 1

               nua = 4*ib+is+2
               nuap = 4*ibp+isp

               nuam = 4*ib+is
               nuapm = 4*ibp+isp+2

             A = 0.5d0 * (g(nua, nuap, 0, ir) -
     $          g(nuapm,nuam, 0, irm) )

             g(nua, nuap, 0, ir) = A
             g(nuapm, nuam, 0, irm) = -A

             A = 0.5d0 * (g_mtau(nua, nuap, 0, ir) -
     $          g_mtau(nuapm,nuam, 0, irm) )

             g_mtau(nua, nuap, 0, ir) = A
             g_mtau(nuapm, nuam, 0, irm) = -A


         enddo
         enddo
         enddo
        enddo 

       enddo
       enddo
       enddo
      endif

c Project out triplet symmetry

      if (rank .eq. 0) then
       do ix = 0, lcx1
        do iy = 0, lcy1
         do iz = 0, lcz1
           ir = ix + iy*lcx + iz*lcy*lcx


           do ib = 0, nb-1
            do is = 0, 1
             do ibp = 0, nb-1
               do isp = 0, 1

               nua = 4*ib+is
               nuap = 4*ibp+isp+2

               nua2 = 4*ib+isp
               nua2p = 4*ibp + is+2

               S = 0.5d0 * (g(nua, nuap, 0, ir) +
     $          g(nua2,nua2p, 0, ir) )

               g(nua, nuap, 0, ir) = S
               g(nua2, nua2p, 0, ir) = S

               S = 0.5d0 * (g_mtau(nua, nuap, 0, ir) +
     $          g_mtau(nua2,nua2p, 0, ir) )

               g_mtau(nua, nuap, 0, ir) = S
               g_mtau(nua2, nua2p, 0, ir) = S

               nua = 4*ib+is+2
               nuap = 4*ibp+isp

               nua2 = 4*ib+isp+2
               nua2p = 4*ibp + is

               S = 0.5d0 * (g(nua, nuap, 0, ir) +
     $          g(nua2,nua2p, 0, ir) )


               g(nua, nuap, 0, ir) = S
               g(nua2, nua2p, 0, ir) = S


               S = 0.5d0 * (g_mtau(nua, nuap, 0, ir) +
     $          g_mtau(nua2,nua2p, 0, ir) )


               g_mtau(nua, nuap, 0, ir) = S
               g_mtau(nua2, nua2p, 0, ir) = S

            enddo
            enddo
           enddo
           enddo
           enddo
          enddo
 
         enddo

       endif

      return
      end
