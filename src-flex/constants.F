      implicit none

c Lattice and frequency grid parameters. (all must be 2^n!)
      integer nb, m  
      integer lcx, lcy, lcz, nc 
      integer llx, lly, llz, nl

      parameter (nb = 2)

      parameter (m =  128) 
      parameter (lcx = 4)
      parameter (lcy = 2)
      parameter (lcz = 1)
      parameter (nc = lcx*lcy*lcz)  ! number of cluster sites

      parameter (llx = 16)     ! lattice dimension >=lcx 
      parameter (lly = 32)
      parameter (llz = 1)
      parameter (nl = llx*lly*llz)

c MPI related constants.

      integer np, mp, ncp, max_message

      parameter (np = 4)
      parameter (mp = m/np)  ! must be an integer
      parameter (ncp = nc/np)  ! must be an integer

c Convenience constants

      integer mp1, m1
      integer llx1, lly1, llz1
      integer ncp1, nc1
      integer lcx1, lcy1, lcz1

      parameter (mp1 = mp - 1)
      parameter (m1 = m - 1)

      parameter (llx1 = llx - 1)
      parameter (lly1 = lly - 1)
      parameter (llz1 = llz - 1)

      parameter (lcx1 = lcx - 1)
      parameter (lcy1 = lcy - 1)
      parameter (lcz1 = lcz - 1)


      parameter (nc1 = nc - 1)
      parameter (ncp1 = ncp - 1)

c Mathematical and physical constants

      REAL kb

      parameter( kb = 8.61734315d-05)  ! eV/K
