      implicit none

c Lattice and frequency grid parameters. (all must be 2^n!)
      integer nb, m  
      integer lcx, lcy, lcz, nc 
      integer llx, lly, llz, nl

      parameter (nb = 3)

      parameter (m =  256) 
      parameter (lcx = 8)
      parameter (lcy = 8)
      parameter (lcz = 1)
      parameter (nc = lcx*lcy*lcz)  ! number of cluster sites

      parameter (llx = 64)     ! lattice dimension >=lcx 
      parameter (lly = 64)
      parameter (llz = 1)
      parameter (nl = llx*lly*llz)

c MPI related constants.

      integer np, mp, ncp, max_message

      parameter (np = 1)
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

      REAL kb, mub, gs

      parameter( kb = 8.61734315d-05)  ! eV/K
      parameter( mub = 5.788381755d-5) ! eV/T
      parameter( gs = 2.002319)
