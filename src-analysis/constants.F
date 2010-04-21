      implicit none

c Lattice and frequency grid parameters. (all must be 2^n!)
      integer nb  
      integer lcx, lcy, lcz, nc 
      integer llx, lly, llz, nl
      integer n_w_max, n_pade_max

      parameter (nb = 3)

      parameter (lcx = 4)
      parameter (lcy = 4)
      parameter (lcz = 1)
      parameter (nc = lcx*lcy*lcz)  ! number of cluster sites

      parameter (llx = 32)     ! lattice dimension >=lcx 
      parameter (lly = 32)
      parameter (llz = 1)
      parameter (nl = llx*lly*llz)

      parameter (n_w_max = 1024)
      parameter (n_pade_max = 10)

c Convenience constants

      integer llx1, lly1, llz1
      integer nc1
      integer lcx1, lcy1, lcz1

      parameter (llx1 = llx - 1)
      parameter (lly1 = lly - 1)
      parameter (llz1 = llz - 1)

      parameter (lcx1 = lcx - 1)
      parameter (lcy1 = lcy - 1)
      parameter (lcz1 = lcz - 1)


      parameter (nc1 = nc - 1)

c Mathematical and physical constants

      DOUBLE PRECISION kb
      parameter( kb = 8.61734315d-05)  ! eV/K
