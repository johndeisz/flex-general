      implicit none

c Lattice and frequency grid parameters. (all must be 2^n!)
      integer nb  
      integer lcx, lcy, lcz, nc 
      integer llx, lly, llz, nl
      integer n_pade_max

      parameter (nb = 3)

      parameter (lcx = 4)
      parameter (lcy = 4)
      parameter (lcz = 1)
      parameter (nc = lcx*lcy*lcz)  ! number of cluster sites

      parameter (llx = 64)     ! lattice dimension >=lcx 
      parameter (lly = 64)
      parameter (llz = 1)
      parameter (nl = llx*lly*llz)

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

      DOUBLE PRECISION gs, kb, pi
      parameter( kb = 8.61734315d-05)  ! eV/K
      parameter( pi = 3.14159265358979323)
      parameter( gs = 2.002319)

