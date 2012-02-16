      DOUBLE COMPLEX function pade_prime(z, n_pade, p, q)

#include "constants.F"

      DOUBLE COMPLEX z
      INTEGER n_pade 
      DOUBLE COMPLEX p(0:n_pade_max), q(0:n_pade_max)

      DOUBLE COMPLEX Pz, Qz, Pzp, Qzp 
      INTEGER i

      Pz  = dcmplx(0.0d0,0.0d0)
      Qz  = z**(n_pade + 1)

      Pzp = dcmplx(0.0d0, 0.0d0)
      Qzp = float(n_pade+1)*z**n_pade

      do i = 0, n_pade
        Pz = Pz + p(i) * z**i
        Pzp = Pzp + float(i) * p(i) * z**(i-1)
        Qz = Qz + q(i) * z**i
        Qzp = Qzp + float(i) * q(i) * z**(i-1)
      enddo

      pade_prime = (Qz*Pzp - Pz*Qzp) / Qz**2 

      return 
      end
