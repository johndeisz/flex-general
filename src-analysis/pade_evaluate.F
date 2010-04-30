      DOUBLE COMPLEX function pade_evaluate(z, n_pade, p, q)

#include "constants.F"

      DOUBLE COMPLEX z
      INTEGER n_pade 
      DOUBLE COMPLEX p(0:n_pade_max), q(0:n_pade_max)

      DOUBLE COMPLEX num, den, result
      INTEGER i

      num  = dcmplx(0.0d0,0.0d0)
      den  = z**(n_pade + 1)

      do i = 0, n_pade
        num = num + p(i) * z**i
        den = den + q(i) * z**i
      enddo

      pade_evaluate = num / den

      return 
      end
