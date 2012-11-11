MODULE pade_eval

#include "../convert.F90"

CONTAINS
  COMPLEX function pade_evaluate(z, n_pade, p, q)

    USE CONSTANTS
    IMPLICIT NONE

    COMPLEX z
    INTEGER n_pade 
    COMPLEX p(0:m/2-1), q(0:m/2-1)
    
    COMPLEX num, den, result
    INTEGER i

    num  = cmplx(0.0d0,0.0d0)
    den  = z**(n_pade + 1)

    do i = 0, n_pade
       num = num + p(i) * z**i
       den = den + q(i) * z**i
    enddo

    pade_evaluate = num / den
    
    return 
  end function pade_evaluate

END MODULE pade_eval
