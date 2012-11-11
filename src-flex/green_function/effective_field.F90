#include "../convert.F90"

subroutine effective_field(iteration, v_pert, h, h_pert, prfld, & 
     prfld_pert, v_pert_eff, h_eff, prfld_eff)

  USE CONSTANTS
  IMPLICIT NONE

  INTEGER iteration
  REAL v_pert(0:nb-1)
  REAL, dimension (0:nb-1,1:3) :: h, h_pert
  REAL prfld, prfld_pert

  REAL, dimension (0:nb-1,1:3) :: h_eff
  REAL prfld_eff
  REAL v_pert_eff(0:nb-1)

  if (iteration .le. 10) then

     h_eff = h + h_pert * float(10-iteration) / 10.0d0
     prfld_eff = prfld + prfld_pert * float(10-iteration) / 10.0d0
     v_pert_eff = v_pert * float(10-iteration) / 10.0d0

  else 

     h_eff = h
     prfld_eff = prfld
     v_pert_eff = 0.0d0
        
  endif
  
  return
end subroutine effective_field
