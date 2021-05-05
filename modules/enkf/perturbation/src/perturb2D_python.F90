subroutine generate_perturbation(xdim, ydim, synforc,randfld, previous_perturbation_exist) ! modify synforc,randfld and return them
use p_pseudo2d_fld
implicit none 
!
integer, intent(in):: xdim, ydim, previous_perturbation_exist 
! previous_perturbation_exist indicates whether previous perturbations exist (=1 exist, =0 inexist)
real*8, intent(inout) :: synforc(xdim*ydim,4), randfld(xdim*ydim,4) !

call p_pseudo2D_fld_sub(xdim, ydim, synforc,randfld, previous_perturbation_exist)
end subroutine generate_perturbation