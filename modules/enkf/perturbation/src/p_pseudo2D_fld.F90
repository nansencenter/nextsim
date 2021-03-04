module p_pseudo2D_fld

use, intrinsic :: ISO_C_BINDING

contains
  subroutine p_pseudo2D_fld_sub(xdim, ydim, synforc,randfld, previous_perturbation_exist) bind(C,NAME='p_pseudo2D_fld_sub') ! modify synforc,randfld and return them

    ! Generates a random 2D field on a 100x100 grid. 
    use mod_random_forcing
    use mod_pseudo
    use m_set_random_seed
    implicit none 
    !
    integer(c_int), intent(in):: xdim, ydim, previous_perturbation_exist 
    ! previous_perturbation_exist indicates whether previous perturbations exist (=1 exist, =0 inexist)
    real(c_double) :: synforc(xdim*ydim,4), randfld(xdim*ydim,4) !
  
    call limits_randf(xdim,ydim)  ! read in setting from pseudo2D.nml
    call init_fvars               ! init field variables
    call init_rand_update(synforc,randfld,previous_perturbation_exist) ! core subroutine, xdim, ydim are set for idm,jdm in the routines
  end subroutine p_pseudo2D_fld_sub

end module p_pseudo2D_fld
