module p_pseudo2D_fld

use, intrinsic :: ISO_C_BINDING

contains
  subroutine p_pseudo2D_fld_sub(xdim, ydim,synforc01,randfld01, perturbation_count) bind(C,NAME='p_pseudo2D_fld_sub')

    ! Generates a random 2D field on a 100x100 grid. 
    use mod_random_forcing
    use mod_pseudo
    use m_set_random_seed
    implicit none 
    !
    integer(c_int), intent(in):: xdim, ydim, perturbation_count ! perturbation_count is only used to decide whether generate ran
    !real(c_double), intent(inout):: synforc01(xdim*ydim, 2), randfld01(xdim*ydim, 10)   ! variables need to be saved in memory
    real*8 :: synforc01(xdim*ydim, 2), randfld01(xdim*ydim, 10) !
  
    call limits_randf(xdim,ydim)  ! read in setting from pseudo2D.nml
    call init_fvars               ! init field variables
    call init_rand_update(synforc01,randfld01, perturbation_count) ! core routine, xdim, ydim are set for idm,jdm in the routines

 ! calcualte synforc01, randfld01 in the routine
    ! output synforc00, synforc01, randfld01.
  end subroutine p_pseudo2D_fld_sub

end module p_pseudo2D_fld
