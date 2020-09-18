module p_pseudo2D_fld

use, intrinsic :: ISO_C_BINDING

contains
  subroutine p_pseudo2D_fld_sub(xdim, ydim, synforc01, randfld01) bind(C,NAME='p_pseudo2D_fld_sub')

    ! Generates a random 2D field on a 100x100 grid. 
    use mod_random_forcing
    use mod_pseudo
    use m_set_random_seed
    implicit none 
    !
    integer(c_int), intent(in):: xdim, ydim
    real(c_float),public, intent(inout):: synforc01(2,xdim*ydim), randfld01(10,xdim*ydim)   ! variables need to be saved in memory
    !real(c_float),public, intent(out):: synforc00(2,xdim*ydim)
    real(c_float),public:: synforc00(2,xdim*ydim), randfld00(10,xdim*ydim)
    !
    
    call limits_randf  ! read in setting from pseudo2D.nml
    call init_fvars    ! init field variables
    call init_rand_update(xdim,ydim) ! core routine, xdim, ydim are set for idm,jdm in the routines
    ! todo: give values of synforc01, randfld01 to variables by synforc_rd, randfld_rd
    ! 
    ! synforc00=synforc01, randfld00= randfld01 
    ! calcualte synforc01, randfld01 in the routine
    ! output synforc00, synforc01, randfld01.
  end subroutine p_pseudo2D_fld_sub

end module p_pseudo2D_fld
