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
    real(c_double), intent(inout):: synforc01(xdim*ydim, 2), randfld01(xdim*ydim, 10)   ! variables need to be saved in memory
    !real(c_double),public, intent(out):: synforc00(xdim*ydim, 2)
    !real(c_double) :: synforc00(xdim*ydim, 2), randfld00(xdim*ydim, 10)
    !
    integer ix, jy, id
    do jy=1,ydim
      do ix=1,xdim
        id = ix + (jy-1)*xdim
        print*,ix,jy,id,synforc01(id,1),synforc01(id,2)
        synforc01(id,1)=-1.*synforc01(id,1)
        synforc01(id,2)=-1.*synforc01(id,2)
       enddo
    enddo    
    call limits_randf(xdim,ydim)  ! read in setting from pseudo2D.nml
    call init_fvars    ! init field variables
    return
    call init_rand_update(synforc01,randfld01) ! core routine, xdim, ydim are set for idm,jdm in the routines
    ! calcualte synforc01, randfld01 in the routine
    ! output synforc00, synforc01, randfld01.
    ! check if it is necessary to return synforc00
  end subroutine p_pseudo2D_fld_sub

end module p_pseudo2D_fld
