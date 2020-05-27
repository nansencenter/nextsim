module p_pseudo2D_fld

use, intrinsic :: ISO_C_BINDING


contains
      subroutine p_pseudo2D_fld_sub() bind(C,NAME='p_pseudo2D_fld_sub')

        ! Generates a random 2D field on a 100x100 grid. 
        use mod_random_forcing
        use mod_pseudo
        use m_set_random_seed
        implicit none 
        
        call limits_randf  ! read in setting from pseudo2D.nml
        call init_fvars    ! init field variables
        call init_rand_update ! core routine and output perturbed field
      end subroutine p_pseudo2D_fld_sub

end module p_pseudo2D_fld
