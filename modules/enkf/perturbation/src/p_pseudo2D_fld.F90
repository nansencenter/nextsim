module p_pseudo2D_fld

use, intrinsic :: ISO_C_BINDING


contains
      subroutine p_pseudo2D_fld_sub() bind(C,NAME='p_pseudo2D_fld_sub')

        ! Generates a random 2D field on a 100x100 grid. 
        use mod_random_forcing
        use mod_pseudo
        use m_set_random_seed
        implicit none 


        integer i,j,n1,n2

        real current_time

        !-- CHeCK: idm, jdm to be read from namelist file --!
        !-- CHeCK: remove the computations in this routine \
        !-- they should  be called from mod_random_forcing --!
        integer, parameter :: idm=360
        integer, parameter :: jdm=360

        integer k,l,iostat

        call init_fvars
        call init_rand_update

      end subroutine p_pseudo2D_fld_sub

end module p_pseudo2D_fld

