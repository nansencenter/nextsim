module p_pseudo2D_fld

use, intrinsic :: ISO_C_BINDING


contains
      subroutine p_pseudo2D_fld_sub() bind(C,NAME='p_pseudo2D_fld_sub')

        ! Generates a random 2D field on a 100x100 grid. 
        use mod_random_forcing
        use mod_pseudo
        use m_set_random_seed
        implicit none 


        character(len=80) outfile

        real rv,rh, alp,bet
        integer i,j,n1,n2

        !-- CHeCK: idm, jdm to be read from namelist file --!
        !-- CHeCK: remove the computations in this routine \
        !-- they should  be called from mod_random_forcing --!
        integer, parameter :: idm=360
        integer, parameter :: jdm=360

        integer k,l,iostat
        real, allocatable, dimension(:,:) :: fld2d,ranfld,accranfld

        ! Initialize
        allocate(fld2d    (idm,jdm))
        allocate(ranfld   (idm,jdm))
        allocate(accranfld(idm,jdm))

        rv=2.0e00  ! vertical correlation range (in nb of layers)
        rh=200.00/3.0   ! horizontal decorrelation scale (number of grid cells)

        !  Generates the vertical correlation of the ensemble.
        alp=exp(-1.0/rv)
        bet=sqrt(1.0-alp**2) ! keeps var(enstmp)=var(ensmem)

        ! FFT dimensions are powers of 2 for best performance 
        n1=2**(ceiling(log(float(idm))/log(2.)))
        n2=2**(ceiling(log(float(jdm))/log(2.)))
        ! set file names
        outfile='pseudo2D_fields.dat'


        call init_fvars
        call init_rand_update
        call set_random_seed

        open (unit=801,file=trim(outfile),status='replace')


        call rand_update

        do i=1,idm 
        do j=1,jdm 
           write(801,'(2i6,f6.2)') i,j,ranfld(i,j)
        end do  
        end do  

        close(801)

        deallocate(accranfld,fld2d,ranfld)
      end subroutine p_pseudo2D_fld_sub

end module p_pseudo2D_fld

