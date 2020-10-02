module mod_random_forcing
! ----------------------------------------------------------------------------
! ToDO:
! - generate grid from ASR and ECMWF forcing files
! - compute scpx and scpy from spheredist8 in m_grid_to_hycom.F
! ----------------------------------------------------------------------------
! ----------------------------------------------------------------------------
! -- init_rand_update - Initializes the random module, reads infile2.in - sets
!                       array sizes of derived types and initializes variables
!                       ran1 and ran. Also sets random number seeed, and sets
!                       the FFT sizes in mod_pseudo
! -- rand_update      - main routine - updates the nondimensional perturbations
!                       contained in variable "ran", and applies the
!                       dimensionalized version of this to the different forcing
!                       fields contained in "forcing_fields". Note that this
!                       routine is (for now) applied only when hf forcing is
!                       enabled (a check is done against forcing update times
!                       "rdtime" given by mod_forcing_nersc and yrflag). Also,
!                       synoptic flags from "mod_forcing_nersc" are not taken
!                       into account.
!
! + Various private routines:
! -- limits_randf     - reads infile2.in and sets forcing variances as well as
!                       spatial and temporal correlation lengths (private)
! -- set_random_seed2 - Sets the random number seed based on date and time
! -- ranfields          Produces a set of nondimensional perturbation fields in
!                       a forcing_fields type
! -- calc_forc_update   Creates a dimensional forcing_field type
! -- assign_force       sets  forcing_field to a constant
! -- assign_vars        Sets variances to a constant
! -- ran_update_ran1    Updates nondimensional forcing_field "ran" using
!                       input variances and nondimensional stochastic forcing
!                       "ran1" - produces a red timeseries specified by input
!                       "alpha"
! -- init_ran(ran)      Allocates variables in the forcing_fields type
! ----------------------------------------------------------------------------
! KAL - 1st offline version 20080830
! NB - to use this routine mod_grid must be initialized (get_grid), as well
!      as the mod_forcing module (init_forcing_nersc). This is needed to
!      allocate grid and forcing variables
! ----------------------------------------------------------------------------

   implicit none
   private

   logical, save :: randf ! Switches on/off random forcing
   logical, save :: debug ! Switches on/off verbose
   logical       :: file_exists ! Check if file exists
   real   , save :: rf_hradius  ! Horizontal decorr length for rand forc [m]
   real   , save :: rf_tradius  ! Temporal decorr length for rand forc
   real   , save :: rh          ! Horizontal decorr length for rand forc [grid cells]
   real   , save :: rv
   integer, save :: rf_prsflg=2 ! initial value
   integer, save :: xdim, ydim
   integer, save :: idm=360, jdm=360
   character(150) :: iopath = '/docker_io' ! modified in namelist
   real,parameter     :: airdns  =  1.2
   real, parameter    :: radian  = 57.2957795
   real, parameter    :: pi      =  3.1415926536
   real, parameter    :: radtodeg= 57.2957795
   real, allocatable, dimension(:,:)  :: synuwind, synvwind, synwndspd,         &
                                         synairtmp, synrelhum, synprecip,       &
                                         synclouds, syntaux, syntauy, synvapmix,&
                                         synradflx, synshwflx, synslp, synssr
   
  ! Random forcing variables:
   type forcing_fields
      real,pointer ::  slp    (:,:) !  Sea level pressure
      real,pointer ::  taux   (:,:) !  wind stress in x direction
      real,pointer ::  tauy   (:,:) !  wind stress in y direction
      real,pointer ::  wndspd (:,:) !  wind speed (tke source)
      real,pointer ::  airtmp (:,:) !  pseudo air temperature
      real,pointer ::  relhum (:,:) !  relative humidity
      real,pointer ::  clouds (:,:) !  cloud cover
      real,pointer ::  precip (:,:) !  precipitation
      real,pointer ::  sss    (:,:) !  SSS for relax
      real,pointer ::  sst    (:,:) !  SST for relax

      real,pointer ::  uwind  (:,:) !  u-component of wind
      real,pointer ::  vwind  (:,:) !  v-component of wind
      real,pointer ::  tauxice(:,:) !  ice stress on water in x dir
      real,pointer ::  tauyice(:,:) !  ice stress on water in y dir
   end type forcing_fields

   type forcing_variances
      real slp
      real taux
      real tauy
      real wndspd
      real airtmp
      real relhum
      real clouds
      real precip
      real sss
      real sst
   end type forcing_variances


   ! Variable containing forcing variances
   type(forcing_variances), save :: vars

   ! Variable containing random forcing
   type(forcing_fields)   , save :: ran , ran1


   interface assignment(=)
      module procedure assign_force
      module procedure assign_vars
   end interface

   interface sqrt
      module procedure var_sqrt
   end interface

   public :: randf, init_rand_update, rand_update,init_fvars,limits_randf

contains
      
   subroutine init_rand_update(synforc01,randfld01,perturbation_count)
      use mod_pseudo
      implicit none
      integer:: perturbation_count
      real :: dx
      real*8, dimension(idm*jdm, 2) :: synforc00, synforc01
      real*8, dimension(idm*jdm,10) :: randfld00, randfld01
      if(.not.randf) then
        write(*,'("randf option switched off in pseudo2D.nml,no perturbation will be applied")')
        return
      end if
      write(*,'("pseudo-random forcing is active for ensemble generation")')
      dx=30  !scpx(idm/2,jdm/2)
      if(debug) print*,'typical model grid scale ', dx
      rh=rf_hradius/dx     ! Decorrelation length is rh grid cells
      rv=rf_tradius        ! Temporal decorrelation scale (days)

      if(debug) print*, "initialized init_ran"
      call init_ran(ran)
      call init_ran(ran1)
      ran=0.
      ran1=0.
      ! Init
      call set_random_seed2
      ! Init fft dimensions in mod_pseudo
      call initfftdim(idm,jdm)

      !-- CHeCK: a conditional here to check ranfld_next exists \
      !-- IF exists, load and move to ranfld_prev, IF NOT run ranfields(ran,rh) --!
      if (perturbation_count>0) then 
         if (debug) print*,  'set perturbations as previous one'
         call load_randfld_synforc(synforc01,randfld01)
         randfld00 = randfld01
         !synforc00 = synforc01   ! todo: ensure synforc00 is not used to addPerturbation 
      else
         if (debug) print*, 'generating initial random field...'
         call ranfields(ran,rh)
         call rand_update('00',randfld00, synforc00)     
      end if
      call rand_update('01',randfld01, synforc01)
   end subroutine

!----------------------------------
   subroutine set_random_seed2
      ! Sets a random seed based on the wall clock time
      implicit none

      integer , dimension(8)::val
      integer cnt
      integer sze
      integer, allocatable, dimension(:):: pt  ! keeps random seed

      call DATE_AND_TIME(values=val)
      call RANDOM_SEED(size=sze)
      allocate( pt(sze))
      call RANDOM_SEED         ! Init - assumes seed is set in some way based on clock, date etc. (not specified in f ortran standard)
      call RANDOM_SEED(GET=pt) ! Retrieve initialized seed
      pt = pt * (val(8)-500)  ! val(8) is milliseconds - this randomizes stuff if random_seed is nut updated often e nough
      call RANDOM_SEED(put=pt)
      deallocate( pt)
   end subroutine set_random_seed2

!c --- Initialize FFT dimensions used in pseudo routines
   subroutine initfftdim(nx,ny)
      use mod_pseudo
      implicit none
      integer, intent(in) :: nx,ny
      fnx = ceiling(log(float(nx))/log(2.))
      fnx = 2**fnx
      fny = ceiling(log(float(ny))/log(2.))
      fny = 2**fny
      if (debug) write(*,'("Fourier transform dimensions ",2i6)')  fnx,fny
   end subroutine

   subroutine limits_randf(xdim, ydim)
      implicit none
      real, parameter :: version2=1.2     ! version of limits routine

      logical :: ex
      integer :: seed, prsflg
      integer, intent(in):: xdim,ydim
      real    :: fversion
      real    :: vslp, vtaux, vtauy, vwndspd, vclouds
      real    :: vairtmp, vprecip, vrelhum, scorr, tcorr
      character(80) :: nmlfile, cwd

      namelist /setup/ iopath, debug
      namelist /pseudo2D/ randf, seed, &
                          vslp, vtaux, vtauy, vwndspd, &
                          vclouds, vairtmp, vprecip, vrelhum, &
                          scorr, tcorr, prsflg

!      randf        = .true.
!      seed         = 11
!      vars%slp     =  10.0
!      vars%taux    =  1.e-3
!      vars%tauy    =  1.e-3
!      vars%wndspd  =  0.64
!      vars%clouds  =  5.e-3
!      vars%airtmp  =  9.0
!      vars%precip  =  1.0
!      vars%relhum  =  1.0
!      rf_hradius   =  500
!      rf_tradius   =  2.0
!      rf_prsflg    =  2

      if (rf_prsflg>2 .or. rf_prsflg<0) then
            print*, 'Pressure flag must be between 0 and 2'
         stop '(limits:limits_randf)'
      end if


      nmlfile ='pseudo2D.nml'    ! name of general configuration namelist file
      open (99,file=nmlfile, status='old', action='read')
      read (99,NML=setup)
      read (99,NML=pseudo2D)
      close(99)

      idm         = xdim   ! remove xdim, ydim from pseudo2D.nml
      jdm         = ydim   
      vars%slp    = vslp
      vars%taux   = vtaux
      vars%tauy   = vtauy
      vars%wndspd = vwndspd
      vars%clouds = vclouds
      vars%airtmp = vairtmp
      vars%precip = vprecip
      vars%relhum = vrelhum
      rf_hradius  = scorr
      rf_tradius  = tcorr
      rf_prsflg   = prsflg

   end subroutine limits_randf


! --- This routine updates the random forcing component, according to
! --- a simple correlation progression with target variance specified
! --- By forcing_variances. At the end of the routine, if applicable,
! --- the random forcing is added to the forcing fields.
   subroutine rand_update(time_index, randfld, synforc)
      implicit none

      ! rt       -- Information on time (mod_year_info)
      ! ran      -- Nondimensional random perturbations
      ! vars     -- Variances of fields ( Real pert = ran * vars)
      ! lrestart -- Special actions are taken if this is a restart
      !type(forcing_fields)    , intent(inout) :: ran
      !type(forcing_variances) , intent(in)    :: vars

      integer :: ix,jy
      real :: alpha, autocorr, nsteps, wspd, mtime
      character(2) :: time_index
      logical, save :: first=.true.

      real, dimension(idm,jdm) :: modlon, modlat

      logical, parameter :: randf_test=.true.
      logical, save :: lfirst=.true.
      real, parameter :: rhoa = 1.2 , cdfac = 0.0012

      real :: cd_new,w4,wfact,wndfac, fcor
      real, dimension(1:idm,1:jdm) :: dpresx,dpresy
      real :: ucor, vcor, ueq,veq,wcor
      real :: wprsfac, minscpx

      real, parameter :: wlat=60.
      integer i,j
      real*8, save :: rdtime=8.d0/24.d0    ! Time step of forcing update
      real*8  :: randfld(idm*jdm,10), synforc(idm*jdm,2)

      ! Autocorrelation between two times "tcorr"
      !KAL - quite high? - autocorr = 0.95
      autocorr = exp(-1.0)

      ! Number of "forcing" steps rdtime between autocorrelation
      ! decay. Autocorrelation time is tcorr
      nsteps = rv/rdtime

      ! This alpha will guarantee autocorrelation "autocorr"
      ! over the time "rv"
      ! rv -> infinity , nsteps -> infinity, alpha -> 1
      ! rv -> 0        , nsteps -> 0       , alpha -> 0 (when 1>autocorr>0)
      alpha=autocorr**(1/nsteps)

      !write(lp,*) 'Rand_update -- Random forcing field update'
      ! Add new random forcing field to the newly read
      ! fields from ecmwf or ncep (:,:,4)
      !ran1=sqrt(vars)*ran

      call calc_forc_update(ran1,ran,sqrt(vars)) ! note: 3ran is only used in subroutine ran_update_ran1

      if (rf_prsflg .eq. 1 .or. rf_prsflg .eq.2 ) then
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! rf_prsflag=1 : wind perturbations calculated from slp, using coriolis
   !                parameter at 40 deg N
   ! rf_prsflag=2 : wind perturbations calculated from slp, using coriolis
   !                parameter at 40 deg N, but limited by the setting of
   !                windspeed, to account for the horizontal scale of pert.
   !                As far as the wind is concerned, this is the same as
   !                reducing pressure perturbations
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! grid size min max
   !      minscpx=minval(scpx)
         minscpx=30
         wprsfac=1.
         ! flag used in prsflg=2
         if (rf_prsflg==2) then
            fcor=2*sin(40./radtodeg)*2*pi/86400; ! Constant

            ! typical pressure gradient
            wprsfac=100.*sqrt(vars%slp)/(rh*minscpx)

            ! results in this typical wind magnitude
            wprsfac=wprsfac/fcor

            ! but should give wind according to vars%wndspd
            ! this is a correction factor for that
            wprsfac=sqrt(vars%wndspd)/(3*wprsfac)
         end if

         dpresx=0.
         dpresy=0.
         do jy=2,jdm
         do ix=2,idm
         ! Pressure gradient. Coversion from mBar to Pa
            dpresx(ix,jy) = &
               100.*(ran1%slp(ix,jy) - ran1%slp(ix-1,jy))/ minscpx
            dpresx(ix,jy)=dpresx(ix,jy)*wprsfac

            dpresy(ix,jy) = &
               100.*(ran1%slp(ix,jy) - ran1%slp(ix,jy-1))/ minscpx
            dpresy(ix,jy)=dpresy(ix,jy)*wprsfac
         end do
         end do


         do jy=1,jdm
         do ix=1,idm
            
            ! Coriolis balance (at 40 deg)
            !fcor=2*sin(max(abs(plat(ix,jy)),20.)/radtodeg)*2*pi/86400;
            fcor=2*sin(40./radtodeg)*2*pi/86400; ! Constant
            fcor=fcor*rhoa
            vcor= dpresx(ix,jy) / (fcor)
            ucor=-dpresy(ix,jy) / (fcor)

            ! In the equatorial band u,v are aligned with the
            ! pressure gradients. Here we use the coriolis
            ! factor above to set it up (to limit the speeds)
            ueq=-dpresx(ix,jy) / abs(fcor)
            veq=-dpresy(ix,jy) / abs(fcor)

            ! Weighting between coriiolis/equator solution
   !         wcor=sin( &
   !            min(abs(plat(ix,jy)),wlat) / wlat * pi * 0.5)
            wcor=sin(wlat) / wlat * pi * 0.5

            synuwind(ix,jy) = wcor*ucor + (1.-wcor)*ueq
            synvwind(ix,jy) = wcor*vcor + (1.-wcor)*veq

            synwndspd(ix,jy) = sqrt(  &
               synuwind(ix,jy)**2 + synvwind(ix,jy)**2)

            ! The rest use uncorrelated fields
            synairtmp(ix,jy) = ran1%airtmp(ix,jy)
            synrelhum(ix,jy) = ran1%relhum(ix,jy)
            synslp   (ix,jy) = ran1%slp   (ix,jy)
            synprecip(ix,jy) = ran1%precip(ix,jy) 
            synrelhum(ix,jy) = min(max(synrelhum(ix,jy),0.0),1.0)
            synprecip(ix,jy) = max(synprecip(ix,jy),0.0)
            synwndspd(ix,jy) = max(synwndspd(ix,jy),0.0) 
         end do
         end do
   ! todo
   !      synuwind = synuwind - compute_mean(synuwind)
   !      synvwind = synvwind - compute_mean(synvwind)

         ! Drag
         do jy=2,jdm-1
         do ix=2,idm-1
            wndfac=(1.+sign(1.,synwndspd(ix,jy)-11.))*.5
            cd_new=(0.49+0.065*synwndspd(ix,jy))*1.0e-3*wndfac+cdfac*(1.-wndfac)

            w4    =.25*( &
               synvwind(ix-1,jy+1)+synvwind(ix,jy+1)+  &
               synvwind(ix-1,jy  )+synvwind(ix,jy  ))
            wfact=sqrt( synuwind(ix,jy)*synuwind(ix,jy)+w4*w4)* airdns*cd_new
            syntaux(ix,jy)=synuwind(ix,jy)*wfact

            w4   =.25*( &
               synuwind(ix  ,jy-1)+synuwind(ix+1,jy-1)+ &
               synuwind(ix+1,jy  )+synuwind(ix  ,jy  ))
            wfact=sqrt( &
                  synvwind(ix,jy)*synvwind(ix,jy)+w4*w4)* airdns*cd_new
            syntauy(ix,jy)=synvwind(ix,jy)*wfact
         end do
         end do

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! rf_prsflag=0 : wind and slp are uncorrelated
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else  ! rf_prsflg .eq. 0
         do jy=2,jdm-1
         do ix=2,idm-1
         !if (ip(ix,jy)==1) then
            syntaux  (ix,jy) = syntaux  (ix,jy) + ran1%taux(ix,jy)
            syntauy  (ix,jy) = syntauy  (ix,jy) + ran1%tauy(ix,jy)

            ! KAL -- winds are nonlinear functions of tau and mainly
            ! KAL -- used for sea ice
            wspd = sqrt(syntaux(ix,jy)**2 + syntauy(ix,jy)**2)
            wspd = max(sqrt(wspd / (cdfac*rhoa)),0.1)
            synuwind(ix,jy) = syntaux (ix,jy) / (wspd*cdfac*rhoa)
            synvwind(ix,jy) = syntauy (ix,jy) / (wspd*cdfac*rhoa)
            synairtmp(ix,jy) = synairtmp(ix,jy)+ran1%airtmp(ix,jy)
            synwndspd(ix,jy) = synwndspd(ix,jy)+ran1%wndspd(ix,jy)
            synrelhum(ix,jy) = synrelhum(ix,jy)+ran1%relhum(ix,jy)
            synprecip(ix,jy) = synprecip(ix,jy)+ran1%precip(ix,jy)
            synrelhum(ix,jy) = min(max(synrelhum(ix,jy),0.0),1.0)
            synprecip(ix,jy) = max(synprecip(ix,jy),0.0)
            synwndspd(ix,jy) = max(synwndspd(ix,jy),0.0)
         !end if
         end do
         end do
      end if ! rf_prsflg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! end if rf_prsflag=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write output files  -- Spatial field dumped on first run （Instead, spatial fields are saved in variables）
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      
      ! ran1 is new random forcing. ran is nondimensional
      ! "Brownian increment".
      call ranfields(ran1,rh)
      !ran= alpha*ran + sqrt(1-alpha*alpha)* ran
      call ran_update_ran1(ran,ran1,alpha)
      call save_randfld_synforc(randfld, synforc) ! save fields to variables, final file output is moved to function exportWindPerturbation
      !call synforc_wr(time_index)
      !call randfld_wr(time_index)
      
   end subroutine rand_update
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine ranfields(ranfld,scorr)
      use mod_pseudo
      implicit none

      type(forcing_fields)    , intent(inout) :: ranfld
      real,                     intent(in)    :: scorr
      real, dimension(idm,jdm) :: gtmp, tmp

      ranfld=0.
      call pseudo2D(ranfld%slp,idm,jdm,1,scorr,fnx,fny)
      call pseudo2D(ranfld%taux,idm,jdm,1,scorr,fnx,fny)
      call pseudo2D(ranfld%tauy,idm,jdm,1,scorr,fnx,fny)
      call pseudo2D(ranfld%wndspd,idm,jdm,1,scorr,fnx,fny)
      call pseudo2D(ranfld%airtmp,idm,jdm,1,scorr,fnx,fny)
      call pseudo2D(ranfld%relhum,idm,jdm,1,scorr,fnx,fny)
      call pseudo2D(ranfld%clouds,idm,jdm,1,scorr,fnx,fny)
      call pseudo2D(ranfld%precip,idm,jdm,1,scorr,fnx,fny)
   end subroutine ranfields

!------------------------------------------
   subroutine calc_forc_update(A,B,C)
      type(forcing_fields), intent(inout) :: A
      type(forcing_fields), intent(in   ) :: B
      type(forcing_variances), intent(in) :: C
      integer :: i,j
      do j=1,jdm
      do i=1,idm
         A%slp   (i,j)=C%slp    * B%slp   (i,j)
         A%taux  (i,j)=C%taux   * B%taux  (i,j)
         A%tauy  (i,j)=C%tauy   * B%tauy  (i,j)
         A%wndspd(i,j)=C%wndspd * B%wndspd(i,j)
         A%airtmp(i,j)=C%airtmp * B%airtmp(i,j)
         A%relhum(i,j)=C%relhum * B%relhum(i,j)
         A%clouds(i,j)=C%clouds * B%clouds(i,j)
         A%precip(i,j)=C%precip * B%precip(i,j)
         A%sss   (i,j)=C%sss    * B%sss   (i,j)
         A%sst   (i,j)=C%sst    * B%sst   (i,j)
      end do
      end do
   end subroutine calc_forc_update

!--------------------------------------------------
   subroutine ran_update_ran1(ran,ran1,alpha)
      type(forcing_fields), intent(inout) :: ran
      type(forcing_fields), intent(   in) :: ran1
      real                , intent(   in) :: alpha

      integer :: ix,jy

      do jy=1,jdm
      do ix=1,idm
         ran%slp   (ix,jy)=alpha*ran%slp   (ix,jy) + sqrt(1-alpha*alpha)*ran1%slp   (ix,jy)
         ran%taux  (ix,jy)=alpha*ran%taux  (ix,jy) + sqrt(1-alpha*alpha)*ran1%taux  (ix,jy)
         ran%tauy  (ix,jy)=alpha*ran%tauy  (ix,jy) + sqrt(1-alpha*alpha)*ran1%tauy  (ix,jy)
         ran%wndspd(ix,jy)=alpha*ran%wndspd(ix,jy) + sqrt(1-alpha*alpha)*ran1%wndspd(ix,jy)
         ran%airtmp(ix,jy)=alpha*ran%airtmp(ix,jy) + sqrt(1-alpha*alpha)*ran1%airtmp(ix,jy)
         ran%relhum(ix,jy)=alpha*ran%relhum(ix,jy) + sqrt(1-alpha*alpha)*ran1%relhum(ix,jy)
         ran%clouds(ix,jy)=alpha*ran%clouds(ix,jy) + sqrt(1-alpha*alpha)*ran1%clouds(ix,jy)
         ran%precip(ix,jy)=alpha*ran%precip(ix,jy) + sqrt(1-alpha*alpha)*ran1%precip(ix,jy)
         ran%sss   (ix,jy)=alpha*ran%sss   (ix,jy) + sqrt(1-alpha*alpha)*ran1%sss   (ix,jy)
         ran%sst   (ix,jy)=alpha*ran%sst   (ix,jy) + sqrt(1-alpha*alpha)*ran1%sst   (ix,jy)
      end do
      end do
   end subroutine

!------------------------------------------
   subroutine init_ran(ran)
   implicit none

      type(forcing_fields), intent(inout) :: ran

      allocate(ran%slp    (idm,jdm))
      allocate(ran%taux   (idm,jdm))
      allocate(ran%tauy   (idm,jdm))
      allocate(ran%wndspd (idm,jdm))
      allocate(ran%airtmp (idm,jdm))
      allocate(ran%relhum (idm,jdm))
      allocate(ran%clouds (idm,jdm))
      allocate(ran%precip (idm,jdm))
      allocate(ran%sss    (idm,jdm))
      allocate(ran%sst    (idm,jdm))
      allocate(ran%uwind  (idm,jdm))
      allocate(ran%vwind  (idm,jdm))
      allocate(ran%tauxice(idm,jdm))
      allocate(ran%tauyice(idm,jdm))
   end subroutine

!------------------------------------------
   subroutine init_fvars()
   implicit none
     ! Allocate fields - some are not used...
   
     IF( .NOT. ALLOCATED( synuwind  ) ) allocate(synuwind (idm,jdm))
     IF( .NOT. ALLOCATED( synvwind  ) ) allocate(synvwind (idm,jdm))
     IF( .NOT. ALLOCATED( synwndspd ) ) allocate(synwndspd(idm,jdm))
     IF( .NOT. ALLOCATED( syntaux   ) ) allocate(syntaux  (idm,jdm))
     IF( .NOT. ALLOCATED( syntauy   ) ) allocate(syntauy  (idm,jdm))
     IF( .NOT. ALLOCATED( synvapmix ) ) allocate(synvapmix(idm,jdm))
     IF( .NOT. ALLOCATED( synairtmp ) ) allocate(synairtmp(idm,jdm))
     IF( .NOT. ALLOCATED( synrelhum ) ) allocate(synrelhum(idm,jdm))
     IF( .NOT. ALLOCATED( synprecip ) ) allocate(synprecip(idm,jdm))
     IF( .NOT. ALLOCATED( synclouds ) ) allocate(synclouds(idm,jdm))
     IF( .NOT. ALLOCATED( synradflx ) ) allocate(synradflx(idm,jdm))
     IF( .NOT. ALLOCATED( synshwflx ) ) allocate(synshwflx(idm,jdm))
     IF( .NOT. ALLOCATED( synslp    ) ) allocate(synslp   (idm,jdm))
     synuwind (:,:)=0.
     synvwind (:,:)=0.
     synwndspd(:,:)=0.
     syntaux  (:,:)=0.
     syntauy  (:,:)=0.
     synvapmix(:,:)=0.
     synairtmp(:,:)=0.
     synrelhum(:,:)=0.
     synprecip(:,:)=0.
     synclouds(:,:)=0.
     synradflx(:,:)=0.
     synshwflx(:,:)=0.
     synslp   (:,:)=0.
   end subroutine

!------------------------------------------
   subroutine load_randfld_synforc(randfld, synforc) ! todo, check if it is necessary to read synforc
      integer :: ix,jy,id
      real*8  :: randfld(idm*jdm, 10), synforc(idm*jdm,2)
      do jy=1,jdm
      do ix=1,idm
         id = (jy-1)*idm + ix
         ran%slp(ix,jy)    = randfld(id,1)
         ran%taux(ix,jy)   = randfld(id,2)
         ran%tauy(ix,jy)   = randfld(id,3)
         ran%wndspd(ix,jy) = randfld(id,4)
         ran%airtmp(ix,jy) = randfld(id,5)
         ran%relhum(ix,jy) = randfld(id,6)
         ran%clouds(ix,jy) = randfld(id,7)
         ran%precip(ix,jy) = randfld(id,8)
         ran%sss(ix,jy)    = randfld(id,9)
         ran%sst(ix,jy)    = randfld(id,10)
      
         synuwind(ix,jy) = synforc(id,1)
         synvwind(ix,jy) = synforc(id,2)! since perturbtions is only applied to wind speeds, the unused variables are commented.
         !synairtmp(ix,jy), synslp(ix,jy), &
         !synprecip(ix,jy), synrelhum(ix,jy) /)
      end do 
      end do 
   end subroutine

   !------------------------------------------
   subroutine save_randfld_synforc(randfld, synforc)
      integer :: ix,jy,id
      real*8  :: randfld(idm*jdm,10), synforc(idm*jdm,2)
      do jy=1,jdm
      do ix=1,idm
         id = (jy-1)*idm + ix
         randfld(id,:) = (/ ran%slp(ix,jy),ran%taux(ix,jy),ran%tauy(ix,jy), &
         ran%wndspd(ix,jy),ran%airtmp(ix,jy),ran%relhum(ix,jy), &
         ran%clouds(ix,jy),ran%precip(ix,jy),ran%sss(ix,jy),ran%sst(ix,jy) /)              
         
         synforc(id,:) = (/ synuwind(ix,jy), synvwind(ix,jy) /) ! since perturbtions is only applied to wind speeds, the unused variables are commented.
         !synairtmp(ix,jy), synslp(ix,jy), &
         !synprecip(ix,jy), synrelhum(ix,jy) /)
         !print*,'save_r',ix,jy,synforc(id,1), synforc(id,2)
      end do 
      end do 
   end subroutine

   ! following functions are not used.
   !------------------------------------------
   real function compute_mean(mat)
      real, dimension(:,:), intent(in) :: mat
      compute_mean = sum(mat)/size(mat)
   end function

   !------------------------------------------
   function var_sqrt(A)
      type(forcing_variances) var_sqrt
      type(forcing_variances), intent(in) :: A

      var_sqrt%slp   = sqrt(A%slp   )
      var_sqrt%taux  = sqrt(A%taux  )
      var_sqrt%tauy  = sqrt(A%tauy  )
      var_sqrt%wndspd= sqrt(A%wndspd)
      var_sqrt%airtmp= sqrt(A%airtmp)
      var_sqrt%relhum= sqrt(A%relhum)
      var_sqrt%clouds= sqrt(A%clouds)
      var_sqrt%precip= sqrt(A%precip)
      var_sqrt%sss   = sqrt(A%sss   )
      var_sqrt%sst   = sqrt(A%sst   )
   end function var_sqrt

!------------------------------------------
   subroutine assign_force(A,r)
      type(forcing_fields), intent(out) :: A
      real, intent(in) :: r

      integer :: i,j

      do j=1,jdm
      do i=1,idm
         A%slp    (i,j) = r
         A%taux   (i,j) = r
         A%tauy   (i,j) = r
         A%wndspd (i,j) = r
         A%airtmp (i,j) = r
         A%relhum (i,j) = r
         A%clouds (i,j) = r
         A%precip (i,j) = r
         A%sss    (i,j) = r
         A%sst    (i,j) = r
         A%uwind  (i,j) = r
         A%vwind  (i,j) = r
         A%tauxice(i,j) = r
         A%tauyice(i,j) = r
      end do
      end do
   end subroutine assign_force

!------------------------------------------
   subroutine assign_vars(A,r)
      type(forcing_variances), intent(out) :: A
      real, intent(in) :: r
      A%slp    = r
      A%taux   = r
      A%tauy   = r
      A%wndspd = r
      A%airtmp = r
      A%relhum = r
      A%clouds = r
      A%precip = r
      A%sss    = r
      A%sst    = r
   end subroutine assign_vars

   subroutine randfld_rd(time_index)

      character(2)  :: time_index
      character(150) :: filename
      integer :: ix,jy

      filename = trim(iopath)//'/randfld.'//time_index

      if (debug) print*, 'reading randfile ', filename

      open(10, file=filename, status="old", action="read")

      do jy=1,jdm
        do ix=1,idm
        read(10, '(10e14.3)') &
                   ran%slp(ix,jy), ran%taux(ix,jy), ran%tauy(ix,jy), &
                   ran%wndspd(ix,jy), ran%airtmp(ix,jy), ran%relhum(ix,jy), &
                   ran%clouds(ix,jy), ran%precip(ix,jy), ran%sss(ix,jy), ran%sst(ix,jy)
        end do !ix
      end do !jy
      close(10)
      if (debug) print*, 'read from file'

   end subroutine

   subroutine synforc_rd(time_index)

         character(2)  :: time_index
         character(150) :: filename
         integer       :: ix,jy, xx, yy

         filename = trim(iopath)//'/synforc.'//time_index

         if (debug) print*, 'reading synforc ', filename

         open(11, file=filename, status="old", action="read")

         do jy=1,jdm
         do ix=1,idm
            read(11,'(6e14.3)') &
                     synuwind(ix,jy), synvwind(ix,jy), &
                     synairtmp(ix,jy), synslp(ix,jy), &
                     synprecip(ix,jy), synrelhum(ix,jy)
         end do !ix
         end do !jy

         close(11)

   end subroutine

   ! subroutine randfld_wr(time_index)

   !       character(2)  :: time_index
   !       character(150) :: filename
   !       integer :: ix,jy

   !       filename = trim(iopath)//'/randfld.'//time_index


   !       if (debug) print*, 'writing randfile ', filename

   !       open(12,file=filename,status='replace')

   !       do jy=1,jdm
   !       do ix=1,idm
   !          write(12,'(10e14.3)') &
   !                   ran%slp(ix,jy), ran%taux(ix,jy), ran%tauy(ix,jy), &
   !                   ran%wndspd(ix,jy), ran%airtmp(ix,jy), ran%relhum(ix,jy), &
   !                   ran%clouds(ix,jy), ran%precip(ix,jy), ran%sss(ix,jy), ran%sst(ix,jy)
   !       end do !ix
   !       end do !jy

   !       close(12)

   ! end subroutine

   ! subroutine synforc_wr(time_index)

   !       character(2)  :: time_index
   !       character(150) :: filename
   !       integer       :: ix,jy
   !       filename = trim(iopath)//'/synforc.'//time_index

   !       if (debug) print*, 'writing synforc ', filename

   !       open(13,file=filename, status='replace')

   !       do jy=1,jdm
   !       do ix=1,idm
   !          write(13,'(2i5,6e14.3)') ix,jy,  &
   !                   synuwind(ix,jy), synvwind(ix,jy), &
   !                   synairtmp(ix,jy), synslp(ix,jy), &
   !                   synprecip(ix,jy), synrelhum(ix,jy)
   !       end do !ix
   !       end do !jy

   !       close(13)

   ! end subroutine

end module mod_random_forcing
