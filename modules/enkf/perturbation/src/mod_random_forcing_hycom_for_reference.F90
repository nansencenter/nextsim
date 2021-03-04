      module mod_random_forcing 
#define FFTW 
      use mod_xc
      implicit none
c --- shield everything
      private
      type forcing_fields
         character(len=5) tforce
         real slp    (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) !  Sea level pressure
         real taux   (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) !  wind stress in x direction
         real tauy   (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) !  wind stress in y direction
         real wndspd (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) !  wind speed (tke source)
         real airtmp (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) !  pseudo air temperature
         real relhum (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) !  relative humidity
         real clouds (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) !  cloud cover
         real precip (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) !  precipitation
         real sss    (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) !  SSS for relax 
         real sst    (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) !  SST for relax 
         real uwind  (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) !  u-component of wind
         real vwind  (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) !  v-component of wind
         real tauxice(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) !  ice stress on water in x dir
         real tauyice(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) !  ice stress on water in y dir
      end type forcing_fields
c
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
c
c --- Operators
      interface operator(*)
         module procedure variances_forcing_fields_mult
      end interface
      interface assignment(=)
         module procedure assign_force
         module procedure assign_vars
      end interface
      interface sqrt
         module procedure var_sqrt
      end interface
c
      integer, save :: fnx, fny    ! fourier transform dimensions
      integer, save, public :: rf_prsflg=0 ! Pressure flag for random forcing
      logical, save, public :: randf       ! Switches on/off random forcing
      real   , save, public :: rf_hradius  ! Horizontal decorr length for rand forc [m]
      real   , save, public :: rf_tradius  ! Temporal decorr length for rand forc
      real   , save :: scorr       ! Horizontal decorr length for rand forc [grid cells]
      real   , save :: tcorr
c
c --- These will hold the forcing fields (dim and nondimensional)
      type(forcing_fields)   , save :: ran, ran1
c
c --- These will hold the forcing variances (input set in m_limits) 
      type(forcing_variances), save, public :: vars
c
c --- Last update of the rf forcing 
      real*8,save  :: dtimerf_last
c
c --- public entities
      public :: init_random_forcing, rand_update, ran
c --- TODO: random forcing now only works with inline option
      contains 

      subroutine init_random_forcing()
      use mod_xc
      use mod_forcing_nersc
      implicit none
      real :: dx
      include 'common_blocks.h'
      call xceget(dx,scpx,itdm/2,jtdm/2)
      if (mnproc==1)print *,'random forcing: typical model grid scale ',
     &   dx
      scorr=rf_hradius/dx     ! Decorrelation length is scorr grid cells
      tcorr=rf_tradius        ! Temporal decorrelation scale (days)
      ran =0.
      ran1=0.
      call set_random_seed2()
      call initfftdim(itdm,jtdm)
      call ranfields(ran)
      end subroutine
c
c --- Initialize FFT dimensions used in pseudo routines
      subroutine initfftdim(nx,ny)
      implicit none
      integer, intent(in) :: nx,ny
      fnx = ceiling(log(float(nx))/log(2.))
      fnx = 2**fnx
      fny = ceiling(log(float(ny))/log(2.))
      fny = 2**fny
      if (mnproc==1) then
         write(lp,'("Fourier transform dimensions ",2i6)')  fnx,fny
         call flush(lp)
      end if
      end subroutine
c
c
c --- Sets a random seed based on the wall clock time
      subroutine set_random_seed2
      use mod_xc, only : mnproc,xcstop, lp
      implicit none 
      integer , dimension(8)::val
      integer cnt
      integer sze
c --- Arrays for random seed
      integer, allocatable, dimension(:):: pt  
      real   , allocatable, dimension(:):: rpt
c
      call DATE_AND_TIME(values=val)
      call RANDOM_SEED(size=sze)
      allocate( pt(sze)) ; allocate(rpt(sze))
c --- Init - assumes seed is set in some way based on clock, 
c --- date etc. (not specified in fortran standard). Sometimes
c --- this initial seed is just set every second 
      call RANDOM_SEED   
c --- Retrieve initialized seed. val(8) is milliseconds - 
      call RANDOM_SEED(GET=pt) 
c --- this randomizes stuff if random_seed is not updated often 
c --- enough. synchronize seed across tasks (needed if pseudo 
c --- is paralellized some day)
      rpt = pt * (val(8)-500)  
      call xcastr(rpt,1)     
      pt=int(rpt)
#if defined MP_TEST_DUMP
      if (mnproc==1) then
      print *,'**********************************************'
      print *,'**********************************************'
      print *,'**********************************************'
      print *,'*** Warning *** Warning *** Warning **********'
      print *,'*** Random Seed set to 1 for MP_TEST_DUMP  ***'
      print *,'**********************************************'
      print *,'**********************************************'
      print *,'**********************************************'
      print *,'**********************************************'
      end if
      pt=1
#endif
      call RANDOM_SEED(put=pt)
      deallocate( pt)
      deallocate(rpt)
      end subroutine set_random_seed2
c
c
c
C --- This routine updates the random forcing component, according to 
C --- a simple correlation progression with target variance specified 
C --- by forcing_variances. At the end of the routine, the random 
c --- forcing is added to the forcing fields kept in index 2. This
c --- routine is only applicable to high-frequency forcing(yrflag >=2).
      subroutine rand_update(dtime,slt,init)
c --- TODO keep ran on first update (no update)
      use mod_xc
      use mod_year_info, only: year_info
      use mod_forcing_nersc
      use mod_hycom_nersc, only:imem
      implicit none
c
      real*8,  intent(in)    :: dtime
      integer, intent(in)    :: slt
      logical, intent(in)    :: init
c
      integer :: ix,jy
      real :: alpha, autocorr, nsteps, wspd
      real*8  :: dtimerf
      real, dimension(itdm,jtdm) :: modlon, modlat, gairtmp, gtaux, 
     &   gdepths , guwnd
      logical, parameter :: randf_test=.false.
      logical, save :: lfirst=.true.
      real, parameter :: rhoa = 1.2 , cdfac = 0.0012
      real :: cd_new,w4,wfact,wndfac, fcor
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: dpresx,dpresy
      real :: ucor, vcor, ueq,veq,wcor
      real :: wprsfac, minscpx, maxscpx
      real, parameter :: radtodeg=57.2957795
      real, parameter :: wlat=15.
      integer :: itst, jtst
      include  'common_blocks.h'
c
      if (yrflag<=2) then 
         if (mnproc==1) 
     &   write(lp,*) 'Random forcing only for yrflag=3'
         call xcstop('(rand_update)')
         stop '(rand_update)'
      end if
c
c --- This is the target time for the random forcing update
      if (slt==1) then
         hf_rdtime=6.d0/24.d0
         dtimerf=floor(dtime/hf_rdtime)*hf_rdtime
         dtimerf_last=0.
      else if (slt==2) then
         dtimerf=(floor(dtime/hf_rdtime)+1)*hf_rdtime
      else
         call xcstop('(rand_update called for slot > 2')
         stop '(rand_update called for slot > 2'
      end if
c
c --- On init we calculate rf regardless
      if (init) then
         continue
c --- This time has already been processed
c      else if (abs(dtimerf_last-dtimefrc_last)< 1e-4) then
      else if (dtime.le.dtimerf_last) then
         return
c --- rf forcing update is in order
      else
         continue
      end if
      dtimerf_last=dtimerf
c
c
c --- Autocorrelation between two times "tcorr"
c --- KAL - quite high? - autocorr = 0.95
c --- LB - Still too high      autocorr = 0.75
      autocorr = exp(-1.0)  ! 0.37

c
c --- Number of "forcing" steps rdtime between autocorrelation
c --- decay. Autocorrelation time is tcorr
c --- hf_rdtime is time step of a hf forcing read - from
c --- mod_forcing_nersc
      nsteps = tcorr/hf_rdtime
c
c --- This alpha will guarantee autocorrelation "autocorr"
c --- over the time "tcorr" 
c --- tcorr -> infinity , nsteps -> infinity, alpha -> 1 
c --- tcorr -> 0        , nsteps -> 0       , alpha -> 0 (when 1>autocorr>0)
      alpha=autocorr**(1/nsteps)    ! alpha = exp(-hf_rdtime/tcorr)
c
      if (mnproc==1) write(lp,*)
     &   'Rand_update -- Random forcing field update.'
c
c --- Calculate dimensionalized random perturbations
      ran1=sqrt(vars)*ran
c
c --- Apply dimensionalized random forcing to high-frequency forcing
c --- fields
      if (rf_prsflg .eq. 1 .or. rf_prsflg .eq.2 ) then
c --- ----------------------------------------------------------------
c --- rf_prsflag=1 : wind perturbations calculated from slp, using coriolis
c ---                parameter at 40 deg N
c --- rf_prsflag=2 : wind perturbations calculated from slp, using coriolis
c ---                parameter at 40 deg N, but limited by the setting of
c ---                windspeed, to account for the horizontal scale of pert.
c ---                As far as the wind is concerned, this is the same as 
c ---                reducing pressure perturbations
c --- ----------------------------------------------------------------
c 
      call xctilr(ran1%slp(1-nbdy,1-nbdy),1,1, 6,6, halo_ps)
c 
c --- grid size min max
      minscpx=minval(scpx,ip==1)
      maxscpx=minval(scpx,ip==1)
      call xcminr(minscpx)
      call xcmaxr(maxscpx)
      wprsfac=1.

c ---  flag used in prsflg=2
      if (rf_prsflg==2) then

c ---    Constant at 40 N
         fcor=2*sin(40./radtodeg)*2*pi/86400;
c
c ---    typical pressure gradient
         wprsfac=100.*sqrt(vars%slp)/(scorr*minscpx)
c
c ---    results in this typical wind magnitude
         wprsfac=wprsfac/fcor
c
c ---    but should give wind according to vars%wndspd
c ---    this is a correction factor for that
         wprsfac=sqrt(vars%wndspd)/wprsfac
      end if
c
c
c --- Calculate pressure gradient
      dpresx=0.
      dpresy=0.
!$OMP PARALLEL DO PRIVATE (ix,jy)
!$OMP&SCHEDULE(STATIC,jblk)
      do jy=1-nbdy+1,jj+nbdy-1
      do ix=1-nbdy+1,ii+nbdy-1

c ---    Pressure gradient. Coversion from mBar to Pa
c ---    Note that this is in u-point, but is used for a v-point
c ---    later 
         if (iu(ix,jy)==1)  then
            dpresx(ix,jy) =
     &         100.*(ran1%slp(ix,jy) - ran1%slp(ix-1,jy))/
     &         scux(ix,jy)
            dpresx(ix,jy)=dpresx(ix,jy)*wprsfac
         else
            dpresx(ix,jy)=0. ! should be extrapolated
         end if
c
c ---    Note that this is in v-point, but is used for a u-point
c ---    later
         if (iv(ix,jy)==1)  then
            dpresy(ix,jy) =
     &         100.*(ran1%slp(ix,jy) - ran1%slp(ix,jy-1))/
     &         scvy(ix,jy)
            dpresy(ix,jy)=dpresy(ix,jy)*wprsfac
         else
         dpresy(ix,jy)=0. ! should be extrapolated
         end if
      end do
      end do
!$OMP END PARALLEL DO



c --- Calculate winds from pressure gradient
!$OMP PARALLEL DO PRIVATE (ix,jy,fcor,ucor,vcor, 
!$OMP&                     ueq,veq,wcor)             
!$OMP&SCHEDULE(STATIC,jblk)
      do jy=1-nbdy+1,jj+nbdy-1
      do ix=1-nbdy+1,ii+nbdy-1
c
c ---    Coriolis balance (at 40 deg)
cKAL     fcor=2*sin(max(abs(plat(ix,jy)),20.)/radtodeg)*2*pi/86400;
         fcor=2*sin(40./radtodeg)*2*pi/86400; ! Constant 
         fcor=fcor*sign(1.,plat(ix,jy))*rhoa
         vcor= dpresx(ix,jy) / (fcor)
         ucor=-dpresy(ix,jy) / (fcor)
c
c ---    In the equatorial band u,v are aligned with the
c ---    pressure gradients. Here we use the coriolis
c ---    factor above to set it up (to limit the speeds)
         ueq=-dpresx(ix,jy) / abs(fcor)
         veq=-dpresy(ix,jy) / abs(fcor)
c
c ---    Weighting between coriiolis/equator solution
         wcor=sin(
     &      min(abs(plat(ix,jy)),wlat) / wlat
     &      *pi*0.5)
c ---    Perturbed winds
         uwind (ix,jy,slt) = uwind (ix,jy,slt)+wcor*ucor + (1.-wcor)*ueq
         vwind (ix,jy,slt) = vwind (ix,jy,slt)+wcor*vcor + (1.-wcor)*veq
         wndspd(ix,jy,slt) = sqrt(uwind(ix,jy,slt)**2 +
     &                            vwind(ix,jy,slt)**2)
c
c ---    The rest use fields independent of slp 
         airtmp(ix,jy,slt) = airtmp(ix,jy,slt)+ran1%airtmp(ix,jy)
         relhum(ix,jy,slt) = relhum(ix,jy,slt)+ran1%relhum(ix,jy)
         slp   (ix,jy,slt) = slp   (ix,jy,slt)+ran1%slp   (ix,jy)
         precip(ix,jy,slt) = precip(ix,jy,slt) ! lognormal precip
     &      *exp(ran1%precip(ix,jy) - 0.5*vars%precip**2) 
         clouds(ix,jy,slt) = clouds(ix,jy,slt)+ran1%clouds(ix,jy)
         relhum(ix,jy,slt) = min(max(relhum(ix,jy,slt),0.0),1.0)
         clouds(ix,jy,slt) = min(max(clouds(ix,jy,slt),0.0),1.0)
         wndspd(ix,jy,slt) = max(wndspd(ix,jy,slt),0.0)
      end do
      end do
!$OMP END PARALLEL DO
c
c
c ---    New drag - Computed directly from winds now
!$OMP PARALLEL DO PRIVATE (ix,jy,wndfac,cd_new,w4,wfact) 
!$OMP&SCHEDULE(STATIC,jblk)
      do jy=1-nbdy+1,jj+nbdy-1
      do ix=1-nbdy+1,ii+nbdy-1
         wndfac=(1.+sign(1.,wndspd(ix,jy,slt)-11.))*.5
         cd_new=(0.49+0.065*wndspd(ix,jy,slt))*1.0e-3*wndfac+cdfac*
     &          (1.-wndfac)
c
         w4    =.25*(vwind(ix-1,jy+1,slt)+vwind(ix,jy+1,slt)+ 
     &               vwind(ix-1,jy  ,slt)+vwind(ix,jy  ,slt))
         wfact=sqrt(uwind(ix,jy,slt)*uwind(ix,jy,slt)+w4*w4)*
     &         airdns*cd_new
         taux(ix,jy,slt)=uwind(ix,jy,slt)*wfact
c
         w4   =.25*(uwind(ix  ,jy-1,slt)+uwind(ix+1,jy-1,slt)+
     &              uwind(ix+1,jy  ,slt)+uwind(ix  ,jy  ,slt))
         wfact=sqrt(vwind(ix,jy,slt)*vwind(ix,jy,slt)+w4*w4)*
     &         airdns*cd_new
         tauy(ix,jy,slt)=vwind(ix,jy,slt)*wfact
      end do
      end do
!$OMP END PARALLEL DO
c
c

c ---    ---------------------------------------------------------
c ---    rf_prsflag=0 : wind and slp are uncorrelated
c ---    ---------------------------------------------------------
      else  ! rf_prsflg .eq. 0
c
c
c
!$OMP PARALLEL DO PRIVATE (ix,jy,wspd) SCHEDULE(STATIC,jblk)
         do jy=1-nbdy,jj+nbdy
         do ix=1-nbdy,ii+nbdy
         !if (ip(ix,jy)==1) then
            taux  (ix,jy,slt) = taux  (ix,jy,slt) + ran1%taux(ix,jy)
            tauy  (ix,jy,slt) = tauy  (ix,jy,slt) + ran1%tauy(ix,jy)

c ---       winds are nonlinear functions of tau and mainly
c ---       used for sea ice
            wspd = sqrt(taux(ix,jy,slt)**2 + tauy(ix,jy,slt)**2)
            wspd = max(sqrt(wspd / (cdfac*rhoa)),0.1)
            uwind(ix,jy,slt) = taux (ix,jy,slt) / (wspd*cdfac*rhoa)
            vwind(ix,jy,slt) = tauy (ix,jy,slt) / (wspd*cdfac*rhoa)
c
            airtmp(ix,jy,slt) = airtmp(ix,jy,slt)+ran1%airtmp(ix,jy)
            wndspd(ix,jy,slt) = wndspd(ix,jy,slt)+ran1%wndspd(ix,jy)
            relhum(ix,jy,slt) = relhum(ix,jy,slt)+ran1%relhum(ix,jy)
            precip(ix,jy,slt) = precip(ix,jy,slt) ! Lognormal 
     &         *exp(ran1%precip(ix,jy) - 0.5*vars%precip**2) 
            relhum(ix,jy,slt) = min(max(relhum(ix,jy,slt),0.0),1.0)
            wndspd(ix,jy,slt) = max(wndspd(ix,jy,slt),0.0)
         !end if
         end do
         end do
      end if ! rf_prsflg
c --- ---------------------------------------------------------
c --- end if rf_prsflag=0
c --- ---------------------------------------------------------
c
c --- Synchronize tiles of new forcing
      call xctilr(uwind (1-nbdy,1-nbdy,slt),1,   1, 6,6, halo_uv)
      call xctilr(vwind (1-nbdy,1-nbdy,slt),1,   1, 6,6, halo_vv)
      call xctilr(taux  (1-nbdy,1-nbdy,slt),1,   1, 6,6, halo_uv)
      call xctilr(tauy  (1-nbdy,1-nbdy,slt),1,   1, 6,6, halo_vv)
      call xctilr(airtmp(1-nbdy,1-nbdy,slt),1,   1, 6,6, halo_ps)
      call xctilr(wndspd(1-nbdy,1-nbdy,slt),1,   1, 6,6, halo_ps)
      call xctilr(relhum(1-nbdy,1-nbdy,slt),1,   1, 6,6, halo_ps)
      call xctilr(precip(1-nbdy,1-nbdy,slt),1,   1, 6,6, halo_ps)
      call xctilr(slp   (1-nbdy,1-nbdy,slt),1,   1, 6,6, halo_ps)
c
c --- ran1 is new random nondimensional random forcing
      call ranfields(ran1)
c
c --- This calculates nondimensional random forcing ran as a AR(1) 
c --- process. ran= alpha*ran + sqrt(1-alpha*alpha)* ran 
      call ran_update_ran1(ran,ran1,alpha)
c
c --- ---------------------------------------------------------
c --- Diagnostics section -- Spatial field dumped on first run
c --- ---------------------------------------------------------
      if (lfirst) then
         ! Test for ranfld .... only done for first pass
         lfirst=.false.
c
         call xcaget(modlon,plon,0)
         call xcaget(modlat,plat,0)
         call xcaget(gdepths,depths,0)
         call xcaget(gairtmp,ran1%airtmp*sqrt(vars%airtmp),0)
         call xcaget(gtaux  ,ran1%taux  *sqrt(vars%taux  ),0)
         call xcaget(guwnd  ,uwind(:,:,slt)                 ,0)
         if (mnproc==1 .and. imem==1) then
            open(10,file='ranfld.dat',status='replace')
            do jy=1,jtdm
            do ix=1,itdm
               write(10,'(2i5,6e14.3)') ix,jy, 
     &         modlon(ix,jy),modlat(ix,jy),gdepths(ix,jy),
     &         gairtmp(ix,jy), gtaux(ix,jy),guwnd(ix,jy)
            end do
            end do
            close(10)
         end if
      end if
c --- ---------------------------------------------------------
c --- Diagnostics section -- time series at test point
c --- ---------------------------------------------------------
      if (randf_test) then
         itst=304
         jtst=796
         if (i0+1<=itst .and. i0+ii>=itst .and. 
     &       j0+1<=jtst .and. j0+jj>=jtst ) then
            open(10,file='rantest.dat',status='unknown',
     &         position='append')
            write(10,'(f18.6,5f14.3)') 
     &         dtimerf,
     &         ran%taux  (itst-i0,jtst-j0)*sqrt(vars%taux),
     &         ran%tauy  (itst-i0,jtst-j0)*sqrt(vars%tauy),
     &         ran%airtmp(itst-i0,jtst-j0)*sqrt(vars%airtmp),
     &         ran%taux(itst-i0,jtst-j0),
     &         ran%tauy(itst-i0,jtst-j0)
            close(10)
         end if
         call inline_diag('hfperturb',dtimerf,slt)
      end if
c --- ---------------------------------------------------------
c --- Diagnostics section -- time series at test point
c --- ---------------------------------------------------------
      end subroutine rand_update
c
c
c --- Calculate a set of nondimensional random forcing fields.
c --- This routine can be sep up in MPI - but so far this doesnt
c --- seem to be a major problem
      subroutine ranfields(ranfld)
      implicit none
c
      type(forcing_fields)    , intent(inout) :: ranfld 
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: tmp
      real, dimension(itdm,jtdm) :: gtmp
c
      ranfld=0.
      call pseudo2D(gtmp,itdm,jtdm,1,scorr,fnx,fny)
      call xcaput(gtmp, ranfld%slp   ,1) ! Node 1 dumps
c
      call pseudo2D(gtmp,itdm,jtdm,1,scorr,fnx,fny)
      call xcaput(gtmp, ranfld%taux  ,1) ! Node 1 dumps
c
      call pseudo2D(gtmp,itdm,jtdm,1,scorr,fnx,fny)
      call xcaput(gtmp, ranfld%tauy  ,1) ! Node 1 dumps
c
      call pseudo2D(gtmp,itdm,jtdm,1,scorr,fnx,fny)
      call xcaput(gtmp, ranfld%wndspd,1) ! Node 1 dumps
c
      call pseudo2D(gtmp,itdm,jtdm,1,scorr,fnx,fny)
      call xcaput(gtmp, ranfld%airtmp,1) ! Node 1 dumps
c
      call pseudo2D(gtmp,itdm,jtdm,1,scorr,fnx,fny)
      call xcaput(gtmp, ranfld%relhum,1) ! Node 1 dumps
c
      call pseudo2D(gtmp,itdm,jtdm,1,scorr,fnx,fny)
      call xcaput(gtmp, ranfld%clouds,1) ! Node 1 dumps
c
      call pseudo2D(gtmp,itdm,jtdm,1,scorr,fnx,fny)
      call xcaput(gtmp, ranfld%precip,1) ! Node 1 dumps
c
      end subroutine ranfields






c --- This routine calculates the pseudo random filds using
c --- the procedure outlined in Evensen (1994) \cite{eve94a}.
c --- TODO - remove legacy cruft
c --- TODO - doublecheck FFTW approach
      subroutine pseudo2D(Amat,nx,ny,lde,rh,n1,n2)
c
      implicit none
      integer, intent(in) :: nx,ny           ! horizontal dimensions
      integer, intent(in) :: lde             ! number of fields to create
      real, intent(out)   :: Amat(nx,ny,lde) ! generated random fields
      real, intent(in)    :: rh              ! Horizontal decorrelation length
      integer, intent(in) :: n1,n2           ! horizontal dimensions in fft grid
c --- Its machine definition time ...
#ifdef AIX
      integer :: maxnn,s1
      integer,save :: naux1,naux2,naux3
      double precision, dimension(:), allocatable,save :: aux1,aux2,aux3
#endif
!#ifdef SGI Commented because not working with Vilje Ve
!      real, allocatable, dimension(:), save :: coeff
!#endif
#if defined(FFTW)
      integer*8, save :: plan
      real*8 :: fftwy(n1,n2)
      include 'fftw3.f'
#endif
c
c --- saving rh used in preprocessing.  Allow for new call to
c --- zeroin if rh changes.
      real, save ::  rh_save=0.0  
c
      real, save :: sigma,sigma2
      real, save :: c
      integer, save :: n1_save=0
      integer, save :: n2_save=0
c
      integer l,p,j,n,m,i
      real kappa2,lambda2,kappa,lambda
      real pi2,deltak,sum,scale
      real a1,b1,tol,fval
c
      real, allocatable    :: fampl(:,:,:)
      real, allocatable    :: phi(:,:)
      real, allocatable    :: y(:,:)   ! Physical field
      complex*16, allocatable :: x(:,:)   ! Fourier amplitudes
      complex, allocatable :: xx(:,:)   ! complex physical field
c
      real, parameter :: dx=1.0
      real, parameter :: pi=3.141592653589
c
      if (lde < 1) then
         call xcstop('pseudo2D: error lde < 1')
         stop 'pseudo2D: error lde < 1'
      end if
      if (rh <= 0.0) then
         call xcstop('pseudo2D: error, rh <= 0.0')
         stop 'pseudo2D: error, rh <= 0.0'
      end if
      allocate(fampl(0:n1/2,-n2/2:n2/2,2))
      allocate(phi(0:n1/2,-n2/2:n2/2))
      allocate(y(0:n1+1,0:n2-1))
      allocate(x(0:n1/2,0:n2-1))
c --- machine tests. Changed to tell at compile-time
#ifndef AIX
#ifndef SGI
#if !defined(FFTW)
!      if (mnproc==1) 
!      print *,'ranfield is only running on the following machines:'
!      print *,'   AIX having essl'
!      print *,'   SGI'
!      print *,'   Any system having FFTW'
!      end if
!      call xcstop ('(pseudo2D)')
!      stop  '(pseudo2D)'
#error No fft method in mod_random_forcing
#endif
#endif
#endif
c
      pi2=2.0*pi
      deltak=pi2**2/(float(n1*n2)*dx*dx)
      kappa=pi2/(float(n1)*dx)
      kappa2=kappa**2
      lambda=pi2/(float(n2)*dx)
      lambda2=lambda**2
      scale=1.0
c
      if (rh /= rh_save .or. n1 /= n1_save .or. n2 /= n2_save) then
         rh_save=rh
         n1_save=n1
         n2_save=n2
!#ifdef SGI commented because not working with Ve, Vilje
!         if(allocated(coeff)) deallocate(coeff)
!         allocate( coeff((n1+15) + 2*(n2+15)) )
!         call dzfft2dui(n1,n2,coeff)
!#endif
#ifdef AIX
         maxnn   = max(n1/2,n2)
         if (maxnn<=2048) then
            naux1= 42000 
         else
            naux1= ceiling(40000+1.64*n1+2.28*n2)
         end if
         if (n1 <= 4096 ) then
            naux2 = 20000
         else if (n1 > 4096 ) then 
            naux2 = ceiling(20000+1.14*n1)
         end if
         if ( n2 > 252) then 
            s1 = min(64, 1+n1/2)
            naux2 = naux2 + ceiling((2*n2+256)*(2.28+s1))
         end if
         naux3=1
         allocate(aux1(naux1))
         allocate(aux2(naux2))
         allocate(aux3(naux3))
         call dcrft2(1,x,n1/2+1,y,n1+2,n1,n2,-1,scale,
     &            aux1,naux1,aux2,naux2,aux3,naux3)
#endif
#if defined(FFTW)
         if (mnproc==1) then
            print *,'Using FFTW for fourier transform'
         end if
#endif
         rh_save=rh
         if (mnproc==1) print *,'Ranfield: Solving for sigma',rh,dx
         a1=0.1e-07
         b1=0.1e-06
         tol=0.1e-10
         call zeroin(func2D,sigma,a1,b1,tol,rh,dx,fval,n1,n2)
c
         sigma2=sigma**2
         sum=0.0
         do p=-n2/2+1,n2/2
         do l=-n1/2+1,n1/2
            sum=sum+exp(
     &         -2.0*(kappa2*float(l*l)+lambda2*float(p*p))/sigma2
     &                 )
         enddo
         enddo
         c=sqrt(1.0/(deltak*sum))
c
         if (mnproc==1) print *,'Ranfield: sigma  ',sigma
         if (mnproc==1) print *,'Ranfield: c=     ',c
      endif ! if rh /= rh_save ....

      do j=1,lde

c ---    Calculating the random wave phases
         call random_number(phi)
         phi=pi2*phi
c
c ---    Calculating the wave amplitues
         do p=-n2/2,n2/2
         do l=0,n1/2 
            fampl(l,p,1)=
     &         exp(-(kappa2*float(l*l)+lambda2*float(p*p))/sigma2)*
     &         cos(phi(l,p))*sqrt(deltak)*c
            fampl(l,p,2)=
     &         exp(-(kappa2*float(l*l)+lambda2*float(p*p))/sigma2)*
     &         sin(phi(l,p))*sqrt(deltak)*c
         enddo
         enddo
         fampl(0,0,2)=0.0
c ---    Put into fft arrays
         do p=0,n2/2-1
            x(:,p)=cmplx(fampl(:,p,1),fampl(:,p,2))
         enddo
         do p=n2/2,n2-1
            x(:,p)=cmplx(fampl(:,-n2+p,1),fampl(:,-n2+p,2))
         enddo
c---     FFT for different machines
!#ifdef SGI commented because not working on Ve Vilje
!         call zdfft2du(-1,n1,n2,x,n1+2,coeff)
!         y=reshape(transfer(x,(/0.0 /) ),(/n1+2,n2/))
!#endif
#ifdef AIX
         call dcrft2(0,x,n1/2+1,y,n1+2,n1,n2,-1,scale,
     &            aux1,naux1,aux2,naux2,aux3,naux3)
#endif
#if defined(FFTW)
         call dfftw_plan_dft_c2r_2d(plan,n1,n2,x,fftwy,FFTW_ESTIMATE)
         call dfftw_execute(plan)
         call dfftw_destroy_plan(plan)
         y(0:n1-1 ,0:n2-1)=fftwy(1:n1,1:n2)
         y(n1:n1+1,0:n2-1)=fftwy(1:2 ,1:n2)
#endif
c
         do m=1,ny
         do i=1,nx
            Amat(i,m,j)=y(i-1,m-1)
         enddo
         enddo
      enddo ! do 1..lde
      deallocate(fampl, phi, y, x)
      end subroutine pseudo2D


      subroutine zeroin(func,zeropkt,ax,bx,tol,length,dx,fval,n1,n2)
      implicit none
c --- Finds zero of function f.
c --- A zero of the function  $func(x,length,dx,n1,n2)$ is computed 
c --- in the interval $[ax,bx]$. Zeroin returns a zero $x$ in the
c --- given interval to within a tolerance  $4*macheps*abs(x) + tol$,
c --- where macheps is the relative machine precision.
c
c --- This function subprogram is a slightly  modified  translation  
c --- of the algol 60 procedure  zero  given in  richard brent,
c --- algorithms for c --- minimization without derivatives, prentice
c --- hall, inc. (1973).
c
      real, external :: func
c
      integer n1,n2, icorr
      real zeropkt,length,dx
c --- ax=left endpoint of initial interval, bx=right endpoint 
c --- of initial interval tol=desired length of the interval of uncertainty of the
      real ax  
      real bx  
      real tol  
      real  a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
      real  abs,sign,fval

c --- compute eps, the relative machine precision

#ifdef DEBUG
      if (mnproc==1) write(*,*)'A=',ax,bx,tol,length,dx
#endif
      icorr=0
c
      eps = 1.0
 10   eps = eps/2.0
      tol1 = 1.0 + eps
      if (tol1 .gt. 1.0) go to 10
c
c --- initialization
 77   a = ax
      b = bx
      fa = func(a,length,dx,n1,n2)
      fb = func(b,length,dx,n1,n2)
c

      if (fa*fb.gt.0.0) then
#ifdef DEBUG
         if (mnproc==1) then
            write(*,*)'fa=',fa
            write(*,*)'fb=',fb
            write(*,*)'fa*fb =',fa*fb,'is greater than zero'
         end if
#endif
         ax=0.1*ax
         bx=10.0*bx
         icorr=icorr+1
         if (icorr < 20) then
            goto 77
         else
            if (mnproc==1) write(*,'(2(a,g13.5))')
     &         'zeroin: No convergence, ax=',ax,' bx=',bx
            call xcstop('(zeroin)')
            stop  '(zeroin)'
         endif 
      endif ! if fa*fb gt ....
c
c --- begin step

 20   c = a
      fc = fa
      d = b - a
      e = d
 30   if (abs(fc) .ge. abs(fb)) go to 40
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
c
c --- convergence test
 40   tol1 = 2.0*eps*abs(b) + 0.5*tol
      xm = .5*(c - b)
      if (abs(xm) .le. tol1) go to 90
      if (fb .eq. 0.0) go to 90
c
c --- is bisection necessary
      if (abs(e) .lt. tol1) go to 70
      if (abs(fa) .le. abs(fb)) go to 70
c
c --- is quadratic interpolation possible
      if (a .ne. c) go to 50
c
c --- linear interpolation
      s = fb/fa
      p = 2.0*xm*s
      q = 1.0 - s
      go to 60
c
c --- inverse quadratic interpolation
 50   q = fa/fc
      r = fb/fc
      s = fb/fa
      p = s*(2.0*xm*q*(q - r) - (b - a)*(r - 1.0))
      q = (q - 1.0)*(r - 1.0)*(s - 1.0)
c
c --- adjust signs
 60   if (p .gt. 0.0) q = -q
      p = abs(p)
c
c --- is interpolation acceptable
      if ((2.0*p) .ge. (3.0*xm*q - abs(tol1*q))) go to 70
      if (p .ge. abs(0.5*e*q)) go to 70
      e = d
      d = p/q
      go to 80
c
c --- bisection
 70   d = xm
      e = d
c
c --- complete step
 80   a = b
      fa = fb
      if (abs(d) .gt. tol1) b = b + d
      if (abs(d) .le. tol1) b = b + sign(tol1, xm)
      fb = func(b,length,dx,n1,n2)
      if ((fb*(fc/abs(fc))) .gt. 0.0) go to 20
      go to 30
c
c --- done
 90   zeropkt = b
      fval=func(b,length,dx,n1,n2)
      end subroutine zeroin




      real function func2D(sigma,length,dx,n1,n2)
c --- Function used to calculate $sigma$ and $c$.
      implicit none
      real sum1,sum2,sigma,length
      real sigma2,pi2,kappa,kappa2,lambda,lambda2,dx
      integer l,p,n1,n2
      real, parameter :: pi=3.141592653589
c
      sigma2=sigma**2
c
      pi2=2.0*pi
      kappa=pi2/(float(n1)*dx)
      kappa2=kappa**2
c
      lambda=pi2/(float(n2)*dx)
      lambda2=lambda**2
c

c --- Calculate sum1
      sum1=0.0
      do p=-n2/2+1,n2/2
      do l=-n1/2+1,n1/2
         sum1=sum1+
     &      exp(-2.0*(kappa2*float(l*l)+lambda2*float(p*p))/sigma2)
     &      *cos(kappa*float(l)*length)
      enddo
      enddo
c
c --- Calculate sum2
      sum2=0.0
      do p=-n2/2+1,n2/2
      do l=-n1/2+1,n1/2
         sum2=sum2+
     &      exp(-2.0*(kappa2*float(l*l)+lambda2*float(p*p))/sigma2)
      enddo
      enddo
c
      func2D = sum1/sum2 - exp(-1.0)
      end function func2D
c
c --- -----------------------------------------------------------------
c --- -----------------------------------------------------------------
c --- -----Operator and auxillary routines for variances---------------
c --- -----------------------------------------------------------------
c --- -----------------------------------------------------------------
      function variances_forcing_fields_mult(A,B)
      implicit none
      type(forcing_fields) variances_forcing_fields_mult
      type(forcing_variances), intent(in) :: A
      type(forcing_fields), intent(in) :: B
      integer :: i,j
      do j=1-nbdy,jj+nbdy
      do i=1-nbdy,ii+nbdy
        variances_forcing_fields_mult%slp   (i,j)=A%slp   *B%slp   (i,j)
        variances_forcing_fields_mult%taux  (i,j)=A%taux  *B%taux  (i,j)
        variances_forcing_fields_mult%tauy  (i,j)=A%tauy  *B%tauy  (i,j)
        variances_forcing_fields_mult%wndspd(i,j)=A%wndspd*B%wndspd(i,j)
        variances_forcing_fields_mult%airtmp(i,j)=A%airtmp*B%airtmp(i,j)
        variances_forcing_fields_mult%relhum(i,j)=A%relhum*B%relhum(i,j)
        variances_forcing_fields_mult%clouds(i,j)=A%clouds*B%clouds(i,j)
        variances_forcing_fields_mult%precip(i,j)=A%precip*B%precip(i,j)
        variances_forcing_fields_mult%sss   (i,j)=A%sss   *B%sss   (i,j)
        variances_forcing_fields_mult%sst   (i,j)=A%sst   *B%sst   (i,j)
      end do
      end do
      end function variances_forcing_fields_mult
c
      function var_sqrt(A)
      implicit none
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
c
      subroutine assign_force(A,r)
      implicit none
      type(forcing_fields), intent(out) :: A
      real, intent(in) :: r
      integer :: i,j
      do j=1-nbdy,jj+nbdy
      do i=1-nbdy,ii+nbdy
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
c
      subroutine assign_vars(A,r)
      implicit none
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
c
      subroutine ran_update_ran1(ran,ran1,alpha)
      implicit none
      type(forcing_fields), intent(inout) :: ran
      type(forcing_fields), intent(   in) :: ran1
      real                , intent(   in) :: alpha
      integer :: ix,jy
      do jy=1-nbdy,jj+nbdy
      do ix=1-nbdy,ii+nbdy
         ran%slp   (ix,jy)=alpha*ran%slp   (ix,jy) + 
     &      sqrt(1-alpha*alpha)*ran1%slp   (ix,jy)
         ran%taux  (ix,jy)=alpha*ran%taux  (ix,jy) + 
     &      sqrt(1-alpha*alpha)*ran1%taux  (ix,jy)
         ran%tauy  (ix,jy)=alpha*ran%tauy  (ix,jy) +  
     &      sqrt(1-alpha*alpha)*ran1%tauy  (ix,jy)
         ran%wndspd(ix,jy)=alpha*ran%wndspd(ix,jy) + 
     &      sqrt(1-alpha*alpha)*ran1%wndspd(ix,jy)
         ran%airtmp(ix,jy)=alpha*ran%airtmp(ix,jy) + 
     &      sqrt(1-alpha*alpha)*ran1%airtmp(ix,jy)
         ran%relhum(ix,jy)=alpha*ran%relhum(ix,jy) +
     &      sqrt(1-alpha*alpha)*ran1%relhum(ix,jy)
         ran%clouds(ix,jy)=alpha*ran%clouds(ix,jy) +
     &      sqrt(1-alpha*alpha)*ran1%clouds(ix,jy)
         ran%precip(ix,jy)=alpha*ran%precip(ix,jy) + 
     &      sqrt(1-alpha*alpha)*ran1%precip(ix,jy)
         ran%sss   (ix,jy)=alpha*ran%sss   (ix,jy) + 
     &      sqrt(1-alpha*alpha)*ran1%sss   (ix,jy)
         ran%sst   (ix,jy)=alpha*ran%sst   (ix,jy) + 
     &      sqrt(1-alpha*alpha)*ran1%sst   (ix,jy)
      end do
      end do
      end subroutine
c
      end module mod_random_forcing
