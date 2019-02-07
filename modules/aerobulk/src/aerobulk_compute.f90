! Rehashed aerobulk_compute for neXtSIM
! AeroBulk / 2015 / L. Brodeau (brodeau@gmail.com)
! https://sourceforge.net/p/aerobulk

MODULE mod_aerobulk_compute_nextsim

   use iso_c_binding

   USE mod_const       !: physical constants
   USE mod_thermo      !: thermodynamics functions

   USE mod_blk_coare   !: COAREv3   algorithm
   USE mod_blk_ncar    !: Large & Yeager algorithm
   USE mod_blk_ecmwf   !: following ECMWF doc...

   IMPLICIT NONE

   PUBLIC :: aerobulk_nextsim_skin, aerobulk_nextsim_no_skin

   PRIVATE :: aerobulk_compute_nextsim

CONTAINS

  subroutine aerobulk_nextsim_skin( calgo, zt, zu, sst, t_zt, &
                                    q_zt, U_zu, slp,          &
                                    QL, QH, Cd_rho_U,         &
                                    l, m,                     &
                                    rad_sw, rad_lw ) bind(c)

     INTEGER,                           INTENT(in)  :: l, m
     CHARACTER(c_char), DIMENSION(l+1), INTENT(in)  :: calgo
     REAL(c_double),                    INTENT(in)  :: zt, zu
     REAL(c_double),    DIMENSION(m+1), INTENT(in)  :: sst, t_zt, q_zt, U_zu, slp
     REAL(c_double),    DIMENSION(m+1), INTENT(out) :: QL, QH, Cd_rho_U
     REAL(c_double),    DIMENSION(m+1), INTENT(in)  :: rad_sw, rad_lw

     ! Locals
     character(len=l) :: calgo_fort
     integer :: i

     ! Loop over the input string to convert to Fortran string
     do i=1,l
       calgo_fort(i:i) = calgo(i)
     enddo

     ! Do init
     jpi = m; jpj = 1;

     ! Call the actual routine
     call aerobulk_compute_nextsim( calgo_fort, zt, zu, sst, t_zt, &
                                    q_zt, U_zu, slp,               &
                                    QL, QH, Cd_rho_U,              &
                                    m,                             &
                                    rad_sw=rad_sw, rad_lw=rad_lw )

  end subroutine aerobulk_nextsim_skin

  subroutine aerobulk_nextsim_no_skin( calgo, zt, zu, sst, t_zt, &
                                       q_zt, U_zu, slp,          &
                                       QL, QH, Cd_rho_U,         &
                                       l, m) bind(c)

     INTEGER,                           INTENT(in)  :: l, m
     CHARACTER(c_char), DIMENSION(l+1), INTENT(in)  :: calgo
     REAL(c_double),                    INTENT(in)  :: zt, zu
     REAL(c_double),    DIMENSION(m,1), INTENT(in)  :: sst, t_zt, q_zt, U_zu, slp
     REAL(c_double),    DIMENSION(m,1), INTENT(out) :: QL, QH, Cd_rho_U

     ! Locals
     character(len=l) :: calgo_fort
     integer :: i

     ! Loop over the input string to convert to Fortran string
     do i=1,l
       calgo_fort(i:i) = calgo(i)
     enddo

     ! Do init
     jpj = m; jpi = 1;

     ! Call the actual routine
     call aerobulk_compute_nextsim( calgo_fort, zt, zu, sst, t_zt, &
                                    q_zt, U_zu, slp,               &
                                    QL, QH, Cd_rho_U,              &
                                    m)

  end subroutine aerobulk_nextsim_no_skin


  SUBROUTINE aerobulk_compute_nextsim( calgo, zt, zu, sst, t_zt, &
                                       q_zt, U_zu, slp,          &
                                       QL, QH, Cd_rho_U,         &
                                       m,                        &
                                       rad_sw, rad_lw )

     !!
     !!******************************
     !! 2015: L. Brodeau (brodeau@gmail.com)
     !!  => all constants taken from mod_thermo and mod_const must
     !!     be done or used from NEMO constant bank...    ... vkarmn ... grav ...
     !!******************************
     !!
     !!======================================================================================
     !!
     !! INPUT :
     !! -------
     !!    *  calgo: what bulk algorithm to use => 'coare'/'coare35'/'ncar'/'ecmwf'
     !!    *  zt   : height for temperature and spec. hum. of air           [m]
     !!    *  zu   : height for wind (10m = traditional anemometric height  [m]
     !!    *  sst  : bulk SST                                               [K]
     !!    *  t_zt : air temperature at zt                                  [K]
     !!    *  q_zt : specific humidity of air at zt                         [kg/kg]
     !!    *  U_zu : wind speed at zu                                       [m/s]
     !!    *  slp  : mean sea-level pressure                                [Pa] ~101000 Pa
     !!
     !! OPTIONAL INPUT (will trigger l_use_skin=TRUE if present!):
     !! ---------------
     !!    *  rad_sw : downwelling shortwave radiation at the surface (>0)  [W/m^2]
     !!    *  rad_lw : downwelling longwave radiation at the surface  (>0)  [W/m^2]
     !!
     !! OUTPUT :
     !! --------
     !!    *  QL       : Latent heat flux                                   [W/m^2]
     !!    *  QH       : Sensible heat flux                                 [W/m^2]
     !!    *  Cd_rho_U : Drag parameter (Cd*XRHO*XUblk)                     [N/(m/s)]
     !!
     !!============================================================================
     !!
     !! I/O ARGUMENTS:
     CHARACTER(len=*),         INTENT(in)  :: calgo
     REAL(wp),                 INTENT(in)  :: zt, zu
     INTEGER,                  INTENT(in)  :: m
     REAL(wp), DIMENSION(m,1), INTENT(in)  :: sst, t_zt, q_zt, U_zu, slp
     REAL(wp), DIMENSION(m,1), INTENT(out) :: QL, QH, Cd_rho_U
     REAL(wp), DIMENSION(m,1), INTENT(in), OPTIONAL :: rad_sw, rad_lw

     REAL(wp), DIMENSION(m,1) ::  &
        &   XSSQ,              & !: Specific humidiyt at the air-sea interface
        &   Cd, Ch, Ce,        & !: bulk transfer coefficients
        &  XTzt,               & !: potential temperature at zt meters
        &  XTzu, XQzu,         & !: potential temperature and specific humidity at zu meters
        &  Ts, qs,             & !:
        &   XUblk,             & !: Bulk scalar wind speed (U_zu corrected for low wind and unstable conditions)
        &   XRHO                 !: density of air

     LOGICAL :: l_use_skin = .FALSE.
     !!------------------------------------------------------------------------------

     ! Cool skin ?
     IF( PRESENT(rad_sw) .AND. PRESENT(rad_lw) ) THEN
        IF((TRIM(calgo) == 'coare').OR.(TRIM(calgo) == 'coare35').OR.(TRIM(calgo) == 'ecmwf')) THEN
           l_use_skin = .TRUE.
           ! PRINT *, ' *** Will use the cool-skin warm-layer scheme of ', TRIM(calgo(1:5)), '!'
        END IF
     END IF

     !! Computing specific humidity at saturation at sea surface temperature :
     XSSQ (:,:) = 0.98*q_sat(sst, slp)

     !! Approximate potential temperarure at zt meters high:
     XTzt = t_zt + gamma_moist(t_zt, q_zt)*zt

     !! Mind that TURB_COARE and TURB_ECMWF will modify SST and SSQ if their
     !! respective Cool Skin Warm Layer parameterization is used
     Ts = sst ; qs = XSSQ


     SELECT CASE(TRIM(calgo))
        !!
     CASE('coare')
        IF( l_use_skin ) THEN
           CALL TURB_COARE ( '3.0', zt, zu, Ts, XTzt, qs, q_zt, U_zu,  &
              &              Cd, Ch, Ce, XTzu, XQzu, XUblk,            &
              &              rad_sw=rad_sw, rad_lw=rad_lw, slp=slp )
        ELSE
           CALL TURB_COARE ( '3.0', zt, zu, Ts, XTzt, qs, q_zt, U_zu,  &
              &              Cd, Ch, Ce, XTzu, XQzu, XUblk )
        END IF
        !!
     CASE('coare35')
        IF( l_use_skin ) THEN
           CALL TURB_COARE ( '3.5', zt, zu, Ts, XTzt, qs, q_zt, U_zu, &
              &              Cd, Ch, Ce, XTzu, XQzu, XUblk,           &
              &              rad_sw=rad_sw, rad_lw=rad_lw, slp=slp )
        ELSE
           CALL TURB_COARE ('3.5', zt, zu, Ts, XTzt, qs, q_zt, U_zu,  &
              &              Cd, Ch, Ce, XTzu, XQzu, XUblk )
        END IF
        !!
        !!
     CASE('ncar')
        CALL TURB_NCAR( zt, zu, Ts, XTzt, qs, q_zt, U_zu, &
           &            Cd, Ch, Ce, XTzu, XQzu, XUblk)
        !!
        !!
     CASE('ecmwf')
        IF( l_use_skin ) THEN
           CALL TURB_ECMWF ( zt, zu, Ts, XTzt, qs, q_zt, U_zu,   &
              &              Cd, Ch, Ce, XTzu, XQzu, XUblk,      &
              &              rad_sw=rad_sw, rad_lw=rad_lw, slp=slp )
        ELSE
           CALL TURB_ECMWF ( zt, zu, Ts, XTzt, qs, q_zt, U_zu,   &
              &              Cd, Ch, Ce, XTzu, XQzu, XUblk)
        END IF
        !!
        !!
     CASE DEFAULT
        write(6,*) 'ERROR: mod_aerobulk_compute.f90 => bulk algorithm ', trim(calgo), ' is unknown!!!'
        STOP
     END SELECT


     !! Need the air density at zu m, so using t and q corrected at zu m:
     XRHO = rho_air(XTzu, XQzu, slp)
     QH   = slp - XRHO*grav*zu      ! QH used as temporary array!
     XRHO = rho_air(XTzu, XQzu, QH)


     !! *** Wind stress ***
     Cd_rho_U = Cd*XRHO*XUblk


     !! *** Latent and Sensible heat fluxes ***
     QL = Ce*XRHO*Lvap(Ts)     * (XQzu - qs) * XUblk
     QH = Ch*XRHO*cp_air(XQzu) * (XTzu - Ts) * XUblk

  END SUBROUTINE aerobulk_compute_nextsim

end module mod_aerobulk_compute_nextsim
