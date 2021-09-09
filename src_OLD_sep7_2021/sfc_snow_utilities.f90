!-----------------------------------------------------------------------------
!+ Utility module 'sfc_snow_utilities'
!------------------------------------------------------------------------------

MODULE sfc_snow_utilities

!------------------------------------------------------------------------------
!
! Description:
!  This module contains all utility subroutines used in the snow model.
!
!  Note : For now we have this in one file which is huge/long. It might make
!  sense to split this in more modules. 
!         
!
! Current Code Owner: MeteoSwiss, Sascha Bellaire
!   phone:
!   fax:
!   email:  sascha.bellaire@meteoswiss.com
!
!
! History:
!   Version    Date       Name
!   ---------- ---------- ----
!   1.0       yyy/mm/dd Sascha Bellaire
!   Initial release, i.e. first stable version
!
!
! Code Description:
!   Language: Fortran 90.
!   Software Standards: "European Standards for Writing and
!   Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! =============================================================================
! + Begin loading modules
! =============================================================================

#ifdef __COSMO__
USE kind_parameters, ONLY :   &
#elif __ICON__
USE mo_kind, ONLY:     &
#endif
    vpp           ! KIND-type parameter for real variables (variable precisionphysics)

USE sfc_snow_data

! Physical and other constants
USE data_constants, ONLY: &
    t0_melt_wp => t0_melt

use data_modelconfig, only: ke_snow

!use fxtr_definition, only: pc_std_p0
!use tsa_lmparam, only: qvsat
use data_runcontrol, only: ntstep

! =============================================================================
! - End loading modules
! =============================================================================

IMPLICIT NONE

PUBLIC           ! All constants and variables in this module are public

CONTAINS

!==============================================================================
! + Begine subroutine new_snow_density
! =============================================================================

SUBROUTINE new_snow_density(rho_hn, t_a, uv, in_qv)

  !$acc routine seq

! ------------------------------------------------------------------------------
!
! Description: This one estimates the new snow density based on meteorological
! conditions
!
!
! Method:  Multi-linear regression model
!
!
!                            (REFERENCE)
!
! -----------------------------------------------------------------------------
! ------------------------------------------------------------------------------
! + Begin declarations
! ------------------------------------------------------------------------------

  ! List of INOUT parameter:

   REAL    (KIND = vpp), INTENT(INOUT) :: &

    rho_hn

  ! List of IN parameter:
   REAL    (KIND = vpp), INTENT(IN) :: &

    t_a                                , & ! air temperature             (K)
    uv                                 , & ! wind speed                  (m s-1)
    in_qv

  ! List of local arrays, vectors and scalars

   REAL    (KIND = vpp), PARAMETER ::  &

!    rh     = 0.8        , & ! relative humidity                          (-)

    beta01 = 3.28       , & ! coefficients
    beta1  = 0.03       , & !
    beta02 = -0.36      , &
    beta2  = -0.75      , &
    beta3  = 0.3

   REAL    (KIND = vpp) ::  &

    vw              , & ! wind speed limited to 2 m s-1                  (m/s)
    t_c             , & ! air temperature converted to degrees Celsius   (Celsius)
    arg             , & ! argument for power low
    in_rh,qv_sat

! ------------------------------------------------------------------------------
! - End declarations
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Start subroutine arguments
! ------------------------------------------------------------------------------

  ! Lower boundary of wind speed - 2 m/s
   vw = MAX(2.0_vpp, uv)

  ! Temperature- Kelvin to Celsius
   t_c = t_a - 273.15_vpp

  ! in RH
  ! qv_sat = qvsat (t_a, REAL(pc_std_p0,vpp) ) 

   in_rh =  0.8 !in_qv / qv_sat

  ! ... power law argument
  IF(t_c .LT. -14.0_vpp) THEN

   arg = beta01 + beta1*t_c + beta2*ASIN(SQRT(in_rh)) + beta3*LOG10(vw)

  ELSE

   arg = beta01 + beta1*t_c + beta2*ASIN(SQRT(in_rh)) + beta3*LOG10(vw) + beta02

  ENDIF

  ! Caclualte new snow density
  rho_hn = 10.0_vpp**arg

!  rho_hn = 100.0_vpp

  write(*,*) 'new snow: ',ntstep,vw,in_rh,in_qv,qv_sat,rho_hn

  ! Limit new snow density
  rho_hn = MAX(100.0_vpp, MIN(rho_hn, 250.0_vpp))

! ------------------------------------------------------------------------------
! - End subroutine
! ------------------------------------------------------------------------------

END SUBROUTINE new_snow_density




!==============================================================================
! + Begin subroutine new_snow_layers
! =============================================================================

SUBROUTINE new_snow_layers(top,  hn, dzm_sn_now                  , &        
                           t_sn_now, t                           , &
                           theta_i_now, theta_w_now, theta_a_now , &
                           rho_sn, rho_hn)

  !$acc routine seq

! ------------------------------------------------------------------------------
!
! Description: Assigns properties to new snow layers
!
!
! Method: Nothing specifc see below ;)
!
!
! -----------------------------------------------------------------------------
! ------------------------------------------------------------------------------
! + Begin declarations
! ------------------------------------------------------------------------------


  ! List of INOUT parameter:
  REAL (KIND=vpp),    INTENT(INOUT)    :: &

    dzm_sn_now(:)         , &         ! height of main level
    t_sn_now(:)           , &         ! temperature
    theta_i_now(:)        , &         ! volumetric ice content
    theta_w_now(:)        , &         ! dito for water
    theta_a_now(:)        , &         ! dito for air
    rho_sn(:)                         ! density

  INTEGER,     INTENT(IN)           :: &

    top                               ! current index of top layer

  ! List of IN parameter:
  REAL (KIND=vpp), INTENT(IN)        :: &

    rho_hn                 , &        ! new snow density
    hn                     , &        ! new snow amount
    t                               ! air temperature

  ! List od local integers, vectors and arrays

  INTEGER         :: &

    ksn                                ! loop vector snow layers

  ! physical and other  constants
  REAL    (KIND = vpp) :: &
   
     t0_melt                            ! melting temperature            (K)


! ------------------------------------------------------------------------------
! - End declarations
! ------------------------------------------------------------------------------

  ! Double precision equivalents of module variables
  t0_melt   = REAL(t0_melt_wp  , vpp)

! ------------------------------------------------------------------------------
! + Start subroutine arguments
! ------------------------------------------------------------------------------

 ! ----------------
 ! Assign layer properties
 ! ----------------

!  IF(top .EQ. 1) THEN ! very first snow layer

      dzm_sn_now(top)    = hn                                             ! layer height

      rho_sn(top)        = rho_hn                                         ! snow density

      theta_i_now(top)   = rho_hn / rho_i                                 ! volumetric ice content
      theta_w_now(top)   = 0.0_vpp                                        ! volumetric water content. NOTE: New snow is always dry !!!
      theta_a_now(top)   = 1.0_vpp - theta_i_now(top) - theta_w_now(top)  ! volumetric air content

      t_sn_now(top)      = MIN(t, t0_melt)                     ! new snow temperature equals air temperature, but not larger melting temperature

!  ELSE

!      dzm_sn_now(top)   =  hn

!      rho_sn(top)       = rho_hn

!      theta_i_now(top) = rho_hn / rho_i
!      theta_w_now(top) = 0.0_vpp                                          ! new snow is always dry
!      theta_a_now(top) = 1.0_vpp - theta_i_now(top) - theta_w_now(top)

!      t_sn_now(top)    = MIN(t, t0_melt)

!  ENDIF

 ! ----------------
 ! Calculate a few things
 ! ----------------

!  ! Layer depth from layer height
!  DO ksn = top, 1, -1

!    zm_sn((top+1) - ksn) = hm_sn(ksn)

!  ENDDO

! ------------------------------------------------------------------------------
! - End subroutine
! ------------------------------------------------------------------------------

END SUBROUTINE new_snow_layers








!==============================================================================
! + Begin subroutine update
! =============================================================================

SUBROUTINE update(top, hm_sn, zm_sn, dzm_sn_now, m_sn      , &
                  theta_i_now, theta_w_now, theta_a_now    , &
                  rho_sn, hcap_sn, hcon_sn)
                   

  !$acc routine seq

! ------------------------------------------------------------------------------
!
! Description: Updates depending fields.
!
!
! Method: Nothing specifc see below ;)
!
!
! -----------------------------------------------------------------------------
! ------------------------------------------------------------------------------
! + Begin declarations
! ------------------------------------------------------------------------------

  ! List of INOUT parameter:
  REAL (KIND=vpp),    INTENT(INOUT)    :: &

    hm_sn(:)              , &         ! height of main level
    zm_sn(:)              , &         ! depth of main level
    rho_sn(:)             , &         ! density
    hcon_sn(:)            , &         ! heat conductivity
    hcap_sn(:)            , &         ! heat capacity
    m_sn(:)               , &         ! snow layer mass
    theta_a_now(:)                    ! volumetric air content


  ! List of IN parameter:
  REAL (KIND=vpp),    INTENT(IN)    :: &
    theta_i_now(:)        , &         ! volumetric ice content
    theta_w_now(:)        , &         ! dito for water
 
    dzm_sn_now(:)                     ! layer thickness


  INTEGER,     INTENT(IN)           :: &

    top                               ! current index of top layer! List of INOUT parameter:

  ! List of local scalars and array's
  INTEGER         :: &

    ksn                                ! loop vector snow layers


! ------------------------------------------------------------------------------
! - End declarations
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Begin subroutine arguments
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Update fields
! ------------------------------------------------------------------------------

  ! -------------------------
  ! Height of snow (main) levels
  ! -------------------------

   DO ksn = 1, top

     IF(ksn .EQ. 1) THEN

       hm_sn(ksn) = dzm_sn_now(ksn)

      ELSE

        hm_sn(ksn) = hm_sn(ksn-1) + dzm_sn_now(ksn)

      ENDIF

   ENDDO


  ! --------------------------
  ! Depth of snow layer (main) levels
  ! -------------------------

   DO ksn = top, 1, -1

     zm_sn( (top+1) - ksn) = hm_sn(ksn) !invert height vector
   
   ENDDO


  ! --------------------------
  ! Volumetric air content
  ! -------------------------

   DO ksn = top, 1, -1

     theta_a_now(ksn) = MAX(0.0_vpp, 1.0_vpp - theta_i_now(ksn) - theta_w_now(ksn))

   ENDDO

  ! --------------------------
  ! Snow layer density
  ! -------------------------

   DO ksn = 1, top, 1

        IF(theta_i_now(ksn) .EQ. 0.0_vpp) THEN

          rho_sn(ksn) = 0.0_vpp

        ELSE

          rho_sn(ksn) = theta_i_now(ksn)*rho_i + theta_w_now(ksn)*rho_w

        ENDIF

   ENDDO


  ! --------------------------
  ! Heat capacity
  ! --------------------------

   DO ksn = 1, top, 1

       IF(rho_sn(ksn) .LT. eps) THEN

         hcap_sn(ksn) = 0.0_vpp

       ELSE

         hcap_sn(ksn) = (  rho_a   * theta_a_now(ksn) * specific_heat_air     &
                         + rho_i   * theta_i_now(ksn) * specific_heat_ice     &
                         + rho_w   * theta_w_now(ksn) * specific_heat_water)  &
                         / rho_sn(ksn)

       ENDIF

   ENDDO


  ! --------------------------
  ! Heat conductivity
  ! --------------------------

   DO ksn = 1, top, 1

      hcon_sn(ksn) = 2.22_vpp * EXP(1.88_vpp * LOG(rho_sn(ksn)/rho_i))

   ENDDO


  ! --------------------------
  ! Snow layer mass
  ! --------------------------

   DO ksn = 1, top, 1

     IF(ksn .EQ. 1) THEN

        m_sn(ksn) = zm_sn(ksn) * ((theta_i_now(ksn) * rho_i) + (theta_w_now(ksn) * rho_w))

      ELSE

        m_sn(ksn) = ABS(zm_sn(ksn) - zm_sn(ksn-1))  * ((theta_i_now(ksn) * rho_i) + (theta_w_now(ksn) * rho_w))

      ENDIF

   ENDDO



! ------------------------------------------------------------------------------
! - End subroutine
! ------------------------------------------------------------------------------

END SUBROUTINE update




!==============================================================================
! + Begin subroutine calc_tch
! =============================================================================

SUBROUTINE calc_tch(t0, t1, tch, uv, q1, h_snow, p)


  !$acc routine seq

! ----------------------------------------------------------------------------
!
! Description: Calculates the turbulent transfer coefficient
!
!
! Method: Schlögl, S., Lehning, M., Nishimura, K., Huwald, H., Cullen, N. J.,  and Mott, R., 2017.
!                 How do stability corrections perform in the stable boundary  layer over snow?,
!                 Bound.-Layer Meteor., https://doi.org/10.1007/s10546-017-0262-1.
!
! -----------------------------------------------------------------------------
! ------------------------------------------------------------------------------
! + Begin declarations
! ------------------------------------------------------------------------------

  ! List of INOUT parameter:
  REAL (KIND=vpp), INTENT(INOUT)        :: &

    tch                     ! transfer coefficient for heat


  ! List of IN parameter:
  REAL (KIND=vpp), INTENT(IN)        :: &

    t0                  , & ! optimized snow surface temperature
    uv                  , & ! wind velcocity
    t1                  , & ! first level temperature (~ 10m)
    q1                  , & ! first level specific humidity
    h_snow              , & ! snow height
    p                       ! surface pressure


  ! Local arrays, vectors and scalars
  REAL    (KIND = vpp) ::  &

    z_ref                 , & ! reference height of meteo values
    z0                    , & ! roughness length

    phi_s                 , & ! virtuell temperature at ground surface
    phi_ref               , & ! virtuell temperature at first level

    psi_m                 , & ! stability correction for momentum
    psi_s                 , & ! stability correction for scalars

    u_star                , & ! friction velocity
    t_star                , & ! temperature scale

    zeta                  , & ! stability paramter

    B                     , & ! buoyancy contribution
    S                     , & ! shear contribution

    E_s                   , & ! saturation water pressure
    e_v                   , & ! water vapour pressure

    q0                        ! specific humidity at surface


  REAL    (KIND = vpp), PARAMETER ::  &

    a1         = -65.35_vpp         , & ! coefficients
    b1         = 0.0017_vpp         , &

    a2         = -813.21_vpp        , &
    b2         = -0.0014_vpp        , &

    c1         = 0.0_vpp            , &
    c2         = 0.0_vpp

! ------------------------------------------------------------------------------
! - End declarations
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Begin subroutine arguments
! ------------------------------------------------------------------------------

! -------------------------
! + Some pre-calculations
! -------------------------

 ! Assign of few parameters
   z_ref    = height_meteo_values !- h_snow
   z0       = roughness_length

 ! Separation of bulk Richardson in
   B = (t0-t1) / (0.5_vpp * (t1+t0))      ! buoyancy
   S = (z_ref*g) /uv**2                      ! shear contribution

 ! Specific humidity at snow surface (asume saturation)
   E_s = 6.112_vpp * EXP( (22.46_vpp * (t0-273.15_vpp)) / (272.62_vpp + t0) )  ! Saturation vapour pressure via Magnus Equation
   e_v = 100.0_vpp * E_s                                          ! Water vapour pressure from relative humidty
   q0  = 0.622 * (e_v/p)                                          ! Specific humidity

 ! Stability correction
   psi_m = a1*B + b1*S + c1   ! for momentum
   psi_s = a2*B + b2*S + c2   ! and scalar

 ! Virtuell temperature at ground (snow) surface
   phi_s   = t0   * (1.0_vpp + rvd_m_o * q0)

 ! Virtuel temperature first level
   phi_ref = t1   * (1.0_vpp + rvd_m_o * q1)

 ! Temperature scale
   t_star  = karman * (phi_s - phi_ref) / (LOG( z_ref / z0) - psi_s)

 ! Stability paramter
   zeta    = (-karman * z_ref * g * t_star) / (phi_s * uv**2)

! -------------------------
! + Calculate transfer coefficient
! -------------------------

   tch = karman**2 / ( (LOG(z_ref/z0) - psi_m*zeta) * (LOG(z_ref/z0) - psi_s*zeta) )

 ! NEUTRAL
 !   tch = karman**2 /  ( LOG(z_ref/z0) * LOG(z_ref/z0) )

! ------------------------------------------------------------------------------
! - End subroutine
! ------------------------------------------------------------------------------

END SUBROUTINE calc_tch




!==============================================================================
! + Begin subroutine turb_fluxes
! =============================================================================

SUBROUTINE turb_flux(sh, lh, p, t1, t0, q1, tch, uv)


  !$acc routine seq

! ----------------------------------------------------------------------------
!
! Description: Calculates the turbulent fluxes latent and sensible
!
!
! Method: Bulk method with stability function by Schloegl et al.
!
! -----------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Begin declarations
! ------------------------------------------------------------------------------

  ! List of INTENT(INOUT) parameters
  REAL (KIND = vpp), INTENT(INOUT) :: &

    sh                   , & ! sensible heat flux
    lh                       ! latent heat flux

  ! List of INTENT(IN) parameters
  REAL (KIND = vpp), INTENT(IN)  :: &
  
    p                    , & ! surface pressure        [PA]
    t0                   , & ! surface temperature     [K]
    t1                   , & ! air temperature         [K]
    q1                   , & ! specific humidty air    [kg/kg]
    tch                  , & ! transfer coefficient    []
    uv                       ! wind speed              [m/s]


  ! Local arrays, vectors and scalars
  REAL    (KIND = vpp) ::  &

    E_s                   , & ! saturation water pressure
    e_v                   , & ! water vapour pressure

    q0                    , & ! specific humidity at surface

    t_v                   , & ! virtuell temperature  

    rho_atm                   ! density of atmosphere


! ------------------------------------------------------------------------------
! - End declarations
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Begin subroutine arguments
! ------------------------------------------------------------------------------

 ! -------------------------
 ! + Some pre-calculations
 ! -------------------------

   ! Specific humidity at surface assuming saturation
   E_s = 6.112_vpp * EXP( (22.46_vpp * (t0-273.15_vpp)) / (272.62_vpp + t0) )  ! Saturation vapour pressure via Magnus Equation
   e_v = 100.0_vpp * E_s                                          ! Water vapour pressure from relative humidty
   q0  = 0.622 * (e_v/p)

   ! Virtuell temperature 
   t_v = t1     * (1.0_vpp + rvd_m_o * q1)

   ! Density of atmosphere
   rho_atm = 0.9_vpp !p / ( gas_constant_air * t_v)

 ! -------------------------
 ! + Calculate turbulent fluxes
 ! -------------------------

   ! Latent heat flux
   lh = tch * uv * rho_atm * lh_v * (q1 - q0)
   lh = MIN(150.0_vpp,MAX(-150.0_vpp,lh))


   ! Sensible heat flux
   sh = tch * uv * rho_atm * specific_heat_air * (t1 - t0)
   sh = MIN(250.0_vpp,MAX(-250.0_vpp,sh))


!   write(*,*) 'turb_flux: ',tch,uv,rho_atm,t1,t0,q1,q0,lh,sh
! ==============================================================================
! - End subroutine
! ==============================================================================

END SUBROUTINE turb_flux





!==============================================================================
! + Begin subroutine rad_flux
! =============================================================================

SUBROUTINE rad_flux(swdir_s, swdifd_s, swdifu_s, lwd_s, lwu_s, &
                    alpha_sn, t, t_sn_sfc)


  !$acc routine seq

! ----------------------------------------------------------------------------
!
! Description: Calculates radiative fluxes
!
!
! Method: See below for details ... ;)
!
! -----------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Begin declarations
! ------------------------------------------------------------------------------

  ! List of INTENT(INOUT) parameters
  REAL (KIND = vpp), INTENT(INOUT) :: &
    swdifu_s            , &    ! upwelling diffuse shortwave radiation, i.e. reflected  [W/m**2]
    alpha_sn            , &    ! albedo snow surface
    lwu_s                      ! long-wave radiation  - upward              [w/m**2]

  ! List of INTENT(IN) paramters
  REAL (KIND = vpp), INTENT(IN)  :: &
    swdir_s             , &    ! short-wave radiation - direct down         [W/m**2]
    swdifd_s            , &    ! short-wave radiation - diffuse down        [W/m**2]
    lwd_s               , &    ! long-wave radiation  - down                [W/m**2]
    t                   , &    ! air temoerature                            [K]
    t_sn_sfc                   ! snow surface temperature                   [K]

 REAL (KIND=vpp), PARAMETER            :: &

    a1 = 0.9_vpp        , &    ! coefficients for the albedo parameterization
    a0 = 0.6_vpp

! ------------------------------------------------------------------------------
! - End declarations
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Begin subroutine arguments
! ------------------------------------------------------------------------------

! -------------------------
! Calculate albedo
! -------------------------

  ! UKMO GCM
  IF(t_sn_sfc - 273.15_vpp .LE. -2.0_vpp) THEN

    alpha_sn = a1

  ELSE

    alpha_sn  = a1 - 0.15_vpp*(a1-a0)*((MIN(t_sn_sfc, t0_melt) - 273.15_vpp) + 2.0_vpp)

  ENDIF

! -------------------------
! Calculate upward short-wave radiation
! -------------------------

  swdifu_s = swdir_s * alpha_sn

! -------------------------
! Calculate upward long-wave radiation
! -------------------------

  lwu_s = sigma * (1.0_vpp - snow_Ctalb) * t_sn_sfc**4

! ==============================================================================
! - End subroutine
! ==============================================================================

END SUBROUTINE rad_flux




!==============================================================================
! + Begin subroutine rad_flux
! =============================================================================

SUBROUTINE solve_1d_heat(dzm_sn, t_sn, t_sn_sfc                 , &
                         for_sn, swabs_sn, zm_sn                , &
                         hcap_sn, hcon_sn, hdif_sn, rho_sn      , &
                         top, t_so_now, zmls, zroc, zalam       , &
                         ke_soil, dt               )


  !$acc routine seq

! ----------------------------------------------------------------------------
!
! Description: This is the implicit solver for the 1D heat equation
!
! Method: First combine snow an soil parameters, second set up a tri-diagonal
!         matrix, third solve using Thomson Alrgorithm.
!
! -----------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Begin declarations
! ------------------------------------------------------------------------------

  ! INTENT(INOUT)
  REAL (KIND=vpp), INTENT(INOUT)     ::  &

    t_sn_sfc                  ! snow surface temperature          (K)


   REAL (KIND=vpp), INTENT(INOUT)   ::   &

    dzm_sn(:)            , &  ! layer thickness

    t_sn(:)              , &  ! layer temperature                 (K)
    t_so_now(:)          , &

    hdif_sn(:)                ! heat diffusivity

  ! INTENT(IN)
   REAL (KIND=vpp), INTENT(IN)   ::   &

    rho_sn(:)             , & ! layer density
    hcon_sn(:)            , & ! heat conductivity
    hcap_sn(:)            , & ! heat capacity
    zm_sn(:)              , & ! main level depth (snow)
    dt                        ! integration timestep

   INTEGER   , INTENT(IN)      ::  &

    top                   , & ! top layer index
    ke_soil                   ! number of soil layers. NOTE without surface
                              ! layer and climatologic layer
   REAL (KIND=vpp), INTENT(IN)   ::   &

    zmls(:)               , &
    zalam(:)              , &
    zroc(:)

   REAL (KIND=vpp), INTENT(IN)        ::  &

    for_sn               , &    ! atmosperic forcing                                                  (W m-2)
    swabs_sn(:)                ! absorbed short-wave radiation in each layer                         (W m-2)

  ! Local arrays, scalar and/or vectors

  INTEGER              ::  &

    ksn                  , &     ! loop index for snow layers
    counter                      ! counter to help with indexing


  REAL (KIND=vpp)       ::  &

    dt_sub                , &  ! sub time step
    dlw_u_sn              , &  ! derivative of upwelling longwave radiation for snow on the ground    (W m-2)
    dz_up                 , &  ! thickness above the layer of interest                                (m)
    dz_low                , &  ! thickness below the layer of interest                                (m)
    beta

  REAL (KIND=vpp), DIMENSION(-top+1:ke_soil+1) :: &

    zm                   , &  ! depth of main levels
    zh                   , &  ! depth of half levels

    dz                   , &  ! layer thickness

    hcon                 , & ! heat conductivity
    hcap                 , & ! heat capacity
    hdif                 , & ! heat diffusivity

    rho                  , &  ! density
    tmp                  , &

    t                    , &  ! temperature of the snow/substrate column

    sw_abs               , &  ! absorbed short-wave radiaton in each layer

    alpha                , &  ! utility variables for building and solving the tri-diagonal matrix
    gamma                , &  !

    a                    , &  !
    b                    , &  !
    c                    , &  !
    d                    , &  !
    e                         ! final snow layer temperature


! ------------------------------------------------------------------------------
! - End declarations
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Begin subroutine arguments
! ------------------------------------------------------------------------------

  IF(top .GE. 1) THEN ! enough snow/layers on the ground

    ! ---------------------------------------------
    ! Assign snow layer properties to local vectors
    ! ----------------------------------------------

    ! Assign default values
    ! ------------------------

      zm(:)   = 0.0_vpp

      hcon(:) = 0.0_vpp
      hcap(:) = 0.0_vpp

      rho(:)  = 0.0_vpp

      t(:)    = 0.0_vpp

      sw_abs  = 0.0_vpp

      counter = 1

   ! Sort layers
   ! ---------------------------

   DO ksn = -top+1, ke_soil+1, 1

     IF(ksn .LE. 0) THEN  ! snow layers

         IF(ksn .EQ. -top+1) THEN

           zm(ksn)     = zm_sn(top)

           hcon(ksn)   = hcon_sn(top)
           hcap(ksn)   = hcap_sn(top)

!           rho(ksn)    = rho_sn(top)

           t(ksn)      = t_sn(top)

           sw_abs(ksn) = swabs_sn(top)


         ELSE

           zm(ksn)     = zm_sn(top-counter)

           hcon(ksn)   = hcon_sn(top-counter)
           hcap(ksn)   = hcap_sn(top-counter)

!           rho(ksn)    = rho_sn(top-counter)

           t(ksn)      = t_sn(top-counter)

           sw_abs(ksn) = swabs_sn(top-counter)

           counter = counter + 1

         ENDIF

     ELSE  ! soil layers

           zm(ksn)     = zm_sn(1) + zmls(ksn)

           hcon(ksn)   = zalam(ksn)
           hcap(ksn)   = zroc(ksn)

!           rho(ksn)    = rho_sn(1)

           t(ksn)      = t_so_now(ksn)

           sw_abs(ksn) = 0.0_vpp

     ENDIF

  ENDDO ! end of snow layers


  ! ------------------------------------------------------------
  ! Some precalculations ...
  ! ------------------------------------------------------------

  ! Derivative of emitted long wave radiation
  ! ----------------------

    dlw_u_sn =  -4.0_vpp * sigma * (1.0_vpp - snow_Ctalb) * t(-top+1)**3


  ! Calculate factors (diffusion) for the linear equations for ...
  ! ----------------------


    ! Loop over all layers ...
    DO ksn = -top+1, ke_soil, 1

      IF(ksn .LE. ke_soil-1) THEN ! ... except

        alpha(ksn)   = dt / hcap(ksn)

        hdif(ksn) = hcon(ksn) * (t(ksn+1) - t(ksn)) /  (zm(ksn+1) - zm(ksn))

      ELSE !(ksn == substrate) ! ... bottom layer

         alpha(ksn) = dt / hcap(ksn)

         hdif(ksn) = 0.0_vpp

      ENDIF
  
  END DO

  ! ------------------------------------------------------------
  ! Setup tridiagonal matrix for set of linear equations for each layer ...
  ! ------------------------------------------------------------

  DO ksn = -top+1, ke_soil


    IF(ksn .EQ. -top+1) THEN ! ... TOP LAYER

    dz_low = zm(ksn+1) - zm(ksn)

      a(ksn)     = 0.0_vpp
      b(ksn)     = 1 + (1.0_vpp - cn) * alpha(ksn) * hcon(ksn)/dz_low - alpha(ksn) * dlw_u_sn
      c(ksn)     = -   (1.0_vpp - cn) * alpha(ksn) * hcon(ksn)/dz_low

      d(ksn)     = t(ksn) + alpha(ksn) * (for_sn - dlw_u_sn*t(ksn) + cn*hdif(ksn))

    ELSEIF (ksn .LE. ke_soil-1) THEN ! ... INNER LAYERS

     dz_up  = zm(ksn)   - zm(ksn-1)
     dz_low = zm(ksn+1) - zm(ksn)

        a(ksn) = -         (1.0_vpp - cn)   * alpha(ksn) *  hcon(ksn-1)/dz_up
        b(ksn) = 1.0_vpp + (1.0_vpp - cn)   * alpha(ksn) * (hcon(ksn)  /dz_low + hcon(ksn-1)/dz_up)
        c(ksn) = -         (1.0_vpp - cn)   * alpha(ksn) *  hcon(ksn)  /dz_low

        d(ksn) = t(ksn) + cn*alpha(ksn) * (hdif(ksn) - hdif(ksn-1)) + alpha(ksn)*sw_abs(ksn)


    ELSEIF (ksn .EQ. ke_soil) THEN ! BOTTOM LAYERS

     dz_up = zm(ksn)   - zm(ksn-1)

       a(ksn) = -         (1.0_vpp - cn) * alpha(ksn) * hcon(ksn-1)/dz_up
       b(ksn) = 1.0_vpp + (1.0_vpp - cn) * alpha(ksn) * hcon(ksn-1)/dz_up
       c(ksn) = 0.0_vpp

       !d(bot_idx) = t(bot_idx) - cn*alpha(bot_idx-1) + alpha(bot_idx)*hdif(bot_idx)
       d(ksn)     = t(ksn) - cn*alpha(ksn)*hdif(ksn-1) + alpha(ksn)*hdif(ksn)

    ENDIF



  ENDDO

  ! ------------------------------------------------------------
  ! Solve the system - Thomas Algorithm
  ! ------------------------------------------------------------

   beta = b(-top+1)

    ! Forward substitution

    DO ksn = -top+1, ke_soil, 1

      IF(ksn .GE. -top+1) THEN

        IF(ksn .EQ. -top+1) THEN

          e(ksn) = d(ksn) / beta

        ELSE

          gamma(ksn) = c(ksn-1) / beta
          beta       = b(ksn) - a(ksn) * gamma(ksn)
          e(ksn)     = (d(ksn) - a(ksn) * e(ksn-1)) / beta

        ENDIF

      ENDIF

    ENDDO

    ! Backward substitution

    DO ksn = ke_soil-1, -top+1, -1

      IF(ksn .GE. -top+1) THEN

        e(ksn) = e(ksn) - gamma(ksn+1) * e(ksn+1)

      ENDIF

    ENDDO

  ! ------------------------------------------------------------
  ! Do some updating required for the next sections
  ! ------------------------------------------------------------

   ! Snow Surface Temperature
    t_sn_sfc = e(-top+1)

    ! Snow layer temperature

    counter = 1

    DO ksn = top,1,-1

      IF(ksn .EQ. top) THEN

       t_sn(ksn)    = e(-top+1)

      ELSE

       t_sn(ksn)    = e(-top+1+counter)

       counter = counter + 1

      ENDIF

    ENDDO


   ! Soil layer temperature ; without 'surface' and climatologic layer





 ELSE ! only a single layer

       ! for now exit code

 ENDIF




! ==============================================================================
! - End subroutine
! ==============================================================================

END SUBROUTINE solve_1d_heat




!==============================================================================
! + Begin subroutine init_phase change
! =============================================================================

SUBROUTINE init_phase_change(t_sn, theta_i, theta_w, top  , & 
                             melt_flag, freeze_flag)

  !$acc routine seq

! ----------------------------------------------------------------------------
!
! Description: Does some initializations for phase change calculations
!
! Method: See below for details ...
!
! -----------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Begin declarations
! ------------------------------------------------------------------------------

  LOGICAL, INTENT(INOUT)                     :: &

   melt_flag             , &
   freeze_flag

  REAL (KIND=vpp), INTENT(IN)                     :: &

   t_sn(:)                  , &
   theta_i(:)               , &
   theta_w(:)

  INTEGER, INTENT(IN)                     :: &

    top     

  INTEGER           :: &

    ksn

! ------------------------------------------------------------------------------
! - End declarations
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Start subroutine arguments
! ------------------------------------------------------------------------------

   !--------------------
   !Initiate melt/freeze flags
   !--------------------

    ! Initiate
    melt_flag   = .FALSE.
    freeze_flag = .FALSE.

    DO ksn = top, 1, -1

     ! Metling
      IF(t_sn(ksn) .GT. t0_melt .AND. theta_i(ksn) .GT. eps2) THEN

        melt_flag = .TRUE.

      ENDIF

     ! Freezing
      IF(t_sn(ksn) .LT. t0_melt .AND. theta_w(ksn) .GT. (theta_r + eps2)) THEN

        freeze_flag = .TRUE.

      ENDIF

    ENDDO  !end of snow layers

! ==============================================================================
! - End subroutine
! ==============================================================================

END SUBROUTINE init_phase_change




!==============================================================================
! + Begin subroutine melt_snow
! =============================================================================

SUBROUTINE melt_snow(zm_sn, dzm_sn, t_sn, hcap_sn,      &
                     theta_i, theta_w, rho_sn,          &
                     top)

  !$acc routine seq

! ----------------------------------------------------------------------------
!
! Description: Check's for melt potential and adjusts properties accordingly
!
! Method: See below for details ...
!
! -----------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Begin declarations
! ------------------------------------------------------------------------------

  ! List of INOUT parameter:
   REAL    (KIND = vpp),  INTENT(INOUT) :: &

    t_sn(:)           , & ! snow temperature

    theta_i(:)        , & ! new values of volumetric ice content
    theta_w(:)            ! new values of volumetric water content


   ! List of INOUT parameter:
   REAL    (KIND = vpp),  INTENT(IN) :: &

    zm_sn(:)          , & ! main level depth
    dzm_sn(:)         , & ! main level height

    rho_sn(:)         , & ! snow layer density                          (kg/m3)

    hcap_sn(:)               ! heat capacity


  ! List of IN parameter:
   INTEGER,   INTENT(IN)    :: &

    top


  ! List of local scalar and paramters
  REAL    (KIND = vpp) ::  &

    dT             , & ! difference between current temperature and melting temperature
    A              , & ! coefficient A (see notes below)
    dtheta_i       , & ! change in volumetric ice content
    dtheta_w       , & ! change in volumetric water content

    q_mf           , & !
    q_rest


  INTEGER  :: &

    ksn

! ------------------------------------------------------------------------------
! - End declarations
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
! + Start subroutine arguments
! ------------------------------------------------------------------------------

 DO ksn = top, 1, -1

    ! --------------------
    ! Now see if any melting is going on -- this implies that (1) the
    ! temperature of the element
    ! is above/equal the melting temperature (2) there is something to melt and
    ! (3) there is enough
    ! room to place the meltwater ...
    ! --------------------

    IF(t_sn(ksn) .GE. t0_melt .AND. theta_i(ksn) .GT. 0.0_vpp .AND. theta_w(ksn) .LT. theta_s) THEN

     ! difference dT between actual melting temperature and layer temperature
     dT = 0.0_vpp
     dT = t0_melt - t_sn(ksn)

    ! --------------------
    ! Now we take into account that there might be some extra energy that could
    ! not
    ! be used by the element above because of complete melting
    ! --------------------

      dT = dT - (q_rest / (hcap_sn(ksn) * rho_sn(ksn) * dzm_sn(ksn)))

      ! ONly do it when there is real potential to melt
      IF(dT .LT. 0.0_vpp) THEN

        ! --------------------
        ! Determine the DECREASE in ice content and the INCREASE of water
        ! content
        ! Adapt A to compute mass changes
        ! --------------------

         A = (hcap_sn(ksn) * rho_sn(ksn)) / (rho_i * lh_f)

         dtheta_i = A * dT
         dtheta_w = - (rho_i / rho_w) * dtheta_i

        ! --------------------
        ! It could happen that there is enough energy available to melt more ice
        ! than is present.
        ! You can only melt so much ice as is there ....
        ! --------------------

         IF( (theta_i(ksn) + dtheta_i) .LT. 0.0_vpp) THEN

           dtheta_i = - theta_i(ksn)
           dtheta_w = - (rho_i / rho_w) * dtheta_i

           dT = dtheta_i / A

         ENDIF

       ! --------------------
        ! It could also be that you are trying to produce more water than is
        ! allowed.
        ! --------------------

          IF( (theta_w(ksn) + dtheta_w) .GT. theta_s) THEN

            dtheta_w = theta_s - theta_w(ksn)
            dtheta_i = - (rho_w/rho_i) *dtheta_w

            dT = dtheta_i / A

          ENDIF

        ! --------------------
        ! Reset/Recalculate properties
        ! --------------------

         ! Layer temperature and transfered energy
         t_sn(ksn) = t_sn(ksn) + dT

         IF(t_sn(ksn) .LE. t0_melt) THEN ! if melting occured it can only be at melt point

           q_rest     = 0.0_vpp
           t_sn(ksn)  = t0_melt

         ELSE

           q_rest = hcap_sn(ksn) * rho_sn(ksn) * dzm_sn(ksn) * (t_sn(ksn) - t0_melt)
           t_sn(ksn)   = t0_melt

         ENDIF

         ! Volumetric freezing power
         q_mf     = q_mf + ((dtheta_i * rho_i * lh_f) / dt)
         dtheta_w = dtheta_w

         ! Contents of ice, water and air
         theta_i(ksn) = theta_i(ksn) + dtheta_i
         theta_w(ksn) = theta_w(ksn) + dtheta_w


      ENDIF ! deltaT check


     ENDIF ! end melt check


  ENDDO ! loop over snow layers

! ==============================================================================
! - End subroutine
! ==============================================================================

END SUBROUTINE melt_snow




!==============================================================================
! + Begin subroutine melt_snow
! =============================================================================

SUBROUTINE freeze_snow(zm_sn, dzm_sn, t_sn, hcap_sn,         &
                       theta_i, theta_w, theta_a, rho_sn,    &
                       top)

  !$acc routine seq

! ----------------------------------------------------------------------------
!
! Description: Does re-freeze the water and adjusts  properties accordingly
!
! Method: See below for details ...
!
! -----------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Begin declarations
! ------------------------------------------------------------------------------

 ! List of INOUT parameter:
   REAL    (KIND = vpp),  INTENT(INOUT) :: &

    t_sn(:)           , & ! snow temperature

    theta_i(:)        , & ! new values of volumetric ice content
    theta_w(:)        , & ! new values of volumetric water content
    theta_a(:)            ! new values of volumetric air content

   ! List of INOUT parameter:
   REAL    (KIND = vpp),  INTENT(IN) :: &

    zm_sn(:)          , & ! main level depth
    dzm_sn(:)         , & ! main level height

    rho_sn(:)         , & ! snow layer density                          (kg/m3)

    hcap_sn(:)               ! heat capacity

  ! List of IN parameter:
   INTEGER,   INTENT(IN)    :: &

    top

  ! List of local scalar and paramters
  REAL    (KIND = vpp) ::  &

    dT             , & ! difference between current temperature and melting temperature
    A              , & ! coefficient A (see notes below)
    dtheta_i       , & ! change in volumetric ice content
    dtheta_w       , & ! change in volumetric water content

    q_mf           , & !
    q_rest


  INTEGER  :: &

    ksn

! ------------------------------------------------------------------------------
! - End declarations
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
! + Start subroutine arguments
! ------------------------------------------------------------------------------


  DO  ksn = top, 1, -1

    ! Freezing within the snowpack can occur if (1) the temperature of the
    ! element is below freezing and if water is present to be refrozen
    IF(t_sn(ksn) .LT. t0_melt .AND. theta_w(ksn) .GT. theta_r) THEN

      ! difference dT between actual layer temperature and freezing temperature
      dT = 0.0_vpp
      dT = t0_melt - t_sn(ksn)

      ! Adapt A to compute mass change
      A = (hcap_sn(ksn) * rho_sn(ksn)) / (rho_i * lh_f)

      ! Compute change in volumetric contenst
      dtheta_i = A * dT
      dtheta_w = -(rho_i/rho_w) * dtheta_i

      ! Make sure that there is enough water to refreeze
      IF( (theta_w(ksn) + dtheta_w) .LT. theta_r) THEN

        dtheta_w = -ABS(theta_w(ksn) - theta_r)
        dtheta_i = -(rho_w / rho_i) * dtheta_w
        dT       = dtheta_i / A

      ENDIF

        ! See if the layer is pure ice
        IF( (theta_i(ksn) + theta_r + dtheta_i) .GE. 1.0_vpp) THEN

          dtheta_w = - ABS(theta_w(ksn) - theta_r)
          dtheta_i = - (rho_w/rho_i) * dtheta_w

          theta_i(ksn)  = 1.0_vpp
          theta_w(ksn)  = theta_r
          theta_a(ksn)  = 0.0_vpp

        ELSE

          theta_i(ksn) = theta_i(ksn) + dtheta_i
          theta_w(ksn) = theta_w(ksn) + dtheta_w
          theta_a(ksn) = MAX(0.0_vpp, 1.0_vpp - theta_i(ksn) - theta_w(ksn))

        ENDIF

    ! --------------------
    ! Reset/Recalculate properties
    ! --------------------

       ! Set some limits
       IF(theta_w(ksn) .GE. 1.0_vpp) THEN

         theta_w(ksn) = 1.0_vpp

       ENDIF

       !Compute the volumetric refreezing power
       q_mf     = q_mf + ((dtheta_i * rho_i * lh_f) / dt)
       dtheta_w = dtheta_w

       t_sn(ksn)     = t_sn(ksn) + dT

    ENDIF

  ENDDO



! ==============================================================================
! - End subroutine
! ==============================================================================

END SUBROUTINE freeze_snow



!==============================================================================
! + Begin subroutine phase_change_surface
! =============================================================================

SUBROUTINE phase_change_surface(top, lh, t_sn_sfc, rho_sn, dt, &
                                dzm_sn, theta_i, theta_w )

  !$acc routine seq

! ----------------------------------------------------------------------------
!
! Description: This subroutine computes phase changes of matter (solid to gas,
!              gas to solid and gas to liquid) and adjusts layer properties
!              accoridngly)
!
! Method: If there is snow and laten heat flux is directed towards the surface
!         add mass and update volumetric contents accoridngly. Otherwise
!         substract mass.
!
! -----------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Begin declarations
! ------------------------------------------------------------------------------

  ! List of INOUT parameter:
  REAL    (KIND = vpp), INTENT(INOUT) :: &

    lh                         , & ! latent heat flux (positive towards surface
    dzm_sn(:)                  , & ! layer thickness
    theta_i(:)                 , & ! volumetric ice content
    theta_w(:)                     ! dito for water

  ! List of IN paramter
  REAL  (KIND = vpp), INTENT(IN)   :: &

    t_sn_sfc                   , & ! snow surface temperature
    rho_sn(:)                  , & ! layer density
    dt                             ! integration time step

  INTEGER           , INTENT(IN) :: &

    top             ! index of top snow layer

  ! List of local scalars and arrays

  REAL (KIND = vpp)        :: &

    dL                         , & ! change of layer thickness
    dM                         , & ! change of layer mass
    M                          , & ! initial mass and volmetric content (water of ice)
    hoar                       , & ! hoar mass

    dzm_sn_old                     ! old value of layer thickness

  INTEGER     :: &

    ksn

! ------------------------------------------------------------------------------
! - End declarations
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Start subroutine arguments
! ------------------------------------------------------------------------------

  ! Initiate some values
   dL   = 0.0_vpp
   dM   = 0.0_vpp
   M    = 0.0_vpp
   hoar = 0.0_vpp

  ! Latent heat flux towards the surface - mass gain
  ! ---------------------

  IF(lh .GT. eps2) THEN ! add mass


    IF(t_sn_sfc .LT. t0_melt) THEN  ! add ice

      dM = lh * dt/lh_s     ! calculate mass change
      lh = 0.0_vpp           ! reset latent heat flux, i.e. energy was used

      hoar = dM

      ! Adjust layer properties accordingly, keep snow density constant

      dzm_sn_old = dzm_sn(top)

      dL = dM/rho_sn(top)

      dzm_sn(top) = dzm_sn(top) + dL

      theta_i(top) = theta_i(top) * (dzm_sn_old/dzm_sn(top))
      theta_i(top) = theta_i(top) + (dM/rho_i*dzm_sn(top))

      theta_w(top) = theta_w(top) * (dzm_sn_old/dzm_sn(top))


    ELSE ! add water

      dM = lh *dt/lh_v      ! calculate mass change
      lh = 0.0_vpp           ! reset latent heat, i.e. energy was used

      theta_w(top) = theta_w(top) + dM/(rho_w*dzm_sn(top))   ! update volumetric water content

    ENDIF
 ! Latent heat flux away from the surface - mass loss
  ! ---------------------
  ELSE
   IF(lh .LT. (-1.0_vpp*eps2)) THEN ! additional check in case lh ist super small, but sligtly positive
     ! Not sure how acc statements need to look for this construct
     DO WHILE (lh .LT. (-1.0_vpp*eps2)) ! while energy is available
       DO ksn = top, 1, -1 ! loop through snow layers
         IF(theta_w(ksn) .GT. eps) THEN ! there is water, i.e. evaporate first
           ! Calculate mass change
           dM = lh * dt/lh_v
            M = theta_w(ksn) * rho_w * dzm_sn(ksn)
             ! Check that you only take the available amount of water
             IF(-dM .GE. M) THEN
               dM = -M
               theta_w(ksn) = theta_w(ksn) + dM/(rho_w*dzm_sn(ksn))
             ELSE
               theta_w(ksn) = theta_w(ksn) + dM/(rho_w*dzm_sn(ksn))
             ENDIF
           lh = lh - dM*lh_v/dt ! update energy used
         ELSEIF (theta_i(ksn) .GT. eps) THEN ! there is no water then sublimate ice matrix
           dM = lh * dt/lh_s
            M = theta_i(ksn) * rho_i * dzm_sn(ksn)
            IF(-dM .GT. M) THEN ! all ice can be sublimated
              dM = -M
              theta_i(ksn) = 0.0_vpp
              dzm_sn(ksn)  = 0.0_vpp
            ELSE
              dzm_sn_old = dzm_sn(ksn)
              dL = dM/rho_sn(ksn)
              dzm_sn(ksn) = dzm_sn(ksn) + dL
              theta_i(ksn) = theta_i(ksn) * (dzm_sn_old/dzm_sn(ksn))
              theta_i(ksn) = theta_i(ksn) + (dM/rho_i*dzm_sn(ksn))
              theta_w(ksn) = theta_w(ksn) * (dzm_sn_old/dzm_sn(ksn))
            ENDIF
           lh = lh - dM*lh_v/dt ! update energy used
         ENDIF
       ENDDO ! end of ksn

       ! MASSIVE HACK HERE
       IF(lh .LT. (-1.0_vpp*eps2)) THEN ! there is still energy left, which should technically be used by the soil layer for now let's erase it
         lh = 0.0_vpp
       ENDIF
     ENDDO ! end of while
   ENDIF
  ENDIF


! ==============================================================================
! - End subroutine
! ==============================================================================

END SUBROUTINE phase_change_surface





!==============================================================================
! + Begin subroutine transport_water
! =============================================================================

SUBROUTINE transport_water(zm_sn, dzm_sn, t_sn, rho_sn, theta_i, theta_a, &
                           theta_w, hcap_sn, top, runoff)

  !$acc routine seq

! ----------------------------------------------------------------------------
!
! Description: This subroutine transports the water down through the snow cover
!
! Method: Percolate the water down. Excess water that cannot be retained in the
!         next (lower) layer is stored in excess_water and moved down the domain.
!         Water that cannot be retained by the  lowest (snow/soil interface)
!         layer will be stored to not loose any mass.
!
! -----------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Begin declarations
! ------------------------------------------------------------------------------

  REAL    (KIND = vpp), INTENT(INOUT) :: &

    rho_sn(:)         , & ! snow layer density             (kg/m3)

    theta_i(:)        , & ! volumetric ice content         (m3/m3)
    theta_w(:)        , & ! volumetric water content       (m3/m3)
    theta_a(:)            ! volumentric air/void content   (m3/m3)

  REAL    (KIND = vpp), INTENT(INOUT) :: &

    runoff                ! melt water runoff

  ! List of IN parameter:
  REAL    (KIND = vpp), INTENT(IN) :: &

    zm_sn(:)          , & ! layer thickness                (m)
    dzm_sn(:)         , & ! layer thickness                (m)
    t_sn(:)           , & ! snow layer temperature         (K)
    hcap_sn(:)            ! heat capacity

  INTEGER           , INTENT(INOUT) :: &

    top             ! index of top snow layer

  ! List of local array's vectors and scalars
  REAL    (KIND = vpp) :: &

    frac_rho        , &    ! fraction of density water to ice
    limit_theta_i   , &    ! abc-formula
    w_up            , &    ! water content of upper layer
    w_low           , &    ! water content of lower layer
    dz_up           , &    ! layer thickness of upper layer
    dz_low          , &    ! dito for lower layer
    dtheta_w_up     , &    ! available water in upper layer
    dtheta_w_low    , &    ! available water in lower layer
    dtheta_w_low_x  , &    ! backup variable for dtheta_w_low
    theta_w_bot     , &    ! volumetric water content of bottom snow layer
    excess_water           ! excess water

  REAL    (KIND = vpp), DIMENSION(1:top) :: &

    w_res           , &    ! effective residual water content
    res_wat_cont    , &    ! potential residual water content
    dtheta_w             !  additional storage capacity (water) due to refreezing

   INTEGER        ::  &

     ksn                  ! dito for snow layer

! ------------------------------------------------------------------------------
! - End declarations
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Start subroutine arguments
! ------------------------------------------------------------------------------

   ! + Initiations
   ! ----------------------------------------

   ! Estimate resiudal water (RWC) content by volume; Coleou and Lesaffre, 1998,
   ! Ann. Glaciol., 26, 64-68

   ! Experimental range:
   !  - density unsoaked: 235 to 580
   !  - density soaked: 328 to 589 kg m-3
   !  - RWC by Mass 0.049 to 0.029
   ! NOTE: That function will limit range to 0.0264 to 0.08 RWC by Vol

   frac_rho      = rho_w/rho_i

   limit_theta_i =   1.0_vpp - frac_rho * ((1.0_vpp + 0.0165_vpp * frac_rho)                   &
                   - SQRT((1.0_vpp + 0.0165_vpp * frac_rho)*(1.0_vpp + 0.0165_vpp * frac_rho)  &
                   - 4.0_vpp * frac_rho * 0.0264_vpp)) / (2.0_vpp * frac_rho)


   ! + Calculate a few properties
   !--------------------

   ! Loop over number of snow layers
  
   DO ksn = 1, top, 1
   
     ! Determine the additional storage capacity (of water) due to refreezing
     ! --------------------------
      dtheta_w(ksn) = hcap_sn(ksn) * rho_sn(ksn) / lh_f / rho_w * MAX(0.0_vpp, (t0_melt - t_sn(ksn)))

     ! --------------------------
     ! Estimate resiudal water (RWC) content by volume; Coleou and Lesaffre, 1998,
     ! Ann. Glaciol., 26, 64-68
     ! --------------------------

     IF(theta_i(ksn) .GT. limit_theta_i) THEN

       ! This case is the limiting case where:
       ! theta_i + (theta_r * (rho_w/rho_i)) >= 1.0
       ! In that case, set the residual water content equal to the pore space

       res_wat_cont(ksn) = (1.0_vpp - theta_i(ksn)) * (rho_i/rho_w)

     ELSE

       IF(theta_i(ksn) .GT. 0.23_vpp) THEN

         res_wat_cont(ksn) = 0.0264_vpp + 0.0099_vpp * (1.0_vpp - theta_i(ksn)) / theta_i(ksn)

       ELSE

         res_wat_cont(ksn) = 0.08_vpp - 0.1023_vpp * (theta_i(ksn) - 0.03_vpp)

       ENDIF

     ENDIF

     ! Limit residual water content
     res_wat_cont(ksn) = MIN(res_wat_cont(ksn), 0.08_vpp)  ! NOTE: only needed in case of theta_i < 0.03


     ! Effective residual water content
     w_res(ksn) = MIN(1.0_vpp - theta_i(ksn) * rho_i/rho_w, res_wat_cont(ksn) + dtheta_w(ksn))
     w_res(ksn) = MAX(0.0_vpp, w_res(ksn))

  END DO

   ! Now start moving the water and adjust properties accordingly
   ! ----------------------------------

   DO ksn = top, 2, -1

     ! Reset excess water to zero
     IF(theta_i(ksn) .LT. eps) THEN ! no more ice in this layer only residual water which needs to be moved

       excess_water = theta_w(ksn) ! add residual water to excess water and ...
       theta_w(ksn) = 0.0_vpp       ! reset water content of said layer

!       top = top - 1

     ELSE

       excess_water = 0.0_vpp

     ENDIF

   ! water content of the upper layer
    w_up = theta_w(ksn)

     IF(ksn .EQ. top .AND. w_up .GT. 0.0_vpp .AND. w_up .LE. w_res(ksn)) THEN

       ! In that case you need to update the volumetric air content and the
       ! density of the top element
       ! as it may have caught some rain! Only top element should be considered,
       ! as when rain would have
       ! infiltrated lower elements as well, w_up > w_res.

       theta_a(ksn) = MAX(0.0_vpp, 1.0_vpp - theta_w(ksn) - theta_i(ksn))
       rho_sn(ksn)  = (theta_i(ksn) * rho_i) + (theta_w(ksn) * rho_w)

         ! we want positive densities
         IF(rho_sn(ksn) .LT. 0.0_vpp) THEN
           rho_sn(ksn) = 0.0_vpp
         ENDIF

     ENDIF


      IF(w_up .GT. w_res(ksn) .OR. excess_water .GT. 0.0_vpp) THEN

         ! ...water is being transfered
         dz_up       = dzm_sn(ksn)
         dz_low      = dzm_sn(ksn-1)

         w_low       = theta_w(ksn-1)

         dtheta_w_up = MAX(0.0_vpp, w_up - w_res(ksn))


         IF(dtheta_w_up .GT. 0.0_vpp .OR. excess_water .GT. 0.0_vpp) THEN

         ! ... dtheta_w_low is determined by also taking excess_water into
         ! account. Maybe excess_water can be stored in this layer.
         dtheta_w_low = dtheta_w_up * (dzm_sn(ksn)/dzm_sn(ksn-1)) + (excess_water/dzm_sn(ksn-1))

           ! now check whether there is enough air left - you might not be able
           ! to move the water
           ! or/and water may refreeze and expand specifically, you might want
           ! to create a water table over ice
           IF( (dtheta_w_low + w_low) .GT. (rho_i/rho_w * (1.0_vpp - theta_i(ksn-1))) ) THEN

             ! Deal with excess water ... Look how much you can leave in the
             ! lower layer (ksn-1).
                 ! If you have too much water even for the lower element (more
                 ! melt or rain per time
                 ! step than can be kept in this element), water is transferred
                 ! to excess_water.
                ! excess_water moves the water downward, trying to insert the
                ! water in lower layer.

             dtheta_w_low_x = dtheta_w_low    ! make backup
             dtheta_w_low   = MAX(0.0_vpp, (rho_i/rho_w * (1.0_vpp - theta_i(ksn-1)) - w_low))

             ! All the water that could not be stored in lower layer is
             ! considered excess_water.
             excess_water = (dtheta_w_low_x - dtheta_w_low) * dzm_sn(ksn-1)

           ELSE

             excess_water = 0.0_vpp

           ENDIF

           ! update volumetric contents, masses and density
           theta_w(ksn)   = w_up  - dtheta_w_up
           theta_w(ksn-1) = w_low + dtheta_w_low

         ENDIF ! end positive water movement

       ENDIF ! end if( W_upper > Wres )

  ENDDO ! loop snow layers


   ! ==========================
   ! Special treatment for lowermost snow layer, i.e runoff at bottom layer
   ! Note: If we want ponding on surface layer, i.e. glacier ice, sea ice, rock
   ! etc. we need to do it here, e.g. if itype = ...
   ! ==========================

   theta_w_bot = theta_w(1)

     IF (theta_w_bot .GT. w_res(1)) THEN

      ! Adjust dependent values accordingly
      theta_w(1) = w_res(1)
      theta_a(1) = 1.0_vpp - theta_w(1) - theta_i(1)

      ! Put all excess water of bottom layer  in runoff, i.e. move it out of the snow cover
      ! Note: if one comments this out you get ponding
      runoff              = runoff + excess_water +  dzm_sn(1) * (theta_w_bot - w_res(1))

     ENDIF

! ==============================================================================
! - End subroutine
! ==============================================================================

END SUBROUTINE transport_water




!==============================================================================
! + Begin subroutine aggregate_layers
! =============================================================================

SUBROUTINE aggregate_layers(top, hm_sn, zm_sn, dzm_sn, m_sn, theta_i, theta_w, theta_a, t_sn, &
                           rho_sn, hcap_sn, hcon_sn, hdif_sn )

  !$acc routine seq

! ----------------------------------------------------------------------------
!
! Description: This subroutine aggregates layers if one is smaller then a
!              pre-defined threshold. 
!
! Method:
!
! -----------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Begin declarations
! ------------------------------------------------------------------------------

 ! List of INOUT parameter:

   REAL    (KIND = vpp), INTENT(INOUT) :: &

    zm_sn(:)               , &  ! snow layer depth
    hm_sn(:)               , &  ! snow layer height
    dzm_sn(:)              , &  ! layer thickness

    hcap_sn(:)             , &  ! heat capacity
    hcon_sn(:)             , &  ! heat conductivity
    hdif_sn(:)             , &  ! heat diffusivity

    m_sn(:)                , &  ! mass
    rho_sn(:)              , &  ! density

    theta_i(:)             , &  ! volumetric ice content
    theta_w(:)             , &  ! volumetric water content
    theta_a(:)             , &  ! volumetric air content

    t_sn(:)                     ! temperature


   INTEGER , INTENT(INOUT)            :: &

    top                         ! top layer index

  ! List of local arrays, vectors and scalars
   INTEGER                         :: &

    ksn                    , &  ! lopp index for snow layer
    i, j                        ! utillity loop indicies

! ------------------------------------------------------------------------------
! - End declarations
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
! + Start subroutine arguments
! ------------------------------------------------------------------------------


    ! First loop: eliminiate 'empty' layers
    ! -------------------

     DO ksn = top, 1, -1

        IF(theta_i(ksn) .LT. 0.01_vpp) THEN

           ! ... and reset roperties
           hm_sn(ksn)   = 0.0_vpp
           zm_sn(ksn)   = 0.0_vpp
           dzm_sn(ksn)  = 0.0_vpp

           theta_i(ksn) = 0.0_vpp
           theta_w(ksn) = 0.0_vpp
           theta_a(ksn) = 0.0_vpp

           rho_sn(ksn)  = 0.0_vpp
           t_sn(ksn)    = 0.0_vpp
           m_sn(ksn)    = 0.0_vpp

           hcap_sn(ksn) = 0.0_vpp
           hcon_sn(ksn) = 0.0_vpp
           hdif_sn(ksn) = 0.0_vpp

           ! Update number of snow layers ...
           top         = top - 1

        ENDIF

      ENDDO


    ! Second loop: aggregate layers if possible
    ! -------------------

      DO ksn = top, 2, -1

        IF(dzm_sn(ksn) .LT. min_height_layer .OR. theta_i(ksn) .LT. 0.01_vpp) THEN ! layer is quite thin

          ! -------------------
          ! Aggregate with lower layer - Adjust properties
          ! -------------------

           dzm_sn(ksn-1)  = dzm_sn(ksn) + dzm_sn(ksn-1)        ! layer thickness

           theta_i(ksn-1) = (theta_i(ksn)*dzm_sn(ksn)                                      &
                          + theta_i(ksn-1)*dzm_sn(ksn-1)) / (dzm_sn(ksn) + dzm_sn(ksn-1))  ! volumetric ice content

           theta_w(ksn-1) = (theta_w(ksn)*dzm_sn(ksn)                                      &
                          + theta_w(ksn-1)*dzm_sn(ksn-1)) / (dzm_sn(ksn) + dzm_sn(ksn-1))  ! volumetric water content


          ! ... and reset roperties
           hm_sn(ksn)   = 0.0_vpp
           zm_sn(ksn)   = 0.0_vpp
           dzm_sn(ksn)  = 0.0_vpp

           theta_i(ksn) = 0.0_vpp
           theta_w(ksn) = 0.0_vpp
           theta_a(ksn) = 0.0_vpp

           rho_sn(ksn)  = 0.0_vpp
           t_sn(ksn)    = 0.0_vpp
           m_sn(ksn)    = 0.0_vpp

           hcap_sn(ksn) = 0.0_vpp
           hcon_sn(ksn) = 0.0_vpp
           hdif_sn(ksn) = 0.0_vpp

           ! Update number of snow layers ...
            top         = top - 1

        ENDIF

      ENDDO


    ! Third loop: it might be that subsurface layers melted, if so shift layer boundaries accordingly
    ! -------------------

     DO i = 1, ke_snow, 1

       IF(dzm_sn(i) .LT. min_height_layer .OR. theta_i(i) .LT. 0.01_vpp) THEN

         DO j = i+1, ke_snow, 1

           dzm_sn(j-1)  = dzm_sn(j)

           theta_i(j-1) = theta_i(j)
           theta_w(j-1) = theta_w(j)

           rho_sn(j-1)   = rho_sn(j)

         ENDDO

       ENDIF

     ENDDO


       ! Sepecial treatment of top layer
       IF(dzm_sn(top) .LT. min_height_layer .OR. theta_i(top) .LT. 0.01_vpp) THEN

         top = top - 1

       ENDIF


    ! Fourth loop:  re-calculate layer depth from thickness
    ! -------------------

!     DO ksn = ke_snow, 1, -1
!
!       IF(ksn .EQ. ke_snow) THEN
!
!         zm_sn(ksn) = dzm_sn(ksn)
!
!       ELSE
!
!         zm_sn(ksn) = zm_sn(ksn+1) + dzm_sn(ksn)
!
!       ENDIF
!
!     ENDDO



! ==============================================================================
! - End subroutine
! ==============================================================================

END SUBROUTINE aggregate_layers




!==============================================================================
! + Begin subroutine clear_layers
! =============================================================================

SUBROUTINE clear_layers(top, hm_sn, zm_sn, dzm_sn, m_sn, theta_i, theta_i_old    , &
                        theta_w, theta_a, t_sn, rho_sn, hcap_sn, hcon_sn, hdif_sn  )

  !$acc routine seq

! ----------------------------------------------------------------------------
!
! Description: This clear the vectors from artificial layer info, i.e. sets it
!              to 0.0_vpp
!
! Method:
!
! -----------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Begin declarations
! ------------------------------------------------------------------------------

  ! List of INOUT parameter:

   REAL    (KIND = vpp), INTENT(INOUT) :: &

    zm_sn(:)               , &  ! snow layer depth
    hm_sn(:)               , &  ! snow layer height
    dzm_sn(:)              , &  ! layer thickness

    hcap_sn(:)             , &  ! heat capacity
    hcon_sn(:)             , &  ! heat conductivity
    hdif_sn(:)             , &  ! heat diffusivity

    m_sn(:)                , &  ! mass
    rho_sn(:)              , &  ! density

    theta_i(:)             , &  ! volumetric ice content
    theta_i_old(:)         , &  ! volumetric ice content (old values before phase change etc.)
    theta_w(:)             , &  ! volumetric water content
    theta_a(:)             , &  ! volumetric air content

    t_sn(:)                     ! temperature


   INTEGER , INTENT(IN)            :: &

    top                         ! top layer index

  ! List of local arrays, vectors and scalars
   INTEGER                         :: &

    ksn

! ------------------------------------------------------------------------------
! - End declarations
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Start subroutine arguments
! ------------------------------------------------------------------------------

 ! + Clear vector
 ! -------------------------

   DO ksn = top+1, ke_snow, 1

     hm_sn(ksn)   = 0.0_vpp
     zm_sn(ksn)   = 0.0_vpp
     dzm_sn(ksn)  = 0.0_vpp

     theta_i(ksn) = 0.0_vpp
     theta_w(ksn) = 0.0_vpp
     theta_a(ksn) = 0.0_vpp

     theta_i_old(ksn) = 0.0_vpp

     rho_sn(ksn)  = 0.0_vpp
     t_sn(ksn)    = 0.0_vpp
     m_sn(ksn)    = 0.0_vpp

     hcap_sn(ksn) = 0.0_vpp
     hcon_sn(ksn) = 0.0_vpp
     hdif_sn(ksn) = 0.0_vpp

   ENDDO

! ==============================================================================
! - End subroutine
! ==============================================================================

END SUBROUTINE clear_layers





!==============================================================================
! + Begin subroutine calc_set_rate
! =============================================================================

SUBROUTINE calc_set_rate(h_snow, zm_sn, dzm_sn, t_sn, m_sn, theta_i, theta_w, &
                         theta_a, theta_i_old, rho_sn, top, melt_flag, dt)

  !$acc routine seq

! ----------------------------------------------------------------------------
!
! Description: Calculates settling rates due to ice loss and overburden stress
!
! Method: Andersen (1976) or Vionnett et al. (2012)
!
! -----------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Begin declarations
! ------------------------------------------------------------------------------

  ! List of INOUT parameter:
   REAL    (KIND = vpp),  INTENT(INOUT) :: &

     zm_sn(:)         , &   ! current main layer depth
     dzm_sn(:)        , &   ! current main layer thickness

     t_sn(:)          , &   ! snow layer temperature
     m_sn(:)          , &   ! snow layer mass

     theta_i(:)       , &   ! volumetric ice content
     theta_w(:)       , &   ! volumetric water content
     theta_a(:)       , &   ! volumetric air content

     rho_sn(:)              ! snow layer density

   REAL    (KIND = vpp), INTENT(INOUT) :: &

     h_snow              ! snow height

  ! List of IN parameter:
   REAL    (KIND = vpp), INTENT(IN) :: &

     theta_i_old(:)    , &  ! old volumetric ice content
     dt                     ! integration time step

   INTEGER, INTENT(IN) :: &

     top


   LOGICAL, INTENT(IN) :: &

     melt_flag

  ! List of local arrays, vectors and scalars

   INTEGER                  :: &

     ksn                                     ! loop index for snow layers

    ! After Andresen (1976)
   REAL(KIND=vpp), PARAMETER ::  &

     c2       = 23.0E-3_vpp              , &  ! Coefficinet [m3/kg]
     eta0     = 9.0E5_vpp                , &  ! The Viscosity Coefficient Eta0  [kg-s/m2]
     c_factor = 0.08_vpp                      ! snow compaction overburden exponential factor (1/K)

   REAL    (KIND = vpp) :: &

     rate_1                              , & ! settling rate for ice loss (s**-1)
     rate_2                              , & ! overburden stress          (s**-1)

     overburden                          , & ! overburden load

     tot_rate                            , & ! total settling rates       (s**-1)

     ddz                                     ! change of layer thickness

   REAL    (KIND = vpp)                          :: &

     dz_old                              , & ! old layer thickness before settling     (m)
     dz_new                                  ! new layer thickness after settling      (m)

   ! After Vionnet (2012)
   REAL(KIND=vpp), PARAMETER ::  &

     a_eta = 0.1_vpp                     , &   !  default  0.1_vpp
     b_eta = 0.023_vpp                   , &   !           0.023_vpp
     c_eta = 250.0_vpp                   , &   !           250.0_vpp
     eta_0 = 7.62237E6_vpp                     !           7.62237E6_vpp

     REAL    (KIND = vpp) :: &

     f1                              , &
     f2                              , &
     eta


! ------------------------------------------------------------------------------
! - End declarations
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Start subroutine arguments
! ------------------------------------------------------------------------------

 ! ------------------------------------------
 ! + Compaction/Settling
 ! -------------------------------------------

     ! Initalizations
     ddz          = 0.0_vpp   ! total change in layer thickness
     overburden   = 0.0_vpp   ! overburden stress (weight of the overlying layers times g)

     DO ksn = top, 1, -1

       ! Reset values for settling rates
       rate_1    = 0.0_vpp      ! due to ice loss
       rate_2    = 0.0_vpp      ! due to overburden stress

       tot_rate  = 0.0_vpp      ! total settliing rate

       dz_old    = 0.0_vpp      ! old layer thickness (before settling) ...
       dz_new    = 0.0_vpp      ! ... and after settling

       !  ... due to loss of ice
       !-------------------

 !      IF(melt_flag) THEN

         IF(theta_i(ksn) <  theta_i_old(ksn)) THEN ! melting ocured in that layers

           rate_1 = - 1.0_vpp/dt * MAX(0.0_vpp,  (theta_i_old(ksn) - theta_i(ksn)) / theta_i_old(ksn) )

         ENDIF

 !      ENDIF

       !  ... due to overburden stress
       !-------------------

      ! Vionnett (2012)
       f1   = 1.0_vpp / (1.0_vpp + 60.0_vpp * ((theta_w(ksn)*rho_w*dzm_sn(ksn)) / (rho_w * dzm_sn(ksn))))
       f2   = 1.0_vpp

       eta = f1 * f2 * eta_0 * (rho_sn(ksn)/c_eta) * exp(a_eta*(t0_melt - t_sn(ksn)) + b_eta*rho_sn(ksn))

       rate_2 = - (overburden + m_sn(ksn)*g/2.0_vpp) / eta


       ! Andersen (1976)
!        rate_2 = -(overburden + m_sn(ksn)*g/2.0_vpp)*EXP(-c_factor * (t_sn(ksn) - t0_melt) - c2*rho_sn(ksn)) / eta0

       ! increase overburden stress NOTE: How to deal with slope angle should be overburden = m*g*cos(alpha)
       overburden = overburden + (m_sn(ksn) * g)

      ! ... calculate change ...
       !-------------------

       ! ... of all (ice loss, overburden, destructive) settling rates (1/s)
       tot_rate = (rate_1*dt) + (rate_2*dt)

       ! ... of layer thickness, i.e. sum over all layers (m)
       ddz          = ddz + MAX(-1.0_vpp * dzm_sn(ksn), dzm_sn(ksn) * tot_rate)

       dz_old = dzm_sn(ksn)
       dz_new = dz_old + MAX(-1.0_vpp * dzm_sn(ksn), dzm_sn(ksn) * tot_rate)


       ! ... volumetric contents
       theta_i(ksn) = MAX(0.0_vpp, theta_i(ksn) * (dz_old / dz_new))    ! ice content
       theta_w(ksn) = MAX(0.0_vpp, theta_w(ksn) * (dz_old / dz_new))    ! water content

       ! ... of layer thickness (m)
       dzm_sn(ksn) = dz_new

  ENDDO

  ! -------------------
  ! + Re-calculate layer depth from layer thickness
  ! -------------------

    DO ksn = top, 1, -1

      IF(ksn .EQ. top) THEN

        zm_sn(ksn) = dzm_sn(ksn)

      ELSE

        zm_sn(ksn) = zm_sn(ksn+1) + dzm_sn(ksn)

      ENDIF

    ENDDO

    ! Calculate change in snow depth
    h_snow = zm_sn(1)


! ==============================================================================
! - End subroutine
! ==============================================================================

END SUBROUTINE calc_set_rate



! ==============================================================================
! + DONE, DONE!!!!
! =============================================================================

END MODULE sfc_snow_utilities



















