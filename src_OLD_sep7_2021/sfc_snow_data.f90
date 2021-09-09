!-----------------------------------------------------------------------------
!+ Data module 'sfc_snow_data'
!------------------------------------------------------------------------------

MODULE sfc_snow_data

!------------------------------------------------------------------------------
!
! Description:
!  This module declares and initializes all parametric scalar and array
!  data which are used in the snow model
!
!
!  Note : Most of this is probably going to sfc_terra_data.f90 or defined as
!  namelist option.
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
!
! =============================================================================
! + Begin loading modules
! =============================================================================

#ifdef __COSMO__
USE sfc_terra_data, ONLY  : &
         eps => eps_div 

USE kind_parameters, ONLY :   &
#elif __ICON__
USE mo_kind, ONLY:     &
#endif
    vpp           ! KIND-type parameter for real variables (variable precisionphysics)

! =============================================================================
! - End loading modules
! =============================================================================

IMPLICIT NONE

PUBLIC           ! All constants and variables in this module are public

! =============================================================================
! + Start declarations
! =============================================================================

! Global (i.e. public) declarations

! 1. Model configurations
! -------------------------------------------------

!  INTEGER                    :: &

!    ke_snow    = 10     , &        ! maximum number of layers, i.e. vertical extent of fields   (-)
!    n_substrate = 1                  ! number of substrate layer, i.e. layers below the snow
                                     ! NOTE: This is ke_soil!!!!!
  
!    !$acc declare create (ke_snow, n_substrate)


  REAL  (KIND=vpp)           ::  &

    min_height_layer = 0.001_vpp       , &  ! minimum layer thickness
    max_height_layer = 0.01_vpp            ! maximum height/thickness of new snow layers                 (m)

    !$acc declare create (min_height_layer, max_height_layer)


  REAL (KIND=vpp), PARAMETER        ::    &

    height_meteo_values = 5.5_vpp       , & ! height of meteo values above ground soil                     (m)
                                             ! NOTE: This should be the replaced
                                             ! by the first atmospheric level.
                                             ! see calc_turb()
    roughness_length    = 0.002_vpp          ! roughness length                                             (m)



! 2. Utility constants
! -------------------------------------------------

 REAL  (KIND=vpp), PARAMETER     ::      &

   eps2   = eps*eps                                          ! and an even smaller number


! 3. Physical constants  NOTE: a bunch of these already exists globally, but
!                              might have slithly different values so double check and replace at some point
! -------------------------------------------------


  REAL  (KIND=vpp), PARAMETER           ::  &

    t0_melt                     = 273.15_vpp                               , & ! melting temperature            (K)
   
    rho_i                       = 917.0_vpp                                , & ! density of ice                 (kg/m**3)
    rho_a                       = 1.1_vpp                                  , & ! density of air                 (kg/m**3)
    rho_w                       = 1000.0_vpp                               , & ! density of water               (kg/m**3)


    specific_heat_air           = 1004.67_vpp                              , & ! it is what it is
    specific_heat_ice           = 2100.0_vpp                               , & ! same here
    specific_heat_water         = 4190.0_vpp                               , & ! and here

    gas_constant_air            = 287.05_vpp                               , & ! gas constant for dry air                    (J kg-1 K-1)
    gas_constant_water_vapor    = 461.51_vpp                               , & ! gas constant for water vapor

    rvd_m_o  =  gas_constant_water_vapor/gas_constant_air - 1.0_vpp        , & ! should be obvious

    lh_v                        = 2.501E6_vpp                              , & ! latent heat of vapourization                (J kg-1)
    lh_f                        = 0.334E6_vpp                              , & ! latent heat of fusion                       (J kg-1)
    lh_s                        = 2.835E6_vpp                              , & ! latent heat of sublimation                  (J kg-1)

    sigma                       = 5.6697E-8_vpp                            , & !  Stefan-Boltzmann constant
    snow_Ctalb                  = 0.004_vpp                                , & !  thermal albedo

    cn                          = 0.5_vpp                                  , & ! Cranck-Nicholson Factor

    karman                      = 0.4_vpp                                  , & ! Karman constant
    g                           = 9.80665_vpp                              , & ! accelaration due to gravity

    theta_s                     = 1.0_vpp                                  , & ! Saturated Water Content, for now we say 1.0
    theta_r                     = 0.0_vpp                                      ! Minimum amount of liquid water that will remain



! ==============================================================================


! Allocate work arrays for the multi-layer snow cover scheme (SNOWPOLINO)

#ifdef ALLOC_WKARR

! Local arrays defined as allocatables here:
! ------------------------------------------

   REAL (KIND=vpp), ALLOCATABLE  :: &

     t_sn_sfc   (:)    , &  ! snow surface temperature                    [K]
     sh_sn      (:)    , &  ! sensible heat flux                          [W/m**2]
     lh_sn      (:)    , &  ! latent heat flux                            [W/m**2]
     tch_sn     (:)    , &  ! transfer coefficient
     alpha_sn   (:)    , &  ! snow surface albdeo                         [-]
     swnet_sn   (:)    , &  ! net short wave radiation                    [W/m**2]
     for_sn     (:)    , &  ! total atmospheric forcing at snow surface   [W/m**2]
     runoff_sn  (:)    , &  ! total melt runoff
     swe_sn     (:)         ! local S.W.E calculation
  
   INTEGER, ALLOCATABLE          :: &

     top        (:)         ! index of first (top) snow layer  [-]

   REAL(KIND=vpp), ALLOCATABLE ::  &

     hm_sn       (:,:)  , & ! height (from bottom) of snow layers (main level)     (m)
     zm_sn       (:,:)  , & ! depth (from top) of snow layers (main level)         (m)

     rho_sn      (:,:)  , & ! density of snow layers (kg/m**3)
     m_sn        (:,:)  , & ! layer mass (kg)
     hcap_sn     (:,:)  , & ! snow layer heat capacity
     hcon_sn     (:,:)  , & !            heat conductivity
     hdif_sn     (:,:)  , & !            heat diffusion

     swabs_sn    (:,:)  , & ! absorbed short wave radiation

     zm_sn_old   (:,:)  , & ! old value of layer depth
     theta_i_old (:,:)  , & !              volumetric ice content
     theta_w_old (:,:)      !              volumetric water content

    LOGICAL, ALLOCATABLE  :: &

     melt_flag   (:)    , & ! flag indicating melting; IF THEN TRUE
     freeze_flag (:)        !      indicating freezing


    REAL(KIND=vpp), ALLOCATABLE     :: &

     tmp_sn      (:)      ! temporary array

!US all these variables are declared, allocated and deallocated in sfc_terra_data
!   REAL(KIND=vpp), ALLOCATABLE     :: &

!    ziw_fr      (:,:)  , &  ! fractional ice content of soil layer
!    zlw_fr      (:,:)  , &  ! fractional liqu. water content of soil layer
!    zw_fr       (:,:)  , &  ! fractional total water content of soil layers
!    zroc        (:,:)  , &  ! heat capacity of soil layers
!    zrocg       (:,:)  , &  ! total volumetric heat capacity of soil
!    zrocg_soil  (:,:)  , &  !

!    zalam       (:,:)  , &
!    zalamtmp    (:,:)  , &
!    hzalam      (:,:)  , &

!    zfcap       (:,:)  , &
!    zporv       (:,:)  , &
!    zpwp        (:,:)  , &
!    zdlam       (:) 

! ==============================================================================

CONTAINS

! ==============================================================================

SUBROUTINE snow_wkarr_alloc (ke_soil, ke_snow, nproma, istat)

  INTEGER, INTENT(IN)  :: ke_soil, ke_snow, nproma
  INTEGER, INTENT(OUT) :: istat

  istat = 0


  ALLOCATE ( t_sn_sfc    (nproma)    , &  ! snow surface temperature            [K]
             sh_sn       (nproma)    , &  ! sensible heat flux                          [W/m**2]
             lh_sn       (nproma)    , &  ! latent heat flux                            [W/m**2]
             tch_sn      (nproma)    , &  ! transfer coefficient
             alpha_sn    (nproma)    , &  ! snow surface albdeo                         [-]
             swnet_sn    (nproma)    , &  ! net short wave radiation                    [W/m**2]
             for_sn      (nproma)    , &  ! total atmospheric forcing at snow surface   [W/m**2]
             runoff_sn   (nproma)    , &  ! total melt runoff
             swe_sn      (nproma)    , &  ! local S.W.E calculation
       STAT=istat)
  !$acc enter data async                                               &
  !$acc create(t_sn_sfc, sh_sn, lh_sn, tch_sn,swe_sn)                  &
  !$acc create(alpha_sn, swnet_sn, for_sn, runoff_sn)


  ALLOCATE ( top         (nproma)    , &  ! index of first (top) snow layer  [-]
       STAT=istat)
  !$acc enter data async create(top)


  ALLOCATE ( hm_sn       (nproma,ke_snow)  , & ! height (from bottom) of snow layers (main level)     (m)
             zm_sn       (nproma,ke_snow)  , & ! depth (from top) of snow layers (main level)         (m)
             rho_sn      (nproma,ke_snow)  , & ! density of snow layers (kg/m**3)
             m_sn        (nproma,ke_snow)  , & ! layer mass (kg)
             hcap_sn     (nproma,ke_snow)  , & ! snow layer heat capacity
             hcon_sn     (nproma,ke_snow)  , & !            heat conductivity
             hdif_sn     (nproma,ke_snow)  , & !            heat diffusion
             swabs_sn    (nproma,ke_snow)  , & ! absorbed short wave radiation
             zm_sn_old   (nproma,ke_snow)  , & ! old value of layer depth
             theta_i_old (nproma,ke_snow)  , & !              volumetric ice content
             theta_w_old (nproma,ke_snow)  , & !              volumetric water content
       STAT=istat)
  !$acc enter data async create(hm_sn, zm_sn, rho_sn, m_sn, hcap_sn)   &
  !$acc create(hcon_sn, hdif_sn, swabs_sn, zm_sn_old, theta_i_old)     &
  !$acc create(theta_w_old)


  ALLOCATE ( melt_flag   (nproma)    , &  ! flag indicating melting; IF THEN TRUE
             freeze_flag (nproma)    , &  !      indicating freezing
       STAT=istat)
  !$acc enter data async create(melt_flag, freeze_flag)


  ALLOCATE ( tmp_sn     (ke_snow)  , & ! temporary array
       STAT=istat)
  !$acc enter data async create(tmp_sn)

 
!ALLOCATE (  ziw_fr      (nproma,ke_soil+1)  , & ! fractional ice content of soil layer
!            zlw_fr      (nproma,ke_soil+1)  , & ! fractional liquid water content of soil layer
!            zw_fr       (nproma,ke_soil+1)  , & ! fractional total water content of soil layers
!            zroc        (nproma,ke_soil+1)  , & !
!            zrocg       (nproma,ke_soil+1)  , & !
!            zrocg_soil  (nproma,ke_soil+1)  , & !
!           
!            zalam       (nproma,ke_soil  )  , & ! heat capacity
!            zalamtmp    (nproma,ke_soil  )  , & ! heat capacity
!            hzalam      (nproma,ke_soil+1)  , & ! heat capacity

!            zfcap       (nproma,ke_soil+1)  , &
!            zporv       (nproma,ke_soil+1)  , &
!            zpwp        (nproma,ke_soil+1)  , &
!            zdlam       (nproma)            , &       
!      STAT=istat)
! !$acc enter data async                                              &
! !$acc create(ziw_fr, zlw_fr, zw_fr, zroc, zrocg)                    &
! !$acc create(zrocg_soil, zalam, zalamtmp, hzalam, zfcap, zporv)     &
! !$acc create(zpwp, zdlam)


END SUBROUTINE snow_wkarr_alloc

!==============================================================================
!==============================================================================

SUBROUTINE snow_wkarr_dealloc (istat)

  INTEGER, INTENT(OUT) :: istat

  istat = 0


  !$acc exit data async                                               &
  !$acc delete(t_sn_sfc, sh_sn, lh_sn)                                &
  !$acc delete(tch_sn, alpha_sn, swnet_sn, for_sn, runoff_sn,swe_sn)  &
  !$acc delete(top, hm_sn, zm_sn, rho_sn, m_sn, hcap_sn, hcon_sn)     &
  !$acc delete(hdif_sn, swabs_sn, zm_sn_old, theta_i_old)             &
  !$acc delete(theta_w_old, melt_flag, freeze_flag, tmp_sn)
! !$acc delete(ziw_fr, zlw_fr, zw_fr, zroc, zrocg, zrocg_soil)        &
! !$acc delete(zalam, zalamtmp, hzalam, zfcap, zporv, zpwp, zdlam)

  DEALLOCATE (t_sn_sfc, sh_sn, lh_sn, tch_sn, alpha_sn, swnet_sn     , &
              for_sn, runoff_sn,swe_sn                                      , &
       STAT=istat)

  DEALLOCATE (top                                                    , &
       STAT=istat)

  DEALLOCATE (hm_sn, zm_sn, rho_sn, m_sn, hcap_sn, hcon_sn, hdif_sn  , &
              swabs_sn, zm_sn_old, theta_i_old, theta_w_old          , &
       STAT=istat)

  DEALLOCATE (melt_flag, freeze_flag                                 , &
       STAT=istat)

  DEALLOCATE (tmp_sn                                                 , &
       STAT=istat)

 
! DEALLOCATE (ziw_fr, zlw_fr, zw_fr, zroc, zrocg, zrocg_soil         , &
!             zalam, zalamtmp, hzalam, zfcap, zporv, zpwp, zdlam     , &
!      STAT=istat)

END SUBROUTINE snow_wkarr_dealloc

!==============================================================================
#endif
      !ALLOC_WKARR

!==============================================================================
! - End declarations
! =============================================================================


END MODULE sfc_snow_data



































