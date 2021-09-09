!------------------------------------------------------------------------------
!+ Source module  "sfc_snow"
!------------------------------------------------------------------------------

MODULE sfc_snow

!------------------------------------------------------------------------------
!
! Description:
!   This module contains the multi-layer snow cover scheme (SNOWPOLINO), called within
!   sfc_terra.f90 and calls additionally required subroutines from
!   sfc_snow_utilities.f90. 
!
!
! Note:
!   For now we have this in one file which is huge/long. It might make
!   sense to split this in more modules.
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

! General modules
#ifdef __COSMO__
USE kind_parameters, ONLY :   &
#elif __ICON__
USE mo_kind, ONLY:     &
#endif
    vpp           ! KIND-type parameter for real variables (variable precisionphysics)

!
#ifdef __COSMO__
USE data_runcontrol, ONLY: ntstep, lsnow, iid_dbg, ibl_b, isc_b
USE data_parallel,   ONLY: my_cart_id
#endif

USE sfc_terra_data

! Additional snow scheme modules
USE sfc_snow_data           ! contains all relevant data and setting
USE sfc_snow_utilities      ! contains all required subroutines (see list in module header)


! =============================================================================
! - End loading modules
! =============================================================================

  IMPLICIT NONE

  PUBLIC           ! All constants and variables in this module are public

  CONTAINS

! =============================================================================
!  + Begin subroutine  - snowpolino 
! =============================================================================

  SUBROUTINE snowpolino (                  & 

                        ! general
                        nvec             , & ! array dimensions
                        ivstart          , & ! start index
                        ivend            , & ! end index
                        iblock           , & ! number of block
                        ke_soil          , & ! number of soil layers without lowermost (climat.) soil layer
                        dt               , & ! time step
                        nclass_gscp      , & ! number of hydrometeor classes of grid scale microphysics

                        ! meteorological
                        u                , & ! zonal wind speed
                        v                , & ! meridional wind speed
                        t                , & ! temperature (air)
                        qv               , & ! specific water vapor content                 (kg/kg)
                        ps               , & ! surface pressure                             ( Pa  )
                        prr_con          , & ! precipitation rate ofrain, convective        (kg/m2*s)
                        prs_con          , & ! precipitation rate of snow, convective       (kg/m2*s)
                        conv_frac        , & ! convective area fraction
                        prr_gsp          , & ! precipitation rate of rain, grid-scale       (kg/m2*s)
                        prs_gsp          , & ! precipitation rate of snow, grid-scale       (kg/m2*s)
                        prg_gsp          , & ! precipitation rate of graupel, grid-scale    (kg/m2*s)
#ifdef TWOMOM_SB
                        prh_gsp          , & ! precipitation rate of hail, grid-scale       (kg/m2*s)
#endif
                        sobs             , & ! net shortwave radiation ( W/m2 )
                        thbs             , & ! net longwave radiation  ( W/m2 )
                        !swdir_s          , & ! direct short wave radiation at ground         (W/m2)
                        !swdifd_s         , & ! diffuse short wave radiation at ground        (W/m2)
                        !swdifu_s         , & ! diffuse/upward short wave radiation at ground (W/m2)
                        !lwd_s            , & ! direct/icoming short wave radiation at ground (W/m2)
                        !lwu_s            , & ! emitted long wave radiation at ground         (W/m2)

                        ! surface
                        t_snow_now       , & ! temperature of the snow-surface              (  K  )
                        t_snow_new       , & ! temperature of the snow-surface              (  K  )
                        zshfl_snow       , & ! sensible heat flux snow/air interface        ( W/m2)
                        zlhfl_snow       , & ! latent   heat flux snow/air interface        ( W/m2)

                        ! snow
                        t_sn_now         , & ! snow temperature (main level)                (  K  )
                        t_sn_new         , & ! snow temperature (main level)                (  K  )
                        theta_i_now      , & ! volumetric ice content                       (  -  )
                        theta_i_new      , & ! volumetric ice content                       (  -  )
                        theta_w_now      , & ! volumetric water content                     (  -  )
                        theta_w_new      , & ! volumetric water content                     (  -  )
                        theta_a_now      , & ! volumetric air content                       (  -  )
                        theta_a_new      , & ! volumetric air content                       (  -  )
                        dzm_sn_now       , & ! snow layer thickness                         (  -  )
                        dzm_sn_new       , & ! snow layer thickness                         (  -  )
                        hn_sn_now        , & ! new snow amounts                             (  m  )
                        hn_sn_new        , & ! new snow amounts                             (  m  )           

                        top_sn_now       , & ! index of the first (top) snow layer          (  -  )
                        top_sn_new       , & ! index of the first (top) snow layer          (  -  )

                        h_snow           , & ! snow height                                  (m)

                        w_snow_now       , & ! water content of snow                        (m H2O)
                        w_snow_new       , & ! water content of snow                        (m H2O)
 
                        rho_snow_new     , & ! equivalent density of snow layer           ( kg / m3) 

                        ! soil
                        t_so_now         , & ! soil temperature (main level)                (  K  )
                        t_so_new         , & ! soil temperature (main level)                (  K  )
                        w_so_now         , & ! total water conent (ice + liquid water)      (m H20)
                        w_so_new         , & ! total water conent (ice + liquid water)      (m H20)
                        w_so_ice_now     , & ! ice content                                  (m H20)
                        w_so_ice_new     , & ! ice content                                  (m H20)
                        zmls             , & ! processing soil level structure
                        plcov            , & ! fraction of plant cover                        --
                        rootdp           , & ! depth of the roots                           ( m  )
                        soiltyp_subs       ) ! type of the soil (keys 0-9)                    --

        
! =============================================================================
!
! Description: This subroutine contains and calls all required code for the new
!     multi-layer snow cover scheme developed during the course of the priority task
!     projects SAINT (Snow Cover Atmosphere INTeraction). 
!
! Method:
!
! Dependencies:
!
! Refernces:
!
! Notes:
!
! =============================================================================

! -----------------------------------------------------------------------------
! + Begin declarations
! -----------------------------------------------------------------------------

! ------------------------
! + Global
! -----------------------

  ! generell
  INTEGER, INTENT(IN)  ::  &
                  nvec,              & ! array dimensions
                  ivstart,           & ! start index
                  ivend,             & ! end index
                  iblock,            & ! number of block
                  ke_soil              ! number of soil layers

  REAL    (KIND = vpp), INTENT(IN)  ::  &
                  dt                   ! time step

  INTEGER, INTENT(IN)  :: &
                  nclass_gscp          ! number of hydrometeor classes of grid scale microphysics

  ! meteorological
  REAL    (KIND = vpp), DIMENSION(nvec), INTENT(IN) :: &
                  u                , & ! zonal wind speed                              ( m/s )
                  v                , & ! meridional wind speed                         ( m/s )
                  t                , & ! temperature                                   (  K  )
                  qv               , & ! specific water vapor content                  (kg/kg)
                  ps               , & ! surface pressure                              ( Pa  )
                  prr_con          , & ! precipitation rate of rain, convective        (kg/m2*s)
                  prs_con          , & ! precipitation rate of snow, convective        (kg/m2*s)
                  conv_frac        , & ! convective area fraction
                  prr_gsp          , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
                  prs_gsp          , & ! precipitation rate of snow, grid-scale        (kg/m2*s)
#ifdef TWOMOM_SB
                  prh_gsp          , & ! precipitation rate of hail, grid-scale        (kg/m2*s)
#endif
                  prg_gsp              ! precipitation rate of graupel, grid-scale     (kg/m2*s)

   REAL    (KIND = vpp), DIMENSION(nvec), INTENT(IN) :: &
                  sobs             , & ! net shortwave radiation (W/m2)
                  thbs                 ! net longwave radiation  (W/m2)

!  REAL    (KIND = vpp), DIMENSION(nvec), INTENT(INOUT) :: &
!                  swdir_s          , & ! direct short wave radiation at ground         (W/m2)
!                  swdifd_s         , & ! diffuse short wave radiation at ground        (W/m2)
!                  swdifu_s         , & ! diffuse/upward short wave radiation at ground (W/m2)
!                  lwd_s            , & ! direct/icoming short wave radiation at ground (W/m2)
!                  lwu_s                ! emitted long wave radiation at ground         (W/m2)

  ! surface
  REAL    (KIND = vpp), DIMENSION(nvec), INTENT(INOUT) :: &
                  t_snow_now           ! temperature of the snow-surface (K)

  REAL    (KIND = vpp), DIMENSION(nvec), INTENT(OUT) :: &
                  t_snow_new           !

  REAL    (KIND = vpp), DIMENSION(nvec), INTENT(OUT) :: &
                  zshfl_snow       , & ! sensible heat flux snow/air interface         (W/m2)
                  zlhfl_snow           ! latent   heat flux snow/air interface         (W/m2)

  ! snow
   REAL    (KIND = vpp)    , DIMENSION(nvec), INTENT(INOUT)  :: &
                  h_snow                  ! snow height


   REAL    (KIND = vpp)    , DIMENSION(nvec,1:n_layers), INTENT(INOUT) :: &
                  t_sn_now          , &   ! snow temperature (main level)              (K)
                  theta_i_now       , &   ! volumetric ice content                     (-)
                  theta_w_now       , &   ! water ice content                          (-)
                  theta_a_now       , &   ! air ice content                            (-)
                  dzm_sn_now              ! layer thickness between main levels        (m)

   REAL    (KIND = vpp)    , DIMENSION(nvec,1:n_layers), INTENT(OUT) :: &
                  t_sn_new          , &   ! snow temperature (main level)              (K)
                  theta_i_new       , &   ! volumetric ice content                     (-)
                  theta_w_new       , &   ! water ice content                          (-)
                  theta_a_new       , &   ! air ice content                            (-)
                  dzm_sn_new              ! layer thickness between main levels        (m)

   REAL    (KIND = vpp)    , DIMENSION(nvec), INTENT(INOUT) :: &
                  hn_sn_now              ! new snow amounts (storage)                  (m)
   REAL    (KIND = vpp)    , DIMENSION(nvec), INTENT(OUT)   :: &
                  hn_sn_new

   REAL    (KIND = vpp)    , DIMENSION(nvec), INTENT(INOUT) :: &
                  top_sn_now              ! index of the first (top) layer index       (-)
   REAL    (KIND = vpp)    , DIMENSION(nvec), INTENT(OUT)   :: &  
                  top_sn_new

   REAL    (KIND = vpp)    , DIMENSION(nvec), INTENT(INOUT) :: &
                  w_snow_now              ! water content (snow water equivialent) of snow cover        (-)
   REAL    (KIND = vpp)    , DIMENSION(nvec), INTENT(OUT)   :: &
                  w_snow_new

   REAL    (KIND = vpp)    , DIMENSION(nvec), INTENT(INOUT) :: &
                  rho_snow_new              ! water content (snow water equivialent) of snow cover        (-)


   ! soil
   REAL    (KIND = vpp), DIMENSION(nvec,0:ke_soil+1), INTENT(INOUT) :: &
                  t_so_now             ! soil temperature (main level)                (  K  )

   REAL    (KIND = vpp), DIMENSION(nvec,0:ke_soil+1), INTENT(OUT) :: &
                  t_so_new             ! soil temperature (main level)                (  K  )
 
   REAL    (KIND = vpp), DIMENSION(nvec,ke_soil+1), INTENT(INOUT) :: &
                  w_so_now         , & ! total water conent (ice + liquid water)      (m H20)
                  w_so_ice_now         ! ice content                                  (m H20)

   REAL    (KIND = vpp), DIMENSION(nvec,ke_soil+1), INTENT(OUT) :: &
                  w_so_new         , & ! total water conent (ice + liquid water)      (m H20)
                  w_so_ice_new         ! ice content                                  (m H20)
      
   REAL    (KIND = vpp), DIMENSION(ke_soil+1), INTENT(IN) :: &
                  zmls                 ! processing soil level structure

   REAL    (KIND = vpp), DIMENSION(nvec), INTENT(IN) :: &
                  plcov            , & ! fraction of plant cover                         --
                  rootdp               ! depth of the roots                            ( m  )

   INTEGER, DIMENSION(nvec), INTENT(IN) :: &
                  soiltyp_subs         ! type of the soil (keys 0-9)                     --

! ------------------------
! + Local
! ------------------------

   INTEGER      :: &

     ! Indices
     i               , & ! loop index in x-drection
     ksn             , & ! loop index for snow layers
     kso             , & ! loop index for soil layers
     mstyp           , & ! soil type index
     l_top           , &
     counter

   REAL(KIND=vpp)  :: &

     ! Timestep parameters
     zdt           , & ! integration time-step [s]

     ! Meteorological paramters

     zsnow_rate    , & ! rate of snow fall               [kg/m**2 s]
     zrain_rate    , & ! rate of rain fall               [kg/m**2 s]

     zuv           , & ! wind speed                      [m/s]

     l_t0_melt
     ! Snow paramters

   REAL (KIND=vpp)  :: &

     rho_hn     , &    ! new snow density                             [kg/m**3]
     k_ext         ! density dependend extinction coefficient
!     swe_sn   (nvec)   ! local S.W.E calculation               


#ifndef ALLOC_WKARR

   REAL (KIND=vpp), DIMENSION(nvec)  :: &

     t_sn_sfc     , &  ! snow surface temperature                    [K]
     sh_sn        , &  ! sensible heat flux                          [W/m**2]
     lh_sn        , &  ! latent heat flux                            [W/m**2] 
     tch_sn       , &  ! transfer coefficient               
     alpha_sn     , &  ! snow surface albdeo                         [-]
     swnet_sn     , &  ! net short wave radiation                    [W/m**2]
     for_sn       , &  ! total atmospheric forcing at snow surface   [W/m**2]
     runoff_sn    , &  ! total melt runoff
     swe_sn            ! local S.W.E calculation

   INTEGER, DIMENSION(nvec)          :: &

     top              ! index of first (top) snow layer  [-] 

   REAL   (KIND=vpp) ::  &

     hm_sn       (nvec,1:n_layers)  , & ! height (from bottom) of snow layers (main level)     (m)      
     zm_sn       (nvec,1:n_layers)  , & ! depth (from top) of snow layers (main level)         (m)

     rho_sn      (nvec,1:n_layers)  , & ! density of snow layers                               (kg/m**3)
     m_sn        (nvec,1:n_layers)  , & ! layer mass                                           (kg)
     hcap_sn     (nvec,1:n_layers)  , & ! snow layer heat capacity
     hcon_sn     (nvec,1:n_layers)  , & !            heat conductivity
     hdif_sn     (nvec,1:n_layers)  , & !            heat diffusion

     swabs_sn    (nvec,1:n_layers)  , & ! absorbed short wave radiation

     zm_sn_old   (nvec,1:n_layers)  , & ! old value of layer depth
     theta_i_old (nvec,1:n_layers)  , & !              volumetric ice content
     theta_w_old (nvec,1:n_layers)      !              volumetric water content

    LOGICAL, DIMENSION(nvec)  :: &

      melt_flag     , & ! flag indicating melting; IF THEN TRUE
      freeze_flag       !      indicating freezing

    REAL    (KIND=vpp) :: &

      tmp_sn (n_layers)         ! temporary vector      


    REAL    (KIND = vpp)   :: &

      ziw_fr      (nvec,ke_soil+1)   , & ! fractional ice content of soil layer
      zlw_fr      (nvec,ke_soil+1)   , & ! fractional liqu. water content of soil layer
      zw_fr       (nvec,ke_soil+1)   , & ! fractional total water content of soil layers
      zroc        (nvec,ke_soil+1)   , & ! heat capacity of soil layers
      zrocg       (nvec,ke_soil+1)   , & ! total volumetric heat capacity of soil
      zrocg_soil  (nvec,ke_soil+1)   , & !

      zalam       (nvec,ke_soil)     , & ! heat conductivity
      zalamtmp    (nvec,ke_soil)     , & ! heat conductivity
      hzalam      (nvec,ke_soil+1)   , & ! heat conductivity (auxilary variable)

      zfcap    (nvec,ke_soil+1)      , & ! field capacity of soil
      zporv    (nvec,ke_soil+1)      , & ! pore volume (fraction of volume)
      zpwp     (nvec,ke_soil+1)      , & ! plant wilting point (fraction of volume)   
      zdlam    (nvec)                    ! heat conductivity parameterilting point  (fraction of volume)
     
#endif

    REAL     (KIND = vpp)  :: &

      zzz                                ! utility variable

    REAL     (KIND = vpp)  :: &

      zwqg           , & ! mean of fcap and pwp
      z4wdpv             ! 4*zwqg/porv

  REAL (KIND=vpp)       ::  &

    dlw_u_sn              , &  ! derivative of upwelling longwave radiation for snow on the ground    (W m-2)
    dz_up                 , &  ! thickness above the layer of interest                                (m)
    dz_low                , &  ! thickness below the layer of interest                                (m)
    beta

  REAL (KIND=vpp), DIMENSION(-n_layers+1:ke_soil+1) :: &

    zm                   , &  ! depth of main levels

    hcon                 , & ! heat conductivity
    hcap                 , & ! heat capacity
    hdif                 , & ! heat diffusivity

    rho                  , &  ! density

    t_sol                , &  ! temperature of the snow/substrate column

    sw_abs               , &  ! absorbed short-wave radiaton in each layer

    alpha                , &  ! utility variables for building and solving the tri-diagonal matrix
    gamma_sol            , &  !

    a                    , &  !
    b                    , &  !
    c                    , &  !
    d                    , &  !
    e                         ! final snow layer temperature

  REAL    (KIND = vpp) ::  &

    dT_sub_melt             , & ! difference between current temperature and melting temperature
    A_sub_melt              , & ! coefficient A_sub_melt (see notes below)
    dtheta_i       , & ! change in volumetric ice content
    dtheta_w       , & ! change in volumetric water content
    q_mf           , & !
    q_rest

  REAL    (KIND = vpp) :: &

    frac_rho        , &    ! fraction of density water to ice
    limit_theta_i   , &    ! abc-formula
    w_up            , &    ! water content of upper layer
    w_low           , &    ! water content of lower layer
    dtheta_w_up     , &    ! available water in upper layer
    dtheta_w_low    , &    ! available water in lower layer
    dtheta_w_low_x  , &    ! backup variable for dtheta_w_low
    theta_w_bot     , &    ! volumetric water content of bottom snow layer
    excess_water           ! excess water

  REAL    (KIND = vpp), DIMENSION(1:n_layers) :: &

    w_res           , &    ! effective residual water content
    res_wat_cont    , &    ! potential residual water content
    dtheta_w_sub_wtr       ! additional storage capacity due to refreezing

  REAL (KIND = vpp)        :: &

    dL                         , & ! change of layer thickness
    dM                         , & ! change of layer mass
    M                          , & ! initial mass and volmetric content (water of ice)
    hoar                       , & ! hoar mass
    dzm_sn_old                     ! old value of layer thickness


    ! After Andresen (1976)
   REAL(KIND=vpp), PARAMETER ::  &

     c2       = 23.0E-3_vpp              , &  ! Coefficinet [m3/kg]
     eta0     = 9.0E5_vpp                , &  ! The Viscosity Coefficient Eta0 [kg-s/m2]
     c_factor = 0.08_vpp                      ! snow compaction overburden exponential factor (1/K)

   REAL    (KIND = vpp) :: &

     rate_1                              , & ! settling rate for ice loss (s**-1)
     rate_2                              , & ! overburden stress (s**-1)
     overburden                          , & ! overburden load
     tot_rate                            , & ! total settling rates (s**-1)
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

   INTEGER(KIND=vpp) :: j,lay

#ifdef __ICON__
  INTEGER :: my_cart_id
#endif

  INTEGER :: my_thrd_id, mcid, mtid, mbid, mvid, k

  LOGICAL ::  ldebug = .FALSE.

  INTEGER, PARAMETER :: N_IMP = 9999999 !125765 !125761

  real (kind = vpp) :: tmp
! ------------------------------------------------------------------------------
! - End declarations
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Begin Section NULLA: Initialisation(s)
! ------------------------------------------------------------------------------

  ! ------------------------
  ! Intiate local scalars
  ! -----------------------

!   write(*,*) ntstep,prs_gsp,prr_gsp,prs_con,prr_con

   ! Set the debug organizational variables
#ifdef __ICON__
   my_cart_id = get_my_global_mpi_id()
#ifdef _OPENMP
   my_thrd_id = omp_get_thread_num()
#endif
#endif
   mcid = iid_dbg
   mbid = ibl_b
   mvid = isc_b
   mtid =   0

! Just do some checkout prints:
  IF (ldebug) THEN
    IF (iblock == mbid .AND. my_cart_id == mcid) THEN
#ifdef _OPENMP
     IF (my_thrd_id == mtid) THEN
#endif
      WRITE(*,'(A,3I5)'   ) 'SFC-DIAGNOSIS terra start:   ', n_layers, ke_soil
#ifdef _OPENMP
     ENDIF
#endif
    ENDIF

    IF (iblock == mbid .AND. my_cart_id == mcid) THEN
    !$acc update host(whatever is printed)
    ENDIF

    DO i = ivstart, ivend
      IF (i== mvid .AND. iblock == mbid .AND. my_cart_id == mcid) THEN
#ifdef _OPENMP
       IF (my_thrd_id == mtid) THEN
#endif
        WRITE(*,'(A,2I5)'  ) 'SFC-DIAGNOSIS snow:   iblock = ', iblock, i
        WRITE(*,'(A      )') ' External Parameters:  '
        WRITE(*,'(A,I28  )') '   soiltyp          :  ', soiltyp_subs(i)
        WRITE(*,'(A,F28.16)') '   plcov            :  ', plcov       (i)
        WRITE(*,'(A,F28.16)') '   rootdp           :  ', rootdp      (i)

        WRITE(*,'(A      )') ' Single level parameters:'
        WRITE(*,'(A,F28.16)') '   u     ke         :  ', u           (i)
        WRITE(*,'(A,F28.16)') '   v     ke         :  ', v           (i)
        WRITE(*,'(A,F28.16)') '   t     ke         :  ', t           (i)
        WRITE(*,'(A,F28.16)') '   qv    ke         :  ', qv          (i)
        WRITE(*,'(A,F28.16)') '   ps               :  ', ps          (i)
        WRITE(*,'(A,F28.16)') '   prr_con          :  ', prr_con     (i)
        WRITE(*,'(A,F28.16)') '   prs_con          :  ', prs_con     (i)
        WRITE(*,'(A,F28.16)') '   conv_frac        :  ', conv_frac   (i)
        WRITE(*,'(A,F28.16)') '   prr_gsp          :  ', prr_gsp     (i)
        WRITE(*,'(A,F28.16)') '   prs_gsp          :  ', prs_gsp     (i)
        WRITE(*,'(A,F28.16)') '   prg_gsp          :  ', prg_gsp     (i)
#ifdef TWOMOM_SB
        WRITE(*,'(A,F28.16)') '   prh_gsp          :  ', prh_gsp     (i)
#endif
        WRITE(*,'(A,F28.16)') '   sobs             :  ', sobs        (i)
        WRITE(*,'(A,F28.16)') '   thbs             :  ', thbs        (i)
!        WRITE(*,'(A,F28.16)') '   swdir_s          :  ', swdir_s     (i)
!        WRITE(*,'(A,F28.16)') '   swdifd_s         :  ', swdifd_s    (i)
!        WRITE(*,'(A,F28.16)') '   swdifu_s         :  ', swdifu_s    (i)
!        WRITE(*,'(A,F28.16)') '   lwd_s            :  ', lwd_s       (i)
!        WRITE(*,'(A,F28.16)') '   lwu_s            :  ', lwu_s       (i)
        WRITE(*,'(A,F28.16)') '   t_snow_now       :  ', t_snow_now  (i)
        WRITE(*,'(A,F28.16)') '   h_snow           :  ', h_snow      (i)
        WRITE(*,'(A,F28.16)') '   w_snow_now       :  ', w_snow_now  (i)
        WRITE(*,'(A,F28.16)') '   hn_sn_now        :  ', hn_sn_now   (i)
        WRITE(*,'(A,F28.16)') '   top_sn_now       :  ', top_sn_now  (i)

        WRITE(*,'(A      )') ' Multi level parameters:'
do k = 1, n_layers
        WRITE(*,'(A,I1,A,F28.16)') '   t_sn_now   (',k,')   :  ', t_sn_now    (i,k)
enddo
do k = 1, n_layers
        WRITE(*,'(A,I1,A,F28.16)') '   theta_i_now(',k,')   :  ', theta_i_now (i,k)
enddo
do k = 1, n_layers
        WRITE(*,'(A,I1,A,F28.16)') '   theta_w_now(',k,')   :  ', theta_w_now (i,k)
enddo
do k = 1, n_layers
        WRITE(*,'(A,I1,A,F28.16)') '   theta_a_now(',k,')   :  ', theta_a_now (i,k)
enddo
do k = 1, n_layers
        WRITE(*,'(A,I1,A,F28.16)') '   dzm_sn_now (',k,')   :  ', dzm_sn_now  (i,k)
enddo
do k = 0, ke_soil+1
        WRITE(*,'(A,I1,A,F28.16)') '   t_so_now(',k,')      :  ', t_so_now    (i,k)
enddo
do k = 1, ke_soil+1
        WRITE(*,'(A,I1,A,F28.16)') '   w_so_now(',k,')      :  ', w_so_now    (i,k)
enddo
do k = 1, ke_soil+1
        WRITE(*,'(A,I1,A,F28.16)') '   w_so_ico_now(',k,')  :  ', w_so_ice_now(i,k)
enddo
do k = 1, ke_soil+1
        WRITE(*,'(A,I1,A,F28.16)') '   zmls        (',k,')  :  ', zmls        (k)
enddo
#ifdef _OPENMP
       ENDIF
#endif
      ENDIF
    ENDDO
  ENDIF


   hdif_sn = 0.0_vpp
   ! time step for soil variables
   zdt      = dt

   !Subroutine arguments

   !$acc data                                                                   &
   !$acc present(u, v, t, qv, ps, prr_con, prs_con, conv_frac)                  &
   !$acc present(prr_gsp, prs_gsp, prg_gsp)                                     &
#ifdef TWOMOM_SB
   !$acc present(prh_gsp)                                                       &
#endif
!   !$acc present(swdir_s, swdifd_s, swdifu_s, lwd_s, lwu_s)                     &
   !$acc present(t_snow_now, t_snow_new, zshfl_snow, zlhfl_snow)                &
   !$acc present(t_sn_now, t_sn_new, theta_i_now, theta_i_new)                  &
   !$acc present(theta_w_now, theta_w_new, theta_a_new, theta_a_now)            &
   !$acc present(dzm_sn_now, dzm_sn_new, hn_sn_now, hn_sn_new)                  &
   !$acc present(top_sn_now, top_sn_new, h_snow, t_so_now, t_so_new)            &
   !$acc present(w_so_ice_now, w_so_ice_new, w_so_now, w_so_new)                &         
   !$acc present(zmls, plcov, rootdp, soiltyp_subs)                             & !SB: Not sure if global constant fields needed here
   !$acc present(w_snow_new)

   ! Local arrays
   !$acc data                                                                   &
   !$acc present (t_sn_sfc, sh_sn, lh_sn, tch_sn, alpha_sn, swnet_sn, for_sn,swe_sn)   &
   !$acc present (runoff_sn, top, hm_sn, zm_sn, rho_sn, m_sn, hcap_sn, hcon_sn) &
   !$acc present (hdif_sn, swabs_sn, zm_sn_old, theta_i_old, theta_w_old)       &
   !$acc present (melt_flag, freeze_flag, tmp_sn)                               &
   !$acc present (ziw_fr, zlw_fr, zw_fr, zroc, zrocg, zrocg_soil)               &
   !$acc present (zalam, zalamtmp, hzalam, zfcap, zporv, zpwp, zdlam)


  ! ------------------------
  ! Intiate local fields
  ! -----------------------

   !$acc parallel async default(none)
   !$acc loop gang vector
   DO i = ivstart, ivend

     top(i) = NINT(top_sn_now(i)) 

   ENDDO
   !$acc end parallel

   !write(*,*) 'snowpolino start: ',top(1)


   !$acc parallel async default(none)
   !$acc loop gang vector
   DO i = ivstart, ivend

     t_sn_sfc(i) = t0_melt

     lh_sn(i)    = 0.0_vpp
     sh_sn(i)    = 0.0_vpp

     tch_sn(i)   = 0.0_vpp

     alpha_sn(i) = 0.85_vpp

     swnet_sn(i) = 0.0_vpp

     for_sn(i)   = 0.0_vpp

     swe_sn(i)   = 0.0_vpp

   ENDDO
   !$acc end parallel




   !$acc parallel async default(none)
   !$acc loop gang vector 
     DO i = ivstart, ivend
       !$acc loop seq
       DO ksn = 1, n_layers

       hm_sn   (i,ksn) = 0.0_vpp
       zm_sn   (i,ksn) = 0.0_vpp
       rho_sn  (i,ksn) = 0.0_vpp
       m_sn    (i,ksn) = 0.0_vpp
       hcon_sn (i,ksn) = 0.0_vpp
       hcap_sn (i,ksn) = 0.0_vpp

       swabs_sn(i,ksn) = 0.0_vpp

     ENDDO
   ENDDO
   !$acc end parallel

   !$acc parallel async
   !$acc loop gang vector private(l_top)
   DO i = ivstart, ivend

     l_top = top(i)
     IF(l_top .GE. 1) THEN

        ! Update dependent variables
        !CALL update(top(i), hm_sn(i,:), zm_sn(i,:), dzm_sn_now(i,:), m_sn(i,:)     , &
        !theta_i_now(i,:), theta_w_now(i,:), theta_a_now(i,:), rho_sn(i,:)          , &
        !hcap_sn(i,:), hcon_sn(i,:))

  ! -------------------------
  ! Height of snow (main) levels
  ! -------------------------
        !$acc loop seq
        do ksn=1,l_top
           if(ksn .eq. 1) then
              hm_sn(i,ksn) = dzm_sn_now(i,ksn)
           else
              hm_sn(i,ksn) = hm_sn(i,ksn-1) + dzm_sn_now(i,ksn)
           endif
        enddo  

  ! --------------------------
  ! Depth of snow layer (main) levels
  ! -------------------------
        !$acc loop seq
        do ksn = l_top,1,-1
           zm_sn(i, (l_top+1) - ksn ) = hm_sn(i,ksn) !invert height vector
           theta_a_now(i,ksn) = max(0.0_vpp,1.0_vpp - theta_i_now(i,ksn) - theta_w_now(i,ksn))
        enddo

  ! --------------------------
  ! Snow layer density
  ! -------------------------
         !$acc loop seq
         do ksn = 1, l_top, 1
            if(theta_i_now(i,ksn) .eq. 0.0_vpp) then
               rho_sn(i,ksn) = 0.0_vpp
            else
               rho_sn(i,ksn) = theta_i_now(i,ksn)*rho_i + theta_w_now(i,ksn)*rho_w
            endif
         enddo
         
  ! --------------------------
  ! Heat capacity
  ! --------------------------
   !$acc loop seq
   DO ksn = 1, l_top, 1
       IF(rho_sn(i,ksn) .LT. eps) THEN
         hcap_sn(i,ksn) = 0.0_vpp
       ELSE
         hcap_sn(i,ksn) = (  rho_a   * theta_a_now(i,ksn) * specific_heat_air     &
                           + rho_i   * theta_i_now(i,ksn) * specific_heat_ice     &
                           + rho_w   * theta_w_now(i,ksn) * specific_heat_water)  &
                            / rho_sn(i,ksn)
       ENDIF
   ENDDO

  ! --------------------------
  ! Heat conductivity
  ! --------------------------
   !$acc loop seq
   DO ksn = 1, l_top, 1
      hcon_sn(i,ksn) = 2.22_vpp * EXP(1.88_vpp * LOG(rho_sn(i,ksn)/rho_i))
   ENDDO

  ! --------------------------
  ! Snow layer mass
  ! --------------------------
   !$acc loop seq
   DO ksn = 1, l_top, 1
     IF(ksn .EQ. 1) THEN
        m_sn(i,ksn) = zm_sn(i,ksn) * ((theta_i_now(i,ksn) * rho_i) + (theta_w_now(i,ksn) * rho_w))
      ELSE
        m_sn(i,ksn) = ABS(zm_sn(i,ksn) - zm_sn(i,ksn-1))  * ((theta_i_now(i,ksn) * rho_i) + (theta_w_now(i,ksn) * rho_w))
      ENDIF
   ENDDO




      ENDIF

   ENDDO
   !$acc end parallel


   !write(*,*) 'start: ',ntstep, 

! ------------------------------------------------------------------------------
! - End Section NULLA
! ------------------------------------------------------------------------------

!    ! Do some printing to check whether or not values look ok 
!    DO i  = ivstart, ivend
!       print *, top_sn_now(i),  NINT(top_sn_now(i)), t_sn_now(i,1) 
!       print *, zm_sn(i,:)
!       print *, hm_sn(i,:)
!       print *, top(i)
!    ENDDO

! ------------------------------------------------------------------------------
! + Begin Section I: New snow
! ------------------------------------------------------------------------------


  !$acc parallel async
  !$acc loop gang vector private (zsnow_rate, zuv, tmp_sn, l_top,l_t0_melt)
  DO i = ivstart, ivend

    l_top = top(i)  
  ! ----------------------
  ! Calculate precipitation rate solid/liquid
  ! ----------------------

    IF ( nclass_gscp >= 2000 ) THEN
      ! only possible when running 2-moment microphysics
#ifdef TWOMOM_SB
      zsnow_rate = prs_gsp(i)+prs_con(i)+prg_gsp(i)+prh_gsp(i) ! [kg/m**2 s]
#endif
    ELSEIF ( nclass_gscp >= 6 ) THEN
      zsnow_rate = prs_gsp(i)+prs_con(i)+prg_gsp(i)            ! [kg/m**2 s]
    ELSE
      zsnow_rate = prs_gsp(i)+prs_con(i)                       ! [kg/m**2 s]
    ENDIF

      zrain_rate = prr_gsp(i)+prr_con(i)  ! [kg/m**2 s]

      !write (*,*) 'new_snow: ',zsnow_rate
  ! ----------------------
  ! Calculate wind speed
  ! ----------------------

    zuv        = SQRT ( u(i)**2 + v(i)**2 )

  ! ----------------------
  ! Calculate new snow amounts
  ! ----------------------



    IF(zsnow_rate .GT. 0.0_vpp) THEN ! there should be solid precipitation

      rho_hn = 0.0_vpp
  
      ! Calcualte potential new snow density
      CALL new_snow_density(rho_hn, t(i), zuv)      

      ! Calculate new snow amounts
       hn_sn_now(i) =  hn_sn_now(i) + ( (zsnow_rate * zdt/rho_w) * rho_w/rho_hn )

       !write(*,*) 'new snow: ', zsnow_rate,zdt,rho_hn, ( (zsnow_rate * zdt/rho_w) * rho_w/rho_hn )

       !write(*,*) hn_sn_now,rho_hn,zsnow_rate,zdt,rho_w
      ! ----------------------
      ! Check if there is enough snow available for new layers, if so take  action
      ! ----------------------

      IF(hn_sn_now(i) .GE. max_height_layer) THEN  ! enough snow for a full layer

         !write(*,*) 'ASFAFS'

        IF(l_top .LT. n_layers) THEN ! snow cover thinner than n_layers*max_snow_height


            !write(*,*) 'GGGGGGGGGGGG' 

             DO WHILE (hn_sn_now(i) .GT. max_height_layer .AND. l_top .LT. n_layers)

             ! Limit hn to max_height_layers
             hn_sn_now(i)  = MIN(hn_sn_now(i), max_height_layer)

             ! Get layer indices
             l_top = l_top + 1  ! new index of top layer

             !write(*,*) ' OH MY A: ', l_top
             ! Assign properties to new snow layers
             !CALL new_snow_layers(top(i), hn_sn_now(i),dzm_sn_now(i,:), t_sn_now(i,:), t(i),  &
             !                   theta_i_now(i,:), theta_w_now(i,:), theta_a_now(i,:),               &
             !                   rho_sn(i,:), rho_hn)

             l_t0_melt              = real(t0_melt,vpp)
             dzm_sn_now(i,l_top)  = hn_sn_now(i)
             rho_sn(i,l_top)      = rho_hn
             theta_i_now(i,l_top) = rho_hn / rho_i
             theta_w_now(i,l_top) = 0.0_vpp
             theta_a_now(i,l_top) = 1.0_vpp - theta_i_now(i,l_top) - theta_w_now(i,l_top)
             t_sn_now(i,l_top)    = min(t(i),l_t0_melt)

              !!! replacing Update function call
              ! -------------------------
              ! Height of snow (main) levels
              ! -------------------------
                    !$acc loop seq
                    do ksn=1,l_top
                       if(ksn .eq. 1) then
                          hm_sn(i,ksn) = dzm_sn_now(i,ksn)
                       else
                          hm_sn(i,ksn) = hm_sn(i,ksn-1) + dzm_sn_now(i,ksn)
                       endif
                    enddo  
      
              ! --------------------------
              ! Depth of snow layer (main) levels
              ! -------------------------
                    !$acc loop seq
                    do ksn = l_top,1,-1
                       zm_sn(i, (l_top+1) - ksn ) = hm_sn(i,ksn) !invert height vector
                       theta_a_now(i,ksn) = max(0.0_vpp,1.0_vpp - theta_i_now(i,ksn) - theta_w_now(i,ksn))
                    enddo
      
              ! --------------------------
              ! Snow layer density
              ! -------------------------
                     !$acc loop seq
                     do ksn = 1, l_top, 1
                        if(theta_i_now(i,ksn) .eq. 0.0_vpp) then
                           rho_sn(i,ksn) = 0.0_vpp
                        else
                           rho_sn(i,ksn) = theta_i_now(i,ksn)*rho_i + theta_w_now(i,ksn)*rho_w
                        endif
                     enddo
                     
              ! --------------------------
              ! Heat capacity
              ! --------------------------
               !$acc loop seq
               DO ksn = 1, l_top, 1
                   IF(rho_sn(i,ksn) .LT. eps) THEN
                     hcap_sn(i,ksn) = 0.0_vpp
                   ELSE
                     hcap_sn(i,ksn) = (  rho_a   * theta_a_now(i,ksn) * specific_heat_air     &
                                       + rho_i   * theta_i_now(i,ksn) * specific_heat_ice     &
                                       + rho_w   * theta_w_now(i,ksn) * specific_heat_water)  &
                                        / rho_sn(i,ksn)
                   ENDIF
               ENDDO
      
              ! --------------------------
              ! Heat conductivity
              ! --------------------------
               !$acc loop seq
               DO ksn = 1, l_top, 1
                  hcon_sn(i,ksn) = 2.22_vpp * EXP(1.88_vpp * LOG(rho_sn(i,ksn)/rho_i))
               ENDDO
      
              ! --------------------------
              ! Snow layer mass
              ! --------------------------
               !$acc loop seq
               DO ksn = 1, l_top, 1
                 IF(ksn .EQ. 1) THEN
                    m_sn(i,ksn) = zm_sn(i,ksn) * ((theta_i_now(i,ksn) * rho_i) + (theta_w_now(i,ksn) * rho_w))
                  ELSE
                    m_sn(i,ksn) = ABS(zm_sn(i,ksn) - zm_sn(i,ksn-1))  * ((theta_i_now(i,ksn) * rho_i) + (theta_w_now(i,ksn) * rho_w))
                  ENDIF
               ENDDO

             ! Update dependent variables
             !CALL update(top(i), hm_sn(i,:), zm_sn(i,:), dzm_sn_now(i,:), m_sn(i,:)                 , &
             !            theta_i_now(i,:), theta_w_now(i,:), theta_a_now(i,:), rho_sn(i,:)          , &
             !            hcap_sn(i,:), hcon_sn(i,:))

             ! Reset new snow storage
             hn_sn_now(i) = hn_sn_now(i) - max_height_layer

             ! Update a few values NOTE this could maybe also happen in update()
             h_snow(i)   = hm_sn(i,l_top)
             t_sn_sfc(i) = t_sn_now(i,l_top)

             ENDDO


        ELSE ! snow cover thicker than n_layers*max_snow_height


           ! --------------------
           ! Merge/Add layers
           ! --------------------

           DO WHILE (hn_sn_now(i) .GT. max_height_layer)

            ! Limit hn to max_height_layers
            hn_sn_now(i)  = MIN(hn_sn_now(i), max_height_layer)

            ! Merge bottom two layers

            theta_i_now(i,1) =   (theta_i_now(i,2)*dzm_sn_now(i,2)                                        &
                               +  theta_i_now(i,1)*dzm_sn_now(i,1)) / (dzm_sn_now(i,2) + dzm_sn_now(i,1))    ! volumetric ice content

            theta_w_now(i,1) =   (theta_w_now(i,2)*dzm_sn_now(i,2)                                        &
                               +  theta_w_now(i,1)*dzm_sn_now(i,1)) / (dzm_sn_now(i,2) + dzm_sn_now(i,1))    ! volumetric water content

            t_sn_now(i,1)    =   (t_sn_now(i,2)*dzm_sn_now(i,2)                                           &
                               +  t_sn_now(i,1)*dzm_sn_now(i,1))    / (dzm_sn_now(i,2) + dzm_sn_now(i,1))    ! snow layer temperature


            dzm_sn_now(i,1)  = dzm_sn_now(i,2) + dzm_sn_now(i,1)                                             ! layer thickness
            
            ! Update top layer index
            l_top = l_top - 1
            !write(*,*) ' OH MY B: ', l_top

             ! ------------------
             ! Move layers down
             ! -----------------

               ! Layer thickness
               ! ----------------------
               tmp_sn(1:n_layers) = dzm_sn_now(i,1:n_layers)

               DO ksn = n_layers, 3, -1
                 dzm_sn_now(i,ksn-1)  = tmp_sn(ksn)
               ENDDO


               ! Volumetric ice content
               ! ----------------------
               tmp_sn(1:n_layers) = theta_i_now(i,1:n_layers)

               DO ksn = n_layers, 3, -1
                 theta_i_now(i,ksn-1)  = tmp_sn(ksn)
               ENDDO


               ! Volumetric water content
               ! ----------------------
               tmp_sn(1:n_layers) = theta_w_now(i,1:n_layers)

               DO ksn = n_layers, 3, -1
                 theta_w_now(i,ksn-1)  = tmp_sn(ksn)
               ENDDO


               ! Layer temperature
               ! ----------------------
               tmp_sn(1:n_layers) = t_sn_now(i,1:n_layers)

               DO ksn = n_layers, 3, -1
                 t_sn_now(i,ksn-1)  = tmp_sn(ksn)
               ENDDO

             ! -------------------
             ! Add new snow layer and assign properties
             ! -------------------

               ! limit top layer index it can/should only be n_layers here
               l_top = MIN(l_top + 1 , n_layers)                         
               !write(*,*) ' OH MY C: ', l_top

               ! Assign properties
               dzm_sn_now(i,l_top)   = hn_sn_now(i)

               rho_sn(i,l_top)       = rho_hn

               theta_i_now(i,l_top)  = rho_hn / rho_i
               theta_w_now(i,l_top)  = 0.0_vpp           ! new snow is always dry
               theta_a_now(i,l_top)  = 1.0_vpp - theta_i_now(i,l_top) - theta_w_now(i,l_top)

               t_sn_now(i,l_top)     = MIN(t(i), t0_melt)

                ! -------------------------
                ! Height of snow (main) levels
                ! -------------------------
                      !$acc loop seq
                      do ksn=1,l_top
                         if(ksn .eq. 1) then
                            hm_sn(i,ksn) = dzm_sn_now(i,ksn)
                         else
                            hm_sn(i,ksn) = hm_sn(i,ksn-1) + dzm_sn_now(i,ksn)
                         endif
                      enddo  
       
                ! --------------------------
                ! Depth of snow layer (main) levels
                ! -------------------------
                      !$acc loop seq
                      do ksn = l_top,1,-1
                         zm_sn(i, (l_top+1) - ksn ) = hm_sn(i,ksn) !invert height vector
                         theta_a_now(i,ksn) = max(0.0_vpp,1.0_vpp - theta_i_now(i,ksn) - theta_w_now(i,ksn))
                      enddo
       
                ! --------------------------
                ! Snow layer density
                ! -------------------------
                       !$acc loop seq
                       do ksn = 1, l_top, 1
                          if(theta_i_now(i,ksn) .eq. 0.0_vpp) then
                             rho_sn(i,ksn) = 0.0_vpp
                          else
                             rho_sn(i,ksn) = theta_i_now(i,ksn)*rho_i + theta_w_now(i,ksn)*rho_w
                          endif
                       enddo
                       
                ! --------------------------
                ! Heat capacity
                ! --------------------------
                 !$acc loop seq
                 DO ksn = 1, l_top, 1
                     IF(rho_sn(i,ksn) .LT. eps) THEN
                       hcap_sn(i,ksn) = 0.0_vpp
                     ELSE
                       hcap_sn(i,ksn) = (  rho_a   * theta_a_now(i,ksn) * specific_heat_air     &
                                         + rho_i   * theta_i_now(i,ksn) * specific_heat_ice     &
                                         + rho_w   * theta_w_now(i,ksn) * specific_heat_water)  &
                                          / rho_sn(i,ksn)
                     ENDIF
                 ENDDO
       
                ! --------------------------
                ! Heat conductivity
                ! --------------------------
                 !$acc loop seq
                 DO ksn = 1, l_top, 1
                    hcon_sn(i,ksn) = 2.22_vpp * EXP(1.88_vpp * LOG(rho_sn(i,ksn)/rho_i))
                 ENDDO
       
                ! --------------------------
                ! Snow layer mass
                ! --------------------------
                 !$acc loop seq
                 DO ksn = 1, l_top, 1
                   IF(ksn .EQ. 1) THEN
                      m_sn(i,ksn) = zm_sn(i,ksn) * ((theta_i_now(i,ksn) * rho_i) + (theta_w_now(i,ksn) * rho_w))
                    ELSE
                      m_sn(i,ksn) = ABS(zm_sn(i,ksn) - zm_sn(i,ksn-1))  * ((theta_i_now(i,ksn) * rho_i) + (theta_w_now(i,ksn) * rho_w))
                    ENDIF
                 ENDDO

               ! Update dependent variables
               !CALL update(top(i), hm_sn(i,:), zm_sn(i,:), dzm_sn_now(i,:), m_sn(i,:)                 , &
               !            theta_i_now(i,:), theta_w_now(i,:), theta_a_now(i,:), rho_sn(i,:)          , &
               !            hcap_sn(i,:), hcon_sn(i,:))

            ! Reset new snow storage
            hn_sn_now(i) = hn_sn_now(i) - max_height_layer

            ! Update a few values NOTE this could maybe also happen in update()
            h_snow(i)   = hm_sn(i,l_top)
            t_sn_sfc(i) = t_sn_now(i,l_top)

           ENDDO ! end do while



         ENDIF ! end thin or thick snow cover


      ENDIF ! end full layer check


    ELSE ! no snowfall check if bottom layer can be splitted

     IF(l_top .LT. n_layers .AND. dzm_sn_now(i,1) .GT. (2.0_vpp*max_height_layer)) THEN ! bottoom layer can be split


       DO WHILE (l_top .LT. n_layers .AND. dzm_sn_now(i,1) .GT. (2.0_vpp*max_height_layer))

         ! -------------------
         ! Move layers up
         ! -------------------

          ! Layer thickness
          tmp_sn(1:n_layers) = dzm_sn_now(i,1:n_layers)

          !$acc loop seq
          DO ksn = 2, n_layers-1, 1
            dzm_sn_now(i,ksn+1) = tmp_sn(ksn)
          ENDDO

          ! Volumetric ice content
          tmp_sn(1:n_layers) = theta_i_now(i,1:n_layers)

          !$acc loop seq
          DO ksn = 2, n_layers-1, 1
            theta_i_now(i,ksn+1) = tmp_sn(ksn)
          ENDDO

          ! Volumetric water content
          tmp_sn(1:n_layers) = theta_w_now(i,1:n_layers)

          !$acc loop seq
          DO ksn = 2, n_layers-1, 1
            theta_w_now(i,ksn+1) = tmp_sn(ksn)
          ENDDO

          ! Layer temperature
          tmp_sn(1:n_layers) = t_sn_now(i,1:n_layers)

          !$acc loop seq
          DO ksn = 2, n_layers-1, 1
            t_sn_now(i,ksn+1) = tmp_sn(ksn)
          ENDDO

         ! ---------------------
         ! Split bottom layer
         ! --------------------

          dzm_sn_now(i,1)  = dzm_sn_now(i,1) - max_height_layer
          dzm_sn_now(i,2)  = max_height_layer

          theta_i_now(i,2) = theta_i_now(i,1)
          theta_w_now(i,2) = theta_w_now(i,1)

          t_sn_now(i,2)    = t_sn_now(i,1)

         ! Update dependent variables
         l_top = l_top + 1
!         write(*,*) ' OH MY D: ', l_top


              !!! replacing Update function call
              ! -------------------------
              ! Height of snow (main) levels
              ! -------------------------
                    !$acc loop seq
                    do ksn=1,l_top
                       if(ksn .eq. 1) then
                          hm_sn(i,ksn) = dzm_sn_now(i,ksn)
                       else
                          hm_sn(i,ksn) = hm_sn(i,ksn-1) + dzm_sn_now(i,ksn)
                       endif
                    enddo  
      
              ! --------------------------
              ! Depth of snow layer (main) levels
              ! -------------------------
                    !$acc loop seq
                    do ksn = l_top,1,-1
                       zm_sn(i, (l_top+1) - ksn ) = hm_sn(i,ksn) !invert height vector
                       theta_a_now(i,ksn) = max(0.0_vpp,1.0_vpp - theta_i_now(i,ksn) - theta_w_now(i,ksn))
                    enddo
      
              ! --------------------------
              ! Snow layer density
              ! -------------------------
                     !$acc loop seq
                     do ksn = 1, l_top, 1
                        if(theta_i_now(i,ksn) .eq. 0.0_vpp) then
                           rho_sn(i,ksn) = 0.0_vpp
                        else
                           rho_sn(i,ksn) = theta_i_now(i,ksn)*rho_i + theta_w_now(i,ksn)*rho_w
                        endif
                     enddo
                     
              ! --------------------------
              ! Heat capacity
              ! --------------------------
               !$acc loop seq
               DO ksn = 1, l_top, 1
                   IF(rho_sn(i,ksn) .LT. eps) THEN
                     hcap_sn(i,ksn) = 0.0_vpp
                   ELSE
                     hcap_sn(i,ksn) = (  rho_a   * theta_a_now(i,ksn) * specific_heat_air     &
                                       + rho_i   * theta_i_now(i,ksn) * specific_heat_ice     &
                                       + rho_w   * theta_w_now(i,ksn) * specific_heat_water)  &
                                        / rho_sn(i,ksn)
                   ENDIF
               ENDDO
      
              ! --------------------------
              ! Heat conductivity
              ! --------------------------
               !$acc loop seq
               DO ksn = 1, l_top, 1
                  hcon_sn(i,ksn) = 2.22_vpp * EXP(1.88_vpp * LOG(rho_sn(i,ksn)/rho_i))
               ENDDO
      
              ! --------------------------
              ! Snow layer mass
              ! --------------------------
               !$acc loop seq
               DO ksn = 1, l_top, 1
                 IF(ksn .EQ. 1) THEN
                    m_sn(i,ksn) = zm_sn(i,ksn) * ((theta_i_now(i,ksn) * rho_i) + (theta_w_now(i,ksn) * rho_w))
                  ELSE
                    m_sn(i,ksn) = ABS(zm_sn(i,ksn) - zm_sn(i,ksn-1))  * ((theta_i_now(i,ksn) * rho_i) + (theta_w_now(i,ksn) * rho_w))
                  ENDIF
               ENDDO

         !CALL update(top(i), hm_sn(i,:), zm_sn(i,:), dzm_sn_now(i,:), m_sn(i,:)         , &
         !            theta_i_now(i,:), theta_w_now(i,:), theta_a_now(i,:), rho_sn(i,:)  , &
         !            hcap_sn(i,:), hcon_sn(i,:))

         ! Update a few values NOTE this could maybe also happen in update()
         h_snow(i)   = hm_sn(i,l_top)
         t_sn_sfc(i) = t_sn_now(i,l_top)

       ENDDO ! end do while

       !if(ntstep .eq. N_IMP) then
       !     write(*,*) 'I AM HERE - split bottom layer'
       !endif


     ENDIF ! split bottom layer


    ENDIF ! snowfall or no snow that was the question

  top(i) = l_top      
  ENDDO
  !$acc end parallel

! ------------------------------------------------------------------------------
! - End Section I: New snow
! ------------------------------------------------------------------------------
!if(ntstep .eq. 165052) then
!   
!   write(*,*) 'end section I: ',h_snow(1),hm_sn(1,top(1))
!
!endif
if(ntstep .eq. N_IMP) then
   do ksn = 1,top(1)
     write(*,*) 'end section I: ',ksn,theta_i_now(1,ksn),theta_a_now(1,ksn),theta_w_now(1,ksn),t_sn_now(1,ksn)

   enddo
endif



! =============================================================================
! + Begin Section II: Atmospheric forcing & Heat Equation 
! =============================================================================



  ! --------------------
  ! Some updating
  ! --------------------

    !$acc parallel async default(none)
    !$acc loop gang vector
    DO i = ivstart, ivend

      IF(top(i) .GE. 1) THEN

        t_sn_sfc(i) = t_sn_now(i,top(i))

      ENDIF


    ENDDO
    !$acc end parallel


  ! ----------------------------------------------------------------------------
  ! + Surface fluxes
  ! ---------------------------------------------------------------------------
  
    !$acc parallel async default(none)
    !$acc loop gang vector private (zuv)
    DO i = ivstart, ivend

      IF(top(i) .GE. 1 ) THEN   ! snow on the ground

          zuv        = SQRT ( u(i)**2 + v(i)**2 )        ! wind speed
      
        ! Calculate transfer coefficient
        ! ---------------  

          CALL calc_tch(t_sn_sfc(i), t(i), tch_sn(i),  zuv, qv(i), h_snow(i), ps(i))


        ! Calculate turbulent fluxes
        ! ---------------

          CALL turb_flux(sh_sn(i), lh_sn(i), ps(i), t(i), t_sn_sfc(i), qv(i), tch_sn(i), zuv)


        ! Calculate radiative fluxes
        ! ---------------

        !  CALL rad_flux(swdir_s(i), swdifd_s(i), swdifu_s(i), lwd_s(i), lwu_s(i), alpha_sn(i), t(i), t_sn_sfc(i))

      ENDIF

    ENDDO
    !$acc end parallel

  ! ----------------------------------------------------------------------------
  ! + Atmospheric forcing
  ! ----------------------------------------------------------------------------


          !$acc parallel async default(none)
          !$acc loop gang vector
          DO i = ivstart, ivend

            IF(top(i) .GE. 1) THEN ! snow on the ground

              ! Calculate net short wave radiation
              ! -----------------------------------

              swnet_sn(i) = sobs(i) !(swdir_s(i) + swdifd_s(i))  - swdifu_s(i)

              ! Distribute absorbed short wave radiation across layers
              ! -----------------------------------

              IF(top(i) .GE. 1) THEN ! sufficient snow cover

                !$acc loop seq
                DO ksn = top(i), 1, -1

                   k_ext = rho_sn(i,ksn) / 3.0_vpp + 50.0_vpp

                   swabs_sn(i,ksn) = swnet_sn(i) * (1.0_vpp - exp(-k_ext * dzm_sn_now(i,ksn)))

                   swnet_sn(i) = swnet_sn(i) - swabs_sn(i,ksn)

                   ! Put the remaining energy into the bottom layer
                   IF(swnet_sn(i) .LT. 5E-4_vpp) THEN

                     swabs_sn(i,1) = swnet_sn(i)

                   ENDIF

                 ENDDO

               ENDIF

              ! Calculate total atmospheric forcing
              ! ------------------------------------
 
!              for_sn(i)  = swabs_sn(i,top(i)) + (lwd_s(i) - lwu_s(i)) + lh_sn(i) + sh_sn(i)
               for_sn(i)  = swabs_sn(i,top(i)) + thbs(i) + lh_sn(i) + sh_sn(i)
       
             ENDIF

           ENDDO
          !$acc end parallel

          !if( top (1) .ge. 1 ) then
          !  write(*,*) 'eb: ',ntstep,for_sn(1),swabs_sn(1,top(1)),thbs(1),lh_sn(1),sh_sn(1) 
          !else
          !  write(*,*) 'eb: ',ntstep,for_sn(1),0.0,thbs(1),lh_sn(1),sh_sn(1)
          !endif 

  ! ----------------------------------------------------------------------------
  ! + Solve the heat equation
  ! ----------------------------------------------------------------------------

     ! ---------------------
     ! Prepare soil properties
     ! ---------------------

     ! Initiations
     ! ------------------

     !$acc parallel async
     !$acc loop gang vector private(mstyp)
     DO i = ivstart, ivend
     
       mstyp     = soiltyp_subs(i)        ! soil type

       zrocg     (i,:) = crhoc (mstyp)              ! heat capacity
       zrocg_soil(i,:) = crhoc (mstyp)              ! heat capacity
       zalam     (i,:) = cala0 (mstyp)              ! heat conductivity parameter
       
       zporv     (i,:) = cporv (mstyp)              ! pore volume
       zfcap     (i,:) = cfcap (mstyp)              ! field capacity     
       zpwp      (i,:) = cpwp  (mstyp)              ! plant wilting point
       zdlam     (i)   = cala1 (mstyp)-cala0(mstyp) ! heat conductivity parameter
     ENDDO
     !$acc end parallel




    
     ! Heat conductivity; NOTE: This is currently only itype_heatcond == 1
     ! -------------------

     !$acc parallel async
     !$acc loop gang vector collapse(2)
     DO kso = 1, ke_soil
       DO i = ivstart, ivend
         zwqg         = 0.5_vpp*(zfcap(i,kso) + zpwp(i,kso))
         z4wdpv       = 4.0_vpp*zwqg/zporv(i,kso)
         ! heat conductivity
         zalamtmp(i,kso) =              zdlam(i)                         &
                               * (0.25_vpp + 0.30_vpp*zdlam(i)           &
                               / (1.0_vpp+0.75_vpp*zdlam(i)))            &
                               * MIN (z4wdpv, 1.0_vpp + (z4wdpv-1.0_vpp) &
                               *(1.0_vpp+0.35_vpp*zdlam(i))              &
                               /(1.0_vpp+1.95_vpp*zdlam(i)))
       ENDDO
     ENDDO
      !$acc end parallel


      !$acc parallel async
      !$acc loop gang vector collapse(2)
      DO kso = 1, ke_soil
        DO i = ivstart, ivend
          zalam(i,kso) = zalam(i,kso) + zalamtmp(i,kso)
          hzalam(i,kso) = zalam(i,kso)
        ENDDO
      ENDDO
      !$acc end parallel





     ! Heat capacity
     ! -------------------

     !$acc parallel async
     !$acc loop gang vector
     DO i = ivstart, ivend
       zw_fr  (i,ke_soil+1)  = w_so_now(i,ke_soil+1)/zdzhs(ke_soil+1)
     ENDDO
     !$acc end parallel

     ! REORDER
     !$acc parallel async
     !$acc loop gang vector collapse(2)
     DO kso   = 1, ke_soil
       DO i = ivstart, ivend
         zw_fr   (i,kso)     = w_so_now(i,kso)/zdzhs(kso)
       ENDDO
     ENDDO
     !$acc end parallel



     !$acc parallel async
     !$acc loop gang vector collapse(2) private(zzz)
     DO kso = 1, ke_soil
       DO i = ivstart, ivend

         ! Scale soil heat capacity with organic fraction -> Chadburn et al., 2015
         IF (zmls(kso) < rootdp(i)) THEN
           zzz = plcov(i)*(rootdp(i)-zmls(kso))/rootdp(i)
           zrocg(i,kso)=(1.0_vpp-zzz)*zrocg_soil(i,kso)+zzz*0.58E+06_vpp
         END IF

       ENDDO
     ENDDO
     !$acc end parallel

     
     !$acc parallel async
     !$acc loop seq
     DO   kso = 1,ke_soil+1
       !$acc loop gang vector
       DO i = ivstart, ivend
         ziw_fr(i,kso) = w_so_ice_now(i,kso)/zdzhs(kso)                            ! ice frac.

         zlw_fr(i,kso) = zw_fr(i,kso) - ziw_fr(i,kso)                              ! liquid water frac.

         zroc(i,kso)   = zrocg(i,kso) + rho_w*zlw_fr(i,kso)*chc_w +          &     ! soil  heat capacity
                                        rho_w*ziw_fr(i,kso)*chc_i
        END DO
      END DO      !soil layers
      !$acc end parallel


     ! --------------------
     ! Call the solver
     ! ---------------------

     !$acc parallel async
     !$acc loop gang vector private(zm,hcon,hcap,hdif,rho,t_sol,sw_abs,alpha,gamma_sol,a,b,c,d,e,l_top,counter)
     DO i = ivstart, ivend

       l_top = top(i)

       IF(l_top .GT. 1) THEN  !!!snow on the ground


         !CALL solve_1d_heat(dzm_sn_now(i,:), t_sn_now(i,:), t_sn_sfc(i), for_sn(i), swabs_sn(i,:)  , &
         !                   zm_sn(i,:), hcap_sn(i,:), hcon_sn(i,:), hdif_sn(i,:), rho_sn(i,:)      , &
         !                   top(i), t_so_now(i,:), zmls(:), zroc(i,:), zalam(i,:)                  , & 
         !                   ke_soil, zdt   )   

         !! inlining the heat solver

        zm    = 0.0_vpp
        hcon  = 0.0_vpp
        hcap  = 0.0_vpp
        rho   = 0.0_vpp
        t_sol = 0.0_vpp
        sw_abs  = 0.0_vpp
        counter = 1
  
        do ksn = -l_top+1, ke_soil+1, 1

               IF(ksn .LE. 0) THEN  ! snow layers
                   IF(ksn .EQ. -l_top+1) THEN
                        zm(ksn)     = zm_sn(i,l_top)
                        hcon(ksn)   = hcon_sn(i,l_top)
                        hcap(ksn)   = hcap_sn(i,l_top)
                        t_sol(ksn)  = t_sn_now(i,l_top)
                        sw_abs(ksn) = swabs_sn(i,l_top)
                   ELSE
                        zm(ksn)     = zm_sn(i,l_top-counter)
                        hcon(ksn)   = hcon_sn(i,l_top-counter)
                        hcap(ksn)   = hcap_sn(i,l_top-counter)
                        t_sol(ksn)  = t_sn_now(i,l_top-counter)
                        sw_abs(ksn) = swabs_sn(i,l_top-counter)
                        counter = counter + 1
                   ENDIF

               ELSE  ! soil layers
                        zm(ksn)     = zm_sn(i,1) + zmls(ksn)
                        !hcon(ksn)   = zalam(i,ksn)
                        hcap(ksn)   = zroc(i,ksn)
                        t_sol(ksn)  = t_so_now(i,ksn)
                        sw_abs(ksn) = 0.0_vpp
               ENDIF

        ENDDO ! end of snow layers

        hcon(1:ke_soil) = zalam(i,1:ke_soil)
  ! ------------------------------------------------------------
  ! Some precalculations ...
  ! ------------------------------------------------------------

  ! Derivative of emitted long wave radiation
  ! ----------------------

    dlw_u_sn =  -4.0_vpp * sigma * (1.0_vpp - snow_Ctalb) * t_sol(-l_top+1)**3


  ! Calculate factors (diffusion) for the linear equations for ...
  ! ----------------------

    ! Loop over all layers ...
    DO ksn = -l_top+1, ke_soil, 1
      IF(ksn .LE. ke_soil-1) THEN ! ... except
        alpha(ksn)   = dt / hcap(ksn)
        hdif(ksn) = hcon(ksn) * (t_sol(ksn+1) - t_sol(ksn)) /  (zm(ksn+1) - zm(ksn))
      ELSE !(ksn == substrate) ! ... bottom layer
         alpha(ksn) = dt / hcap(ksn)
         hdif(ksn) = 0.0_vpp
      ENDIF
    END DO


 ! ------------------------------------------------------------
  ! Setup tridiagonal matrix for set of linear equations for each layer ...
  ! ------------------------------------------------------------

  DO ksn = -l_top+1, ke_soil

    IF(ksn .EQ. -l_top+1) THEN ! ... TOP LAYER

    dz_low = zm(ksn+1) - zm(ksn)

      a(ksn)     = 0.0_vpp
      b(ksn)     = 1 + (1.0_vpp - cn) * alpha(ksn) * hcon(ksn)/dz_low - alpha(ksn) * dlw_u_sn
      c(ksn)     = -   (1.0_vpp - cn) * alpha(ksn) * hcon(ksn)/dz_low

      d(ksn)     = t_sol(ksn) + alpha(ksn) * (for_sn(i) - dlw_u_sn*t_sol(ksn) + cn*hdif(ksn))

    ELSEIF (ksn .LE. ke_soil-1) THEN ! ... INNER LAYERS

     dz_up  = zm(ksn)   - zm(ksn-1)
     dz_low = zm(ksn+1) - zm(ksn)

        a(ksn) = -         (1.0_vpp - cn)   * alpha(ksn) *  hcon(ksn-1)/dz_up
        b(ksn) = 1.0_vpp + (1.0_vpp - cn)   * alpha(ksn) * (hcon(ksn)  /dz_low + hcon(ksn-1)/dz_up)
        c(ksn) = -         (1.0_vpp - cn)   * alpha(ksn) *  hcon(ksn)  /dz_low

        d(ksn) = t_sol(ksn) + cn*alpha(ksn) * (hdif(ksn) - hdif(ksn-1)) + alpha(ksn)*sw_abs(ksn)


    ELSEIF (ksn .EQ. ke_soil) THEN ! BOTTOM LAYERS

     dz_up = zm(ksn)   - zm(ksn-1)

       a(ksn) = -         (1.0_vpp - cn) * alpha(ksn) * hcon(ksn-1)/dz_up
       b(ksn) = 1.0_vpp + (1.0_vpp - cn) * alpha(ksn) * hcon(ksn-1)/dz_up
       c(ksn) = 0.0_vpp

       !d(bot_idx) = t(bot_idx) - cn*alpha(bot_idx-1) + alpha(bot_idx)*hdif(bot_idx)
       d(ksn)     = t_sol(ksn) - cn*alpha(ksn)*hdif(ksn-1) + alpha(ksn)*hdif(ksn)

    ENDIF

  ENDDO

  ! ------------------------------------------------------------
  ! Solve the system - Thomas Algorithm
  ! ------------------------------------------------------------

   beta = b(-l_top+1)
    ! Forward substitution

    DO ksn = -l_top+1, ke_soil, 1
      IF(ksn .GE. -l_top+1) THEN
        IF(ksn .EQ. -l_top+1) THEN
          e(ksn) = d(ksn) / beta
        ELSE
          gamma_sol(ksn) = c(ksn-1) / beta
          beta       = b(ksn) - a(ksn) * gamma_sol(ksn)
          e(ksn)     = (d(ksn) - a(ksn) * e(ksn-1)) / beta
        ENDIF
      ENDIF
    ENDDO

    ! Backward substitution

    DO ksn = ke_soil-1, -l_top+1, -1
      IF(ksn .GE. -l_top+1) THEN
        e(ksn) = e(ksn) - gamma_sol(ksn+1) * e(ksn+1)
      ENDIF
    ENDDO

   ! ------------------------------------------------------------
   ! Do some updating required for the next sections
   ! ------------------------------------------------------------
   ! Snow Surface Temperature

    t_sn_sfc(i) = e(-l_top+1)
    ! Snow layer temperature
    counter = 1
    DO ksn = l_top,1,-1
      IF(ksn .EQ. l_top) THEN
       t_sn_now(i,ksn)    = e(-l_top+1)
      ELSE
       t_sn_now(i,ksn)    = e(-l_top+1+counter)
       counter = counter + 1
      ENDIF
    ENDDO

       ENDIF

     ENDDO
     !$acc end parallel

! =============================================================================
! + End Section II: Atmospheric forcing & Heat Equation
! =============================================================================

!if(ntstep .eq. 165052) then
!   
!   write(*,*) 'end section II: ',h_snow(1),hm_sn(1,top(1)) !h_snow(1)
!
!endif
if(ntstep .eq. N_IMP) then
   do ksn = 1,top(1)
     write(*,*) 'end section II: ',ksn,theta_i_now(1,ksn),theta_a_now(1,ksn),theta_w_now(1,ksn),t_sn_now(1,ksn)

   enddo
endif



! =============================================================================
! + Begin Section III: Phase Changes
! =============================================================================



    ! Keep old values
    ! -----------------

    !$acc parallel async default(none)
    !$acc loop seq
    DO ksn = 1, n_layers, 1
      !$acc loop gang vector
       DO i = ivstart, ivend

         zm_sn_old(i,ksn)   = zm_sn(i,ksn)
    
         theta_i_old(i,ksn) = theta_i_now(i,ksn)
         theta_w_old(i,ksn) = theta_w_now(i,ksn)

       ENDDO
    ENDDO
    !$acc end parallel
    
    ! Initiate phase change
    ! -------------------

     !$acc parallel async
     !$acc loop gang vector private(l_top)
     DO i = ivstart, ivend
       l_top = top(i)
       IF(l_top .GE. 1) THEN ! snow on the ground

         !CALL init_phase_change(t_sn_now(i,:)                       , & 
         !                       theta_i_now(i,:), theta_w_now(i,:)  , &
         !                       top(i)                              , &
         !                       melt_flag(i), freeze_flag(i)) 

         ! Initiate
         melt_flag(i)   = .FALSE.
         freeze_flag(i) = .FALSE.

         DO ksn = l_top, 1, -1
            ! Metling
            IF(t_sn_now(i,ksn) .GT. t0_melt .AND. theta_i_now(i,ksn) .GT. eps2) THEN
                 melt_flag(i) = .TRUE.
            ENDIF
            ! Freezing
            IF(t_sn_now(i,ksn) .LT. t0_melt .AND. theta_w_now(i,ksn) .GT. (theta_r + eps2)) THEN
                 freeze_flag(i) = .TRUE.
            ENDIF
         ENDDO  !end of snow layers

       ENDIF

     ENDDO
     !$acc end parallel 


    ! Check for melting
    ! --------------------
    
     !$acc parallel async
     !$acc loop gang vector private(l_top)
     DO i = ivstart, ivend
       l_top = top(i)
       IF(l_top .GE. 1) THEN ! snow on the ground
         IF(melt_flag(i)) THEN ! melting

           !CALL melt_snow(zm_sn(i,:), dzm_sn_now(i,:), t_sn_now(i,:), hcap_sn(i,:),  &
           !               theta_i_now(i,:), theta_w_now(i,:), rho_sn(i,:),           &
           !               top(i))

           q_rest = 0.0_vpp
           DO ksn = l_top, 1, -1

                ! --------------------
                ! Now see if any melting is going on -- this implies that (1) the
                ! temperature of the element
                ! is above/equal the melting temperature (2) there is something to melt and
                ! (3) there is enough
                ! room to place the meltwater ...
                ! --------------------
                
                IF(t_sn_now(i,ksn) .GE. t0_melt .AND. theta_i_now(i,ksn) .GT. 0.0_vpp .AND. theta_w_now(i,ksn) .LT. theta_s) THEN
             
                 ! difference dT between actual melting temperature and layer temperature
                 dT_sub_melt = 0.0_vpp
                 dT_sub_melt = t0_melt - t_sn_now(i,ksn)
             
                ! --------------------
                ! Now we take into account that there might be some extra energy that could
                ! not
                ! be used by the element above because of complete melting
                ! --------------------
             
                  dT_sub_melt = dT_sub_melt - (q_rest / (hcap_sn(i,ksn) * rho_sn(i,ksn) * dzm_sn_now(i,ksn)))
                  ! ONly do it when there is real potential to melt
                  IF(dT_sub_melt .LT. 0.0_vpp) THEN
                    ! --------------------
                    ! Determine the DECREASE in ice content and the INCREASE of water
                    ! content
                    ! Adapt A_sub_melt to compute mass changes
                    ! --------------------
                     A_sub_melt = (hcap_sn(i,ksn) * rho_sn(i,ksn)) / (rho_i * lh_f)
                     dtheta_i = A_sub_melt * dT_sub_melt
                     dtheta_w = - (rho_i / rho_w) * dtheta_i
                    ! --------------------
                    ! It could happen that there is enough energy available to melt more ice
                    ! than is present.
                    ! You can only melt so much ice as is there ....
                    ! --------------------
                     IF( (theta_i_now(i,ksn) + dtheta_i) .LT. 0.0_vpp) THEN
                       dtheta_i = - theta_i_now(i,ksn)
                       dtheta_w = - (rho_i / rho_w) * dtheta_i
                       dT_sub_melt = dtheta_i / A_sub_melt
                     ENDIF
  
                    ! --------------------
                    ! It could also be that you are trying to produce more water than is
                    ! allowed.
                    ! --------------------
                      IF( (theta_w_now(i,ksn) + dtheta_w) .GT. theta_s) THEN
                        dtheta_w = theta_s - theta_w_now(i,ksn)
                        dtheta_i = - (rho_w/rho_i) *dtheta_w
                        dT_sub_melt = dtheta_i / A_sub_melt
                      ENDIF
  
                    ! --------------------
                    ! Reset/Recalculate properties
                    ! --------------------
                     ! Layer temperature and transfered energy
                     t_sn_now(i,ksn) = t_sn_now(i,ksn) + dT_sub_melt
                     IF(t_sn_now(i,ksn) .LE. t0_melt) THEN ! if melting occured it can only be at melt point
                       q_rest     = 0.0_vpp
                       t_sn_now(i,ksn)  = t0_melt
                     ELSE
                       q_rest = hcap_sn(i,ksn) * rho_sn(i,ksn) * dzm_sn_now(i,ksn) * (t_sn_now(i,ksn) - t0_melt)
                       t_sn_now(i,ksn)   = t0_melt
                     ENDIF
                     ! Volumetric freezing power
                     q_mf     = q_mf + ((dtheta_i * rho_i * lh_f) / dt)
                     dtheta_w = dtheta_w
                     ! Contents of ice, water and air
                     theta_i_now(i,ksn) = theta_i_now(i,ksn) + dtheta_i
                     theta_w_now(i,ksn) = theta_w_now(i,ksn) + dtheta_w
                  ENDIF ! deltaT check
                 ENDIF ! end melt check
           ENDDO ! loop over snow layers

           ! Update dependent variables
           !CALL update(top(i), hm_sn(i,:), zm_sn(i,:), dzm_sn_now(i,:), m_sn(i,:)          , &
           !            theta_i_now(i,:), theta_w_now(i,:), theta_a_now(i,:), rho_sn(i,:)   , &
           !            hcap_sn(i,:), hcon_sn(i,:))

           ! -------------------------
           ! Height of snow (main) levels
           ! -------------------------
                 !$acc loop seq
                 do ksn=1,l_top
                    if(ksn .eq. 1) then
                       hm_sn(i,ksn) = dzm_sn_now(i,ksn)
                    else
                       hm_sn(i,ksn) = hm_sn(i,ksn-1) + dzm_sn_now(i,ksn)
                    endif
                 enddo  
         
           ! --------------------------
           ! Depth of snow layer (main) levels
           ! -------------------------
                 !$acc loop seq
                 do ksn = l_top,1,-1
                    zm_sn(i, (l_top+1) - ksn ) = hm_sn(i,ksn) !invert height vector
                    theta_a_now(i,ksn) = max(0.0_vpp,1.0_vpp - theta_i_now(i,ksn) - theta_w_now(i,ksn))
                 enddo
         
           ! --------------------------
           ! Snow layer density
           ! -------------------------
                  !$acc loop seq
                  do ksn = 1, l_top, 1
                     if(theta_i_now(i,ksn) .eq. 0.0_vpp) then
                        rho_sn(i,ksn) = 0.0_vpp
                     else
                        rho_sn(i,ksn) = theta_i_now(i,ksn)*rho_i + theta_w_now(i,ksn)*rho_w
                     endif
                  enddo
                  
           ! --------------------------
           ! Heat capacity
           ! --------------------------
            !$acc loop seq
            DO ksn = 1, l_top, 1
                IF(rho_sn(i,ksn) .LT. eps) THEN
                  hcap_sn(i,ksn) = 0.0_vpp
                ELSE
                  hcap_sn(i,ksn) = (  rho_a   * theta_a_now(i,ksn) * specific_heat_air     &
                                    + rho_i   * theta_i_now(i,ksn) * specific_heat_ice     &
                                    + rho_w   * theta_w_now(i,ksn) * specific_heat_water)  &
                                     / rho_sn(i,ksn)
                ENDIF
            ENDDO
         
           ! --------------------------
           ! Heat conductivity
           ! --------------------------
            !$acc loop seq
            DO ksn = 1, l_top, 1
               hcon_sn(i,ksn) = 2.22_vpp * EXP(1.88_vpp * LOG(rho_sn(i,ksn)/rho_i))
            ENDDO
         
           ! --------------------------
           ! Snow layer mass
           ! --------------------------
            !$acc loop seq
            DO ksn = 1, l_top, 1
              IF(ksn .EQ. 1) THEN
                 m_sn(i,ksn) = zm_sn(i,ksn) * ((theta_i_now(i,ksn) * rho_i) + (theta_w_now(i,ksn) * rho_w))
               ELSE
                 m_sn(i,ksn) = ABS(zm_sn(i,ksn) - zm_sn(i,ksn-1))  * ((theta_i_now(i,ksn) * rho_i) + (theta_w_now(i,ksn) * rho_w))
               ENDIF
            ENDDO
         
          ENDIF ! end check in melt or not
       ENDIF ! end check if top .ge. 1
     ENDDO ! end loop over horizontal pixels
     !$acc end parallel


     ! Check for freezing
     ! --------------------

     !$acc parallel async
     !$acc loop gang vector
     DO i = ivstart, ivend
       l_top = top(i)
       IF(l_top .GE. 1) THEN ! snow on the ground
         IF(freeze_flag(i)) THEN ! freeze

           !CALL freeze_snow(zm_sn(i,:), dzm_sn_now(i,:), t_sn_now(i,:), hcap_sn(i,:)  , &
           !                 theta_i_now(i,:), theta_w_now(i,:), theta_a_now(i,:)        , &
           !                 rho_sn(i,:), top(i))
           DO  ksn = l_top, 1, -1
              ! Freezing within the snowpack can occur if (1) the temperature of the
              ! element is below freezing and if water is present to be refrozen
              IF(t_sn_now(i,ksn) .LT. t0_melt .AND. theta_w_now(i,ksn) .GT. theta_r) THEN
                ! difference dT between actual layer temperature and freezing temperature
                dT_sub_melt = 0.0_vpp
                dT_sub_melt = t0_melt - t_sn_now(i,ksn)
                ! Adapt A_sub_melt to compute mass change
                A_sub_melt = (hcap_sn(i,ksn) * rho_sn(i,ksn)) / (rho_i * lh_f)
                ! Compute change in volumetric contenst
                dtheta_i = A_sub_melt * dT_sub_melt
                dtheta_w = -(rho_i/rho_w) * dtheta_i
                ! Make sure that there is enough water to refreeze
                IF( (theta_w_now(i,ksn) + dtheta_w) .LT. theta_r) THEN
                  dtheta_w = -ABS(theta_w_now(i,ksn) - theta_r)
                  dtheta_i = -(rho_w / rho_i) * dtheta_w
                  dT_sub_melt       = dtheta_i / A_sub_melt
                ENDIF
                ! See if the layer is pure ice
                IF( (theta_i_now(i,ksn) + theta_r + dtheta_i) .GE. 1.0_vpp) THEN
                  dtheta_w = - ABS(theta_w_now(i,ksn) - theta_r)
                  dtheta_i = - (rho_w/rho_i) * dtheta_w
                  theta_i_now(i,ksn)  = 1.0_vpp
                  theta_w_now(i,ksn)  = theta_r
                  theta_a_now(i,ksn)  = 0.0_vpp
                ELSE
                  theta_i_now(i,ksn) = theta_i_now(i,ksn) + dtheta_i
                  theta_w_now(i,ksn) = theta_w_now(i,ksn) + dtheta_w
                  theta_a_now(i,ksn) = MAX(0.0_vpp, 1.0_vpp - theta_i_now(i,ksn) - theta_w_now(i,ksn))
                ENDIF
         
                ! --------------------
                ! Reset/Recalculate properties
                ! --------------------
                ! Set some limits
                IF(theta_w_now(i,ksn) .GE. 1.0_vpp) THEN
                  theta_w_now(i,ksn) = 1.0_vpp
                ENDIF
         
                !Compute the volumetric refreezing power
                q_mf     = q_mf + ((dtheta_i * rho_i * lh_f) / dt)
                dtheta_w = dtheta_w
                t_sn_now(i,ksn)     = t_sn_now(i,ksn) + dT_sub_melt
             ENDIF
           ENDDO
 
           ! Update dependent variables
!           CALL update(top(i), hm_sn(i,:), zm_sn(i,:), dzm_sn_now(i,:), m_sn(i,:)          , &
!                       theta_i_now(i,:), theta_w_now(i,:), theta_a_now(i,:), rho_sn(i,:)   , &
!                       hcap_sn(i,:), hcon_sn(i,:))

           ! -------------------------
           ! Height of snow (main) levels
           ! -------------------------
                 !$acc loop seq
                 do ksn=1,l_top
                    if(ksn .eq. 1) then
                       hm_sn(i,ksn) = dzm_sn_now(i,ksn)
                    else
                       hm_sn(i,ksn) = hm_sn(i,ksn-1) + dzm_sn_now(i,ksn)
                    endif
                 enddo  
         
           ! --------------------------
           ! Depth of snow layer (main) levels
           ! -------------------------
                 !$acc loop seq
                 do ksn = l_top,1,-1
                    zm_sn(i, (l_top+1) - ksn ) = hm_sn(i,ksn) !invert height vector
                    theta_a_now(i,ksn) = max(0.0_vpp,1.0_vpp - theta_i_now(i,ksn) - theta_w_now(i,ksn))
                 enddo
         
           ! --------------------------
           ! Snow layer density
           ! -------------------------
                  !$acc loop seq
                  do ksn = 1, l_top, 1
                     if(theta_i_now(i,ksn) .eq. 0.0_vpp) then
                        rho_sn(i,ksn) = 0.0_vpp
                     else
                        rho_sn(i,ksn) = theta_i_now(i,ksn)*rho_i + theta_w_now(i,ksn)*rho_w
                     endif
                  enddo
                  
           ! --------------------------
           ! Heat capacity
           ! --------------------------
            !$acc loop seq
            DO ksn = 1, l_top, 1
                IF(rho_sn(i,ksn) .LT. eps) THEN
                  hcap_sn(i,ksn) = 0.0_vpp
                ELSE
                  hcap_sn(i,ksn) = (  rho_a   * theta_a_now(i,ksn) * specific_heat_air     &
                                    + rho_i   * theta_i_now(i,ksn) * specific_heat_ice     &
                                    + rho_w   * theta_w_now(i,ksn) * specific_heat_water)  &
                                     / rho_sn(i,ksn)
                ENDIF
            ENDDO
         
           ! --------------------------
           ! Heat conductivity
           ! --------------------------
            !$acc loop seq
            DO ksn = 1, l_top, 1
               hcon_sn(i,ksn) = 2.22_vpp * EXP(1.88_vpp * LOG(rho_sn(i,ksn)/rho_i))
            ENDDO
         
           ! --------------------------
           ! Snow layer mass
           ! --------------------------
            !$acc loop seq
            DO ksn = 1, l_top, 1
              IF(ksn .EQ. 1) THEN
                 m_sn(i,ksn) = zm_sn(i,ksn) * ((theta_i_now(i,ksn) * rho_i) + (theta_w_now(i,ksn) * rho_w))
               ELSE
                 m_sn(i,ksn) = ABS(zm_sn(i,ksn) - zm_sn(i,ksn-1))  * ((theta_i_now(i,ksn) * rho_i) + (theta_w_now(i,ksn) * rho_w))
               ENDIF
            ENDDO

          ENDIF ! if freeze
       ENDIF ! if there is snow or not
     ENDDO ! end loop over horizontal pixels
     !$acc end parallel

! =============================================================================
! - End Section III: Phase Changes
! =============================================================================

if(ntstep .eq. N_IMP) then
   do ksn = 1,top(1)
     write(*,*) 'end section III: ',ksn,theta_i_now(1,ksn),theta_a_now(1,ksn),theta_w_now(1,ksn),t_sn_now(1,ksn)

   enddo
endif


!REAL (KIND = vpp)        :: &
!
!    dL                         , & ! change of layer thickness
!    dM                         , & ! change of layer mass
!    M                          , & ! initial mass and volmetric content (water of ice)
!    hoar                       , & ! hoar mass
!    dzm_sn_old                     ! old value of layer thickness

! =============================================================================
! + Begin Section IV: Water Transport
! =============================================================================

  !$acc parallel async
  !$acc loop gang vector private(l_top)
  DO i = ivstart, ivend
    l_top = top(i)
    IF(l_top .GE. 1) THEN ! snow on the ground

      ! Compute Sublimation/Deposition/Condensation/Evaporation
      ! -------------------
      !CALL phase_change_surface(top(i), lh_sn(i), t_sn_sfc(i), rho_sn(i,:), zdt     , &
      !                          dzm_sn_now(i,:), theta_i_now(i,:), theta_w_now(i,:)   )

      ! Initiate some values
       dL   = 0.0_vpp
       dM   = 0.0_vpp
       M    = 0.0_vpp
       hoar = 0.0_vpp
      ! Latent heat flux towards the surface - mass gain
      ! ---------------------
      IF(lh_sn(i) .GT. eps2) THEN ! add mass
        IF(t_sn_sfc(i) .LT. t0_melt) THEN  ! add ice
    
          dM = lh_sn(i) * dt/lh_s     ! calculate mass change
          lh_sn(i) = 0.0_vpp           ! reset latent heat flux, i.e. energy was used
          hoar = dM
    
          ! Adjust layer properties accordingly, keep snow density constant
    
          dzm_sn_old = dzm_sn_now(i,l_top)
          dL = dM/rho_sn(i,l_top)
          dzm_sn_now(i,l_top) = dzm_sn_now(i,l_top) + dL
          theta_i_now(i,l_top) = theta_i_now(i,l_top) * (dzm_sn_old/dzm_sn_now(i,l_top))
          theta_i_now(i,l_top) = theta_i_now(i,l_top) + (dM/(rho_i*dzm_sn_now(i,l_top)))
          theta_w_now(i,l_top) = theta_w_now(i,l_top) * (dzm_sn_old/dzm_sn_now(i,l_top))
                  if (ntstep .eq. N_IMP) then
                      write(*,*) 'update_C: ',i,ksn,theta_i_now(i,ksn)
                  endif

        ELSE ! add water
          dM = lh_sn(i) *dt/lh_v      ! calculate mass change
          lh_sn(i) = 0.0_vpp           ! reset latent heat, i.e. energy was used
          theta_w_now(i,l_top) = theta_w_now(i,l_top) + dM/(rho_w*dzm_sn_now(i,l_top))   ! update volumetric water content
        ENDIF
      ! Latent heat flux away from the surface - mass loss
      ! ---------------------
      ELSE
       IF(lh_sn(i) .LT. (-1.0_vpp*eps2)) THEN ! additional check in case lh ist super small, but sligtly positive
         ! Not sure how acc statements need to look for this construct
         DO WHILE (lh_sn(i) .LT. (-1.0_vpp*eps2)) ! while energy is available
           DO ksn = l_top, 1, -1 ! loop through snow layers
             IF(theta_w_now(i,ksn) .GT. eps) THEN ! there is water, i.e. evaporate first
               ! Calculate mass change
               dM = lh_sn(i) * dt/lh_v
                M = theta_w_now(i,ksn) * rho_w * dzm_sn_now(i,ksn)
                 ! Check that you only take the available amount of water
                 IF(-dM .GE. M) THEN
                   dM = -M
                   theta_w_now(i,ksn) = theta_w_now(i,ksn) + dM/(rho_w*dzm_sn_now(i,ksn))
                 ELSE
                   theta_w_now(i,ksn) = theta_w_now(i,ksn) + dM/(rho_w*dzm_sn_now(i,ksn))
                 ENDIF
               lh_sn(i) = lh_sn(i) - dM*lh_v/dt ! update energy used
             ELSEIF (theta_i_now(i,ksn) .GT. eps) THEN ! there is no water then sublimate ice matrix
               dM = lh_sn(i) * dt/lh_s
                M = theta_i_now(i,ksn) * rho_i * dzm_sn_now(i,ksn)
                IF(-dM .GT. M) THEN ! all ice can be sublimated
                  dM = -M
                  theta_i_now(i,ksn) = 0.0_vpp
                  dzm_sn_now(i,ksn)  = 0.0_vpp
                  if (ntstep .eq. N_IMP) then
                      write(*,*) 'update_B: ',ntstep,i,ksn,theta_i_now(i,ksn)
                  endif

                ELSE
                  dzm_sn_old = dzm_sn_now(i,ksn)
                  dL = dM/rho_sn(i,ksn)
                  dzm_sn_now(i,ksn) = dzm_sn_now(i,ksn) + dL
                  theta_i_now(i,ksn) = theta_i_now(i,ksn) * (dzm_sn_old/dzm_sn_now(i,ksn))
                  theta_i_now(i,ksn) = theta_i_now(i,ksn) + (dM/(rho_i*dzm_sn_now(i,ksn)))
                  theta_w_now(i,ksn) = theta_w_now(i,ksn) * (dzm_sn_old/dzm_sn_now(i,ksn))
                  if (ntstep .eq. N_IMP) then
                      write(*,*) 'update_A: ',i,ksn,theta_i_now(i,ksn)
                  endif
                ENDIF
               lh_sn(i) = lh_sn(i) - dM*lh_v/dt ! update energy used
             ENDIF
           ENDDO ! end of ksn
   
       ! MASSIVE HACK HERE
       IF(lh_sn(i) .LT. (-1.0_vpp*eps2)) THEN ! there is still energy left, which should technically be used by the soil layer for now let's erase it
         lh_sn(i) = 0.0_vpp
       ENDIF
     ENDDO ! end of while
   ENDIF
  ENDIF

if(ntstep .eq. N_IMP) then
   do ksn = 1,top(1)
     write(*,*) 'in between section IV (A): ',ksn,theta_i_now(1,ksn),theta_a_now(1,ksn),theta_w_now(1,ksn)
   enddo
endif

      ! Update dependent variables
      !CALL update(top(i), hm_sn(i,:), zm_sn(i,:), dzm_sn_now(i,:), m_sn(i,:)         , &
      !            theta_i_now(i,:), theta_w_now(i,:), theta_a_now(i,:),rho_sn(i,:)   , &
      !            hcap_sn(i,:), hcon_sn(i,:))

           ! -------------------------
           ! Height of snow (main) levels
           ! -------------------------
                 !$acc loop seq
                 do ksn=1,l_top
                    if(ksn .eq. 1) then
                       hm_sn(i,ksn) = dzm_sn_now(i,ksn)
                    else
                       hm_sn(i,ksn) = hm_sn(i,ksn-1) + dzm_sn_now(i,ksn)
                    endif
                 enddo  
         
           ! --------------------------
           ! Depth of snow layer (main) levels
           ! -------------------------
                 !$acc loop seq
                 do ksn = l_top,1,-1
                    zm_sn(i, (l_top+1) - ksn ) = hm_sn(i,ksn) !invert height vector
                    theta_a_now(i,ksn) = max(0.0_vpp,1.0_vpp - theta_i_now(i,ksn) - theta_w_now(i,ksn))
                 enddo
         
           ! --------------------------
           ! Snow layer density
           ! -------------------------
                  !$acc loop seq
                  do ksn = 1, l_top, 1
                     if(theta_i_now(i,ksn) .eq. 0.0_vpp) then
                        rho_sn(i,ksn) = 0.0_vpp
                     else
                        rho_sn(i,ksn) = theta_i_now(i,ksn)*rho_i + theta_w_now(i,ksn)*rho_w
                     endif
                  enddo
                  
           ! --------------------------
           ! Heat capacity
           ! --------------------------
            !$acc loop seq
            DO ksn = 1, l_top, 1
                IF(rho_sn(i,ksn) .LT. eps) THEN
                  hcap_sn(i,ksn) = 0.0_vpp
                ELSE
                  hcap_sn(i,ksn) = (  rho_a   * theta_a_now(i,ksn) * specific_heat_air     &
                                    + rho_i   * theta_i_now(i,ksn) * specific_heat_ice     &
                                    + rho_w   * theta_w_now(i,ksn) * specific_heat_water)  &
                                     / rho_sn(i,ksn)
                ENDIF
            ENDDO
         
           ! --------------------------
           ! Heat conductivity
           ! --------------------------
            !$acc loop seq
            DO ksn = 1, l_top, 1
               hcon_sn(i,ksn) = 2.22_vpp * EXP(1.88_vpp * LOG(rho_sn(i,ksn)/rho_i))
            ENDDO
         
           ! --------------------------
           ! Snow layer mass
           ! --------------------------
            !$acc loop seq
            DO ksn = 1, l_top, 1
              IF(ksn .EQ. 1) THEN
                 m_sn(i,ksn) = zm_sn(i,ksn) * ((theta_i_now(i,ksn) * rho_i) + (theta_w_now(i,ksn) * rho_w))
               ELSE
                 m_sn(i,ksn) = ABS(zm_sn(i,ksn) - zm_sn(i,ksn-1))  * ((theta_i_now(i,ksn) * rho_i) + (theta_w_now(i,ksn) * rho_w))
               ENDIF
            ENDDO


      ! Rain on snow
      ! -------------------

        ! Note: This still needs to be integrated. Simply use zrain_rate and put
        ! it in the top layer. We also would need to put the energy as an
        ! additional source into the solver for the heat equation or in the
        ! forcing directly.

      ! -------------------
      ! Transport the water
      ! -------------------

      ! CALL transport_water(zm_sn(i,:), dzm_sn_now(i,:), t_sn_now(i,:), rho_sn(i,:),   &
      !                      theta_i_now(i,:), theta_a_now(i,:), theta_w_now(i,:),      &
      !                      hcap_sn(i,:), top(i), runoff_sn(i))
      frac_rho      = rho_w/rho_i
      limit_theta_i =   1.0_vpp - frac_rho * ((1.0_vpp + 0.0165_vpp * frac_rho)                   &
                      - SQRT((1.0_vpp + 0.0165_vpp * frac_rho)*(1.0_vpp + 0.0165_vpp * frac_rho)  &
                      - 4.0_vpp * frac_rho * 0.0264_vpp)) / (2.0_vpp * frac_rho)
   
      ! + Calculate a few properties
      !--------------------
      ! Loop over number of snow layers
      DO ksn = 1, l_top, 1
        ! Determine the additional storage capacity (of water) due to refreezing
        ! --------------------------
         dtheta_w_sub_wtr(ksn) = hcap_sn(i,ksn) * rho_sn(i,ksn) / lh_f / rho_w * MAX(0.0_vpp, (t0_melt - t_sn_now(i,ksn)))
        ! --------------------------
        ! Estimate resiudal water (RWC) content by volume; Coleou and Lesaffre, 1998,
        ! Ann. Glaciol., 26, 64-68
        ! --------------------------
        IF(theta_i_now(i,ksn) .GT. limit_theta_i) THEN
          ! This case is the limiting case where:
          ! theta_i + (theta_r * (rho_w/rho_i)) >= 1.0
          ! In that case, set the residual water content equal to the pore space
          res_wat_cont(ksn) = (1.0_vpp - theta_i_now(i,ksn)) * (rho_i/rho_w)
        ELSE
          IF(theta_i_now(i,ksn) .GT. 0.23_vpp) THEN
            res_wat_cont(ksn) = 0.0264_vpp + 0.0099_vpp * (1.0_vpp - theta_i_now(i,ksn)) / theta_i_now(i,ksn)
          ELSE
            res_wat_cont(ksn) = 0.08_vpp - 0.1023_vpp * (theta_i_now(i,ksn) - 0.03_vpp)
          ENDIF
        ENDIF
        ! Limit residual water content
        res_wat_cont(ksn) = MIN(res_wat_cont(ksn), 0.08_vpp)  ! NOTE: only needed in case of theta_i < 0.03
        ! Effective residual water content
        w_res(ksn) = MIN(1.0_vpp - theta_i_now(i,ksn) * rho_i/rho_w, res_wat_cont(ksn) + dtheta_w_sub_wtr(ksn))
        w_res(ksn) = MAX(0.0_vpp, w_res(ksn))
     END DO
   
      ! Now start moving the water and adjust properties accordingly
      ! ----------------------------------
      DO ksn = l_top, 2, -1
        ! Reset excess water to zero
        IF(theta_i_now(i,ksn) .LT. eps) THEN ! no more ice in this layer only residual water which needs to be moved
          excess_water = theta_w_now(i,ksn) ! add residual water to excess water and ...
          theta_w_now(i,ksn) = 0.0_vpp       ! reset water content of said layer
        ELSE
          excess_water = 0.0_vpp
        ENDIF
      ! water content of the upper layer
       w_up = theta_w_now(i,ksn)
        IF(ksn .EQ. l_top .AND. w_up .GT. 0.0_vpp .AND. w_up .LE. w_res(ksn)) THEN
          ! In that case you need to update the volumetric air content and the
          ! density of the top element
          ! as it may have caught some rain! Only top element should be considered,
          ! as when rain would have
          ! infiltrated lower elements as well, w_up > w_res.
          theta_a_now(i,ksn) = MAX(0.0_vpp, 1.0_vpp - theta_w_now(i,ksn) - theta_i_now(i,ksn))
          rho_sn(i,ksn)  = (theta_i_now(i,ksn) * rho_i) + (theta_w_now(i,ksn) * rho_w)
            ! we want positive densities
            IF(rho_sn(i,ksn) .LT. 0.0_vpp) THEN
              rho_sn(i,ksn) = 0.0_vpp
            ENDIF
        ENDIF
         IF(w_up .GT. w_res(ksn) .OR. excess_water .GT. 0.0_vpp) THEN
            ! ...water is being transfered
            dz_up       = dzm_sn_now(i,ksn)
            dz_low      = dzm_sn_now(i,ksn-1)
            w_low       = theta_w_now(i,ksn-1)
            dtheta_w_up = MAX(0.0_vpp, w_up - w_res(ksn))
            IF(dtheta_w_up .GT. 0.0_vpp .OR. excess_water .GT. 0.0_vpp) THEN
            ! ... dtheta_w_low is determined by also taking excess_water into
            ! account. Maybe excess_water can be stored in this layer.
            dtheta_w_low = dtheta_w_up * (dzm_sn_now(i,ksn)/dzm_sn_now(i,ksn-1)) + (excess_water/dzm_sn_now(i,ksn-1))
              ! now check whether there is enough air left - you might not be able
              ! to move the water
              ! or/and water may refreeze and expand specifically, you might want
              ! to create a water table over ice
              IF( (dtheta_w_low + w_low) .GT. (rho_i/rho_w * (1.0_vpp - theta_i_now(i,ksn-1))) ) THEN
                ! Deal with excess water ... Look how much you can leave in the
                ! lower layer (ksn-1).
                    ! If you have too much water even for the lower element (more
                    ! melt or rain per time
                    ! step than can be kept in this element), water is transferred
                    ! to excess_water.
                   ! excess_water moves the water downward, trying to insert the
                   ! water in lower layer.
                dtheta_w_low_x = dtheta_w_low    ! make backup
                dtheta_w_low   = MAX(0.0_vpp, (rho_i/rho_w * (1.0_vpp - theta_i_now(i,ksn-1)) - w_low))
                ! All the water that could not be stored in lower layer is
                ! considered excess_water.
                excess_water = (dtheta_w_low_x - dtheta_w_low) * dzm_sn_now(i,ksn-1)
              ELSE
                excess_water = 0.0_vpp
              ENDIF
              ! update volumetric contents, masses and density
              theta_w_now(i,ksn)   = w_up  - dtheta_w_up
              theta_w_now(i,ksn-1) = w_low + dtheta_w_low
            ENDIF ! end positive water movement
          ENDIF ! end if( W_upper > Wres )
     ENDDO ! loop snow layers
   
      ! ==========================
      ! Special treatment for lowermost snow layer, i.e runoff at bottom layer
      ! Note: If we want ponding on surface layer, i.e. glacier ice, sea ice, rock
      ! etc. we need to do it here, e.g. if itype = ...
      ! ==========================
      theta_w_bot = theta_w_now(i,1)
        IF (theta_w_bot .GT. w_res(1)) THEN
         ! Adjust dependent values accordingly
         theta_w_now(i,1) = w_res(1)
         theta_a_now(i,1) = 1.0_vpp - theta_w_now(i,1) - theta_i_now(i,1)
         ! Put all excess water of bottom layer  in runoff, i.e. move it out of the snow cover
         ! Note: if one comments this out you get ponding
         runoff_sn(i)              = runoff_sn(i) + excess_water +  dzm_sn_now(i,1) * (theta_w_bot - w_res(1))
        ENDIF
   
       ! Update dependent variables
       !CALL update(top(i), hm_sn(i,:), zm_sn(i,:), dzm_sn_now(i,:), m_sn(i,:)         , &
       !            theta_i_now(i,:), theta_w_now(i,:), theta_a_now(i,:),rho_sn(i,:)   , &
       !            hcap_sn(i,:), hcon_sn(i,:))

           ! -------------------------
           ! Height of snow (main) levels
           ! -------------------------
                 !$acc loop seq
                 do ksn=1,l_top
                    if(ksn .eq. 1) then
                       hm_sn(i,ksn) = dzm_sn_now(i,ksn)
                    else
                       hm_sn(i,ksn) = hm_sn(i,ksn-1) + dzm_sn_now(i,ksn)
                    endif
                 enddo  
         
           ! --------------------------
           ! Depth of snow layer (main) levels
           ! -------------------------
                 !$acc loop seq
                 do ksn = l_top,1,-1
                    zm_sn(i, (l_top+1) - ksn ) = hm_sn(i,ksn) !invert height vector
                    theta_a_now(i,ksn) = max(0.0_vpp,1.0_vpp - theta_i_now(i,ksn) - theta_w_now(i,ksn))
                 enddo
         
           ! --------------------------
           ! Snow layer density
           ! -------------------------
                  !$acc loop seq
                  do ksn = 1, l_top, 1
                     if(theta_i_now(i,ksn) .eq. 0.0_vpp) then
                        rho_sn(i,ksn) = 0.0_vpp
                     else
                        rho_sn(i,ksn) = theta_i_now(i,ksn)*rho_i + theta_w_now(i,ksn)*rho_w
                     endif
                  enddo
                  
           ! --------------------------
           ! Heat capacity
           ! --------------------------
            !$acc loop seq
            DO ksn = 1, l_top, 1
                IF(rho_sn(i,ksn) .LT. eps) THEN
                  hcap_sn(i,ksn) = 0.0_vpp
                ELSE
                  hcap_sn(i,ksn) = (  rho_a   * theta_a_now(i,ksn) * specific_heat_air     &
                                    + rho_i   * theta_i_now(i,ksn) * specific_heat_ice     &
                                    + rho_w   * theta_w_now(i,ksn) * specific_heat_water)  &
                                     / rho_sn(i,ksn)
                ENDIF
            ENDDO
         
           ! --------------------------
           ! Heat conductivity
           ! --------------------------
            !$acc loop seq
            DO ksn = 1, l_top, 1
               hcon_sn(i,ksn) = 2.22_vpp * EXP(1.88_vpp * LOG(rho_sn(i,ksn)/rho_i))
            ENDDO
         
           ! --------------------------
           ! Snow layer mass
           ! --------------------------
            !$acc loop seq
            DO ksn = 1, l_top, 1
              IF(ksn .EQ. 1) THEN
                 m_sn(i,ksn) = zm_sn(i,ksn) * ((theta_i_now(i,ksn) * rho_i) + (theta_w_now(i,ksn) * rho_w))
               ELSE
                 m_sn(i,ksn) = ABS(zm_sn(i,ksn) - zm_sn(i,ksn-1))  * ((theta_i_now(i,ksn) * rho_i) + (theta_w_now(i,ksn) * rho_w))
               ENDIF
            ENDDO

 



    ENDIF

  ENDDO
  !$acc end parallel


! =============================================================================
! - End Section IV: Water transport
! =============================================================================

if(ntstep .eq. N_IMP) then
   do ksn = 1,top(1)
     write(*,*) 'end section IV: ',ksn,theta_i_now(1,ksn),theta_a_now(1,ksn),theta_w_now(1,ksn),t_sn_now(1,ksn)
   enddo
endif


! =============================================================================
! + Begin Section V: Settling
! =============================================================================





  ! Check for layers that have been melted during phase change and take action
  ! -----------------------------------------------

  !$acc parallel async
  !$acc loop gang vector
  DO i = ivstart, ivend
    l_top = top(i)
    IF(l_top .GE. 2) THEN ! snow on the ground

      ! Aggregate profile if feasible
      !CALL aggregate_layers(top(i), hm_sn(i,:), zm_sn(i,:), dzm_sn_now(i,:), m_sn(i,:)            , &
      !                      theta_i_now(i,:), theta_w_now(i,:), theta_a_now(i,:)                  , &
      !                      t_sn_now(i,:), rho_sn(i,:), hcap_sn(i,:), hcon_sn(i,:), hdif_sn(i,:)    )

      ! First loop: eliminiate 'empty' layers
      ! -------------------
  
       DO ksn = l_top, 1, -1
          IF(theta_i_now(i,ksn) .LT. 0.01_vpp) THEN

!             write(*,*) 'OH MY E: ',ntstep,i,ksn,theta_i_now(i,ksn),theta_a_now(i,ksn),theta_w_now(i,ksn),l_top,runoff_sn(1)

             ! ... and reset roperties
             hm_sn(i,ksn)   = 0.0_vpp
             zm_sn(i,ksn)   = 0.0_vpp
             dzm_sn_now(i,ksn)  = 0.0_vpp
  
             theta_i_now(i,ksn) = 0.0_vpp
             theta_w_now(i,ksn) = 0.0_vpp
             theta_a_now(i,ksn) = 0.0_vpp

             rho_sn(i,ksn)  = 0.0_vpp
             t_sn_now(i,ksn)    = 0.0_vpp
             m_sn(i,ksn)    = 0.0_vpp
 
             hcap_sn(i,ksn) = 0.0_vpp
             hcon_sn(i,ksn) = 0.0_vpp
             hdif_sn(i,ksn) = 0.0_vpp
 
            ! Update number of snow layers ...
             l_top         = l_top - 1
             
          ENDIF
        ENDDO

   !do ksn = 1,l_top
   !  write(*,*) 'middle of section V: ',ksn,theta_i_now(1,ksn),theta_a_now(1,ksn),theta_w_now(1,ksn),t_sn_now(1,ksn)

   !enddo

  
      ! Second loop: aggregate layers if possible
      ! -------------------
  
       DO ksn = l_top, 2, -1
         IF(dzm_sn_now(i,ksn) .LT. min_height_layer .OR. theta_i_now(i,ksn) .LT. 0.01_vpp) THEN ! layer is quite thin
  
            ! -------------------
            ! Aggregate with lower layer - Adjust properties
            ! -------------------
  
  
             theta_i_now(i,ksn-1) = (theta_i_now(i,ksn)*dzm_sn_now(i,ksn)                                      &
                            + theta_i_now(i,ksn-1)*dzm_sn_now(i,ksn-1)) / (dzm_sn_now(i,ksn) + dzm_sn_now(i,ksn-1))  ! volumetric ice content
  
             theta_w_now(i,ksn-1) = (theta_w_now(i,ksn)*dzm_sn_now(i,ksn)                                      &
                            + theta_w_now(i,ksn-1)*dzm_sn_now(i,ksn-1)) / (dzm_sn_now(i,ksn) + dzm_sn_now(i,ksn-1))  ! volumetric water content
  
             dzm_sn_now(i,ksn-1)  = dzm_sn_now(i,ksn) + dzm_sn_now(i,ksn-1)        ! layer thickness
 
            ! ... and reset roperties
             hm_sn(i,ksn)   = 0.0_vpp
             zm_sn(i,ksn)   = 0.0_vpp
             dzm_sn_now(i,ksn)  = 0.0_vpp
  
             theta_i_now(i,ksn) = 0.0_vpp
             theta_w_now(i,ksn) = 0.0_vpp
             theta_a_now(i,ksn) = 0.0_vpp
  
             rho_sn(i,ksn)  = 0.0_vpp
             t_sn_now(i,ksn)    = 0.0_vpp
             m_sn(i,ksn)    = 0.0_vpp
  
             hcap_sn(i,ksn) = 0.0_vpp
             hcon_sn(i,ksn) = 0.0_vpp
             hdif_sn(i,ksn) = 0.0_vpp
  
             ! Update number of snow layers ...
              l_top         = l_top - 1

         ENDIF
       ENDDO
    !do ksn = 1,l_top
    ! write(*,*) 'middle of section V B: ',ksn,theta_i_now(1,ksn),theta_a_now(1,ksn),theta_w_now(1,ksn),t_sn_now(1,ksn)
    !enddo
 
!  
!      ! Third loop: it might be that subsurface layers melted, if so shift layer boundaries accordingly
!      ! -------------------
!  
       DO lay = 1, n_layers, 1
         IF(dzm_sn_now(i,lay) .LT. min_height_layer .OR. theta_i_now(i,lay) .LT. 0.01_vpp) THEN
           DO j = lay+1, n_layers, 1  
             dzm_sn_now(i,j-1)  = dzm_sn_now(i,j)
             theta_i_now(i,j-1) = theta_i_now(i,j)
             theta_w_now(i,j-1) = theta_w_now(i,j)
             rho_sn(i,j-1)      = rho_sn(i,j)
          ENDDO
         ENDIF
       ENDDO
  
       ! Sepecial treatment of top layer
       IF(dzm_sn_now(i,l_top) .LT. min_height_layer .OR. theta_i_now(i,l_top) .LT. 0.01_vpp) THEN
           l_top = l_top - 1
       ENDIF
  
  
      ! Clear/clean properties
      !CALL clear_layers(top(i), hm_sn(i,:), zm_sn(i,:), dzm_sn_now(i,:), m_sn(i,:)                , &
      !                  theta_i_now(i,:), theta_i_old(i,:), theta_w_now(i,:), theta_a_now(i,:)    , &
      !                  t_sn_now(i,:), rho_sn(i,:), hcap_sn(i,:), hcon_sn(i,:), hdif_sn(i,:)        )

        DO ksn = l_top+1, n_layers, 1
          hm_sn(i,ksn)   = 0.0_vpp
          zm_sn(i,ksn)   = 0.0_vpp
          dzm_sn_now(i,ksn)  = 0.0_vpp
          theta_i_now(i,ksn) = 0.0_vpp
          theta_w_now(i,ksn) = 0.0_vpp
          theta_a_now(i,ksn) = 0.0_vpp
          theta_i_old(i,ksn) = 0.0_vpp
          rho_sn(i,ksn)  = 0.0_vpp
          t_sn_now(i,ksn)    = 0.0_vpp
          m_sn(i,ksn)    = 0.0_vpp
          hcap_sn(i,ksn) = 0.0_vpp
          hcon_sn(i,ksn) = 0.0_vpp
          hdif_sn(i,ksn) = 0.0_vpp
        ENDDO
     
       ! Update dependent variables
       !CALL update(top(i), hm_sn(i,:), zm_sn(i,:), dzm_sn_now(i,:), m_sn(i,:)                     , &
       !            theta_i_now(i,:), theta_w_now(i,:), theta_a_now(i,:),rho_sn(i,:)               , &
       !            hcap_sn(i,:), hcon_sn(i,:)                                                       )

           ! -------------------------
           ! Height of snow (main) levels
           ! -------------------------
                 !$acc loop seq
                 do ksn=1,l_top
                    if(ksn .eq. 1) then
                       hm_sn(i,ksn) = dzm_sn_now(i,ksn)
                    else
                       hm_sn(i,ksn) = hm_sn(i,ksn-1) + dzm_sn_now(i,ksn)
                    endif
                 enddo  
         
           ! --------------------------
           ! Depth of snow layer (main) levels
           ! -------------------------
                 !$acc loop seq
                 do ksn = l_top,1,-1
                    zm_sn(i, (l_top+1) - ksn ) = hm_sn(i,ksn) !invert height vector
                    theta_a_now(i,ksn) = max(0.0_vpp,1.0_vpp - theta_i_now(i,ksn) - theta_w_now(i,ksn))
                 enddo
         
           ! --------------------------
           ! Snow layer density
           ! -------------------------
                  !$acc loop seq
                  do ksn = 1, l_top, 1
                     if(theta_i_now(i,ksn) .eq. 0.0_vpp) then
                        rho_sn(i,ksn) = 0.0_vpp
                     else
                        rho_sn(i,ksn) = theta_i_now(i,ksn)*rho_i + theta_w_now(i,ksn)*rho_w
                     endif
                  enddo
                  
           ! --------------------------
           ! Heat capacity
           ! --------------------------
            !$acc loop seq
            DO ksn = 1, l_top, 1
                IF(rho_sn(i,ksn) .LT. eps) THEN
                  hcap_sn(i,ksn) = 0.0_vpp
                ELSE
                  hcap_sn(i,ksn) = (  rho_a   * theta_a_now(i,ksn) * specific_heat_air     &
                                    + rho_i   * theta_i_now(i,ksn) * specific_heat_ice     &
                                    + rho_w   * theta_w_now(i,ksn) * specific_heat_water)  &
                                     / rho_sn(i,ksn)
                ENDIF
            ENDDO
         
           ! --------------------------
           ! Heat conductivity
           ! --------------------------
            !$acc loop seq
            DO ksn = 1, l_top, 1
               hcon_sn(i,ksn) = 2.22_vpp * EXP(1.88_vpp * LOG(rho_sn(i,ksn)/rho_i))
            ENDDO
         
           ! --------------------------
           ! Snow layer mass
           ! --------------------------
            !$acc loop seq
            DO ksn = 1, l_top, 1
              IF(ksn .EQ. 1) THEN
                 m_sn(i,ksn) = zm_sn(i,ksn) * ((theta_i_now(i,ksn) * rho_i) + (theta_w_now(i,ksn) * rho_w))
               ELSE
                 m_sn(i,ksn) = ABS(zm_sn(i,ksn) - zm_sn(i,ksn-1))  * ((theta_i_now(i,ksn) * rho_i) + (theta_w_now(i,ksn) * rho_w))
               ENDIF
            ENDDO


    ENDIF
    top(i) = l_top
  ENDDO
  !$acc end parallel


  ! Apply some fail-safes before doing the settling, which is temperature  dependent
  ! ------------------------------

  !$acc parallel async default(none)
  !$acc loop gang vector
  DO i = ivstart, ivend

    t_sn_sfc(i) = MAX(220.0_vpp, MIN(t_sn_sfc(i), t0_melt))

  ENDDO
  !$acc end parallel


  !$acc parallel async default(none)
  !$acc loop gang vector
  DO i = ivstart, ivend

    !DO ksn = top(i), 1, -1
    !
    !write(*,*) ntstep,i,ksn,t_sn_now(i,ksn) != MIN(MAX(220.0_vpp, t_sn_now(i,ksn)), t0_melt)
    !
    !ENDDO    

    !$acc loop seq
    DO ksn = top(i), 1, -1

      t_sn_now(i,ksn) = MIN(MAX(220.0_vpp, t_sn_now(i,ksn)), t0_melt)

    ENDDO    

  ENDDO
  !$acc end parallel


  ! Calculate settling rates
  ! ---------------------------------

  !$acc parallel async
  !$acc loop gang vector
  DO i = ivstart, ivend
    l_top = top(i)
    IF(l_top .GE. 1) THEN ! snow on the ground

      !CALL calc_set_rate(h_snow(i), zm_sn(i,:), dzm_sn_now(i,:), t_sn_now(i,:), m_sn(i,:)        , &
      !                   theta_i_now(i,:), theta_w_now(i,:), theta_a_now(i,:), theta_i_old(i,:)  , &
      !                   rho_sn(i,:), top(i), melt_flag(i), zdt)


     ! Initalizations
     ddz          = 0.0_vpp   ! total change in layer thickness
     overburden   = 0.0_vpp   ! overburden stress (weight of the overlying layers times g)

     DO ksn = l_top, 1, -1

       ! Reset values for settling rates
       rate_1    = 0.0_vpp      ! due to ice loss
       rate_2    = 0.0_vpp      ! due to overburden stress

       tot_rate  = 0.0_vpp      ! total settliing rate

       dz_old    = 0.0_vpp      ! old layer thickness (before settling) ...
       dz_new    = 0.0_vpp      ! ... and after settling

       !  ... due to loss of ice
       !-------------------
 !      IF(melt_flag) THEN
         IF(theta_i_now(i,ksn) <  theta_i_old(i,ksn)) THEN ! melting ocured in that layers
           rate_1 = - 1.0_vpp/dt * MAX(0.0_vpp,  (theta_i_old(i,ksn) - theta_i_now(i,ksn)) / theta_i_old(i,ksn) )
         ENDIF
 !      ENDIF
       !  ... due to overburden stress
       !-------------------
      ! Vionnett (2012)
       f1   = 1.0_vpp / (1.0_vpp + 60.0_vpp * ((theta_w_now(i,ksn)*rho_w*dzm_sn_now(i,ksn)) / (rho_w * dzm_sn_now(i,ksn))))
       f2   = 1.0_vpp

       eta = f1 * f2 * eta_0 * (rho_sn(i,ksn)/c_eta) * exp(a_eta*(t0_melt - t_sn_now(i,ksn)) + b_eta*rho_sn(i,ksn))
       rate_2 = -1.0_vpp * (overburden + (m_sn(i,ksn)*g/2.0_vpp)) / eta

!       if(ntstep .eq. 165052) then

!           write(*,*) i,ksn,f1,f2,overburden,g,eta,eta_0,rho_sn(i,ksn),t_sn_now(i,ksn)
!       endif
       ! Andersen (1976)
!        rate_2 = -(overburden + m_sn(ksn)*g/2.0_vpp)*EXP(-c_factor * (t_sn(ksn) - t0_melt) - c2*rho_sn(ksn)) / eta0

       ! increase overburden stress NOTE: How to deal with slope angle should be overburden = m*g*cos(alpha)
       overburden = overburden + (m_sn(i,ksn) * g)

      ! ... calculate change ...
       !-------------------

       ! ... of all (ice loss, overburden, destructive) settling rates (1/s)
       tot_rate = (rate_1*dt) + (rate_2*dt)

       ! ... of layer thickness, i.e. sum over all layers (m)
       ddz          = ddz + MAX(-1.0_vpp * dzm_sn_now(i,ksn), dzm_sn_now(i,ksn) * tot_rate)

       dz_old = dzm_sn_now(i,ksn)
       dz_new = dz_old + MAX(-1.0_vpp * dzm_sn_now(i,ksn), dzm_sn_now(i,ksn) * tot_rate)

!       if (ntstep .eq. 165052) then
!          write(*,*) i,ksn,dzm_sn_now(i,ksn),dz_old,MAX(-1.0_vpp * dzm_sn_now(i,ksn), dzm_sn_now(i,ksn) * tot_rate),rate_1,rate_2,dt
!       endif

       ! ... volumetric contents
       theta_i_now(i,ksn) = MAX(0.0_vpp, theta_i_now(i,ksn) * (dz_old / dz_new))    ! ice content
       theta_w_now(i,ksn) = MAX(0.0_vpp, theta_w_now(i,ksn) * (dz_old / dz_new))    ! water content

       ! ... of layer thickness (m)
       dzm_sn_now(i,ksn) = dz_new

  ENDDO

  ! -------------------
  ! + Re-calculate layer depth from layer thickness
  ! -------------------

    DO ksn = l_top, 1, -1
      IF(ksn .EQ. l_top) THEN
        zm_sn(i,ksn) = dzm_sn_now(i,ksn)
      ELSE
        zm_sn(i,ksn) = zm_sn(i,ksn+1) + dzm_sn_now(i,ksn)
      ENDIF
    ENDDO

    ! Calculate change in snow depth
    h_snow(i) = zm_sn(i,1)
    ENDIF
  ENDDO
  !$acc end parallel

! =============================================================================
! - End Section V: Settling
! =============================================================================

if(ntstep .eq. N_IMP) then
   do ksn = 1,top(1)
     write(*,*) 'end section V: ',ksn,theta_i_now(1,ksn),theta_a_now(1,ksn),theta_w_now(1,ksn),t_sn_now(1,ksn)

   enddo
endif





! =============================================================================
! + Begin Section VI: Preparations for the next time step 
! =============================================================================



   ! Check for layers that have been melted during phase change
   ! ------------------------------------

    !$acc parallel async
    !$acc loop gang vector
    DO i = ivstart, ivend
      l_top = top(i)

      IF(l_top .GE. 2) THEN ! snow on the ground

!      ! Aggregate profile if feasible
!      !CALL aggregate_layers(top(i), hm_sn(i,:), zm_sn(i,:), dzm_sn_now(i,:), m_sn(i,:)            , &
!      !                      theta_i_now(i,:), theta_w_now(i,:), theta_a_now(i,:)                  , &
!      !                      t_sn_now(i,:), rho_sn(i,:), hcap_sn(i,:), hcon_sn(i,:), hdif_sn(i,:)    )

!      ! First loop: eliminiate 'empty' layers
!      ! -------------------
!       DO ksn = l_top, 1, -1
!          IF(theta_i_now(i,ksn) .LT. 0.01_vpp) THEN
!             ! ... and reset roperties
!             hm_sn(i,ksn)   = 0.0_vpp
!             zm_sn(i,ksn)   = 0.0_vpp
!             dzm_sn_now(i,ksn)  = 0.0_vpp
!  
!             theta_i_now(i,ksn) = 0.0_vpp
!             theta_w_now(i,ksn) = 0.0_vpp
!             theta_a_now(i,ksn) = 0.0_vpp
!  
!             rho_sn(i,ksn)  = 0.0_vpp
!             t_sn_now(i,ksn)    = 0.0_vpp
!             m_sn(i,ksn)    = 0.0_vpp
!  
!             hcap_sn(i,ksn) = 0.0_vpp
!             hcon_sn(i,ksn) = 0.0_vpp
!             hdif_sn(i,ksn) = 0.0_vpp
!  
!             ! Update number of snow layers ...
!             l_top         = l_top - 1
!             write(*,*) 'OH MY G: ', l_top 

!          ENDIF
!        ENDDO
  
      ! Second loop: aggregate layers if possible
      ! -------------------
!        DO ksn = l_top, 2, -1
!          IF(dzm_sn_now(i,ksn) .LT. min_height_layer .OR. theta_i_now(i,ksn) .LT. 0.01_vpp) THEN ! layer is quite thin
!  
!            ! -------------------
!            ! Aggregate with lower layer - Adjust properties
!            ! -------------------
!  
!             dzm_sn_now(i,ksn-1)  = dzm_sn_now(i,ksn) + dzm_sn_now(i,ksn-1)        ! layer thickness
!  
!             theta_i_now(i,ksn-1) = (theta_i_now(i,ksn)*dzm_sn_now(i,ksn)                                      &
!                            + theta_i_now(i,ksn-1)*dzm_sn_now(i,ksn-1)) / (dzm_sn_now(i,ksn) + dzm_sn_now(i,ksn-1))  ! volumetric ice content
!  
!             theta_w_now(i,ksn-1) = (theta_w_now(i,ksn)*dzm_sn_now(i,ksn)                                      &
!                            + theta_w_now(i,ksn-1)*dzm_sn_now(i,ksn-1)) / (dzm_sn_now(i,ksn) + dzm_sn_now(i,ksn-1))  ! volumetric water content
!  
!  
!            ! ... and reset roperties
!             hm_sn(i,ksn)   = 0.0_vpp
!             zm_sn(i,ksn)   = 0.0_vpp
!             dzm_sn_now(i,ksn)  = 0.0_vpp
!  
!             theta_i_now(i,ksn) = 0.0_vpp
!             theta_w_now(i,ksn) = 0.0_vpp
!             theta_a_now(i,ksn) = 0.0_vpp
!  
!             rho_sn(i,ksn)  = 0.0_vpp
!             t_sn_now(i,ksn)    = 0.0_vpp
!             m_sn(i,ksn)    = 0.0_vpp
!  
!             hcap_sn(i,ksn) = 0.0_vpp
!             hcon_sn(i,ksn) = 0.0_vpp
!             hdif_sn(i,ksn) = 0.0_vpp
!  
!             ! Update number of snow layers ...
!              l_top         = l_top - 1
!              write(*,*) 'OH MY K:', l_top 
!          ENDIF
!  
!        ENDDO
  
  
      ! Third loop: it might be that subsurface layers melted, if so shift layer boundaries accordingly
      ! -------------------
  
!       DO lay = 1, n_layers, 1
!         IF(dzm_sn_now(i,lay) .LT. min_height_layer .OR. theta_i_now(i,lay) .LT. 0.01_vpp) THEN
!           DO j = lay+1, n_layers, 1
!  
!             dzm_sn_now(i,j-1)  = dzm_sn_now(i,j)
!             theta_i_now(i,j-1) = theta_i_now(i,j)
!             theta_w_now(i,j-1) = theta_w_now(i,j)
!             rho_sn(i,j-1)      = rho_sn(i,j)
!  
!           ENDDO
!         ENDIF
!       ENDDO
  
         ! Sepecial treatment of top layer
!         IF(dzm_sn_now(i,l_top) .LT. min_height_layer .OR. theta_i_now(i,l_top) .LT. 0.01_vpp) THEN
!           l_top = l_top - 1
!         ENDIF
  
  
      ! Clear/clean properties
      !CALL clear_layers(top(i), hm_sn(i,:), zm_sn(i,:), dzm_sn_now(i,:), m_sn(i,:)                , &
      !                  theta_i_now(i,:), theta_i_old(i,:), theta_w_now(i,:), theta_a_now(i,:)    , &
      !                  t_sn_now(i,:), rho_sn(i,:), hcap_sn(i,:), hcon_sn(i,:), hdif_sn(i,:)        )

!        DO ksn = l_top+1, n_layers, 1
!          hm_sn(i,ksn)   = 0.0_vpp
!          zm_sn(i,ksn)   = 0.0_vpp
!          dzm_sn_now(i,ksn)  = 0.0_vpp
!          theta_i_now(i,ksn) = 0.0_vpp
!          theta_w_now(i,ksn) = 0.0_vpp
!          theta_a_now(i,ksn) = 0.0_vpp
!          theta_i_old(i,ksn) = 0.0_vpp
!          rho_sn(i,ksn)  = 0.0_vpp
!          t_sn_now(i,ksn)    = 0.0_vpp
!          m_sn(i,ksn)    = 0.0_vpp
!          hcap_sn(i,ksn) = 0.0_vpp
!          hcon_sn(i,ksn) = 0.0_vpp
!          hdif_sn(i,ksn) = 0.0_vpp
!        ENDDO
     
       ! Update dependent variables
       !CALL update(top(i), hm_sn(i,:), zm_sn(i,:), dzm_sn_now(i,:), m_sn(i,:)                     , &
       !            theta_i_now(i,:), theta_w_now(i,:), theta_a_now(i,:),rho_sn(i,:)               , &
       !            hcap_sn(i,:), hcon_sn(i,:)                                                       )

           ! -------------------------
           ! Height of snow (main) levels
           ! -------------------------
                 !$acc loop seq
                 do ksn=1,l_top
                    if(ksn .eq. 1) then
                       hm_sn(i,ksn) = dzm_sn_now(i,ksn)
                    else
                       hm_sn(i,ksn) = hm_sn(i,ksn-1) + dzm_sn_now(i,ksn)
                    endif
                 enddo  
         
           ! --------------------------
           ! Depth of snow layer (main) levels
           ! -------------------------
                 !$acc loop seq
                 do ksn = l_top,1,-1
                    zm_sn(i, (l_top+1) - ksn ) = hm_sn(i,ksn) !invert height vector
                    theta_a_now(i,ksn) = max(0.0_vpp,1.0_vpp - theta_i_now(i,ksn) - theta_w_now(i,ksn))
                 enddo
         
           ! --------------------------
           ! Snow layer density
           ! -------------------------
                  !$acc loop seq
                  do ksn = 1, l_top, 1
                     if(theta_i_now(i,ksn) .eq. 0.0_vpp) then
                        rho_sn(i,ksn) = 0.0_vpp
                     else
                        rho_sn(i,ksn) = theta_i_now(i,ksn)*rho_i + theta_w_now(i,ksn)*rho_w
                     endif
                  enddo
                  
           ! --------------------------
           ! Heat capacity
           ! --------------------------
            !$acc loop seq
            DO ksn = 1, l_top, 1
                IF(rho_sn(i,ksn) .LT. eps) THEN
                  hcap_sn(i,ksn) = 0.0_vpp
                ELSE
                  hcap_sn(i,ksn) = (  rho_a   * theta_a_now(i,ksn) * specific_heat_air     &
                                    + rho_i   * theta_i_now(i,ksn) * specific_heat_ice     &
                                    + rho_w   * theta_w_now(i,ksn) * specific_heat_water)  &
                                     / rho_sn(i,ksn)
                ENDIF
            ENDDO
         
           ! --------------------------
           ! Heat conductivity
           ! --------------------------
            !$acc loop seq
            DO ksn = 1, l_top, 1
               hcon_sn(i,ksn) = 2.22_vpp * EXP(1.88_vpp * LOG(rho_sn(i,ksn)/rho_i))
            ENDDO
         
           ! --------------------------
           ! Snow layer mass
           ! --------------------------
            !$acc loop seq
            DO ksn = 1, l_top, 1
              IF(ksn .EQ. 1) THEN
                 m_sn(i,ksn) = zm_sn(i,ksn) * ((theta_i_now(i,ksn) * rho_i) + (theta_w_now(i,ksn) * rho_w))
               ELSE
                 m_sn(i,ksn) = ABS(zm_sn(i,ksn) - zm_sn(i,ksn-1))  * ((theta_i_now(i,ksn) * rho_i) + (theta_w_now(i,ksn) * rho_w))
               ENDIF
            ENDDO


      ENDIF

      top(i) = l_top
    ENDDO
    !$acc end parallel

if(ntstep .eq. N_IMP) then
   do ksn = 1,top(1)
     write(*,*) 'end section VI: ',ksn,theta_i_now(1,ksn),theta_a_now(1,ksn),theta_w_now(1,ksn),t_sn_now(1,ksn)
   enddo
endif

! =============================================================================
! + Begin Section VII: Updating of prognostic variables
! =============================================================================


    ! Some resetting in case of complete melt out
    ! -------------------

    !$acc parallel async default(none)
    !$acc loop gang vector
   ! DO i = ivstart, ivend
  
   !   ! In case snow melted completly also reset indices
   !   IF(h_snow(i) .LT. eps) THEN

   !     top(i)    = 0
   !     h_snow(i) = 0.0_vpp
   !     swe_sn(i)    = 0.0_vpp
   !

   !    ENDIF

   !  ENDDO
     !$acc end parallel


    ! Update snow layer properties
    ! -------------------------------------

    !$acc parallel async default(none)
    !$acc loop seq
    DO ksn = 1, n_layers, 1
      DO i = ivstart, ivend

        t_sn_new(i,ksn)    = t_sn_now(i,ksn)     ! snow layer temperature  

        theta_i_new(i,ksn) = theta_i_now(i,ksn)  ! volumetric ice content
        theta_w_new(i,ksn) = theta_w_now(i,ksn)  ! volumetric water content
        theta_a_new(i,ksn) = theta_a_now(i,ksn)  ! volumetric air content

        dzm_sn_new(i,ksn)  = dzm_sn_now(i,ksn)   ! layer thickness

      END DO
    END DO
    !$acc end parallel


    ! Calculate snow water equivialent
    ! -------------------------------------



    !$acc parallel async default(none)
    !$acc loop seq
    DO ksn = 1, n_layers, 1
      DO i = ivstart, ivend

       swe_sn(i) = swe_sn(i) + ( dzm_sn_now(i,ksn) * rho_sn(i,ksn) ) 

      END DO
    ENDDO
    !$acc end parallel



    ! Update snow height and surface temperature
    ! -----------------------------

    !$acc parallel async default(none)
    !$acc loop gang vector
    DO i = ivstart, ivend
      top_sn_new (i) = REAL(top(i))    ! index of first (top) snow layer
      hn_sn_new  (i) = hn_sn_now(i)    ! new snow storage

      top_sn_now (i) = REAL(top(i))    ! index of first (top) snow layer

      IF(top(i) .GE. 1) THEN
        h_snow     (i) = hm_sn(i,top(i))   ! snow height NOTE: for this line to work
                                        !                   the snow height update needs to be commented
                                        !                   out in sfc_terra.f90
                                        !                   or #ifdef later on
        !write(*,*) 'real: ', h_snow(i)
      ELSE
         h_snow(i) = 0.0_vpp
      ENDIF

      if (h_snow(i) .gt. 0.0) then
        rho_snow_new(i) = swe_sn(i)/h_snow(i)
      !t_snow_new (i) =  t_sn_sfc(i)     ! snow surface temperature
      
        tmp = 0.0_vpp
        do ksn = 1,top(i)
           tmp = tmp + dzm_sn_now(i,ksn) * t_sn_now(i,ksn)
        enddo  
      
      t_snow_new (i) = tmp / h_snow(i) !t_sn_now(i,1) 
      w_snow_new (i) = swe_sn(i)/1000.0_vpp

      else
        !rho_snow_new(i) = 0.0_vpp
      endif
    !write(*,*) 'update: ',top_sn_new(i),hn_sn_new(i),h_snow(i),t_snow_new(i),w_snow_new(i)
    END DO
    !$acc end parallel


    ! Update turbulent fluxes
    ! -----------------------------

    !$acc parallel async default(none)
    !$acc loop gang vector
    DO i = ivstart, ivend

      zshfl_snow(i) = sh_sn(i)
      zlhfl_snow(i) = lh_sn(i)
 !     zqhfl_snow(i) = lh_sn(i)/lh_s
    
    ENDDO
    !$acc end parallel


!    write(*,*) 'snowpolino end: ',top(1)
!    write(*,*) 'hsnow: ', ntstep,h_snow(1),t_so_now(1,1),sobs(1),prs_gsp,prr_gsp
!     write(*,*) 'energy: ',ntstep,h_snow(1),t_so_now(1,1),sobs(1)
!    END IF
!   write(*,*) 'end: ',ntstep,top(1),h_snow(1) 

if(ntstep .eq. N_IMP) then
   do ksn = 1,top(1)
     write(*,*) 'end section VII: ',ksn,theta_i_now(1,ksn),theta_a_now(1,ksn),theta_w_now(1,ksn),t_sn_now(1,ksn)
   enddo
endif


    if( mod(ntstep,120) .eq. 0) then 
       open(22, file = './out/pro_snow.txt')
       write(22,'(A7,   2I12.5,a,/)',    ADVANCE = 'YES') 'idx'         , ntstep, top(1)
       write(22,'(A7,  1000F12.3,a,/)',   ADVANCE = 'YES') 'hm_sn'      , hm_sn(1,:)
       write(22,'(A7,  1000F12.3,a,/)',   ADVANCE = 'YES') 'zm_sn'      , zm_sn(1,:)
       write(22,'(A7,  1000F12.3,a,/)',   ADVANCE = 'YES') 'dz'         , dzm_sn_now(1,:)
       write(22,'(A7,  1000F12.3,a,/)',   ADVANCE = 'YES') 't_sn'       , t_sn_now(1,:)
       write(22,'(A,   1000F12.3,a,/)',   ADVANCE = 'YES') 'theta_i'    , theta_i_now(1,:)
       write(22,'(A,   1000F12.3,a,/)',   ADVANCE = 'YES') 'theta_w'    , theta_w_now(1,:)
       write(22,'(A,   1000F12.3,a,/)',   ADVANCE = 'YES') 'theta_a'    , theta_a_now(1,:)
       write(22,'(A7,  1000F12.3,a,/)',   ADVANCE = 'YES') 'rho_sn'     , rho_sn(1,:)
       write(22,'(A7,  1000F12.3,a,/)',   ADVANCE = 'YES') 'm_sn'       , m_sn(1,:)
       write(22,'(A,   1000F12.3,a,/)',   ADVANCE = 'YES') 'hcap_sn'    , hcap_sn(1,:)
       write(22,'(A,   1000F12.3,a,/)',   ADVANCE = 'YES') 'hcon_sn'    , hcon_sn(1,:)
       write(22,'(A,   1000F12.3,a,/)',   ADVANCE = 'YES') 'hdif_sn'    , hdif_sn(1,:)
    endif

    if( mod(ntstep,120) .eq. 0) then
       open(44, file = './out/met.txt')  
       write(44,'(1I6, 13F18.3,a,/)', ADVANCE = 'YES') ntstep,t(1),t_sn_sfc(1),qv(1),zuv,& 
&                           zrain_rate,sobs(1),thbs(1),lh_sn(1),sh_sn(1),h_snow(1),0.9
       !close(44)
    endif

    if(ntstep .gt. N_IMP) then
      stop
    endif

! ------------------------------------------------------------------------------
! - End Section '?': Update prognostic variables
! ------------------------------------------------------------------------------

!$acc end data
!$acc end data


END SUBROUTINE snowpolino

! =============================================================================
!  - End subroutine  - snowpolino
! =============================================================================


! DONE, DONE!!!
END MODULE sfc_snow

!==============================================================================
! + End of module src_soil_multlay
!==============================================================================







