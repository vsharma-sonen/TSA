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

use data_modelconfig, only: ke_snow
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

                        ! surface
                        t_snow_now       , & ! temperature of the snow-surface              (  K  )
                        t_snow_new       , & ! temperature of the snow-surface              (  K  )
                        zshfl_snow       , & ! sensible heat flux snow/air interface        ( W/m2)
                        zlhfl_snow       , & ! latent   heat flux snow/air interface        ( W/m2)

                        ! snow
                        t_sn_now         , & ! snow temperature (nodes)                (  K  )
                        t_sn_new         , & ! snow temperature (nodes)                (  K  )
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
                        t_so_now         , & 
                        t_so_new           &
                        ) 

        
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


   REAL    (KIND = vpp)    , DIMENSION(nvec,1:ke_snow), INTENT(INOUT) :: &
                  theta_i_now       , &   ! volumetric ice content                     (-)
                  theta_w_now       , &   ! water ice content                          (-)
                  theta_a_now       , &   ! air ice content                            (-)
                  dzm_sn_now              ! layer thickness between main levels        (m)

   REAL    (KIND = vpp)    , DIMENSION(nvec,1:ke_snow), INTENT(OUT) :: &
                  theta_i_new       , &   ! volumetric ice content                     (-)
                  theta_w_new       , &   ! water ice content                          (-)
                  theta_a_new       , &   ! air ice content                            (-)
                  dzm_sn_new              ! layer thickness between main levels        (m)

   REAL    (KIND = vpp)    , DIMENSION(nvec,1:ke_snow+1), INTENT(INOUT) :: t_sn_now ! NODAL TEMPERATURES
 
   REAL    (KIND = vpp)    , DIMENSION(nvec,1:ke_snow+1), INTENT(OUT)   :: t_sn_new ! UPDATED NODAL TEMPERATURES

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

! ------------------------
! + Local
! ------------------------

   INTEGER      :: &

     ! Indices
     i               , & ! loop index in x-drection
     ksn             , & ! loop index for snow layers
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
     swe_sn       , &  ! local S.W.E calculation
     delta_rad         ! delta coefficient for linearizing net longwave radiation
 

   INTEGER, DIMENSION(nvec)          :: &

     top        , &      ! index of first (top) snow layer  [-] 
     merge_point

   REAL   (KIND=vpp) ::  &

     rho_sn      (nvec,1:ke_snow)  , & ! density of snow layers                               (kg/m**3)
     hcap_sn     (nvec,1:ke_snow)  , & ! snow layer heat capacity
     hcon_sn     (nvec,1:ke_snow)  , & !            heat conductivity
     hdif_sn     (nvec,1:ke_snow)  , & !            heat diffusion

     swabs_sn    (nvec,1:ke_snow)  , & ! absorbed short wave radiation

     theta_i_old (nvec,1:ke_snow)  , & !              volumetric ice content
     theta_w_old (nvec,1:ke_snow)      !              volumetric water content

    LOGICAL, DIMENSION(nvec)  :: &

      melt_flag     , & ! flag indicating melting; IF THEN TRUE
      freeze_flag       !      indicating freezing

    REAL    (KIND=vpp) :: &

      tmp_sn (ke_snow)         ! temporary vector      


#endif

  REAL (KIND=vpp)       ::  &
    dlw_u_sn              , &  ! derivative of upwelling longwave radiation for snow on the ground    (W m-2)
    dz_up                 , &  ! thickness above the layer of interest                                (m)
    dz_low                , &  ! thickness below the layer of interest                                (m)
    beta

  REAL (KIND=vpp), DIMENSION(1:ke_snow+1) :: &
    a_matrix                    , &  !
    b_matrix                    , &  !
    c_matrix                    , &  !
    d_matrix

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

  REAL    (KIND = vpp), DIMENSION(1:ke_snow) :: &

    w_res           , &    ! effective residual water content
    res_wat_cont    , &    ! potential residual water content
    dtheta_w_sub_wtr       ! additional storage capacity due to refreezing

  REAL ( KIND = vpp), DIMENSION(nvec,1:ke_snow) :: &
    t_sn_elem


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

   integer(kind=vpp) :: elems_left,l_mergep,l_thinlayer_merge,merge_elems
   integer(kind=vpp) :: num_new_elems

   real(kind=vpp), parameter :: new_snow_elem = 0.01
#ifdef __ICON__
  INTEGER :: my_cart_id
#endif

  INTEGER :: my_thrd_id, mcid, mtid, mbid, mvid, k

  LOGICAL ::  ldebug = .false.

  INTEGER, PARAMETER :: N_IMP = 999999999 !200000 !8659 !125765 !125761

  real (kind = vpp) :: l_hsnow,l_mass,l_lh_sn
  real (kind = vpp) :: tmp,tmp_rate_1,tmp_rate_2

  real (kind = vpp) :: alpha_solver_down,alpha_solver_up
  real (kind = vpp) :: beta_solver_down,beta_solver_up
  real (kind = vpp) :: coeff
  real (kind = vpp) :: t_emiss, emiss
  real (kind = vpp) :: lower_bc
! ------------------------------------------------------------------------------
! - End declarations
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Begin Section NULLA: Initialisation(s)
! ------------------------------------------------------------------------------

! Just do some checkout prints:
  IF (ldebug) THEN
      WRITE(*,'(A,3I5)'   ) 'SFC-DIAGNOSIS snowpolino start:   ', ke_snow, ke_snow

    DO i = ivstart, ivend
        WRITE(*,'(A,2I5)'  ) 'SFC-DIAGNOSIS snow:   iblock = ', iblock, i
        WRITE(*,'(A      )') ' External Parameters:  '

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
        WRITE(*,'(A,F28.16)') '   t_snow_now       :  ', t_snow_now  (i)
        WRITE(*,'(A,F28.16)') '   h_snow           :  ', h_snow      (i)
        WRITE(*,'(A,F28.16)') '   w_snow_now       :  ', w_snow_now  (i)
        WRITE(*,'(A,F28.16)') '   hn_sn_now        :  ', hn_sn_now   (i)
        WRITE(*,'(A,F28.16)') '   top_sn_now       :  ', top_sn_now  (i)

        WRITE(*,'(A      )') ' Multi level parameters:'
        do k = 1, ke_snow+1
             WRITE(*,'(A,I1,A,F28.16)') '   t_sn_now   (',k,')   :  ', t_sn_now    (i,k)
        enddo
        do k = 1, ke_snow
             WRITE(*,'(A,I1,A,F28.16)') '   theta_i_now(',k,')   :  ', theta_i_now (i,k)
        enddo
        do k = 1, ke_snow
             WRITE(*,'(A,I1,A,F28.16)') '   theta_w_now(',k,')   :  ', theta_w_now (i,k)
        enddo
        do k = 1, ke_snow
             WRITE(*,'(A,I1,A,F28.16)') '   theta_a_now(',k,')   :  ', theta_a_now (i,k)
        enddo
        do k = 1, ke_snow
             WRITE(*,'(A,I1,A,F28.16)') '   dzm_sn_now (',k,')   :  ', dzm_sn_now  (i,k)
        enddo
    ENDDO
  ENDIF


   hdif_sn = 0.0_vpp
   zdt      = dt


  ! ------------------------
  ! Intiate local fields
  ! -----------------------

   DO i = ivstart, ivend
     top(i) = NINT(top_sn_now(i)) 
   ENDDO


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

   DO i = ivstart, ivend
      DO ksn = 1, ke_snow
         rho_sn  (i,ksn) = 0.0_vpp
         hcon_sn (i,ksn) = 0.0_vpp
         hcap_sn (i,ksn) = 0.0_vpp
         swabs_sn(i,ksn) = 0.0_vpp
      ENDDO
   ENDDO

   DO i = ivstart, ivend

     l_top = top(i)
     IF(l_top .GE. 1) THEN

         theta_a_now(i,:) = max(0.0_vpp,1.0_vpp - theta_i_now(i,:) - theta_w_now(i,:))

         ! --------------------------
         ! Snow layer density
         ! -------------------------
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
         do ksn = 1, l_top, 1
             if(rho_sn(i,ksn) .LT. eps) then
               hcap_sn(i,ksn) = 0.0_vpp
             else
               hcap_sn(i,ksn) = (  rho_a   * theta_a_now(i,ksn) * specific_heat_air     &
                                 + rho_i   * theta_i_now(i,ksn) * specific_heat_ice     &
                                 + rho_w   * theta_w_now(i,ksn) * specific_heat_water)  &
                                  / rho_sn(i,ksn)
             endif
         enddo

         ! --------------------------
         ! Heat conductivity
         ! --------------------------
         do ksn = 1, l_top, 1
             hcon_sn(i,ksn) = 2.22_vpp * EXP(1.88_vpp * LOG(rho_sn(i,ksn)/rho_i))
         enddo

     ENDIF ! if (l_top .ge. 1) 
   ENDDO ! loop over grid points

! ------------------------------------------------------------------------------
! - End Section Initializations
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Begin Section I: New snow
! ------------------------------------------------------------------------------

  DO i = ivstart, ivend

    l_top = top(i)  
    l_mergep = merge_point(i)
 
    ! VARUN: TO CHECK
    zrain_rate = (prr_gsp(i)+prr_con(i))  ! [kg/m**2 s]
    zsnow_rate = (prs_gsp(i)+prs_con(i))
 
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
      CALL new_snow_density(rho_hn, t(i), zuv, qv(i))      

      ! Calculate new snow amounts
      hn_sn_now(i) =  hn_sn_now(i) + ( (zsnow_rate * zdt/rho_hn) )


      num_new_elems = floor(hn_sn_now(i)/new_snow_elem)

      if(num_new_elems >= 1) then

        hn_sn_now(i) = hn_sn_now(i) - num_new_elems * new_snow_elem

        if( num_new_elems + l_top > ke_snow) then
           merge_elems = num_new_elems + l_top - ke_snow
           do j = 1,merge_elems
            theta_i_now(i,l_mergep) =   (theta_i_now(i,l_mergep+j)*dzm_sn_now(i,l_mergep+j)                                        &
                               +  theta_i_now(i,l_mergep)*dzm_sn_now(i,l_mergep)) / (dzm_sn_now(i,l_mergep+j) + dzm_sn_now(i,l_mergep))    ! volumetric ice content

            theta_w_now(i,l_mergep) =   (theta_w_now(i,l_mergep+j)*dzm_sn_now(i,l_mergep+j)                                        &
                               +  theta_w_now(i,l_mergep)*dzm_sn_now(i,l_mergep)) / (dzm_sn_now(i,l_mergep+j) + dzm_sn_now(i,l_mergep))    ! volumetric ice content

            t_sn_elem(i,l_mergep)   =   (t_sn_elem(i,l_mergep+j)*dzm_sn_now(i,l_mergep+j)                                           &
                               +  t_sn_elem(i,l_mergep)*dzm_sn_now(i,l_mergep))    / (dzm_sn_now(i,l_mergep+j) + dzm_sn_now(i,l_mergep))    ! snow layer temperature

            dzm_sn_now(i,l_mergep)  = dzm_sn_now(i,l_mergep) + dzm_sn_now(i,l_mergep+j)                                             ! layer thickness

           enddo

           elems_left = (l_top - l_mergep - merge_elems)

           ! shift all elems down the array
           theta_i_now(i,l_mergep+1:l_mergep+elems_left) = theta_i_now(i,l_mergep+merge_elems+1:l_top)
           theta_w_now(i,l_mergep+1:l_mergep+elems_left) = theta_w_now(i,l_mergep+merge_elems+1:l_top)
           t_sn_elem  (i,l_mergep+1:l_mergep+elems_left) =   t_sn_elem(i,l_mergep+merge_elems+1:l_top)

           l_top = l_top - merge_elems;


        endif ! if( num_new_elems + l_top > ke_snow)

        do j=1,num_new_elems
   
          ! add new elems
          dzm_sn_now(i,l_top + j) = new_snow_elem
          rho_sn(i,l_top + j)     = rho_hn
          theta_i_now(i,l_top+j) = rho_hn / rho_i
          theta_w_now(i,l_top+j) = 0.0_vpp
          theta_a_now(i,l_top+j) = 1.0_vpp - theta_i_now(i,l_top+j) - theta_w_now(i,l_top+j)
          t_sn_elem(i,l_top+j)    = min(t(i),l_t0_melt)
 
        enddo
          l_top = l_top + num_new_elems

      endif ! if (num_new_elems >= 1)

      ! shift for changing merge point 

      if(dzm_sn_now(i,l_mergep) > 1.0_vpp) then
         l_mergep = l_mergep + 1
      endif

  endif

     top(i) = l_top      
     merge_point(i) = l_mergep
enddo

! ------------------------------------------------------------------------------
! - End Section I: New snow
! ------------------------------------------------------------------------------

!! updating element properties
   DO i = ivstart, ivend

     l_top = top(i)
     IF(l_top .GE. 1) THEN

         theta_a_now(i,:) = max(0.0_vpp,1.0_vpp - theta_i_now(i,:) - theta_w_now(i,:))

         ! --------------------------
         ! Snow layer density
         ! -------------------------
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
         do ksn = 1, l_top, 1
             if(rho_sn(i,ksn) .LT. eps) then
               hcap_sn(i,ksn) = 0.0_vpp
             else
               hcap_sn(i,ksn) = (  rho_a   * theta_a_now(i,ksn) * specific_heat_air     &
                                 + rho_i   * theta_i_now(i,ksn) * specific_heat_ice     &
                                 + rho_w   * theta_w_now(i,ksn) * specific_heat_water)  &
                                  / rho_sn(i,ksn)
             endif
         enddo

         ! --------------------------
         ! Heat conductivity
         ! --------------------------
         do ksn = 1, l_top, 1
             hcon_sn(i,ksn) = 2.22_vpp * EXP(1.88_vpp * LOG(rho_sn(i,ksn)/rho_i))
         enddo

     ENDIF ! if (l_top .ge. 1) 
   ENDDO ! loop over grid points

! =============================================================================
! + Begin Section II: Atmospheric forcing & Heat Equation 
! =============================================================================

  ! --------------------
  ! Some updating
  ! --------------------

    DO i = ivstart, ivend
      l_top = top(i)
      IF(l_top .GE. 1) THEN
        t_sn_sfc(i) = t_sn_now(i,l_top+1)
      ENDIF
    ENDDO


  ! ----------------------------------------------------------------------------
  ! + Surface fluxes
  ! ---------------------------------------------------------------------------
  
    DO i = ivstart, ivend

      IF(top(i) .GE. 1 ) THEN   ! snow on the ground

          zuv        = SQRT ( u(i)*u(i) + v(i)*v(i) )        ! wind speed
      
        ! Calculate transfer coefficient
        ! ---------------  
          CALL calc_tch(t_sn_sfc(i), t(i), tch_sn(i),  zuv, qv(i), h_snow(i), ps(i))

        ! Calculate turbulent fluxes
        ! ---------------
          CALL turb_flux(sh_sn(i), lh_sn(i), ps(i), t(i), t_sn_sfc(i), qv(i), tch_sn(i), zuv)

        ! Calculate radiative fluxes
        ! ---------------

        !  CALL rad_flux(swdir_s(i), swdifd_s(i), swdifu_s(i), lwd_s(i), lwu_s(i), alpha_sn(i), t(i), t_sn_sfc(i))

           emiss   = thbs(i)/(sigma*t(i)*t(i)*t(i)*t(i))
           t_emiss = sqrt(sqrt(emiss)) * t(i)

           delta_rad(i) = sigma * (t_emiss + t_sn_sfc(i)) * (t_emiss * t_emiss + t_sn_sfc(i) * t_sn_sfc(i))


      ENDIF
    l_lh_sn = lh_sn(i)
      !lh_sn(i) = 0.0_vpp
    ENDDO
  ! ----------------------------------------------------------------------------
  ! + Atmospheric forcing
  ! ----------------------------------------------------------------------------


          DO i = ivstart, ivend

            IF(top(i) .GE. 1) THEN ! snow on the ground

              ! Calculate net short wave radiation
              ! -----------------------------------

              swnet_sn(i) = sobs(i) !(swdir_s(i) + swdifd_s(i))  - swdifu_s(i)

              for_sn(i)  = sobs(i) + lh_sn(i) ! not taken into linearized boundary conditions

             ENDIF
           ENDDO

  ! ----------------------------------------------------------------------------
  ! + Solve the heat equation
  ! ----------------------------------------------------------------------------


     DO i = ivstart, ivend

       l_top = top(i)

       lower_bc = t_so_new(i,0)
       IF(l_top .GT. 1) THEN  !!!snow on the ground


       counter = 0
       DO ksn = l_top+1,1,-1
          counter = counter+1

          IF(counter .EQ. 1) THEN ! ... top node ! NEUMANN BC

            alpha_solver_up   = dzm_sn_now(i,ksn-1) / ( 6.0_vpp * dt )
            beta_solver_up    = ( hcon_sn(i,ksn-1) / ( rho_sn(i,ksn-1)*hcap_sn(i,ksn-1) ) ) * (1.0_vpp/dzm_sn_now(i,ksn-1))
        
            a_matrix(counter)     =  0.0_vpp
            b_matrix(counter)     =  2.0_vpp * alpha_solver_up + beta_solver_up + & 
&                                    tch_sn(i) + delta_rad(i)

            c_matrix(counter)     =  alpha_solver_up - beta_solver_up
 
            d_matrix(counter)     =  2.0_vpp * alpha_solver_up * t_sn_now(i,ksn) + alpha_solver_up * t_sn_now(i,ksn-1) + &
&                                    tch_sn(i) * t(i) + &
&                                    delta_rad(i) * t(i)

          ELSEIF (counter .EQ. l_top) THEN ! ... bottom node ! Dirichlet BC

            a_matrix(counter) = 0.0_vpp
            b_matrix(counter) = 1.0_vpp
            c_matrix(counter) = 0.0_vpp
 
            d_matrix(counter) = lower_bc  

          ELSE  ! Middle nodes

            alpha_solver_up  = dzm_sn_now(i,ksn) / ( 6.0_vpp * dt )
            beta_solver_up   = ( hcon_sn(i,ksn) / ( rho_sn(i,ksn)*hcap_sn(i,ksn) ) ) * (1.0_vpp/dzm_sn_now(i,ksn))

            alpha_solver_down = dzm_sn_now(i,ksn-1) / ( 6.0_vpp * dt )
            beta_solver_down  = ( hcon_sn(i,ksn-1) / ( rho_sn(i,ksn-1)*hcap_sn(i,ksn-1) ) ) * (1.0_vpp/dzm_sn_now(i,ksn-1))

            a_matrix(counter) = alpha_solver_up - beta_solver_up 
            b_matrix(counter) = 2.0_vpp * (alpha_solver_up + alpha_solver_down) + beta_solver_up + beta_solver_down 
            c_matrix(counter) = alpha_solver_down - beta_solver_down

            d_matrix(counter) = alpha_solver_up * t_sn_now(i,ksn+1) +  & 
&                                2.0_vpp * (alpha_solver_up + alpha_solver_down) * t_sn_now(i,ksn) + &
&                                alpha_solver_down * t_sn_now(i,ksn-1)

          ENDIF ! if block for splitting between top, middle and bottom nodes
       ENDDO

       ! ------------------------------------------------------------
       ! Solve the system - Thomas Algorithm
       ! ------------------------------------------------------------

       ! step 1: forward elimination
       do ksn=2,l_top
         coeff  = a_matrix(ksn)/b_matrix(ksn-1)
         b_matrix(ksn) = b_matrix(ksn) - coeff * c_matrix(ksn-1)
         d_matrix(ksn) = d_matrix(ksn) - coeff * d_matrix(ksn-1)
      enddo

      ! step 2: back substitution
      t_sn_new(i,l_top) = d_matrix(l_top)/b_matrix(l_top)
  
      do ksn=l_top-1,1,-1
         t_sn_new(i,ksn) = (d_matrix(ksn) - c_matrix(ksn) * t_sn_new(i,ksn))/b_matrix(ksn)
      enddo


       ENDIF ! if block for snow or no snow
     ENDDO

! =============================================================================
! + End Section II: Atmospheric forcing & Heat Equation
! =============================================================================


!! updating element properties
   DO i = ivstart, ivend

     l_top = top(i)
     IF(l_top .GE. 1) THEN

         theta_a_now(i,:) = max(0.0_vpp,1.0_vpp - theta_i_now(i,:) - theta_w_now(i,:))

         ! --------------------------
         ! Snow layer density
         ! -------------------------
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
         do ksn = 1, l_top, 1
             if(rho_sn(i,ksn) .LT. eps) then
               hcap_sn(i,ksn) = 0.0_vpp
             else
               hcap_sn(i,ksn) = (  rho_a   * theta_a_now(i,ksn) * specific_heat_air     &
                                 + rho_i   * theta_i_now(i,ksn) * specific_heat_ice     &
                                 + rho_w   * theta_w_now(i,ksn) * specific_heat_water)  &
                                  / rho_sn(i,ksn)
             endif
         enddo

         ! --------------------------
         ! Heat conductivity
         ! --------------------------
         do ksn = 1, l_top, 1
             hcon_sn(i,ksn) = 2.22_vpp * EXP(1.88_vpp * LOG(rho_sn(i,ksn)/rho_i))
         enddo

     ENDIF ! if (l_top .ge. 1) 
   ENDDO ! loop over grid points

!!! ENDING UPDATING ELEMENT PROPERTIES


! =============================================================================
! + Begin Section III: Phase Changes
! =============================================================================

    ! Keep old values
    ! -----------------

    DO ksn = 1, ke_snow, 1
       DO i = ivstart, ivend
         theta_i_old(i,ksn) = theta_i_now(i,ksn)
         theta_w_old(i,ksn) = theta_w_now(i,ksn)
       ENDDO
    ENDDO
    
    ! Initiate phase change
    ! -------------------

    DO i = ivstart, ivend
      l_top = top(i)
      IF(l_top .GE. 1) THEN ! snow on the ground
        ! Initiate
        melt_flag(i)   = .FALSE.
        freeze_flag(i) = .FALSE.
        DO ksn = l_top, 1, -1
           ! Metling
           IF(t_sn_elem(i,ksn) .GT. t0_melt .AND. theta_i_now(i,ksn) .GT. eps2) THEN
                melt_flag(i) = .TRUE.
           ENDIF
           ! Freezing
           IF(t_sn_elem(i,ksn) .LT. t0_melt .AND. theta_w_now(i,ksn) .GT. (theta_r + eps2)) THEN
                freeze_flag(i) = .TRUE.
           ENDIF
           !write(*,*) 'melt: ',ksn,melt_flag(i),freeze_flag(i)
        ENDDO  !end of snow layers
      ENDIF
    ENDDO


    ! Check for melting
    ! --------------------
    
    DO i = ivstart, ivend
       l_top = top(i)
       IF(l_top .GE. 1) THEN ! snow on the ground
         IF(melt_flag(i)) THEN ! melting

           q_rest = 0.0_vpp
           DO ksn = l_top, 1, -1

                ! --------------------
                ! Now see if any melting is going on -- this implies that (1) the
                ! temperature of the element
                ! is above/equal the melting temperature (2) there is something to melt and
                ! (3) there is enough
                ! room to place the meltwater ...
                ! --------------------
                !write(*,*) 'melting: ',ntstep,ksn,q_rest
               
                IF(t_sn_elem(i,ksn) .GE. t0_melt .AND. theta_i_now(i,ksn) .GT. 0.0_vpp .AND. theta_w_now(i,ksn) .LT. theta_s) THEN
            

 
                  ! difference dT between actual melting temperature and layer temperature
                  dT_sub_melt = 0.0_vpp
                  dT_sub_melt = t0_melt - t_sn_elem(i,ksn)
                  
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
                     dtheta_i   = A_sub_melt * dT_sub_melt
                     dtheta_w   = - (rho_i / rho_w) * dtheta_i
                    ! --------------------
                    ! It could happen that there is enough energy available to melt more ice
                    ! than is present.
                    ! You can only melt so much ice as is there ....
                    ! --------------------
                     IF( (theta_i_now(i,ksn) + dtheta_i) .LT. 0.0_vpp) THEN
                       dtheta_i    = - theta_i_now(i,ksn)
                       dtheta_w    = - (rho_i / rho_w) * dtheta_i
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
                     t_sn_elem(i,ksn) = t_sn_elem(i,ksn) + dT_sub_melt
                     IF(t_sn_elem(i,ksn) .LE. t0_melt) THEN ! if melting occured it can only be at melt point
                       q_rest     = 0.0_vpp
                       t_sn_elem(i,ksn)  = t0_melt
                     ELSE
                       q_rest = hcap_sn(i,ksn) * rho_sn(i,ksn) * dzm_sn_now(i,ksn) * (t_sn_elem(i,ksn) - t0_melt)
                       t_sn_elem(i,ksn)   = t0_melt
                     ENDIF
                     ! Volumetric freezing power
                     q_mf     = q_mf + ((dtheta_i * rho_i * lh_f) / (dt ) )
                     dtheta_w = dtheta_w
                     ! Contents of ice, water and air
                     theta_i_now(i,ksn) = theta_i_now(i,ksn) + dtheta_i
                     theta_w_now(i,ksn) = theta_w_now(i,ksn) + dtheta_w
                  ENDIF ! deltaT check
                  !write(*,*) 'melting: ',ntstep,ksn,dtheta_i,dtheta_w

                 ENDIF ! end melt check
           ENDDO ! loop over snow layers


                 do ksn = l_top,1,-1
                    theta_a_now(i,ksn) = max(0.0_vpp,1.0_vpp - theta_i_now(i,ksn) - theta_w_now(i,ksn))
                 enddo
         
           ! --------------------------
           ! Snow layer density
           ! -------------------------
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
                DO ksn = 1, l_top, 1
                   hcon_sn(i,ksn) = 2.22_vpp * EXP(1.88_vpp * LOG(rho_sn(i,ksn)/rho_i))
                ENDDO
         
        
          ENDIF ! end check in melt or not
       ENDIF ! end check if top .ge. 1
     ENDDO ! end loop over horizontal pixels


     ! Check for freezing
     ! --------------------

     DO i = ivstart, ivend
       l_top = top(i)
       IF(l_top .GE. 1) THEN ! snow on the ground
         IF(freeze_flag(i)) THEN ! freeze

           DO  ksn = l_top, 1, -1
              ! Freezing within the snowpack can occur if (1) the temperature of the
              ! element is below freezing and if water is present to be refrozen
              IF(t_sn_elem(i,ksn) .LT. t0_melt .AND. theta_w_now(i,ksn) .GT. theta_r) THEN
                ! difference dT between actual layer temperature and freezing temperature
                dT_sub_melt = 0.0_vpp
                dT_sub_melt = t0_melt - t_sn_elem(i,ksn)
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
                  theta_i_now(i,ksn)  = 1.0_vpp - theta_r
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
                q_mf     = q_mf + ((dtheta_i * rho_i * lh_f) /( dt ))
                dtheta_w = dtheta_w
                t_sn_elem(i,ksn)     = t_sn_elem(i,ksn) + dT_sub_melt
             ENDIF
           ENDDO
 

           ! --------------------------
           ! Depth of snow layer (main) levels
           ! -------------------------
                 do ksn = l_top,1,-1
                    theta_a_now(i,ksn) = max(0.0_vpp,1.0_vpp - theta_i_now(i,ksn) - theta_w_now(i,ksn))
                 enddo
         
           ! --------------------------
           ! Snow layer density
           ! -------------------------
                  do ksn = 1, l_top, 1
                     if(theta_i_now(i,ksn) .eq. 0.0_vpp) then
                        rho_sn(i,ksn) = 0.0_vpp
                     else
                        rho_sn(i,ksn) = theta_i_now(i,ksn)*rho_i + theta_w_now(i,ksn)*rho_w
                     endif
                  enddo
         
          ENDIF ! if freeze
       ENDIF ! if there is snow or not
     ENDDO ! end loop over horizontal pixels

! =============================================================================
! - End Section III: Phase Changes
! =============================================================================


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
  excess_water=0.0_vpp
  DO i = ivstart, ivend
    l_top = top(i)
    IF(l_top .GE. 1) THEN ! snow on the ground

      ! Initiate some values
       dL   = 0.0_vpp
       dM   = 0.0_vpp
       M    = 0.0_vpp
       hoar = 0.0_vpp
      ! Latent heat flux towards the surface - mass gain
      ! ---------------------
      IF(lh_sn(i) .GT. eps2) THEN ! add mass
        IF(t_sn_sfc(i) .LT. t0_melt) THEN  ! add ice
    
          dM = lh_sn(i) *( dt)/lh_s     ! calculate mass change
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
                      write(*,*) 'update_C: ',iblock,i,ksn,theta_i_now(i,l_top)
                  endif
        ELSE ! add water
          dM = lh_sn(i) *( dt) /lh_v      ! calculate mass change
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
               dM = lh_sn(i) *( dt)/lh_v
                M = theta_w_now(i,ksn) * rho_w * dzm_sn_now(i,ksn)
                 ! Check that you only take the available amount of water
                 IF(-dM .GE. M) THEN
                   dM = -M
                   theta_w_now(i,ksn) = theta_w_now(i,ksn) + dM/(rho_w*dzm_sn_now(i,ksn))
                 ELSE
                   theta_w_now(i,ksn) = theta_w_now(i,ksn) + dM/(rho_w*dzm_sn_now(i,ksn))
                 ENDIF
               lh_sn(i) = lh_sn(i) - dM*lh_v/( dt ) ! update energy used
             ELSEIF (theta_i_now(i,ksn) .GT. eps) THEN ! there is no water then sublimate ice matrix
               dM = lh_sn(i) * ( dt ) /lh_s
                M = theta_i_now(i,ksn) * rho_i * dzm_sn_now(i,ksn)
                IF(-dM .GT. M) THEN ! all ice can be sublimated
                  dM = -M
                  theta_i_now(i,ksn) = 0.0_vpp
                  dzm_sn_now(i,ksn)  = 0.0_vpp
                  if (ntstep .eq. N_IMP) then
                      write(*,*) 'update_B: ',iblock,ntstep,i,ksn,theta_i_now(i,ksn)
                  endif

                ELSE
                  dzm_sn_old = dzm_sn_now(i,ksn)
                  dL = dM/rho_sn(i,ksn)
                  dzm_sn_now(i,ksn) = dzm_sn_now(i,ksn) + dL
                  theta_i_now(i,ksn) = theta_i_now(i,ksn) * (dzm_sn_old/dzm_sn_now(i,ksn))
                  theta_i_now(i,ksn) = theta_i_now(i,ksn) + (dM/(rho_i*dzm_sn_now(i,ksn)))
                  theta_w_now(i,ksn) = theta_w_now(i,ksn) * (dzm_sn_old/dzm_sn_now(i,ksn))
                  if (ntstep .eq. N_IMP) then
                      write(*,*) 'update_A: ',iblock,i,ksn,theta_i_now(i,ksn)
                  endif
                ENDIF
               lh_sn(i) = lh_sn(i) - dM*lh_v/( dt ) ! update energy used
             ENDIF
           ENDDO ! end of ksn
   
       ! MASSIVE HACK HERE
       IF(lh_sn(i) .LT. (-1.0_vpp*eps2)) THEN ! there is still energy left, which should technically be used by the soil layer for now let's erase it
         lh_sn(i) = 0.0_vpp
       ENDIF
     ENDDO ! end of while
   ENDIF
  ENDIF

         
           ! --------------------------
           ! Depth of snow layer (main) levels
           ! -------------------------
                 do ksn = l_top,1,-1
                    theta_a_now(i,ksn) = max(0.0_vpp,1.0_vpp - theta_i_now(i,ksn) - theta_w_now(i,ksn))
                 enddo
         
           ! --------------------------
           ! Snow layer density
           ! -------------------------
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
                 DO ksn = 1, l_top, 1
                    hcon_sn(i,ksn) = 2.22_vpp * EXP(1.88_vpp * LOG(rho_sn(i,ksn)/rho_i))
                 ENDDO
         

      ! -------------------
      ! Transport the water
      ! -------------------

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
         dtheta_w_sub_wtr(ksn) = (hcap_sn(i,ksn) * rho_sn(i,ksn)) * (1.0 / lh_f ) * (1.0 / rho_w)  * MAX(0.0_vpp, (t0_melt - t_sn_elem(i,ksn)))
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
          IF(theta_i_now(i,ksn) .LT. 0.23_vpp) THEN
            res_wat_cont(ksn) = 0.0264_vpp + 0.0099_vpp * (1.0_vpp - theta_i_now(i,ksn)) / theta_i_now(i,ksn)
          ELSE
            res_wat_cont(ksn) = 0.08_vpp - 0.1023_vpp * (theta_i_now(i,ksn) - 0.03_vpp)
            !res_wat_cont(ksn) = 0.01_vpp - 0.1023_vpp * (theta_i_now(i,ksn) - 0.03_vpp)
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
      runoff_sn(i) = 0.0_vpp
      theta_w_bot = theta_w_now(i,1)
        IF (theta_w_bot .GT. w_res(1)) THEN
         ! Adjust dependent values accordingly
         theta_w_now(i,1) = w_res(1)
         theta_a_now(i,1) = 1.0_vpp - theta_w_now(i,1) - theta_i_now(i,1)
         ! Put all excess water of bottom layer  in runoff, i.e. move it out of the snow cover
         ! Note: if one comments this out you get ponding
         runoff_sn(i)              = runoff_sn(i) + excess_water +  dzm_sn_now(i,1) * (theta_w_bot - w_res(1))
        ENDIF
   
           ! --------------------------
           ! Depth of snow layer (main) levels
           ! -------------------------
                 do ksn = l_top,1,-1
                    theta_a_now(i,ksn) = max(0.0_vpp,1.0_vpp - theta_i_now(i,ksn) - theta_w_now(i,ksn))
                 enddo
         
           ! --------------------------
           ! Snow layer density
           ! -------------------------
                  do ksn = 1, l_top, 1
                     if(theta_i_now(i,ksn) .eq. 0.0_vpp) then
                        rho_sn(i,ksn) = 0.0_vpp
                     else
                        rho_sn(i,ksn) = theta_i_now(i,ksn)*rho_i + theta_w_now(i,ksn)*rho_w
                     endif
                  enddo
                  
    ENDIF


  ENDDO ! loop over horizontal pixels i=ivstart,ivend

! =============================================================================
! - End Section IV: Water transport
! =============================================================================

! =============================================================================
! + Begin Section V: Settling
! =============================================================================

  ! Calculate settling rates
  ! ---------------------------------

  DO i = ivstart, ivend
    l_top = top(i)

    IF(l_top .GE. 1) THEN ! snow on the ground

     ! Initalizations
     ddz          = 0.0_vpp   ! total change in layer thickness
     overburden   = 0.0_vpp   ! overburden stress (weight of the overlying layers times g)

     tmp_rate_1 = 0.0_vpp
     tmp_rate_2 = 0.0_vpp

     DO ksn = l_top, 1, -1

       ! Reset values for settling rates
       rate_1    = 0.0_vpp      ! due to ice loss
       rate_2    = 0.0_vpp      ! due to overburden stress

       tot_rate  = 0.0_vpp      ! total settliing rate

       dz_old    = 0.0_vpp      ! old layer thickness (before settling) ...
       dz_new    = 0.0_vpp      ! ... and after settling

       !  ... due to loss of ice
       !-------------------
         IF(theta_i_now(i,ksn) <  theta_i_old(i,ksn)) THEN ! melting ocured in that layers
           rate_1 = (-0.5_vpp/dt) * MAX(0.0_vpp,  (theta_i_old(i,ksn) - theta_i_now(i,ksn)) / theta_i_old(i,ksn) )
           !rate_1 = -0.5_vpp/dt * MAX(0.0_vpp,  (theta_i_old(i,ksn) - theta_i_now(i,ksn)) / theta_i_old(i,ksn) )
         ENDIF
       !  ... due to overburden stress
       !-------------------
       ! Vionnett (2012)
       f1   = 1.0_vpp / (1.0_vpp + 60.0_vpp * ((theta_w_now(i,ksn)*rho_w*dzm_sn_now(i,ksn)) / (rho_w * dzm_sn_now(i,ksn))))
       f2   = 1.0_vpp

       eta = f1 * f2 * eta_0 * (rho_sn(i,ksn)/c_eta) * exp(a_eta*(t0_melt - t_sn_elem(i,ksn)) + b_eta*rho_sn(i,ksn))
       rate_2 = -1.0_vpp * (overburden + ( (dzm_sn_now(i,ksn)*rho_sn(i,ksn))*g/2.0_vpp)) / eta


       ! increase overburden stress NOTE: How to deal with slope angle should be overburden = m*g*cos(alpha)
       overburden = overburden + ( (dzm_sn_now(i,ksn)*rho_sn(i,ksn)) * g)

      ! ... calculate change ...
       !-------------------

       ! ... of all (ice loss, overburden, destructive) settling rates (1/s)
       tot_rate = (rate_1*dt) + (rate_2*dt)


       tmp_rate_1 = tmp_rate_1 + rate_1
       tmp_rate_2 = tmp_rate_2 + rate_2

       ! ... of layer thickness, i.e. sum over all layers (m)
       ddz          = ddz + MAX(-1.0_vpp * dzm_sn_now(i,ksn), dzm_sn_now(i,ksn) * tot_rate)


       dz_old = dzm_sn_now(i,ksn)
       dz_new = dz_old + MAX(-1.0_vpp * dzm_sn_now(i,ksn), dzm_sn_now(i,ksn) * tot_rate)

       ! ... volumetric contents
       theta_i_now(i,ksn) = MAX(0.0_vpp, theta_i_now(i,ksn) * (dz_old / dz_new))    ! ice content
       theta_w_now(i,ksn) = MAX(0.0_vpp, theta_w_now(i,ksn) * (dz_old / dz_new))    ! water content

       ! ... of layer thickness (m)
       dzm_sn_now(i,ksn) = dz_new

  ENDDO

    ENDIF

  ENDDO

! =============================================================================
! - End Section V: Settling
! =============================================================================

  do i = ivstart, ivend

    l_top = top(i)
    l_mergep = merge_point(i)
    l_thinlayer_merge = 0

    if(l_top .GE. 2) THEN ! snow on the ground
       do ksn = l_top, 2, -1
          if(dzm_sn_now(i,ksn) .LT. min_height_layer .OR. theta_i_now(i,ksn) .LT. 0.01_vpp) THEN ! layer is quite thin

            l_thinlayer_merge = l_thinlayer_merge + 1

            ! -------------------
            ! Aggregate with lower layer - Adjust properties
            ! -------------------
            theta_i_now(i,ksn-1) = (theta_i_now(i,ksn)*dzm_sn_now(i,ksn)                                      &
                             + theta_i_now(i,ksn-1)*dzm_sn_now(i,ksn-1)) / (dzm_sn_now(i,ksn) + dzm_sn_now(i,ksn-1))  ! volumetric ice content

            theta_w_now(i,ksn-1) = (theta_w_now(i,ksn)*dzm_sn_now(i,ksn)                                      &
                             + theta_w_now(i,ksn-1)*dzm_sn_now(i,ksn-1)) / (dzm_sn_now(i,ksn) + dzm_sn_now(i,ksn-1))  ! volumetric water content

            t_sn_elem(i,ksn-1) = (t_sn_elem(i,ksn)*dzm_sn_now(i,ksn)                                      &
                                + t_sn_elem(i,ksn-1)*dzm_sn_now(i,ksn-1)) / (dzm_sn_now(i,ksn) + dzm_sn_now(i,ksn-1))  ! volumetric ice content
           
            dzm_sn_now(i,ksn-1)  = dzm_sn_now(i,ksn) + dzm_sn_now(i,ksn-1)        ! layer thickness

            ! ... and reset roperties
            dzm_sn_now(i,ksn)  = 0.0_vpp
            theta_i_now(i,ksn) = 0.0_vpp
            theta_w_now(i,ksn) = 0.0_vpp
            theta_a_now(i,ksn) = 0.0_vpp
            t_sn_elem(i,ksn)    = 0.0_vpp

          endif ! if element to be merged or not

          j=0
          do lay=1,l_top
             if((theta_i_now(i,lay) .gt. 0.01_vpp) ) then
    
                j=j+1
                dzm_sn_now(i,j) = dzm_sn_now(i,lay)
                theta_i_now(i,j) = theta_i_now(i,lay)
                theta_w_now(i,j) = theta_w_now(i,lay)
                theta_a_now(i,j) = theta_a_now(i,lay)
                t_sn_elem(i,j) = t_sn_elem(i,lay) 
 
             endif
          enddo
          l_top = j

          do lay = l_top+1, ke_snow, 1
             dzm_sn_now(i,lay)  = 0.0_vpp
             theta_i_now(i,lay) = 0.0_vpp
             theta_w_now(i,lay) = 0.0_vpp
             theta_a_now(i,lay) = 0.0_vpp
             t_sn_elem(i,lay)    = 0.0_vpp
          enddo

       enddo ! end of loop over elements for merging of thin layers


       ! split merge-point layer: remove one cm from this layer

       if ( l_thinlayer_merge .gt. 0 ) then

          if( (l_mergep .eq. 1) .and. (dzm_sn_now(i,1) < 0.01) ) then
          
          
          else

            ! create space for new layers from merge-point layer split
             dzm_sn_now(l_mergep+1+l_thinlayer_merge:l_top+l_thinlayer_merge,i) =  dzm_sn_now( l_mergep+1:l_top , i )
            theta_i_now(l_mergep+1+l_thinlayer_merge:l_top+l_thinlayer_merge,i) = theta_i_now( l_mergep+1:l_top , i )
            theta_w_now(l_mergep+1+l_thinlayer_merge:l_top+l_thinlayer_merge,i) = theta_w_now( l_mergep+1:l_top , i )
            theta_a_now(l_mergep+1+l_thinlayer_merge:l_top+l_thinlayer_merge,i) = theta_a_now( l_mergep+1:l_top , i )
              t_sn_elem(l_mergep+1+l_thinlayer_merge:l_top+l_thinlayer_merge,i) =   t_sn_elem( l_mergep+1:l_top , i )

            do j=1,l_thinlayer_merge

                 dzm_sn_now(l_mergep+j,i) = new_snow_elem
                theta_i_now(l_mergep+j,i) = theta_i_now(l_mergep,i)
                theta_a_now(l_mergep+j,i) = theta_a_now(l_mergep,i)
                theta_w_now(l_mergep+j,i) = theta_w_now(l_mergep,i)
                  t_sn_elem(l_mergep+j,i) =   t_sn_elem(l_mergep,i) 

            enddo

                  dzm_sn_now(l_mergep,i) = dzm_sn_now(l_mergep,i) - (new_snow_elem) * l_thinlayer_merge
                  l_top = l_top + l_thinlayer_merge
         
          endif

       endif

       if ( (dzm_sn_now(l_mergep,i) < 0.01) .and. (l_mergep > 1) ) then
            l_mergep = l_mergep - 1
       endif

   endif
 enddo

! =============================================================================
! + Begin Section VI: Preparations for the next time step 
! =============================================================================
  DO i = ivstart, ivend
      l_top = top(i)
      if ( l_top .ge. 1) then
           ! update
         
           ! --------------------------
           ! Depth of snow layer (main) levels
           ! -------------------------
                 do ksn = l_top,1,-1
                    theta_a_now(i,ksn) = max(0.0_vpp,1.0_vpp - theta_i_now(i,ksn) - theta_w_now(i,ksn))
                 enddo
         
           ! --------------------------
           ! Snow layer density
           ! -------------------------
                  do ksn = 1, l_top, 1
                     if(theta_i_now(i,ksn) .eq. 0.0_vpp) then
                        rho_sn(i,ksn) = 0.0_vpp
                     else
                        rho_sn(i,ksn) = theta_i_now(i,ksn)*rho_i + theta_w_now(i,ksn)*rho_w
                     endif
                  enddo
         
      ENDIF

      top(i) = l_top
    ENDDO

! =============================================================================
! + Begin Section VII: Updating of prognostic variables
! =============================================================================

    ! Update snow layer properties
    ! -------------------------------------

    DO ksn = 1, ke_snow, 1
      DO i = ivstart, ivend

        t_sn_new(i,ksn)    = t_sn_now(i,ksn)     ! snow layer temperature  

        theta_i_new(i,ksn) = theta_i_now(i,ksn)  ! volumetric ice content
        theta_w_new(i,ksn) = theta_w_now(i,ksn)  ! volumetric water content
        theta_a_new(i,ksn) = theta_a_now(i,ksn)  ! volumetric air content

        dzm_sn_new(i,ksn)  = dzm_sn_now(i,ksn)   ! layer thickness

      END DO
    END DO


    ! Calculate snow water equivialent
    ! -------------------------------------
    DO ksn = 1, ke_snow, 1
      DO i = ivstart, ivend
       swe_sn(i) = swe_sn(i) + ( dzm_sn_now(i,ksn) * rho_sn(i,ksn) ) 
      END DO
    ENDDO

    ! Update snow height and surface temperature
    ! -----------------------------

    DO i = ivstart, ivend
      top_sn_new (i) = REAL(top(i))    ! index of first (top) snow layer
      hn_sn_new  (i) = hn_sn_now(i)    ! new snow storage

      top_sn_now (i) = REAL(top(i))    ! index of first (top) snow layer

      IF(top(i) .GE. 1) THEN

         h_snow(i) = 0.0_vpp
         do ksn = 1, ke_snow,1
            h_snow(i) = h_snow(i) + dzm_sn_now(i,ksn)
         enddo
      ELSE
         h_snow(i) = 0.0_vpp
      ENDIF

      if (h_snow(i) .gt. 0.0) then
        rho_snow_new(i) = swe_sn(i)/h_snow(i)
        t_snow_new (i) = t_sn_now(i,top(i)) 
        w_snow_new (i) = swe_sn(i)/1000.0_vpp
      endif
    END DO


    ! Update turbulent fluxes
    ! -----------------------------

    DO i = ivstart, ivend

      zshfl_snow(i) = sh_sn(i)
      zlhfl_snow(i) = lh_sn(i)
 !     zqhfl_snow(i) = lh_sn(i)/lh_s
    
    ENDDO


    if( mod(ntstep,120) .eq. 0) then 
       open(22, file = './out/pro_snow.txt')
       write(22,'(A7,   2I12.5,a,/)',    ADVANCE = 'YES') 'idx'         , ntstep, top(1)
       write(22,'(A7,  1000F12.3,a,/)',   ADVANCE = 'YES') 'dz'         , dzm_sn_now(1,:)
       write(22,'(A7,  1000F12.3,a,/)',   ADVANCE = 'YES') 't_sn'       , t_sn_now(1,:)
       write(22,'(A,   1000F12.3,a,/)',   ADVANCE = 'YES') 'theta_i'    , theta_i_now(1,:)
       write(22,'(A,   1000F12.3,a,/)',   ADVANCE = 'YES') 'theta_w'    , theta_w_now(1,:)
       write(22,'(A,   1000F12.3,a,/)',   ADVANCE = 'YES') 'theta_a'    , theta_a_now(1,:)
       write(22,'(A7,  1000F12.3,a,/)',   ADVANCE = 'YES') 'rho_sn'     , rho_sn(1,:)
       write(22,'(A,   1000F12.3,a,/)',   ADVANCE = 'YES') 'hcap_sn'    , hcap_sn(1,:)
       write(22,'(A,   1000F12.3,a,/)',   ADVANCE = 'YES') 'hcon_sn'    , hcon_sn(1,:)
       write(22,'(A,   1000F12.3,a,/)',   ADVANCE = 'YES') 'hdif_sn'    , hdif_sn(1,:)
    endif

    write(*,*) 'hsnow: ',ntstep,h_snow(1)

    if( mod(ntstep,30) .eq. 0) then
       open(44, file = './out/met.txt')  
       do i = ivstart,ivend
           l_top = top(i)
           write(44,'(1I6,1I6,1I6,15F18.3,a,/)', ADVANCE = 'YES') ntstep,i,l_top,ddz,t_sn_sfc(i),qv(i),zuv,& 
       &                           zrain_rate,sobs(i),thbs(i),l_lh_sn,sh_sn(i),h_snow(i),t_so_now(i,0),t_so_now(i,ke_soil+1),dzm_sn_new(i,l_top)
       enddo
    !   !close(44)
    endif

    if(ntstep .gt. N_IMP) then
      stop
    endif

! ------------------------------------------------------------------------------
! - End Section '?': Update prognostic variables
! ------------------------------------------------------------------------------



END SUBROUTINE snowpolino

! =============================================================================
!  - End subroutine  - snowpolino
! =============================================================================


! DONE, DONE!!!
END MODULE sfc_snow

!==============================================================================
! + End of module src_soil_multlay
!==============================================================================

