!+ Data module for all global meteorological fields
!------------------------------------------------------------------------------

MODULE data_fields

!------------------------------------------------------------------------------
!
! Description for TSA version
!  This is data_block_fields from COSMO version 5.07 without modifications
!
!------------------------------------------------------------------------------
!
! Description:
!  This module declares all meteorological fields that have to reside in 
!  the long term storage, i.e. that are used in more than one module.
!  Fields included are
!    - constant fields defining the reference atmosphere
!    - external parameter fields
!    - prognostic variables
!    - tendency fields for the prognostic variables
!    - fields for surface values
!    - fields that are computed in the parametrization packages 
!      or in the dynamics
!    - fields for model-output and diagnostics
!    - fields for the boundary values
!
!  All fields are declared as allocatable arrays. They are allocated in the
!  setup of the model and deallocated in the cleanup at the end of the
!  program.
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Ulrich Schaettler
!  Initial release
! 1.3        1998/04/15 Guenther Doms
!  Definition of new tendency arrays for convection
! 1.7        1998/07/16 Guenther Doms
!  Removal of global array 'rrssk'.
! 1.20       1999/01/07 Guenhter Doms
!  Renaming of some global variables.
! 1.24       1999/03/01 Guenther Doms
!  Declaration of a new 3-D prognostic variable (qi).
! 1.30       1999/06/24 Matthias Raschendorfer
!  Declaration of a 3-D prognostic var. (tke) and its tendenci field (tketens).
!  Declaration of a 3-D array (rcld).
!  Declaration of 5 3-D arrays for canopy layers (c_big,c_sml, r_air, t_e,qv_e).
!  Declaration of a 2-D arrays (tfh,tfm) and (h_can,d_pat).
!  Declaration of 4 2-D variables for new convection closures
! 1.33       1999/10/14 Reinhold Hess
!  Declaration of new global 2-D arrays idiv_hum and aevap_s for diagnosis
!  of the model water budget
!  Declaration of a new 2-D array 'sai' for surface area index (M.Raschendorfer)
! 1.34       1999/12/10 Ulrich Schaettler
!  Named all boundary fields with "_bd" (for consistency)
! 1.39       2000/05/03 Ulrich Schaettler
!  Add declaration of a boundary field for w (used for interactive nesting)
! 2.2        2000/08/18 Matthias Raschendorfer
!  Declaration of the 2-D arrays 'eai' and 'tai'.
! 2.4        2001/01/29 Christoph Schraff
!  Declaration of the 2-D array 'prne_con', for humidity balancing at T-nudging.
! 2.8        2001/07/06 Ulrich Schaettler
!  Added new fields for multi-layer soil model and surface fluxes
! 2.11       2001/09/28 Ulrich Schaettler
!  Added new fields for lateral values of cloud ice
! 2.17       2002/05/08 Ulrich Schaettler
!  New fields for Kain-Fritsch convection scheme
! 2.18       2002/07/16 Ulrich Schaettler
!  New fields for specific rain and snow content;
!  included declaration of a1t, a2t from src_leapfrog
! 3.5        2003/09/02 Ulrich Schaettler
!  New fields phi_tot, rla_tot to avoid global communication in the radiation
! 3.6        2003/12/11 Reinhold Schrodin
!  New field freshsnow for new multi-layer soil model added
! 3.7        2004/02/18 Ulrich Schaettler
!  New fields for computing synthetic satellite images (synme5-7, synmsg)
!  New field for storing convective cloud water (clw_con)
!  Renamed alb (alb_rad), idiv_hum (tdiv_hum), phi (rlat), rla (rlon),
!      cphi (crlat), acphir (acrlat), tgphi (tgrlat) (for consistency with GME
! 3.13       2004/12/03 Ulrich Schaettler
!  New fields for graupel scheme (qg, prg_gsp, grau_gsp): Thorsten Reinhardt
!  New fields for 3D turbulence  (tkhm, tkhh): Jochen Foerstner
!  New fields for fast_waves_rk: 1/SQRT(G): sqrtg_r_(s,u,v,w): Jochen Foerstner
!  New fields for external parameters (not yet used: for_e, for_d) R. Schrodin
!  Renamed w_ice to w_so_ice (to be consistent with GME)
! 3.17       2005/12/12 Reinhold Schrodin
!  New fields rho_snow (prognostic snow density) and h_snow (snow height) added
! 3.18       2006/03/03 Ulrich Schaettler
!  New fields for the CLM version              (by CLM Community)
!  New fields for the lake model FLake         (Dmitrii Mironov)
!  New field rh_2m for relative humidity in 2m (Matthias Raschendorfer)
!  New field tinc_lh to gather temperature increments due to latent heat 
!   conversion in microphysics (mainly saturation adjustment; Jochen Foerstner)
! 3.19       2006/04/25 Ulrich Schaettler
!  New field T_S_LAKE to save lake values from Nudging cycle
! 3.21       2006/12/04 Jochen Foerstner, Burkhardt Rockel, Christoph Schraff
!  New variables (LMK): t0, dt0dz, hd_mask_dcoeff, qvt_diff
!  Put declaration of zwcon here, because it is used in several modules
!  New boundary variables: qr_bd, qs_bd, qg_bd
!  Renamed variable sunshhrs, sodwdir to dursun, sodwddm
!  Additional variable introduced for deep atmosphere (Ronny Petrik)
!  Additional fields for Bechtold convection scheme (MeteoSwiss)
!  New fields for time-integrated analysis increments from nudging (section 10)
! V3_23        2007/03/30 Matthias Raschendorfer, Matteo Buzzi, Jochen Foerstner
!  Added field 'tfv' containing the laminar reduction factor for evaporation
!  Added fields for computations of topographic radiation correction (Matteo Buzzi)
!  Added fields for new relaxation (Jochen Foerstner)
! V4_3         2008/02/25 Matthias Raschendorfer
!  Introduction of a 3D diagnostic field 'edr' for the eddy dissipotion rate
! V4_4         2008/07/16 Jan-Peter Schulz
!  Added external parameters (sso_stdh, sso_gamma, sso_theta, sso_sigma) and
!  fields (ut_sso, vt_sso, tt_sso, ustr_sso, vstr_sso, vdis_sso,
!  austr_sso, avstr_sso, avdis_sso) for SSO scheme
! V4_5         2008/09/10 Ulrich Schaettler
!  Moved field for damping coefficients from src_relaxation to data_fields 
!  for global use
! V4_8         2009/02/16 Ulrich Schaettler
!  New fields for convective and dynamical gust in 10m
!  New fields for additional output of radiation values (for CLM)
!  Add global field for reference pressure at half levels (p0hl) (G. Zaengl)
! V4_9         2009/07/16 Ulrich Schaettler
!  New fields for diffusion masks for t-, q- and u-fields
! V4_10        2009/09/11 Ulrich Schaettler
!  New field for snow melt
!  New field for sea-mask (Jan-Peter Schulz for sea-ice model)
! V4_11        2009/11/30 Ekaterina Machulskaya, Juergen Helmert, Lucio Torrisi
!  Additional fields for the snow model (EM)
!  Additional external parameters for aerosol distributions, surface emissivity
!   and stomata resistance;  (JH)
!  Additional fields in output for SMA (LT)
! V4_12        2010/05/11 Michael Baldauf, Ulrich Schaettler
!  New fields dzeta_dlam, dzeta_dphi
!  Renamed hd_mask_dcoeff to hd_mask_dcoeff_p for diffusion of pressure
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_15        2010/11/19 Ulrich Schaettler
!  The field t_s_lake is removed again after adaptations in the SST Analysis
! V4_20        2011/08/31 Ulrich Schaettler
!  tgrlat needs 2 dimensions because of v-point dependence 
!    (reported by Andreas Will)
!  Introduction of additional 3D-arrays 'tket_(conv, sso, hshr)' for additional
!   TKE source by the action of other sub grid scale flow patterns.
!    (by Matthias Raschendorfer)
! V4_22        2012/01/31 Thorsten Reinhardt
!  Additional fields used to update the solar zenith angle in the radiation.
! V4_23        2012/05/10 Michael Baldauf, Oliver Fuhrer, Juergen Helmert, CLM
!  Eliminate fields sqrtg_r_*, dzeta_*, which are now in grid_metrics_utilities
!  Eliminate field qvt_diff, which has been eliminated from the model at all (OF)
!  Eliminated fields for Kain-Fritsch convection
!  Introduction of prescribed surface albedo: new 2D fields alb_dry, alb_sat (CLM)
!     and alb_dif (JH)
!  Introduction of new diagnostic variable for maximum wind speed in 10m height
!  Implementation of time dependent boundary values for aerosol optical depths
!   (Introduction of fields aer_[du,su,ss,bc,or]_bd) (CLM)
! V4_24        2012/06/22 Michael Baldauf
!  Added new field for reference temperature on half levels: t0hl(:,:,:)
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer, Ulrich Blahak
!                         Ekaterina Machulskaya
!  Eliminated fields for the humidity tracers qx 
!    (are now treated with new tracer module)
!  Added fields for the 2-moment scheme: precipitation rate and amount of hail (UB)
!  Added a new field tg_radstep: ground temperature at the last call to radiation (EM)
! V4_26        2012/12/06 Ulrich Schaettler
!  Adapted variable names of multi-layer snow model to corresponding
!   short names for I/O
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler 
!  Declare allocatable fields as pointer for MESSY
!  Global SAVE attribute for all constants
! V4_29        2013/10/04 Ulrich Schaettler, Astrid Kerkweg
!  Added 3 global fields: fr_snow, fr_wi, ustar_fv (which were local to
!   the soil model before): for a better handling in MESSY
! V5_1         2014-11-28 Ulrich Blahak, Oliver Fuhrer
!  Added advective tendency of TKE for itype_turb=3.
!  Added some more fields for 2-moment scheme
!  Replaced ireals by wp (working precision) (OF)
! V5_2         2015-05-21 Ulrich Schaettler
!  Include field cu_evap from convection
! V5_3         2015-10-09 Ulrich Blahak, Ulrich Schaettler
!  Added new fields:
!     ttens_diab for pure diabatic temperature tendency output and
!     ttdiab_conv  for pure diabatic temperature tendency due to convection
!     u_m, v_m: wind speeds on mass points (needed for COSMO-ICON physics)
!  Renamed alb_rad_coarse to alb_rad_rc; added fields for radiation
!  (but still commented, because they are not activated yet)
! V5_4a        2016-05-10 Matthias Raschendorfer
!  Introduced global fields for 
!    turbulent transfer velocities:                   tvm, tvh, tkr, tkred_sfc
!    horizontal wind tendencies from the turbulence scheme:  ut_turb, vt_turb
!    effective length scale of circulation patterns:         l_pat
! V5_4b        2016-07-12 Ulrich Blahak, Ulrich Schaettler, 
!                         Xavier Lapillonne, Jochen Foerstner
!  Added h0noise for idealized surface fluxes for the new blocked transfer 
!   scheme code (UB)
!  Added new fields for COSMO-ICON version of TERRA: w_p, w_s (US)
!  Added integer field isoiltyp, needed in new TERRA version (US)
!  Added field dqvdt_conv for 3d convective moisture convergence (XL)
!  Added fields w_conv, qhfl needed for convection "prepare" task (US)
!  Added fields llakemask, ttens_conv needed by Tiedtke-Bechtold convection (JF)
! V5_4c        2016-10-06 KIT, Ulrich Schaettler
!  Added additional fields for TWOMOM_SB
! V5_4e        2017-03-23 KIT, Ulrich Schaettler
!  Even more fields for TWOMOM_SB
!  Renamed variable qvsflx to qvfl_s (in accordance with ICON)
!  Eliminated h_can, c_big, c_sml, r_air
! V5_4f        2017-09-01 Ulrich Blahak, Ulrich Schaettler
!  Added comp_dbzobs and comp_dbzsim (dBZ-Composits on model grid) for
!   radar forward OPERATOR (UB)
!  Added track output variables for conv. cells:
!   uhel_track, rot_track, vih_track, vii_track, wmax_track, lpi_track, dbz_track (UB)
!  Removed fields for old 2-level soil model (US)
! V5_5         2018-02-23 Ulrich Blahak
!  Added echotop_z and echotop_p (UB)
!  Added mconv_hor_track (UB)
!  Declared comp_dbz[obs,sim] as threedimensional variables now
! V5_6a        2019-05-21 Jan-Peter Schulz
!  Introduce skin temperature approach
! V5_6b        2019-10-16 Pavel Khain, Harel Muskatel, Ulrich Blahak, Aleksei Poliukhov
!  Added new fields for CLOUDRAD: aot, ssa, asy (and _bd); qc_rad_con, qc_rad_adiab, weff_rad
!  Replaced TWOMOM_SB and COSMOART by CLOUDRAD for definition of some fields
! V5_7_SNOW    2021-01-08 Varun Sharma
!  Added variables (t_sn,theta_i,theta_w,theta_a,dzm_sn,hn_sn,top_sn) for SNOWPOLINO
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:
USE kind_parameters, ONLY :   &
    wp           ! KIND-type parameter for real variables

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Global (i.e. public) Declarations:

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------
#ifndef MESSY
  REAL  (KIND=wp),     TARGET, ALLOCATABLE ::           &
#else
  REAL  (KIND=wp),     POINTER             ::           &
#endif
    rho0 (:,:,:),     & ! reference density at the full model levels  (kg/m3)
    dp0  (:,:,:),     & ! reference pressure thickness of layers      ( Pa  )
    p0   (:,:,:),     & ! reference pressure at full levels           ( Pa  )
    p0hl (:,:,:),     & ! reference pressure at half levels           ( Pa  )
    t0hl (:,:,:),     & ! reference temperature at half levels        ( K   )
    dt0dz(:,:,:),     & ! temperature grad. of reference atmosphere   ( K/m )
    t0   (:,:,:),     & ! reference temperature                       ( K   )
    hhl  (:,:,:)        ! geometrical height of half levels           ( m   )

! 2. external parameter fields                                        (unit)
! ----------------------------
#ifndef MESSY
  REAL  (KIND=wp),     TARGET, ALLOCATABLE ::           &
#else
  REAL  (KIND=wp),     POINTER             ::           &
#endif
    hsurf  (:,:),   & ! height of surface topography                  ( m   )
    sso_stdh (:,:), & ! standard deviation of sub-grid scale orography( m   )
    sso_gamma(:,:), & ! anisotropy of sub-grid scale orography          --
    sso_theta(:,:), & ! angle betw. principal axis of orography and E ( rad )
    sso_sigma(:,:), & ! mean slope of sub-grid scale orography          --
    aer_su (:,:),   & ! monthly aerosol climatology sulfate drops     (0 - 1)
    aer_du (:,:),   & ! monthly aerosol climatology total dust        (0 - 1)
    aer_or (:,:),   & ! monthly aerosol climatology organic (water sol.)(0 - 1)
    aer_bc (:,:),   & ! monthly aerosol climatology black carbon      (0 - 1)
    aer_ss (:,:),   & ! monthly aerosol climatology sea salt          (0 - 1)
#ifdef CLOUDRAD
    aot (:,:,:),    & !new Aleksei Poliukhov 2017
    ssa (:,:,:),    & !
    asy (:,:,:),    & !
#endif
    emis_rad (:,:), & ! external thermal emissivity map               (0 - 1)
    rsmin2d(:,:)  , & ! minimum stomata resistance                    ( s/m )
    swi    (:,:,:), & ! soil wetness index                            (0 - 1)
    gz0    (:,:),   & ! surface roughness * g                         (m2/s2)
    fr_land(:,:),   & ! fraction of land in a grid element              --
    soiltyp(:,:),   & ! type of the soil (keys 0-9)                     --
    vio3   (:,:),   & ! vertical integrated ozone contents            (Pa O3)
    hmo3   (:,:),   & ! ozone maximum                                 ( Pa  )
    rlat   (:,:),   & ! geographical latitude                         ( rad )
    rlon   (:,:),   & ! geographical longitude                        ( rad )
    rlattot(:,:),   & ! geographical latitude                         ( rad )
    rlontot(:,:),   & ! geographical longitude                        ( rad )
    fccos  (:,:),   & ! horizontal coriolis-parameter                 ( 1/s )
    fc     (:,:),   & ! coriolis-parameter                            ( 1/s )
    rmy    (:,:,:), & ! Davis-parameter for relaxation (mass, qv, qc)   --
    rmyq   (:,:),   & ! Davis-parameter for relaxation (qr, qs, qg)     --
    hd_mask_dcoeff_p(:,:,:), & ! 3D-domain mask for horizontal diffusion  --
    hd_mask_dcoeff_t(:,:,:), & ! 3D-domain mask for horizontal diffusion  --
    hd_mask_dcoeff_q(:,:,:), & ! 3D-domain mask for horizontal diffusion  --
    hd_mask_dcoeff_u(:,:,:), & ! 3D-domain mask for horizontal diffusion  --
    ofa_hdx(:,:,:), & !
    ofa_hdy(:,:,:), & !
    hd_mask(:,:,:)    ! 3D-domain mask for horizontal diffusion         --

#ifndef MESSY
  INTEGER,             TARGET, ALLOCATABLE ::           &
#else
  INTEGER,             POINTER             ::           &
#endif
    isoiltyp(:,:)     ! type of the soil (keys 0-9)                     --

#ifndef MESSY
  REAL  (KIND=wp),     TARGET, ALLOCATABLE ::           &
#else
  REAL  (KIND=wp),     POINTER             ::           &
#endif
    crlat  (:,:),   & ! cosine of transformed latitude                  --
    acrlat (:,:),   & ! 1 / ( crlat * radius of the earth )           ( 1/m )
    tgrlat (:,:),   & ! tangens of transformed latitude                 --
    aerlan (:,:),   & ! aerosol-distribution on rural areas             --
    aerurb (:,:),   & ! aerosol-distribution on urban areas             --
    aerdes (:,:),   & ! aerosol-distribution on desert areas            --
    aersea (:,:),   & ! aerosol-distribution on the sea                 --
    plcov  (:,:),   & ! fraction of plant cover                         --
    lai    (:,:),   & ! leaf area index of plants                       --
    tai    (:,:),   & ! transpiration area index                        --
    sai    (:,:),   & ! surface area index                              --
    eai    (:,:),   & ! (evaporative) earth area index                  --
    rootdp (:,:),   & ! depth of the roots                            ( m  )
    for_e  (:,:),   & ! ground fraction covered by evergreen forest     --
    for_d  (:,:),   & ! ground fraction covered by deciduous forest     --
    skinc  (:,:),   & ! skin conductivity                             (W m-2 K-1)
!   h_can  (:,:),   & ! hight of the vertically resolved canopy       ( m )
    d_pat  (:,:),   & ! horizontal pattern length scale               ( m )
    l_pat  (:,:),   & ! effective length scale of circulation patterns
    alb_dry(:,:),   & ! surface albedo field for dry soil
    alb_sat(:,:),   & ! surface albedo field for saturated soil
    alb_dif(:,:)      ! diffuse albedo field

  ! also for MESSY, these fields are only allocatable?
! REAL  (KIND=wp),     TARGET, ALLOCATABLE ::           &
!   c_big (:,:,:),  & ! effective drag coefficient of canopy elements
!                     ! larger than or equal to the tubulent length
!                     ! scale                                         (1/m)
!   c_sml (:,:,:),  & ! effective drag coefficient of canopy elements
!                     ! smaller than the tubulent length scale        (1/m)
!   r_air (:,:,:)     ! air containing fraction of a gridbox inside
!                     ! the canopy                                    ( 1 )

  LOGICAL, ALLOCATABLE ::           &
    least_lbdz(:,:), & ! mask for eastern  lateral boundary zone
    lwest_lbdz(:,:), & ! mask for western  lateral boundary zone
    lnorth_lbdz(:,:),& ! mask for northern lateral boundary zone
    lsouth_lbdz(:,:)   ! mask for southern lateral boundary zone

  ! external parameter fields for the lake model FLake
#ifndef MESSY
  REAL  (KIND=wp),     TARGET, ALLOCATABLE ::           &
#else
  REAL  (KIND=wp),     POINTER             ::           &
#endif
    fr_lake(:,:),   & ! lake fraction in a grid element [0,1]         (  -  )
    depth_lk(:,:),  & ! lake depth                                    (  m  )
    fetch_lk(:,:),  & ! wind fetch over lake                          (  m  )
    dp_bs_lk(:,:),  & ! depth of the thermally active layer
                      ! of bottom sediments                           (  m  )
    t_bs_lk (:,:),  & ! climatological temperature at the bottom of
                      ! the thermally active layer of sediments       (  K  )
    gamso_lk(:,:)     ! attenuation coefficient for
                      ! solar radiation in lake water                 ( 1/m )


#ifndef MESSY
  LOGICAL,             TARGET, ALLOCATABLE ::           &
#else
  LOGICAL,             POINTER ::                       &
#endif
    llandmask(:,:), & ! landpoint mask
    llakemask(:,:), & ! lake point mask
    lseamask (:,:)    ! ocean point mask, i.e. water but not lake

! 3. prognostic variables                                             (unit)
! -----------------------

#ifndef MESSY
  REAL  (KIND=wp),     TARGET, ALLOCATABLE ::           &
#else
  REAL  (KIND=wp),     POINTER             ::           &
#endif
    u (:,:,:,:),    & ! zonal wind speed                              ( m/s )
    v (:,:,:,:),    & ! meridional wind speed                         ( m/s )
    w (:,:,:,:),    & ! vertical wind speed (defined on half levels)  ( m/s )
    t (:,:,:,:),    & ! temperature                                   (  k  )
    pp(:,:,:,:),    & ! deviation from the reference pressure         ( pa  )
! needed for COSMO-ICON Physics
    u_m(:,:,:),     & ! zonal wind speed  at mass point               ( m/s )
    v_m(:,:,:),     & ! meridional wind speed at mass point           ( m/s )

! fields of the turbulent scheme defined on half-levels:

   tke(:,:,:,:),    & ! SQRT(2 * turbulent kinetik energy)            ( m/s )
   edr(:,:,:),      & ! eddy dissipation rate of TKE (EDR)            (m2/s3)
   tketens(:,:,:),  & ! non-advective tke-tendency (on half-levels)   ( m/s )
   tket_adv(:,:,:), & ! pure advective tke-tendency (on half-levels)  ( m/s )
   tket_conv(:,:,:),& ! TKE-tendency due to convective buoyancy       ( m2/s3)
   tket_hshr(:,:,:),& ! TKE-tendency due to (sep.) horiz. shear       ( m2/s3)
   tket_sso (:,:,:)   ! TKE-tendency due to SSO wake production       ( m2/s3)

! 4. tendency fields for the prognostic variables                     (unit )
! -----------------------------------------------
!    time tendencies  by diabatic and adiabatic processes
!    without sound-wave terms

#ifndef MESSY
  REAL  (KIND=wp),     ALLOCATABLE ::           &
#else
  REAL  (KIND=wp),     POINTER             ::           &
#endif
    utens (:,:,:),  & ! u-tendency without sound-wave terms           ( m/s2)
    vtens (:,:,:),  & ! v-tendency without sound-wave terms           ( m/s2)
    wtens (:,:,:),  & ! w-tendency without sound-wave terms           ( m/s2)
                      ! (defined on half levels)
    ttens (:,:,:),  & ! t-tendency without sound-wave terms           ( K/s )
    pptens(:,:,:)     ! pp-tendency without sound-wave terms          (Pa/s )


! 5. fields for surface values and soil/canopy model variables        (unit )
! -----------------------------------------------------

#ifndef MESSY
  REAL  (KIND=wp),     TARGET, ALLOCATABLE ::           &
#else
  REAL  (KIND=wp),     POINTER             ::           &
#endif
    ps       (:,:,:),  & ! surface pressure                           ( pa  )
    t_snow   (:,:,:),  & ! temperature of the snow-surface            (  K  )
    t_s      (:,:,:),  & ! temperature of the ground surface (soil)   (  K  )
    t_sk     (:,:,:),  & ! skin temperature                           (  K  )
    t_g      (:,:,:),  & ! weighted surface temperature               (  K  )
    tg_radstep (:,:),  & ! ground temperature
                         ! at the last call of the radiation routine  (  K  )
    qv_s     (:,:,:),  & ! specific water vapor content at the surface(kg/kg)
    t_cl     (:,:),    & ! temperature between medium and lower
                         ! soil layer (climatology)                   (  K  )
    t_so     (:,:,:,:),& ! multi-layer soil temperature               (  K  )
    w_snow   (:,:,:),  & ! water content of snow                      (m H2O)
    w_i      (:,:,:),  & ! water content of interception water        (m H2O)
    w_p      (:,:,:),  & ! water content of interception water pond   (m H2O)
    w_s      (:,:,:),  & ! water content of interception snow water   (m H2O)
    w_so     (:,:,:,:),& ! multi-layer soil moisture                  (m H2O)
    w_so_ice (:,:,:,:),& ! multi-layer soil ice                       (m H2O)
    w_cl     (:,:),    & ! climatological water content               (m H2O)
    t_e      (:,:,:),  & ! surface temperature of the canopy elements (  K  )
    qv_e     (:,:,:),  & ! surface value of qv of the canopy elements (Kg/Kg)
    freshsnow(:,:),    & ! weighting function indicating 'freshness' of snow
    snow_melt(:,:),    & ! snow melt amount                           (kg/m2)
    h_snow   (:,:,:),  & ! snow height                                (  m  )
    rho_snow (:,:,:),  & ! prognostic snow density                    (kg/m3)
    fr_snow  (:,:),    & ! surface fraction covered by snow           (  -  )
    fr_wi    (:,:),    & ! surface fraction covered by interception water (-)
    ustar_fv (:,:),    & ! friction velocity (ustar)                  ( m/s )
    lev_snow     (:,:,:,:),& ! vertical grid in snowpack              (  m  )
    t_snow_mult  (:,:,:,:),& ! temperature of the snow-surface        (  K  )
    dzh_snow_mult(:,:,:,:),& ! layer thickness between half levels in snow   (  m  )
    wliq_snow    (:,:,:,:),& ! liquid water content in snow           (m H2O)
    w_snow_mult  (:,:,:,:),& ! total (liquid + solid) water content of snow  (m H2O)
    rho_snow_mult(:,:,:,:)   ! prognostic snow density                (kg/m3)

  ! fields for prognostic variables of the lake model FLake or ocean
  ! variables
#ifndef MESSY
  REAL  (KIND=wp),     TARGET, ALLOCATABLE ::           &
#else
  REAL  (KIND=wp),     POINTER             ::           &
#endif
    fr_ice(:,:),    & ! ice fraction for ocean/lake surfaces          (  -  )
    t_ice (:,:,:),  & ! temperature of ice/water surface              (  K  )
    h_ice (:,:,:)     ! lake/sea ice thickness                        (  m  )

#ifndef MESSY
  REAL  (KIND=wp),     TARGET, ALLOCATABLE ::           &
#else
  REAL  (KIND=wp),     POINTER             ::           &
#endif
    t_mnw_lk(:,:,:),& ! mean temperature of the water column          (  K  )
    t_wml_lk(:,:,:),& ! mixed-layer temperature                       (  K  )
    t_bot_lk(:,:,:),& ! temperature at the water-bottom sediment
                      ! interface                                     (  K  )
    t_b1_lk (:,:,:),& ! temperature at the bottom of the upper layer
                      ! of the sediments                              (  K  )
    c_t_lk  (:,:,:),& ! shape factor with respect to the
                      ! temperature profile in lake thermocline       (  -  )
    h_ml_lk (:,:,:),& ! thickness of the mixed-layer                  (  m  )
    h_b1_lk (:,:,:)   ! thickness of the upper layer
                      ! of bottom sediments                           (  m  )

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------

#ifndef MESSY
  REAL  (KIND=wp),     TARGET, ALLOCATABLE ::           &
#else
  REAL  (KIND=wp),     POINTER             ::           &
#endif

    ttens_diab(:,:,:), & ! temperature tendency due to pure diabatic processes for output (  K/s )
    tinc_lh(:,:,:), & ! temperature increment due to latent heat      (  K  )

!   air density at present time level (main levels)       
    rho (:,:,:),    & ! total density of air                          (kg/m3)

!   coefficients for turbulent diffusion in the atmosphere
!   (defined on half levels)
                      ! vertical   turbulent diffusion coefficients
    tkvm(:,:,:),    & ! ... for momentum                              (m2/s)
    tkvh(:,:,:),    & ! ... for heat and moisture                     (m2/s)
                      ! horizontal turbulent diffusion coefficients
    tkhm(:,:,:),    & ! ... for momentum                              (m2/s)
    tkhh(:,:,:),    & ! ... for heat and moisture                     (m2/s)

    ! variables needed for turb_prepare and other turbulence issues
    l_hori (:,:),   & ! horizontal grid spacing (location dependent horizontal scale)
    hdiv (:,:,:),   & ! horizontal divergence
    hdef2(:,:,:),   & ! horizontal deformation square
    dwdx (:,:,:),   & ! horizontal gradients of vertical wind
    dwdy (:,:,:),   & ! horizontal gradients of vertical wind

    ! note that ut_turb, vt_turb are defined on the mass grid point!
    ! only needed, if vertical diffusion is computed with turbulence scheme
    ut_turb (:,:,:),& ! u-tendency from vertical diffusion in turbu ( m/s^2)
    vt_turb (:,:,:),& ! v-tendency due to convection                  ( m/s^2)

!   vertical varying implicitness of vertical diffusion
!   (def. at half levels)
    a1t(:),         & !                                               ( -- )
    a2t(:),         & !                                               ( -- )

    ! Rayleigh damping coefficient
    rdcoef(:),      & !

!   turbulence statistics in the atmosphere
!   (defined on full levels)
    rcld(:,:,:),    & ! standard deviation of the saturation deficit    --

!   turbulent transfer coefficients at the surface
    tcm (:,:),      & ! for momentum                                  ( -- )
    tch (:,:),      & ! for heat and moisture                         ( -- )
!   turbulent transfer velocities at the surface
    tvm (:,:),      & ! for momentum                                  ( m/s)
    tvh (:,:),      & ! for heat and moisture                         ( m/s)
    tkr (:,:),      & ! reference surf. diff. coeff. (l*Ustar)        (m2/s)
    tkred_sfc(:,:), & ! reduction factor for minimum diffusion coefficients near the surface
!   turbulent transfer factors for laminar- and roughness-layer transfer
    tfm (:,:),      & ! of momentum                                     --
    tfh (:,:),      & ! of scalars                                      --
    tfv (:,:),      & ! of water vapor compared to heat                 --

!   random noise on idealized surface fluxes
    h0noise(:,:),   &

!   fields from the radiation scheme
    sohr          (:,:,:), & ! rate of solar heating                           ( K/s )
    sotr          (:,:,:), & ! solar transmissivity
    sotr_par      (:,:)  , & ! solar transmissivity, photosynthetic active radiation
!US these fields belong to a development which did never make it into the official
!   version. But if this development will be activated again (see also subroutine
!   radiation_average in radiation_interface) they are needed as global fields
!   Components of photosynthetic active radiation surface flux
!   flpar_s_dir   (:,:)  , & ! direct component (aka parallel component)         ( W/m2)
!   flpar_s_difd  (:,:)  , & ! diffuse downward component                        ( W/m2)
!   flpar_s_difu  (:,:)  , & ! diffuse upward component                          ( W/m2)
    thhr          (:,:,:), & ! rate of thermal heating                           ( K/s )
    sodwddm       (:,:)  , & ! downward direct solar radiative flux / smu0       ( W/m2)
    qc_rad        (:,:,:), & ! subgrid-scale specific cloud liq. water cont.     (kg/kg)
    qi_rad        (:,:,:), & ! subgrid-scale specific ice water                  (kg/kg)
#ifdef CLOUDRAD
    qr_rad        (:,:,:), & ! subgrid-scale specific rain content               (kg/kg)
    qs_rad        (:,:,:), & ! subgrid-scale specific snow content               (kg/kg)
    qg_rad        (:,:,:), & ! subgrid-scale specific graupel content            (kg/kg)
    qc_rad_con    (:,:,:), &
    qc_rad_adiab  (:,:,:), &
    weff_rad      (:,:,:), &
#endif
    clc_sgs       (:,:,:), & ! subgrid-scale stratiform cloud cover                --
    alb_rad       (:,:)  , & ! albedo of the ground                                --
    alb_rad_rc    (:,:)  , & ! albedo of ground on coarse grid (nradcoarse>1)      --
    tg_rc         (:,:)  , & ! ground temperature on the coarse radiation grid
    sobs          (:,:)  , & ! solar radiation at the ground                     ( W/m2)
    thbs          (:,:)  , & ! thermal radiation at the ground                   ( W/m2)
    pabs          (:,:)  , & ! photosynthetic active radiation at the ground     ( W/m2)
    sobt          (:,:)  , & ! solar radiation at upper boundary of atmosphere   ( W/m2)
    thbt          (:,:)  , & ! thermal radiation at upper boundary of atmosphere ( W/m2)
    clch          (:,:)  , & ! cloud cover with high clouds                        --   
    clcm          (:,:)  , & ! cloud cover with medium clouds                      --   
    clcl          (:,:)  , & ! cloud cover with low clouds                         --   
    clct          (:,:)  , & ! total cloud cover                                   --   
    sun_el        (:,:)  , & ! sun elevation angle                               ( deg )
    sun_azi       (:,:)      ! sun azimuth  angle                                ( deg )

#ifndef MESSY
  REAL  (KIND=wp),     TARGET, ALLOCATABLE ::           &
#else
  REAL  (KIND=wp),     POINTER             ::           &
#endif
!   fields from the convection scheme
    clc_con(:,:,:), & ! cloud cover due to convection                   --     
    clw_con(:,:,:), & ! convective cloud liquid water
    prr_con(:,:),   & ! precipitation rate of rain, convective        (kg/m2*s)
    prs_con(:,:),   & ! precipitation rate of snow, convective        (kg/m2*s)
    prne_con(:,:),  & ! precipitation rate, no evaporat., convective  (kg/m2*s)
    bas_con(:,:),   & ! level index of convective cloud base            --
    top_con(:,:),   & ! level index of convective cloud top             --
    tt_conv (:,:,:),& ! temperature tendency due to convection        ( K/s  )
    ttens_conv (:,:,:),& ! temperature tendency for IFS convection ( K/s  )
    ttdiab_conv(:,:,:),& ! pure diabatic temperature tendency due to convection    ( K/s  )
    qvt_conv(:,:,:),& ! humidity    tendency due to convection        ( 1/s  )
    qct_conv(:,:,:),& ! qc tendency due to convection                 ( 1/s  )
    qit_conv(:,:,:),& ! qi tendency due to convection                 ( 1/s  )
    qrt_conv(:,:,:),& ! qr tendency due to convection                 ( 1/s  )
    qst_conv(:,:,:),& ! qs tendency due to convection                 ( 1/s  )
! note that ut_conv, vt_conv are defined on the mass grid point!
    ut_conv (:,:,:),& ! u-tendency due to convection                  ( m/s^2)
    vt_conv (:,:,:),& ! v-tendency due to convection                  ( m/s^2)
    mflx_con(:,:),  & ! cloud base massflux                           (kg/m2*s)
    cape_con(:,:),  & ! convective available energy                   (   J/kg)
    tke_con (:,:),  & ! convective turbulent energy                   (   J/kg)
    qcvg_con(:,:),  & ! moisture convergence for Kuo-type closure     (    1/s)
    w0avg   (:,:,:),& ! running average of w
    cu_evap (:)       ! factor for evaporation of rain

#ifndef MESSY
  REAL  (KIND=wp),     TARGET, ALLOCATABLE ::           &
#else
  REAL  (KIND=wp),     POINTER             ::           &
#endif
!   fields from the grid-scale precipitation scheme
    qrs    (:,:,:), & ! precipitation water content (water loading)   (kg/kg)
    prr_gsp(:,:),   & ! precipitation rate of rain, grid-scale        (kg/m2*s)
    prs_gsp(:,:),   & ! precipitation rate of snow, grid-scale        (kg/m2*s)
    prg_gsp(:,:),   & ! precipitation rate of graupel, grid-scale     (kg/m2*s)
    prh_gsp(:,:)      ! precipitation rate of hail, grid-scale        (kg/m2*s)
                      !   (only for the 2-moment scheme)

#ifndef MESSY
  REAL  (KIND=wp),     TARGET, ALLOCATABLE ::           &
#else
  REAL  (KIND=wp),     POINTER             ::           &
#endif
!   fields from the sub-grid scale orography scheme
! note that ut_sso, vt_sso are defined on the mass grid point!
    ut_sso (:,:,:), & ! u-tendency due to SSO                         ( m/s2)
    vt_sso (:,:,:), & ! v-tendency due to SSO                         ( m/s2)
    tt_sso (:,:,:), & ! temperature tendency due to SSO               ( K/s )
    ustr_sso (:,:), & ! u-stress (surface momentum flux) due to SSO   ( N/m2)
    vstr_sso (:,:), & ! v-stress (surface momentum flux) due to SSO   ( N/m2)
    vdis_sso (:,:), & ! vert. int. dissipation of kin. en. due to SSO ( W/m2)
    austr_sso(:,:), & ! average of ustr_sso                           ( N/m2)
    avstr_sso(:,:), & ! average of vstr_sso                           ( N/m2)
    avdis_sso(:,:)    ! average of vdis_sso                           ( W/m2)

#ifdef CLOUDRAD
  REAL  (KIND=wp),     POINTER             ::           &
   reffc_out(:,:,:),   &   ! effective radius of cloud droplets   ( m)  
   reffi_out(:,:,:),   &   ! effective radius of ice particles    ( m)
   reffr_out(:,:,:),   &   ! effective radius of rain drops       ( m)
   reffs_out(:,:,:),   &   ! effective radius of snow flakes      ( m)
   reffg_out(:,:,:),   &   ! effective radius of ice graupel      ( m)
   odepthw_so (:,:,:), &   ! optical depth of cloud droplets      ( 1)
   odepthw_th (:,:,:), &   ! optical depth of cloud droplets      ( 1)
   odepthi_so (:,:,:), &   ! optical depth of ice particles       ( 1)
   odepthi_th (:,:,:), &   ! optical depth of ice particles       ( 1)
   odepthr_so (:,:,:), &   ! optical depth of rain drops          ( 1)
   odepthr_th (:,:,:), &   ! optical depth of rain drops          ( 1)
   odepths_so (:,:,:), &   ! optical depth of snow flakes         ( 1)
   odepths_th (:,:,:), &   ! optical depth of snow flakes         ( 1)
   odepthg_so (:,:,:), &   ! optical depth of graupel             ( 1)
   odepthg_th (:,:,:)      ! optical depth of graupel             ( 1)
#endif

#ifndef MESSY
  REAL  (KIND=wp),     TARGET, ALLOCATABLE ::           &
#else
  REAL  (KIND=wp),     POINTER             ::           &
#endif
!   fields for the radiation correction scheme
    ! these are actual values
    swdir_s  (:,:), & ! direct comp. of solar radiative flux at surface ( W/m2)
    swdifd_s (:,:), & ! diffuse downward comp. of short wave rad. flux  ( W/m2)
    swdifu_s (:,:), & ! diffuse upward   comp. of short wave rad. flux  ( W/m2)
    swtrdir_s(:,:), & ! direct comp. of solar radiative transm. at surface
    swtrdifd_s(:,:),& ! diffuse downward comp. of short wave rad. transm.
    swtrdifu_s(:,:),& ! diffuse upward   comp. of short wave rad. transm.
    lwd_s    (:,:), & !         downward comp. of long  wave rad. flux  ( W/m2)
    lwu_s    (:,:), & !         upward   comp. of long  wave rad. flux  ( W/m2)

    ! these are accumulated values
    aswdir_s (:,:), & ! direct comp. of solar radiative flux at surface ( W/m2)
    aswdifd_s(:,:), & ! diffuse downward comp. of short wave rad. flux  ( W/m2)
    aswdifu_s(:,:), & ! diffuse upward   comp. of short wave rad. flux  ( W/m2)
    alwd_s   (:,:), & !         downward comp. of long  wave rad. flux  ( W/m2)
    alwu_s   (:,:), & !         upward   comp. of long  wave rad. flux  ( W/m2)

    ! this is the essential correction factor
    swdir_cor(:,:), & ! direct short wave radiation correction factor actual value

    ! these are topographic parameters
    skyview  (:,:), & ! sky view
    slo_asp  (:,:), & ! slope aspect
    slo_ang  (:,:), & ! slope angle
    horizon(:,:,:), & ! horizon

!   fields that are computed in the dynamics and / or physics
    dqvdt     (:,:,:),& ! threedimensional moisture convergence         ( 1/s )
    dqvdt_conv(:,:,:),& ! 3d convective moisture convergence            ( 1/s )
    qvfl_s    (:,:),  & ! surface flux of water vapour                  (kg/m2s)
    dpsdt     (:,:),  & ! tendency of the surface pressure              ( pa/s)
    umfl_s    (:,:),  & ! u-momentum flux (surface)                     ( N/m2)
    vmfl_s    (:,:),  & ! v-momentum flux (surface)                     ( N/m2)
    qmomflux  (:,:),  & ! momentum flux at surface (for FLake)          ( N/m2)
    shfl_s    (:,:),  & ! sensible heat flux (surface)                  ( W/m2)
    lhfl_s    (:,:),  & ! latent heat flux (surface)                    ( W/m2)
    aumfl_s   (:,:),  & ! average u-momentum flux (surface)             ( N/m2)
    avmfl_s   (:,:),  & ! average v-momentum flux (surface)             ( N/m2)
    ashfl_s   (:,:),  & ! average sensible heat flux (surface)          ( W/m2)
    alhfl_s   (:,:),  & ! average latent heat flux (surface)            ( W/m2)
    rstom     (:,:),  & ! stomata resistance                            ( s/m )
    lhfl_bs   (:,:),  & ! latent heat flux from bare soil evap.         ( W/m2)
    lhfl_pl   (:,:,:),& ! latent heat flux from plants                  ( W/m2)
    alhfl_bs  (:,:),  & ! average latent heat flux from bare soil evap. ( W/m2)
    alhfl_pl  (:,:,:),& ! average latent heat flux from plants          ( W/m2)

    ! used in "prepare" for Tiedtke convection (for lconf_avg)
    w_conv    (:,:,:),& ! 3d vertical velocity on full levels for conv
    qhfl      (:,:)     ! surface flux of water vapour for convection

  ! fields used in the Runge-Kutta scheme (only per time step)
  REAL (KIND = wp),     ALLOCATABLE :: &
    wcon(:,:,:),    & ! contravariant vertical velocity
    uadvt(:,:,:),   & ! advective tendency of u
    vadvt(:,:,:),   & ! advective tendency of v
    wadvt(:,:,:),   & ! advective tendency of w
    ppadvt(:,:,:),  & ! advective tendency of pp
    tadvt(:,:,:)      ! advective tendency of t

! 7. fields for model output and diagnostics                          (unit )
! ------------------------------------------

#ifndef MESSY
  REAL  (KIND=wp),     TARGET, ALLOCATABLE ::           &
#else
  REAL  (KIND=wp),     POINTER             ::           &
#endif
    t_2m    (:,:),  & ! temperature in 2m                             (  K  )
    t_2m_av (:,:),  & ! time mean temperature in 2m                   (  K  )
    qv_2m   (:,:),  & ! specific water vapor content in 2m            (kg/kg)
    td_2m   (:,:),  & ! dew-point in 2m                               (  K  )
    td_2m_av(:,:),  & ! time mean dew-point in 2m                     (  K  )
    rh_2m   (:,:),  & ! relative humidity in 2m                       (  %  )
    u_10m   (:,:),  & ! zonal wind in 10m                             ( m/s )
    u_10m_av(:,:),  & ! time_mean zonal wind in 10m                   ( m/s )
    v_10m   (:,:),  & ! meridional wind in 10m                        ( m/s )
    v_10m_av(:,:),  & ! time mean meridional wind in 10m              ( m/s )
    tmin_2m (:,:),  & ! minimum temperature in 2m                     (  K  )
    tmax_2m (:,:),  & ! maximum temperature in 2m                     (  K  )
    vmax_10m(:,:),  & ! maximal wind gust in 10m                      ( m/s )
    vabsmx_10m(:,:),& ! maximal wind speed in 10m                     ( m/s )
    vgust_dyn(:,:), & ! maximal dynamical wind gust in 10m            ( m/s )
    vgust_con(:,:), & ! maximal convective wind gust in 10m           ( m/s )
    asob_s  (:,:),  & ! average solar radiation budget (surface)      ( W/m2)
    athb_s  (:,:),  & ! average thermal radiation budget (surface)    ( W/m2)
    apab_s  (:,:),  & ! average photosynthetic active radiation (sfc) ( W/m2)
    asob_t  (:,:),  & ! average solar radiation budget (model top)    ( W/m2) 
    athb_t  (:,:),  & ! average thermal radiation budget (model top)  ( W/m2)
     sod_t  (:,:),  & ! solar downward radiation at top of atmosphere (     )
    asod_t  (:,:),  & ! averaged solar downward radiation at top      (     )
    dursun  (:,:),  & ! sunshine duration                             (  s  )
    dursun_m(:,:),  & ! maximum possible sunshine duration            (  s  )
    dursun_r(:,:),  & ! relative sunshine duration                    (  s  )
    rain_gsp(:,:),  & ! amount of rain from grid-scale precip. (sum)  (kg/m2)
    snow_gsp(:,:),  & ! amount of snow from grid-scale precip. (sum)  (kg/m2)
    grau_gsp(:,:),  & ! amount of graupel from grid-scale precip. (sum) (kg/m2)
    hail_gsp(:,:),  & ! amount of hail from grid-scale precip. (sum)  (kg/m2)
                      !   (only for the 2-moment scheme)
    rain_con(:,:),  & ! amount of rain from convective precip. (sum)  (kg/m2)
    snow_con(:,:),  & ! amount of snow from convective precip. (sum)  (kg/m2)

    uhel_track(:,:),& ! updraft helicity track (max over last output interval)    (J/kg)
    rot_track(:,:), & ! rotation track (max over last output interval)            (1/s)
    vih_track(:,:), & ! track of total column hydrometeors (max over last output interval) (kg/m2)
    vii_track(:,:), & ! track of total column hydrometeors T < -10 deg C (max over last output interval) (kg/m2)
    wmax_track(:,:),& ! track of max. vertical velocity (last output interval)    (m/s)
    lpi_track(:,:) ,& ! track of max. vertical velocity (last output interval)    (m/s)
    dbz_track(:,:) ,& ! track of column max. radar reflectivity                   (dBZ)
    echotop_z(:,:,:),  & ! track of dBZ echo tops (height) of some thresholds     (m MSL) 
    echotop_p(:,:,:),  & ! track of dBZ echo tops (pressure) of some thresholds   (hPa) 
    mconv_hor_track(:,:) ,& ! track of 1km-layer-averaged horizontal moisture convergence (kg/(m3 s))

#ifdef RADARFWO
    comp_dbzobs(:,:,:), & ! radar reflectivity composite (obs) from the radar operator (dBZ)
    comp_dbzsim(:,:,:), & ! radar reflectivity composite (sim) from the radar operator (dBZ)
#endif

    runoff_s(:,:),  & ! surface water runoff; sum over forecast       (kg/m2)
    runoff_g(:,:),  & ! soil water runoff; sum over forecast          (kg/m2)
    tdiv_hum(:,:),  & ! vertical sum for  divergence of humidity      (kg/m2)
    aevap_s (:,:)     ! accumulated surface moisture flux             (kg/m2)

! 8. fields for the boundary values                                   (unit )
! ---------------------------------

#ifndef MESSY
  REAL  (KIND=wp),     TARGET, ALLOCATABLE ::           &
#else
  REAL  (KIND=wp),     POINTER             ::           &
#endif
    u_bd   (:,:,:,:), & ! boundary field for u                        ( m/s )
    v_bd   (:,:,:,:), & ! boundary field for v                        ( m/s )
    w_bd   (:,:,:,:), & ! boundary field for w                        ( m/s )
    t_bd   (:,:,:,:), & ! boundary field for t                        (  K  )
    pp_bd  (:,:,:,:), & ! boundary field for pp                       (  pa )
    qv_s_bd  (:,:,:), & ! boundary field for qv_s                     (kg/kg)
    t_snow_bd(:,:,:), & ! boundary field for t_snow                   (  K  )
    t_s_bd   (:,:,:), & ! boundary field for t_s                      (  K  )
    w_snow_bd(:,:,:), & ! boundary field for w_snow                   (m H2O)

    ! and for the CLM Version
    hmo3_bd  (:,:,:), & ! boundary field for hmo3                     (m    )
    vio3_bd  (:,:,:), & ! boundary field for vio3                     (pa O3)
    w_cl_bd  (:,:,:), & ! boundary field for w_cl                     (m H2O)
    t_cl_bd  (:,:,:), & ! boundary field for t_cl                     (  K  )
    lai_bd   (:,:,:), & ! boundary field for lai                      ( --  )
    rootdp_bd(:,:,:), & ! boundary field for rootdp                   (m    )
    plcov_bd (:,:,:)    ! boundary field for plcov                    ( --  )

#ifndef MESSY
  REAL  (KIND=wp),     TARGET, ALLOCATABLE ::           &
#else
  REAL  (KIND=wp),     POINTER             ::           &
#endif

#ifdef CLOUDRAD
    aot_bd(:,:,:,:),  & !new Aleksei Poliukhov 2017
    ssa_bd(:,:,:,:),  & !
    asy_bd(:,:,:,:),  & !
#endif
    aer_du_bd(:,:,:), & ! boundary field for dust aerosols            ( --  )
    aer_su_bd(:,:,:), & ! boundary field for sulfate aerosols         ( --  )
    aer_ss_bd(:,:,:), & ! boundary field for sea salt aerosols        ( --  )
    aer_bc_bd(:,:,:), & ! boundary field for black carbon aerosols    ( --  )
    aer_or_bd(:,:,:)    ! boundary field for organic matter aerosols  ( --  )

! 9. fields for the synthetic satellite images
! --------------------------------------------

#ifndef MESSY
  REAL  (KIND=wp),     TARGET, ALLOCATABLE ::           &
#else
  REAL  (KIND=wp),     POINTER             ::           &
#endif
    synme5 (:,:,:)  , & ! Meteosat 5
    synme6 (:,:,:)  , & ! Meteosat 6
    synme7 (:,:,:)  , & ! Meteosat 7
    synmsg (:,:,:)      ! Meteosat Second Generation

! 10. analysis increment fields
! -----------------------------

#ifndef MESSY
  REAL  (KIND=wp),     TARGET, ALLOCATABLE ::           &
#else
  REAL  (KIND=wp),     POINTER             ::           &
#endif
    ff_anai (:,:,:) , & ! wind velocity                               ( m/s )
    dd_anai (:,:,:) , & ! wind direction                              ( rad )
    t_anai  (:,:,:) , & ! temperature                                 (  k  )
    p_anai  (:,:,:) , & ! deviation from the reference pressure       ( Pa  )
    qv_anai (:,:,:) , & ! specific water vapor content                (kg/kg)
    qc_anai (:,:,:)     ! specific cloud water content (via saturation adjustm)
!   fi_anai (:,:,:) , & ! geopotential
!   pmsl_anai (:,:) , & ! mean sea level pressure
!   tqv_anai  (:,:) , & ! wind velocity
!   tqc_anai  (:,:)     ! wind velocity

!==============================================================================
!VS <
! 12. fields for the multi_layer snow cover scheme (SNOWPOLINO)
! -------------------------------

  REAL  (KIND=wp),     TARGET, ALLOCATABLE ::           &

    t_sn     (:,:,:,:) , & ! multi-layer snow temperature               (  K  )

    theta_i  (:,:,:,:) , & ! volumentric ice content                    (  -  )
    theta_w  (:,:,:,:) , & ! volumentric water content                  (  -  )
    theta_a  (:,:,:,:) , & ! volumentric air content                    (  -  )

    dzm_sn   (:,:,:,:) , & ! layer thickness                            (  m  )

    hn_sn    (:,:,:)   , & ! new snow amounts (storage)                 (  m  )

    top_sn   (:,:,:)       ! index of first (top) snow layer            (  -  )

!VS >


END MODULE data_fields
