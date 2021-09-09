!+ Interface module for organizing fields in blocked data structure for TSA
!- Modification of sfc_interface for TERRA
!------------------------------------------------------------------------------

MODULE tsa_sfc_interface

!------------------------------------------------------------------------------
!
! Description:
!  This module contains routines to initialize and execute the terra
!  scheme in blocked data format when running the COSMO-Model. The scheme
!  have to be in this blocked data format (dimensions: nproma, ke(1), (nblock))
!  and the interface takes care that all necessary input and output data for
!  the scheme are transformed from the COSMO-Model data 
!  format (ie,je,ke) to the blocked data format.
!
!  Compared to the other blocked parameterization schemes an additional 
!  copying is necessary here, because all surface schemes are now working
!  only on the special type of grid points they need. This means that
!  TERRA only works on land grid points, FLAKE only on lake points and
!  the SEAICE schemes only on grid points with sea-ice. This is in contrast
!  to the former method, where TERRA for example used the IF (llandmask(i,j))
!  construct for every loop.
!  
!  The following public routines are contained:
!
!   - tsa_sfc_init:
!      => allocates necessary data structures to hold the masks for copying 
!         land and lake grid points to the input data for the blocked schemes
!         (note that the mask for sea-ice is changing during the run time, 
!          because sea ice can melt. Therefore this mask is not allocated
!          for the whole model run time).
!      => initializes some soil variables with reasonable values
!   - tsa_sfc_init_copy (imode, ib, ipend,  ierror, yerror):
!      => copies fields to block structure (imode=1) or from block (imode=2)
!      => substitution for "register" code
!
!   - tsa_sfc_organize (ib, ipend, ierror, yerror)
!      => copies all fields to special structures used in the interfaces
!         to TERRA, FLAKE, SEAICE (and back afterwards)
!      => calls the schemes
!
!   - tsa_sfc_finalize
!      => cleaning up at the end of the program
!
!  And these routines are private: they are only called within this module:
!
!   - sfc_in_wkarr_alloc (nproma, ke_soil, ke_snow)
!   - sfc_in_wkarr_dealloc
!      => to allocate / deallocate local work arrays if -DALLOC_WKARR is set.
!         otherwise, these fields are used as local arrays in the subroutines
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  Ulrich.Schaettler@dwd.de
!
! History:
! Version      Date       Name
! ----------   ---------- ----
! V5_4e        2017-03-23 Ulrich Schaettler
!  Initial release
! V5_4f        2017-09-01 Ulrich Schaettler, Ulrich Blahak, Valentin Clement
!  Set h_snow_new_b also for landpoints to h_snow_t (output from TERRA)
!   and include it in copyFromBlock-copylist in all cases (US)
!  Bug fix in TWOMOM_SB: need prh_gsp_b instead of prg_gsp_b (UB)
!  Added update of host / device, when running on GPUs (VC)
! V5_4g        2017-11-13 Ulrich Schaettler
!  Added new local variables dum1/2fl_s_t to acc present_or_create list
!  Adaptations to latest ICON updates (2017-11-08)
! V5_4h        2017-12-15 Ulrich Schaettler
!  Declared and updated new fields tsnred for GPU
!  Allocate and compute global variables zzhls, zdzhs, zdzms in sfc_init
! V5_5         2018-02-23 Ulrich Schaettler
!  Modifications to run the full block of parameterizations on GPU
!   Moved definition of data section at the beginning
!   Commented update device / host
! V5_5a        2018-06-22 Xavier Lapillonne
!  Ported subroutine sfc_prepare to GPU
! V5_5b        2018-10-29 Xavier Lapillonne
!  Removed update host / device for Flake, but keep it for flake_init
! V5_6a        2019-05-21 Jan-Peter Schulz
!  Introduce skin temperature formulation by Schulz and Vogel (2017)
! V5_6b        2019-10-16 Ulrich Schaettler
!  Adaptations due to Re-unification with ICON
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
! 
! Declarations:
!
! Modules used:
!
USE kind_parameters,  ONLY : wp,vpp

USE data_block_fields,ONLY :                                                &
    ! in and inout variables
    fr_land_b       , fr_lake_b       , depth_lk_b      , fc_b            , &
    isoiltyp_b      , plcov_b         , rootdp_b        , sai_b           , &
    tai_b           , eai_b           , skinc_b         , rsmin2d_b       , &
    u_m_b           , v_m_b           , t_b             , p0_b            , &
    ps_b            , t_s_b           , t_sk_b          , t_snow_b        , &
    t_snow_mult_b   ,                                                       &
    t_g_b           , qv_s_b          , w_snow_b        , rho_snow_b      , &
    rho_snow_mult_b , h_snow_b        , w_i_b           , w_p_b           , &
    w_s_b           , t_so_b          , w_so_b          , w_so_ice_b      , &
    u_10m_b         , v_10m_b         , freshsnow_b     , fr_snow_b       , &
    wliq_snow_b     , w_snow_mult_b   , dzh_snow_mult_b , prr_con_b       , &
    prs_con_b       , prr_gsp_b       , prs_gsp_b       , prg_gsp_b       , &
    prh_gsp_b       , qvfl_s_b        ,                                     &
    tch_b           , tcm_b           , tfv_b           , sobs_b          , &
    thbs_b          , pabs_b          , runoff_s_b      , runoff_g_b      , &
    shfl_s_b        , lhfl_s_b        , lhfl_bs_b       , lhfl_pl_b       , &
    rstom_b         , ptot_b          , pp_b            , qmomflux_b      , &
    fetch_lk_b      , dp_bs_lk_b      , t_bs_lk_b       , gamso_lk_b      , &
    t_mnw_lk_b      , t_mnw_lk_new_b  , t_wml_lk_b      , t_wml_lk_new_b  , &
    t_bot_lk_b      , t_bot_lk_new_b  , t_b1_lk_b       , t_b1_lk_new_b   , &
    c_t_lk_b        , c_t_lk_new_b    , h_ml_lk_b       , h_ml_lk_new_b   , &
    h_b1_lk_b       , h_b1_lk_new_b   , umfl_s_b        , vmfl_s_b        , &
    t_ice_b         , h_ice_b         , sso_sigma_b     , snow_melt_b     , &
    ! output variables
    t_snow_new_b        , h_snow_new_b        , t_ice_new_b , h_ice_new_b , &
    t_s_new_b           , t_sk_new_b          , t_snow_new_b,               &
    t_snow_mult_new_b   ,                                                   &
    w_snow_new_b        , rho_snow_new_b      , rho_snow_mult_new_b ,       &
    w_i_new_b           , w_p_new_b           , w_s_new_b           ,       &
    t_so_new_b          , w_so_new_b          , w_so_ice_new_b      ,       &
    wliq_snow_new_b     , w_snow_mult_new_b   , dzh_snow_mult_new_b ,       &
    t_g_new_b

!VS <
USE data_block_fields,ONLY :                                                &
    ! in and inout variables
    t_sn_b          , theta_i_b       , theta_w_b       , theta_a_b       , &
    dzm_sn_b        , top_sn_b        , hn_sn_b                           , &

    ! output variables                                                   
    t_sn_new_b          , theta_i_new_b       , theta_w_new_b ,             &
    theta_a_new_b       , dzm_sn_new_b        , top_sn_new_b,               &
    hn_sn_new_b
!VS >

USE data_fields


USE data_modelconfig, ONLY : dt, dt2, ke, idt_qv,                           &
                             ke_soil, ke_snow, czmls, czhls, ie, je, idt_qv,&
                             msoilgrib

!USE data_parallel,    ONLY : my_cart_id

USE data_runcontrol,  ONLY : idbg_level, ldebug_soi, lprintdeb_all,         &
                             l2tls, ibot_w_so, itype_turb, itype_tran,      &
                             ntstep, nblock, nproma, nlastproma, llake,     &
                             lseaice, itype_gscp, lmulti_snow, itype_vdif,  &
                             itype_canopy, nnow, nold, nnew, ltur
!VS <
USE data_runcontrol,  ONLY : lsnow
!VS >


 
USE sfc_terra_data,   ONLY : lsoilinit_dfi, zzhls, zdzhs, zdzms, cf_snow
USE sfc_terra,        ONLY : terra
USE sfc_terra_init,   ONLY : terra_init
USE sfc_utilities,    ONLY : diag_snowfrac_tg

#ifdef ALLOC_WKARR
USE sfc_terra_data,   ONLY:  terra_wkarr_alloc,   terra_wkarr_dealloc
!USE sfc_flake_data,   ONLY:  flake_wkarr_alloc,   flake_wkarr_dealloc
!USE sfc_seaice,       ONLY:  seaice_wkarr_alloc,  seaice_wkarr_dealloc
!VS <
USE sfc_snow_data,    ONLY:  snow_wkarr_alloc,    snow_wkarr_dealloc
!VS >
#endif
!VS <
USE sfc_snow,          ONLY: snowpolino
USE sfc_snow_data,     ONLY: lh_s
USE sfc_snow_init,     ONLY: snow_init
!VS >

USE tsa_data,         ONLY : z0, z0_b, qv, qv_b, ymodel, mind_ilon, mind_jlat

!USE environment,       ONLY: model_abort

!==============================================================================

IMPLICIT NONE

!------------------------------------------------------------------------------

PRIVATE

!==============================================================================

!------------------------------------------------------------------------------
! Global (i.e. public) Declarations:
!!DL PUBLIC :: sfc_init, sfc_prepare, sfc_organize, sfc_init_copy, sfc_finalize, &
!!DL           sfcCopyList, mind_landgp, mind_lakegp, mind_sicegp, n_landgp,     &
!!DL           n_lakegp
PUBLIC :: tsa_sfc_init, tsa_sfc_organize, tsa_sfc_init_copy, tsa_sfc_finalize, &
          mind_landgp, n_landgp                       

!==============================================================================
! Module variables
!!DL TYPE(CopylistStruct) :: sfcCopyList

!!DL INTEGER, ALLOCATABLE :: mind_lakegp(:,:), & ! to set up list of lake points
!!DL                            n_lakegp(:)  , & ! number of lake points per block
!!DL                         mind_landgp(:,:), & ! to set up list of land points
!!DL                            n_landgp(:)  , & ! number of land points per block
!!DL                         mind_sicegp(:)      ! number of sea ice points per block
 
INTEGER, ALLOCATABLE :: mind_landgp(:,:), & ! to set up list of land points
                        n_landgp(:)         ! number of land points per block

#ifdef ALLOC_WKARR
INTEGER       , DIMENSION(:), ALLOCATABLE   ::                                           &
    isoiltyp_t

LOGICAL       , DIMENSION(:), ALLOCATABLE   ::                                           &
    liceana_t

REAL (KIND=wp), DIMENSION(:), ALLOCATABLE   ::                                           &
    rootdp_t       , rsmin2d_t      , ps_t           , h_snow_gp_t    , meltrate_t     , &
    w_p_now_t      , w_p_new_t      , w_s_now_t      , w_s_new_t      , u_10m_t        , &
    v_10m_t        , zf_snow_t      , prr_con_t      , prs_con_t      , conv_frac_t    , &
    prr_gsp_t      , prs_gsp_t      , prg_gsp_t      , prh_gsp_t      , sobs_t         , &
    thbs_t         , pabs_t         , sso_sigma_t    , tfv_t          , shfl_snow_t    , &
    lhfl_snow_t    , lhfl_bs_t      , rstom_t        , u_t            , v_t            , &
    t_t            , qv_t           , ptot_t         , ai_uf_t        , alb_red_uf_t   , &
    snow_melt_t    , tsnred_t

! VS <
!REAL (KIND=wp), DIMENSION(:), ALLOCATABLE   ::                                          &
!   swdir_s_t       , swdifd_s_t     , swdifu_s_t     , lwd_s_t        , lwu_s_t

REAL (KIND=wp), DIMENSION(:,:), ALLOCATABLE               ::                            &
    t_sn_now_t     , t_sn_new_t                                                        , &
    theta_i_now_t  , theta_i_new_t                                                     , &
    theta_w_now_t  , theta_w_new_t                                                     , &
    theta_a_now_t  , theta_a_new_t                                                     , &
    dzm_sn_now_t   , dzm_sn_new_t

REAL (KIND=wp), DIMENSION(:), ALLOCATABLE                  ::                           &
    hn_sn_now_t    , hn_sn_new_t

REAL (KIND=wp), DIMENSION(:), ALLOCATABLE                 ::                            &
    top_sn_now_t   , top_sn_new_t
!VS >


REAL (KIND=wp), DIMENSION(:),   ALLOCATABLE             ::                               &
    plcov_t        , sai_t          , tai_t          , eai_t          , skinc_t        , &
    tcm_t          , tch_t          , t_snow_now_t   , t_snow_new_t   , t_s_now_t      , &
    t_s_new_t      , t_sk_now_t     , t_sk_new_t     , t_g_now_t      , qv_s_t         , &
    w_snow_now_t   , w_snow_new_t   ,                                                    &
    rho_snow_now_t , rho_snow_new_t , h_snow_t       , w_i_now_t      , w_i_new_t      , &
    freshsnow_t    , runoff_s_t     , runoff_g_t     , w_imp_t        , w_isa_t        , &
    fr_paved_t     , sa_uf_t        , shfl_s_t       , lhfl_s_t       , qvfl_s_t       , &
    t_g_new_t      , dum1fl_s_t     , dum2fl_s_t

REAL (KIND=wp), DIMENSION(:,:), ALLOCATABLE               ::                             &
    t_so_now_t     , t_so_new_t

REAL (KIND=wp), DIMENSION(:,:), ALLOCATABLE               ::                             &
    w_so_now_t     , w_so_new_t     , w_so_ice_now_t , w_so_ice_new_t

REAL (KIND=wp), DIMENSION(:,:),   ALLOCATABLE             ::                             &
    lhfl_pl_t

REAL (KIND=wp), DIMENSION(:,:), ALLOCATABLE               ::                             &
    t_snow_mult_now_t, t_snow_mult_new_t

REAL (KIND=wp), DIMENSION(:,:), ALLOCATABLE               ::                             &
    rho_snow_mult_now_t, rho_snow_mult_new_t, wliq_snow_now_t    , wliq_snow_new_t    ,  &
    wtot_snow_now_t    , wtot_snow_new_t    , dzh_snow_now_t     , dzh_snow_new_t


!!DL REAL (KIND=wp), DIMENSION(:), ALLOCATABLE   ::                                           &
!!DL     ! lake input fields
!!DL     fr_lake_t          , depth_lk_t         , fetch_lk_t         , dp_bs_lk_t         ,  &
!!DL     t_bs_lk_t          , gamso_lk_t         , fc_t               , qmomflux_t         ,  &
!!DL     ! lake fields for previous time step
!!DL     t_snow_t_p         , h_snow_t_p         , t_ice_t_p          , h_ice_t_p          ,  &
!!DL     t_mnw_lk_t_p       , t_wml_lk_t_p       , t_bot_lk_t_p       , c_t_lk_t_p         ,  &
!!DL     h_ml_lk_t_p        , t_b1_lk_t_p        , h_b1_lk_t_p        , t_s_t_p            ,  &
!!DL     ! lake fields for new time step
!!DL     t_snow_t_n         , h_snow_t_n         , t_ice_t_n          , h_ice_t_n          ,  &
!!DL     t_mnw_lk_t_n       , t_wml_lk_t_n       , t_bot_lk_t_n       , c_t_lk_t_n         ,  &
!!DL     h_ml_lk_t_n        , t_b1_lk_t_n        , h_b1_lk_t_n        , t_s_t_n            ,  &
!!DL     ! and some sea ice fields
!!DL     fr_ice_t           , albsi_t_p          , albsi_t_n

! local fields for the interface 
! these are not activated in COSMO right now, therefore they are local
REAL (KIND=wp), ALLOCATABLE                 ::                                           &
    h_snow_gp_b        (:),   & ! if snow tiles are considered this is needed for snow aging
    conv_frac_b        (:),   & ! convective area fraction
                                     ! could be set depend on tropics or not???
    ! the next variables are output from terra, but not further processed in COSMO
    meltrate_b         (:),   & ! snow melting rate
    shfl_snow_b        (:),   & ! sensible heat flux sfc (snow covered)
    lhfl_snow_b        (:)      ! latent heat flux sfc   (snow covered)
#endif

!==============================================================================
! Module procedures in "tsa_sfc_interface" 
!==============================================================================

CONTAINS 

!==============================================================================
!+ Module procedure "tsa_sfc_init" in "tsa_sfc_interface" 
!------------------------------------------------------------------------------    

SUBROUTINE tsa_sfc_init

!------------------------------------------------------------------------------
!
! Description:
!   initialize the surface schemes
!
!------------------------------------------------------------------------------

! Locals
! ------

INTEGER :: izdebug,    & ! local debug level
           istat,      & ! error variable
           ipend,      & ! end of block
           il, ip, ib, & ! loop indices for blocks
           i,j, kso      ! indices in (i,j) data format

!----------- End of header ----------------------------------------------------

  istat = 0

  ! init debugging message level
  izdebug = idbg_level

  ! GPU allocation and update on all processors for the following fields
  ! is done in routines from acc_global_data

  ! Allocate and compute additional fields for the model geometry
  ALLOCATE (zzhls          (ke_soil+1)       , & ! depth of the half level soil layers in m
            zdzhs          (ke_soil+1)       , & ! layer thickness between half levels
            zdzms          (ke_soil+1)  )        ! distance between main levels

  zzhls(1) = 2.0_wp*czmls(1) ! depth of first half level
  zdzhs(1) = zzhls(1)        ! layer thickness betw. half levels of uppermost layer
  zdzms(1) = czmls(1)        ! layer thickness between soil surface and main level
                             !   of uppermost layer

  DO kso = 2,ke_soil+1
    zzhls(kso) = zzhls(kso-1) + 2.0_wp*(czmls(kso) -zzhls(kso-1))
    zdzhs(kso) = zzhls(kso) - zzhls(kso-1) ! layer thickness betw. half levels
    zdzms(kso) = czmls(kso) - czmls(kso-1) ! layer thickness betw. main levels
  ENDDO

  ! Allocate and compute n_landgp, mind_landgp
  ALLOCATE (   n_landgp(nblock))
  ALLOCATE (mind_landgp(nproma,nblock))

  ! loop over the full domain to search for land points
  n_landgp(:)      =  0
  mind_landgp(:,:) = -1

  DO ib = 1, nblock
    IF (ib == nblock) THEN
      ipend = nlastproma
    ELSE
      ipend = nproma
    ENDIF

    il = 0
    DO ip = 1, ipend
!!DL  IF (ymodel == 'COSMO') THEN
        i = mind_ilon(ip,ib)
        j = mind_jlat(ip,ib)
!!DL  ELSE
!!DL    i=ip
!!DL    j=ib
!!DL  ENDIF
      IF (llandmask(i,j)) THEN
        ! this is a land point
        il = il+1
        mind_landgp(il,ib) = ip
      ENDIF
    ENDDO

    n_landgp(ib) = il
  ENDDO

!######################################################################################
!!DL  IF (llake) THEN
!!DL    ! Allocate and compute n_lakegp, mind_lakegp
!!DL    ALLOCATE (   n_lakegp(nblock))
!!DL    ALLOCATE (mind_lakegp(nproma,nblock))

!!DL    ! loop over the full domain to search for lake points
!!DL    n_lakegp(:)      =  0
!!DL    mind_lakegp(:,:) = -1

!!DL    DO ib = 1, nblock
!!DL      IF (ib == nblock) THEN
!!DL        ipend = nlastproma
!!DL      ELSE
!!DL        ipend = nproma
!!DL      ENDIF

!!DL      il = 0
!!DL      DO ip = 1, ipend
!!DL        i = mind_ilon(ip,ib)
!!DL        j = mind_jlat(ip,ib)
!!DL        IF (depth_lk(i,j) > 0.0_wp) THEN
!!DL          ! this is a lake point
!!DL          il = il+1
!!DL          mind_lakegp(il,ib) = ip
!!DL        ENDIF
!!DL      ENDDO

!!DL      n_lakegp(ib) = il
!!DL    ENDDO
!!DL  ENDIF

!!DL  IF (lseaice) THEN
!!DL    ALLOCATE (mind_sicegp(nproma))
!!DL     ! not necessary to allocate it for every block, because it has to be set up
!!DL     ! in every time step again (sea ice can melt during a time step)
!!DL  ENDIF

  ! Initialize the "new" time levels for certain blocked surface variables, 
  ! which are only defined over land and would remain with the undef-value
  ! over water from the allocation otherwise (and this then destroys the grib
  ! field):
  w_snow_new_b  (:) =   0.0_wp
  h_snow_new_b  (:) =   0.0_wp
  rho_snow_new_b(:) = 250.0_wp
  w_i_new_b     (:) =   0.0_wp
  w_p_new_b     (:) =   0.0_wp
  w_s_new_b     (:) =   0.0_wp

  IF (lmulti_snow) THEN
    w_snow_mult_new_b  (:,:) =   0.0_wp
    wliq_snow_new_b    (:,:) =   0.0_wp
    rho_snow_mult_new_b(:,:) = 250.0_wp
    dzh_snow_mult_new_b(:,:) =   0.0_wp
  ENDIF

!VS <
    ! for the new multi layer snow cover scheme
    ! ---------------------------------------------
   if (lsnow) then
     t_sn_new_b         (:,:) = 0.0_wp
     theta_i_new_b      (:,:) = 0.0_wp
     theta_w_new_b      (:,:) = 0.0_wp
     theta_a_new_b      (:,:) = 0.0_wp
     dzm_sn_new_b       (:,:) = 0.0_wp
     hn_sn_new_b        (:)   = 0.0_wp
     top_sn_new_b       (:)   = 0.0_wp
   endif
!VS >




#ifdef ALLOC_WKARR
  CALL sfc_in_wkarr_alloc(nproma, ke_soil, ke_snow)
  CALL terra_wkarr_alloc (ke_soil, ke_snow, nproma, istat)
!VS <
  IF (lsnow) then
    CALL snow_wkarr_alloc (ke_soil, ke_snow, nproma, istat)
  endif
!VS >

#endif

END SUBROUTINE tsa_sfc_init

!==============================================================================
!==============================================================================
!+ Module procedure "tsa_sfc_init_copy" in "tsa_sfc_interface" 
!+ NEW for TSA - Do the copy to block and from block without "register"
!------------------------------------------------------------------------------    

SUBROUTINE tsa_sfc_init_copy (imode,ib,ipend,ierror,yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   Do all required copies to/from block for TSA
!
!------------------------------------------------------------------------------

! Subroutines arguments
! --------------------

  INTEGER, INTENT(IN)  :: imode   ! mode for copying   
                                  ! =1 : copy to block
                                  ! =2 : copy from block
  INTEGER, INTENT(IN)  :: ib      ! current block index
  INTEGER, INTENT(IN)  :: ipend   ! length of current block
  INTEGER, INTENT(OUT) :: ierror
  CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg

!------------------------------------------------------------------------------

  ! Locals
  INTEGER :: nx
  INTEGER :: i,j,k,ip

!------------------------------------------------------------------------------

! Do the copy to block

  IF (imode == 1) THEN
!!DL  nx = nnew
      nx = nnow
    DO ip = 1, ipend
!DL COSMO or ICON?
      i = mind_ilon(ip,ib)
      j = mind_jlat(ip,ib)

      fr_land_b(ip)    = fr_land(i,j)
      z0_b(ip)         = z0(i,j)                   !DL: new for parturs_newblock
      isoiltyp_b(ip)   = INT (soiltyp(i,j))
      plcov_b(ip)      = plcov(i,j)
      rootdp_b(ip)     = rootdp(i,j)
      sai_b(ip)        = sai(i,j)
      tai_b(ip)        = tai(i,j)
      eai_b(ip)        = eai(i,j)
      skinc_b(ip)      = skinc(i,j)
      rsmin2d_b(ip)    = rsmin2d(i,j)
      sso_sigma_b(ip)  = sso_sigma(i,j)

      u_m_b(ip,ke)  = u(i,j,ke,nx)
      v_m_b(ip,ke)  = v(i,j,ke,nx)
      t_b(ip,ke)    = t(i,j,ke,nx)
      pp_b(ip,ke)   = pp(i,j,ke,nx)
      qv_b(ip,ke)   = qv(i,j,ke,nx)

      p0_b(ip,ke)      = p0(i,j,ke)
      ps_b(ip)            = ps(i,j,nx)
      t_s_b(ip)           = t_s(i,j,nx)
      t_sk_b(ip)          = t_sk(i,j,nx)
      t_snow_b(ip)        = t_snow(i,j,nx)
      t_g_b(ip)           = t_g(i,j,nx)
      qv_s_b(ip)          = qv_s(i,j,nx)
      w_snow_b(ip)        = w_snow(i,j,nx)
      rho_snow_b(ip)      = rho_snow(i,j,nx)
      h_snow_b(ip)        = h_snow(i,j,nx)
      w_i_b(ip)           = w_i(i,j,nx)
      w_p_b(ip)           = w_p(i,j,nx)
      w_s_b(ip)           = w_s(i,j,nx)
      u_10m_b(ip)         = u_10m(i,j)
      v_10m_b(ip)         = v_10m(i,j)
      freshsnow_b(ip)     = freshsnow(i,j)
      fr_snow_b(ip)       = fr_snow(i,j)

      prr_con_b(ip) = prr_con(i,j)
      prr_con_b(ip) = prr_con(i,j)
      prs_con_b(ip) = prs_con(i,j)
      prr_gsp_b(ip) = prr_gsp(i,j)
      prs_gsp_b(ip) = prs_gsp(i,j)
      prg_gsp_b(ip) = prg_gsp(i,j)

      tch_b(ip) = tch(i,j)
      tcm_b(ip) = tcm(i,j)
      tfv_b(ip) = tfv(i,j)

      sobs_b(ip) = sobs(i,j)
      thbs_b(ip) = thbs(i,j)
      pabs_b(ip) = pabs(i,j)

      runoff_s_b(ip)  = runoff_s(i,j)
      runoff_g_b(ip)  = runoff_g(i,j)
      shfl_s_b(ip)    = shfl_s(i,j)
      lhfl_s_b(ip)    = lhfl_s(i,j)
      qvfl_s_b(ip)    = qvfl_s(i,j)
      lhfl_bs_b(ip)   = lhfl_bs(i,j)
      rstom_b(ip)     = rstom(i,j)
      snow_melt_b(ip) = snow_melt(i,j)

! all new time levels
      t_s_new_b(ip)           = t_s(i,j,nnew)
      t_sk_new_b(ip)          = t_sk(i,j,nnew)
      t_g_new_b(ip)           = t_g(i,j,nnew)
      t_snow_new_b(ip)        = t_snow(i,j,nnew)
      w_snow_new_b(ip)          = w_snow(i,j,nnew)
      rho_snow_new_b(ip)      = rho_snow(i,j,nnew)
      h_snow_new_b(ip)        = h_snow(i,j,nnew)
      w_i_new_b(ip)           = w_i(i,j,nnew)
      w_p_new_b(ip)           = w_p(i,j,nnew)
      w_s_new_b(ip)           = w_s(i,j,nnew)     

    END DO

#ifdef TWOMOM_SB
  IF (ALLOCATED(prh_gsp)) THEN
    DO ip = 1, ipend
      i = mind_ilon(ip,ib)
      j = mind_jlat(ip,ib)
      prh_gsp_b(ip) = prh_gsp(i,j)
    ENDDO
  END IF
#endif

    DO ip = 1, ipend
      i = mind_ilon(ip,ib)
      j = mind_jlat(ip,ib)
  
      DO k=0,ke_soil+1
        t_so_b(ip,k)     = t_so(i,j,k,nx)
        t_so_new_b(ip,k) = t_so(i,j,k,nnew)
      ENDDO
      DO k=1,ke_soil
        lhfl_pl_b(ip,k)  = lhfl_pl(i,j,k)
      ENDDO
    ENDDO

    DO ip = 1, ipend
      i = mind_ilon(ip,ib)
      j = mind_jlat(ip,ib)

      DO k=1,ke_soil+1
       w_so_b(ip,k)         = w_so(i,j,k,nx)
       w_so_ice_b(ip,k)     = w_so_ice(i,j,k,nx)
       w_so_new_b(ip,k)     = w_so(i,j,k,nnew)
       w_so_ice_new_b(ip,k) = w_so_ice(i,j,k,nnew)
      ENDDO
    ENDDO

!!  IF (lmulti_snow) THEN
      DO ip = 1, ipend
        i = mind_ilon(ip,ib)
        j = mind_jlat(ip,ib)
        DO k=1,ke_snow
          t_snow_mult_b(ip,k)     = t_snow_mult(i,j,k,nx)
          rho_snow_mult_b(ip,k)   = rho_snow_mult(i,j,k,nx)
          wliq_snow_b(ip,k)       = wliq_snow(i,j,k,nx)
          w_snow_mult_b(ip,k)     = w_snow_mult(i,j,k,nx)
          dzh_snow_mult_b(ip,k)   = dzh_snow_mult(i,j,k,nx)
     
          t_snow_mult_new_b(ip,k)   = t_snow_mult(i,j,k,nnew)
          rho_snow_mult_new_b(ip,k) = rho_snow_mult(i,j,k,nnew)
          wliq_snow_new_b(ip,k)     = wliq_snow(i,j,k,nnew)
          w_snow_mult_new_b(ip,k)   = w_snow_mult(i,j,k,nnew)
          dzh_snow_mult_b(ip,k)     = dzh_snow_mult(i,j,k,nnew)
        ENDDO
      ENDDO
!!  ENDIF

    if(lsnow) then
       do ip = 1,ipend
          i = mind_ilon(ip,ib)
          j = mind_jlat(ip,ib)
          do k=1,ke_snow

             theta_i_b(ip,k) = theta_i(i,j,k,nx)
             theta_w_b(ip,k) = theta_w(i,j,k,nx)
             theta_a_b(ip,k) = theta_a(i,j,k,nx)
             
             dzm_sn_b (ip,k) = dzm_sn(i,j,k,nx)
             t_sn_b(ip,k)    = t_sn(i,j,k,nx)
                        
             theta_i_new_b(ip,k) = theta_i(i,j,k,nnew)
             theta_w_new_b(ip,k) = theta_w(i,j,k,nnew)
             theta_a_new_b(ip,k) = theta_a(i,j,k,nnew)
             
             dzm_sn_new_b (ip,k) = dzm_sn(i,j,k,nnew)
             t_sn_new_b(ip,k)    = t_sn(i,j,k,nnew)

          enddo 

          top_sn_b(ip)  = top_sn(i,j,nx)
          top_sn_new_b(ip)  = top_sn(i,j,nnew)

       enddo
    endif

!------------------------------------------------------------------------------

! Do the copy from block

  ELSEIF (imode == 2) THEN

    nx = nnow

    DO ip = 1, ipend
      i = mind_ilon(ip,ib)
      j = mind_jlat(ip,ib)

      t_s(i,j,nx)            = t_s_b(ip)
      t_s(i,j,nnew)          = t_s_new_b(ip)
      t_sk(i,j,nx)           = t_sk_b(ip)
      t_sk(i,j,nnew)         = t_sk_new_b(ip)
      t_snow(i,j,nx)         = t_snow_b(ip)
      t_snow(i,j,nnew)       = t_snow_new_b(ip)
      t_g(i,j,nx)            = t_g_b(ip)
      t_g(i,j,nnew)          = t_g_new_b(ip)
      qv_s(i,j,nx)           = qv_s_b(ip)
      w_snow(i,j,nx)         = w_snow_b(ip)
      w_snow(i,j,nnew)       = w_snow_new_b(ip)
      rho_snow(i,j,nx)       = rho_snow_b(ip)
      rho_snow(i,j,nnew)     = rho_snow_new_b(ip)
      h_snow(i,j,nx)         = h_snow_b(ip)
      h_snow(i,j,nnew)       = h_snow_new_b(ip)
      w_i(i,j,nx)            = w_i_b(ip)
      w_i(i,j,nnew)          = w_i_new_b(ip)
      w_p(i,j,nx)            = w_p_b(ip)
      w_p(i,j,nnew)          = w_p_new_b(ip)
      w_s(i,j,nx)            = w_s_b(ip)
      w_s(i,j,nnew)          = w_s_new_b(ip)
      freshsnow(i,j)      = freshsnow_b(ip)
      fr_snow(i,j)        = fr_snow_b(ip)
      tch(i,j)               = tch_b(ip)
      tcm(i,j)               = tcm_b(ip)
      tfv(i,j)               = tfv_b(ip)
      runoff_s(i,j)          = runoff_s_b(ip)
      runoff_g(i,j)          = runoff_g_b(ip)
      shfl_s(i,j)            = shfl_s_b(ip)
      lhfl_s(i,j)            = lhfl_s_b(ip)
      qvfl_s(i,j)            = qvfl_s_b(ip)
      lhfl_bs(i,j)           = lhfl_bs_b(ip)
!!DL  lhfl_pl(i,j)           = lhfl_pl_b(ip)
      rstom(i,j)             = rstom_b(ip)
      snow_melt(i,j)         = snow_melt_b(ip)
    END DO

    DO ip = 1,ipend
      i = mind_ilon(ip,ib)
      j = mind_jlat(ip,ib)

      DO k = 0, ke_soil+1
        t_so(i,j,k,nx)     = t_so_b(ip,k)
        t_so(i,j,k,nnew)   = t_so_new_b(ip,k)
      ENDDO
      DO k = 1, ke_soil
        lhfl_pl(i,j,k)     = lhfl_pl_b(ip,k)
      ENDDO
      DO k = 1, ke_soil + 1
        w_so(i,j,k,nx)       = w_so_b(ip,k)
        w_so(i,j,k,nnew)     = w_so_new_b(ip,k)
        w_so_ice(i,j,k,nx)  = w_so_ice_b(ip,k)
        w_so_ice(i,j,k,nnew) = w_so_ice_new_b(ip,k)
      ENDDO
    ENDDO

!! IF (lsnow_mult) THEN
   DO ip = 1,ipend
      i = mind_ilon(ip,ib)
      j = mind_jlat(ip,ib)
      DO k = 1, ke_snow
        t_snow_mult(i,j,k,nx)     = t_snow_mult_b(ip,k)
        t_snow_mult(i,j,k,nnew)   = t_snow_mult_new_b(ip,k)
        rho_snow_mult(i,j,k,nx)   = rho_snow_mult_b(ip,k)
        rho_snow_mult(i,j,k,nnew) = rho_snow_mult_new_b(ip,k)
        wliq_snow(i,j,k,nx)       = wliq_snow_b(ip,k)
        wliq_snow(i,j,k,nnew)     = wliq_snow_new_b(ip,k)
        w_snow_mult(i,j,k,nx)     = w_snow_mult_b(ip,k)
        w_snow_mult(i,j,k,nnew)   = w_snow_mult_new_b(ip,k)
        dzh_snow_mult(i,j,k,nx)   = dzh_snow_mult_b(ip,k)
        dzh_snow_mult(i,j,k,nnew) = dzh_snow_mult_new_b(ip,k)
      ENDDO
   ENDDO
   if(lsnow) then
     do ip = 1,ipend
         i = mind_ilon(ip,ib)
         j = mind_jlat(ip,ib)
         do k = 1,ke_snow
            
            theta_i(i,j,k,nx)   = theta_i_b(ip,k)
            theta_i(i,j,k,nnew) = theta_i_new_b(ip,k)

            theta_w(i,j,k,nx)   = theta_w_b(ip,k)
            theta_w(i,j,k,nnew) = theta_w_new_b(ip,k)

            theta_a(i,j,k,nx)   = theta_a_b(ip,k)
            theta_a(i,j,k,nnew) = theta_a_new_b(ip,k)

            dzm_sn(i,j,k,nx)   = dzm_sn_b(ip,k)
            dzm_sn(i,j,k,nnew) = dzm_sn_new_b(ip,k)

            t_sn(i,j,k,nx)   = t_sn_b(ip,k)
            t_sn(i,j,k,nnew) = t_sn_new_b(ip,k)

         enddo
            top_sn(i,j,nx)   = top_sn_b(ip)
            top_sn(i,j,nnew) = top_sn_new_b(ip)
     enddo
   endif



!! ENDIF

!------------------------------------------------------------------------------
  
  ELSE            ! wrong imode

    WRITE (*,*) "Wrong imode in tsa_sfc_init_copy: ",imode," Has to be 1 or 2!"
    ierror = 222
    yerrmsg = 'Wrong imode in tsa_sfc_init_copy'
    RETURN

!------------------------------------------------------------------------------
  ENDIF
!------------------------------------------------------------------------------

END SUBROUTINE tsa_sfc_init_copy

!==============================================================================
!==============================================================================
!+ Module procedure to perform TERRA                                            
!------------------------------------------------------------------------------    

SUBROUTINE tsa_sfc_organize (ib, ipend, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   This module calls TERRA                                                     
!
!------------------------------------------------------------------------------

! Subroutines arguments
! --------------------

  INTEGER, INTENT(IN)  :: ib                 ! current block index
  INTEGER, INTENT(IN)  :: ipend              ! length of current block
  INTEGER, INTENT(OUT) :: ierror
  CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg

! Locals
! ------

  INTEGER :: izdebug      !local debug level

  INTEGER :: k, ip, il, ifull, n_sicegp, n_gscpclass, itl, nx, izerror, &
!VS <
             i, j, iv
!VS >


  REAL(KIND=wp) :: zdt

!!DL  REAL(KIND=wp), POINTER :: qv_b(:,:)      ! pointer to tracer field

#ifndef ALLOC_WKARR
! Local fields for additional copying
! -----------------------------------

  INTEGER       , DIMENSION(nproma)                         ::                             &
      isoiltyp_t

  LOGICAL       , DIMENSION(nproma)                         ::                             &
      liceana_t

  REAL (KIND=wp), DIMENSION(nproma)                         ::                             &
      rootdp_t       , rsmin2d_t      , ps_t           , h_snow_gp_t    , meltrate_t     , &
      w_p_now_t      , w_p_new_t      , w_s_now_t      , w_s_new_t      , u_10m_t        , &
      v_10m_t        , zf_snow_t      , prr_con_t      , prs_con_t      , conv_frac_t    , &
      prr_gsp_t      , prs_gsp_t      , prg_gsp_t      , prh_gsp_t      , sobs_t         , &
      thbs_t         , pabs_t         , sso_sigma_t    , tfv_t          , shfl_snow_t    , &
      lhfl_snow_t    , lhfl_bs_t      , rstom_t        , u_t            , v_t            , &
      t_t            , qv_t           , ptot_t         , ai_uf_t        , alb_red_uf_t   

  REAL (KIND=wp), DIMENSION(nproma)                         ::                             &
      plcov_t        , sai_t          , tai_t          , eai_t          , skinc_t        , &
      tcm_t          , tch_t          , t_snow_now_t   , t_snow_new_t   , t_s_now_t      , &
      t_s_new_t      , t_sk_now_t     , t_sk_new_t     , t_g_now_t      , qv_s_t         , &
      w_snow_now_t   , w_snow_new_t   ,                                                    &
      rho_snow_now_t , rho_snow_new_t , h_snow_t       , w_i_now_t      , w_i_new_t      , &
      freshsnow_t    , runoff_s_t     , runoff_g_t     , w_imp_t        , w_isa_t        , &
      fr_paved_t     , sa_uf_t        , shfl_s_t       , lhfl_s_t       , qvfl_s_t       , &
      t_g_new_t      , snow_melt_t    , tsnred_t       , dum1fl_s_t     , dum2fl_s_t

  REAL (KIND=wp), DIMENSION(nproma,0:ke_soil+1)             ::                             &
      t_so_now_t     , t_so_new_t
!VS <
  REAL (KIND=vpp), DIMENSION(nproma,1:ke_snow)             ::                            &
      t_sn_now_t     , t_sn_new_t                                                       , &
      theta_i_now_t  , theta_i_new_t                                                    , &
      theta_w_now_t  , theta_w_new_t                                                    , &
      theta_a_now_t  , theta_a_new_t                                                    , &
      dzm_sn_now_t   , dzm_sn_new_t

  REAL (KIND=vpp), DIMENSION(nproma)                        ::                            &
      hn_sn_now_t    , hn_sn_new_t

  REAL (KIND=vpp), DIMENSION(nproma)                        ::                            &
      top_sn_now_t   , top_sn_new_t
!VS >


  REAL (KIND=wp), DIMENSION(nproma,ke_soil+1)               ::                             &
      w_so_now_t     , w_so_new_t     , w_so_ice_now_t , w_so_ice_new_t

  REAL (KIND=wp), DIMENSION(nproma,ke_soil+1)               ::                             &
      !note: in ICON this variable is ke_soil+1, in COSMO only ke_soil
      !      make the variable given to sfc_terra therefore ke_soil+1
      lhfl_pl_t

  REAL (KIND=wp), DIMENSION(nproma,0:ke_snow)               ::                             &
      t_snow_mult_now_t, t_snow_mult_new_t

  REAL (KIND=wp), DIMENSION(nproma,ke_snow)                 ::                             &
      rho_snow_mult_now_t, rho_snow_mult_new_t, wliq_snow_now_t    , wliq_snow_new_t    ,  &
      wtot_snow_now_t    , wtot_snow_new_t    , dzh_snow_now_t     , dzh_snow_new_t


!!DL  REAL (KIND=wp), DIMENSION(nproma)                ::                                   &
!!DL    ! lake input fields
!!DL    fr_lake_t          , depth_lk_t         , fetch_lk_t         , dp_bs_lk_t         , &
!!DL    t_bs_lk_t          , gamso_lk_t         , fc_t               , qmomflux_t         , &
!!DL    ! lake fields for previous time step
!!DL    t_snow_t_p         , h_snow_t_p         , t_ice_t_p          , h_ice_t_p          , &
!!DL    t_mnw_lk_t_p       , t_wml_lk_t_p       , t_bot_lk_t_p       , c_t_lk_t_p         , &
!!DL    h_ml_lk_t_p        , t_b1_lk_t_p        , h_b1_lk_t_p        , t_s_t_p            , &
!!DL    ! lake fields for new time step
!!DL    t_snow_t_n         , h_snow_t_n         , t_ice_t_n          , h_ice_t_n          , &
!!DL    t_mnw_lk_t_n       , t_wml_lk_t_n       , t_bot_lk_t_n       , c_t_lk_t_n         , &
!!DL    h_ml_lk_t_n        , t_b1_lk_t_n        , h_b1_lk_t_n        , t_s_t_n            , &
!!DL    ! and some sea ice fields
!!DL    fr_ice_t           , albsi_t_p          , albsi_t_n

! local fields for the interface 
! these are not activated in COSMO right now, therefore they are local
  REAL (KIND=wp)       ::          &
    h_snow_gp_b        (nproma),   & ! if snow tiles are considered this is needed for snow aging
    conv_frac_b        (nproma),   & ! convective area fraction
                                     ! could be set depend on tropics or not???
    ! the next variables are output from terra, but not further processed in COSMO
    meltrate_b         (nproma),   & ! snow melting rate
    shfl_snow_b        (nproma),   & ! sensible heat flux sfc (snow covered)
    lhfl_snow_b        (nproma)      ! latent heat flux sfc   (snow covered)
#endif

!----------- End of header ----------------------------------------------------

!------------------------------------------------------------------------------
! Begin subroutine sfc_organize
!------------------------------------------------------------------------------
  
  ierror     = 0
  izerror    = 0
  yerrmsg(:) = ' '

  ! init debugging message level
  izdebug = idbg_level

  IF (izdebug > 5 .AND. ib == 1) THEN
!!DL    PRINT *, '  SURFACE SCHEMES: SOIL MODEL '
    PRINT *, '  SURFACE SCHEME: TERRA '
  ENDIF

  ! get the correct time step and select timelevel
  IF (l2tls) THEN
    nx   = nnow
    zdt  = dt
  ELSE
    nx   = nold
    zdt  = dt2
  ENDIF

!!DL   ! Retrieve required tracers in block format
!!DL   CALL trcr_get(izerror, idt_qv, ptr_tlev=nx, ptr_b=qv_b)
!!DL   ierror=MAX(ierror,izerror)

  ! determine n_gscpclass: number of considered hydrometeors
  ! ICON has different values than COSMO!
  SELECT CASE (itype_gscp)
  CASE (1)
    n_gscpclass = 3
  CASE (2)
    n_gscpclass = 4
  CASE (3)
    n_gscpclass = 5
  CASE (4)
    n_gscpclass = 6
  CASE (100:)
    ! only possible for TWOMOM_SB
    n_gscpclass = itype_gscp
  END SELECT

  !------------------------------------------------------------------------------
  ! Start GPU data region
  !------------------------------------------------------------------------------

  !$acc data                                                                         &
  !$acc present (h_snow_new_b   , mind_landgp     , n_landgp       , liceana_t     ) &
  !$acc present (isoiltyp_t     , isoiltyp_b      , rootdp_t       , rootdp_b      ) &
  !$acc present (plcov_t        , plcov_b         , sai_t          , sai_b         ) &
  !$acc present (tai_t          , tai_b           , eai_t          , eai_b         ) &
  !$acc present (skinc_t        , skinc_b                                          ) &
  !$acc present (rsmin2d_t      , rsmin2d_b       , sso_sigma_t    , sso_sigma_b   ) &
  !$acc present (ps_t           , ps_b            , t_snow_now_t   , t_snow_b      ) &
  !$acc present (t_s_now_t      , t_s_b           , t_sk_now_t     , t_sk_b        ) &
  !$acc present (t_g_now_t      , t_g_b                                            ) &
  !$acc present (qv_s_t         , qv_s_b          , w_snow_now_t   , w_snow_b      ) &
  !$acc present (rho_snow_now_t , rho_snow_b      , freshsnow_t    , freshsnow_b   ) &
  !$acc present (h_snow_t       , h_snow_b        , h_snow_gp_t    , h_snow_gp_b   ) &
  !$acc present (w_i_now_t      , w_i_b           , w_p_now_t      , w_p_b         ) &
  !$acc present (w_s_now_t      , w_s_b           , u_10m_t        , u_10m_b       ) &
  !$acc present (v_10m_t        , v_10m_b         , zf_snow_t      , fr_snow_b     ) &
  !$acc present (tch_t          , tch_b           , tcm_t          , tcm_b         ) &
  !$acc present (runoff_s_t     , runoff_s_b      , runoff_g_t     , runoff_g_b    ) &
  !$acc present (snow_melt_t    , snow_melt_b     , prr_con_t      , prr_con_b     ) &
  !$acc present (prs_con_t      , prs_con_b       , conv_frac_t    , conv_frac_b   ) &
  !$acc present (prr_gsp_t      , prr_gsp_b       , prs_gsp_t      , prs_gsp_b     ) &
  !$acc present (prg_gsp_t      , prg_gsp_b       , tfv_t          , tfv_b         ) &
  !$acc present (sobs_t         , sobs_b          , thbs_t         , thbs_b        ) &
  !$acc present (pabs_t         , pabs_b          , u_t            , u_m_b         ) &
  !$acc present (v_t            , v_m_b           , t_t            , t_b           ) &
  !$acc present (qv_t           , qv_b            , ptot_t         , ptot_b        ) &
  !$acc present (p0_b           , pp_b                                             ) &
#ifdef TWOMOM_SB
  !$acc present (prh_gsp_t      , prh_gsp_b                                        ) &
#endif
  !$acc present (lhfl_pl_t          , lhfl_pl_b                                    ) &
  !$acc present (t_so_now_t         , t_so_b                                       ) &
  !$acc present (w_so_now_t         , w_so_b                                       ) &
  !$acc present (w_so_ice_now_t     , w_so_ice_b                                   ) &
  !$acc present (t_snow_mult_now_t  , t_snow_mult_b                                ) &
  !$acc present (rho_snow_mult_now_t, rho_snow_mult_b                              ) &
  !$acc present (wliq_snow_now_t    , wliq_snow_b                                  ) &
  !$acc present (wtot_snow_now_t    , w_snow_mult_b                                ) &
  !$acc present (dzh_snow_now_t     , dzh_snow_mult_b                              ) &
  !$acc present (tsnred_t                                                          ) &

  ! Output from TERRA
  !$acc present (shfl_s_t       , shfl_s_b        , lhfl_s_t       , lhfl_s_b      ) &
  !$acc present (qvfl_s_t       , qvfl_s_b        , meltrate_t     , meltrate_b    ) &
  !$acc present (lhfl_bs_t      , lhfl_bs_b       , rstom_t        , rstom_b       ) &
  !$acc present (rho_snow_new_t , rho_snow_new_b  , t_snow_new_t   , t_snow_new_b  ) &
  !$acc present (w_i_new_t      , w_i_new_b       , w_p_new_t      , w_p_new_b     ) &
  !$acc present (w_snow_new_t   , w_snow_new_b    , t_g_new_t      , t_g_new_b     ) &
  !$acc present (t_s_new_t      , t_s_new_b       , t_sk_new_t     , t_sk_new_b    ) &
  !$acc present (w_s_new_t      , w_s_new_b                                        ) &
  !$acc present (lhfl_snow_t    , lhfl_snow_b     , shfl_snow_t    , shfl_snow_b   ) &
  !$acc present (t_so_new_t     , t_so_new_b      , w_so_new_t     , w_so_new_b    ) &
  !$acc present (w_so_ice_new_t , w_so_ice_new_b                                   ) &

  !$acc present (t_snow_mult_new_t  , t_snow_mult_new_b                            ) &
  !$acc present (rho_snow_mult_new_t, rho_snow_mult_new_b                          ) &
  !$acc present (wliq_snow_new_t    , wliq_snow_new_b                              ) &
  !$acc present (wtot_snow_new_t    , w_snow_mult_new_b                            ) &
  !$acc present (dzh_snow_new_t     , dzh_snow_mult_new_b                          )

  ! init local input fields for TERRA not yet considered in COSMO
  !$acc kernels
  h_snow_gp_b(:) = h_snow_b(:)     ! if not working with snow tiles
  conv_frac_b(:) = 0.05_wp
  !$acc end kernels

  ! and compute ptot_b again
  !$acc kernels
  ptot_b(:,ke) = p0_b(:,ke) + pp_b(:,ke)
  !$acc end kernels

!------------------------------------------------------------------------------
! TERRA 
!------------------------------------------------------------------------------
 
  IF (n_landgp(ib) > 0) THEN

    ! Copy all land points to special variables, if needed as input
    !$acc parallel
    !$acc loop gang vector private (ifull)
    DO il = 1, n_landgp(ib)
      ! grid point index in the full variable
      ifull = mind_landgp(il,ib)

      isoiltyp_t     (il)    = isoiltyp_b      (ifull)
      liceana_t      (il)    = .FALSE.
      rootdp_t       (il)    = rootdp_b        (ifull)
      plcov_t        (il)    = plcov_b         (ifull)
      sai_t          (il)    = sai_b           (ifull)
      tai_t          (il)    = tai_b           (ifull)
      eai_t          (il)    = eai_b           (ifull)
      skinc_t        (il)    = skinc_b         (ifull)
      rsmin2d_t      (il)    = rsmin2d_b       (ifull)
      sso_sigma_t    (il)    = sso_sigma_b     (ifull)

      ps_t           (il)    = ps_b            (ifull)
      t_snow_now_t   (il)    = t_snow_b        (ifull)
      t_s_now_t      (il)    = t_s_b           (ifull)
      t_sk_now_t     (il)    = t_sk_b          (ifull)
      t_g_now_t      (il)    = t_g_b           (ifull)
      qv_s_t         (il)    = qv_s_b          (ifull)
      w_snow_now_t   (il)    = w_snow_b        (ifull)
      rho_snow_now_t (il)    = rho_snow_b      (ifull)
      freshsnow_t    (il)    = freshsnow_b     (ifull)
      h_snow_t       (il)    = h_snow_b        (ifull)
      h_snow_gp_t    (il)    = h_snow_gp_b     (ifull)
      w_i_now_t      (il)    = w_i_b           (ifull)
      w_p_now_t      (il)    = w_p_b           (ifull)
      w_s_now_t      (il)    = w_s_b           (ifull)
      u_10m_t        (il)    = u_10m_b         (ifull)
      v_10m_t        (il)    = v_10m_b         (ifull)
      zf_snow_t      (il)    = fr_snow_b       (ifull)
      tch_t          (il)    = tch_b           (ifull)
      tcm_t          (il)    = tcm_b           (ifull)
      runoff_s_t     (il)    = runoff_s_b      (ifull)
      runoff_g_t     (il)    = runoff_g_b      (ifull)
      snow_melt_t    (il)    = snow_melt_b     (ifull)

      prr_con_t      (il)    = prr_con_b       (ifull)
      prs_con_t      (il)    = prs_con_b       (ifull)
      conv_frac_t    (il)    = conv_frac_b     (ifull)
      prr_gsp_t      (il)    = prr_gsp_b       (ifull)
      prs_gsp_t      (il)    = prs_gsp_b       (ifull)
      prg_gsp_t      (il)    = prg_gsp_b       (ifull)
#ifdef TWOMOM_SB
      prh_gsp_t      (il)    = prh_gsp_b       (ifull)
#endif
      tfv_t          (il)    = tfv_b           (ifull)
      sobs_t         (il)    = sobs_b          (ifull)
      thbs_t         (il)    = thbs_b          (ifull)
      pabs_t         (il)    = pabs_b          (ifull)

      u_t            (il)    = u_m_b           (ifull,ke)
      v_t            (il)    = v_m_b           (ifull,ke)
      t_t            (il)    = t_b             (ifull,ke)
      qv_t           (il)    = qv_b            (ifull,ke)
      ptot_t         (il)    = ptot_b          (ifull,ke)
    ENDDO
    !$acc end parallel

    !$acc parallel
    DO k = 1, ke_soil
      !$acc loop gang vector private (ifull)
      DO il = 1, n_landgp(ib)
        ! grid point index in the full variable
        ifull = mind_landgp(il,ib)

        lhfl_pl_t          (il,k)  = lhfl_pl_b           (ifull,k)
      ENDDO
    ENDDO
    !$acc end parallel

    !$acc parallel
    DO k = 0, ke_soil+1
      !$acc loop gang vector private (ifull)
      DO il = 1, n_landgp(ib)
        ! grid point index in the full variable
        ifull = mind_landgp(il,ib)

        t_so_now_t         (il,k)  = t_so_b              (ifull,k)   ! tl
      ENDDO
    ENDDO
    !$acc end parallel

!VS <
 if(lsnow) then
    ! 3D fields
    !$acc parallel async
    DO k = 1,ke_snow
      !$acc loop gang vector private (ifull)
      DO il = 1, n_landgp(ib)
        ! grid point index in the full variable
        ifull = mind_landgp(il,ib)

        ! Use input field in working precision
        t_sn_now_t    (il,k) = t_sn_b    (ifull,k)
        theta_i_now_t (il,k) = theta_i_b (ifull,k)
        theta_w_now_t (il,k) = theta_w_b (ifull,k)
        theta_a_now_t (il,k) = theta_a_b (ifull,k)
        dzm_sn_now_t  (il,k) = dzm_sn_b  (ifull,k)


      ENDDO
    ENDDO
    !$acc end parallel

    ! 2D fields
    !$acc parallel async
    !$acc loop gang vector private (ifull)
    DO il = 1, n_landgp(ib)
        ! grid point index in the full variable
        ifull = mind_landgp(il,ib)

        top_sn_now_t  (il  ) =  top_sn_b  (ifull  )
        hn_sn_now_t   (il  ) =  hn_sn_b   (ifull  )

!        swdir_s_t     (il  ) =  swdir_s_b  (ifull )
!        swdifd_s_t    (il  ) =  swdifd_s_b (ifull )
!        swdifu_s_t    (il  ) =  swdifu_s_b (ifull )
!        lwd_s_t       (il  ) =  lwd_s_b    (ifull )
!        lwu_s_t       (il  ) =  lwu_s_b    (ifull )


    ENDDO
    !$acc end parallel

  ENDIF
!VS >

    !$acc parallel
    DO k = 1, ke_soil+1
      !$acc loop gang vector private (ifull)
      DO il = 1, n_landgp(ib)
        ! grid point index in the full variable
        ifull = mind_landgp(il,ib)

        w_so_now_t         (il,k)  = w_so_b              (ifull,k)   ! tl
        w_so_ice_now_t     (il,k)  = w_so_ice_b          (ifull,k)   ! tl
      ENDDO
    ENDDO
    !$acc end parallel

!    endif

    IF (lmulti_snow) THEN
      !$acc parallel
      DO k = 0, ke_snow
        !$acc loop gang vector private (ifull)
        DO il = 1, n_landgp(ib)
          ! grid point index in the full variable
          ifull = mind_landgp(il,ib)

          t_snow_mult_now_t  (il,k)  = t_snow_mult_b       (ifull,k)   ! tl
        ENDDO
      ENDDO
      !$acc end parallel

      !$acc parallel
      DO k = 1, ke_snow
        !$acc loop gang vector private (ifull)
        DO il = 1, n_landgp(ib)
          ! grid point index in the full variable
          ifull = mind_landgp(il,ib)

          rho_snow_mult_now_t(il,k)  = rho_snow_mult_b     (ifull,k)   ! tl
          wliq_snow_now_t    (il,k)  = wliq_snow_b         (ifull,k)   ! tl
          wtot_snow_now_t    (il,k)  = w_snow_mult_b       (ifull,k)   ! tl
          dzh_snow_now_t     (il,k)  = dzh_snow_mult_b     (ifull,k)   ! tl
        ENDDO
      ENDDO
      !$acc end parallel
    ENDIF

    ! tsnred is a dummy for COSMO, which has to be set to 0.0
    !$acc kernels
    tsnred_t(:) = 0.0_wp
    !$acc end kernels

    ! Do initialization for TERRA in the first time step
    ! (should be called before the time loop, but what about the blocked data structure?)
    IF ( (ntstep == 0 ) .OR. lsoilinit_dfi) THEN
      lsoilinit_dfi = .FALSE.
      CALL terra_init(                                &
       init_mode         = 1                        , & 
                ! 1: full coldstart is executed
                ! 2: warmstart with full fields for h_snow from snow analysis
                ! 3: warmstart (within assimilation cycle) with analysis increments for h_snow
       nvec              = nproma                       , & !IN array dimensions                     ! nvec
       ivstart           = 1                            , & !IN optional start index                 ! ivstart
       ivend             = n_landgp(ib)                 , & !IN optional end   index                 ! ivend
       iblock            = ib                           , & !IN number of block
       ke_soil           = ke_soil                      , & !IN without lowermost (climat.) soil layer
       ke_snow           = ke_snow                      , & !IN 
       zmls              = czmls              (:)       , & !IN processing soil level structure 
       soiltyp_subs      = isoiltyp_t         (:)       , & !IN type of the soil (keys 0-9)         --    
       rootdp            = rootdp_t           (:)       , & !IN depth of the roots                ( m  )
       plcov             = plcov_t            (:)       , & !IN fraction of plant cover                  --
       t_snow_now        = t_snow_now_t       (:)       , & !INOUT temperature of the snow-surface (  K  )
       t_snow_mult_now   = t_snow_mult_now_t  (:,:)     , & !INOUT temperature of the snow-surface (  K  )
       t_s_now           = t_s_now_t          (:)       , & !INOUT temperature of the ground surface (  K  )
       t_s_new           = t_s_new_t          (:)       , & !OUT temperature of the ground surface   (  K  )
       t_sk_now          = t_sk_now_t         (:)       , & !INOUT skin temperature                  (  K  )
       t_sk_new          = t_sk_new_t         (:)       , & !OUT skin temperature                    (  K  )
       w_snow_now        = w_snow_now_t       (:)       , & !INOUT water content of snow         (m H2O) 
       h_snow            = h_snow_t           (:)       , & !INOUT snow height
       rho_snow_now      = rho_snow_now_t     (:)       , & !IN  snow density                    (kg/m**3)
       rho_snow_mult_now = rho_snow_mult_now_t(:,:)     , & !INOUT snow density               (kg/m**3) 
       t_so_now          = t_so_now_t         (:,:)     , & !INOUT soil temperature (main level)    (  K  )
       t_so_new          = t_so_new_t         (:,:)     , & !OUT soil temperature (main level)      (  K  )
       w_so_now          = w_so_now_t         (:,:)     , & !IN  total water content (ice + liquid water) (m H20)
       w_so_new          = w_so_new_t         (:,:)     , & !OUT total water content (ice + liquid water) (m H20)
       w_so_ice_now      = w_so_ice_now_t     (:,:)     , & !IN  ice content   (m H20)
       w_so_ice_new      = w_so_ice_new_t     (:,:)     , & !OUT ice content   (m H20)
       wliq_snow_now     = wliq_snow_now_t    (:,:)     , & !INOUT liquid water content in the snow     (m H2O)
       wtot_snow_now     = wtot_snow_now_t    (:,:)     , & !INOUT total (liquid+solid) water content of snow  (m H2O)
       dzh_snow_now      = dzh_snow_now_t     (:,:)    )    !INOUT layer thickness between half levels in snow (  m  )

      ! and call diagnosis for snow fraction and computation of t_g_now
      IF (itype_canopy == 1) THEN
        CALL diag_snowfrac_tg(                     &
            ivstart   = 1                        , &
            ivend     = n_landgp(ib)             , & ! start/end indices
            t_snow    = t_snow_now_t  (:)        , & ! snow temp
            t_soiltop = t_s_now_t     (:)        , & ! soil top temp
            w_snow    = w_snow_now_t  (:)        , & ! snow WE
            rho_snow  = rho_snow_now_t(:)        , & ! snow density
            freshsnow = freshsnow_t   (:)        , & ! fresh snow fraction
            sso_sigma = sso_sigma_t              , & ! sso stdev
            tai       = tai_t         (:)        , & ! effective leaf area index
            snowfrac  = zf_snow_t                , & ! OUT: snow cover fraction
            t_g       = t_g_now_t     (:)        )   ! OUT: averaged ground temp
      ELSE IF (itype_canopy == 2) THEN
        CALL diag_snowfrac_tg(                     &
            ivstart   = 1                        , &
            ivend     = n_landgp(ib)             , & ! start/end indices
            t_snow    = t_snow_now_t  (:)        , & ! snow temp
            t_soiltop = t_sk_now_t    (:)        , & ! skin temperature
            w_snow    = w_snow_now_t  (:)        , & ! snow WE
            rho_snow  = rho_snow_now_t(:)        , & ! snow density
            freshsnow = freshsnow_t   (:)        , & ! fresh snow fraction
            sso_sigma = sso_sigma_t              , & ! sso stdev
            tai       = tai_t         (:)        , & ! effective leaf area index
            snowfrac  = zf_snow_t                , & ! OUT: snow cover fraction
            t_g       = t_g_now_t     (:)        )   ! OUT: averaged ground temp
      END IF
    ENDIF

!VS <

IF(lsnow) THEN

   ! Do initializations for the multi-layer snow cover schmeme (SNOWPOLINO)
   IF( (ntstep == 0) ) THEN

     CALL snow_init(                                      &
       ! general
       nvec              = nproma                       , & !IN array dimensions nvec
       ivstart           = 1                            , & !IN optional start index
       ivend             = n_landgp(ib)                 , & !IN optional end index
       ke_soil           = ke_soil                      , & !IN without lowermost (climat.) soil layer      US ??
       dt                = zdt                          , & !IN integration timestep

       ! snow
       t_sn_now          = t_sn_now_t         (:,:)     , & !INOUT snow temperature (main level)             (  K  )
       t_sn_new          = t_sn_new_t         (:,:)     , & !OUT snow temperature (main level)               (  K  )

       theta_i_now       = theta_i_now_t      (:,:)     , & !INOUT volumetric ice content                    (  -  )
       theta_i_new       = theta_i_new_t      (:,:)     , & !OUT   volumetric ice content                    (  -  )

       theta_w_now       = theta_w_now_t      (:,:)     , & !INOUT volumetric water content                  (  -  )
       theta_w_new       = theta_w_new_t      (:,:)     , & !OUT   volumetric water content                  (  -  )

       theta_a_now       = theta_a_now_t      (:,:)     , & !INOUT volumetric air content                    (  -  )
       theta_a_new       = theta_a_new_t      (:,:)     , & !OUT   volumetric air content                    (  -  )

       dzm_sn_now        = dzm_sn_now_t       (:,:)     , & !INOUT snow layer thickness                      (  -  )
       dzm_sn_new        = dzm_sn_new_t       (:,:)     , & !OUT   snow layer thickness                      (  -  )

       hn_sn_now         = hn_sn_now_t        (:)       , & !INOUT new snow amounts (storage)                (  m  )
       hn_sn_new         = hn_sn_new_t        (:)       , & !OUT   new snow amounts (storage)                (  m  )

       top_sn_now        = top_sn_now_t       (:)       , & !INOUT index of the first (top) snow layer       (  -  )
       top_sn_new        = top_sn_new_t       (:)       , & !OUT   index of the first (top) snow layer       (  -  )

       h_snow            = h_snow_t           (:)       , & !INOUT snow height

       ! soil
       t_so_now          = t_so_now_t         (:,:)     , & !INOUT soil temperature (main level)             (  K  )
       t_so_new          = t_so_new_t         (:,:)       ) !OUT soil temperature (main level)               (  K  )

   ENDIF ! ntstep
ENDIF ! lsnow

!VS>



    ! now call the TERRA scheme
    CALL terra (                                          &
       nvec              = nproma                       , & !IN array dimensions                       ! nvec
       ivstart           = 1                            , & !IN optional start index                   ! ivstart
       ivend             = n_landgp(ib)                 , & !IN optional end   index                   ! ivend
       iblock            = ib                           , & !IN number of block
       ke_soil           = ke_soil                      , & !IN without lowermost (climat.) soil layer      US ??
       ke_snow           = ke_snow                      , & !IN 
!!DL   ke_soil_hy        = ibot_w_so                    , & !IN number of hydrological active soil layers
       ke_soil_hy        = ke_soil                      , & !IN number of hydrological active soil layers
       zmls              = czmls  (:)                   , & !IN processing soil level structure 
       icant             = itype_tran                   , & !IN canopy-type
       nclass_gscp       = n_gscpclass                  , & !IN number of hydrometeor classes          ! ?????
       dt                = zdt                          , & !IN 
       soiltyp_subs      = isoiltyp_t         (:)       , & !IN type of the soil (keys 0-9)              --    
       plcov             = plcov_t            (:)       , & !IN fraction of plant cover                  --
       rootdp            = rootdp_t           (:)       , & !IN depth of the roots                     ( m  )
       sai               = sai_t              (:)       , & !IN surface area index                       --
       tai               = tai_t              (:)       , & !IN surface area index                       --
       eai               = eai_t              (:)       , & !IN surface area index                       --
       skinc             = skinc_t            (:)       , & !IN skin conductivity                      (W/m**2/K)
       rsmin2d           = rsmin2d_t          (:)       , & !IN minimum stomata resistance             ( s/m )
  ! for TERRA_URB
!          fr_paved          = fr_paved_t         (:)       , & !IN fraction of paved area
!          sa_uf             = sa_uf_t            (:)       , & !IN total impervious surface-area index
!          ai_uf             = ai_uf_t            (:)       , & !IN surface area index of the urban fabric
!          alb_red_uf        = alb_red_uf_t       (:)       , & !IN albedo reduction factor for the urban fabric

       u                 =  u_t               (:)       , & !IN zonal wind speed
       v                 =  v_t               (:)       , & !IN meridional wind speed 
       t                 =  t_t               (:)       , & !IN temperature                            (  K  )
       qv                =  qv_t              (:)       , & !IN specific water vapor content           (kg/kg)
       ptot              =  ptot_t            (:)       , & !IN base state pressure                    ( Pa  ) 
       ps                =  ps_t              (:)       , & !IN surface pressure                       ( Pa  )

       t_snow_now        = t_snow_now_t       (:)       , & !INOUT temperature of the snow-surface     (  K  )
       t_snow_new        = t_snow_new_t       (:)       , & !OUT temperature of the snow-surface       (  K  )

       t_snow_mult_now   = t_snow_mult_now_t  (:,:)     , & !INOUT temperature of the snow-surface     (  K  )
       t_snow_mult_new   = t_snow_mult_new_t  (:,:)     , & !OUT temperature of the snow-surface       (  K  )

       t_s_now           = t_s_now_t          (:)       , & !INOUT temperature of the ground surface   (  K  )
       t_s_new           = t_s_new_t          (:)       , & !OUT temperature of the ground surface     (  K  )

       t_sk_now          = t_sk_now_t         (:)       , & !INOUT skin temperature                    (  K  )
       t_sk_new          = t_sk_new_t         (:)       , & !OUT skin temperature                      (  K  )

       t_g               = t_g_now_t          (:)       , & !INOUT weighted surface temperature        (  K  )
       qv_s              = qv_s_t             (:)       , & !INOUT specific humidity at the surface    (kg/kg)

       w_snow_now        = w_snow_now_t       (:)       , & !INOUT water content of snow               (m H2O) 
       w_snow_new        = w_snow_new_t       (:)       , & !OUT water content of snow                 (m H2O) 

       rho_snow_now      = rho_snow_now_t     (:)       , & !IN  snow density                          (kg/m**3)
       rho_snow_new      = rho_snow_new_t     (:)       , & !OUT snow density                          (kg/m**3)

       rho_snow_mult_now = rho_snow_mult_now_t(:,:)     , & !INOUT snow density                        (kg/m**3) 
       rho_snow_mult_new = rho_snow_mult_new_t(:,:)     , & !OUT snow density                          (kg/m**3) 

       h_snow            = h_snow_t           (:)       , & !INOUT snow height
       h_snow_gp         = h_snow_gp_t        (:)       , & !IN grid-point averaged snow height
       meltrate          = meltrate_t         (:)       , & !OUT snow melting rate
       tsnred            = tsnred_t           (:)       , & !dummy for COSMO

       w_i_now           = w_i_now_t          (:)       , & !INOUT water content of interception water (m H2O)
       w_i_new           = w_i_new_t          (:)       , & !OUT water content of interception water   (m H2O)

       w_p_now           = w_p_now_t          (:)       , & !INOUT water content of interception water (m H2O)
       w_p_new           = w_p_new_t          (:)       , & !OUT water content of interception water   (m H2O)

       w_s_now           = w_s_now_t          (:)       , & !INOUT water content of interception water (m H2O)
       w_s_new           = w_s_new_t          (:)       , & !OUT water content of interception water   (m H2O)

       t_so_now          = t_so_now_t         (:,:)     , & !INOUT soil temperature (main level)       (  K  )
       t_so_new          = t_so_new_t         (:,:)     , & !OUT soil temperature (main level)         (  K  )

       w_so_now          = w_so_now_t         (:,:)     , & !IN  total water content (ice+liquid)      (m H20)
       w_so_new          = w_so_new_t         (:,:)     , & !OUT total water content (ice+liquid)      (m H20)

       w_so_ice_now      = w_so_ice_now_t     (:,:)     , & !IN  ice content                           (m H20)
       w_so_ice_new      = w_so_ice_new_t     (:,:)     , & !OUT ice content                           (m H20)

       u_10m             = u_10m_t            (:)       , & !IN zonal wind in 10m                      ( m/s )
       v_10m             = v_10m_t            (:)       , & !IN meridional wind in 10m                 ( m/s )
       freshsnow         = freshsnow_t        (:)       , & !INOUT indicator for age of snow in top of snow layer
       zf_snow           = zf_snow_t          (:)       , & !INOUT snow-cover fraction

       wliq_snow_now     = wliq_snow_now_t    (:,:)     , & !INOUT liquid water content in the snow    (m H2O)
       wliq_snow_new     = wliq_snow_new_t    (:,:)     , & !OUT liquid water content in the snow      (m H2O)

       wtot_snow_now     = wtot_snow_now_t    (:,:)     , & !INOUT total (liquid+solid) water content of snow  (m H2O)
       wtot_snow_new     = wtot_snow_new_t    (:,:)     , & !OUT   total (liquid+solid) water content of snow  (m H2O)

       dzh_snow_now      = dzh_snow_now_t     (:,:)     , & !INOUT layer thickness between half levels in snow (  m  )
       dzh_snow_new      = dzh_snow_new_t     (:,:)     , & !OUT layer thickness between half levels in snow   (  m  )

       prr_con           = prr_con_t    (:)       , & !IN precipitation rate of rain, convective       (kg/m2*s)
       prs_con           = prs_con_t    (:)       , & !IN precipitation rate of snow, convective       (kg/m2*s)
       conv_frac         = conv_frac_t  (:)       , & !IN convective area fraction
       prr_gsp           = prr_gsp_t    (:)       , & !IN precipitation rate of rain, grid-scale       (kg/m2*s)
       prs_gsp           = prs_gsp_t    (:)       , & !IN precipitation rate of snow, grid-scale       (kg/m2*s)
       prg_gsp           = prg_gsp_t    (:)       , & !IN precipitation rate of graupel, grid-scale    (kg/m2*s)
#ifdef TWOMOM_SB
       prh_gsp           = prh_gsp_t    (:)       , & !IN precipitation rate of hail, grid-scale       (kg/m2*s)
#endif

       tch               = tch_t        (:)       , & !INOUT turbulent transfer coefficient for heat     ( -- )
       tcm               = tcm_t        (:)       , & !INOUT turbulent transfer coefficient for momentum ( -- )
       tfv               = tfv_t        (:)       , & !INOUT laminar reduction factor for evaporation    ( -- )

       sobs              = sobs_t       (:)       , & !IN solar radiation at the ground               (W/m2)
       thbs              = thbs_t       (:)       , & !IN thermal radiation at the ground             (W/m2)
       pabs              = pabs_t       (:)       , & !IN photosynthetic active radiation             (W/m2)

       runoff_s          = runoff_s_t   (:)       , & !INOUT surface water runoff; sum over forecast  (kg/m2)
       runoff_g          = runoff_g_t   (:)       , & !INOUT soil water runoff; sum over forecast     (kg/m2)
! for TERRA_URB
!      w_imp             = w_imp_t      (:)       , & !INOUT impervious water storage                  --
!      w_isa             = w_isa_t      (:)       , & !INOUT same, multiplied by fr_paved              --

       zshfl_s           = dum1fl_s_t   (:)       , & !OUT sensible heat flux soil/air interface    (W/m2) 
       zlhfl_s           = dum2fl_s_t   (:)       , & !OUT latent   heat flux soil/air interface    (W/m2) 
       zshfl_snow        = shfl_snow_t  (:)       , & !OUT sensible heat flux snow/air interface    (W/m2) 
       zlhfl_snow        = lhfl_snow_t  (:)       , & !OUT latent   heat flux snow/air interface    (W/m2) 
       lhfl_bs           = lhfl_bs_t    (:)       , & !OUT latent heat flux from bare soil evap.    (W/m2)
       lhfl_pl           = lhfl_pl_t    (:,:)     , & !OUT latent heat flux from plants    evap.    (W/m2)
       rstom             = rstom_t      (:)       , & !OUT stomatal resistance                      ( s/m )
! is optional
       zshfl_sfc         = shfl_s_t               , & !OUT sensible heat flux surface interface     (W/m2) 
       zlhfl_sfc         = lhfl_s_t               , & !OUT latent   heat flux surface interface     (W/m2) 
       zqhfl_sfc         = qvfl_s_t            )      !OUT moisture heat flux surface interface     (kg/m2/s) 
       !  this one must be taken by TERRA_URB

    ! and call diagnosis for snow fraction and computation of t_g_new
    IF (itype_canopy == 1) THEN
      CALL diag_snowfrac_tg(                     &
          ivstart   = 1                        , &
          ivend     = n_landgp(ib)             , & ! start/end indices
          t_snow    = t_snow_new_t  (:)        , & ! snow temp
          t_soiltop = t_s_new_t     (:)        , & ! soil top temp
          w_snow    = w_snow_new_t  (:)        , & ! snow WE
          rho_snow  = rho_snow_new_t(:)        , & ! snow density
          freshsnow = freshsnow_t   (:)        , & ! fresh snow fraction
          meltrate  = meltrate_t               , & ! snow melting rate
          sso_sigma = sso_sigma_t              , & ! sso stdev
          tai       = tai_t         (:)        , & ! effective leaf area index
          snowfrac  = zf_snow_t                , & ! OUT: snow cover fraction
          t_g       = t_g_new_t     (:)        )   ! OUT: averaged ground temp
    ELSE IF (itype_canopy == 2) THEN
      CALL diag_snowfrac_tg(                     &
          ivstart   = 1                        , &
          ivend     = n_landgp(ib)             , & ! start/end indices
          t_snow    = t_snow_new_t  (:)        , & ! snow temp
          t_soiltop = t_sk_new_t    (:)        , & ! skin temperature
          w_snow    = w_snow_new_t  (:)        , & ! snow WE
          rho_snow  = rho_snow_new_t(:)        , & ! snow density
          freshsnow = freshsnow_t   (:)        , & ! fresh snow fraction
          meltrate  = meltrate_t               , & ! snow melting rate
          sso_sigma = sso_sigma_t              , & ! sso stdev
          tai       = tai_t         (:)        , & ! effective leaf area index
          snowfrac  = zf_snow_t                , & ! OUT: snow cover fraction
          t_g       = t_g_new_t     (:)        )   ! OUT: averaged ground temp
    END IF

    ! compute snow melt amount:
    !$acc parallel
    !$acc loop gang vector
    DO il = 1, n_landgp(ib)
      snow_melt_t(il) = snow_melt_t(il) + meltrate_t(il) * zdt
    ENDDO
    !$acc end parallel

    ! Copy the land points back to the blocked format

    !$acc parallel
    !$acc loop gang vector private (ifull)
    DO il = 1, n_landgp(ib)
      ! grid point index in the full variable
      ifull = mind_landgp(il,ib)

      t_snow_b        (ifull)  = t_snow_now_t   (il)
      t_snow_new_b    (ifull)  = t_snow_new_t   (il)
      t_s_b           (ifull)  = t_s_now_t      (il)
      t_s_new_b       (ifull)  = t_s_new_t      (il)
      t_sk_b          (ifull)  = t_sk_now_t     (il)
      t_sk_new_b      (ifull)  = t_sk_new_t     (il)
      t_g_new_b       (ifull)  = t_g_new_t      (il)
      qv_s_b          (ifull)  = qv_s_t         (il)
      w_snow_b        (ifull)  = w_snow_now_t   (il)
      w_snow_new_b    (ifull)  = w_snow_new_t   (il)
      rho_snow_b      (ifull)  = rho_snow_now_t (il)
      rho_snow_new_b  (ifull)  = rho_snow_new_t (il)
      h_snow_b        (ifull)  = h_snow_t       (il)
      h_snow_new_b    (ifull)  = h_snow_t       (il)   ! only h_snow is interface to terra???
      fr_snow_b       (ifull)  = zf_snow_t      (il)    
      freshsnow_b     (ifull)  = freshsnow_t    (il)
      meltrate_b      (ifull)  = meltrate_t     (il)
      w_i_b           (ifull)  = w_i_now_t      (il)
      w_i_new_b       (ifull)  = w_i_new_t      (il)
      w_p_b           (ifull)  = w_p_now_t      (il)    
      w_p_new_b       (ifull)  = w_p_new_t      (il)    
      w_s_b           (ifull)  = w_s_now_t      (il)    
      w_s_new_b       (ifull)  = w_s_new_t      (il)    
      tfv_b           (ifull)  = tfv_t          (il)    
      shfl_snow_b     (ifull)  = shfl_snow_t    (il)    
      lhfl_snow_b     (ifull)  = lhfl_snow_t    (il)    
      lhfl_bs_b       (ifull)  = lhfl_bs_t      (il)    
      rstom_b         (ifull)  = rstom_t        (il)    
      snow_melt_b     (ifull)  = snow_melt_t    (il)    

      tch_b           (ifull)  = tch_t          (il)
      tcm_b           (ifull)  = tcm_t          (il)
      runoff_s_b      (ifull)  = runoff_s_t     (il)
      runoff_g_b      (ifull)  = runoff_g_t     (il)
      shfl_s_b        (ifull)  = shfl_s_t       (il)
      lhfl_s_b        (ifull)  = lhfl_s_t       (il)
      qvfl_s_b        (ifull)  = qvfl_s_t       (il)
    ENDDO
    !$acc end parallel

    !$acc parallel
    DO k = 1, ke_soil
      !$acc loop gang vector private (ifull)
      DO il = 1, n_landgp(ib)
        ! grid point index in the full variable
        ifull = mind_landgp(il,ib)

        lhfl_pl_b     (ifull,k)  = lhfl_pl_t    (il,k)  
      ENDDO
    ENDDO
    !$acc end parallel

    !$acc parallel
    DO k = 0, ke_soil+1
      !$acc loop gang vector private (ifull)
      DO il = 1, n_landgp(ib)
        ! grid point index in the full variable
        ifull = mind_landgp(il,ib)

        t_so_b              (ifull,k)  = t_so_now_t         (il,k)
        t_so_new_b          (ifull,k)  = t_so_new_t         (il,k)
      ENDDO
    ENDDO
    !$acc end parallel


!VS <
IF(lsnow) THEN
   ! now call the multi-layer snow cover scheme (SNOWPOLINO)
   CALL snowpolino(                                      &
       ! general
       nvec              = nproma                       , & !IN array dimensions nvec
       ivstart           = 1                            , & !IN optional start index
       ivend             = n_landgp(ib)                 , & !IN optional end index
       iblock            = ib                           , & !IN number of block
       ke_soil           = ke_soil                      , & !IN without lowermost (climat.) soil layer      US ??
       dt                = zdt                          , & !IN integration timestep
       nclass_gscp       = n_gscpclass                  , & !IN number of hydrometeor classes
      ! meteorological
       u                 =  u_t               (:)       , & !IN zonal wind speed
       v                 =  v_t               (:)       , & !IN meridional wind speed
       t                 =  t_t               (:)       , & !IN temperature                                  (    K)
       qv                =  qv_t              (:)       , & !IN specific water vapor content                 (kg/kg)
       ps                =  ps_t              (:)       , & !IN surface pressure                             ( Pa  )
       prr_con           = prr_con_t          (:)       , & !IN precipitation rate ofrain, convective        (kg/m2*s)
       prs_con           = prs_con_t          (:)       , & !IN precipitation rate of snow, convective       (kg/m2*s)
       conv_frac         = conv_frac_t        (:)       , & !IN convective area fraction
       prr_gsp           = prr_gsp_t          (:)       , & !IN precipitation rate of rain, grid-scale       (kg/m2*s)
       prs_gsp           = prs_gsp_t          (:)       , & !IN precipitation rate of snow, grid-scale       (kg/m2*s)
       prg_gsp           = prg_gsp_t          (:)       , & !IN precipitation rate of graupel, grid-scale    (kg/m2*s)
#ifdef TWOMOM_SB
       prh_gsp           = prh_gsp_t          (:)       , & !IN precipitation rate of hail, grid-scale       (kg/m2*s)
#endif
       sobs              = sobs_t       (:)       , & !IN solar radiation at the ground               (W/m2)
       thbs              = thbs_t       (:)       , & !IN thermal radiation at the ground             (W/m2)

!       swdir_s           = swdir_s_t          (:)       , & !
!       swdifd_s          = swdifd_s_t         (:)       , & !
!       swdifu_s          = swdifu_s_t         (:)       , & !
!       lwd_s             = lwd_s_t            (:)       , & !
!       lwu_s             = lwu_s_t            (:)       , & !
       ! surface
       t_snow_now        = t_snow_now_t       (:)       , & !INOUT temperature of the snow-surface           (  K  )
       t_snow_new        = t_snow_new_t       (:)       , & !OUT temperature of the snow-surface             (  K  )
       zshfl_snow        = shfl_snow_t        (:)       , & !OUT sensible heat flux snow/air interface       ( W/m2)
       zlhfl_snow        = lhfl_snow_t        (:)       , & !OUT latent   heat flux snow/air interface       ( W/m2)

       ! snow
       t_sn_now          = t_sn_now_t         (:,:)     , & !INOUT snow temperature (main level)             (  K  )
       t_sn_new          = t_sn_new_t         (:,:)     , & !OUT snow temperature (main level)               (  K  )

       theta_i_now       = theta_i_now_t      (:,:)     , & !INOUT volumetric ice content                    (  -  )
       theta_i_new       = theta_i_new_t      (:,:)     , & !OUT   volumetric ice content                    (  -  )

       theta_w_now       = theta_w_now_t      (:,:)     , & !INOUT volumetric water content                  (  -  )
       theta_w_new       = theta_w_new_t      (:,:)     , & !OUT   volumetric water content                  (  -  )

       theta_a_now       = theta_a_now_t      (:,:)     , & !INOUT volumetric air content                    (  -  )
       theta_a_new       = theta_a_new_t      (:,:)     , & !OUT   volumetric air content                    (  -  )

       dzm_sn_now        = dzm_sn_now_t       (:,:)     , & !INOUT snow layer thickness                      (  -  )
       dzm_sn_new        = dzm_sn_new_t       (:,:)     , & !OUT   snow layer thickness                      (  -  )

       hn_sn_now         = hn_sn_now_t        (:)       , & !INOUT new snow amounts (storage)                (  m  )
       hn_sn_new         = hn_sn_new_t        (:)       , & !OUT   new snow amounts (storage)                (  m  )

       top_sn_now        = top_sn_now_t       (:)       , & !INOUT index of the first (top) snow layer       (  -  )
       top_sn_new        = top_sn_new_t       (:)       , & !OUT   index of the first (top) snow layer       (  -  )

       h_snow            = h_snow_t           (:)       , & !INOUT snow height

       w_snow_now        = w_snow_now_t       (:)       , & !INOUT water content of snow                     (m H2O)
       w_snow_new        = w_snow_new_t       (:)       , & !OUT water content  of snow                      (m H2O)

       rho_snow_new      = rho_snow_new_t     (:)       , & !INOUT density of snow

       ! soil
       t_so_now          = t_so_now_t         (:,:)     , & !INOUT soil temperature (main level)             (  K  )
       t_so_new          = t_so_new_t         (:,:)     , & !OUT soil temperature (main level)               (  K  )

       w_so_now          = w_so_now_t         (:,:)     , & !IN  total water      content (ice+liquid)      (m H20)
       w_so_new          = w_so_new_t         (:,:)     , & !OUT total water      content (ice+liquid)      (m H20)

       w_so_ice_now      = w_so_ice_now_t     (:,:)     , & !IN  ice content                                 (m H20)
       w_so_ice_new      = w_so_ice_new_t     (:,:)     , & !OUT ice content                                 (m H20)

       zmls              = czmls              (:)       , & !IN processing soil level structure

       plcov             = plcov_t            (:)       , & !IN fraction of plant cover                        --
       rootdp            = rootdp_t           (:)       , & !IN depth of the roots                           ( m  )

       soiltyp_subs      = isoiltyp_t         (:)         ) !IN type of the soil       (keys 0-9)              --

    !write(*,*) 'HOLY SNOW: ', h_snow_t 

    ! and call diagnosis for snow fraction and computation of t_g_new
    IF (itype_canopy == 1) THEN
      CALL diag_snowfrac_tg(                     &
          ivstart   = 1                        , &
          ivend     = n_landgp(ib)             , & ! start/end indices
          t_snow    = t_snow_new_t  (:)        , & ! snow temp
          t_soiltop = t_s_new_t     (:)        , & ! soil top temp
          w_snow    = w_snow_new_t  (:)        , & ! snow WE
          rho_snow  = rho_snow_new_t(:)        , & ! snow density
          freshsnow = freshsnow_t   (:)        , & ! fresh snow fraction
          meltrate  = meltrate_t               , & ! snow melting rate
          sso_sigma = sso_sigma_t              , & ! sso stdev
          tai       = tai_t         (:)        , & ! effective leaf area index
          snowfrac  = zf_snow_t                , & ! OUT: snow cover fraction
          t_g       = t_g_new_t     (:)        )   ! OUT: averaged ground temp
    ELSE IF (itype_canopy == 2) THEN
      CALL diag_snowfrac_tg(                     &
          ivstart   = 1                        , &
          ivend     = n_landgp(ib)             , & ! start/end indices
          t_snow    = t_snow_new_t  (:)        , & ! snow temp
          t_soiltop = t_sk_new_t    (:)        , & ! skin temperature
          w_snow    = w_snow_new_t  (:)        , & ! snow WE
          rho_snow  = rho_snow_new_t(:)        , & ! snow density
          freshsnow = freshsnow_t   (:)        , & ! fresh snow fraction
          meltrate  = meltrate_t               , & ! snow melting rate
          sso_sigma = sso_sigma_t              , & ! sso stdev
          tai       = tai_t         (:)        , & ! effective leaf area index
          snowfrac  = zf_snow_t                , & ! OUT: snow cover fraction
          t_g       = t_g_new_t     (:)        )   ! OUT: averaged ground temp
    END IF

    ! Update fractional snow cover; Note idiag_snowfrac=1 only
    ! --------------------

    !$acc parallel async
    !$acc loop gang vector
    DO il = 1, n_landgp(ib)

      zf_snow_t(il) = MIN(1.0_vpp, w_snow_new_t(il)/cf_snow) 

    ENDDO
    !$acc end parallel


    ! Update surface fluxes
    ! ---------------------

    !$acc parallel async
    !$acc loop gang vector
    DO il = 1, n_landgp(ib)

      shfl_s_t(il) = dum1fl_s_t(il)*(1.0_vpp - zf_snow_t(il)) + shfl_snow_t(il)*zf_snow_t(il) 
      lhfl_s_t(il) = dum2fl_s_t(il)*(1.0_vpp - zf_snow_t(il)) + lhfl_snow_t(il)*zf_snow_t(il)
     
    ENDDO
    !$acc end parallel
   

   ! Copy back to block structure
   ! -----------------------

   !2-D fields
  
    !$acc parallel async
    !$acc loop gang vector private (ifull)
    DO il = 1, n_landgp(ib)
      ! grid point index in the full variable
      ifull = mind_landgp(il,ib)

      shfl_s_b        (ifull)  = shfl_s_t       (il)
      lhfl_s_b        (ifull)  = lhfl_s_t       (il)
      h_snow_b        (ifull)  = h_snow_t       (il)
      h_snow_new_b    (ifull)  = h_snow_t       (il)

    ENDDO
    !$acc end parallel


   ! 3D fields
    !$acc parallel async
    DO k = 1, ke_snow
      !$acc loop gang vector private (ifull)
      DO il = 1, n_landgp(ib)
        ! grid point index in the full variable
        ifull = mind_landgp(il,ib)

        t_sn_b     (ifull,k)       = t_sn_now_t (il,k)
        t_sn_new_b (ifull,k)       = t_sn_new_t (il,k)

        theta_i_b     (ifull,k) = theta_i_now_t (il,k)
        theta_i_new_b (ifull,k) = theta_i_new_t (il,k)

        theta_w_b     (ifull,k) = theta_w_now_t (il,k)
        theta_w_new_b (ifull,k) = theta_w_new_t (il,k)

        theta_a_b     (ifull,k) = theta_a_now_t (il,k)
        theta_a_new_b (ifull,k) = theta_a_new_t (il,k)

        dzm_sn_b      (ifull,k)  = dzm_sn_now_t (il,k)
        dzm_sn_new_b  (ifull,k)  = dzm_sn_new_t (il,k)

      ENDDO
    ENDDO
    !$acc end parallel

   ! more 2-D fields

    !$acc parallel async
    !$acc loop gang vector private (ifull)
    DO il = 1, n_landgp(ib)
      ! grid point index in the full variable
      ifull = mind_landgp(il,ib)

        top_sn_b      (ifull  ) = top_sn_now_t (il  )
        top_sn_new_b  (ifull  ) = top_sn_new_t (il  )

        hn_sn_b      (ifull  )  = hn_sn_now_t (il  )
        hn_sn_new_b  (ifull  )  = hn_sn_new_t (il  )

        t_snow_new_b (ifull) = t_snow_new_t (il)
        w_snow_new_b (ifull) = w_snow_new_t (il)
        rho_snow_new_b(ifull) = rho_snow_new_t (il)

    ENDDO
    !$acc end parallel

ENDIF ! end SNOWPOLINO

!VS >

    !$acc parallel
    DO k = 1, ke_soil+1
      !$acc loop gang vector private (ifull)
      DO il = 1, n_landgp(ib)
        ! grid point index in the full variable
        ifull = mind_landgp(il,ib)

        w_so_b              (ifull,k)  = w_so_now_t         (il,k)
        w_so_new_b          (ifull,k)  = w_so_new_t         (il,k)
        w_so_ice_b          (ifull,k)  = w_so_ice_now_t     (il,k)
        w_so_ice_new_b      (ifull,k)  = w_so_ice_new_t     (il,k)
      ENDDO
    ENDDO
    !$acc end parallel

    IF (lmulti_snow) THEN
      !$acc parallel
      DO k = 0, ke_snow
        !$acc loop gang vector private (ifull)
        DO il = 1, n_landgp(ib)
         ! grid point index in the full variable
          ifull = mind_landgp(il,ib)

          t_snow_mult_b       (ifull,k)  = t_snow_mult_now_t  (il,k)
          t_snow_mult_new_b   (ifull,k)  = t_snow_mult_new_t  (il,k)
        ENDDO
      ENDDO
      !$acc end parallel

      !$acc parallel
      DO k = 1, ke_snow
        !$acc loop gang vector private (ifull)
        DO il = 1, n_landgp(ib)
          ! grid point index in the full variable
          ifull = mind_landgp(il,ib)

          rho_snow_mult_b     (ifull,k)  = rho_snow_mult_now_t(il,k)
          rho_snow_mult_new_b (ifull,k)  = rho_snow_mult_new_t(il,k)
          wliq_snow_b         (ifull,k)  = wliq_snow_now_t    (il,k)
          wliq_snow_new_b     (ifull,k)  = wliq_snow_new_t    (il,k)
          w_snow_mult_b       (ifull,k)  = wtot_snow_now_t    (il,k)
          w_snow_mult_new_b   (ifull,k)  = wtot_snow_new_t    (il,k)
          dzh_snow_mult_b     (ifull,k)  = dzh_snow_now_t     (il,k)
          dzh_snow_mult_new_b (ifull,k)  = dzh_snow_new_t     (il,k)
        ENDDO
      ENDDO
    !$acc end parallel
    ENDIF

  ENDIF ! n_landgp(ib) > 0

!write(*,*) h_snow_t
!------------------------------------------------------------------------------
! End of module procedure tsa_sfc_organize
!------------------------------------------------------------------------------

END SUBROUTINE tsa_sfc_organize
       
!==============================================================================
!==============================================================================
!+ Module procedure "tsa_sfc_finalize" to finalize calls to TERRA
!------------------------------------------------------------------------------

SUBROUTINE tsa_sfc_finalize

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine handles the finalization of TERRA
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

! Local parameters:
! ----------------

INTEGER :: izerror

!------------------------------------------------------------------------------
! Begin Subroutine sfc_finalize
!------------------------------------------------------------------------------

  izerror = 0
  ! Deallocate the working arrays
!!DL  DEALLOCATE (czmls, czhls, msoilgrib)
  DEALLOCATE (zzhls, zdzhs, zdzms)

!!DL
  IF (ymodel == 'COSMO') THEN
  DEALLOCATE (   n_landgp)
  DEALLOCATE (mind_landgp)
  ENDIF

!!DL   IF (llake) THEN
!!DL     DEALLOCATE (   n_lakegp)
!!DL     DEALLOCATE (mind_lakegp)
!!DL   ENDIF

!!DL   IF (lseaice) THEN
!!DL     DEALLOCATE (mind_sicegp)
!!DL   ENDIF

  !$acc exit data delete (czmls, zzhls, zdzhs, zdzms)
!!DL   !$acc exit data delete (mind_landgp, n_landgp, mind_lakegp, n_lakegp, mind_sicegp)
  !$acc exit data delete (mind_landgp, n_landgp)

#ifdef ALLOC_WKARR
  CALL terra_wkarr_dealloc (izerror)
  CALL sfc_in_wkarr_dealloc

!!DL   IF (llake) THEN
!!DL     CALL flake_wkarr_dealloc (izerror)
!!DL   ENDIF
!!DL   IF (lseaice) THEN
!!DL     CALL seaice_wkarr_dealloc (izerror)
!!DL   ENDIF
  IF (lsnow) then
    CALL snow_wkarr_dealloc (izerror)
  endif

!US  IF (izerror /= 0) THEN
!US    CALL model_abort(my_cart_id, izerror, 'surface wkarr_dealloc', 'sfc_finalize')
!US  ENDIF
#endif

!------------------------------------------------------------------------------
! End of module procedure tsa_sfc_finalize
!------------------------------------------------------------------------------

END SUBROUTINE tsa_sfc_finalize

!==============================================================================
#ifdef ALLOC_WKARR
!==============================================================================
!+ Module procedure for allocating working arrays for the surface schemes
!------------------------------------------------------------------------------

SUBROUTINE sfc_in_wkarr_alloc (nproma, ke_soil, ke_snow)

!------------------------------------------------------------------------------
!
! Purpose:   Allocation of the working arrays required in the interface
!            for the surface schemes
!
!------------------------------------------------------------------------------

! Declarations:
INTEGER, INTENT(IN) :: &
  nproma, ke_soil, ke_snow

INTEGER :: izl, ist

! End of header
!==============================================================================

  izl = 0
  ist = 0

  ALLOCATE ( isoiltyp_t         (nproma), STAT=izl ); isoiltyp_t    = 0       ; ist=ist+izl
  ALLOCATE ( liceana_t          (nproma), STAT=izl ); liceana_t     = .FALSE. ; ist=ist+izl
  ALLOCATE ( plcov_t            (nproma), STAT=izl ); plcov_t       = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( rootdp_t           (nproma), STAT=izl ); rootdp_t      = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( sai_t              (nproma), STAT=izl ); sai_t         = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( tai_t              (nproma), STAT=izl ); tai_t         = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( eai_t              (nproma), STAT=izl ); eai_t         = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( skinc_t            (nproma), STAT=izl ); skinc_t       = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( rsmin2d_t          (nproma), STAT=izl ); rsmin2d_t     = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( sso_sigma_t        (nproma), STAT=izl ); sso_sigma_t   = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( ps_t               (nproma), STAT=izl ); ps_t          = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( t_snow_now_t       (nproma), STAT=izl ); t_snow_now_t  = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( t_snow_new_t       (nproma), STAT=izl ); t_snow_new_t  = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( t_s_now_t          (nproma), STAT=izl ); t_s_now_t     = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( t_s_new_t          (nproma), STAT=izl ); t_s_new_t     = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( t_sk_now_t         (nproma), STAT=izl ); t_sk_now_t    = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( t_sk_new_t         (nproma), STAT=izl ); t_sk_new_t    = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( t_g_now_t          (nproma), STAT=izl ); t_g_now_t     = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( t_g_new_t          (nproma), STAT=izl ); t_g_new_t     = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( qv_s_t             (nproma), STAT=izl ); qv_s_t        = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( w_snow_now_t       (nproma), STAT=izl ); w_snow_now_t  = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( w_snow_new_t       (nproma), STAT=izl ); w_snow_new_t  = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( rho_snow_now_t     (nproma), STAT=izl ); rho_snow_now_t= 0.0_wp  ; ist=ist+izl
  ALLOCATE ( rho_snow_new_t     (nproma), STAT=izl ); rho_snow_new_t= 0.0_wp  ; ist=ist+izl
  ALLOCATE ( h_snow_t           (nproma), STAT=izl ); h_snow_t      = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( h_snow_gp_t        (nproma), STAT=izl ); h_snow_gp_t   = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( meltrate_t         (nproma), STAT=izl ); meltrate_t    = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( snow_melt_t        (nproma), STAT=izl ); snow_melt_t   = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( tsnred_t           (nproma), STAT=izl ); tsnred_t      = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( w_i_now_t          (nproma), STAT=izl ); w_i_now_t     = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( w_i_new_t          (nproma), STAT=izl ); w_i_new_t     = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( w_p_now_t          (nproma), STAT=izl ); w_p_now_t     = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( w_p_new_t          (nproma), STAT=izl ); w_p_new_t     = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( w_s_now_t          (nproma), STAT=izl ); w_s_now_t     = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( w_s_new_t          (nproma), STAT=izl ); w_s_new_t     = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( u_10m_t            (nproma), STAT=izl ); u_10m_t       = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( v_10m_t            (nproma), STAT=izl ); v_10m_t       = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( freshsnow_t        (nproma), STAT=izl ); freshsnow_t   = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( zf_snow_t          (nproma), STAT=izl ); zf_snow_t     = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( prr_con_t          (nproma), STAT=izl ); prr_con_t     = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( prs_con_t          (nproma), STAT=izl ); prs_con_t     = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( conv_frac_t        (nproma), STAT=izl ); conv_frac_t   = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( prr_gsp_t          (nproma), STAT=izl ); prr_gsp_t     = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( prs_gsp_t          (nproma), STAT=izl ); prs_gsp_t     = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( prg_gsp_t          (nproma), STAT=izl ); prg_gsp_t     = 0.0_wp  ; ist=ist+izl
#ifdef TWOMOM_SB
  ALLOCATE ( prh_gsp_t          (nproma), STAT=izl ); prh_gsp_t     = 0.0_wp  ; ist=ist+izl
#endif
  ALLOCATE ( tch_t              (nproma), STAT=izl ); tch_t         = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( tcm_t              (nproma), STAT=izl ); tcm_t         = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( tfv_t              (nproma), STAT=izl ); tfv_t         = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( sobs_t             (nproma), STAT=izl ); sobs_t        = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( thbs_t             (nproma), STAT=izl ); thbs_t        = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( pabs_t             (nproma), STAT=izl ); pabs_t        = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( runoff_s_t         (nproma), STAT=izl ); runoff_s_t    = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( runoff_g_t         (nproma), STAT=izl ); runoff_g_t    = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( shfl_s_t           (nproma), STAT=izl ); shfl_s_t      = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( lhfl_s_t           (nproma), STAT=izl ); lhfl_s_t      = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( qvfl_s_t           (nproma), STAT=izl ); qvfl_s_t      = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( dum1fl_s_t         (nproma), STAT=izl ); dum1fl_s_t    = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( dum2fl_s_t         (nproma), STAT=izl ); dum2fl_s_t    = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( shfl_snow_t        (nproma), STAT=izl ); shfl_snow_t   = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( lhfl_snow_t        (nproma), STAT=izl ); lhfl_snow_t   = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( lhfl_bs_t          (nproma), STAT=izl ); lhfl_bs_t     = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( rstom_t            (nproma), STAT=izl ); rstom_t       = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( u_t                (nproma), STAT=izl ); u_t           = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( v_t                (nproma), STAT=izl ); v_t           = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( t_t                (nproma), STAT=izl ); t_t           = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( qv_t               (nproma), STAT=izl ); qv_t          = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( ptot_t             (nproma), STAT=izl ); ptot_t        = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( fr_paved_t         (nproma), STAT=izl ); fr_paved_t    = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( sa_uf_t            (nproma), STAT=izl ); sa_uf_t       = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( ai_uf_t            (nproma), STAT=izl ); ai_uf_t       = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( alb_red_uf_t       (nproma), STAT=izl ); alb_red_uf_t  = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( w_imp_t            (nproma), STAT=izl ); w_imp_t       = 0.0_wp  ; ist=ist+izl
  ALLOCATE ( w_isa_t            (nproma), STAT=izl ); w_isa_t       = 0.0_wp  ; ist=ist+izl

  ALLOCATE ( t_so_now_t         (nproma,0:ke_soil+1) ,STAT=izl ); t_so_now_t        = 0.0_wp; ist=ist+izl
  ALLOCATE ( t_so_new_t         (nproma,0:ke_soil+1) ,STAT=izl ); t_so_new_t        = 0.0_wp; ist=ist+izl

  ALLOCATE ( w_so_now_t         (nproma,  ke_soil+1) ,STAT=izl ); w_so_now_t        = 0.0_wp; ist=ist+izl
  ALLOCATE ( w_so_new_t         (nproma,  ke_soil+1) ,STAT=izl ); w_so_new_t        = 0.0_wp; ist=ist+izl
  ALLOCATE ( w_so_ice_now_t     (nproma,  ke_soil+1) ,STAT=izl ); w_so_ice_now_t    = 0.0_wp; ist=ist+izl
  ALLOCATE ( w_so_ice_new_t     (nproma,  ke_soil+1) ,STAT=izl ); w_so_ice_new_t    = 0.0_wp; ist=ist+izl

  ALLOCATE ( lhfl_pl_t          (nproma,  ke_soil+1) ,STAT=izl ); lhfl_pl_t         = 0.0_wp; ist=ist+izl

  ALLOCATE ( t_snow_mult_now_t  (nproma,0:ke_snow) ,STAT=izl ); t_snow_mult_now_t   = 0.0_wp; ist=ist+izl
  ALLOCATE ( t_snow_mult_new_t  (nproma,0:ke_snow) ,STAT=izl ); t_snow_mult_new_t   = 0.0_wp; ist=ist+izl

  ALLOCATE ( rho_snow_mult_now_t(nproma,  ke_snow) ,STAT=izl ); rho_snow_mult_now_t = 0.0_wp; ist=ist+izl
  ALLOCATE ( rho_snow_mult_new_t(nproma,  ke_snow) ,STAT=izl ); rho_snow_mult_new_t = 0.0_wp; ist=ist+izl
  ALLOCATE ( wliq_snow_now_t    (nproma,  ke_snow) ,STAT=izl ); wliq_snow_now_t     = 0.0_wp; ist=ist+izl
  ALLOCATE ( wliq_snow_new_t    (nproma,  ke_snow) ,STAT=izl ); wliq_snow_new_t     = 0.0_wp; ist=ist+izl
  ALLOCATE ( wtot_snow_now_t    (nproma,  ke_snow) ,STAT=izl ); wtot_snow_now_t     = 0.0_wp; ist=ist+izl
  ALLOCATE ( wtot_snow_new_t    (nproma,  ke_snow) ,STAT=izl ); wtot_snow_new_t     = 0.0_wp; ist=ist+izl
  ALLOCATE ( dzh_snow_now_t     (nproma,  ke_snow) ,STAT=izl ); dzh_snow_now_t      = 0.0_wp; ist=ist+izl
  ALLOCATE ( dzh_snow_new_t     (nproma,  ke_snow) ,STAT=izl ); dzh_snow_new_t      = 0.0_wp; ist=ist+izl

!VS <
  ALLOCATE ( t_sn_now_t         (nproma,1:ke_snow) ,STAT=izl ); t_sn_now_t    = 0.0_wp; ist=ist+izl
  ALLOCATE ( t_sn_new_t         (nproma,1:ke_snow) ,STAT=izl ); t_sn_new_t    = 0.0_wp; ist=ist+izl
  ALLOCATE ( theta_i_now_t      (nproma,1:ke_snow) ,STAT=izl ); theta_i_now_t = 0.0_wp; ist=ist+izl
  ALLOCATE ( theta_i_new_t      (nproma,1:ke_snow) ,STAT=izl ); theta_i_new_t = 0.0_wp; ist=ist+izl
  ALLOCATE ( theta_w_now_t      (nproma,1:ke_snow) ,STAT=izl ); theta_w_now_t = 0.0_wp; ist=ist+izl
  ALLOCATE ( theta_w_new_t      (nproma,1:ke_snow) ,STAT=izl ); theta_w_new_t = 0.0_wp; ist=ist+izl
  ALLOCATE ( theta_a_now_t      (nproma,1:ke_snow) ,STAT=izl ); theta_a_now_t = 0.0_wp; ist=ist+izl
  ALLOCATE ( theta_a_new_t      (nproma,1:ke_snow) ,STAT=izl ); theta_a_new_t = 0.0_wp; ist=ist+izl
  ALLOCATE ( dzm_sn_now_t       (nproma,1:ke_snow) ,STAT=izl ); dzm_sn_now_t  = 0.0_wp; ist=ist+izl
  ALLOCATE ( dzm_sn_new_t       (nproma,1:ke_snow) ,STAT=izl ); dzm_sn_new_t  = 0.0_wp; ist=ist+izl
  ALLOCATE ( hn_sn_now_t        (nproma          ) ,STAT=izl ); hn_sn_now_t   = 0.0_wp; ist=ist+izl
  ALLOCATE ( hn_sn_new_t        (nproma          ) ,STAT=izl ); hn_sn_new_t   = 0.0_wp; ist=ist+izl
  ALLOCATE ( top_sn_now_t       (nproma          ) ,STAT=izl ); top_sn_now_t  = 0.0_wp; ist=ist+izl
  ALLOCATE ( top_sn_new_t       (nproma          ) ,STAT=izl ); top_sn_new_t  = 0.0_wp; ist=ist+izl
!VS >


  ! allocate all these arrays for the device
  !$acc enter data                                                               &
  !$acc create ( isoiltyp_t    , plcov_t       , rootdp_t      , sai_t         ) &
  !$acc create ( tai_t         , eai_t         , skinc_t       , rsmin2d_t     ) &
  !$acc create ( sso_sigma_t   , liceana_t                                     ) &
  !$acc create ( ps_t          , t_snow_now_t  , t_snow_new_t  , t_s_now_t     ) &
  !$acc create ( t_s_new_t     , t_sk_now_t    , t_sk_new_t    , t_g_now_t     ) &
  !$acc create ( t_g_new_t     , qv_s_t                                        ) &
  !$acc create ( w_snow_now_t  , w_snow_new_t  , rho_snow_now_t, rho_snow_new_t) &
  !$acc create ( h_snow_t      , h_snow_gp_t   , meltrate_t    , snow_melt_t   ) &
  !$acc create ( tsnred_t      , w_i_now_t     , w_i_new_t     , w_p_now_t     ) &
  !$acc create ( w_p_new_t     , w_s_now_t     , w_s_new_t     , u_10m_t       ) &
  !$acc create ( v_10m_t       , freshsnow_t   , zf_snow_t     , prr_con_t     ) &
  !$acc create ( prs_con_t     , conv_frac_t   , prr_gsp_t     , prs_gsp_t     ) &
  !$acc create ( prg_gsp_t     , tch_t         , tcm_t         , tfv_t         ) &
  !$acc create ( sobs_t        , thbs_t        , pabs_t        , runoff_s_t    ) &
  !$acc create ( runoff_g_t    , shfl_s_t      , lhfl_s_t      , qvfl_s_t      ) &
  !$acc create ( dum1fl_s_t    , dum2fl_s_t    , shfl_snow_t   , lhfl_snow_t   ) &
  !$acc create ( lhfl_bs_t     , rstom_t       , u_t           , v_t           ) &
  !$acc create ( t_t           , qv_t          , ptot_t        , fr_paved_t    ) &
  !$acc create ( sa_uf_t       , ai_uf_t       , alb_red_uf_t  , w_imp_t       ) &
  !$acc create ( w_isa_t       , t_so_now_t    , t_so_new_t    , lhfl_pl_t     ) &
  !$acc create ( w_so_now_t    , w_so_new_t    , w_so_ice_now_t, w_so_ice_new_t) &

  !$acc create ( t_snow_mult_now_t  , t_snow_mult_new_t  , rho_snow_mult_now_t ) &
  !$acc create ( rho_snow_mult_new_t, wliq_snow_now_t    , wliq_snow_new_t     ) &
  !$acc create ( wtot_snow_now_t    , wtot_snow_new_t    , dzh_snow_now_t      ) &
  !$acc create ( dzh_snow_new_t                                                ) &

!VS <
!  !$acc create (swdir_s_t      , swdifd_s_t      , swdifu_s_t, lwd_s_t, lwu_s_t) &
  !$acc create (t_sn_now_t     , t_sn_new_t    , theta_i_now_t , theta_i_new_t ) &
  !$acc create (theta_w_now_t  , theta_w_new_t , theta_a_now_t , theta_a_new_t ) &
  !$acc create (dzm_sn_now_t   , dzm_sn_new_t  , top_sn_now_t  , top_sn_new_t  ) &
  !$acc create (hn_sn_now_t    , hn_sn_new_t                                   )
!VS >



#ifdef TWOMOM_SB
  !$acc create ( prh_gsp_t     )
#endif

!!DL   IF (llake .OR. lseaice) THEN
!!DL     ! lake input fields
!!DL     ALLOCATE ( fr_lake_t    (nproma), STAT=izl ); fr_lake_t    = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( depth_lk_t   (nproma), STAT=izl ); depth_lk_t   = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( fetch_lk_t   (nproma), STAT=izl ); fetch_lk_t   = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( dp_bs_lk_t   (nproma), STAT=izl ); dp_bs_lk_t   = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( t_bs_lk_t    (nproma), STAT=izl ); t_bs_lk_t    = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( gamso_lk_t   (nproma), STAT=izl ); gamso_lk_t   = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( fc_t         (nproma), STAT=izl ); fc_t         = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( qmomflux_t   (nproma), STAT=izl ); qmomflux_t   = 0.0_wp  ; ist=ist+izl

!!DL     ! lake fields for previous time step
!!DL     ALLOCATE ( t_snow_t_p   (nproma), STAT=izl ); t_snow_t_p   = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( h_snow_t_p   (nproma), STAT=izl ); h_snow_t_p   = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( t_ice_t_p    (nproma), STAT=izl ); t_ice_t_p    = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( h_ice_t_p    (nproma), STAT=izl ); h_ice_t_p    = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( t_mnw_lk_t_p (nproma), STAT=izl ); t_mnw_lk_t_p = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( t_wml_lk_t_p (nproma), STAT=izl ); t_wml_lk_t_p = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( t_bot_lk_t_p (nproma), STAT=izl ); t_bot_lk_t_p = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( c_t_lk_t_p   (nproma), STAT=izl ); c_t_lk_t_p   = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( h_ml_lk_t_p  (nproma), STAT=izl ); h_ml_lk_t_p  = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( t_b1_lk_t_p  (nproma), STAT=izl ); t_b1_lk_t_p  = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( h_b1_lk_t_p  (nproma), STAT=izl ); h_b1_lk_t_p  = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( t_s_t_p      (nproma), STAT=izl ); t_s_t_p      = 0.0_wp  ; ist=ist+izl

!!DL     ! lake fields for new time step
!!DL     ALLOCATE ( t_snow_t_n   (nproma), STAT=izl ); t_snow_t_n   = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( h_snow_t_n   (nproma), STAT=izl ); h_snow_t_n   = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( t_ice_t_n    (nproma), STAT=izl ); t_ice_t_n    = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( h_ice_t_n    (nproma), STAT=izl ); h_ice_t_n    = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( t_mnw_lk_t_n (nproma), STAT=izl ); t_mnw_lk_t_n = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( t_wml_lk_t_n (nproma), STAT=izl ); t_wml_lk_t_n = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( t_bot_lk_t_n (nproma), STAT=izl ); t_bot_lk_t_n = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( c_t_lk_t_n   (nproma), STAT=izl ); c_t_lk_t_n   = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( h_ml_lk_t_n  (nproma), STAT=izl ); h_ml_lk_t_n  = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( t_b1_lk_t_n  (nproma), STAT=izl ); t_b1_lk_t_n  = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( h_b1_lk_t_n  (nproma), STAT=izl ); h_b1_lk_t_n  = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( t_s_t_n      (nproma), STAT=izl ); t_s_t_n      = 0.0_wp  ; ist=ist+izl

!!DL     ! and some sea ice fields
!!DL     ALLOCATE ( fr_ice_t     (nproma), STAT=izl ); fr_ice_t     = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( albsi_t_p    (nproma), STAT=izl ); albsi_t_p    = 0.0_wp  ; ist=ist+izl
!!DL     ALLOCATE ( albsi_t_n    (nproma), STAT=izl ); albsi_t_n    = 0.0_wp  ; ist=ist+izl

!!DL     ! allocate all these arrays for the device
!!DL     !$acc enter data                                                               &
!!DL     ! lake input fields
!!DL     !$acc create ( fr_lake_t    , depth_lk_t   , fetch_lk_t   , dp_bs_lk_t   )     &
!!DL     !$acc create ( t_bs_lk_t    , gamso_lk_t   , fc_t         , qmomflux_t   )     &

!!DL     ! lake fields for previous time step
!!DL     !$acc create ( t_snow_t_p   , h_snow_t_p   , t_ice_t_p    , h_ice_t_p    )     &
!!DL     !$acc create ( t_mnw_lk_t_p , t_wml_lk_t_p , t_bot_lk_t_p , c_t_lk_t_p   )     &
!!DL     !$acc create ( h_ml_lk_t_p  , t_b1_lk_t_p  , h_b1_lk_t_p  , t_s_t_p      )     &

!!DL     ! lake fields for new time step
!!DL     !$acc create ( t_snow_t_n   , h_snow_t_n   , t_ice_t_n    , h_ice_t_n    )     &
!!DL     !$acc create ( t_mnw_lk_t_n , t_wml_lk_t_n , t_bot_lk_t_n , c_t_lk_t_n   )     &
!!DL     !$acc create ( h_ml_lk_t_n  , t_b1_lk_t_n  , h_b1_lk_t_n  , t_s_t_n      )     &
!!DL 
!!DL     ! and some sea ice fields
!!DL     !$acc create ( fr_ice_t     , albsi_t_n    , albsi_t_p    )

!!DL   ENDIF


  ! local fields for the interface for fields which are not activated in COSMO right now
  ALLOCATE ( h_snow_gp_b (nproma), STAT=izl); h_snow_gp_b = 0.0_wp; ist=ist+izl
  ALLOCATE ( conv_frac_b (nproma), STAT=izl); conv_frac_b = 0.0_wp; ist=ist+izl

  ! the next variables are output from terra, but not further processed in COSMO
  ALLOCATE ( meltrate_b  (nproma), STAT=izl); meltrate_b  = 0.0_wp; ist=ist+izl
  ALLOCATE ( shfl_snow_b (nproma), STAT=izl); shfl_snow_b = 0.0_wp; ist=ist+izl
  ALLOCATE ( lhfl_snow_b (nproma), STAT=izl); lhfl_snow_b = 0.0_wp; ist=ist+izl

  ! allocate all these arrays for the device
  !$acc enter data                                                               &
  !$acc create ( h_snow_gp_b, conv_frac_b, meltrate_b, shfl_snow_b, lhfl_snow_b )

!US  IF (ist /= 0) THEN
!US    CALL model_abort(my_cart_id, ist, 'Allocation of the interface working arrays failed', &
!US                     'sfc_in_wkarr_alloc')
!US  ENDIF

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE sfc_in_wkarr_alloc

!==============================================================================
!==============================================================================
!+ Module procedure for deallocating working arrays for the surface interface
!------------------------------------------------------------------------------

SUBROUTINE sfc_in_wkarr_dealloc

!------------------------------------------------------------------------------
!
! Purpose:   Allocation of the working arrays required in the interface
!            for the surface schemes
!
!------------------------------------------------------------------------------

! Declarations:

INTEGER :: izl, ist

! End of header
!==============================================================================

  izl = 0
  ist = 0

  DEALLOCATE ( isoiltyp_t         , STAT=izl ); ist=ist+izl
  DEALLOCATE ( liceana_t          , STAT=izl ); ist=ist+izl
  DEALLOCATE ( plcov_t            , STAT=izl ); ist=ist+izl
  DEALLOCATE ( rootdp_t           , STAT=izl ); ist=ist+izl
  DEALLOCATE ( sai_t              , STAT=izl ); ist=ist+izl
  DEALLOCATE ( tai_t              , STAT=izl ); ist=ist+izl
  DEALLOCATE ( eai_t              , STAT=izl ); ist=ist+izl
  DEALLOCATE ( skinc_t            , STAT=izl ); ist=ist+izl
  DEALLOCATE ( rsmin2d_t          , STAT=izl ); ist=ist+izl
  DEALLOCATE ( sso_sigma_t        , STAT=izl ); ist=ist+izl
  DEALLOCATE ( ps_t               , STAT=izl ); ist=ist+izl
  DEALLOCATE ( t_snow_now_t       , STAT=izl ); ist=ist+izl
  DEALLOCATE ( t_snow_new_t       , STAT=izl ); ist=ist+izl
  DEALLOCATE ( t_s_now_t          , STAT=izl ); ist=ist+izl
  DEALLOCATE ( t_s_new_t          , STAT=izl ); ist=ist+izl
  DEALLOCATE ( t_sk_now_t         , STAT=izl ); ist=ist+izl
  DEALLOCATE ( t_sk_new_t         , STAT=izl ); ist=ist+izl
  DEALLOCATE ( t_g_now_t          , STAT=izl ); ist=ist+izl
  DEALLOCATE ( t_g_new_t          , STAT=izl ); ist=ist+izl
  DEALLOCATE ( qv_s_t             , STAT=izl ); ist=ist+izl
  DEALLOCATE ( w_snow_now_t       , STAT=izl ); ist=ist+izl
  DEALLOCATE ( w_snow_new_t       , STAT=izl ); ist=ist+izl
  DEALLOCATE ( rho_snow_now_t     , STAT=izl ); ist=ist+izl
  DEALLOCATE ( rho_snow_new_t     , STAT=izl ); ist=ist+izl
  DEALLOCATE ( h_snow_t           , STAT=izl ); ist=ist+izl
  DEALLOCATE ( h_snow_gp_t        , STAT=izl ); ist=ist+izl
  DEALLOCATE ( meltrate_t         , STAT=izl ); ist=ist+izl
  DEALLOCATE ( snow_melt_t        , STAT=izl ); ist=ist+izl
  DEALLOCATE ( tsnred_t           , STAT=izl ); ist=ist+izl
  DEALLOCATE ( w_i_now_t          , STAT=izl ); ist=ist+izl
  DEALLOCATE ( w_i_new_t          , STAT=izl ); ist=ist+izl
  DEALLOCATE ( w_p_now_t          , STAT=izl ); ist=ist+izl
  DEALLOCATE ( w_p_new_t          , STAT=izl ); ist=ist+izl
  DEALLOCATE ( w_s_now_t          , STAT=izl ); ist=ist+izl
  DEALLOCATE ( w_s_new_t          , STAT=izl ); ist=ist+izl
  DEALLOCATE ( u_10m_t            , STAT=izl ); ist=ist+izl
  DEALLOCATE ( v_10m_t            , STAT=izl ); ist=ist+izl
  DEALLOCATE ( freshsnow_t        , STAT=izl ); ist=ist+izl
  DEALLOCATE ( zf_snow_t          , STAT=izl ); ist=ist+izl
  DEALLOCATE ( prr_con_t          , STAT=izl ); ist=ist+izl
  DEALLOCATE ( prs_con_t          , STAT=izl ); ist=ist+izl
  DEALLOCATE ( conv_frac_t        , STAT=izl ); ist=ist+izl
  DEALLOCATE ( prr_gsp_t          , STAT=izl ); ist=ist+izl
  DEALLOCATE ( prs_gsp_t          , STAT=izl ); ist=ist+izl
  DEALLOCATE ( prg_gsp_t          , STAT=izl ); ist=ist+izl
#ifdef TWOMOM_SB
  DEALLOCATE ( prh_gsp_t          , STAT=izl ); ist=ist+izl
#endif
  DEALLOCATE ( tch_t              , STAT=izl ); ist=ist+izl
  DEALLOCATE ( tcm_t              , STAT=izl ); ist=ist+izl
  DEALLOCATE ( tfv_t              , STAT=izl ); ist=ist+izl
  DEALLOCATE ( sobs_t             , STAT=izl ); ist=ist+izl
  DEALLOCATE ( thbs_t             , STAT=izl ); ist=ist+izl
  DEALLOCATE ( pabs_t             , STAT=izl ); ist=ist+izl
  DEALLOCATE ( runoff_s_t         , STAT=izl ); ist=ist+izl
  DEALLOCATE ( runoff_g_t         , STAT=izl ); ist=ist+izl
  DEALLOCATE ( shfl_s_t           , STAT=izl ); ist=ist+izl
  DEALLOCATE ( lhfl_s_t           , STAT=izl ); ist=ist+izl
  DEALLOCATE ( qvfl_s_t           , STAT=izl ); ist=ist+izl
  DEALLOCATE ( dum1fl_s_t         , STAT=izl ); ist=ist+izl
  DEALLOCATE ( dum2fl_s_t         , STAT=izl ); ist=ist+izl
  DEALLOCATE ( shfl_snow_t        , STAT=izl ); ist=ist+izl
  DEALLOCATE ( lhfl_snow_t        , STAT=izl ); ist=ist+izl
  DEALLOCATE ( lhfl_bs_t          , STAT=izl ); ist=ist+izl
  DEALLOCATE ( rstom_t            , STAT=izl ); ist=ist+izl
  DEALLOCATE ( u_t                , STAT=izl ); ist=ist+izl
  DEALLOCATE ( v_t                , STAT=izl ); ist=ist+izl
  DEALLOCATE ( t_t                , STAT=izl ); ist=ist+izl
  DEALLOCATE ( qv_t               , STAT=izl ); ist=ist+izl
  DEALLOCATE ( ptot_t             , STAT=izl ); ist=ist+izl
  DEALLOCATE ( fr_paved_t         , STAT=izl ); ist=ist+izl
  DEALLOCATE ( sa_uf_t            , STAT=izl ); ist=ist+izl
  DEALLOCATE ( ai_uf_t            , STAT=izl ); ist=ist+izl
  DEALLOCATE ( alb_red_uf_t       , STAT=izl ); ist=ist+izl
  DEALLOCATE ( w_imp_t            , STAT=izl ); ist=ist+izl
  DEALLOCATE ( w_isa_t            , STAT=izl ); ist=ist+izl

  DEALLOCATE ( t_so_now_t         , STAT=izl ); ist=ist+izl
  DEALLOCATE ( t_so_new_t         , STAT=izl ); ist=ist+izl

  DEALLOCATE ( w_so_now_t         , STAT=izl ); ist=ist+izl
  DEALLOCATE ( w_so_new_t         , STAT=izl ); ist=ist+izl
  DEALLOCATE ( w_so_ice_now_t     , STAT=izl ); ist=ist+izl
  DEALLOCATE ( w_so_ice_new_t     , STAT=izl ); ist=ist+izl

  DEALLOCATE ( lhfl_pl_t          , STAT=izl ); ist=ist+izl

  DEALLOCATE ( t_snow_mult_now_t  , STAT=izl ); ist=ist+izl
  DEALLOCATE ( t_snow_mult_new_t  , STAT=izl ); ist=ist+izl

  DEALLOCATE ( rho_snow_mult_now_t, STAT=izl ); ist=ist+izl
  DEALLOCATE ( rho_snow_mult_new_t, STAT=izl ); ist=ist+izl
  DEALLOCATE ( wliq_snow_now_t    , STAT=izl ); ist=ist+izl
  DEALLOCATE ( wliq_snow_new_t    , STAT=izl ); ist=ist+izl
  DEALLOCATE ( wtot_snow_now_t    , STAT=izl ); ist=ist+izl
  DEALLOCATE ( wtot_snow_new_t    , STAT=izl ); ist=ist+izl
  DEALLOCATE ( dzh_snow_now_t     , STAT=izl ); ist=ist+izl
  DEALLOCATE ( dzh_snow_new_t     , STAT=izl ); ist=ist+izl

!VS <
!  DEALLOCATE ( swdir_s_t          , STAT=izl ); ist=ist+izl
!  DEALLOCATE ( swdifd_s_t         , STAT=izl ); ist=ist+izl
!  DEALLOCATE ( swdifu_s_t         , STAT=izl ); ist=ist+izl
!  DEALLOCATE ( lwd_s_t            , STAT=izl ); ist=ist+izl
!  DEALLOCATE ( lwu_s_t            , STAT=izl ); ist=ist+izl
  DEALLOCATE ( t_sn_now_t         , STAT=izl ); ist=ist+izl
  DEALLOCATE ( t_sn_new_t         , STAT=izl ); ist=ist+izl
  DEALLOCATE ( theta_i_now_t      , STAT=izl ); ist=ist+izl
  DEALLOCATE ( theta_i_new_t      , STAT=izl ); ist=ist+izl
  DEALLOCATE ( theta_w_now_t      , STAT=izl ); ist=ist+izl
  DEALLOCATE ( theta_w_new_t      , STAT=izl ); ist=ist+izl
  DEALLOCATE ( theta_a_now_t      , STAT=izl ); ist=ist+izl
  DEALLOCATE ( theta_a_new_t      , STAT=izl ); ist=ist+izl
  DEALLOCATE ( dzm_sn_now_t       , STAT=izl ); ist=ist+izl
  DEALLOCATE ( dzm_sn_new_t       , STAT=izl ); ist=ist+izl
  DEALLOCATE ( hn_sn_now_t        , STAT=izl ); ist=ist+izl
  DEALLOCATE ( hn_sn_new_t        , STAT=izl ); ist=ist+izl
  DEALLOCATE ( top_sn_now_t       , STAT=izl ); ist=ist+izl
  DEALLOCATE ( top_sn_new_t       , STAT=izl ); ist=ist+izl
!VS >



  ! deallocate all these arrays for the device
  !$acc exit  data                                                               &
  !$acc delete ( isoiltyp_t    , plcov_t       , rootdp_t      , sai_t         ) &
  !$acc delete ( tai_t         , eai_t         , skinc_t       , rsmin2d_t     ) &
  !$acc delete ( sso_sigma_t   , liceana_t                                     ) &
  !$acc delete ( ps_t          , t_snow_now_t  , t_snow_new_t  , t_s_now_t     ) &
  !$acc delete ( t_s_new_t     , t_sk_now_t    , t_sk_new_t    , t_g_now_t     ) &
  !$acc delete ( t_g_new_t     , qv_s_t                                        ) &
  !$acc delete ( w_snow_now_t  , w_snow_new_t  , rho_snow_now_t, rho_snow_new_t) &
  !$acc delete ( h_snow_t      , h_snow_gp_t   , meltrate_t    , snow_melt_t   ) &
  !$acc delete ( tsnred_t      , w_i_now_t     , w_i_new_t     , w_p_now_t     ) &
  !$acc delete ( w_p_new_t     , w_s_now_t     , w_s_new_t     , u_10m_t       ) &
  !$acc delete ( v_10m_t       , freshsnow_t   , zf_snow_t     , prr_con_t     ) &
  !$acc delete ( prs_con_t     , conv_frac_t   , prr_gsp_t     , prs_gsp_t     ) &
  !$acc delete ( prg_gsp_t     , tch_t         , tcm_t         , tfv_t         ) &
  !$acc delete ( sobs_t        , thbs_t        , pabs_t        , runoff_s_t    ) &
  !$acc delete ( runoff_g_t    , shfl_s_t      , lhfl_s_t      , qvfl_s_t      ) &
  !$acc delete ( dum1fl_s_t    , dum2fl_s_t    , shfl_snow_t   , lhfl_snow_t   ) &
  !$acc delete ( lhfl_bs_t     , rstom_t       , u_t           , v_t           ) &
  !$acc delete ( t_t           , qv_t          , ptot_t        , fr_paved_t    ) &
  !$acc delete ( sa_uf_t       , ai_uf_t       , alb_red_uf_t  , w_imp_t       ) &
  !$acc delete ( w_isa_t       , t_so_now_t    , t_so_new_t    , lhfl_pl_t     ) &
  !$acc delete ( w_so_now_t    , w_so_new_t    , w_so_ice_now_t, w_so_ice_new_t) &

  !$acc delete ( t_snow_mult_now_t  , t_snow_mult_new_t  , rho_snow_mult_now_t ) &
  !$acc delete ( rho_snow_mult_new_t, wliq_snow_now_t    , wliq_snow_new_t     ) &
  !$acc delete ( wtot_snow_now_t    , wtot_snow_new_t    , dzh_snow_now_t      ) &
  !$acc delete ( dzh_snow_new_t                                                ) &

!VS <
!  !$acc delete (swdir_s_t,       swdifd_s_t,    swdifu_s_t,   lwd_s_t,  lwu_s_t) &
  !$acc delete (t_sn_now_t     , t_sn_new_t    , theta_i_now_t , theta_i_new_t ) &
  !$acc delete (theta_w_now_t  , theta_w_new_t , theta_a_now_t , theta_a_new_t ) &
  !$acc delete (dzm_sn_now_t   , dzm_sn_new_t  , top_sn_now_t  , top_sn_new_t  ) &
  !$acc delete (hn_sn_now_t    , hn_sn_new_t                                   )
!VS >


#ifdef TWOMOM_SB
  !$acc delete ( prh_gsp_t     )
#endif

!!DL   IF (llake .OR. lseaice) THEN
!!DL     ! lake input fields
!!DL     DEALLOCATE ( fr_lake_t    , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( depth_lk_t   , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( fetch_lk_t   , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( dp_bs_lk_t   , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( t_bs_lk_t    , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( gamso_lk_t   , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( fc_t         , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( qmomflux_t   , STAT=izl ); ist=ist+izl

!!DL     ! lake fields for previous time step
!!DL     DEALLOCATE ( t_snow_t_p   , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( h_snow_t_p   , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( t_ice_t_p    , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( h_ice_t_p    , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( t_mnw_lk_t_p , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( t_wml_lk_t_p , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( t_bot_lk_t_p , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( c_t_lk_t_p   , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( h_ml_lk_t_p  , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( t_b1_lk_t_p  , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( h_b1_lk_t_p  , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( t_s_t_p      , STAT=izl ); ist=ist+izl

!!DL     ! lake fields for new time step
!!DL     DEALLOCATE ( t_snow_t_n   , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( h_snow_t_n   , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( t_ice_t_n    , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( h_ice_t_n    , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( t_mnw_lk_t_n , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( t_wml_lk_t_n , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( t_bot_lk_t_n , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( c_t_lk_t_n   , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( h_ml_lk_t_n  , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( t_b1_lk_t_n  , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( h_b1_lk_t_n  , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( t_s_t_n      , STAT=izl ); ist=ist+izl

!!DL     ! and some sea ice fields
!!DL     DEALLOCATE ( fr_ice_t     , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( albsi_t_p    , STAT=izl ); ist=ist+izl
!!DL     DEALLOCATE ( albsi_t_n    , STAT=izl ); ist=ist+izl

!!DL     ! deallocate all these arrays for the device
!!DL     !$acc exit  data                                                               &
!!DL     ! lake input fields
!!DL     !$acc delete ( fr_lake_t    , depth_lk_t   , fetch_lk_t   , dp_bs_lk_t   )     &
!!DL     !$acc delete ( t_bs_lk_t    , gamso_lk_t   , fc_t         , qmomflux_t   )     &

!!DL     ! lake fields for previous time step
!!DL     !$acc delete ( t_snow_t_p   , h_snow_t_p   , t_ice_t_p    , h_ice_t_p    )     &
!!DL     !$acc delete ( t_mnw_lk_t_p , t_wml_lk_t_p , t_bot_lk_t_p , c_t_lk_t_p   )     &
!!DL     !$acc delete ( h_ml_lk_t_p  , t_b1_lk_t_p  , h_b1_lk_t_p  , t_s_t_p      )     &

!!DL     ! lake fields for new time step
!!DL     !$acc delete ( t_snow_t_n   , h_snow_t_n   , t_ice_t_n    , h_ice_t_n    )     &
!!DL     !$acc delete ( t_mnw_lk_t_n , t_wml_lk_t_n , t_bot_lk_t_n , c_t_lk_t_n   )     &
!!DL     !$acc delete ( h_ml_lk_t_n  , t_b1_lk_t_n  , h_b1_lk_t_n  , t_s_t_n      )     &

!!DL     ! and some sea ice fields
!!DL     !$acc delete ( fr_ice_t     , albsi_t_p    , albsi_t_n    )

!!DL   ENDIF


  ! local fields for the interface for fields which are not activated in COSMO right now
  DEALLOCATE ( h_snow_gp_b , STAT=izl); ist=ist+izl
  DEALLOCATE ( conv_frac_b , STAT=izl); ist=ist+izl

  ! the next variables are output from terra, but not further processed in COSMO
  DEALLOCATE ( meltrate_b  , STAT=izl); ist=ist+izl
  DEALLOCATE ( shfl_snow_b , STAT=izl); ist=ist+izl
  DEALLOCATE ( lhfl_snow_b , STAT=izl); ist=ist+izl

  ! deallocate all these arrays for the device
  !$acc exit  data                                                               &
  !$acc delete ( h_snow_gp_b, conv_frac_b, meltrate_b, shfl_snow_b, lhfl_snow_b )

!US  IF (ist /= 0) THEN
!US    CALL model_abort(my_cart_id, ist, 'Deallocation of the interface working arrays failed', &
!US                     'sfc_in_wkarr_dealloc')
!US  ENDIF

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE sfc_in_wkarr_dealloc

!==============================================================================
#endif
!==============================================================================

END MODULE tsa_sfc_interface
