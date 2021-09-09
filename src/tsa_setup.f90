! ---------------------------------------------------------------------
! Routines to allocate, initialize and de-allocate the COSMO environment
! for the Terra-Module
!==============================================================================

MODULE tsa_setup

!==============================================================================
!
! Description:
!  This file contains subroutines that handle the setup of TSA:
!   read the namelist variables, allocate variables, initialize variables and 
!   de-allocate variables
!  required for the TSA version of th soil module terra.
!
!
!  Currently included:
!
!    - read_namelist
!      reads all the namelist variables to define the configuration
!
!    - allocate_fields:
!      allocates variables required for TSA
!
!NEW - allocate_block_fields:
!      allocates variables required for TSA (block structure)
!
!    - init_variables:
!      initialize variables required for TSA
!
!    - clean_up:
!      de-allocates variables required for TSA
!
!
! Current Code Owner: DWD, Juergen Helmert
!  phone:  +49  69  8062 2704
!  fax:    +49  69  8062 3721
!  email:  Juergen.Helmert@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! V1.0       2003/10/01 Felix Ament, University of Bonn
!  Initial version
! V2.0       2006/08/01 Felix Ament, MeteoSwiss
!  Felix Ament, MeteoSwiss, August 2006
! V4.13      2010/08/01 Guy de Morsier, MeteoSwiss (GDM)
!  Revised code to correspond to FieldExtra and COSMO version 4.13
! V5.01      2015-12-01 Yiftach Ziv, IMS (XYZ)
!  Arranged code to adhere to coding standards.
!  Removed subroutine constants, all constants available from data_constants
! V5.07      2020-07-15 Doerte Liermann, Ulrich Schaettler
!  Added subroutine to allocate block variables
!  General refactoring: renamed original module terra_lmenv to tsa_setup
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!======================================================

! USE declarations
USE kind_parameters,       ONLY:                                               &
    wp           ! KIND-type parameter for real variables

USE tsa_data, ONLY:                                                            &
    mind_ilon, mind_jlat, czhls_const, qv, qv_bd, prr_gsp_bd, prs_gsp_bd,      &
    so_down_bd, th_down_bd, pabs_bd, ps_bd, uuid, astore, store,               &
    runoff, rlat_geo, rlon_geo, ydate_ini, ydate_end,                          &
    z0_b, qv_b, zflmg, nstore, nastore,                                        &
    nout_interval, outdir, outprefix, ntype_output, lconstout,                 &
    radofiledir, metfiledir, metfileprefix,                                    &
    ntype_atminput, ntype_radinput, ntype_raininput, rundef,                   &
    w_g0, w_i0, t_soil0, w_snow0, t_cl0, w_cl0,                                &
    soilinitdir, soilinitprefix, soiltyp_const,                                &
    ke_model, ngp, ngpi, ngrid, ymodel, yform_write, yform_read,               &
    ntstep_max, tincr_max, constfilename, nsub_x, nsub_y, which_subreg,        &
    itype_hydparam,                                                            &
    lrel_in        , lvol_in        , lvegadapt      , lhomosoil      ,        &
    lhomoinit      , lconstout      , lrootadapt     , lcheck         ,        &
    lmulti_in      ,                  lz0local       , lconstvegalb   ,        &
    lext_monthly   , linfil_revised , lhourly_data   , ldestaggeruv   ,        &
    lpar           , lgettcl        , lcalc          ,                         &
    lai_const      , lai_min_const  , lai_max_const  , plcov_const    ,        &
    plcov_min_const, plcov_max_const, rootdp_const   , plcov_mn       ,        &
    plcov_mx       , rootdp0        , lai_mn         , lai_mx         ,        &
    z0             , rstom_mx       , rstom_mn       , vegalb         ,        &
    z0_const       , rstom_mn_const , rstom_mx_const , vegalb_const   ,        &
    dz             , dz_u           , nrecmax        , hsurface       ,        &
    rain_fac       , lat


USE data_io,               ONLY:                                               &
    lana_rho_snow

USE data_fields,           ONLY:                                               &
    soiltyp      , plcov        , rootdp       , sai          , tai          , &
    eai          , lai          , llandmask    , skinc        , sso_sigma    , &
    rsmin2d      , rstom        , rlat         , rlon         , hsurf        , &
    fr_wi        , fr_snow      , fr_land      ,                               &
    u            , v            , t            , pp           , ps           , &
    u_bd         ,                t_bd         ,                               &
    p0           ,                                                             &
    t_s          , t_sk         , t_g          , t_snow       , qv_s         , &
    w_snow       , w_i          , w_s          , w_p          ,                &
    prr_gsp      , prs_gsp      , prg_gsp      , prr_con      , prs_con      , &
    tch          , tcm          , tfh          , tfv          , sobs         , &
    thbs         , pabs         , runoff_s     , runoff_g     , qvfl_s       , &
    lhfl_s       , shfl_s       , t_2m         , u_10m        , v_10m        , &
    h_snow       , rho_snow     , wliq_snow    , freshsnow    , lhfl_pl      , &
    lhfl_bs      , t_so         , w_so         , w_so_ice     , t_cl         , &
    snow_melt    ,                                                             &
    w_snow_mult  , dzh_snow_mult, t_snow_mult  , rho_snow_mult, &
!VS < 
    theta_i,theta_w,theta_a,dzm_sn,t_sn,top_sn
!VS >

USE data_block_fields, ONLY:                                                   &
    hsurf_b        , isoiltyp_b     , fr_land_b      , fr_snow_b      ,        &
    plcov_b        , rootdp_b       , sso_sigma_b    , rsmin2d_b      ,        &
    rstom_b        , skinc_b        , sai_b          , eai_b          ,        &
    tai_b          ,                                                           &
    u_m_b          , v_m_b          , t_b            ,                         &
    p0_b           , pp_b           , ps_b           , ptot_b         ,        &
    t_g_b          , t_g_new_b      , t_s_b          , t_s_new_b      ,        &
    t_sk_b         , t_sk_new_b     , t_snow_b       , t_snow_new_b   ,        &
    t_so_b         , t_so_new_b     ,                                          &
    w_i_b          , w_i_new_b      , w_p_b          , w_p_new_b      ,        &
    w_s_b          , w_s_new_b      ,                                          &
    qv_s_b         , qv_s_new_b     , freshsnow_b    , h_snow_b       ,        &
    h_snow_new_b   , rho_snow_b     , rho_snow_new_b , runoff_g_b     ,        &
    runoff_s_b     , snow_melt_b    ,                                          &
    w_snow_b       , w_snow_new_b   , w_so_b         , w_so_new_b     ,        &
    w_so_ice_b     , w_so_ice_new_b , wliq_snow_b    , wliq_snow_new_b,        &
    prr_gsp_b      , prs_gsp_b      , prg_gsp_b      , prh_gsp_b      ,        &
    prr_con_b      , prs_con_b      ,                                          &
    qvfl_s_b       , shfl_s_b       , lhfl_bs_b      , lhfl_pl_b      ,        &
    lhfl_s_b       , pabs_b         , sobs_b         , thbs_b         ,        &
    tch_b          , tcm_b          , tfv_b          , u_10m_b        ,        &
    v_10m_b        ,                                                           &
    w_snow_mult_b      , w_snow_mult_new_b  , t_snow_mult_b      ,             &
    t_snow_mult_new_b  , dzh_snow_mult_b    , dzh_snow_mult_new_b,             &
    rho_snow_mult_b    , rho_snow_mult_new_b

USE data_modelconfig,      ONLY:                                               &
    ie, je, ke, dt, dt2, ke_soil, ke_snow, czmls,                              &
    istartpar, iendpar, jstartpar, jendpar,                                    &
    pollon, pollat, polgam, dlon, dlat, startlon, startlat, czmls

USE data_runcontrol,       ONLY:                                               &
    nnew, nnow, nblock, nproma, nlastproma, lstomata, llake, lseaice,          &
    itype_trvg, itype_tran, l2tls, itype_canopy, itype_evsl, itype_heatcond,   &
    itype_hydbound, itype_mire, itype_root, lmelt, lmelt_var, lmulti_snow,     &
    itype_calendar

USE sfc_terra_data,        ONLY:                                               &
    crsmin, idiag_snowfrac

USE utilities,             ONLY:                                               &
    phirot2phi, rla2rlarot, rlarot2rla, get_utc_date

USE tsa_lmparam,           ONLY:                                               &
    gen_area_indices

USE tsa_gribio,            ONLY:                                               &
    get_grib_info_icon

USE support_datetime,      ONLY:                                               &
    date_delta

USE fxtr_definition,       ONLY:                                               &
    kind_idate

USE data_block_fields,     ONLY:                                               &
       t_sn_b ,    t_sn_new_b ,                                                &
    theta_i_b , theta_i_new_b ,                                                &
    theta_w_b , theta_w_new_b ,                                                &
    theta_a_b , theta_a_new_b ,                                                &
     dzm_sn_b ,  dzm_sn_new_b ,                                                &
      hn_sn_b ,   hn_sn_new_b ,                                                &
     top_sn_b ,  top_sn_new_b                                     
       
USE data_runcontrol, ONLY : lsnow,lsnow_coldstart


!==============================================================================

IMPLICIT NONE

!==============================================================================

CONTAINS

!==============================================================================

!----------------------------------------------------------------------
SUBROUTINE read_namelist
!----------------------------------------------------------------------

!-------------------------------------------------------------------------------
!
! Description:
!   Sets definitions for reading of the namelist of TSA
!
!-------------------------------------------------------------------------------

  IMPLICIT NONE

  ! Local variables for subregion selection and time management

  INTEGER :: &
    i,j    , & ! loop variables
    cnt                                      ! subregions counter
!DL
  ! Local variables for definition of nblock and nlastproma
  ! uuid for ICON horizontal grid
  INTEGER :: nij                                      ! number of grid points
  INTEGER :: isc, jsc, iec, jec, ipend, ib, ip
  INTEGER :: ierr

  LOGICAL                               :: isfound
  CHARACTER(LEN=1)                      :: uuid_in(16)
!DL

!PK>
  INTEGER, ALLOCATABLE :: &
    ie_sub(:,:)                          , & ! index variable for sub region division
    je_sub(:,:)                              ! index variable for sub region division

  REAL (KIND=wp), ALLOCATABLE           :: &
    startlon_sub(:,:)                    , & ! first longitude of sub region
    startlat_sub(:,:)                        ! first latitude of sub region
!PK<
  INTEGER (KIND=kind_idate)             :: &
    initial_date                         , & ! start date of run
    final_date                               ! end date of run
  CHARACTER (LEN=28)                    :: yactdate2   ! date in long format (dummy)
  INTEGER                               :: doy         ! day of year (dummy)
  REAL (KIND=wp)                        :: acthour      ! current hour (dummy)

! Namelists Variables:
!-------------------------
  namelist /RUN_TERRA/                     & ! 
    dt,                                    & ! 
    ntstep_max,                            & ! 
    ke_soil,                               & ! 
    lvegadapt,                             & ! 
    ie,je,                                 & ! 
    ydate_ini,                             & ! initial and final dates of simulation (yyyymmddhh)
    ydate_end,                             & !   ydate_end added following JT (XYZ)
!!DL
    ymodel,                                & ! COSMO or ICON
    ngrid,                                 & ! number of used horizontal grid
    ngpi,                                  & ! generatingProcessIdentifier - defining the ICON domain
    yform_read,                            & ! 'apix' for eccodes, 'grb1' for libDWD (only GRIB1)
    nproma,                                & ! number of grid points per block -> block structure
!!DL
    lrootadapt,                            & ! 
    dlon,                                  & ! 
    dlat,                                  & ! 
    startlon,                              & ! 
    startlat,                              & ! 
    pollon,                                & ! 
    pollat,                                & ! 
    polgam,                                & ! 
!PK>
    nsub_x,                                & ! 
    nsub_y,                                & ! 
    which_subreg,                          & ! 
!PK<
    lcheck,                                & ! 
    lz0local,                              & ! 
    lconstvegalb,                          & ! 

! GDM>
    lstomata,                              & ! 
    itype_evsl,                            & ! 
    itype_hydbound,                        & ! 
    itype_hydparam,                        & ! 
    itype_root,                            & ! 
!DL
    itype_canopy,                          & !
    itype_mire,                            & !
    idiag_snowfrac,                        & !
!DL
    czhls_const,                           & ! 
    nrecmax,                               & ! 
    lmelt,                                 & ! 
    lmelt_var,                             & ! 
    itype_heatcond,                        & ! 
    linfil_revised,                        & ! 
    lmulti_snow,                           & !
    lsnow,                                 & !
    lsnow_coldstart,                       & ! 
    ke_snow,                               & ! 
!   itype_hydcond,                         & ! 
!DL kexpdec,                         
    crsmin,                                & ! 
! GDM<
    lcalc                                    ! 
!-------------------------

  namelist /EXTPARA/                       & ! 
       lat,                                & ! 
       hsurface,                           & ! 
       constfilename,                      & ! 
       lhomosoil,                          & ! 
       soiltyp_const,                      & ! 
       plcov_const,                        & ! 
       rootdp_const,                       & ! 
       lai_const,                          & ! 
       lai_min_const,                      & ! 
       lai_max_const,                      & ! 
       plcov_min_const,                    & ! 
       plcov_max_const,                    & ! 
       z0_const,                           & ! 
       rstom_mn_const,                     & ! 
       rstom_mx_const,                     & ! 
       vegalb_const,                       & !    
       lext_monthly,                       & ! 
       lgettcl                               ! 

!-------------------------
  namelist /SOILINIT/                      &
       lhomoinit,                          & ! 
       t_soil0,                            & ! 
       t_cl0,                              & ! 
       w_snow0,                            & ! 
       w_g0,                               & ! 
       w_i0,                               & ! 
       w_cl0,                              & ! 
       lrel_in,                            & ! 
       soilinitdir,                        & ! 
       soilinitprefix,                     & ! 
       lvol_in,                            & ! 
       lmulti_in                             ! 

!-------------------------
  namelist /METFORCING/                    &        
       ntype_atminput,                     & ! 
       ntype_raininput,                    & ! 
       ntype_radinput,                     & ! 
       metfiledir,                         & ! 
       metfileprefix,                      & ! 
       radofiledir,                        & ! 
       rain_fac,                           & ! 
       lhourly_data,                       & ! 
       tincr_max,                          & ! 
       ldestaggeruv,                       & ! 
       lpar,                               & ! 
       dz,                                 & ! 
       dz_u,                               & ! 
       ke_model                              ! 

!-------------------------
  namelist /OUTPUT/                        &
       ntype_output,                       & ! 
       yform_write,                        & ! DL
       outdir,                             & ! 
       outprefix,                          & ! 
       nout_interval,                      & ! 
       lconstout                             ! 


  ! Setting default values
!-------------------------
  ydate_ini='2010010100'    
  ydate_end='0000000000'
  nrecmax=-1                 
  lcheck=.TRUE.
  lvegadapt=.TRUE.
!DL set lhomoinit and lhomosoil to .FALSE. 
!   set lmulti_in to .TRUE. because this is the standard
! lhomoinit=.TRUE.
! lhomosoil=.TRUE.
! lmulti_in=.FALSE.
  lhomoinit=.FALSE.
  lhomosoil=.FALSE.
  lmulti_in=.TRUE.
  ie=1
  je=1
  dlon=1
  dlat=1  

!PK>  Default values for rotation:
  pollon=-170
  pollat=32.5
  polgam=0.0
  startlon=1
  startlat=1
  nsub_x=1
  nsub_y=1
  which_subreg=1
!PK<
  lat=52
  hsurface=0.0
  lrel_in=.FALSE.
  lvol_in=.FALSE.
  dz=2.0
  dz_u=10.0
  ke_model=-1
!GDM>
  lmulti_snow=.FALSE.
  ke_snow=2
!GDM<
!VS<
  lsnow = .FALSE.
  lsnow_coldstart = .FALSE.
!VS >
  lconstout=.TRUE.
  rain_fac=1.0
  tincr_max=0
  ntype_atminput=1
  ntype_raininput=1
  ntype_radinput=1
  ntype_output=1
  dt=60
  ntstep_max=0
  ke_soil=7
  czhls_const(1:8)=(/ 0.01, 0.03, 0.09, 0.27, 0.81, 2.43, 7.29, 21.87 /)
  nout_interval=60  
!DL>
  yform_read ='apix'
  yform_write='api2'
  ymodel='COSMO'
  ngrid=44
  ngpi=11
  nproma=1
  itype_canopy=1   ! canopy energetically not represented (no skin temp.), 2: skin temp. present
  itype_mire=0     ! no mire, =1: Approach from Alla Yurova et al., 2014
  idiag_snowfrac=1 !
!DL<
!GDM>
  itype_evsl=2     ! bare soil evaporation after BATS
  itype_heatcond=1 ! no soil moisture dependent heat conductivity
  itype_hydbound=1 ! free drainage for lower hydrological boundary cond.
  itype_hydparam=1 ! soil moisture drainage and diffusion parametr.
  itype_root=1     ! uniform root density distribution
!GDM<
!PK>
!  itype_hydcond=0
!DLkexpdec=2.0
  crsmin=150.0
!PK<
  lconstvegalb=.TRUE.
  lstomata=.TRUE.
  lz0local=.FALSE.
  lext_monthly=.FALSE.
  lgettcl=.FALSE.
  lmelt=.FALSE.
  lmelt_var=.FALSE.
  linfil_revised=.FALSE.
  lcalc=.TRUE.
  lrootadapt=.TRUE.
  lhourly_data=.FALSE.  
  ldestaggeruv=.FALSE.  ! per default, don't unstagger wind velocities
  lpar=.FALSE.          ! per default, compute PAR=0.5*SOBS
  radofiledir=''

! Reading Namelist Files
  OPEN(99,FILE='INPUT_TERRA',FORM='FORMATTED')
  READ(99,RUN_TERRA)
  READ(99,EXTPARA)
  READ(99,SOILINIT)
  READ(99,METFORCING)
  READ(99,OUTPUT)
  CLOSE(99)

!PK>
  ! Modify startlon,startlat & ie,je to relate to sub-region
  !   and not entire domain if computation by subregions is selected 
  IF ( nsub_x*nsub_y>1 ) THEN
    ! Allocate temporary matrices to hold information about subregions
    ALLOCATE(ie_sub(nsub_x,nsub_y))
    ALLOCATE(je_sub(nsub_x,nsub_y))
    ALLOCATE(startlon_sub(nsub_x,nsub_y))
    ALLOCATE(startlat_sub(nsub_x,nsub_y))
    ! Compute number of x and y values in each subregion
    ie_sub(:,:) = FLOOR(REAL(ie)/REAL(nsub_x))            !value of each row
    ie_sub(nsub_x,:) = ie - SUM( ie_sub(1:(nsub_x-1),1) ) !value of last row
    je_sub(:,:) = FLOOR(REAL(je)/REAL(nsub_y))            !value of each coloumn
    je_sub(:,nsub_y) = je -  SUM( je_sub(1,1:(nsub_y-1)) )!value of last coloumn
    ! Compute startlon and startlat values for each subregion
    startlon_sub(1,:) = startlon
    DO i=2,nsub_x
      startlon_sub(i,:) = startlon_sub(i-1,1) + ie_sub(i-1,1) * dlon
    ENDDO
    startlat_sub(:,1) = startlat
    DO j=2,nsub_y
      startlat_sub(:,j) = startlat_sub(1,j-1) + je_sub(1,j-1) * dlat
    ENDDO
    IF (which_subreg>nsub_x*nsub_y) THEN
       WRITE(6,*) "No. of subregion too large!"
       STOP 
    ENDIF
    ! Choose subregion and set startlon, startlat, ie, je for the specific sub-region
    cnt=0 !counter for regions
    DO i=1,nsub_x
      DO j=1,nsub_y
        cnt=cnt+1
        IF (cnt==which_subreg) THEN
           startlon=startlon_sub(i,j)
           startlat=startlat_sub(i,j)
           ie=ie_sub(i,j)
           je=je_sub(i,j)
        ENDIF
      ENDDO
    ENDDO
    ! Deallocate temporary matrices
    DEALLOCATE(ie_sub,je_sub,startlon_sub,startlat_sub)
  ENDIF

  ! if end date is given, re-compute ntstep_max from it. caution: this will overwrite ntstep_max!
  READ(ydate_ini,*) initial_date
  READ(ydate_end,*) final_date
  IF (final_date > initial_date) THEN !may overwrite ntstep_max !
    ntstep_max = (date_delta(final_date*100,initial_date)) / ( dt/60 )
  ELSE  !if final date not given, compute it
    CALL get_utc_date(ntstep_max, ydate_ini, dt, itype_calendar, ydate_end, yactdate2, doy, acthour)
  ENDIF

! Screen output of (potentially modified) model configuration
!DL uncommented
  WRITE(*,*) "Model configuration:"
  WRITE(*,*) "----------------------------"
  IF (ymodel == 'COSMO') THEN
    WRITE(*,*) "pollat, pollon, polgam = ", pollat, pollon, polgam
    WRITE(*,*) "Subregion ", which_subreg, " of ", nsub_x*nsub_y, " regions."
    WRITE(*,*) "startlon, ie = ", startlon, ie
    WRITE(*,*) "startlat, je = ", startlat, je
  ELSE
    WRITE(*,*) "Number of ICON grid = ", ngrid
    WRITE(*,*) "Generating Process Identifier (domain) = ", ngpi
    WRITE(*,*) "Subregion ", which_subreg, " of ", nsub_x*nsub_y, " regions."
  ENDIF
  WRITE(*,*) "Initial and final date: ", ydate_ini, " - ", ydate_end
  WRITE(*,*) "Time step and ntstep_max: ", dt, ntstep_max
  WRITE(*,*) "----------------------------"
  WRITE(*,*)
!PK<
!DL uncommented

!!DL##############################################################################
! Some constraints for ICON, Get number of grid points for ICON
! Define block structure (nblock and nlastproma)
  IF (TRIM(ymodel) == 'ICON') THEN
    yform_read='apix'
    yform_write='api2'
    CALL get_grib_info_icon(constfilename,ngrid,ngpi,nij,uuid_in,isfound)
    IF (isfound) THEN
      ! Set uuid of horizontal ICON grid for later checks (defined in data_terra_standalone)
      uuid = uuid_in
      ! Set number of grid points for later use
      ngp  = nij
    ELSE
      WRITE (6,*) 'File with external parameter not found:',constfilename
      STOP 'No extpar file!!'
    ENDIF
  ELSE
    nij = ie * je
  ENDIF

  PRINT *, "nij, ngrid, ngpi, uuid :", nij, ngrid, ngpi !, uuid

! Compute nblock, nlastproma for COSMO and ICON
  IF (MOD(nij,nproma) == 0) THEN
    nblock     = INT (nij/nproma)
    nlastproma = nproma
  ELSE
    nblock     = INT(nij/nproma) + 1
    nlastproma = MOD(nij,nproma)
  ENDIF

  PRINT *, "nblock, nlastproma: ", nblock, nlastproma


! From src_block_fields
  ! Allocate and initialize mind_ilon, mind_jlat
  ! compute the index-arrays

  ! Allocation
  ALLOCATE (mind_ilon(nproma, nblock), mind_jlat(nproma, nblock), STAT=ierr)
  IF (ierr>0) THEN
    WRITE(*,*) 'ERROR: Index arrays allocation failed'
    STOP ' Allocation of mind_ilon,mind_jlat failed!'
  END IF

  IF (ymodel == 'COSMO') THEN
    ! Set index
    isc = 1
    jsc = 1
    iec = ie
    jec = je
    i = isc - 1
    j = jsc

    ipend = nproma       ! compute domain for all but last block is 1:nproma
    DO ib = 1, nblock
      IF (ib == nblock) ipend = nlastproma          ! last block 1:nlastproma
      DO ip = 1, ipend
        i = i + 1
        IF (i > iec) THEN
          j = j + 1
          i = isc
        ENDIF
        IF (j > jec) THEN
          ! Error: this must not happen, because
          ! nproma*nblock <= (iec-isc+1)*(jec-jsc+1)
          WRITE(*,*) ' ERROR: inconsistent sizes nproma*nblock, nlastproma, nij :', &
               nproma*nblock, nlastproma, nij
          STOP
        ENDIF
  
        mind_ilon (ip,ib) = i
        mind_jlat (ip,ib) = j
      ENDDO
    ENDDO

  ELSE     ! ICON
    ipend = nproma
    DO ib = 1, nblock
      IF (ib == nblock) ipend = nlastproma          ! last block 1:nlastproma
      DO ip = 1, ipend
        mind_ilon (ip,ib) = ip
        mind_jlat (ip,ib) = ib
      ENDDO
    ENDDO

  ENDIF    ! COSMO model
  
!!DL##############################################################################

  ! set multi_layer scheme to be TRUE
! lmulti_layer=.TRUE.
!GDM>
  IF (lstomata) THEN
    WRITE(*,*) "lstomata is now set to .TRUE. => rsmin2d in ALLOCATE_fields"
  ENDIF
!GDM<
  ke=1
  ! WRITE YUSPECIF files containing all parameters.
  ! this is useful to compare runs as not all parameters
  ! need necessarily be in the namelists and default values
  ! can change with time.
  OPEN(99,FILE='YUSPECIF',FORM='FORMATTED')
  WRITE(99,RUN_TERRA)
  WRITE(99,EXTPARA)
  WRITE(99,SOILINIT)
  WRITE(99,METFORCING)
  WRITE(99,OUTPUT)
  CLOSE(99)

END SUBROUTINE read_namelist
!----------------------------------------------------------------------




SUBROUTINE allocate_fields
!-------------------------------------------------------------------------------
! Description:
!   Allocates variables for TSA
!-------------------------------------------------------------------------------

  ! Allocate COSMO fields
  INTEGER, PARAMETER :: ntlev=2
  INTEGER            :: i, j

  IF(ymodel == 'COSMO') THEN
    i = ie
    j = je
  ELSE
    i = nproma
    j = nblock
  ENDIF

  ALLOCATE(                   &
       p0(i,j,ke:ke)      , &   
       soiltyp(i,j)       , &
       lai(i,j)           , &
       plcov(i,j)         , &   
       rootdp(i,j)        , &   
       sai(i,j)           , &   
       tai(i,j)           , &   
       eai(i,j)           , &   
       llandmask(i,j)     , &   
       u(i,j,ke:ke,ntlev) , &   
       v(i,j,ke:ke,ntlev) , &   
       t(i,j,ke:ke,ntlev) , &   
       qv(i,j,ke:ke,ntlev), &
!      qv(i,j,ke:ke)      , &   
       pp(i,j,ke:ke,ntlev), &   
       ps(i,j,ntlev)      , &       
       t_snow(i,j,ntlev)  , &
       t_s(i,j,ntlev)     , &   
!DL
       t_sk(i,j,ntlev)    , &
       skinc(i,j)         , &
       sso_sigma(i,j)     , &
!DL
       t_g(i,j,ntlev)     , &    
       u_10m(i,j)         , &    
       v_10m(i,j)         , &    
       t_2m(i,j)          , &    
       qv_s(i,j,ntlev)    , &    
       w_snow(i,j,ntlev)  , &   
       w_i(i,j,ntlev)     , &       
!DL
       w_p(i,j,ntlev)     , &
       w_s(i,j,ntlev)     , &
!DL
       prr_con(i,j)       , &    
       prs_con(i,j)       , &    
       prr_gsp(i,j)       , &    
       prs_gsp(i,j)       , &    
!DL
       prg_gsp(i,j)       , &
!DL
       tch(i,j)           , &     
       tcm(i,j)           , &     
       tfh (i,j)          , &    
       tfv (i,j)          , &    
       sobs(i,j)          , &     
       thbs(i,j)          , &     
       pabs(i,j)          , &     
       runoff_s(i,j)      , &   
       runoff_g(i,j)      , &
       lhfl_s(i,j)        , &
       shfl_s(i,j)        , &
       lhfl_bs(i,j)       , &
       qvfl_s(i,j),         & ! surface flux of water vapour                  (kg/m2s)
       h_snow(i,j,ntlev)  , &
       rho_snow(i,j,ntlev) )

  ALLOCATE(                   &
       u_bd(i,j,ke:ke,2)  , & ! from here meteo_fields
       t_bd(i,j,ke:ke,2)  , &
       qv_bd(i,j,ke:ke,2) , &
!      qv_bd(i,j,ke:ke)   , &
       ps_bd(i,j,2)       , &
       prr_gsp_bd(i,j,2)  , &
       prs_gsp_bd(i,j,2)  , &
       so_down_bd(i,j,2)  , &
       th_down_bd(i,j,2)  , &
       pabs_bd(i,j,2) )

  ALLOCATE(                             &
          t_so(i,j,0:ke_soil+1,ntlev)  , &
          w_so(i,j,ke_soil+1,ntlev)    , &
          w_so_ice(i,j,ke_soil+1,ntlev), &
          t_cl(i,j) )
  ALLOCATE(                              &
           wliq_snow(i,j,ke_snow,ntlev)    , &
           w_snow_mult(i,j,ke_snow,ntlev)  , &
           dzh_snow_mult(i,j,ke_snow,ntlev), &
           t_snow_mult(i,j,ke_snow,ntlev)  , &
           rho_snow_mult(i,j,ke_snow,ntlev) )


  if ( lsnow ) then
     allocate( theta_i(i,j,ke_snow,ntlev) , &
               theta_w(i,j,ke_snow,ntlev) , &
               theta_a(i,j,ke_snow,ntlev) , &
               dzm_sn(i,j,ke_snow,ntlev) , &
               t_sn(i,j,ke_snow,ntlev) , &
               top_sn(i,j,ntlev) )

  endif

     



  ALLOCATE(rsmin2d(i,j))
  ALLOCATE( lhfl_pl(i,j,ke_soil), &
               rstom(i,j)          , &
               freshsnow(i,j) )
  ALLOCATE(czmls(ke_soil+1))

  ALLOCATE(                  &
       plcov_mx(i,j)     , &
       plcov_mn(i,j)     , &
       rootdp0(i,j)      , &
       lai_mn(i,j)       , &
       lai_mx(i,j)       , &
       z0(i,j)           , &
       rlon(i,j)         , &
       rlat(i,j)         , &
       hsurf(i,j)        , &
       fr_wi(i,j)        , &
       fr_snow(i,j)      , &
       fr_land(i,j)      , &
       snow_melt(i,j)    , &
! XYZ> adding geographical lat and long following JT of GUF
       rlon_geo(i,j)     , &
       rlat_geo(i,j) )
!<XYZ

  ALLOCATE(                   &
       astore(i,j,nastore), &
       store(i,j,nstore) )
  astore=0.0
  ALLOCATE(runoff(i,j,ke_soil+1))
  runoff=0.0
  ALLOCATE(zflmg(i,j,ke_soil+1))

  ALLOCATE(                   &
       vegalb(i,j)        , &
       rstom_mn(i,j)      , &
       rstom_mx(i,j) )

END SUBROUTINE allocate_fields

!-------------------------------------------------------------------------------

SUBROUTINE allocate_block_fields 

!-------------------------------------------------------------------------------
! Description:
!   Allocates block variables for TSA  
!   - modification of "block_fields_allocate" (src_block_fields_org.f90)
!-------------------------------------------------------------------------------

! Locals:
! -------
  INTEGER :: ist
  INTEGER :: izl
  REAL(KIND=wp), PARAMETER :: &
 
    r_init_val=-99999999_wp   ! It is an error to use a block field
                              ! without assigning it a value before.
                              ! All fields are therefore initialized to
                              ! non realistic value
  INTEGER      , PARAMETER :: &
    i_init_val=-99999999      ! non-realistic integer values#
!-------------------------------------------------------------------------------

  ist=0
  izl=0

PRINT *, 'ALLOCATE_BLOCK: nblock, nproma:', nblock, nproma

  ! Allocate COSMO/ICON  block fields

  ! constant fields for the reference atmosphere
  ALLOCATE(p0_b           (nproma,ke:ke)  , STAT=izl); p0_b           = r_init_val; ist=ist+izl

  ! external parameters, needed in several parameterizations
  ALLOCATE(hsurf_b        (nproma)     , STAT=izl); hsurf_b        = r_init_val; ist=ist+izl
  ALLOCATE(fr_land_b      (nproma)     , STAT=izl); fr_land_b      = r_init_val; ist=ist+izl
  ALLOCATE(z0_b           (nproma)     , STAT=izl); z0_b           = r_init_val; ist=ist+izl
! ALLOCATE(llandmask_b    (nproma)     , STAT=izl); llandmask_b    = .FALSE.   ; ist=ist+izl
  ALLOCATE(isoiltyp_b     (nproma)     , STAT=izl); isoiltyp_b     = i_init_val; ist=ist+izl
  ALLOCATE(plcov_b        (nproma)     , STAT=izl); plcov_b        = r_init_val; ist=ist+izl
  ALLOCATE(rootdp_b       (nproma)     , STAT=izl); rootdp_b       = r_init_val; ist=ist+izl
  ALLOCATE(sai_b          (nproma)     , STAT=izl); sai_b          = r_init_val; ist=ist+izl
  ALLOCATE(tai_b          (nproma)     , STAT=izl); tai_b          = r_init_val; ist=ist+izl
  ALLOCATE(eai_b          (nproma)     , STAT=izl); eai_b          = r_init_val; ist=ist+izl
  ALLOCATE(skinc_b        (nproma)     , STAT=izl); skinc_b        = r_init_val; ist=ist+izl
  ALLOCATE(rsmin2d_b      (nproma)     , STAT=izl); rsmin2d_b      = r_init_val; ist=ist+izl
  ALLOCATE(sso_sigma_b    (nproma)     , STAT=izl); sso_sigma_b    = r_init_val; ist=ist+izl
  
  ! dynamical variables
  ALLOCATE(ptot_b         (nproma,ke:ke)  , STAT=izl); ptot_b         = r_init_val; ist=ist+izl
! ALLOCATE(epr_b          (nproma,ke:ke)  , STAT=izl); epr_b          = r_init_val; ist=ist+izl
  ALLOCATE(pp_b           (nproma,ke:ke)  , STAT=izl); pp_b           = r_init_val; ist=ist+izl
  ALLOCATE(t_b            (nproma,ke:ke)  , STAT=izl); t_b            = r_init_val; ist=ist+izl
  ALLOCATE(qv_b           (nproma,ke:ke)  , STAT=izl); qv_b           = r_init_val; ist=ist+izl
  ALLOCATE(u_m_b          (nproma,ke:ke)  , STAT=izl); u_m_b          = r_init_val; ist=ist+izl
  ALLOCATE(v_m_b          (nproma,ke:ke)  , STAT=izl); v_m_b          = r_init_val; ist=ist+izl

  ! fields for surface values
  ALLOCATE(ps_b           (nproma)     , STAT=izl); ps_b           = r_init_val; ist=ist+izl
  ALLOCATE(t_s_b          (nproma)     , STAT=izl); t_s_b          = r_init_val; ist=ist+izl
  ALLOCATE(t_s_new_b      (nproma)     , STAT=izl); t_s_new_b      = r_init_val; ist=ist+izl
  ALLOCATE(t_sk_b         (nproma)     , STAT=izl); t_sk_b          = r_init_val; ist=ist+izl
  ALLOCATE(t_sk_new_b     (nproma)     , STAT=izl); t_sk_new_b      = r_init_val; ist=ist+izl
  ALLOCATE(t_g_b          (nproma)     , STAT=izl); t_g_b          = r_init_val; ist=ist+izl
  ALLOCATE(t_g_new_b      (nproma)     , STAT=izl); t_g_new_b      = r_init_val; ist=ist+izl
  ALLOCATE(qv_s_b         (nproma)     , STAT=izl); qv_s_b         = r_init_val; ist=ist+izl
  ALLOCATE(qv_s_new_b     (nproma)     , STAT=izl); qv_s_new_b     = r_init_val; ist=ist+izl

  ! from microphysics
! ALLOCATE(tinc_lh_b      (nproma,ke:ke)  , STAT=izl); tinc_lh_b      = r_init_val; ist=ist+izl
! ALLOCATE(qrs_b          (nproma,ke:ke)  , STAT=izl); qrs_b          = r_init_val; ist=ist+izl
  ALLOCATE(prr_con_b      (nproma)     , STAT=izl); prr_con_b      = r_init_val; ist=ist+izl
  ALLOCATE(prs_con_b      (nproma)     , STAT=izl); prs_con_b      = r_init_val; ist=ist+izl
  ALLOCATE(prr_gsp_b      (nproma)     , STAT=izl); prr_gsp_b      = r_init_val; ist=ist+izl
  ALLOCATE(prs_gsp_b      (nproma)     , STAT=izl); prs_gsp_b      = r_init_val; ist=ist+izl
  ALLOCATE(prg_gsp_b      (nproma)     , STAT=izl); prg_gsp_b      = r_init_val; ist=ist+izl
  ALLOCATE(prh_gsp_b      (nproma)     , STAT=izl); prh_gsp_b      = r_init_val; ist=ist+izl

  ! from radiation scheme (needed in surface TERRA)
  ALLOCATE(sobs_b         (nproma)     , STAT=izl); sobs_b         = r_init_val; ist=ist+izl
  ALLOCATE(thbs_b         (nproma)     , STAT=izl); thbs_b         = r_init_val; ist=ist+izl
  ALLOCATE(pabs_b         (nproma)     , STAT=izl); pabs_b         = r_init_val; ist=ist+izl

  ! additional variables
! ALLOCATE(t_2m_b         (nproma)     , STAT=izl); t_2m_b         = r_init_val; ist=ist+izl
! ALLOCATE(qv_2m_b        (nproma)     , STAT=izl); qv_2m_b        = r_init_val; ist=ist+izl
! ALLOCATE(td_2m_b        (nproma)     , STAT=izl); td_2m_b        = r_init_val; ist=ist+izl
! ALLOCATE(rh_2m_b        (nproma)     , STAT=izl); rh_2m_b        = r_init_val; ist=ist+izl
  ALLOCATE(u_10m_b        (nproma)     , STAT=izl); u_10m_b        = r_init_val; ist=ist+izl
  ALLOCATE(v_10m_b        (nproma)     , STAT=izl); v_10m_b        = r_init_val; ist=ist+izl
  ALLOCATE(shfl_s_b       (nproma)     , STAT=izl); shfl_s_b       = r_init_val; ist=ist+izl
  ALLOCATE(lhfl_s_b       (nproma)     , STAT=izl); lhfl_s_b       = r_init_val; ist=ist+izl
  ALLOCATE(qvfl_s_b       (nproma)     , STAT=izl); qvfl_s_b       = r_init_val; ist=ist+izl

  ! from the surface schemes
  ALLOCATE(t_snow_b       (nproma)     , STAT=izl); t_snow_b       = r_init_val; ist=ist+izl
  ALLOCATE(t_snow_new_b   (nproma)     , STAT=izl); t_snow_new_b   = r_init_val; ist=ist+izl
  ALLOCATE(w_snow_b       (nproma)     , STAT=izl); w_snow_b       = r_init_val; ist=ist+izl
  ALLOCATE(w_snow_new_b   (nproma)     , STAT=izl); w_snow_new_b   = r_init_val; ist=ist+izl
  ALLOCATE(h_snow_b       (nproma)     , STAT=izl); h_snow_b       = r_init_val; ist=ist+izl
  ALLOCATE(h_snow_new_b   (nproma)     , STAT=izl); h_snow_new_b   = r_init_val; ist=ist+izl
  ALLOCATE(rho_snow_b     (nproma)     , STAT=izl); rho_snow_b     = r_init_val; ist=ist+izl
  ALLOCATE(rho_snow_new_b (nproma)     , STAT=izl); rho_snow_new_b = r_init_val; ist=ist+izl
  ALLOCATE(freshsnow_b    (nproma)     , STAT=izl); freshsnow_b    = r_init_val; ist=ist+izl
  ALLOCATE(fr_snow_b      (nproma)     , STAT=izl); fr_snow_b      = r_init_val; ist=ist+izl
  ALLOCATE(w_i_b          (nproma)     , STAT=izl); w_i_b          = r_init_val; ist=ist+izl
  ALLOCATE(w_i_new_b      (nproma)     , STAT=izl); w_i_new_b      = r_init_val; ist=ist+izl
  ALLOCATE(w_p_b          (nproma)     , STAT=izl); w_p_b          = r_init_val; ist=ist+izl
  ALLOCATE(w_p_new_b      (nproma)     , STAT=izl); w_p_new_b      = r_init_val; ist=ist+izl
  ALLOCATE(w_s_b          (nproma)     , STAT=izl); w_s_b          = r_init_val; ist=ist+izl
  ALLOCATE(w_s_new_b      (nproma)     , STAT=izl); w_s_new_b      = r_init_val; ist=ist+izl

  ALLOCATE(t_so_b         (nproma,0:ke_soil+1), STAT=izl); t_so_b         = r_init_val; ist=ist+izl
  ALLOCATE(t_so_new_b     (nproma,0:ke_soil+1), STAT=izl); t_so_new_b     = r_init_val; ist=ist+izl
  ALLOCATE(w_so_b         (nproma,  ke_soil+1), STAT=izl); w_so_b         = r_init_val; ist=ist+izl
  ALLOCATE(w_so_new_b     (nproma,  ke_soil+1), STAT=izl); w_so_new_b     = r_init_val; ist=ist+izl
  ALLOCATE(w_so_ice_b     (nproma,  ke_soil+1), STAT=izl); w_so_ice_b     = r_init_val; ist=ist+izl
  ALLOCATE(w_so_ice_new_b (nproma,  ke_soil+1), STAT=izl); w_so_ice_new_b = r_init_val; ist=ist+izl

  ALLOCATE(t_snow_mult_b      (nproma,0:ke_snow), STAT=izl); t_snow_mult_b      = r_init_val; ist=ist+izl
  ALLOCATE(t_snow_mult_new_b  (nproma,0:ke_snow), STAT=izl); t_snow_mult_new_b  = r_init_val; ist=ist+izl
  ALLOCATE(w_snow_mult_b      (nproma,  ke_snow), STAT=izl); w_snow_mult_b      = r_init_val; ist=ist+izl
  ALLOCATE(w_snow_mult_new_b  (nproma,  ke_snow), STAT=izl); w_snow_mult_new_b  = r_init_val; ist=ist+izl
  ALLOCATE(wliq_snow_b        (nproma,  ke_snow), STAT=izl); wliq_snow_b        = r_init_val; ist=ist+izl
  ALLOCATE(wliq_snow_new_b    (nproma,  ke_snow), STAT=izl); wliq_snow_new_b    = r_init_val; ist=ist+izl
  ALLOCATE(rho_snow_mult_b    (nproma,  ke_snow), STAT=izl); rho_snow_mult_b    = r_init_val; ist=ist+izl
  ALLOCATE(rho_snow_mult_new_b(nproma,  ke_snow), STAT=izl); rho_snow_mult_new_b= r_init_val; ist=ist+izl
  ALLOCATE(dzh_snow_mult_b    (nproma,  ke_snow), STAT=izl); dzh_snow_mult_b    = r_init_val; ist=ist+izl
  ALLOCATE(dzh_snow_mult_new_b(nproma,  ke_snow), STAT=izl); dzh_snow_mult_new_b= r_init_val; ist=ist+izl

  ALLOCATE(runoff_s_b     (nproma)     , STAT=izl); runoff_s_b     = r_init_val; ist=ist+izl
  ALLOCATE(runoff_g_b     (nproma)     , STAT=izl); runoff_g_b     = r_init_val; ist=ist+izl
  ALLOCATE(lhfl_bs_b      (nproma)     , STAT=izl); lhfl_bs_b      = r_init_val; ist=ist+izl
  ALLOCATE(lhfl_pl_b   (nproma,ke_soil), STAT=izl); lhfl_pl_b      = r_init_val; ist=ist+izl
  ALLOCATE(rstom_b        (nproma)     , STAT=izl); rstom_b        = r_init_val; ist=ist+izl
  ALLOCATE(snow_melt_b    (nproma)     , STAT=izl); snow_melt_b    = r_init_val; ist=ist+izl

  ! from turbulence scheme
  ALLOCATE(tcm_b          (nproma)     , STAT=izl); tcm_b          = r_init_val; ist=ist+izl
  ALLOCATE(tch_b          (nproma)     , STAT=izl); tch_b          = r_init_val; ist=ist+izl
  ALLOCATE(tfv_b          (nproma)     , STAT=izl); tfv_b          = r_init_val; ist=ist+izl

  !VS < 
  ! for the multi-layer snow cover scheme (SNOWPOLINO)
  if( lsnow ) then
     ALLOCATE(t_sn_b         (nproma,1:ke_snow), STAT=izl); t_sn_b            = r_init_val; ist=ist+izl
     ALLOCATE(t_sn_new_b     (nproma,1:ke_snow), STAT=izl); t_sn_new_b        = r_init_val; ist=ist+izl
   
     ALLOCATE(theta_i_b      (nproma,1:ke_snow), STAT=izl); theta_i_b         = r_init_val; ist=ist+izl
     ALLOCATE(theta_i_new_b  (nproma,1:ke_snow), STAT=izl); theta_i_new_b     = r_init_val; ist=ist+izl
   
     ALLOCATE(theta_w_b      (nproma,1:ke_snow), STAT=izl); theta_w_b         = r_init_val; ist=ist+izl
     ALLOCATE(theta_w_new_b  (nproma,1:ke_snow), STAT=izl); theta_w_new_b     = r_init_val; ist=ist+izl
   
     ALLOCATE(theta_a_b      (nproma,1:ke_snow), STAT=izl); theta_a_b         = r_init_val; ist=ist+izl
     ALLOCATE(theta_a_new_b  (nproma,1:ke_snow), STAT=izl); theta_a_new_b     = r_init_val; ist=ist+izl
   
     ALLOCATE(dzm_sn_b       (nproma,1:ke_snow), STAT=izl); dzm_sn_b          = r_init_val; ist=ist+izl
     ALLOCATE(dzm_sn_new_b   (nproma,1:ke_snow), STAT=izl); dzm_sn_new_b      = r_init_val; ist=ist+izl
   
     ALLOCATE(hn_sn_b        (nproma           ), STAT=izl); hn_sn_b           = r_init_val; ist=ist+izl
     ALLOCATE(hn_sn_new_b    (nproma           ), STAT=izl); hn_sn_new_b       = r_init_val; ist=ist+izl
   
     ALLOCATE(top_sn_b       (nproma           ), STAT=izl); top_sn_b          = r_init_val; ist=ist+izl
     ALLOCATE(top_sn_new_b   (nproma           ), STAT=izl); top_sn_new_b      = r_init_val; ist=ist+izl
  endif
!VS >


  IF (ist /= 0) THEN
     WRITE (*,*) ist, " fields in block structure could not be allocated!! STOP!! "
     STOP 'ALLOCATING block fields'
  ENDIF


END SUBROUTINE allocate_block_fields

!-------------------------------------------------------------------------------

SUBROUTINE init_variables()
!-------------------------------------------------------------------------------
! Description:
!   initialize variables for TSA
!-------------------------------------------------------------------------------

  ! initiate COSMO fields and values
  INTEGER :: i,j,k, istart, jstart

  istart=1
  istartpar=1
  iendpar=ie
  IF ( ymodel /= 'COSMO') iendpar=nproma
  jstart=1
  jstartpar=1
  jendpar=je
  IF ( ymodel /= 'COSMO') jendpar=nblock

  dt2=2*dt 
  
  ke=1
  nnow=1
  nnew=2
!GDM>
  itype_trvg=2 ! BATS version of vegetation transpiration
  itype_tran=1 ! new surface-atmosphere transfer
  lstomata=.FALSE.
!GDM< END
!DL
  llake=.FALSE.
  lseaice=.FALSE.
!DL
  l2tls=.TRUE.
  tfh=1.0
  tfv=1.0
!US
  itype_calendar = 0   ! Gregorian Calendar


!DL Next part only for COSMO but needed for tsa_sfc_interface ...
  IF (ymodel == 'COSMO') THEN
! set new rlat & rlon coordinates:
  DO i=1,ie  
     rlon(i,:)=startlon+(i-1)*dlon
  ENDDO
  DO j=1,je
     rlat(:,j)=startlat+(j-1)*dlat
  ENDDO

!XYZ>
! Compute true lon/lat coordinate values following JT of GUF
  DO i = 1 , ie
    DO j = 1 , je
      rlon_geo(i,j) = rlarot2rla( rlat(i,j), rlon(i,j), pollat, pollon, polgam )
      rlat_geo(i,j) = phirot2phi( rlat(i,j), rlon(i,j), pollat, pollon, polgam )
    ENDDO
  ENDDO
!XYZ<
  ENDIF            ! COSMO model
!!DL

!XYZ> cancelled subroutine constants.
!      calling subroutine set_constants from data_constants instead.
!!DL    CALL set_constants
!XYZ<

  czmls(1)=0.5*czhls_const(1)
  DO k=2,ke_soil+1
    czmls(k)=0.5*(czhls_const(k-1)+czhls_const(k))
  ENDDO

  lana_rho_snow=.TRUE.

  !CALL gen_area_indices

  ! superfluous fields
  v=0.0
  pp=0
  prr_con=0.0 
  prs_con=0.0
  prs_gsp=0.0    
!DL ?? prr_gsp=0.0, not necessary as prr_gsp will ALWAYS be defined 
  w_p=0.0
  w_s=0.0
!DL

  runoff_s=0.0
  runoff_g=0.0

END SUBROUTINE init_variables

!-------------------------------------------------------------------------------

SUBROUTINE clean_up

!-------------------------------------------------------------------------------
! Description:
!   de-allocates variables required for TSA
!-------------------------------------------------------------------------------

  INTEGER :: ist, izl

  ist=0
  izl=0

  ! Deallocate COSMO fields
  DEALLOCATE(                                     &
       p0, soiltyp, plcov, rootdp, sai, tai, eai, &   
       llandmask, u, v, t, qv, pp, ps           , &
       t_snow, t_s, t_g, qv_s                   , &
       t_sk, skinc, sso_sigma                   , &
       w_snow, w_i                              , &
       prr_con, prs_con, prr_gsp, prs_gsp       , &    
       tch, tcm, tfh , tfv                      , &    
       sobs, thbs, pabs                         , &
       runoff_s, runoff_g,lai, lhfl_bs          , &
       lhfl_s,shfl_s,t_2m,u_10m,v_10m           , &
       h_snow,rho_snow )

  DEALLOCATE(                                     &
       t_so,w_so,czmls,w_so_ice )

  DEALLOCATE(                                     &
       wliq_snow,w_snow_mult,dzh_snow_mult   ,    &
       t_snow_mult,rho_snow_mult )
  DEALLOCATE(rsmin2d)

  DEALLOCATE(lhfl_pl,rstom,freshsnow)
  DEALLOCATE(astore,store)
  DEALLOCATE(runoff)
  DEALLOCATE(zflmg)
  DEALLOCATE(vegalb,rstom_mn,rstom_mx)

!DL
  ! Deallocate also the block fields
  ! constant fields for the reference atmosphere
  DEALLOCATE(p0_b           , STAT=izl); ist=ist+izl

  ! external parameters, needed in several parameterizations
  DEALLOCATE(hsurf_b        , STAT=izl); ist=ist+izl
  DEALLOCATE(fr_land_b      , STAT=izl); ist=ist+izl
! DEALLOCATE(llandmask_b    , STAT=izl); ist=ist+izl
  DEALLOCATE(isoiltyp_b     , STAT=izl); ist=ist+izl
  DEALLOCATE(plcov_b        , STAT=izl); ist=ist+izl
  DEALLOCATE(rootdp_b       , STAT=izl); ist=ist+izl
  DEALLOCATE(sai_b          , STAT=izl); ist=ist+izl
  DEALLOCATE(tai_b          , STAT=izl); ist=ist+izl
  DEALLOCATE(eai_b          , STAT=izl); ist=ist+izl
  DEALLOCATE(skinc_b        , STAT=izl); ist=ist+izl
  DEALLOCATE(rsmin2d_b      , STAT=izl); ist=ist+izl
  DEALLOCATE(sso_sigma_b    , STAT=izl); ist=ist+izl

  ! dynamical variables
  DEALLOCATE(ptot_b         , STAT=izl); ist=ist+izl
! DEALLOCATE(epr_b          , STAT=izl); ist=ist+izl
  DEALLOCATE(pp_b           , STAT=izl); ist=ist+izl
  DEALLOCATE(t_b            , STAT=izl); ist=ist+izl
  DEALLOCATE(qv_b           , STAT=izl); ist=ist+izl
  DEALLOCATE(u_m_b          , STAT=izl); ist=ist+izl
  DEALLOCATE(v_m_b          , STAT=izl); ist=ist+izl

  ! fields for surface values
  DEALLOCATE(ps_b           , STAT=izl); ist=ist+izl
  DEALLOCATE(t_s_b          , STAT=izl); ist=ist+izl
  DEALLOCATE(t_s_new_b      , STAT=izl); ist=ist+izl
  DEALLOCATE(t_sk_b         , STAT=izl); ist=ist+izl
  DEALLOCATE(t_sk_new_b     , STAT=izl); ist=ist+izl
  DEALLOCATE(t_g_b          , STAT=izl); ist=ist+izl
  DEALLOCATE(t_g_new_b      , STAT=izl); ist=ist+izl
  DEALLOCATE(qv_s_b         , STAT=izl); ist=ist+izl
  DEALLOCATE(qv_s_new_b     , STAT=izl); ist=ist+izl

  ! from microphysics
! DEALLOCATE(tinc_lh_b      , STAT=izl); ist=ist+izl
! DEALLOCATE(qrs_b          , STAT=izl); ist=ist+izl
  DEALLOCATE(prr_con_b      , STAT=izl); ist=ist+izl
  DEALLOCATE(prs_con_b      , STAT=izl); ist=ist+izl
  DEALLOCATE(prr_gsp_b      , STAT=izl); ist=ist+izl
  DEALLOCATE(prs_gsp_b      , STAT=izl); ist=ist+izl
  DEALLOCATE(prg_gsp_b      , STAT=izl); ist=ist+izl
!!  DEALLOCATE(prh_gsp_b      , STAT=izl); ist=ist+izl

  ! from radiation scheme (needed in surface TERRA)
  DEALLOCATE(sobs_b         , STAT=izl); ist=ist+izl
  DEALLOCATE(thbs_b         , STAT=izl); ist=ist+izl
  DEALLOCATE(pabs_b         , STAT=izl); ist=ist+izl

  ! additional variables
! DEALLOCATE(t_2m_b         , STAT=izl); ist=ist+izl
! DEALLOCATE(qv_2m_b        , STAT=izl); ist=ist+izl
! DEALLOCATE(td_2m_b        , STAT=izl); ist=ist+izl
! DEALLOCATE(rh_2m_b        , STAT=izl); ist=ist+izl
  DEALLOCATE(u_10m_b        , STAT=izl); ist=ist+izl
  DEALLOCATE(v_10m_b        , STAT=izl); ist=ist+izl
  DEALLOCATE(shfl_s_b       , STAT=izl); ist=ist+izl
  DEALLOCATE(lhfl_s_b       , STAT=izl); ist=ist+izl
  DEALLOCATE(qvfl_s_b       , STAT=izl); ist=ist+izl

  ! from the surface schemes
  DEALLOCATE(t_snow_b       , STAT=izl); ist=ist+izl
  DEALLOCATE(t_snow_new_b   , STAT=izl); ist=ist+izl
  DEALLOCATE(w_snow_b       , STAT=izl); ist=ist+izl
  DEALLOCATE(w_snow_new_b   , STAT=izl); ist=ist+izl
  DEALLOCATE(h_snow_b       , STAT=izl); ist=ist+izl
  DEALLOCATE(h_snow_new_b   , STAT=izl); ist=ist+izl
  DEALLOCATE(rho_snow_b     , STAT=izl); ist=ist+izl
  DEALLOCATE(rho_snow_new_b , STAT=izl); ist=ist+izl
  DEALLOCATE(freshsnow_b    , STAT=izl); ist=ist+izl
  DEALLOCATE(fr_snow_b      , STAT=izl); ist=ist+izl
  DEALLOCATE(w_i_b          , STAT=izl); ist=ist+izl
  DEALLOCATE(w_i_new_b      , STAT=izl); ist=ist+izl
  DEALLOCATE(w_p_b          , STAT=izl); ist=ist+izl
  DEALLOCATE(w_p_new_b      , STAT=izl); ist=ist+izl
  DEALLOCATE(w_s_b          , STAT=izl); ist=ist+izl
  DEALLOCATE(w_s_new_b      , STAT=izl); ist=ist+izl

  DEALLOCATE(t_so_b         , STAT=izl); ist=ist+izl
  DEALLOCATE(t_so_new_b     , STAT=izl); ist=ist+izl
  DEALLOCATE(w_so_b         , STAT=izl); ist=ist+izl
  DEALLOCATE(w_so_new_b     , STAT=izl); ist=ist+izl
  DEALLOCATE(w_so_ice_b     , STAT=izl); ist=ist+izl
  DEALLOCATE(w_so_ice_new_b , STAT=izl); ist=ist+izl

  DEALLOCATE(t_snow_mult_b      , STAT=izl); ist=ist+izl
  DEALLOCATE(t_snow_mult_new_b  , STAT=izl); ist=ist+izl
  DEALLOCATE(w_snow_mult_b      , STAT=izl); ist=ist+izl
  DEALLOCATE(w_snow_mult_new_b  , STAT=izl); ist=ist+izl
  DEALLOCATE(wliq_snow_b        , STAT=izl); ist=ist+izl
  DEALLOCATE(wliq_snow_new_b    , STAT=izl); ist=ist+izl
  DEALLOCATE(rho_snow_mult_b    , STAT=izl); ist=ist+izl
  DEALLOCATE(rho_snow_mult_new_b, STAT=izl); ist=ist+izl
  DEALLOCATE(dzh_snow_mult_b    , STAT=izl); ist=ist+izl
  DEALLOCATE(dzh_snow_mult_new_b, STAT=izl); ist=ist+izl

  DEALLOCATE(runoff_s_b     , STAT=izl); ist=ist+izl
  DEALLOCATE(runoff_g_b     , STAT=izl); ist=ist+izl
  DEALLOCATE(lhfl_bs_b      , STAT=izl); ist=ist+izl
  DEALLOCATE(lhfl_pl_b      , STAT=izl); ist=ist+izl
  DEALLOCATE(rstom_b        , STAT=izl); ist=ist+izl
  DEALLOCATE(snow_melt_b    , STAT=izl); ist=ist+izl

  ! from turbulence scheme
  DEALLOCATE(tcm_b          , STAT=izl); ist=ist+izl
  DEALLOCATE(tch_b          , STAT=izl); ist=ist+izl
  DEALLOCATE(tfv_b          , STAT=izl); ist=ist+izl

  IF (ist /= 0) THEN
     WRITE (*,*) ist, " fields in block structure could not be deallocated!! STOP!! "
     STOP 'DEALLOCATING block fields'
  ENDIF


!DL

END SUBROUTINE clean_up

!==============================================================================

END MODULE tsa_setup
