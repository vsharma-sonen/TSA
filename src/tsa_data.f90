!+ Source module "tsa_data" - Defining extra variables needed for Terra Stand Alone (TSA)
!------------------------------------------------------------------------------

MODULE tsa_data

!------------------------------------------------------------------------------
!
! Description:
!  The Module tsa_data defines extra variables needed for Terra Stand Alone (TSA)
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
! V3.0       2008/02/01 Felix Ament, MeteoSwiss
!  with contributions of
!  - Gerd Vogel and Jürgen Helmert (DWD): Brooks and
!    Coorey paramterization of drainage and diffusion,
!    revised formulation of infiltration, root density 
!    distribution
!  - Alexander Block (BTU): soil heat conductivity 
!    dependant of actual soil moisture content
!  - Gerd Vogel, Eric Jäger (ETHZ), Catherine Meissner (FZK) and others:
!    carefull testing and bug removal
! V4.13      2010/08/01 Guy de Morsier, MeteoSwiss (GDM)
!  Revised code to correspond to FieldExtra and COSMO version 4.13
! V5.01      2015-12-01 Yiftach Ziv, IMS (XYZ)
!  Replaced ireals by wp (working precision) (XYZ)
!  Arranged code to adhere to coding standards.
!
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!======================================================
!
! Declarations:
!
! Modules used:
!------------------------------------------------------------------------------

USE kind_parameters

!==============================================================================

IMPLICIT NONE

!==============================================================================

! 1. Parameters
!-------------------------------------------------------
 REAL(KIND=wp), PARAMETER :: &
     rundef=-9999.0_wp

 INTEGER, PARAMETER ::       &
     ke_soil_max=40        , & ! maximum number of soil layers 
     nastore=8             , & ! number additional 2d-fields (temporally averaged), 
                               ! which are stored just for output purposes 
     nstore=5                  ! see above, but not averaged


! 2. Namelist variables
!-------------------------------------------------------

 INTEGER                  :: &
     ntstep_max            , &  ! total number of time steps
     tincr_max             , &  ! max. gap in input data, expressed in hours
     nout_interval         , &  ! number of time steps between outputs  
     soiltyp_const         , &  ! soil type (IF lhomosoil=.true.)
     ntype_atminput        , &  ! data source of atmospheric screen level variables
                                ! =1: LM analysis at lowest model level
                                ! =2: LM analysis data at 2m (temperature, humidity) 
                                !     10m (wind)
                                ! =3: Synop analysis by DWD 
     ntype_raininput       , &  ! data source of precipitation data
                                ! =1: LM analysis 
                                ! =2: RADOLAN data by DWD
     ntype_radinput        , &  ! data source of radiative surface fluxes
                                ! =1: LM analysis including t_g
                                ! =2: LM analysis including t_s,t_snow,w_snow and
                                !        without apab_s 
     ntype_output          , &  ! extent of output:
                                ! =1: soil moisture profiles ONLY
                                ! =2: soil moisture and temperature profiles
                                ! =3: soil moisture and temperature profiles +
                                !     forcing data + some diagnostic fields
     lmgrib_aveint         , &  ! 
     nrecmax               , &  ! Maximum record number in binary forcing file
                                ! If the mor than nrecmax are read, reading starts
                                !  at the first record (looping).
     itype_hydparam        , &  ! type of soil moisture drainage and diffusion parameterization
                                ! =1: standard
                                ! =3: Brooks and Coorey + DWD soil types
                                ! =5: Brooks and Coorey + USDA soil types
!DL
     ngrid                 , &  ! number of ICON grid
     ngpi                  , &  ! generating process identifier for ICON domain
     ngp                   , &  ! number of gridpoints of ICON
!DL
     ke_model                   ! specify lowest model level (useful, IF input files
                                !  contain several levels and are unordered)

!DL
 CHARACTER (LEN=4)        :: &
     yform_read            , &  ! input  format: grb1(=libDWD), apix(=eccodes)
     yform_write                ! output format: grb1, api1, api2
                                ! grb1=libDWD,GRIB1/api1=eccodes,GRIB1/api2=eccodes,GRIB2
 CHARACTER (LEN=5)        :: &
     ymodel                     ! distinguish between COSMO and ICON
!DL


 CHARACTER (LEN=100)      :: &
     outdir                , &  ! output directory
     outprefix             , &  ! prefix common to all output files    
     metfiledir            , &  ! directory of meteorological forcing data
     metfileprefix         , &  ! prefix common to all meteorological forcing files
     radofiledir           , &  ! directory containing RADOLAN files
     constfilename         , &  ! gribfile with time constant surface parameters
                                !  
     soilinitdir           , &  ! directory containing initial soil conditions
     soilinitprefix             ! prefix common to all initial soil condition files

!    Initial conditions (soil parameters) in case of homoginous soil lhomoinit=.true. (lhomosoil=.true.)  
 REAL (KIND=wp) ::           &
     t_soil0(ke_soil_max)  , & ! soil temperature
     t_cl0                 , & ! climatological (constant) soil temperature 
                               !    lowest level 
     w_snow0               , & ! water content of snow storage 
     w_g0(ke_soil_max)     , & ! soil moisture
     w_i0                  , & ! water content of interception store
     w_cl0                 , & ! soil moisture of lowest layer 
     plcov_const           , & ! fixed plant cover (lvegadapt=.false.)
     plcov_min_const       , & ! minimum and maximum value of plant cover
     plcov_max_const       , & ! throughout a year 
     rootdp_const          , & ! root depth
     lai_const             , & ! fixed leaf area index (lvegadapt=.false.)
     lai_min_const         , & ! minimum and maximum value of leaf area index
     lai_max_const         , & !    throughout a year 
     z0_const              , & ! roughness length
     rstom_mn_const        , & ! minimum and maximum value of stomatal
     rstom_mx_const        , & !   resistance
     vegalb_const              ! albedo of vegetation

REAL (KIND=wp) ::                &
     czhls_const(ke_soil_max+1), & ! depth of soil layer boundaries
     lat                       , & ! latitude (in case of lhomosoil=.true.)
     hsurface                  , & ! surface elevation (in case of lhomosoil=.true.)
     dz                        , & ! reference level height of temperature and humidity
     dz_u                      , & ! reference level height of windspeed
     rain_fac                      ! factor to convert rain rates into LM units 
                                   !  (kg/s/m**2 = mm/s)
!DL
!    uuid of ICON grid for comparison
 CHARACTER(LEN=1)   ::  uuid(16)
!DL

 CHARACTER (LEN=14) ::  &
     ydate_ini        , & ! initial date of simulation (yyyymmddhh) 
     ydate_end            ! final date of simulation (yyyymmddhh)


 LOGICAL ::          &
     lrel_in       , & ! homogenoeus input of soil moisture relative to
                       ! field capacity (lhomoinit=.true.)
     lvol_in       , & ! homogenoeus input of soil moisture in vol. %
                       ! (lhomoinit=.true.)
     lvegadapt     , & ! simulating the anual cycle of vegetation ?   
     lhomosoil     , & ! homogenoeus soil parameters ?
     lhomoinit     , & ! homogenous initial soil conditions ?   
     lconstout     , & ! including constant fields in the first outputfile?
     lrootadapt    , & ! varying the root depth according to anaual cycle?
     lcheck        , & ! some checks
     lmulti_in     , & ! specifying intial soil conditions.true: multi layer data 
                       ! is provided; false: two layer data is available
     lbug_u10m     , & ! emulates 10m wind bug
     lz0local      , & ! GRIB only: reads local roughness length (ee=101, tab=202) 
                       ! without contribution due to subgrid-scale orographie
     lconstvegalb  , & ! spatially varying plant albedo (ee=213, tab=202)
     lext_monthly  , & ! experimental! external parameters are read every new month
     linfil_revised, & ! revised parameterization of infiltration
     lhourly_data  , & ! GRIB input contains hourly data (precip, radiation)
     ldestaggeruv  , & ! destagger velocities on input from external grib file
     lpar          , & ! try to read PAR from forcing data?
     lgettcl       , & ! read T_CL from external parameter file
     lcalc             ! debug option: lcalc=false ==> terra and parturs are
                       !  not invoked

  INTEGER ::         &
    nsub_x         , & ! number of sub-regions in zonal direction
    nsub_y         , & ! number of sub-regions in meridional direction
    which_subreg       ! current selected sub-region

! 3. Standard LM variables
!-------------------------------------------------------
 INTEGER                  :: ierror   ! error status variable
!>XYZ: LEN changed from 80 to 255 to correspond to src_soil_multilay.f90
 CHARACTER (LEN=255)      :: yerror
!<XYZ
 CHARACTER (LEN=80)       :: dummy



! 4. Additional fields
!-------------------------------------------------------

  !!DL
  ! "mind" fields formerly defined in src_block_fields
  ! Fields for index-reorganization
  INTEGER,  ALLOCATABLE :: &
    mind_ilon(:,:)       , & ! keeps the i (longitude) index
    mind_jlat(:,:)           ! keeps the j (latitude)  index
                             ! for every grid point in the (nproma,nblock) ordering
  !!DL

  ! fields, which are not declared in data_fields.f90
  REAL (KIND=wp), ALLOCATABLE :: &
    ! needed, as there is no tracer structure in TSA
    qv      (:,:,:,:)    , & ! specific water vapor content                  (kg/kg)
    qv_bd   (:,:,:,:)    , & ! boundary field for qv                       (kg/kg)

    !US but rlat, rlon are the geographical lat/lon in COSMO? Is this different here?
    rlat_geo(:,:)        , & ! geographical latitude                         ( rad )
    rlon_geo(:,:)        , & ! geographical longitude                        ( rad )

    plcov_mn(:,:)        , & ! min of plant cover
    plcov_mx(:,:)        , & ! max of plant cover
    rootdp0 (:,:)        , & ! maximum rootdepth
    lai_mn  (:,:)        , & ! min of leaf area index
    lai_mx  (:,:)        , & ! max of leaf area index
    z0      (:,:)        , & ! roughness length
    rstom_mx(:,:)        , &
    rstom_mn(:,:)        , &
    vegalb  (:,:)

  ! fields, which are not declared in data_block_fields.f90
  REAL (KIND=wp), ALLOCATABLE ::   &
    qv_b    (:,:)       , & ! needed, as there is no tracer structure in TSA
    z0_b      (:)           ! surface roughness                             (  m )

! Additional meteorological boundary fields
 REAL (KIND=wp), ALLOCATABLE :: &
     ps_bd(:,:,:)             , & ! boundary field of surface pressure
     prr_gsp_bd(:,:,:)        , & ! boundary field of rain
     prs_gsp_bd(:,:,:)        , & ! --- """ --- of snow
     so_down_bd(:,:,:)        , & ! --- """ --- of downwelling solar radiation
     th_down_bd(:,:,:)        , & ! --- """ --- of downwelling thermal radiation
     pabs_bd(:,:,:)               ! --- """ --- of downwelling PAR

 REAL (KIND=wp), ALLOCATABLE :: &
     astore(:,:,:)            , & !
     store(:,:,:)             , & !
     zflmg(:,:,:)             , & !
     runoff(:,:,:)            , & !
     sobs2(:,:)               , & !
     thbs2(:,:)                   !

END MODULE tsa_data
