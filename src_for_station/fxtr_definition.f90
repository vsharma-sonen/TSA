!+ Program fieldextra, specific support: global definitions
!+ from ~bettems/projects/fieldextra_v10_1_0
!
!****************************************************************************
MODULE fxtr_definition
!=============================================================================
!
! Global definitions
! > User defined kinds
! > User defined types & operator overloading
! > Size of INT
! > Size of global repository and of namelist cache
! > String length
! > Dimensions of static arrays
! > Dimensions used in DWD GRIB1 library
! > Physical, mathematical and other constants
! > Symbolic constants
! > Fortran units
! > Thresholds and scalings for formating output
!
!
! Current Code Owner
! ------------------
! Jean-Marie Bettems
! Swiss Meteorological Institute, Zurich, Switzerland
! jean-marie.bettems@meteoswiss.ch
!
!-----------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE



! Public entities
!================
! User defined kinds
PUBLIC :: kind_idate, kind_ilarge, kind_opkey
! User defined types
! ... misc
PUBLIC :: ty_r_fmt
! ... external resources
PUBLIC :: ty_fld_dict, ty_loc, ty_loc_cat, ty_reg, ty_kal, ty_adaboost, ty_ab_features
! ... mesh characteristics
PUBLIC :: ty_scan_mode, ty_hc_desc, ty_hc_cache, ty_vc_desc
! ... field characteristics
PUBLIC :: ty_fld_product, ty_fld_orig, ty_tr_unit, ty_fld_epsid, ty_fld_prob
! ... operator characteristics
PUBLIC :: ty_logical_cond, ty_compare_op, ty_conv_op, ty_regrid_op
! ... namelist specifications
PUBLIC :: ty_out_spec, ty_out_mode, ty_fld_spec, ty_fld_spec_root, ty_out_dred
! ... data storage
PUBLIC :: ty_1d_field, ty_2d_field, ty_ic_data, ty_gb_field, ty_calc_spec, ty_calc_outdata, ty_out_store
! Operator overloading
PUBLIC :: OPERATOR(==), OPERATOR(/=)
! Size of integer
PUBLIC :: int_in_byte
! Size of global repository and of namelist cache
PUBLIC :: grep_size, nl_cache_size
! String length
PUBLIC :: ln_exp, ln_str, ln_loc, ln_reg, ln_buf
! Dimensions of static arrays
PUBLIC :: mx_ifile, mx_ij, mx_time, mx_tlag
PUBLIC :: mx_ofile, mx_field, mx_lev, mx_atm_lev, mx_svat_lev, mx_member, mx_probability, mx_gp
PUBLIC :: mx_inout
PUBLIC :: mx_location, mx_loc_cat, mx_loc_addfile
PUBLIC :: mx_set_reg
PUBLIC :: mx_vc_coef, mx_grids
PUBLIC :: mx_fld_dict, mx_fld_key, mx_ada_set
PUBLIC :: mx_model_spec, mx_precip_component
PUBLIC :: mx_size_values, mx_iteration, mx_alias, mx_token, mx_tags
PUBLIC :: mx_op_coef, mx_conv_op
! Dimensions used in DWD GRIB1 library
PUBLIC :: gb_nrbit, gb_ldims, gb_ldsup, gb_lpds, gb_lgds, gb_lbms, gb_lbds, gb_lbmp
! Physical, mathematical and other constants
PUBLIC :: pc_g, pc_emissivity_surface, pc_boltzman_cst, pc_r_d, pc_r_v, pc_rdv, pc_o_rdv
PUBLIC :: pc_t0, pc_b1, pc_b2w, pc_b2i, pc_b3, pc_b4w, pc_b4i, pc_lh_v, pc_lh_s, pc_cp_d 
PUBLIC :: pc_rvd_o, pc_rdocp, pc_std_lapse, pc_std_t0, pc_std_p0, pc_omega, pc_earth_rad
PUBLIC :: pc_cel_to_fahr1, pc_cel_to_fahr2
PUBLIC :: mc_rad_to_deg, mc_deg_to_rad, mc_eps, mc_pi
PUBLIC :: oc_tab
! Symbolic constants
PUBLIC :: sc_sv_n, sc_sv_w, sc_sv_s, sc_sv_e, sc_sv_all
PUBLIC :: sc_hor_n, sc_hor_w, sc_hor_s, sc_hor_e
PUBLIC :: sc_subset_none, sc_subset_gplist, sc_subset_loclist, sc_subset_domain, sc_subset_loccat
PUBLIC :: sc_hc_type_unknown, sc_hc_type_latlon, sc_hc_type_rotlatlon
PUBLIC :: sc_hc_type_swiss, sc_hc_type_gboaga_w, sc_hc_type_gboaga_e
PUBLIC :: sc_earth_spherical, sc_earth_spheroidal, sc_earth_unknown
PUBLIC :: sc_vc_notset, sc_vc_unsafe, sc_vc_safe
PUBLIC :: sc_vc_type_unknown, sc_vc_type_pressure, sc_vc_type_height, sc_vc_type_sleve
PUBLIC :: sc_vc_coding_unknown, sc_vc_coding_cosmo_old1, sc_vc_coding_cosmo_old2, sc_vc_coding_cosmo_new
PUBLIC :: sc_model_unknown, sc_model_cosmo, sc_model_ifs, sc_model_na
PUBLIC :: sc_pcat_unknown, sc_pcat_determinist, sc_pcat_epsmemb, sc_pcat_epsprob, sc_pcat_nbhprob
PUBLIC :: sc_pcat_epsmean, sc_pcat_epsmedian, sc_pcat_observation
PUBLIC :: sc_suite_unknown, sc_suite_standard, sc_suite_dwd_eps, sc_suite_cosmo_leps
PUBLIC :: sc_ltyp_unknown, sc_ltyp_single, sc_ltyp_multilev, sc_ltyp_multilay
PUBLIC :: sc_vref_unknown, sc_vref_geog, sc_vref_native
PUBLIC :: rundef, iundef, cundef
! Fortran units
PUBLIC :: unit_ascii, unit_template, unit_cache, unit_location, unit_kalman, unit_region
PUBLIC :: unit_diagnostic, unit_adaboost, unit_dictionary, unit_nl
! Thresholds and scalings for formating output
PUBLIC :: out_scale_hvc, out_precip_min, out_cloud_min, out_rsnow_min



! Parameters and kinds
!=====================
! Number of bytes used by integer storage
INTEGER, PARAMETER      :: int_in_byte = BIT_SIZE(1)/8

! Size of namelist cache
INTEGER, PARAMETER      :: nl_cache_size = 5000

! KIND parameter for storage of large integers
INTEGER, PARAMETER      :: kind_ilarge = SELECTED_INT_KIND(12)
! KIND Parameters for storage of date (yyyymmddhhmm)
INTEGER, PARAMETER      :: kind_idate = SELECTED_INT_KIND(12)

! Parameters and kind for global repository
! > global repository is managed in fxtr_storage
! > size of global repository grep_size is expressed as a power of two
! > for each operator a unique integer identifier is derived from the index
!   of the different elements defining the operator; the kind of the identifier
!   is chosen to allow values as large as grep_size**3 (i.e. kind is
!   ceiling[LOG10(2)*3*grep_size_power])
INTEGER, PARAMETER      :: grep_size_power = 12
INTEGER, PARAMETER      :: grep_size = 2**grep_size_power
INTEGER, PARAMETER      :: kind_opkey = SELECTED_INT_KIND(grep_size_power)


! String size (name convention ln_*)
!-----------------------------------
INTEGER, PARAMETER      :: ln_exp      = 3          ! Length of experiment tag strings
INTEGER, PARAMETER      :: ln_loc      = 8          ! Length of location%tag string
INTEGER, PARAMETER      :: ln_reg      = 20         ! Length of location%region string
INTEGER, PARAMETER      :: ln_str      = 384        ! Length of standard string
INTEGER, PARAMETER      :: ln_buf      = 2000000    ! Length of buffer string


! Dimensions of static arrays (name convention: mx_*)
!----------------------------------------------------
! 1. Sizing of input data 
INTEGER, PARAMETER      :: mx_ifile       = 5000      ! Max. # of input files
INTEGER, PARAMETER      :: mx_ij          = 1200*1200 ! Max. size of horizontal mesh
INTEGER, PARAMETER      :: mx_time        = 750       ! Max. # of validation dates
INTEGER, PARAMETER      :: mx_tlag        = 25        ! Max. # of times in tlag list
! 2. Sizing of output data
INTEGER, PARAMETER      :: mx_ofile       = 2500      ! Max. # of output files
INTEGER, PARAMETER      :: mx_field       = 1200      ! Max. # of collected fields for each output
INTEGER, PARAMETER      :: mx_atm_lev     = 100       ! Max. # of collected levels for each atm. field
INTEGER, PARAMETER      :: mx_svat_lev    = 10        ! Max. # of collected levels for each SVAT field
INTEGER, PARAMETER      :: mx_lev         = MAX(mx_atm_lev, mx_svat_lev)
INTEGER, PARAMETER      :: mx_probability = 20        ! Max. # of collected probability intervals 
INTEGER, PARAMETER      :: mx_member      = 50        ! Max. # of collected EPS members
INTEGER, PARAMETER      :: mx_gp          = 2500      ! Max. # locations via loclist, i/jlist, loccat
! 3. Max. number of (compressed) input/output files pairs
INTEGER, PARAMETER      :: mx_inout       = 10000     
! 4. Sizing of data%values array
!    Used to put a soft limit on memory allocation (normally most of the memory
!    is allocated for data%values). On a machine with 32 bits REAL, a value of
!    500000000 corresponds to a limit of 2 GB.
INTEGER(KIND=kind_ilarge), PARAMETER                &
                        :: mx_size_values = 750000000

! Max. # auxiliary file
! * Field dictionary
INTEGER, PARAMETER      :: mx_fld_dict  = 10    
! * AdaBoost coefficients sets
INTEGER, PARAMETER      :: mx_ada_set  = 50    

! Max. # vertical coordinates coefficients
INTEGER, PARAMETER      :: mx_vc_coef = 1000    

! Max. # keys in each field dictionary
INTEGER, PARAMETER      :: mx_fld_key = 2000    

! Max. # different horizontal meshes
INTEGER, PARAMETER      :: mx_grids = 10

! Max. # field name aliases (ignore tri)
INTEGER, PARAMETER      :: mx_alias = 10

! Max. # of processing iterations
! (the input routines in fxtr_control must be adapted when increasing this value)
INTEGER, PARAMETER      :: mx_iteration = 4

! Max. # tokens in a logical expression
INTEGER, PARAMETER      :: mx_token = 10

! Max. # tags in list of tags to use
INTEGER, PARAMETER      :: mx_tags = 10

! Max. # coefficients defining comparison operator (see ty_compare_op)
INTEGER, PARAMETER      :: mx_op_coef = 10

! Max. # convolution operators kept in cache
INTEGER, PARAMETER      :: mx_conv_op = 100

! Max. # of model specifications in namelist file
INTEGER, PARAMETER      :: mx_model_spec = 20
! Max. # of precipitation components
INTEGER, PARAMETER      :: mx_precip_component = 10

! Related to locations definitions
! Maximum number of locations in location_list file
INTEGER, PARAMETER      :: mx_location = 10000
! Largest allowed category in location_list file (0 <= cat <= mx_loc_cat)
INTEGER, PARAMETER      :: mx_loc_cat = 20
! Maximum number of files defining additional location characteristics
INTEGER, PARAMETER      :: mx_loc_addfile = 20

! Related to regions definitions
! Maximum number of region set in region_list file
INTEGER, PARAMETER      :: mx_set_reg = 10


! Dimensions used in DWD GRIB1 library (name convention: gb_*)
!-------------------------------------------------------------
INTEGER, PARAMETER     :: gb_nrbit =     16    ! Packing rate, in number of bits per value
INTEGER, PARAMETER     :: gb_ldims =     20    ! Size of dims vector
INTEGER, PARAMETER     :: gb_ldsup =     73    ! Size of dsup vector
INTEGER, PARAMETER     :: gb_lpds  =    321    ! Size of product definition block
INTEGER, PARAMETER     :: gb_lgds  =    626    ! Size of grid description block
INTEGER, PARAMETER     :: gb_lbms  =      3    ! Size of bitmap data description block
INTEGER, PARAMETER     :: gb_lbds  =     11    ! Size of binary data description block
INTEGER, PARAMETER     :: gb_lbmp  = 100000    ! Size of bitmap


! Physical, mathematical and other constants
!-------------------------------------------
! Physical constants (name convention: pc_*)
! acceleration due to gravity
REAL, PARAMETER        :: pc_g     = 9.80665
! earth rotation velocity
REAL, PARAMETER        :: pc_omega = 7.292E-5 ! [rad/s]
! earth radius
REAL, PARAMETER        :: pc_earth_rad = 6371.229E3
! for radiation processes
REAL, PARAMETER        :: pc_emissivity_surface = 0.996
REAL, PARAMETER        :: pc_boltzman_cst       = 5.6697E-8 
! for thermodynamic functions (values taken from LM v3.2)
REAL, PARAMETER        :: pc_t0    = 273.15     
REAL, PARAMETER        :: pc_r_d   = 287.05    ! Gas constant for dry air
REAL, PARAMETER        :: pc_r_v   = 461.51    ! Gas constant for water vapor
REAL, PARAMETER        :: pc_rdv   = pc_r_d / pc_r_v
REAL, PARAMETER        :: pc_o_rdv = 1. - pc_rdv
REAL, PARAMETER        :: pc_b1    = 610.78
REAL, PARAMETER        :: pc_b2w   = 17.2693882
REAL, PARAMETER        :: pc_b2i   = 21.8745584
REAL, PARAMETER        :: pc_b3    = 273.16
REAL, PARAMETER        :: pc_b4w   = 35.86
REAL, PARAMETER        :: pc_b4i   = 7.66
REAL, PARAMETER        :: pc_lh_v  = 2.501E6  ! latent heat of vapourization
REAL, PARAMETER        :: pc_lh_s  = 2.835E6  ! latent heat of sublimation
REAL, PARAMETER        :: pc_cp_d  = 1005.0   ! specific heat of dry air at constant pressure
REAL, PARAMETER        :: pc_rvd_o = (pc_r_v / pc_r_d) - 1.0
REAL, PARAMETER        :: pc_rdocp = pc_r_d/pc_cp_d 
! for transforming Celsius in Fahrenheit
REAL, PARAMETER        :: pc_cel_to_fahr1  = 1.8
REAL, PARAMETER        :: pc_cel_to_fahr2  = 32.0
! for standard atmosphere
REAL, PARAMETER        :: pc_std_lapse = 0.0065   ! reference temperature lapse rate [K/m]
REAL, PARAMETER        :: pc_std_t0    = 288.15   ! reference temperature at sea level [K]
REAL, PARAMETER        :: pc_std_p0    = 101325.  ! reference pressure at sea level [Pa]

! Mathematical constants (name convention: mc_*)
REAL, PARAMETER        :: mc_pi            = 3.141592653
REAL, PARAMETER        :: mc_rad_to_deg    = 57.2957795
REAL, PARAMETER        :: mc_deg_to_rad    = 0.0174532925
REAL, PARAMETER        :: mc_eps           = 1.E-10

! Other constants (name convention: oc_*)
CHARACTER(LEN=1), PARAMETER :: oc_tab = ACHAR(9)     ! ASCII code for tabulator


! Symbolic constants (name convention: sc_*)
!-------------------------------------------
!   for model type
INTEGER, PARAMETER     :: sc_model_unknown      = -1   ! unknown
INTEGER, PARAMETER     :: sc_model_cosmo        = 0    ! COSMO model
INTEGER, PARAMETER     :: sc_model_ifs          = 1    ! IFS model
INTEGER, PARAMETER     :: sc_model_na           = 100  ! not available (observations)
!   for product category
INTEGER, PARAMETER     :: sc_pcat_unknown       = -1   ! unknown
INTEGER, PARAMETER     :: sc_pcat_determinist   = 0    ! model
INTEGER, PARAMETER     :: sc_pcat_epsmemb       = 1    !
INTEGER, PARAMETER     :: sc_pcat_epsprob       = 2    !
INTEGER, PARAMETER     :: sc_pcat_epsmean       = 3    !
INTEGER, PARAMETER     :: sc_pcat_epsmedian     = 4    !
INTEGER, PARAMETER     :: sc_pcat_nbhprob       = 5    !
INTEGER, PARAMETER     :: sc_pcat_observation   = 100  ! observation
!   for production suite
INTEGER, PARAMETER     :: sc_suite_unknown      = -1   ! unknown
INTEGER, PARAMETER     :: sc_suite_standard     = 0    ! standard
INTEGER, PARAMETER     :: sc_suite_dwd_eps      = 1    ! DWD EPS 
INTEGER, PARAMETER     :: sc_suite_cosmo_leps   = 2    ! COSMO LEPS
!   for grid point subset
INTEGER, PARAMETER     :: sc_subset_none    = 0
INTEGER, PARAMETER     :: sc_subset_gplist  = 1
INTEGER, PARAMETER     :: sc_subset_loclist = 2
INTEGER, PARAMETER     :: sc_subset_loccat  = 3
INTEGER, PARAMETER     :: sc_subset_domain  = 4
!   for type of hcoord 
INTEGER, PARAMETER     :: sc_hc_type_unknown   = -1
INTEGER, PARAMETER     :: sc_hc_type_latlon    = 0
INTEGER, PARAMETER     :: sc_hc_type_rotlatlon = 10
INTEGER, PARAMETER     :: sc_hc_type_swiss     = 200
INTEGER, PARAMETER     :: sc_hc_type_gboaga_w  = 201
INTEGER, PARAMETER     :: sc_hc_type_gboaga_e  = 202
!   for reference system for vector fields (see GRIB code table 7)
INTEGER, PARAMETER     :: sc_vref_unknown       = -1   ! unknown
INTEGER, PARAMETER     :: sc_vref_geog          = 0    ! easterly and northerly directions
INTEGER, PARAMETER     :: sc_vref_native        = 1    ! native (increasing i and j co-ordinates)
!   for earth model (GRIB code table 7)
INTEGER, PARAMETER     :: sc_earth_unknown = -1        ! unknown
INTEGER, PARAMETER     :: sc_earth_spherical = 0       ! sphere of radius 6367.47km
INTEGER, PARAMETER     :: sc_earth_spheroidal = 1      ! oblate spheroidal
!   for origin of vcoord information
INTEGER, PARAMETER     :: sc_vc_notset = 0
INTEGER, PARAMETER     :: sc_vc_unsafe = 1
INTEGER, PARAMETER     :: sc_vc_safe   = 2
!   for coding of vcoord 
INTEGER, PARAMETER     :: sc_vc_coding_unknown    = -1
INTEGER, PARAMETER     :: sc_vc_coding_cosmo_old1 = 1
INTEGER, PARAMETER     :: sc_vc_coding_cosmo_old2 = 2
INTEGER, PARAMETER     :: sc_vc_coding_cosmo_new  = 3
!   for type of vcoord (values must not be changed, because used in called LM routine)
INTEGER, PARAMETER     :: sc_vc_type_unknown  = 0
INTEGER, PARAMETER     :: sc_vc_type_pressure = 1
INTEGER, PARAMETER     :: sc_vc_type_height   = 2
INTEGER, PARAMETER     :: sc_vc_type_sleve    = 3
!   for type of levels
INTEGER, PARAMETER     :: sc_ltyp_unknown  = -1
INTEGER, PARAMETER     :: sc_ltyp_single   = 0
INTEGER, PARAMETER     :: sc_ltyp_multilev = 1
INTEGER, PARAMETER     :: sc_ltyp_multilay = 2
!   for horizon, in 4 directions
INTEGER, PARAMETER     :: sc_hor_s  = 1
INTEGER, PARAMETER     :: sc_hor_w  = 2
INTEGER, PARAMETER     :: sc_hor_n  = 3
INTEGER, PARAMETER     :: sc_hor_e  = 4
!   for skyview directions, total and in 4 directions
INTEGER, PARAMETER     :: sc_sv_s   = 1
INTEGER, PARAMETER     :: sc_sv_w   = 2
INTEGER, PARAMETER     :: sc_sv_n   = 3
INTEGER, PARAMETER     :: sc_sv_e   = 4
INTEGER, PARAMETER     :: sc_sv_all = 5

! "Undefined" flags
REAL, PARAMETER        :: rundef = -99999.
INTEGER, PARAMETER     :: iundef = -99999
CHARACTER(LEN=*), PARAMETER :: cundef = "???"   


! Fortran units (name convention: unit_*)
!----------------------------------------
INTEGER, PARAMETER     :: unit_ascii      = 10  ! For ascii output
INTEGER, PARAMETER     :: unit_template   = 11  ! For ascii template file
INTEGER, PARAMETER     :: unit_cache      = 12  ! For data cache file
INTEGER, PARAMETER     :: unit_location   = 13  ! For locations list
INTEGER, PARAMETER     :: unit_kalman     = 14  ! For kalman coefficients
INTEGER, PARAMETER     :: unit_region     = 15  ! For definition of regions
INTEGER, PARAMETER     :: unit_adaboost   = 16  ! For adaboost coefficients
INTEGER, PARAMETER     :: unit_dictionary = 17  ! For field dictionary
INTEGER, PARAMETER     :: unit_diagnostic = 18  ! For diagnostic
INTEGER, PARAMETER     :: unit_nl         = 19  ! For namelist


! Thresholds and scalings for formating output (name convention: out_*)
!----------------------------------------------------------------------
! Scaling factor for height based v.coord.
REAL, PARAMETER        :: out_scale_hvc = 0.0001 
! Minimum precip. value to calculate precipitation percentage
REAL, PARAMETER        :: out_precip_min = 0.05
! Minimum specific cloud (qi+qc) to compute cloud ice ratio
REAL, PARAMETER        :: out_cloud_min = 1.E-10
! Minimum snow density [kg/m**3] to compute h_snow
REAL, PARAMETER        :: out_rsnow_min = 1.0



! Types definitions  (name convention: ty_*)
!===========================================

! Misc
!-----

! To store information about Fortran format descriptor (for real numbers)
TYPE         :: ty_r_fmt
  INTEGER                                      :: rep           ! repeat count
  CHARACTER(LEN=2)                             :: desc          ! edit descriptor
  INTEGER                                      :: w             ! field width    
  INTEGER                                      :: d             ! number of digits
  CHARACTER(LEN=10)                            :: stripped      ! format w/o repeat count
  INTEGER                                      :: ierr          ! error flag
END TYPE ty_r_fmt


! External resources
!-------------------

! To store field dictionary
TYPE         :: ty_fld_dict
  CHARACTER(LEN=ln_str)                        :: name
  INTEGER                                      :: nbr_key
  CHARACTER(LEN=ln_str), DIMENSION(mx_fld_key) :: key
  LOGICAL, DIMENSION(mx_fld_key)               :: is_constant
  INTEGER, DIMENSION(mx_fld_key)               :: order_vec_component
  CHARACTER(LEN=ln_str), DIMENSION(mx_fld_key) :: ass_vec_component
  LOGICAL, DIMENSION(mx_fld_key)               :: is_slevel
  LOGICAL, DIMENSION(mx_fld_key)               :: is_svat_field
  CHARACTER(LEN=ln_str), DIMENSION(mx_fld_key) :: ass_mlev_field
  CHARACTER(LEN=ln_str), DIMENSION(mx_fld_key) :: parent_field
  LOGICAL, DIMENSION(mx_fld_key)               :: mtlev_parent
  INTEGER, DIMENSION(mx_fld_key)               :: grib1_table_version
  INTEGER, DIMENSION(mx_fld_key)               :: grib1_parameter
  INTEGER, DIMENSION(mx_fld_key)               :: grib1_level_type
  INTEGER, DIMENSION(mx_fld_key)               :: grib1_level
  INTEGER, DIMENSION(mx_fld_key)               :: grib1_trange_type
  INTEGER, DIMENSION(mx_fld_key)               :: grib2_discipline
  INTEGER, DIMENSION(mx_fld_key)               :: grib2_parameter_cat
  INTEGER, DIMENSION(mx_fld_key)               :: grib2_parameter_numb
  INTEGER, DIMENSION(mx_fld_key)               :: grib2_stat_proc
  INTEGER, DIMENSION(mx_fld_key)               :: grib2_type_ffs
  INTEGER, DIMENSION(mx_fld_key)               :: grib2_factor_ffs
  INTEGER, DIMENSION(mx_fld_key)               :: grib2_value_ffs
  INTEGER, DIMENSION(mx_fld_key)               :: grib2_type_sfs
  INTEGER, DIMENSION(mx_fld_key)               :: grib2_factor_sfs
  INTEGER, DIMENSION(mx_fld_key)               :: grib2_value_sfs
  INTEGER, DIMENSION(mx_fld_key)               :: grib2_type_gen_proc
  CHARACTER(LEN=ln_str), DIMENSION(mx_fld_key) :: netcdf_units
  CHARACTER(LEN=ln_str), DIMENSION(mx_fld_key) :: netcdf_name
END TYPE ty_fld_dict


! To store definition of locations 
TYPE          :: ty_loc
  CHARACTER(LEN=30)                            :: description ! Description of the location 
  CHARACTER(LEN=ln_loc)                        :: tag         ! Unique string defining location
  CHARACTER(LEN=ln_reg)                        :: region      ! Description of region
  INTEGER                                      :: flag        ! Flag (dd,ff,NA,NA,NA,NA,NA,NA)
  INTEGER                                      :: index       ! Station id 
  INTEGER                                      :: category    ! Station category (user defined)
  INTEGER                                      :: altitude    ! Station altitude
  REAL                                         :: latitude    ! Station geog. latitude
  REAL                                         :: longitude   ! Station geog. longitude
  ! derived quantities
  INTEGER                                      :: ic_idx      ! Index of incore field to derive gpi/gpj
  INTEGER                                      :: gpi         ! (i,j)-coordinates and height of 
  INTEGER                                      :: gpj         ! associated grid point, or 0 if 
  INTEGER                                      :: gph         ! calculated by fieldextra
  INTEGER                                      :: regset_idx  ! Index of associated region set
  INTEGER                                      :: reg_idx     ! Index of associated region
  ! additional characteristics
  REAL, DIMENSION(4)                           :: horizon     ! horizon factors (S,W,N,E)
  REAL, DIMENSION(5)                           :: skyview     ! skyview factors (S,W,N,E,ALL)
END TYPE ty_loc


! To store information about location categories
TYPE          :: ty_loc_cat
  INTEGER                                       :: nbr_locations  ! Number of locations in category
  CHARACTER(LEN=ln_loc), DIMENSION(mx_location) :: tag_list       ! List of locations in category
END TYPE ty_loc_cat


! To store definitions of regions
TYPE          :: ty_reg
  INTEGER                                      :: nbr_summits  ! Nbr. of summits defining the region
  REAL, DIMENSION(:), POINTER                  :: sum_lat      ! Summit latitude
  REAL, DIMENSION(:), POINTER                  :: sum_lon      ! Summit longitude
  REAL, DIMENSION(:), POINTER                  :: sum_gpi_r    ! Summit i-coord. (real number!)
  REAL, DIMENSION(:), POINTER                  :: sum_gpj_r    ! Summit j-coord. (real number!)
  INTEGER                                      :: nbr_gp       ! Nbr. of grid points inside of region
  INTEGER, DIMENSION(4)                        :: bbox         ! Region bounding box (imin,jmin,imax,jmax)
  INTEGER, DIMENSION(:), POINTER               :: gpi          ! i-coord. of points inside region
  INTEGER, DIMENSION(:), POINTER               :: gpj          ! j-coord. of points inside region
  INTEGER                                      :: ic_idx       ! Index of incore field to derive gpi/gpj
  CHARACTER(LEN=ln_loc)                        :: tag          ! Region tag
END TYPE ty_reg


! To store Kalman coefficients
TYPE          :: ty_kal
  INTEGER                                      :: year           !  
  INTEGER                                      :: month          ! --> validation date 
  INTEGER                                      :: day            !____
  REAL, DIMENSION(mx_gp,mx_time)               :: coef1_t2m      ! 
  REAL, DIMENSION(mx_gp,mx_time)               :: coef2_t2m      ! --> kal. coef. for T_2M
  LOGICAL, DIMENSION(mx_gp,mx_time)            :: coef_t2m_set   !____ 
  REAL, DIMENSION(mx_gp,mx_time)               :: coef1_td2m     ! 
  REAL, DIMENSION(mx_gp,mx_time)               :: coef2_td2m     ! --> kal. coef. for TD_2M
  LOGICAL, DIMENSION(mx_gp,mx_time)            :: coef_td2m_set  !____ 
END TYPE ty_kal


! To store adaboost coefficents
TYPE          :: ty_adaboost
  CHARACTER(LEN=ln_str), POINTER, DIMENSION(:) :: outfieldnames      ! name of the out field
  INTEGER                                      :: nbaseclassifiers   ! number of base classifiers
  CHARACTER(LEN=ln_str)                        :: baseclassifiername ! name of base classifier
  INTEGER                                      :: nfeatures          ! number of features
  CHARACTER(LEN=ln_str), POINTER, DIMENSION(:) :: featurenames       ! name of the features
  INTEGER, POINTER, DIMENSION(:)               :: levels             ! model level of the features
  ! param. for base classifiers: 
  INTEGER, POINTER, DIMENSION(:)               :: features           !   concerning feature
  REAL, POINTER, DIMENSION(:)                  :: thresholds         !   threshold value
  REAL, POINTER, DIMENSION(:)                  :: signs              !   side of the event
  REAL, POINTER, DIMENSION(:)                  :: weights            !   voting weight
END TYPE ty_adaboost


! Used in adaboost to store mapping between feature names, indices and field data
! Meaning of type:  1: 1 dim., constant over time
!                   2: 1 dim., constant over grid point
!                   3: 2 dim.
!                   4: incore field
TYPE          :: ty_ab_features
  CHARACTER(LEN=ln_str)                        :: name        ! name of the feature/field
  INTEGER                                      :: level       ! field level
  INTEGER                                      :: type        ! type/dimensionality of the data
  REAL, POINTER, DIMENSION(:,:)                :: data        ! pointer to data
  REAL, POINTER, DIMENSION(:)                  :: ic_data     ! pointer to incore data
END TYPE ty_ab_features


! Mesh characteristics
!---------------------

! To store information about horizontal coordinates
! ... data scanning mode
TYPE         :: ty_scan_mode
  LOGICAL                          :: valid             ! valid definition of scanning mode ?
  LOGICAL                          :: defered           ! definition of scanning mode defered ?
  LOGICAL                          :: positive_i        ! scanning mode in +i direction ?
  LOGICAL                          :: positive_j        ! scanning mode in +j direction ?
  LOGICAL                          :: consecutive_i     ! consecutive i ?
END TYPE ty_scan_mode

! ... grid description
!     (1) coordinate system (i,j) and units:
!         type = sc_hc_type_latlon (geog. lat/lon)
!           i-axis: West-East direction, j-axis: South-North direction, units: degree
!         type = sc_hc_type_rotlatlon (rotated lat/lon)
!           i-axis: rot. West-East direction, j-axis: rot. South-North direction, units: degree
!         type = sc_hc_type_swiss (metric swiss)
!           i-axis: West-East direction, j-axis: South-North direction, units: meters
!         type = sc_hc_type_gboaga_w, sc_hc_type_gboaga_e (metric italian)
!           i-axis: West-East direction, j-axis: South-North direction, units: meters
!     (2) values are set to undef for information not relevant for the considered type
TYPE         :: ty_hc_desc
  INTEGER                          :: type              ! Data representation type (sc_hc_type_*)
  TYPE(ty_scan_mode)               :: scan_mode         ! Scanning mode
  INTEGER                          :: earth_model       ! Earth model (sc_earth_*)
  REAL                             :: rotation_pollat   ! Grid rotation - North pole latitude (deg)
  REAL                             :: rotation_pollon   !               - North pole longitude (deg)
  REAL                             :: rotation_gamma    !               - Rotation angle (deg)
  INTEGER                          :: ni                ! Nbr. of points in i-direction 
  INTEGER                          :: nj                ! Nbr. of points in j-direction 
  INTEGER                          :: n_points          ! Total nbr. of points (i.e. ni*nj)
  REAL                             :: di                ! Increment in i direction 
  REAL                             :: dj                ! Increment in j direction 
  REAL                             :: swi               ! Start corner, i-coord. 
  REAL                             :: swj               ! Start corner, j-coord. 
  REAL                             :: nei               ! End corner, i-coord. 
  REAL                             :: nej               ! End corner, j-coord. 
  INTEGER                          :: vref              ! Reference system for vector fields (sc_vref_*)
END TYPE ty_hc_desc

! ... for caching some information
TYPE         :: ty_hc_cache
  TYPE(ty_scan_mode)               :: scan_mode         
  INTEGER                          :: ni                
  INTEGER                          :: nj                
  INTEGER, DIMENSION(:,:), POINTER :: look_up_table     ! (i,j) --> position in data array
  INTEGER, DIMENSION(:), POINTER   :: ij2data           ! position in (i,j) --> position in data
END TYPE ty_hc_cache


! To store information about vertical coordinates
TYPE         :: ty_vc_desc
  INTEGER                          :: coding            ! Coding of v. coord. (sc_vc_coding_*)
  INTEGER                          :: type              ! Type of v. coord. (sc_vc_type_*)
  INTEGER                          :: n_layer           ! Number of layers (ke)
  INTEGER                          :: lowest_layer      ! Index of lowest layer
  REAL, DIMENSION(4)               :: ref               ! Reference values (e.g. p0)
  REAL, DIMENSION(3)               :: ref_opt           ! Additional reference values (e.g. sleve)
  INTEGER                          :: n_coef            ! Nbr. of coefficients
  REAL, DIMENSION(mx_vc_coef)      :: coef              ! Coefficients defining v. coord.
END TYPE ty_vc_desc


! Field characteristics
!----------------------

! To store product specification
TYPE          :: ty_fld_product
  INTEGER                            :: model           ! Model (sc_model_*)
  INTEGER                            :: pcat            ! Product category (sc_pcat_*)
  INTEGER                            :: suite           ! Suite (sc_suite_*)
END TYPE ty_fld_product


! To store specification of field origin
!   Sources used to set this information:
!   > derived from GRIB (genproc_local, assimilation)
!   > derived from GRIB, can be overwritten by user (experiment, genproc, center, subcenter)
!   > defined by user (model_tag)
TYPE          :: ty_fld_orig
  INTEGER                            :: center           ! Originating center
  INTEGER                            :: subcenter        !             subcenter
  INTEGER                            :: genproc          ! Generating process, pds 4
  INTEGER                            :: genproc_local    !                     local coding 
  INTEGER                            :: assimilation     !                     from assimilation?
                                                         !                     (yes,no,undefined)
  INTEGER                            :: model_tag        ! Model tag (indexed string)
  INTEGER                            :: experiment       ! Experiment tag
END TYPE ty_fld_orig


! To store information about time range unit
TYPE         :: ty_tr_unit
  INTEGER                            :: conversion_factor   ! conversion factor [minutes]
  INTEGER                            :: tr_unit             ! time range unit (GRIB table 4) 
END TYPE ty_tr_unit


! To store specification on EPS member
TYPE          :: ty_fld_epsid
  INTEGER                            :: member          ! member identificator
  INTEGER                            :: nbr_members     ! total number of members
  INTEGER                            :: d_base_date     ! Driving EPS: base date
  INTEGER                            :: d_base_time     ! Driving EPS: base time
  INTEGER                            :: d_nbr_members   ! Driving EPS: number of used members
  INTEGER                            :: d_member        ! Driving EPS: member id
  INTEGER                            :: d_size_cluster  ! Driving EPS: size of cluster
END TYPE ty_fld_epsid


! To store specification of EPS probability field
TYPE          :: ty_fld_prob
  REAL                               :: threshold_low   ! lower treshold 
  REAL                               :: threshold_high  ! higher treshold 
  INTEGER                            :: nbr_eps_members ! nbr. associated EPS members
  LOGICAL                            :: weighted_eps    ! is product derived from weighted eps
END TYPE ty_fld_prob


! Operator characteristics
!-------------------------

! To store characteristics of logical condition
TYPE          :: ty_logical_cond
  INTEGER                                      :: n_token     ! Number of tokens
  CHARACTER(LEN=ln_str), DIMENSION(mx_token)   :: field       ! Token: field
  CHARACTER(LEN=ln_str), DIMENSION(mx_token)   :: operator    !        operator
  REAL, DIMENSION(mx_token)                    :: threshold   !        threshold
END TYPE ty_logical_cond


! To store characteristics of comparison operator
TYPE          :: ty_compare_op
  CHARACTER(LEN=ln_str)                        :: name        ! Operator name
  REAL, DIMENSION(mx_op_coef)                  :: coef        !          coefficients
END TYPE ty_compare_op


! To store characteristics of convolution operator
TYPE          :: ty_conv_op
  CHARACTER(LEN=ln_str)                        :: mask_name       ! Mask name
  INTEGER                                      :: mask_radius     !      radius
  CHARACTER(LEN=ln_str)                        :: mask_weight     !      weight
  INTEGER                                      :: mask_weight_arg !      weight argument
  INTEGER                                      :: op_radius       ! Operator max. radius
  INTEGER                                      :: op_size         !          effective size
  INTEGER, DIMENSION(:,:), POINTER             :: op_mask         !          (i,j) of mask
  REAL, DIMENSION(:), POINTER                  :: op_weight       !          weight
END TYPE ty_conv_op


! To store characteristics of re-gridding operator
TYPE          :: ty_regrid_op
  ! Target grid
  CHARACTER(LEN=ln_str)                        :: associated_tag        ! associated in-core field
  CHARACTER(LEN=ln_str)                        :: grid_type             ! type
  REAL                                         :: start_x               ! start x
  REAL                                         :: start_y               ! start y
  REAL                                         :: d_x                   ! delta x
  REAL                                         :: d_y                   ! delta y
  INTEGER                                      :: n_x                   ! nbr. points in x direction
  INTEGER                                      :: n_y                   ! nbr. points in y direction
  REAL                                         :: pollon                ! lon. of rotated N pole
  REAL                                         :: pollat                ! lat. of rotated N pole
  ! Interpolation alg.
  CHARACTER(LEN=ln_str)                        :: interpolation_method  ! method
  CHARACTER(LEN=ln_str)                        :: neighbourhood_type    ! type of neighbourhood
  REAL                                         :: neighbourhood_radius  ! radius of neighbourhood
  REAL                                         :: exp_scale             ! scale for exp weight
END TYPE ty_regrid_op


! Namelist specifications
!------------------------

! To store specification of output type
! This information specifies characteristics of output data
TYPE          :: ty_out_spec
  CHARACTER(LEN=ln_str)                        :: type   
  CHARACTER(LEN=ln_str)                        :: template   
  CHARACTER(LEN=ln_str)                        :: nomatch_action   
  TYPE(ty_scan_mode)                           :: scanning_mode 
  CHARACTER(LEN=ln_str)                        :: format 
  CHARACTER(LEN=1)                             :: separator
  CHARACTER(LEN=ln_str)                        :: cheader 
  LOGICAL                                      :: fusedt
  REAL                                         :: undefcode
  LOGICAL                                      :: depreciated
  INTEGER                                      :: verbosity
  LOGICAL                                      :: noundef
  LOGICAL                                      :: dupcst
  LOGICAL                                      :: dwhhack 
  CHARACTER(LEN=ln_str)                        :: probinfo
  CHARACTER(LEN=ln_str)                        :: epsminfo
  CHARACTER(LEN=ln_str)                        :: mminfo
  CHARACTER(LEN=ln_str)                        :: text1 
  CHARACTER(LEN=ln_str)                        :: text2 
END TYPE ty_out_spec


! To store specification of output mode
! This information specifies the modus operandi of some algorithms used in the production of the output
TYPE          :: ty_out_mode
  ! ... model tag (for log)
  CHARACTER(LEN=ln_str)                                  :: spec_tag
  ! ... precipitation components
  CHARACTER(LEN=ln_str), DIMENSION(mx_precip_component)  :: precip_all
  CHARACTER(LEN=ln_str), DIMENSION(mx_precip_component)  :: precip_snow
  CHARACTER(LEN=ln_str), DIMENSION(mx_precip_component)  :: precip_rain
  CHARACTER(LEN=ln_str), DIMENSION(mx_precip_component)  :: precip_convective
  CHARACTER(LEN=ln_str), DIMENSION(mx_precip_component)  :: precip_gridscale
  ! ... algorithm for relative humidity
  LOGICAL                                                :: rh_clipped       ! clip relhum to 100%
  ! ... algorithm for eps probability
  LOGICAL                                                :: epsprob_weighted ! weighted probability 
  ! ... algorithm for neighbourhood prob.
  INTEGER                                                :: nbhprob_shape    ! shape 
  INTEGER                                                :: nbhprob_rxy      ! spatial radius
  INTEGER                                                :: nbhprob_rt       ! temporal radius
  INTEGER                                                :: nbhprob_rfade    ! radius of fading zone
END TYPE ty_out_mode


! To store specification of fields to extract or generate
TYPE          :: ty_fld_spec
  INTEGER                                    :: key             ! Unique key (refers to field_key in storage)
  TYPE(ty_fld_product)                       :: product         ! Product associated with field
  INTEGER                                    :: name            ! Name of field (indexed)
  INTEGER                                    :: tag             ! Tag of field (indexed)
  LOGICAL                                    :: ignore_tri      ! Ignore tri when extracting field
  INTEGER                                    :: pspec           ! Product specification (indexed)
  INTEGER                                    :: ltype           ! Level type of field (indexed)
  INTEGER                                    :: level           ! Level of field
  TYPE(ty_fld_epsid)                         :: epsid           ! EPS member identificator
  TYPE(ty_fld_prob)                          :: prob            ! Probability thresholds
  INTEGER, DIMENSION(mx_tags)                :: usetag          ! Tags of allowed parent fields (indexed)
  INTEGER                                    :: useoperator     ! Operator to generate field (indexed)
  LOGICAL                                    :: multi_tlev      ! Touch multiple time levels?
  INTEGER                                    :: set_genproc     ! Set generating process
  INTEGER(KIND=kind_idate)                   :: set_refdate     ! Set reference date
  INTEGER                                    :: set_startperiod ! Set start time of period of validity
  INTEGER                                    :: set_leadtime    ! Set end time of period of validity
  REAL                                       :: scale           ! Scaling factor
  REAL                                       :: offset          ! Offset factor
  INTEGER                                    :: vop             ! Vertical operator (indexed)
  INTEGER                                    :: vop_usetag      !   tag associated fields (indexed)
  INTEGER                                    :: vop_lev         !   target levels (indexed)
  INTEGER                                    :: vop_nlev        !   number of levels
  INTEGER                                    :: top             ! Time operator (indexed)
  INTEGER                                    :: pop             ! Point operator (indexed)
  INTEGER                                    :: hop             ! Lateral operator (indexed)
  INTEGER                                    :: hop_mask        !   associated mask (indexed)
  INTEGER                                    :: rop             ! Re-gridding operator (indexed)
  INTEGER                                    :: merge_with      ! Merge operator (indexed)
  INTEGER                                    :: merge_mask      !   associated mask (indexed)
  INTEGER                                    :: compare_with    ! Compare operator (indexed)
  INTEGER                                    :: compare_fct     !   associated function (indexed)
  INTEGER                                    :: compare_mask    !   associated mask (indexed)
  REAL                                       :: undefcode       ! Code for undefined values in output
END TYPE ty_fld_spec

! To avoid allocating huge sparse arrays, a 2-level structure is used to store
! information of fields to generate:  the root level is indexed by the file,
! the 2nd level by the field and the iteration
TYPE          :: ty_fld_spec_root
  TYPE(ty_fld_spec), DIMENSION(:,:), POINTER :: field => NULL()
END TYPE ty_fld_spec_root


! To store specification on data reduction
TYPE          :: ty_out_dred
  INTEGER                                    :: iteration
  INTEGER                                    :: subset_type
END TYPE ty_out_dred


! Data storage
!-------------

! To create array of field pointers (1d data array)
TYPE          :: ty_1d_field
  REAL, POINTER, DIMENSION(:)                :: data => NULL()
  INTEGER, DIMENSION(2)                      :: level
  INTEGER                                    :: ltype
  INTEGER                                    :: epsm
END TYPE ty_1d_field


! To create array of field pointers (2d data array)
TYPE          :: ty_2d_field
  REAL, POINTER, DIMENSION(:,:)              :: data => NULL()
  INTEGER, DIMENSION(2)                      :: level
  INTEGER                                    :: ltype
  INTEGER                                    :: epsm
END TYPE ty_2d_field


! To store information about INCORE data
TYPE         :: ty_ic_data
  INTEGER                                    :: nbr_field     ! Total number of in-core fields
  CHARACTER(LEN=ln_str), DIMENSION(mx_field) :: tags          ! Tags associated to in-core fields
  INTEGER                                    :: hsurf         ! Index of HSURF
  INTEGER                                    :: zsurf         !       of ZSURF
  INTEGER                                    :: gem           !       of GEM
  INTEGER                                    :: fr_land       !       of FR_LAND
  INTEGER                                    :: rlat          !       of RLAT
  INTEGER                                    :: rlon          !       of RLON
  INTEGER                                    :: swiss_we      !       of SWISS_WE
  INTEGER                                    :: swiss_sn      !       of SWISS_SN
  INTEGER                                    :: boagaw_we     !       of BOAGAW_WE
  INTEGER                                    :: boagaw_sn     !       of BOAGAW_SN
  INTEGER                                    :: boagae_we     !       of BOAGAE_WE
  INTEGER                                    :: boagae_sn     !       of BOAGAE_SN
  INTEGER,DIMENSION(mx_lev)                  :: hfl           !       of HFL
  INTEGER,DIMENSION(mx_lev)                  :: hhl           !       of HHL
END TYPE ty_ic_data


! To store GRIB field and associated meta-informations
TYPE         :: ty_gb_field
  ! Data 
  REAL,DIMENSION(:),POINTER                  :: data           ! Expended data set
  ! Meta-data
  INTEGER                                    :: ednr           ! GRIB edition number
  TYPE(ty_fld_product)                       :: product        ! Product
  TYPE(ty_fld_orig)                          :: origin         ! Origin of field
  CHARACTER(LEN=ln_str), DIMENSION(2)        :: name           ! Name of field and user defined tag
  CHARACTER(LEN=ln_str), DIMENSION(mx_alias) :: alias          ! Aliases (ignore tri)
  INTEGER                                    :: dictionary     ! Active dictionary
  INTEGER                                    :: tbl            ! GRIB table version number
  INTEGER                                    :: parameter      ! GRIB indicator of parameter
  INTEGER                                    :: ltype          ! Level type
  INTEGER, DIMENSION(2)                      :: level          ! Level
  INTEGER, DIMENSION(5)                      :: trange         ! Time range
  INTEGER, DIMENSION(8)                      :: reftime        ! Reference date and time
  TYPE(ty_fld_prob)                          :: prob           ! Probability thresholds
  TYPE(ty_fld_epsid)                         :: epsid          ! Identificator of EPS field
  TYPE(ty_hc_desc)                           :: hcoord         ! Description of horizontal coordinates
  TYPE(ty_vc_desc)                           :: vcoord         ! Description of vertical coordinates
  ! Usage
  LOGICAL                                    :: point_incore   ! Track association with incore data
END TYPE ty_gb_field 


! To store information used when computing new fields
! ... Specifications of field to compute
TYPE    ::  ty_calc_spec  
  TYPE(ty_fld_product)                         :: product     ! Product
  CHARACTER(LEN=ln_str)                        :: name        ! name
  INTEGER                                      :: dictionary  ! Active dictionary
  INTEGER                                      :: ltype       ! level type
  INTEGER                                      :: level       ! level
  TYPE(ty_fld_prob)                            :: prob        ! probability thresholds
  TYPE(ty_fld_epsid)                           :: epsid       ! EPS id
  CHARACTER(LEN=ln_str), DIMENSION(:), POINTER :: usetag      ! tags accepted for parent fields
  CHARACTER(LEN=ln_str)                        :: useoperator ! operator to use
  TYPE(ty_out_mode)                            :: mode        ! mode of production
END TYPE ty_calc_spec
! ... Output data and associated meta information
TYPE     :: ty_calc_outdata     
  REAL, DIMENSION(:,:), POINTER                :: values      ! Output field values
  LOGICAL, DIMENSION(:,:), POINTER             :: flag        !   field flag
  TYPE(ty_fld_prob), POINTER                   :: prob        !   field probability info
  TYPE(ty_fld_orig), POINTER                   :: origin      !   field origin
  INTEGER, DIMENSION(:), POINTER               :: level       !   field level
  INTEGER, DIMENSION(:,:), POINTER             :: trange      !   time range
  INTEGER, POINTER                             :: vref        !   type of reference system
  REAL, DIMENSION(:), POINTER                  :: ngrid       !   associated native grid
  INTEGER, POINTER                             :: spec        !   additional field specifications
END TYPE ty_calc_outdata


! To store data and meta-information associated with each output
!
! -> a set of 2 dimensional fields defined on a the same base grid are 
!    associated with each output. Staggered grids are supported by
!    keeping track of the shift between the base grid and the field
!    native grid. Each 2d field is defined on the same set of grid
!    points - which can be the whole model domain, a sub-domain or 
!    any collection of points -  and is available for a set of
!    validation dates.
! -> fields values are collected in values(:,:,:) array, where:
!       first dim. is for location index 
!       second dim. is for field index 
!       third dim. is for validation date index 
!    The number of locations is given by nbr_gp, the number of fields
!    by nbr_field and the number of validation dates by nbr_vdate. The
!    dimensions of the array is at least nbr_gp x nbr_field x nbr_vdate,
!    but may be larger.
!    Initialization status for a specific field and date is tracked
!    with the defined(:,:) array, where the first dim. is for the
!    field index and the second dim. for the validation date index.
! -> It is also possible to set a flag associated with each element
!    of the value array, by using the flag(:,:,:) array; the 
!    meaning of this flag is context dependent (e.g. it is used
!    to keep track of the type of temperature used in the computation 
!    of t_2m_s). flag(:,:,:) is only allocated when this information
!    is used in the final output (e.g. PSA_TABLE); in this case, the
!    value of flag_is_required is et to true.
! -> any multi-levels field can be reconstructed with the information 
!    of level_idx. The array level_idx(:,1:nbr_level(:)) points
!    to the set of associated levels, sorted in ascending order.
! -> any EPS field can be reconstructed with the information of 
!    eps_member_idx. The array eps_member_idx(:,1:nbr_eps_member(:)) 
!    points to the set of associated members, sorted in ascending order.
! -> characteristics of fields are documented in:
!      field_key, field_is_passive, field_product, field_origin,
!      field_name, field_dictionary, field_trange, 
!      field_ltype, field_level, field_epsid, field_prob,
!      field_vref, field_ngrid, field_undefcode
!      Comments: 
!        field_key: unique field identifier, used to track fields
!                   within the generate iterations
!        field_is_passive: passive fields can not be used as parent
!                          fields and are automatically transfered
!                          from one iteration to the next
!        field_name(1,:): field name,  field_name(2,:): user defined tag
!        field_dictionary: field dictionary used to interpret field_name
!        field_vref: type of reference system for 2d vector field,
!                    value according to bit 5 of GRIB code table 7
!        field_ngrid: transformation to compute field native grid
!                     from model base grid defined by grid_hcoord
!                     (1): value of i-shift  (2): value of j-shift
!                     (3): status (0 means identity transformation)
!        field_undefcode: code used in output for undefined values
! -> field transformations are documented in:
!      field_scale, field_offset
!      field_hop, field_hop_mask, field_top, field_pop, 
!      field_vop, field_vop_nlev, field_vop_lev, field_vop_usetag 
! -> validation dates are documented in:
!      pout_active, field_vdate
!      Comments:
!        pout_active: set to false when associated validation date is 
!                     not transferred into output file
!        field_vdate: validation date, field_vdate(1:nbr_vdate) is
!                     sorted in ascending order
! -> locations are documented in:
!      gp_coord, gp_idx, gp_lat, gp_lon, gp_domain
!      loc_used, loc_idx, loc_tag, loc_h, loc_region.
!      Comments:
!        gp_domain: specification of sub-domain for output
!                (1): imin    (2): jmin  
!                (3): imax    (4): jmax
!                (5): iincr   (6): jincr
!        gp_coord: coordinates of grid point n in input data array
!                (1,n): coordinate i (west-east); 
!                (2,n): coordinate j (south-north);
!                (3,n): position in input data array 
!        gp_idx: position of grid point (i,j) in value array 
!                (only set when gp_domain is defined!)
!        gp_lat,
!        gp_lon: geographical coordinates of grid point [millidegree]
!        loc_idx: index in locations array of associated additional
!                 information.
! -> All fields must be associated with the same horizontal mesh, described by
!      grid_hcoord
!    All fields must not necessarily be associated with the same vertical mesh;
!    this is only the case when some fields are defined on model levels. The
!    common vertical mesh, if it is defined, is decribed by
!      grid_vcoord
!    The status of grid_vcoord information is kept in
!      grid_vcoord_set
! -> when defined, a common data origin is described by
!      center, subcenter, model_tag and experiment 
! -> when defined, a common  reference date is described by
!      base_date
! -> the associated output file is univoquely defined by ofile_idx,
!    its index in out_paths array. 
!    Additional information on output file is available:
!      output file name (ofile_name),
!      bogus output, i.e. error occured while building file (ofile_bogus),
!      completness of collected values (ofile_complete),
!      first write in output file (ofile_firstwrite),
!      use file name postfix as long as data is not complete (ofile_usepostfix),
!      write in append mode (ofile_replace)
!      indices of input files contributing to collected values (ifile_used)
! -> order of fields appearance in output file follows field_idx(:) ,
!    last_field_sorted is used to keep track of already sorted fields
!
TYPE                                            :: ty_out_store
  INTEGER                                          :: center
  INTEGER                                          :: subcenter
  CHARACTER(LEN=ln_str)                            :: model_tag
  INTEGER                                          :: experiment
  INTEGER(KIND=kind_idate)                         :: base_date
  INTEGER                                          :: ofile_idx
  CHARACTER(LEN=ln_str)                            :: ofile_name
  LOGICAL                                          :: ofile_bogus
  LOGICAL                                          :: ofile_complete
  LOGICAL                                          :: ofile_firstwrite
  LOGICAL                                          :: ofile_usepostfix
  LOGICAL                                          :: ofile_replace
  LOGICAL, DIMENSION(mx_ifile)                     :: ifile_used
  TYPE(ty_vc_desc)                                 :: grid_vcoord
  INTEGER                                          :: grid_vcoord_set
  TYPE(ty_hc_desc)                                 :: grid_hcoord
  INTEGER                                          :: nbr_gp
  INTEGER, DIMENSION(:,:), POINTER                 :: gp_coord
  INTEGER, DIMENSION(:,:), POINTER                 :: gp_idx
  REAL, DIMENSION(:), POINTER                      :: gp_lat
  REAL, DIMENSION(:), POINTER                      :: gp_lon
  INTEGER, DIMENSION(6)                            :: gp_domain
  LOGICAL                                          :: loc_used
  INTEGER, DIMENSION(mx_gp)                        :: loc_idx
  CHARACTER(LEN=ln_loc), DIMENSION(mx_gp)          :: loc_tag
  CHARACTER(LEN=ln_reg), DIMENSION(mx_gp)          :: loc_region
  INTEGER, DIMENSION(mx_gp)                        :: loc_h
  INTEGER                                          :: nbr_vdate
  INTEGER(KIND=kind_idate), DIMENSION(mx_time)     :: field_vdate
  LOGICAL, DIMENSION(mx_time)                      :: pout_active
  INTEGER                                          :: nbr_field
  INTEGER, DIMENSION(:), POINTER                   :: field_idx
  INTEGER                                          :: last_field_sorted
  INTEGER, DIMENSION(:), POINTER                   :: field_key
  LOGICAL, DIMENSION(:), POINTER                   :: field_is_passive
  TYPE(ty_fld_product), DIMENSION(:), POINTER      :: field_product
  TYPE(ty_fld_orig), DIMENSION(:), POINTER         :: field_origin
  CHARACTER(LEN=ln_str), DIMENSION(:,:), POINTER   :: field_name
  INTEGER, DIMENSION(:), POINTER                   :: field_dictionary
  INTEGER, DIMENSION(:,:,:), POINTER               :: field_trange
  INTEGER, DIMENSION(:), POINTER                   :: field_ltype
  INTEGER, DIMENSION(:,:), POINTER                 :: field_level
  TYPE(ty_fld_epsid), DIMENSION(:), POINTER        :: field_epsid
  TYPE(ty_fld_prob), DIMENSION(:), POINTER         :: field_prob
  INTEGER, DIMENSION(:), POINTER                   :: field_vref
  REAL, DIMENSION(:,:), POINTER                    :: field_ngrid
  REAL, DIMENSION(:), POINTER                      :: field_undefcode
  REAL, DIMENSION(:), POINTER                      :: field_scale
  REAL, DIMENSION(:), POINTER                      :: field_offset
  INTEGER, DIMENSION(:), POINTER                   :: field_hop
  INTEGER, DIMENSION(:), POINTER                   :: field_hop_mask
  INTEGER, DIMENSION(:), POINTER                   :: field_vop
  INTEGER, DIMENSION(:), POINTER                   :: field_vop_nlev
  INTEGER, DIMENSION(:), POINTER                   :: field_vop_lev
  INTEGER, DIMENSION(:), POINTER                   :: field_vop_usetag
  INTEGER, DIMENSION(:), POINTER                   :: field_top
  INTEGER, DIMENSION(:), POINTER                   :: field_pop
  INTEGER, DIMENSION(:), POINTER                   :: nbr_level
  INTEGER, DIMENSION(:,:), POINTER                 :: level_idx
  INTEGER, DIMENSION(:), POINTER                   :: nbr_eps_member
  INTEGER, DIMENSION(:,:), POINTER                 :: eps_member_idx
  REAL, DIMENSION(:,:,:), POINTER                  :: values
  LOGICAL, DIMENSION(:,:), POINTER                 :: defined
  LOGICAL, DIMENSION(:,:,:), POINTER               :: flag
  LOGICAL                                          :: flag_is_required
END TYPE ty_out_store



! Interfaces
!===========
INTERFACE OPERATOR(==)
  MODULE PROCEDURE are_scan_mode_equal
END INTERFACE

INTERFACE OPERATOR(/=)
  MODULE PROCEDURE do_scan_mode_differ
END INTERFACE



CONTAINS

                     !*******************************!
                     ! ==> Operators overloading <== !
                     !*******************************!



!+****************************************************************************
FUNCTION are_scan_mode_equal(scan_mode_1, scan_mode_2)
!=============================================================================
!
! Overload '==' for ty_scan_mode
!
!------------------------------------------------------------------------------
  ! Dummy arguments
  TYPE(ty_scan_mode), INTENT(IN)     :: scan_mode_1
  TYPE(ty_scan_mode), INTENT(IN)     :: scan_mode_2
  LOGICAL                            :: are_scan_mode_equal


  IF ( scan_mode_1%valid         .EQV. scan_mode_2%valid         .AND.    &
       scan_mode_1%defered       .EQV. scan_mode_2%defered       .AND.    &
       scan_mode_1%positive_i    .EQV. scan_mode_2%positive_i    .AND.    &
       scan_mode_1%positive_j    .EQV. scan_mode_2%positive_j    .AND.    &
       scan_mode_1%consecutive_i .EQV. scan_mode_2%consecutive_i     ) THEN
    are_scan_mode_equal = .TRUE.
  ELSE
    are_scan_mode_equal = .FALSE.
  ENDIF

END FUNCTION are_scan_mode_equal


!+****************************************************************************
FUNCTION do_scan_mode_differ(scan_mode_1, scan_mode_2)
!=============================================================================
!
! Overload '/=' for ty_scan_mode
!
!------------------------------------------------------------------------------
  ! Dummy arguments
  TYPE(ty_scan_mode), INTENT(IN)     :: scan_mode_1
  TYPE(ty_scan_mode), INTENT(IN)     :: scan_mode_2
  LOGICAL                            :: do_scan_mode_differ


  do_scan_mode_differ = .TRUE.
  IF ( are_scan_mode_equal(scan_mode_1, scan_mode_2) ) do_scan_mode_differ = .FALSE.


END FUNCTION do_scan_mode_differ



END MODULE fxtr_definition
