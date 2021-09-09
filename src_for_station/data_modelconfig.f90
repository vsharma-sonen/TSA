!+ Data module for the configuration of spatial and time discretization
!------------------------------------------------------------------------------

MODULE data_modelconfig

!------------------------------------------------------------------------------
!
! Description for TSA version
!  This is data_modelconfig from COSMO version 5.07 without modifications
!
!------------------------------------------------------------------------------
!
! Description:
!  This module contains all variables (scalars and arrays) concerned with the
!  configuration of the horizontal and vertical grid, the time discretization
!  and the reference atmosphere. These are
!    - the vertical coordinate parameters and related variables
!    - horizontal and vertical sizes of the arrays and related variables
!    - start- and end-indices for the computations in the horizontal layers
!    - constants for the horizontal rotated grid and related variables
!    - variables for the time discretization and related variables
!
!  The arrays are declared as allocatable arrays and are allocated in the
!  setup and deallocated in the cleanup of the model. The scalar variables
!  are either part of a NAMELIST (read in the setup) or are initialized
!  in the routine calculate_consts.
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
! 1.5        1998/06/29 Guenther Doms
!  Definition of a new global array hhlr (height of half levels ref. to z=0)
! 1.30       1999/06/24 Matthias Raschendofer
!  Definition of a new INTEGER-parameter kcm used for allocation of canopy fields.
! 1.34       1999/12/10 Ulrich Schaettler
!  Added variables vhmx_vol and vhmx_cfl
! 1.39       2000/05/03 Ulrich Schaettler
!  Renamed variables concerned with latitude and longitude from phi, rla to
!  lat and lon. Included variables for specifying layer index corresponding
!  to certain pressures (has been in data_diagnostics before).
! 2.8        2001/07/06 Ulrich Schaettler
!  Added new variables for multi-layer soil model
! 2.11       2001/09/28 Ulrich Schaettler
!  Renamed a variable for multi-layer soil model
! 2.14       2002/02/15 Ulrich Schaettler
!  New variables for the smooth level vertical coordinate (SLEVE)
! 3.6        2003/12/11 Reinhold Schrodin
!  Added (and changed) some variables for multi-layer soil model
! 3.14       2005/01/25 Ulrich Schaettler
!  Introduced variable ke_rd: lowest level with Rayleigh damping
! 3.18       2006/03/03 Ulrich Schaettler
!  Editorial Changes
! 3.19       2006/04/25 Ulrich Schaettler
!  New variable to store the depth of soil half levels (for NetCDF I/O)
! 3.21       2006/12/04 Jochen Foerstner, Klaus Stephan, Burkhardt Rockel
!  New variables betagw, beta2sw, beta2gw for treatment of sound and gravity waves
!  New index variables for LM mainlevels on 950, 700 HPa
!  Introduced new Namelist parameter polgam
! V4_5         2008/09/10 Guenther Zaengl
!  Added parameters for new reference atmosphere
! V4_11        2009/11/30 Ekaterina Machulskaya, Guenther Zaengl
!  Introduced ke_snow for number of layers in snow model (EM)
!  Introduced alphaass for Asselin filter weight (GZ)
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_24        2012/06/22 Michael Baldauf
!  Introduced variables dx_min, dy_min, dz_min
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  MESSy interface introduced
!  lalloc_XXX added: these switches avoid the necessary distinction
!                    between IF ALLOCATED and IF ASSOCIATED
! V4_28        2013/07/12 Ulrich Schaettler
!  Introduced endlon_tot, endlat_tot as global variables
!  Removed vertical grid and reference atmosphere variables and put them
!   to vgrid_refatm_utils
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Introduced index variables for the COSMO humidity tracers
! V5_1         2014-11-28 Ulrich Schaettler, Oliver Fuhrer
!  Introduced 1D field hhl_prof(:) for providing a special hhl profile to 
!   all subdomains for calculations formerly refering to hhlr
!  Replaced ireals by wp (working precision) (OF)
! V5_4f        2017-09-01 Ulrich Blahak, Ulrich Schaettler
!  Added lcalc_XXXtrack: switches to avoid unnecessary computations
!   of cell track variables, which are not in the output list (UB)
!  Removed obsolete switches lmulti_layer and associated lalloc_w_g3x (US)
! V5_5         2018-02-23 Ulrich Blahak
!  Added lcalc_echotopp and lcalc_echotopz to trigger computation of
!   dBZ echo tops in hPa and/or in meters MSL (UB)
!  Added lcalc_mconvtrack to trigger computation of horizontal moisture convergence track
! V5_5b        2018-10-29 Ulrich Schaettler
!  Introduced global logical lcalc_sometracks to switch off track calculation 
!   completedly (set in organize_data)
! V5_6b        2019-10-16 Pavel Khain
!  Introduced indices for tracers used by ICON-ART or MACC aerosol mixing ratios
! V5_7         2020-02-21 Astrid Kerkweg, Ulrich Schaettler
!  Introduce the index kstart_moist to switch off saturation adjustment in the
!    higher stratosphere (important for high atmosphere version only) (AK)
!  Removed iintegers (US)
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
USE kind_parameters, ONLY :   &
    wp           ! KIND-type parameter for real variables

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Global (i.e. public) Declarations:

! 1. vertical coordinate parameters and related variables
! -------------------------------------------------------

  REAL  (KIND=wp), ALLOCATABLE ::           &
    hhl_prof(:)     ! a special hhl-profile, which is synchronized over all 
                    ! subdomains. It will be allocated from 0:ke1, 
                    ! hhl_prof(0) contains an approximate zflat-value

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------

  INTEGER  ::     &
    ! number of grid points for the total domain
    ie_tot,       & ! number of grid points in zonal direction
    je_tot,       & ! number of grid points in meridional direction
    ke_tot,       & ! number of grid points in vertical direction
    nlandpoints_tot,  & ! number of land points in the grid

    ! number of grid points for this domain
    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ke,           & ! number of grid points in vertical direction
    ke_soil,      & ! number of layers in the multi-layer soil model
    ke_snow,      & ! number of layers in multi-layer snow model
    nlandpoints,  & ! number of land points in the grid

    ke1,          & ! KE+1
    ke_rd,        & ! lowest level with Rayleigh-damping
    ieje,         & ! IE*JE
    iejeke,       & ! IE*JE*KE
    ieke,         & ! IE*KE
    ie_max,       & ! Max. of ie on all processors
    je_max,       & ! Max. of je on all processors

    ! start- and end-indices for vertically bounded regions
    kcm             ! index of the lowest model layer, higher than the canopy
 
  ! adopted from ICON for COSMO high atmosphere version
  ! the saturation pressure parametrisations are not valid anymore in the higher
  ! atmosphere. Therefore a start index is defined from which on downward 
  ! saturation adjustment is calculated. By default this value is 1 

  INTEGER  :: kstart_moist = 1 


! 2a. Variables for the new multi-layer soil model
! --------------------------------------------------------------------

  REAL (KIND=wp), ALLOCATABLE       ::           &
    czmls(:),     & ! depth of the main soil layers in meters
    czhls(:)        ! depth of the half soil layers in meters

  INTEGER,        ALLOCATABLE       ::           &
    msoilgrib(:)    ! grib coded depth of main soil levels in centimeters
                    ! (careful: the first level will be coded with 1,
                    !           but is in the depth of 0.5 cm!)

! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the 
!    forecast for the prognostic variables in a horizontal layer.
!    Note, that the indices for the wind-speeds u and v differ from 
!    the other ones because of the use of the staggered Arakawa-C-grid.
!    
  INTEGER  ::     &
!   zonal direction
    istart,       & ! start index for the forecast of w, t, qd, qw and pp
    iend,         & ! end index for the forecast of w, t, qd, qw and pp
    istartu,      & ! start index for the forecast of u
    iendu,        & ! end index for the forecast of u
    istartv,      & ! start index for the forecast of v
    iendv,        & ! end index for the forecast of v
    istartpar,    & ! start index for computations in the parallel program
    iendpar,      & ! end index for computations in the parallel program

!   meridional direction
    jstart,       & ! start index for the forecast of w, t, qd, qw and pp
    jend,         & ! end index for the forecast of w, t, qd, qw and pp
    jstartu,      & ! start index for the forecast of u
    jendu,        & ! end index for the forecast of u
    jstartv,      & ! start index for the forecast of v
    jendv,        & ! end index for the forecast of v
    jstartpar,    & ! start index for computations in the parallel program
    jendpar         ! end index for computations in the parallel program

! 4. constants for the horizontal rotated grid and related variables
! ------------------------------------------------------------------

  REAL  (KIND=wp)                  ::           &
    pollon,       & ! longitude of the rotated north pole (in degrees, E>0)
    pollat,       & ! latitude of the rotated north pole (in degrees, N>0)
    polgam,       & ! angle between the north poles of the systems
    dlon,         & ! grid point distance in zonal direction (in degrees)
    dlat,         & ! grid point distance in meridional direction (in degrees)
    startlon_tot, & ! transformed longitude of the lower left grid point
                    ! of the total domain (in degrees, E>0)
    startlat_tot, & ! transformed latitude of the lower left grid point
                    ! of the total domain (in degrees, N>0)
    endlon_tot,   & ! transformed longitude of the upper right grid point
                    ! of the total domain (in degrees, E>0)
    endlat_tot,   & ! transformed latitude of the upper right grid point
                    ! of the total domain (in degrees, N>0)
    startlon,     & ! transformed longitude of the lower left grid point
                    ! of this subdomain (in degrees, E>0)
    startlat,     & ! transformed latitude of the lower left grid point
                    ! of this subdomain (in degrees, N>0)
    eddlon,       & ! 1 / dlon
    eddlat,       & ! 1 / dlat
    edadlat,      & ! 1 / (radius of the earth * dlat)
    dlonddlat,    & ! dlon / dlat
    dlatddlon,    & ! dlat / dlon
    degrad,       & ! factor for transforming degree to rad
    raddeg,       & ! factor for transforming rad to degree
    dx_min,       & !
    dy_min,       & !
    dz_min          !

! 5. variables for the time discretization and related variables
! --------------------------------------------------------------

  REAL  (KIND=wp)                  ::           &
    dt,           & ! long time-step
    ed2dt,        & ! 1 / (2 * dt)
    dt2,          & ! 2 * dt
    dtdeh,        & ! dt / 3600 seconds
    epsass,       & ! eps for the Asselin-filter
    alphaass,     & ! weight for Williams 2009 modification (0.5 < alpha <= 1)
    betasw,       & ! beta-variable for treatment of soundwaves    in w
    betagw,       & ! beta-variable for treatment of gravity-waves in w
    beta2sw,      & ! beta-variable for treatment of soundwaves    in p*, T*
    beta2gw,      & ! beta-variable for treatment of gravity-waves in p*, T*
    vhmx_vol,     & ! maximum absolute horizontal wind in total model domain
    vhmx_cfl        ! maximum absolute horizontal wind velocity from CFL

  INTEGER  ::     &
    nehddt          ! 3600 seconds / dt

! 7. Layer index corresponding to a specified pressure
!    and other organizational variables
! ----------------------------------------------------

  INTEGER  ::     &
    klv950,       & ! k index of the LM-mainlevel, on 950 HPa
    klv850,       & ! k index of the LM-mainlevel, on 850 HPa
    klv800,       & ! k index of the LM-mainlevel, on 800 HPa
    klv700,       & ! k index of the LM-mainlevel, on 700 HPa
    klv500,       & ! k index of the LM-mainlevel, on 500 HPa
    klv400,       & ! k index of the LM-mainlevel, on 400 HPa
    klv300          ! k index of the LM-mainlevel, on 300 HPa

! 8. Organizational variables to handle the COSMO humidity tracers
!    and other related things
! ----------------------------------------------------

  INTEGER, SAVE ::   &
    ! these are the indices used for the COSMO humidity tracers
    ! (including the 2-moment scheme)
    idt_qv  = -99,   &
    idt_qc  = -99,   &
    idt_qi  = -99,   &
    idt_qr  = -99,   &
    idt_qs  = -99,   &
    idt_qg  = -99,   &
    idt_qh  = -99,   &
    idt_qni = -99,   &
    idt_qnc = -99,   &
    idt_qnr = -99,   &
    idt_qns = -99,   &
    idt_qng = -99,   &
    idt_qnh = -99

  INTEGER, SAVE ::   &
    idt_ncn = -99

  INTEGER, PRIVATE  :: i

  INTEGER, SAVE :: &
    ! these are the indices used for the ICON-ART aerosol mixing ratios:
    idt_icon_art(1:3) = (/ (-99, i=1,3) /), &
    ! these are the indices used for the MACC aerosol mixing ratios:
    idt_aermr(1:11) = (/ (-99, i=1,11) /)

  INTEGER                        ::  &
    ntrcr_radprog                   !< Number of prognostic tracers for radiation

  ! DEFINE LOGICAL to avoid differentiation of " IF ALLOCATED" and
  ! "IF ASSOCIATED"  between pure COSMO and COSMO/MESSY setup
  LOGICAL ::                       &
    lalloc_h_ice     = .FALSE.,    &
    lalloc_t_cl      = .FALSE.,    &
    lalloc_t_s_bd    = .FALSE.,    &
    lalloc_prr_gsp   = .FALSE.,    &
    lalloc_prr_con   = .FALSE.,    &
    lalloc_prs_gsp   = .FALSE.,    &
    lalloc_prs_con   = .FALSE.,    &
    lalloc_prg_gsp   = .FALSE.,    &
    lalloc_tke       = .FALSE.,    &
    lalloc_prh_gsp   = .FALSE.,    &  
    lalloc_tinc_lh   = .FALSE.,    & 
    lalloc_t_ice     = .FALSE.,    &
    lalloc_pertstoph = .FALSE.


  LOGICAL ::                       &
    lcalc_lpitrack   = .FALSE.,    &
    lcalc_uheltrack  = .FALSE.,    &
    lcalc_rottrack   = .FALSE.,    &
    lcalc_dbztrack   = .FALSE.,    &
    lcalc_viltrack   = .FALSE.,    &
    lcalc_echotopz   = .FALSE.,    &
    lcalc_echotopp   = .FALSE.,    &
    lcalc_mconvtrack = .FALSE.,    &
    lcalc_sometracks = .FALSE.

!=======================================================================
END MODULE data_modelconfig
