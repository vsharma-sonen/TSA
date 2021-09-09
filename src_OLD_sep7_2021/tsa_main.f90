!+ Source program "terra" - main program of Terra Stand Alone (TSA)
!------------------------------------------------------------------------------

PROGRAM tsa_main

!------------------------------------------------------------------------------
!
! Description:
!  This is a new version of TSA with the "block structure" of COSMO / ICON
!
!  The program terra is the main program of the StandAlone version of the
!  soil module terra which is in sfc_terra, extracted from the Local-Model of 
!  the COSMO Consortium.
!
!  The program organizes all data, environment and parameters needed for
!  a standalone execution of the terra module.
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
!  - Gerd Vogel and Juergen Helmert (DWD): Brooks and
!    Coorey paramterization of drainage and diffusion,
!    revised formulation of infiltration, root density 
!    distribution
!  - Alexander Block (BTU): soil heat conductivity 
!    dependant of actual soil moisture content
!  - Gerd Vogel, Eric Jaeger (ETHZ), Catherine Meissner (FZK) and others:
!    carefull testing and bug removal
! V4.13      2010/08/01 Guy de Morsier, MeteoSwiss (GDM)
!  Revised code to correspond to FieldExtra and COSMO version 4.13
! V5.01      2015-12-01 Yiftach Ziv, IMS (XYZ)
!  Revised code to correspond to COSMO version 5.01
!  Eliminated use of tracers, since TSA cannot operate with them.
!  Arranged code to adhere to coding standards.
! V5.03      2016-04-01 Yiftach Ziv, IMS (XYZ)
!  Added rlat_geo, rlon_geo to represent geographical lat,lon
!   of non-rotated system to be used in subroutine vegadapt
!   following Julian Todter (JT) of GUF
! V5.07      2020/02/21     Doerte Liermann, Ulrich Schaettler
!  Replaced data_parameters by kind_parameters and removed iintegers
!  Added additional variables for block structure
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
!  A) data modules from COSMO
!------------------------------------------------------------------------------

USE kind_parameters, ONLY :   &
    wp                          ! KIND-type parameter for real variables

!------------------------------------------------------------------------------

USE data_block_fields, ONLY:                                                 &
    fr_land_b, u_m_b, v_m_b, t_b, qv_s_b, t_g_b, tcm_b, tch_b

!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &


! 2. physical constants and related variables
! -------------------------------------------
    r_d,          & ! gas constant for dry air
    rvd_m_o,      & ! r_v/r_d - 1
    cp_d,         & ! specific heat of dry air at constant pressure
    rdocp,        & ! r_d / cp_d
    lh_v            ! latent heat of vapourization

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------
    p0         ,    & ! base state pressure                           (Pa) 

! 2. external parameter fields                                        (unit)
! ----------------------------
    plcov      ,    & ! fraction of plant cover                         --
    rootdp     ,    & ! depth of the roots                            ( m  )
    lai        ,    & ! leaf area index of plants                       --
    llandmask  ,    & ! landpoint mask                                  --

! 3. prognostic variables                                             (unit)
! -----------------------
    u          ,    & ! zonal wind speed                              ( m/s )
    v          ,    & ! meridional wind speed                         ( m/s )
    t          ,    & ! temperature                                   (  k  )

! 5. fields for surface values and soil model variables               (unit )
! -----------------------------------------------------
    ps        ,     & ! surface pressure                              ( pa  )
    t_s       ,     & ! temperature of the ground surface             (  K  )
    t_g       ,     & ! weighted surface temperature                  (  K  )
    qv_s      ,     & ! specific humidity at the surface              (kg/kg)
    prr_gsp     ,   & ! precipitation rate of rain, grid-scale        (kg/m2*s)
    prs_gsp     ,   & ! precipitation rate of snow, grid-scale        (kg/m2*s)
    tch         ,   & ! turbulent transfer coefficient for heat       ( -- )
    sobs        ,   & ! solar radiation at the ground                 ( W/m2)
    thbs        ,   & ! thermal radiation at the ground               ( W/m2)

! 7. fields for model output and diagnostics                          (unit )
! ---------------------------------------------------------------
    runoff_s    ,   & ! surface water runoff; sum over forecast       (kg/m2)
    runoff_g    ,   & ! soil water runoff; sum over forecast          (kg/m2)
    shfl_s      ,   &
    lhfl_s      ,   &
    lhfl_bs           ! average latent heat flux from bare soil evap. ( W/m2)

! end of data_fields

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------
    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ke              ! number of grid points in vertical direction

!------------------------------------------------------------------------------

USE data_parallel,      ONLY:  &
    my_cart_id        ! rank of this subdomain in the global communicator

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

    ntstep,       & ! actual time step
    nnow  ,       & ! corresponds to ntstep 
    nnew  ,       & ! corresponds to ntstep + 1
    nproma,       & ! block size for physical parameterizations
    nlastproma,   & ! size of last block
    nblock          ! number of blocks

!------------------------------------------------------------------------------
!   B) Data module, utilities and routines for TSA
!------------------------------------------------------------------------------

USE tsa_data      ! use all variables - especially designated for TSA

USE tsa_setup,          ONLY:                                                 &
    read_namelist, allocate_fields, allocate_block_fields, init_variables,    &
    clean_up

USE tsa_lmparam,        ONLY:                                                 &
    parturs_newblock, gen_area_indices, vegadapt, stats, model_abort

USE tsa_input,          ONLY:                                                 &
    read_const_fields, read_const_fields_icon,                                &
    read_initial_fields, read_initial_fields_icon,                            &
    read_metforc

USE tsa_output,         ONLY:                                                 &
    grbout, grbout_eccodes, binout, ascout
 
USE tsa_sfc_interface,  ONLY:                                                 &
    tsa_sfc_init, tsa_sfc_init_copy, tsa_sfc_finalize, tsa_sfc_organize

!------------------------------------------------------------------------------

!==============================================================================

IMPLICIT NONE

!==============================================================================

!------------------------------------------------------------------------------
! Declarations
!------------------------------------------------------------------------------

  CHARACTER (LEN=250)       :: yzerrmsg

  ! Local scalars:
  ! -------------

  INTEGER        ::  &
  !   Indices
    nt             , & ! time-level for integration
    k              , & ! loop index for snow layers
    i              , & ! loop index in x-direction              
    j              , & ! loop index in y-direction              
    ipend          , & ! 
    im1, jm1       , & ! i-1, j-1
    ib             , & ! loop index for block structure: ib=1,nblock
    izerror            ! error status variable

  !   Date handling
  INTEGER        ::  &
    jahrestag      , & !
    jahrestag_old  , & !
    month          , & !
    month_old      , & !
    hsurf_size     , & !
    runoff_size        !

  REAL (KIND=wp)           ::  &
    zst              , & ! dummy argument for Stmt. function
    zsge             , & ! dummy argument for Stmt. function
    zsp              , & ! dummy argument for Stmt. function
    zuv              , & ! wind velocity in lowest atmospheric layer   (m/s)
    zexcoeff         , & !
    zrho_atm             ! air density of lowest atmospheric layer     (kg/m**3)

  REAL (KIND=wp), ALLOCATABLE, DIMENSION(:,:) ::    &
    tzuv             , & !
    tzrho                !

!==============================================================================
! Start of Main program 
!==============================================================================

!========================================================
! Section 1: Initializations
!========================================================

! Reading Namelist configuration
  CALL read_namelist()                ! from terra_io

! Allocation of all required LM fields
  CALL allocate_fields()

! Read time constant 
  CALL init_variables()               ! from terra_lmenv

!DL
! Allocation of all required LM fields in block form
  CALL allocate_block_fields()
!DL

  SELECT CASE (TRIM(ymodel))
  CASE ('COSMO')
    IF (lext_monthly) THEN
       CALL read_const_fields(1)        ! from terra_io
    ELSE 
       CALL read_const_fields()         ! from terra_io
    ENDIF

    CALL read_initial_fields()          ! from terra_io

  CASE ('ICON')
    IF (lext_monthly) THEN
       CALL read_const_fields_icon(1)
    ELSE
       CALL read_const_fields_icon()
    ENDIF

    CALL read_initial_fields_icon()
  END SELECT

! read gen_area_indices
  CALL gen_area_indices                 ! from tsa_lmparam

  IF (ntype_output<=3) THEN
!DL GRIB output routine depending on yform_write
    IF (yform_write == 'grb1') THEN
      CALL grbout(0)                   ! from terra_io
    ELSE
      CALL grbout_eccodes(0)           ! from terra_io
    ENDIF
  ENDIF
  ! allocate tzuv and tzrho
  IF (ymodel == 'COSMO') THEN
    ALLOCATE(tzuv(ie,je), tzrho(ie,je))
  ELSE
    ALLOCATE(tzuv(nproma,nblock), tzrho(nproma,nblock))
  ENDIF

! initialize old month and year
  jahrestag_old=-1
  month_old=-1

! initialize sensible and latent heat fluxes at surface
  shfl_s (:,:) = 0.0_wp
  lhfl_bs(:,:) = 0.0_wp
  lhfl_s (:,:) = 0.0_wp

!DL initialize block structure
  CALL tsa_sfc_init()

!========================================================
! Section 2: Time steping
!========================================================

!!DO nt = 0, ntstep_max-1
  DO ntstep = 0, ntstep_max-1

!!  ntstep=nt
    WRITE(6,*) "Step",ntstep


    !--------------------------------------------------------
    ! Section 2.1: Reading meteorological forcing
    !--------------------------------------------------------

     CALL read_metforc(ntstep,jahrestag,month)                   ! from terra_io


    !--------------------------------------------------------
    ! Section 2.2: Update plant parameters LAI, PLCOV and ROOTDP
    !--------------------------------------------------------

    IF ( lvegadapt .AND. (jahrestag/=jahrestag_old)) THEN    
      CALL vegadapt(jahrestag,hsurface,lat,             &   ! from terra_lmparam
           lai_mn,lai_mx,plcov_mn,plcov_mx,rootdp0,     &
           lai,plcov,rootdp)
      WRITE(*,*) "Vegetation adapted for ",jahrestag
      CALL gen_area_indices
      jahrestag_old=jahrestag
    ENDIF
    IF (lext_monthly .AND. (month/=month_old)) THEN
      CALL read_const_fields(month)
      CALL gen_area_indices
      month_old=month
    ENDIF

    IF (lcalc) THEN

      !--------------------------------------------------------
      ! Section 2.3: Loop over blocks
      !--------------------------------------------------------

      DO ib = 1, nblock
        IF (ib == nblock) THEN
           ipend = nlastproma
        ELSE
           ipend = nproma
        ENDIF

        CALL tsa_sfc_init_copy  (1, ib, ipend,  ierror, yerror)
        IF (ierror /= 0) THEN
           CALL model_abort (my_cart_id, 100+ierror, yerror, 'tsa_sfc_init_copyToBlock')
        ENDIF

        ! Calculate transfer coefficients
        CALL parturs_newblock (zh=dz_u, iblock=ib, nvec=ipend, &
                               fr_land = fr_land_b(:),         &
                               z0      = z0_b(:),              &
                               u       = u_m_b(:,ke),          &
                               v       = v_m_b(:,ke),          &
                               t       = t_b(:,ke),            &
                               qv      = qv_b(:,ke),           &
                               qv_s    = qv_s_b(:),            &
                               t_g     = t_g_b(:),             &
                               tcm     = tcm_b(:),             &
                               tch     = tch_b(:)        )

        ! Call soil module TERRA (within tsa_sfc_organize)
        CALL tsa_sfc_organize (ib, ipend,  ierror, yerror)
        IF (ierror /= 0) THEN
           CALL model_abort (my_cart_id, 100+ierror, yerror, 'tsa_sfc_organize')
        ENDIF

        CALL tsa_sfc_init_copy  (2, ib, ipend,  ierror, yerror)
        IF (ierror /= 0) THEN
           CALL model_abort (my_cart_id, 100+ierror, yerror, 'tsa_sfc_init_copyFromBlock')
        ENDIF

      ENDDO

      !--------------------------------------------------------
      ! Section 2.5: Diagnostics and Output
      !--------------------------------------------------------

      ! Diagnose sensible and latent heat flux over land surfaces
      IF (ymodel == 'COSMO') THEN
        DO j=1,je
          jm1 = MAX( 1, j-1 )
          DO i=1,ie 
            IF (llandmask(i,j)) THEN
              im1 = MAX( 1, i-1)
              zuv        = 0.5*SQRT ( (u(i,j,ke,nnew) + u(im1,j,ke,nnew))**2      &
                                     +(v(i,j,ke,nnew) + v(i,jm1,ke,nnew))**2 ) 
              zrho_atm = ps(i,j,ke) / (r_d*  t_g(i,j,nnow)*(1+rvd_m_o*qv_s(i,j,nnow)))
              zexcoeff= tch(i,j)* zuv * zrho_atm
              shfl_s(i,j)= -cp_d * zexcoeff * &
                   (t(i,j,ke,nnew)*( (ps(i,j,nnew)/p0(i,j,ke))**rdocp ) - t_g(i,j,nnew) )
              lhfl_s(i,j)= -lh_v * zexcoeff * &
                   (qv(i,j,ke,nnew)-qv_s(i,j,nnew))
!                  (qv(i,j,ke)-qv_s(i,j,nnew))           ! in tracer form - commented
              tzuv(i,j) = zuv
              tzrho(i,j) = zrho_atm
            ENDIF
          ENDDO
        ENDDO
      ELSE
!DL  ......................ICON............................................
        DO j = 1, nblock
          IF (j == nblock) THEN
            ipend = nlastproma
          ELSE
            ipend = nproma
          ENDIF
          DO i = 1, ipend
            IF (llandmask(i,j)) THEN
                zuv        = 0.5*SQRT ( u(i,j,ke,nnew)**2 + v(i,j,ke,nnew)**2 )
                zrho_atm = ps(i,j,ke) / (r_d*  t_g(i,j,nnow)*(1+rvd_m_o*qv_s(i,j,nnow)))
                zexcoeff= tch(i,j)* zuv * zrho_atm
                shfl_s(i,j)= -cp_d * zexcoeff * &
                     (t(i,j,ke,nnew)*( (ps(i,j,nnew)/p0(i,j,ke))**rdocp ) - t_g(i,j,nnew) )
                lhfl_s(i,j)= -lh_v * zexcoeff * &
                     (qv(i,j,ke,nnew)-qv_s(i,j,nnew))
!                    (qv(i,j,ke)-qv_s(i,j,nnew))           ! in tracer form - commented
                tzuv(i,j) = zuv
                tzrho(i,j) = zrho_atm
             ENDIF
           ENDDO
         ENDDO
      ENDIF

!!DL IF ( ntstep == ntstep_max-1  ) THEN
!      PRINT *, 'Diagnostic, ntstep = ', ntstep
!      CALL stats('qv      ',qv(:,:,ke,nnew),llandmask,UNIT=6)
!     CALL stats('qv      ',qv(:,:,ke),     llandmask,UNIT=6)    ! XYZ> in TSA qv is 4 dimensional
!      CALL stats('qv_s    ',qv_s(:,:,nnew), llandmask,UNIT=6)
!      CALL stats('tch     ',tch(:,:),       llandmask,UNIT=6)
!      CALL stats('zuv     ',tzuv(:,:),      llandmask,UNIT=6)
!      CALL stats('zrho_atm',tzrho(:,:),     llandmask,UNIT=6)
!      CALL stats('lhfl_s  ',lhfl_s(:,:),    llandmask,UNIT=6)
!      CALL stats('shfl_s  ',shfl_s(:,:),    llandmask,UNIT=6)
!      CALL stats('lhfl_bs ',lhfl_bs(:,:),   llandmask,UNIT=6)
!DL some additional print outs - test
!      CALL stats('t       ',t(:,:,ke,nnew), llandmask,UNIT=6)
!      CALL stats('ps      ',ps(:,:,nnew)  , llandmask,UNIT=6)
!      CALL stats('p0      ',p0(:,:,ke)    , llandmask,UNIT=6)
!      CALL stats('t_g     ',t_g(:,:,nnew) , llandmask,UNIT=6)
       

!!DL ENDIF
    
    ENDIF ! lcalc

    ! Sum up all fields that are exported as time averaged values
    astore(:,:,1)=astore(:,:,1)+ lhfl_s
    astore(:,:,2)=astore(:,:,2)+ shfl_s
    astore(:,:,3)=astore(:,:,3)+ prr_gsp 
    astore(:,:,4)=astore(:,:,4)+ prs_gsp
    astore(:,:,5)=astore(:,:,5)+ sobs
    astore(:,:,6)=astore(:,:,6)+ thbs
    astore(:,:,7)=astore(:,:,7)+ runoff_s
    astore(:,:,8)=astore(:,:,8)+ runoff_g

    ! Reset runoff
    runoff_s=0.0
    runoff_g=0.0
   
    ! Output
    IF (MODULO(ntstep+1,nout_interval)==0) THEN
      IF (ntype_output<=3) THEN
!DL     GRIB output routine depending on yform_write
        IF (yform_write == 'grb1') THEN
          CALL grbout(ntstep+1)
        ELSE
          CALL grbout_eccodes(ntstep+1)
        ENDIF 
      ELSEIF (ntype_output==4) THEN
         CALL binout(ntstep+1)
      ELSEIF (ntype_output==5) THEN
         CALL ascout(ntstep+1)
      ELSE
         WRITE(6,*) "Unknown output type!"
         STOP
      ENDIF
      runoff=0.0
   ENDIF

    !--------------------------------------------------------
    ! Section 2.6: Finalize time loop
    !--------------------------------------------------------
    nnow=nnew
    nnew=nnew+1
    IF (nnew>2) THEN
      nnew=1
    ENDIF

  ENDDO ! End of time steping loop

!========================================================
! Section 3: Clean up
!========================================================

!!CALL WRITE_ctl()
!!IF (lconstout) CALL WRITE_const_ctl(nconstout)

  CALL clean_up
!CALL tsa_sfc_finalize

WRITE(6,*) "Program finished!"

!########################################################################################
!INCLUDE 'terra_lmenv.f90'    ! subroutines to establish the fields and variables 
!                             ! needed by the program

END PROGRAM tsa_main

