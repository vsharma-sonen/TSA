!==============================================================================
! Subroutines for Input of Terra Stand Alone (TSA)
!==============================================================================

MODULE  tsa_input

!------------------------------------------------------------------------------
!
! Description:
!  This file contains subroutines for input operation of the stand alone 
!  version of the soil module TERRA.
!
!  Currently included: Subroutines for input:
!
!    - read_metforc
!      reads the meteorological forcing data for TSA
!
!    - read_statbin
!      reads data from bin files for TSA
!
!    - read_const_fields
!      reads constant fields and assigns them to TSA
!
!    - read_const_fields_icon (NEW)
!      reads constant icon fields and assigns them to TSA
!
!    - read_initial_fields
!      reads initial conditions and assigns them to TSA
!
!    - read_initial_fields_icon (NEW)
!      reads initial icon fields and assigns them to TSA
!
!    - read_lmgrib
!      reads the local model grib files and assigns fields to TSA
!
!    - read_icongrib (NEW)
!      reads the icon LAM (?) grib files and assigns fields to TSA
!
!    - read_radolan
!      reads RADOLAN files for TSA
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
! V5.01      2015-12-01 Yiftach Ziv & Pavel Khain, IMS 
!  Added itype_hydcond, crsmin, kexpdec to namelist to better control
!    soil parameters, following Jurgen Helmert of DWD
!  Added functionality of sub-regions to improve computing
!    efficiency following Julian Todter(JT) of GUF  (PK)
!  Updated precision to working precision (wp)  (XYZ)
!  Arranged code to adhere to coding standards. (XYZ)
! V5.07      2020-07-15 Doerte Liermann, Ulrich Schaettler
!  Added eccodes to read/write grib files
!  General refactoring: this is now only a module for input
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!======================================================

! USE declarations
USE kind_parameters,       ONLY:                                            &
    wp           ! KIND-type parameter for real variables

USE tsa_data,              ONLY:                                            &
    czhls_const, qv, qv_bd, prr_gsp_bd, prs_gsp_bd,                         &
    so_down_bd, th_down_bd, pabs_bd, ps_bd, uuid,                           &
    w_g0, w_i0, t_soil0, w_snow0, t_cl0, w_cl0,                             &
    rlat_geo, rlon_geo, ydate_ini,                                          &
    radofiledir, metfiledir, metfileprefix,                                 &
    ntype_atminput, ntype_radinput, ntype_raininput, rundef,                &
    soilinitdir, soilinitprefix, soiltyp_const,                             &
    ke_model, ngp, ngpi, ngrid, ymodel, yform_write, yform_read,            &
    ntstep_max, tincr_max, constfilename,                                   &
    lrel_in        , lvol_in        , lvegadapt      , lhomosoil      ,     &
    lhomoinit      , lcheck         , lz0local       , lconstvegalb   ,     &
    lhourly_data   , ldestaggeruv   ,                                       &
    lpar           , lgettcl        , lcalc          ,                      &
    lai_const      , lai_min_const  , lai_max_const  , plcov_const    ,     &
    plcov_min_const, plcov_max_const, rootdp_const   , plcov_mn       ,     &
    plcov_mx       , rootdp0        , lai_mn         , lai_mx         ,     &
    z0             , rstom_mx       , rstom_mn       , vegalb         ,     &
    z0_const       , rstom_mn_const , rstom_mx_const , vegalb_const   ,     &
    dz             , dz_u           , nrecmax        , hsurface

USE data_io,               ONLY:                                            &
    intgribf, intgribc, irealgrib, iblock, lfd, lds, iwlength,              &
    iwlength

USE data_fields,           ONLY:                                            &
    u, t, ps, p0, prr_gsp, prs_gsp, u_10m, v_10m, t_2m, t_s, t_so, w_so,    &
    w_snow, rho_snow, freshsnow, t_snow, t_g, w_i, qv_s, tch, tcm,          &
    pabs, sobs, thbs, llandmask, lai, plcov, rootdp, soiltyp, hsurf,        &
    fr_land, rlon, rlat, t_cl, w_cl, w_so_ice, rsmin2d,                     &
    u_bd, t_bd

USE data_modelconfig,      ONLY:                                            &
    pollon, pollat, polgam, dlon, dlat, startlon, startlat, ke_soil, czmls, &
    ie, je, ke, dt, ke_snow
    
USE data_runcontrol,       ONLY:                                            &
    nnew, nnow, nblock, nproma, nlastproma, lstomata, itype_canopy,         &
    itype_evsl, itype_heatcond, itype_hydbound, itype_mire, itype_root,     &
    lmelt, lmelt_var, lmulti_snow, ntstep, itype_calendar

USE data_constants,        ONLY:                                            &
    sigma, b1, b2w, b3, b4w, r_d, g, o_m_rdv, rdv

USE sfc_terra_data,        ONLY:                                            &
    cf_snow, csalb_p, ctalb, cporv, crsmax, crsmin, idiag_snowfrac

USE meteo_utilities,       ONLY:                                            &
    tgcom

USE utilities,             ONLY:                                            &
    phirot2phi, phi2phirot, rla2rlarot, rlarot2rla, get_utc_date

USE tsa_gribio,            ONLY:                                            &
    read_grib_eccodes, read_grib, init_griblist, clear_griblist,            &
    add2griblist, get_gribinfo, get_grib_info, get_grib_info_icon

USE tsa_interpol,          ONLY:                                            &
    general_interpol2d

USE tsa_lmparam,           ONLY:                                            &
    qvsat, calc_albedo, near_surface, stats

USE support_datetime,      ONLY:                                            &
    date_delta, date_sum

USE fxtr_definition,       ONLY:                                            &
    kind_idate   , & !
    pc_std_p0        ! 101325 - standard slp for terra_io

!==============================================================================

IMPLICIT NONE

!==============================================================================

CONTAINS

!======================================================================
! I) subroutines for input: 
!======================================================================

!----------------------------------------------------------------------
SUBROUTINE read_metforc(n,doy,month)
!----------------------------------------------------------------------
! Description:
!
! Main routine to control input of atmospheric, radiation and
! precipitation data from external files. Following steps
! are executed:
! 1) generate date string
! 2) read initial data if first timestep
! 3) read data at next input interval
! 4) temporal interpolation for current simulation time
! 5) compute diagnostics (p0, T_2m, U_10m, net radiation)
!----------------------------------------------------------------------

!>XYZ commented - no tracers
!USE src_tracer, ONLY : trcr_get, trcr_errorstr
!<XYZ
  IMPLICIT NONE

! Subroutine arguments:
  INTEGER, INTENT(IN) :: n          ! number of actual time step
  INTEGER, INTENT(OUT):: doy,month  ! day of year, month
! Local variables
  LOGICAL                              :: &
    lnew_lmana,                           &
    lnew_rado,                            &
    lnew_rad                              
  LOGICAL                              :: &
    allow_bc_gap,                         &
    lt2m,                                 &
    lu10m
  INTEGER                              :: &
    idummy,                               &
    ntemp,                                &
    nave_new,                             &
    nave_new2,                            &
    nave_old2,                            &
    i, j, iend
  INTEGER(KIND=kind_idate)             :: actual_date
  REAL(KIND=wp)                        :: &
    acthour,                              &
    rdummy,                               &
    fac,                                  &
    fac1
  REAL(KIND=wp)                        :: &
    also,                                 &
    alth
  CHARACTER (LEN=14)                   :: &
    yactdate1,                            &
    yactdate1_next,                       &
    yeffdate
  CHARACTER (LEN=28)                   :: yactdate2
  CHARACTER (LEN=2)                    :: min


  INTEGER,  SAVE :: n1, n2, nrain1, nrain2
  INTEGER,  SAVE :: nave_old
  INTEGER,  SAVE :: delta_bc
  INTEGER,  SAVE :: ncount_rain=0, nrec=0
  INTEGER (KIND=kind_idate), SAVE :: oldbc_date, newbc_date


  INTEGER             :: izerror                        ! error STATUS variable

  CHARACTER (LEN=80)  :: yerrmsg

!XYZ> no use of tracers in TSA, qv is 4 dimentional
!!! Tracer pointers
!!!----------------
!!!  REAL (KIND=wp),     POINTER :: &
!!!    qv  (:,:,:) => NULL(),       &                    ! QV at tlev=nx
!!!    qv_bd  (:,:,:) => NULL()



!!!Retrieve the required microphysics tracers
! izerror = 0
! yerrmsg = '       '
! CALL trcr_get(izerror, idt_qv, ptr_tlev = n, ptr = qv)
! IF (izerror /= 0) THEN
!   yerrmsg = trcr_errorstr(izerror)
!   PRINT *, 'trcr_get ERROR: ', izerror, yerrmsg
! ENDIF
!XYZ<



  !=========================================================
  ! Section 1: Initialization
  !=========================================================

  lnew_lmana=.FALSE.       !
  lnew_rado =.FALSE.       !
  lnew_rad=.FALSE.         !

  allow_bc_gap=.FALSE.
  IF (ntype_atminput<=2) THEN
     allow_bc_gap = .TRUE.
  ENDIF

  CALL get_utc_date (n, ydate_ini, dt, itype_calendar,     &
                     yactdate1, yactdate2, doy,acthour) 
  READ(yactdate1(5:6),*) month
  !US with the new get_utc_date from COSMO, yactdate1 has length 14
  !   therefore read only the first 10 digits
  READ(yactdate1(1:10),*) actual_date

  !=========================================================
  ! Section 2: First Subroutine Call
  !=========================================================

  nnew = 1
  nnow = 1
  IF ( n == 0 ) THEN

     ! Initial setting of time levels
     n1=1
     n2=1
     ! Initial setting of previous and next boundary conditions
     oldbc_date = actual_date
     newbc_date = actual_date

     ! Reading of LM-Grib Files or ICON grib fields according to ymodel
     IF ((ntype_atminput<=3).OR.(ntype_raininput==1)) THEN
        IF (ymodel == 'COSMO') THEN
          WRITE(6,*) 'Reading LM-Gribfile at ', yactdate1(1:10)
          CALL read_lmgrib(yactdate1(1:10), n2, yeffdate, nave_old, look_forward=allow_bc_gap)
        ELSE
          WRITE(6,*) 'Reading ICON-Gribfile at ', yactdate1(1:10)
          CALL read_icongrib(yactdate1(1:10), n2, yeffdate, nave_old, look_forward=allow_bc_gap)
        ENDIF
        IF ( yactdate1(1:10) /= yeffdate ) THEN
           WRITE(6,*) 'Gap detected; using LM-Gribfile at ', yeffdate
        ENDIF
     ENDIF

     ! Reading of binary files for homogenous meteorological forcing
     IF (ntype_atminput==4) THEN
        CALL read_statascii_smet(nrec,1)
        nrec=nrec+1
        lnew_rado=.TRUE.
     ENDIF

     ! Reading rain data from RADOLAN files
     IF (ntype_raininput==2) THEN
        nrain1=1
        nrain2=2
        WRITE(6,*)'Reading Radolanfile at ',yactdate1(1:10)//'50'
        CALL read_radolan(yactdate1(1:10)//'50',nrain2)
        lnew_rado=.TRUE.
     ENDIF

  ENDIF

  !=========================================================
  ! Section 3: Reading new data if necessary
  !=========================================================

  !IF ( actual_date >= newbc_date ) THEN             ! next bc required?

     ! next hour
     CALL get_utc_date (1, yactdate1, 31.0_wp, itype_calendar,    & 
          yactdate1_next, yactdate2, idummy, rdummy) 

     ! exchanging the time levels
    ! ntemp=n2
    ! n2=n1
    ! n1=ntemp

     ! Reading of LM-Grib Files
     IF ((ntype_atminput<=3).OR.(ntype_raininput==1)) THEN
        IF ( ymodel == 'COSMO' ) THEN
          WRITE(6,*) 'Reading LM-Gribfile at ', yactdate1_next(1:10)
          CALL read_lmgrib(yactdate1_next(1:10), n2, yeffdate, nave_new, look_forward=allow_bc_gap)
        ELSE
          WRITE(6,*) 'Reading ICON-Gribfile at ', yactdate1_next(1:10)
          CALL read_icongrib(yactdate1_next(1:10), n2, yeffdate, nave_new, look_forward=allow_bc_gap)
        ENDIF
        lnew_lmana=.TRUE.
        IF ( yactdate1_next(1:10) /= yeffdate ) THEN
           yactdate1_next(1:10) = yeffdate 
           WRITE(6,*) 'Gap detected; using LM/ICON-Gribfile at ', yactdate1_next(1:10)
        ENDIF
     ENDIF

     ! Reading of binary files for homogenous meteorological forcing
     IF (ntype_atminput==4) THEN
        CALL read_statascii_smet(nrec,1)
        nrec=nrec+1
        IF ((nrecmax>0).AND.(nrecmax<=nrec)) THEN
           nrec=nrec-nrecmax
        ENDIF
        lnew_rado=.TRUE.
        lnew_rad=.TRUE.
        nrain2=1
     ENDIF

     oldbc_date = newbc_date
     READ(yactdate1_next(1:10),*) newbc_date
     delta_bc = date_delta(newbc_date*100,oldbc_date*100)

  !ENDIF

  ! Reading next RADOLAN FILE
  ! Radolan data is valid at xx.50UTC
  IF (ntype_raininput==2) THEN
     fac= (acthour - INT(acthour))*3600.
     IF ((fac>=3000).AND.(fac-dt<3000)) THEN  
        CALL get_utc_date (1,yactdate1,3601.0_wp, itype_calendar, & 
             yactdate1_next,yactdate2,idummy,rdummy) 
        ntemp=nrain2
        nrain2=nrain1
        nrain1=ntemp
        WRITE(6,*) 'Reading RADOLAN FILE at ',yactdate1_next(1:10)//'50'
        CALL read_radolan(yactdate1_next(1:10)//'50',nrain2)
        lnew_rado=.TRUE.
     ENDIF
  ENDIF



  !=========================================================
  ! Section 5: Temporal Interpolation
  !=========================================================
!  fac = (date_delta(actual_date*100,oldbc_date*100)    &
!        + (acthour-INT(acthour))*60.0_wp) / REAL(delta_bc)  
  fac1=1.0 !-fac

  ! Atmospheric variables at the reference level ! VARUN 
  IF (ntype_atminput<=4) THEN                   
     u(:,:,ke,1) = fac1*u_bd(:,:,ke,1) !+ fac*u_bd(:,:,ke,n2)
     t(:,:,ke,1) = fac1*t_bd(:,:,ke,1) !+ fac*t_bd(:,:,ke,n2)
     qv(:,:,ke,1)= fac1*qv_bd(:,:,ke,1) !+ fac*qv_bd(:,:,ke,n2)
!    qv(:,:,ke)= fac1*qv_bd(:,:,ke) + fac*qv_bd(:,:,ke)              !XYZ> uncommented qv non-tracer is 4 dimentional
     ps(:,:,1)  = fac1*ps_bd(:,:,1) !+ fac*ps_bd(:,:,n2)
  ELSE
     WRITE(6,*) "Invalid ntype_atminput! (Time interpolation)"
     STOP
  ENDIF

  ! Rain data
  IF (ntype_raininput==1) THEN     
     IF (lnew_lmana) THEN         
        IF (lhourly_data) THEN
           prr_gsp(:,:)= prr_gsp_bd(:,:,n2)
           prs_gsp(:,:)= prs_gsp_bd(:,:,n2)
        ELSE IF (nave_new>nave_old) THEN
           prr_gsp(:,:)= (prr_gsp_bd(:,:,n2)-prr_gsp_bd(:,:,n1)) &
           / FLOAT(nave_new-nave_old)
           prs_gsp(:,:)= (prs_gsp_bd(:,:,n2)-prs_gsp_bd(:,:,n1)) &
           / FLOAT(nave_new-nave_old)
        ELSE
           prr_gsp(:,:)= prr_gsp_bd(:,:,n2)/FLOAT(nave_new)
           prs_gsp(:,:)= prs_gsp_bd(:,:,n2)/FLOAT(nave_new)
        ENDIF
     ENDIF
  ELSEIF ((ntype_raininput==2).OR.(ntype_raininput==4)) THEN
     IF (lnew_rado) THEN
        WHERE (t(:,:,ke,nnew)>273.15) 
           prr_gsp(:,:)=prr_gsp_bd(:,:,1)
           prs_gsp(:,:)=0.0_wp
        ELSEWHERE
           prr_gsp(:,:)=0.0_wp
           prs_gsp(:,:)=prr_gsp_bd(:,:,1)
        END WHERE
     ENDIF
  ELSE
     WRITE(6,*) "Invalid ntype_raininput! (Time interpolation)"
  ENDIF

  ! Radiation
  IF (ntype_radinput<=2) THEN
     IF (lnew_lmana) THEN
        IF (lhourly_data) THEN
           so_down_bd(:,:,n1) = so_down_bd(:,:,n2)
           th_down_bd(:,:,n1) = th_down_bd(:,:,n2)
           pabs = pabs_bd(:,:,n2)
        ELSE IF (nave_new>nave_old) THEN
           so_down_bd(:,:,n1) = &
                (nave_new*so_down_bd(:,:,n2) - nave_old*so_down_bd(:,:,n1)) &
                / FLOAT(nave_new-nave_old)
           th_down_bd(:,:,n1) = &
                (nave_new*th_down_bd(:,:,n2) - nave_old*th_down_bd(:,:,n1)) &
                / FLOAT(nave_new-nave_old)
           pabs = (nave_new*pabs_bd(:,:,n2) - nave_old*pabs_bd(:,:,n1))     &
                  / FLOAT(nave_new-nave_old)
        ELSE
           so_down_bd(:,:,n1) = so_down_bd(:,:,n2)
           th_down_bd(:,:,n1) = th_down_bd(:,:,n2)
           pabs = pabs_bd(:,:,n2)
        ENDIF
     ENDIF
  ELSEIF (ntype_radinput==4) THEN
     IF (lnew_rad) THEN
        so_down_bd(:,:,1)=so_down_bd(:,:,1)
        th_down_bd(:,:,1)=th_down_bd(:,:,1)
     ENDIF
  ELSE
     WRITE(6,*) "Invalid ntype_radinput! (Time interpolation)"
     STOP
  ENDIF

  IF (lnew_lmana) THEN
     nave_old=nave_new
  ENDIF

  IF (lcalc) THEN

     !=========================================================
     ! Section 6: Conversions
     !=========================================================

     ! Copy time levels
     !u(:,:,ke,nnow)  = u(:,:,ke,nnew)
     !t(:,:,ke,nnow)  = t(:,:,ke,nnew) 
     !qv(:,:,ke,nnow) = qv(:,:,ke,nnew)
     !ps(:,:,nnow)    = ps(:,:,nnew)
     p0(:,:,ke)      = ps(:,:,1) * EXP( -g*dz/(r_d*t(:,:,ke,1)) )

     ! diagnosis of net radiation from down-welling fluxes
     IF (ymodel == 'COSMO') THEN
       DO i=1,ie
          DO j=1,je
             IF (llandmask(i,j)) THEN
                CALL calc_albedo(i,j,1,also,alth)    
                ! write(*,*) 'albedo: ',ntstep, also 
                sobs(i,j)= so_down_bd(i,j,1)*(1.0_wp-also)
                thbs(i,j)= th_down_bd(i,j,1)-(1.0_wp-alth)*sigma*t_g(i,j,1)**4
             ENDIF
          ENDDO
       ENDDO
     ELSE
!DL  ICON
       iend = nproma
       DO j=1,nblock
          IF(j==nblock) iend = nlastproma
          DO i=1,iend
             IF (llandmask(i,j)) THEN
                CALL calc_albedo(i,j,nnow,also,alth)     
                write(*,*) 'albedo: ',ntstep, also 
                sobs(i,j)= so_down_bd(i,j,n1)*(1.0_wp-also)
                thbs(i,j)= th_down_bd(i,j,n1)-(1.0_wp-alth)*sigma*t_g(i,j,nnow)**4
             ENDIF
          ENDDO
       ENDDO
     ENDIF

     IF (ntype_radinput==4) THEN
        pabs=0.5_wp*sobs
     ENDIF
        
     ! diagnosis of T_2m and u_10m
     IF ( ABS(dz_u-10.0_wp) < 1.0E-6_wp ) THEN
        u_10m(:,:) = u(:,:,ke,nnew)
        lu10m=.FALSE.
     ELSE IF ( dz_u>10.0_wp) THEN
        lu10m=.TRUE.
     ELSE
        WRITE(6,*) "dz_u lower than 10m is not allowed!"
     ENDIF

     IF ( ABS(dz-2.0_wp) < 1.0E-6_wp ) THEN
        t_2m(:,:)  = t(:,:,ke,nnew) 
        lt2m=.FALSE.
     ELSE IF ( dz > 2.0_wp ) THEN
        lt2m=.TRUE.
     ELSE
        WRITE(6,*) "dz lower than 2m is not allowed!"
     ENDIF

     IF ( lt2m .OR. lu10m) THEN
       CALL near_surface(nnew,lt2m,lu10m)
     ENDIF

     ! check
     IF (lcheck) THEN
        OPEN(UNIT=444,FILE='YUDRDAT',FORM='FORMATTED',STATUS='UNKNOWN',POSITION='APPEND')
        WRITE(444,*)
        WRITE(444,*) 'Meteorological forcing at nt=',ntstep
        CALL stats("u",u(:,:,ke,nnow),llandmask,UNIT=444)
        CALL stats("t",t(:,:,ke,nnow),llandmask,UNIT=444)
        CALL stats("qv",qv(:,:,ke,nnow),llandmask,UNIT=444)
!       CALL stats("qv",qv(:,:,ke),llandmask,UNIT=444)         !XYZ> commented
        CALL stats("ps",ps(:,:,nnow),llandmask,UNIT=444)
        CALL stats("sobs",sobs,llandmask,UNIT=444)
        CALL stats("thbs",thbs,llandmask,UNIT=444)
        CALL stats("pabs",pabs,llandmask,UNIT=444)
        CALL stats("prr_gsp",prr_gsp,llandmask,UNIT=444)
        CALL stats("prs_gsp",prs_gsp,llandmask,UNIT=444)
        CALL stats("t_2m",t_2m,llandmask,UNIT=444)
        CALL stats("u_10m",u_10m,llandmask,UNIT=444)
        CLOSE(UNIT=444)
     ENDIF

  ENDIF   !lcalc

END SUBROUTINE read_metforc
!----------------------------------------------------------------------



SUBROUTINE read_statbin(nr,nx)
!----------------------------------------------------------------------
! Description:
!
! Read binary file and assign fields to corresponding TERRA fields 
!----------------------------------------------------------------------

 !USE fxtr_definition, ONLY : pc_std_p0

  IMPLICIT NONE

! Subroutine arguments:
  INTEGER, INTENT(IN) :: &
       nr, &               ! record number
       nx                  ! time level
! Local variables
  INTEGER, PARAMETER :: nvarmax=6
  REAL (KIND=wp) :: rvec(nvarmax)
!  REAL (KIND=4) :: rvec(nvarmax)

!Declaration of STATEMENT-FUNCTIONS
!--------------------
  REAL (KIND=wp) :: fpvsw, zst, fqvs, zsge, zsp
  fpvsw(zst)       = b1 * EXP( b2w * (zst - b3) / (zst - b4w) )
  fqvs (zsge, zsp) = rdv * zsge / ( zsp - o_m_rdv * zsge )

  ! READ one record to rvec
  OPEN (90,FILE=TRIM(metfiledir)//'/'//TRIM(metfileprefix),ACCESS='direct',RECL=nvarmax*4)
  READ (90,REC=nr+1) rvec
  CLOSE (90)

  ! assign rvec to corresponding TERRA fields 
  u_bd(:,:,ke,nx)=rvec(3)
  t_bd(:,:,ke,nx)=rvec(1)+273.15
  qv_bd(:,:,ke,nx)= MIN( 1.0_wp , rvec(2)/100.0_wp )                    &                    
        * fqvs ( fpvsw ( rvec(1)+273.15_wp ), REAL(pc_std_p0,wp) )
! qv_bd(:,:,ke)= min( 1.0 , rvec(2)/100.0 ) * &                        !XYZ> commented
!      fqvs ( fpvsw ( REAL(rvec(1)+273.15,wp) ), REAL(pc_std_p0,wp) )
  ps_bd(:,:,nx)=REAL(pc_std_p0,wp)
  IF (rvec(1)>0) THEN
     prr_gsp_bd(:,:,nx)=rvec(4)/3600.0_wp
     prs_gsp_bd(:,:,nx)=0.0
  ELSE
     prr_gsp_bd(:,:,nx)=0.0
     prs_gsp_bd(:,:,nx)=rvec(4)/3600.0_wp
  ENDIF
  so_down_bd(:,:,nx)=rvec(5)
  th_down_bd(:,:,nx)=rvec(6)
  pabs_bd(:,:,nx)=0.0_wp   

END SUBROUTINE read_statbin
!----------------------------------------------------------------------

SUBROUTINE read_statascii(nr,nx)
!----------------------------------------------------------------------
! Description:
!
! Read binary file and assign fields to corresponding TERRA fields 
!----------------------------------------------------------------------

 !USE fxtr_definition, ONLY : pc_std_p0

  IMPLICIT NONE

! Subroutine arguments:
  INTEGER, INTENT(IN) :: &
       nr, &               ! record number
       nx                  ! time level
! Local variables
  INTEGER, PARAMETER :: nvarmax=6
  REAL (KIND=wp) :: rvec(nvarmax)
!  REAL (KIND=4) :: rvec(nvarmax)
  real (kind=wp) :: a(15)
  integer, parameter :: out_unit = 90
!Declaration of STATEMENT-FUNCTIONS
!--------------------
  REAL (KIND=wp) :: fpvsw, zst, fqvs, zsge, zsp
  fpvsw(zst)       = b1 * EXP( b2w * (zst - b3) / (zst - b4w) ) ! saturation vapor pressure
  fqvs (zsge, zsp) = rdv * zsge / ( zsp - o_m_rdv * zsge ) ! saturation mixing ratio

  ! READ one record to rvec
!  OPEN (90,FILE=TRIM(metfiledir)//'/'//TRIM(metfileprefix),ACCESS='direct',RECL=nvarmax*4)
!  READ (90,REC=nr+1) rvec
!  CLOSE (90)

  if(nr == 0) then
    OPEN (out_unit,FILE='./input/met_data.txt')
  endif
    read (out_unit,*) a(:)

!  write(*,*) ' IN STATS_ASCII: ',nr,nx
!  write(*,*) a 
  ! RVEC(1) = TA
  ! RVEC(2) = relative humidity
  ! RVEC(3) = wind speed
  ! RVEC(4) = precip in mm ( or mm/h ?? )
  ! RVEC(5) = Downwelling SHORTWAVE
  ! RVEC(6) = Downwelling LONGWAVE
    
  ! assign rvec to corresponding TERRA fields 
  u_bd(:,:,ke,nx)= a(6)
  t_bd(:,:,ke,nx)= a(1)
  qv_bd(:,:,ke,nx)= MIN( 1.0_wp , a(2) )                    &                    
        * fqvs ( fpvsw ( a(1) ), REAL(pc_std_p0,wp) )
! qv_bd(:,:,ke)= min( 1.0 , rvec(2)/100.0 ) * &                        !XYZ> commented
!      fqvs ( fpvsw ( REAL(rvec(1)+273.15,wp) ), REAL(pc_std_p0,wp) )
  
  ps_bd(:,:,nx)=101000.0_wp !REAL(pc_std_p0,wp)
  !IF (rvec(1)>0) THEN
     prr_gsp_bd(:,:,nx) =  a(11)/900.0_wp !rvec(4)/3600.0_wp
     prs_gsp_bd(:,:,nx) =  0.0
  !ELSE
  !   prr_gsp_bd(:,:,nx)=0.0
  !   prs_gsp_bd(:,:,nx)=0.001_wp !rvec(4)/3600.0_wp
  !ENDIF
  so_down_bd(:,:,nx) = a(9)
  th_down_bd(:,:,nx) = a(10)
  pabs_bd(:,:,nx)=0.0_wp   

END SUBROUTINE read_statascii
!----------------------------------------------------------------------

SUBROUTINE read_statascii_smet(nr,nx)
!----------------------------------------------------------------------
! Description:
!
! Read binary file and assign fields to corresponding TERRA fields 
!----------------------------------------------------------------------

 !USE fxtr_definition, ONLY : pc_std_p0

  IMPLICIT NONE

! Subroutine arguments:
  INTEGER, INTENT(IN) :: &
       nr, &               ! record number
       nx                  ! time level
! Local variables
  INTEGER, PARAMETER :: nvarmax=6
  REAL (KIND=wp) :: rvec(nvarmax)
!  REAL (KIND=4) :: rvec(nvarmax)
  real (kind=wp) :: a(7)
  integer, parameter :: out_unit = 90
  integer :: i
!Declaration of STATEMENT-FUNCTIONS
!--------------------
  REAL (KIND=wp) :: fpvsw, zst, fqvs, zsge, zsp
  fpvsw(zst)       = b1 * EXP( b2w * (zst - b3) / (zst - b4w) ) ! saturation vapor pressure
  fqvs (zsge, zsp) = rdv * zsge / ( zsp - o_m_rdv * zsge ) ! saturation mixing ratio

  ! READ one record to rvec
!  OPEN (90,FILE=TRIM(metfiledir)//'/'//TRIM(metfileprefix),ACCESS='direct',RECL=nvarmax*4)
!  READ (90,REC=nr+1) rvec
!  CLOSE (90)

  if(nr == 0) then
    OPEN (out_unit,FILE='./input/input_meteo.txt')
    do i = 1, 16,1
       read(out_unit,*)
    enddo
  endif
    read (out_unit,*) a(:)

!  write(*,*) ' IN STATS_ASCII: ',nr,nx
!  write(*,*) a 
  ! RVEC(1) = TA
  ! RVEC(2) = relative humidity
  ! RVEC(3) = wind speed
  ! RVEC(4) = precip in mm ( or mm/h ?? )
  ! RVEC(5) = Downwelling SHORTWAVE
  ! RVEC(6) = Downwelling LONGWAVE
    
  ! assign rvec to corresponding TERRA fields 
  u_bd(:,:,ke,nx)= a(3)
  t_bd(:,:,ke,nx)= a(1)
  qv_bd(:,:,ke,nx)= (MIN( 100.0_wp , a(2) ) / 100.0_wp)                    &                    
        * fqvs ( fpvsw ( a(1) ), REAL(pc_std_p0,wp) )
! qv_bd(:,:,ke)= min( 1.0 , rvec(2)/100.0 ) * &                        !XYZ> commented
!      fqvs ( fpvsw ( REAL(rvec(1)+273.15,wp) ), REAL(pc_std_p0,wp) )
  
  ps_bd(:,:,nx)=101000.0_wp !REAL(pc_std_p0,wp)
  !IF (rvec(1)>0) THEN
     prr_gsp_bd(:,:,nx) =  a(6)/30.0_wp !rvec(4)/3600.0_wp
     prs_gsp_bd(:,:,nx) =  0.0
  !ELSE
  !   prr_gsp_bd(:,:,nx)=0.0
  !   prs_gsp_bd(:,:,nx)=0.001_wp !rvec(4)/3600.0_wp
  !ENDIF
  so_down_bd(:,:,nx) = a(4)
  th_down_bd(:,:,nx) = a(5)
  pabs_bd(:,:,nx)=0.0_wp   

END SUBROUTINE read_statascii_smet
!----------------------------------------------------------------------


SUBROUTINE read_const_fields(month)
!----------------------------------------------------------------------
! Description:
!
! Read assign fields to corresponding TERRA fields
!  from external files or from periodic vegetation assumption 
!----------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
  INTEGER, INTENT(IN), OPTIONAL :: &
     month

! Local variables

  LOGICAL                                        :: &
       lcrop                                      , & ! 
       is_found                                       ! 
  INTEGER                                        :: &
     ie_in                                        , & ! 
     je_in                                        , & ! 
     i                                            , & ! 
     j                                            , & ! 
     i0_crop                                      , & ! 
     j0_crop                                      , & ! 
     isoil                                        , & ! 
     ilev                                             ! 
!DL
  INTEGER                                        :: &
     imiss
  REAL (KIND=wp)                                 :: &
     r
  REAL (KIND=wp)                                 :: &
       rlat0_in                                   , & ! 
       rlon0_in                                   , & ! 
       dlat_in                                    , & ! 
       dlon_in                                    , & ! 
       pollat_in                                  , & ! 
       pollon_in                                      ! 
  REAL (KIND=wp), ALLOCATABLE                    :: &
       rlon_in(:,:)                               , & ! 
       rlat_in(:,:)                               , & ! 
       hsurf_in(:,:)                              , & ! 
       fr_land_in(:,:)                            , & ! 
       soiltyp_in(:,:)                            , & ! 
       z0_in(:,:)                                 , & ! 
       plcov_mn_in(:,:)                           , & ! 
       plcov_mx_in(:,:)                           , & ! 
       plcov_in(:,:)                              , & ! 
       lai_mn_in(:,:)                             , & ! 
       lai_mx_in(:,:)                             , & ! 
       lai_in(:,:)                                , & ! 
       vegalb_in(:,:)                             , & ! 
       rstom_mn_in(:,:)                           , & ! 
       rstom_mx_in(:,:)                           , & ! 
       rootdp_in(:,:)                             , & ! 
       t_cl_in(:,:)                                   ! 
  !    psand_in(:,:)                              , & ! 
  !    psand(:,:)                                 , & ! 
  !    pclay_in(:,:)                              , & ! 
  !    pclay(:,:)                                     ! 
  LOGICAL, ALLOCATABLE                           :: &
       mask1(:,:)                                     !
  CHARACTER (LEN=100)                            :: constfilename2

!!
  ! homogeneous soil/surface paramters
  IF (lhomosoil) THEN                   
     soiltyp  = soiltyp_const
     rootdp   = rootdp_const
     rootdp0  = rootdp_const
     plcov    = plcov_const    
     plcov_mn = plcov_min_const   
     plcov_mx = plcov_max_const           
     lai      = lai_const
     lai_mn   = lai_min_const
     lai_mx   = lai_max_const
     z0       = z0_const
!GDM>
     llandmask=.TRUE. ! GDM added for call to stats
     IF (.NOT. lstomata) THEN
!GDM<
        rstom_mn=crsmin
        rstom_mx=crsmax
     ELSE
        rstom_mn=rstom_mn_const
        rstom_mx=rstom_mx_const
     ENDIF
     IF (lconstvegalb) THEN
        vegalb=csalb_p
     ELSE
        vegalb=vegalb_const
     ENDIF

  ! in CASE soil/surface paramters provided by a grib-FILE
  ELSE                           

     IF (PRESENT(month)) THEN
        WRITE(constfilename2,'(A,A,I2.2,A)') TRIM(constfilename),'_2006',month,'0100'
     ELSE
        constfilename2=constfilename
     ENDIF


     ! determine size and position of input files
     IF (yform_read == 'apix') THEN
     CALL get_grib_info(constfilename2,                &
          ie_in,je_in,rlat0_in,rlon0_in,              &
          dlat_in,dlon_in,pollat_in,pollon_in,is_found)
     ELSE
     CALL get_gribinfo(constfilename2,                &
          ie_in,je_in,rlat0_in,rlon0_in,              &
          dlat_in,dlon_in,pollat_in,pollon_in,is_found)
     ENDIF

     IF ( .NOT. is_found ) THEN
        WRITE(6,*) "File ",TRIM(constfilename2)," does not exits!"," LOC A"
        STOP
     ENDIF

     ! allocate input fields
     ALLOCATE (                    &
          rlon_in(ie_in,je_in)   , &
          rlat_in(ie_in,je_in)   , &
          hsurf_in(ie_in,je_in)  , &  
          fr_land_in(ie_in,je_in), &
          soiltyp_in(ie_in,je_in), &
          z0_in(ie_in,je_in)     , &
          rootdp_in(ie_in,je_in) , &
          mask1(ie_in,je_in) )
     rlon_in=rundef
     rlat_in=rundef
     hsurf_in=rundef
     fr_land_in=rundef
     soiltyp_in=rundef
     z0_in=rundef
     rootdp_in=rundef
     IF (lvegadapt) THEN
        ALLOCATE (                     &
             plcov_mn_in(ie_in,je_in), &
             plcov_mx_in(ie_in,je_in), &
             lai_mn_in(ie_in,je_in)  , &
             lai_mx_in(ie_in,je_in) )
        plcov_mn_in=rundef
        plcov_mx_in=rundef
        lai_mn_in=rundef
        lai_mx_in=rundef
     ELSE
        ALLOCATE (                     &
             plcov_in(ie_in,je_in)   , &
             lai_in(ie_in,je_in) )
        plcov_in=rundef
        lai_in=rundef
     ENDIF
     IF (lstomata) THEN 
        ALLOCATE (                     &
             rstom_mn_in(ie_in,je_in), &
             rstom_mx_in(ie_in,je_in) )
        rstom_mn_in=rundef
        rstom_mx_in=rundef
     ENDIF
     IF (.NOT. lconstvegalb) THEN 
        ALLOCATE (vegalb_in(ie_in,je_in))
        vegalb_in=rundef
     ENDIF
     IF (lgettcl) THEN
        ALLOCATE (t_cl_in(ie_in,je_in))
        t_cl_in=rundef
     ENDIF
     !    IF (itype_hydparam==4) THEN
     !       ALLOCATE(pclay_in(ie_in,je_in), &
     !            pclay(ie,je), & 
     !            psand_in(ie_in,je_in), &
     !            psand(ie,je))
     !       pclay=rundef
     !       psand_in=rundef
     !       psand=rundef
     !    ENDIF

     IF (yform_read == 'apix') THEN
     ! READ data from gribfile
     CALL init_griblist(99,ie_in,je_in,1)
     CALL add2griblist(rlon_in,     'RLON','surface')
     CALL add2griblist(rlat_in,     'RLAT','surface')
     CALL add2griblist(hsurf_in,    'HSURF','surface')
     CALL add2griblist(fr_land_in,  'FR_LAND','surface')
     CALL add2griblist(soiltyp_in,  'SOILTYP','surface')
     IF (lz0local) THEN
        WRITE(6,*) "Z0LOC not implemented yet!"
!! DL   CALL add2griblist(z0_in,    82,250,1)   ! GDM> new numbers instead of 100/202 (SKYVIEW)
     ELSE
        CALL add2griblist(z0_in,    'Z0','surface')
     ENDIF
     IF (lvegadapt) THEN
        CALL add2griblist(plcov_mn_in, 'PLCOV_MN','surface')
        CALL add2griblist(plcov_mx_in, 'PLCOV_MX','surface')
        CALL add2griblist(lai_mn_in,   'LAI_MN','surface')
        CALL add2griblist(lai_mx_in,   'LAI_MX','surface')
     ELSE
        CALL add2griblist(plcov_in, 'PLCOV','surface')
        CALL add2griblist(lai_in,   'LAI','surface')
     ENDIF
     CALL add2griblist(rootdp_in,   'ROOTDP','surface')
     IF (lstomata) THEN
        WRITE(6,*) "lstomata=true not implementd yet!"
!! DL   CALL add2griblist(rstom_mn_in, 212,201,1) ! GDM> new numbers instead of 84/202 (AER_SO4)
! GDM>  CALL add2griblist(rstom_mx_in,  85,202,1)
     ENDIF
     IF (.NOT.lconstvegalb) THEN
        WRITE(6,*) "lconstvegalb=false not implementd yet!"
!! DL   CALL add2griblist(vegalb_in,   213,201,1) ! GDM> new numbers instead of 86/202 (AER_DUST)
     ENDIF
!    IF (lgettcl) THEN
!       IF (lmulti_layer) THEN
!          WRITE(6,*) "lgettcl / lmulti_layer =true not implemented yet!"
!          ilev=INT(czmls(ke_soil+1)*100+0.99)
!! GRIB1  CALL add2griblist(t_cl_in(:,:),'T_SO','depthBelowLand',0,ilev)
!! GRIB2  CALL add2griblist(t_cl_in(:,:),'T_2M_CL','depthBelowLand',0,ilev)
!! DL     CALL add2griblist(t_cl_in(:,:),197,201,111,0,ilev)
!       ELSE
!          CALL add2griblist(t_cl_in(:,:),'T_CL_LM','depthBelowLand',0,41)
!       ENDIF
!    ENDIF
     !    IF (itype_hydparam==4) THEN
     !        CALL add2griblist(psand_in,  98,202,1)
     !        CALL add2griblist(pclay_in,  97,202,1)
     !     ENDIF

     CALL read_grib_eccodes(constfilename2)
     CALL clear_griblist()

     ELSE

     ! READ data from gribfile
     CALL init_griblist(99,ie_in,je_in,1)  
     CALL add2griblist(rlon_in,    115,202,1)
     CALL add2griblist(rlat_in,    114,202,1)
     CALL add2griblist(hsurf_in,     8,  2,1)
     CALL add2griblist(fr_land_in,  81,  2,1)
     CALL add2griblist(soiltyp_in,  57,202,1)
     IF (lz0local) THEN
        CALL add2griblist(z0_in,    82,250,1)   ! GDM> new numbers instead of 100/202 (SKYVIEW)
     ELSE
        CALL add2griblist(z0_in,    83,  2,1)
     ENDIF
     IF (lvegadapt) THEN
        CALL add2griblist(plcov_mn_in, 68,202,1)
        CALL add2griblist(plcov_mx_in, 67,202,1)
        CALL add2griblist(lai_mn_in,   70,202,1)
        CALL add2griblist(lai_mx_in,   69,202,1)
     ELSE
        CALL add2griblist(plcov_in, 87,  2,1)
        CALL add2griblist(lai_in,   61,202,1)
     ENDIF
     CALL add2griblist(rootdp_in,   62,202,1)
     IF (lstomata) THEN 
        CALL add2griblist(rstom_mn_in, 212,201,1) ! GDM> new numbers instead of 84/202 (AER_SO4)
! GDM>  CALL add2griblist(rstom_mx_in,  85,202,1)
     ENDIF
     IF (.NOT.lconstvegalb) THEN 
        CALL add2griblist(vegalb_in,   213,201,1) ! GDM> new numbers instead of 86/202 (AER_DUST)
     ENDIF
!    IF (lgettcl) THEN
!       IF (lmulti_layer) THEN
!          ilev=INT(czmls(ke_soil+1)*100+0.99)
!          CALL add2griblist(t_cl_in(:,:),197,201,111,0,ilev)
!       ELSE
!          CALL add2griblist(t_cl_in(:,:),85,2,111,0,41)
!       ENDIF
!    ENDIF
     !    IF (itype_hydparam==4) THEN
     !        CALL add2griblist(psand_in,  98,202,1)
     !        CALL add2griblist(pclay_in,  97,202,1)
     !     ENDIF
 
     CALL read_grib(constfilename2)
     CALL clear_griblist()              
     ENDIF
!
!DL ##############################################################################
!   Check, if FR_LAND, SOILTYP, HSURF etc are available
    imiss = 0
    IF (ALL(fr_land_in==rundef)) THEN
       imiss = imiss + 1
       WRITE(6,*) "FR_LAND is missing!"
    ELSEIF (ALL(soiltyp_in==rundef)) THEN
       imiss = imiss + 1
       WRITE(6,*) "SOILTYP is missing!"
    ELSEIF (ALL(hsurf_in==rundef)) THEN
       imiss = imiss + 1
       WRITE(6,*) "HSURF is missing!"
    ELSEIF (ALL(z0_in==rundef)) THEN
       imiss = imiss + 1
       WRITE(6,*) "z0 is missing!"
    ELSEIF (ALL(rootdp_in==rundef)) THEN
       imiss = imiss + 1
       WRITE(6,*) "ROOTDP is missing!"
    ENDIF
    IF(lvegadapt) THEN
       IF (ALL(plcov_mn_in==rundef)) THEN
          imiss = imiss + 1
          WRITE(6,*) "PLCOV_MN is missing!"
       ELSEIF (ALL(plcov_mx_in==rundef)) THEN
          imiss = imiss + 1
          WRITE(6,*) "PLCOV_MX is missing!"
       ELSEIF (ALL(lai_mn_in==rundef)) THEN
          imiss = imiss + 1
          WRITE(6,*) "LAI_MN is missing!"
       ELSEIF (ALL(lai_mx_in==rundef)) THEN
          imiss = imiss + 1
          WRITE(6,*) "LAI_MX is missing!"
       ENDIF 
    ELSE
       IF (ALL(plcov_in==rundef)) THEN
          imiss = imiss + 1
          WRITE(6,*) "PLCOV is missing!"
       ELSEIF (ALL(lai_in==rundef)) THEN
          imiss = imiss + 1
          WRITE(6,*) "LAI is missing!"
       ENDIF
    ENDIF
    IF (imiss > 0) THEN
       WRITE(6,*) "!!! There are ",imiss," fields missing in read_const_fields!!!"
       STOP 'read_const_fields'
    ENDIF
!DL ##############################################################################

     ! YUCHKDAT
     IF (lcheck) THEN
        OPEN(UNIT=444,FILE='YUCHKDAT',FORM='FORMATTED',STATUS='UNKNOWN')
        WRITE(444,*)
        WRITE(444,*) 'Reading ',TRIM(constfilename2)
        CALL stats("rlon",rlon_in,UNIT=444)
        CALL stats("rlat",rlat_in,UNIT=444)
        CALL stats("hsurf",hsurf_in,UNIT=444)
        CALL stats("fr_land",fr_land_in,UNIT=444)
        CALL stats("soiltyp",soiltyp_in,UNIT=444)
        CALL stats("z0",z0_in,UNIT=444)
        IF (lvegadapt) THEN
           CALL stats("plcov_mn",plcov_mn_in,UNIT=444)
           CALL stats("plcov_mx",plcov_mx_in,UNIT=444)
           CALL stats("lai_mn",lai_mn_in,UNIT=444)
           CALL stats("lai_mx",lai_mx_in,UNIT=444)
        ELSE
           CALL stats("plcov",plcov_in,UNIT=444)
           CALL stats("lai",lai_in,UNIT=444)
        ENDIF
        CALL stats("rootdp",rootdp_in,UNIT=444)
        IF (lstomata) THEN  ! GDM CHANGED
           CALL stats("rstom_mn",rstom_mn_in,UNIT=444)
           CALL stats("rstom_mx",rstom_mx_in,UNIT=444)
        ENDIF
        IF (.NOT.lconstvegalb) THEN
           CALL stats("vegalb",vegalb_in,UNIT=444)
        ENDIF
        IF (lgettcl) THEN
           CALL stats("t_cl",t_cl_in,UNIT=444)
        ENDIF
        CLOSE(UNIT=444)
     ENDIF

     ! tell user what we are doing with external parameters
     lcrop=.FALSE.
     IF ( (ABS(dlon_in-dlon)<1.0E-3_wp)          &
          .AND. (ABS(dlat_in-dlat)<1.0E-3_wp) )  THEN
        i0_crop=NINT((startlon-rlon0_in)/dlon+1)
        j0_crop=NINT((startlat-rlat0_in)/dlat+1)
        IF ((i0_crop>=1)                                                &
            .AND. (j0_crop>=1)                                          &
            .AND. (ABS((i0_crop-1)*dlon-startlon+rlon0_in)<1.0E-3_wp)   &
            .AND. (ABS((j0_crop-1)*dlat-startlat+rlat0_in)<1.0E-3_wp) ) THEN
           lcrop=.TRUE.
        ENDIF
     ENDIF
     IF (lcrop) THEN
        IF (i0_crop==1 .AND. j0_crop==1 .AND. ie==ie_in .AND. je==je_in) THEN
           WRITE(6,*) "Perfect match of external parameters!"
        ELSE
           WRITE(6,*) "Croping external parameters!"
        ENDIF
     ELSE
        WRITE(6,*) "Interpolating external parameters!"
     ENDIF

     ! convert geographical coordinates into rotated LM coordinates
     IF (.NOT.lcrop) THEN
!DL Check, if rlon_in and rlat_in are available
        IF (ALL(rlon_in==rundef) .OR. ALL(rlat_in==rundef)) THEN
           WRITE(6,*) "rlon_in and/or rlat_in are missing! Stop in read_const_fields!"
           STOP
        ENDIF
        DO i=1,ie_in
           DO j=1,je_in
              r=phi2phirot(rlat_in(i,j),rlon_in(i,j),pollat,pollon)
              rlon_in(i,j)=rla2rlarot(rlat_in(i,j),rlon_in(i,j),pollat,pollon,0.0_wp)   ! XYZ> added polgam=0 for compatability with utilities
              rlat_in(i,j)=r
           ENDDO
        ENDDO
     ENDIF

     ! Spatial interpolation
     IF (lcrop) THEN
        fr_land=fr_land_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1) 
        hsurf=hsurf_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1) 
        soiltyp=soiltyp_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1) 
        z0=z0_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1) 
        IF (lvegadapt) THEN
           plcov_mn=plcov_mn_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1) 
           plcov_mx=plcov_mx_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1) 
           lai_mn=lai_mn_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1) 
           lai_mx=lai_mx_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1) 
        ELSE
           plcov=plcov_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1) 
           lai=lai_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1) 
        ENDIF
        rootdp=rootdp_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1) 
        llandmask=fr_land>=0.5     !SZB
        IF (lstomata) THEN 
           rstom_mn=rstom_mn_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)
! GDM      rstom_mx=rstom_mx_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)
           rsmin2d =rstom_mn  ! GDM CHANGED
           rstom_mx=rstom_mn  ! GDM CHANGED
        ELSE
           rstom_mn=crsmin
           rstom_mx=crsmax
        ENDIF
        IF (.NOT.lconstvegalb) THEN 
           vegalb=vegalb_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)
        ELSE
           vegalb=csalb_p
        ENDIF
        IF (lgettcl) THEN
           t_cl=t_cl_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)
        ENDIF
        !IF (itype_hydparam==4) THEN
        !   psand=psand_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)*100.0 
        !   pclay=pclay_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)*100.0 
        !ENDIF
     ELSE
        CALL general_interpol2d(rlat_in,rlon_in,fr_land_in,     &
             startlat,startlon,dlat,dlon,fr_land,               &
             2,2,0.51_wp,1.0_wp)
        mask1=(fr_land_in>=0.5_wp)      !SZB
        llandmask=(fr_land>=0.5_wp)     !SZB
        CALL general_interpol2d(rlat_in,rlon_in,hsurf_in,       &
             startlat,startlon,dlat,dlon,hsurf,                 &
             2,2,0.51_wp,1.0_wp)
        CALL general_interpol2d(rlat_in,rlon_in,soiltyp_in,     &
             startlat,startlon,dlat,dlon,soiltyp,               &
             2,4,0.51_wp,1.0_wp,mask1,llandmask)
        CALL general_interpol2d(rlat_in,rlon_in,z0_in,          &
             startlat,startlon,dlat,dlon,z0,                    &
             5,2,0.51_wp,1.0_wp)
        IF (lvegadapt) THEN
           CALL general_interpol2d(rlat_in,rlon_in,plcov_mn_in, &
                startlat,startlon,dlat,dlon,plcov_mn,           &
                2,2,0.51_wp,1.0_wp,mask1,llandmask)
           CALL general_interpol2d(rlat_in,rlon_in,plcov_mx_in, &
                startlat,startlon,dlat,dlon,plcov_mx,           &
                2,2,0.51_wp,1.0_wp,mask1,llandmask)
           CALL general_interpol2d(rlat_in,rlon_in,lai_mn_in,   &
                startlat,startlon,dlat,dlon,lai_mn,             &
                2,2,0.51_wp,1.0_wp,mask1,llandmask)
           CALL general_interpol2d(rlat_in,rlon_in,lai_mx_in,   &
                startlat,startlon,dlat,dlon,lai_mx,             &
                2,2,0.51_wp,1.0_wp,mask1,llandmask)
        ELSE
           CALL general_interpol2d(rlat_in,rlon_in,plcov_in,    &
                startlat,startlon,dlat,dlon,plcov,              &
                2,2,0.51_wp,1.0_wp,mask1,llandmask)
           CALL general_interpol2d(rlat_in,rlon_in,lai_in,      &
                startlat,startlon,dlat,dlon,lai,                &
                2,2,0.51_wp,1.0_wp,mask1,llandmask)
        ENDIF
        CALL general_interpol2d(rlat_in,rlon_in,rootdp_in,      &
             startlat,startlon,dlat,dlon,rootdp,                &
             2,2,0.51_wp,1.0_wp,mask1,llandmask)
        IF (lstomata) THEN 
           CALL general_interpol2d(rlat_in,rlon_in,rstom_mn_in, &
                startlat,startlon,dlat,dlon,rstom_mn,           &
                2,2,0.51_wp,1.0_wp,mask1,llandmask)
! GDM>     CALL general_interpol2d(rlat_in,rlon_in,rstom_mx_in, &
!               startlat,startlon,dlat,dlon,rstom_mx, &
!               2,2,0.51_wp,1.0_wp,mask1,llandmask)
           rsmin2d =rstom_mn  ! GDM CHANGED
           rstom_mx=rstom_mn  ! GDM CHANGED
        ELSE
           rstom_mn=crsmin
           rstom_mx=crsmax
        ENDIF
        IF (.NOT.lconstvegalb) THEN 
           CALL general_interpol2d(rlat_in,rlon_in,vegalb_in,   &
                startlat,startlon,dlat,dlon,vegalb,             &
                2,2,0.51_wp,1.0_wp,mask1,llandmask)
        ELSE
           vegalb=csalb_p
        ENDIF
        IF (lgettcl) THEN
           CALL general_interpol2d(rlat_in,rlon_in,t_cl_in,     &
                startlat,startlon,dlat,dlon,t_cl,               &
                2,2,0.51_wp,1.0_wp,mask1,llandmask)
        ENDIF
        !IF (itype_hydparam==4) THEN
        !   CALL general_interpol2d(rlat_in,rlon_in,psand_in, &
        !        startlat,startlon,dlat,dlon,psand, &
        !        2,2,0.51_wp,1.0_wp,mask1,llandmask)
        !   CALL general_interpol2d(rlat_in,rlon_in,pclay_in, &
        !        startlat,startlon,dlat,dlon,pclay, &
        !        2,2,0.51_wp,1.0_wp,mask1,llandmask)
        !   psand=psand*100.0
        !   pclay=pclay*100.0
        !ENDIF

        ! ensure consistency between fr_land (llandmask) and soiltyp
        WHERE (fr_land<0.5_wp)     !SZB
           soiltyp=9
        END WHERE
        WHERE (fr_land>=0.5_wp .AND. NINT(soiltyp)>=9)     !SZB
           fr_land=0.0_wp
        END WHERE
        ! final land-sea mask according to soil-type map;
        ! all point with llandmask=.FALSE. will not be considered for TERRA simulations 
        llandmask= fr_land>=0.5_wp   !SZB
     ENDIF

     ! set initial values for variables exhibiting an annual cycle
     IF (lvegadapt) THEN
        plcov=plcov_mx
        rootdp0=rootdp                                                     
        lai=lai_mx
     ENDIF
     ! Convert plcov from % to fraction, if needed (ecoclimap)
     IF (MAXVAL(plcov)>2) THEN
        plcov=plcov/100.0_wp
     ENDIF

     ! Input fields are no longer needed
     DEALLOCATE(rlon_in,rlat_in,hsurf_in,fr_land_in,soiltyp_in,z0_in,rootdp_in,mask1)
     IF (lvegadapt) THEN
        DEALLOCATE(plcov_mn_in,plcov_mx_in,lai_mn_in,lai_mx_in)
     ELSE
        DEALLOCATE(plcov_in,lai_in)
     ENDIF
     IF (lstomata) THEN 
        DEALLOCATE(rstom_mn_in,rstom_mx_in)
     ENDIF
     IF (.NOT.lconstvegalb) THEN 
        DEALLOCATE(vegalb_in)
     ENDIF 
     IF (lgettcl) THEN
        DEALLOCATE(t_cl_in)
     ENDIF
     !IF (itype_hydparam==4) THEN
     !   DEALLOCATE(pclay_in,psand_in)
     !ENDIF

  ENDIF

  ! Fill rlat and rlon field
  DO i=1,ie
     DO j=1,je
        rlat(i,j)=phirot2phi(startlat+(j-1)*dlat, & !SZB
         startlon+(i-1)*dlon,pollat,pollon,0.0_wp)  ! XYZ> added polgam=0 for compatability with utilities
        rlon(i,j)=rlarot2rla(startlat+(j-1)*dlat, & !SZB
         startlon+(i-1)*dlon,pollat,pollon,0.0_wp)  ! XYZ> added polgam=0 for compatability with utilities
     ENDDO
  ENDDO

  ! Diagnostic Ouput on the screen
  IF (lcheck) THEN
     CALL stats("rlon",rlon,llandmask)
     CALL stats("rlat",rlat,llandmask)
     CALL stats("hsurf",hsurf,llandmask)
     CALL stats("fr_land",fr_land,llandmask)
     CALL stats("soiltyp",REAL(soiltyp,wp),llandmask)
     CALL stats("z0",z0,llandmask)
     CALL stats("plcov_mn",plcov_mn,llandmask)
     CALL stats("plcov_mx",plcov_mx,llandmask)
     CALL stats("plcov",plcov,llandmask)
     CALL stats("lai_mn",lai_mn,llandmask)
     CALL stats("lai_mx",lai_mx,llandmask)
     CALL stats("lai",lai,llandmask)
     CALL stats("rootdp",rootdp,llandmask)
     CALL stats("rstom_mn",rstom_mn,llandmask)
     CALL stats("rstom_mx",rstom_mx,llandmask)
     IF (lstomata) CALL stats("rsmin2d",rsmin2d,llandmask) ! GDM CHANGED
     CALL stats("vegalb",vegalb,llandmask)
  ENDIF

END SUBROUTINE read_const_fields
!----------------------------------------------------------------------


SUBROUTINE read_const_fields_icon (month)
!----------------------------------------------------------------------
! Description:
!
!  ICON: Read assign fields to corresponding TERRA fields
!  from external files or from periodic vegetation assumption 
!----------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
  INTEGER, INTENT(IN), OPTIONAL                  :: &
     month

! Local variables

  LOGICAL                                        :: &
       is_found                                       ! 
  INTEGER                                        :: &
       numberofpoints                             ,&
       ilev                                           ! 
   
 CHARACTER (LEN=1)                               :: &
       uuid_in(16)
  INTEGER                                        :: &
       imiss, i, j, ii
  REAL (KIND=wp), ALLOCATABLE                    :: &
       rlon_in(:,:)                               , & ! 
       rlat_in(:,:)                               , & ! 
       hsurf_in(:,:)                              , & ! 
       fr_land_in(:,:)                            , & ! 
       soiltyp_in(:,:)                            , & ! 
       z0_in(:,:)                                 , & ! 
       plcov_mn_in(:,:)                           , & ! 
       plcov_mx_in(:,:)                           , & ! 
       plcov_in(:,:)                              , & ! 
       lai_mn_in(:,:)                             , & ! 
       lai_mx_in(:,:)                             , & ! 
       lai_in(:,:)                                , & ! 
       vegalb_in(:,:)                             , & ! 
       rstom_mn_in(:,:)                           , & ! 
       rstom_mx_in(:,:)                           , & ! 
       rootdp_in(:,:)                             , & ! 
       t_cl_in(:,:)                                   ! 
  !    psand_in(:,:)                              , & ! 
  !    psand(:,:)                                 , & ! 
  !    pclay_in(:,:)                              , & ! 
  !    pclay(:,:)                                     ! 
  CHARACTER (LEN=100)                            :: constfilename2

!!
  ! homogeneous soil/surface paramters
  IF (lhomosoil) THEN                   
     soiltyp  = soiltyp_const
     rootdp   = rootdp_const
     rootdp0  = rootdp_const
     plcov    = plcov_const    
     plcov_mn = plcov_min_const   
     plcov_mx = plcov_max_const           
     lai      = lai_const
     lai_mn   = lai_min_const
     lai_mx   = lai_max_const
     z0       = z0_const
!GDM>
     llandmask=.TRUE. ! GDM added for call to stats
     IF (.NOT. lstomata) THEN
!GDM<
        rstom_mn=crsmin
        rstom_mx=crsmax
     ELSE
        rstom_mn=rstom_mn_const
        rstom_mx=rstom_mx_const
     ENDIF
     IF (lconstvegalb) THEN
        vegalb=csalb_p
     ELSE
        vegalb=vegalb_const
     ENDIF

  ! in CASE soil/surface paramters provided by a grib-FILE
  ELSE                           

     IF (PRESENT(month)) THEN
        WRITE(constfilename2,'(A,A,I2.2,A)') TRIM(constfilename),'_2006',month,'0100'
     ELSE
        constfilename2=constfilename
     ENDIF


     ! determine size and position of input files
     CALL get_grib_info_icon (constfilename2,                &
          ngrid,ngpi,numberofpoints,uuid_in,is_found)

     IF ( .NOT. is_found ) THEN
        WRITE(6,*) "File ",TRIM(constfilename2)," does not exits!"," LOC B"

        STOP 
     ENDIF
     IF (ANY(uuid_in(:) /= uuid(:))) THEN
        WRITE(6,*) "uuid_in in read_const_field_icon does not match"
        STOP 'WRONG UUID'
     ENDIF

     ! allocate input fields
     ALLOCATE (                    &
          rlon_in( nproma,nblock)   , &
          rlat_in( nproma,nblock)   , &
          hsurf_in( nproma,nblock)  , &  
          fr_land_in( nproma,nblock), &
          soiltyp_in( nproma,nblock), &
          z0_in( nproma,nblock)     , &
          rootdp_in( nproma,nblock) )
     rlon_in=rundef
     rlat_in=rundef
     hsurf_in=rundef
     fr_land_in=rundef
     soiltyp_in=rundef
     z0_in=rundef
     rootdp_in=rundef
     IF (lvegadapt) THEN
        ALLOCATE (                         &
             plcov_mn_in( nproma,nblock) , &
             plcov_mx_in( nproma,nblock) , &
             lai_mn_in( nproma,nblock)   , &
             lai_mx_in( nproma,nblock)   )
        plcov_mn_in=rundef
        plcov_mx_in=rundef
        lai_mn_in=rundef
        lai_mx_in=rundef
     ELSE
        ALLOCATE (                     &
             plcov_in( nproma,nblock)   , &
             lai_in( nproma,nblock) )
        plcov_in=rundef
        lai_in=rundef
     ENDIF
     IF (lstomata) THEN 
        ALLOCATE (                     &
             rstom_mn_in( nproma,nblock), &
             rstom_mx_in( nproma,nblock) )
        rstom_mn_in=rundef
        rstom_mx_in=rundef
     ENDIF
     IF (.NOT. lconstvegalb) THEN 
        ALLOCATE (vegalb_in( nproma,nblock))
        vegalb_in=rundef
     ENDIF
     IF (lgettcl) THEN
        ALLOCATE (t_cl_in( nproma,nblock))
        t_cl_in=rundef
     ENDIF
     !    IF (itype_hydparam==4) THEN
     !       ALLOCATE(pclay_in(ie_in,je_in), &
     !            pclay(ie,je), & 
     !            psand_in(ie_in,je_in), &
     !            psand(ie,je))
     !       pclay=rundef
     !       psand_in=rundef
     !       psand=rundef
     !    ENDIF

     ! READ data from gribfile
     ! Design list of required fields
!!DL CALL init_griblist(99,ie_in,je_in,1)
     CALL init_griblist(99,nproma,nblock,1)
     CALL add2griblist(rlon_in,     'CLON','surface')
     CALL add2griblist(rlat_in,     'CLAT','surface')
     CALL add2griblist(hsurf_in,    'HSURF','surface')
     CALL add2griblist(fr_land_in,  'FR_LAND','surface')
     CALL add2griblist(soiltyp_in,  'SOILTYP','surface')
     IF (lz0local) THEN
        WRITE(6,*) "Z0LOC not implemented yet!"
!! DL   CALL add2griblist(z0_in,    82,250,1)   ! GDM> new numbers instead of 100/202 (SKYVIEW)
     ELSE
        CALL add2griblist(z0_in,    'Z0','surface')
     ENDIF
     IF (lvegadapt) THEN
        CALL add2griblist(plcov_mn_in, 'PLCOV_MN','surface')
        CALL add2griblist(plcov_mx_in, 'PLCOV_MX','surface')
        CALL add2griblist(lai_mn_in,   'LAI_MN','surface')
        CALL add2griblist(lai_mx_in,   'LAI_MX','surface')
     ELSE
        CALL add2griblist(plcov_in, 'PLCOV','surface')
        CALL add2griblist(lai_in,   'LAI','surface')
     ENDIF
     CALL add2griblist(rootdp_in,   'ROOTDP','surface')
     IF (lstomata) THEN
        WRITE(6,*) "lstomata=true not implementd yet!"
!! DL   CALL add2griblist(rstom_mn_in, 212,201,1) ! GDM> new numbers instead of 84/202 (AER_SO4)
! GDM>  CALL add2griblist(rstom_mx_in,  85,202,1)
     ENDIF
     IF (.NOT.lconstvegalb) THEN
        WRITE(6,*) "lconstvegalb=false not implementd yet!"
!! DL   CALL add2griblist(vegalb_in,   213,201,1) ! GDM> new numbers instead of 86/202 (AER_DUST)
     ENDIF
!    IF (lgettcl) THEN
!       IF (lmulti_layer) THEN
!          WRITE(6,*) "lgettcl / lmulti_layer =true not implemented yet!"
!          ilev=INT(czmls(ke_soil+1)*100+0.99)
!! GRIB1  CALL add2griblist(t_cl_in(:,:),'T_SO','depthBelowLand',0,ilev)
!! GRIB2  CALL add2griblist(t_cl_in(:,:),'T_2M_CL','depthBelowLand',0,ilev)
!! DL     CALL add2griblist(t_cl_in(:,:),197,201,111,0,ilev)
!       ELSE
!          CALL add2griblist(t_cl_in(:,:),'T_CL_LM','depthBelowLand',0,41)
!       ENDIF
!    ENDIF
     !    IF (itype_hydparam==4) THEN
     !        CALL add2griblist(psand_in,  98,202,1)
     !        CALL add2griblist(pclay_in,  97,202,1)
     !     ENDIF

!-------------------------------------------
     CALL read_grib_eccodes(constfilename2)
!-------------------------------------------
     CALL clear_griblist()

!
!DL ##############################################################################
!   Check, if FR_LAND, SOILTYP, HSURF etc are available
    imiss = 0
    IF (ALL(fr_land_in==rundef)) THEN
       imiss = imiss + 1
       WRITE(6,*) "FR_LAND is missing!"
    ELSEIF (ALL(soiltyp_in==rundef)) THEN
       imiss = imiss + 1
       WRITE(6,*) "SOILTYP is missing!"
    ELSEIF (ALL(hsurf_in==rundef)) THEN
       imiss = imiss + 1
       WRITE(6,*) "HSURF is missing!"
    ELSEIF (ALL(z0_in==rundef)) THEN
       imiss = imiss + 1
       WRITE(6,*) "z0 is missing!"
    ELSEIF (ALL(rootdp_in==rundef)) THEN
       imiss = imiss + 1
       WRITE(6,*) "ROOTDP is missing!"
    ENDIF
    IF(lvegadapt) THEN
       IF (ALL(plcov_mn_in==rundef)) THEN
!!DL      imiss = imiss + 1
          WRITE(6,*) "PLCOV_MN is missing!"
       ELSEIF (ALL(plcov_mx_in==rundef)) THEN
          imiss = imiss + 1
          WRITE(6,*) "PLCOV_MX is missing!"
       ELSEIF (ALL(lai_mn_in==rundef)) THEN
!!DL      imiss = imiss + 1
          WRITE(6,*) "LAI_MN is missing!"
       ELSEIF (ALL(lai_mx_in==rundef)) THEN
          imiss = imiss + 1
          WRITE(6,*) "LAI_MX is missing!"
       ENDIF 
    ELSE
       IF (ALL(plcov_in==rundef)) THEN
          imiss = imiss + 1
          WRITE(6,*) "PLCOV is missing!"
       ELSEIF (ALL(lai_in==rundef)) THEN
          imiss = imiss + 1
          WRITE(6,*) "LAI is missing!"
       ENDIF
    ENDIF
    IF (imiss > 0) THEN
       WRITE(6,*) "!!! There are ",imiss," fields missing in read_const_fields_icon!!!"
       STOP 'read_const_fields_icon'
    ENDIF
!DL ##############################################################################

     ! YUCHKDAT
     IF (lcheck) THEN
        OPEN(UNIT=444,FILE='YUCHKDAT',FORM='FORMATTED',STATUS='UNKNOWN')
        WRITE(444,*)
        WRITE(444,*) 'Reading ',TRIM(constfilename2)
        CALL stats("rlon",rlon_in,UNIT=444)
        CALL stats("rlat",rlat_in,UNIT=444)
        CALL stats("hsurf",hsurf_in,UNIT=444)
        CALL stats("fr_land",fr_land_in,UNIT=444)
        CALL stats("soiltyp",soiltyp_in,UNIT=444)
        CALL stats("z0",z0_in,UNIT=444)
        IF (lvegadapt) THEN
           CALL stats("plcov_mn",plcov_mn_in,UNIT=444)
           CALL stats("plcov_mx",plcov_mx_in,UNIT=444)
           CALL stats("lai_mn",lai_mn_in,UNIT=444)
           CALL stats("lai_mx",lai_mx_in,UNIT=444)
        ELSE
           CALL stats("plcov",plcov_in,UNIT=444)
           CALL stats("lai",lai_in,UNIT=444)
        ENDIF
        CALL stats("rootdp",rootdp_in,UNIT=444)
        IF (lstomata) THEN  ! GDM CHANGED
           CALL stats("rstom_mn",rstom_mn_in,UNIT=444)
           CALL stats("rstom_mx",rstom_mx_in,UNIT=444)
        ENDIF
        IF (.NOT.lconstvegalb) THEN
           CALL stats("vegalb",vegalb_in,UNIT=444)
        ENDIF
        IF (lgettcl) THEN
           CALL stats("t_cl",t_cl_in,UNIT=444)
        ENDIF
        CLOSE(UNIT=444)
     ENDIF


     rlon_geo    = rlon_in
     rlat_geo    = rlat_in
     rlon        = rlon_in
     rlat        = rlat_in
     fr_land     = fr_land_in
     hsurf       = hsurf_in
     soiltyp     = soiltyp_in
     z0          = z0_in
     IF (lvegadapt) THEN
        plcov_mn = plcov_mn_in
        plcov_mx = plcov_mx_in
        lai_mn   = lai_mn_in
        lai_mx   = lai_mx_in
     ELSE
        plcov    = plcov_in
        lai      = lai_in
     ENDIF
     rootdp      = rootdp_in
     llandmask   = fr_land>=0.5     !SZB
     IF (lstomata) THEN 
        rstom_mn = rstom_mn_in
        rsmin2d  = rstom_mn
        rstom_mx = rstom_mn
     ELSE
        rstom_mn = crsmin
        rstom_mx = crsmax
     ENDIF
     IF (.NOT.lconstvegalb) THEN 
        vegalb   = vegalb_in
     ELSE
        vegalb   = csalb_p
     ENDIF
     IF (lgettcl) THEN
        t_cl     = t_cl_in
     ENDIF

     ! set initial values for variables exhibiting an annual cycle
     IF (lvegadapt) THEN
        plcov    = plcov_mx
        rootdp0  = rootdp                                                     
        lai      = lai_mx
     ENDIF

     ! Convert plcov from % to fraction, if needed (ecoclimap)
     IF (MAXVAL(plcov)>2) THEN
        plcov=plcov/100.0_wp
     ENDIF

     ! Input fields are no longer needed
     DEALLOCATE(rlon_in,rlat_in,hsurf_in,fr_land_in,soiltyp_in,z0_in,rootdp_in)
     IF (lvegadapt) THEN
        DEALLOCATE(plcov_mn_in,plcov_mx_in,lai_mn_in,lai_mx_in)
     ELSE
        DEALLOCATE(plcov_in,lai_in)
     ENDIF
     IF (lstomata) THEN 
        DEALLOCATE(rstom_mn_in,rstom_mx_in)
     ENDIF
     IF (.NOT.lconstvegalb) THEN 
        DEALLOCATE(vegalb_in)
     ENDIF 
     IF (lgettcl) THEN
        DEALLOCATE(t_cl_in)
     ENDIF
     !IF (itype_hydparam==4) THEN
     !   DEALLOCATE(pclay_in,psand_in)
     !ENDIF

     IF (ALL (plcov_mn == rundef) .OR. ALL(lai_mn == rundef) ) THEN
        lvegadapt =.FALSE.
        WRITE (6,*) "lvegadapt is set to .FALSE. because of missing plcov_mn/lai_mn!!"
     ENDIF

  ENDIF

  ! Diagnostic Ouput on the screen
  IF (lcheck) THEN
     CALL stats("rlon",rlon,llandmask)
     CALL stats("rlat",rlat,llandmask)
     CALL stats("hsurf",hsurf,llandmask)
     CALL stats("fr_land",fr_land,llandmask)
     CALL stats("soiltyp",REAL(soiltyp,wp),llandmask)
     CALL stats("z0",z0,llandmask)
     CALL stats("plcov_mn",plcov_mn,llandmask)
     CALL stats("plcov_mx",plcov_mx,llandmask)
     CALL stats("plcov",plcov,llandmask)
     CALL stats("lai_mn",lai_mn,llandmask)
     CALL stats("lai_mx",lai_mx,llandmask)
     CALL stats("lai",lai,llandmask)
     CALL stats("rootdp",rootdp,llandmask)
     CALL stats("rstom_mn",rstom_mn,llandmask)
     CALL stats("rstom_mx",rstom_mx,llandmask)
     IF (lstomata) CALL stats("rsmin2d",rsmin2d,llandmask) ! GDM CHANGED
     CALL stats("vegalb",vegalb,llandmask)
  ENDIF

END SUBROUTINE read_const_fields_icon
!----------------------------------------------------------------------


SUBROUTINE read_initial_fields
!----------------------------------------------------------------------
! Description:
!
! Read initial fields and assign them to corresponding TERRA fields 
!----------------------------------------------------------------------


!  USE fxtr_definition, ONLY : pc_std_p0

  IMPLICIT NONE

  ! Parameters
  INTEGER,                  PARAMETER          :: ke_soil_old=3
  REAL    (KIND=wp),        PARAMETER          :: &
       zsm_level_old(ke_soil_old+1)             = &
         (/ 0.00_wp, 0.10_wp, 1.00_wp, 1.90_wp/), &
          zt_level_old(ke_soil_old)             = &
         (/ 0.00_wp, 0.09_wp, 0.49_wp/)

  ! Local Variables
  LOGICAL                                      :: &
       lcrop                                    , & !
       is_found                                     !
  REAL (KIND=wp)                               :: &
       dt                                       , & !
       fac                                          !
  CHARACTER (LEN=100)                          :: &
       filename                                 , & ! 
       filename1                                , & ! 
       filename2                                    !
  CHARACTER (LEN=2)                            :: str
  INTEGER                                      :: &
       ie_in                                    , & !
       je_in                                    , & !
       kso                                      , & !
       ilev                                     , & !
       idest                                    , & !
       k                                        , & !
       k1                                       , & !
       i                                        , & !
       j                                        , & !
       i0_crop                                  , & !
       j0_crop                                  , & !
       imin                                     , & !
       imax                                     , & !
       jmin                                     , & !
       jmax                                     , & !
       cvalid                                   , & !
       ir                                       , & !
!DL for tests
       icount                                   , & !
       i2                                       , & !
       j2                                           !
  REAL (KIND=wp)                               :: &
       rlat0_in                                 , & !
       rlon0_in                                 , & !
       dlat_in                                  , & !
       dlon_in                                  , & !
       pollon_in                                , & !
       pollat_in                                , & !
       rlat_geo                                 , & !
       rlon_geo                                 , & !
       rlat_rot_in                              , & !
       rlon_rot_in                                  !
  REAL (KIND=wp), TARGET, ALLOCATABLE          :: &
       t_snow_in(:,:)                           , & !
       w_snow_in(:,:)                           , & !
       rho_snow_in(:,:)                         , & !
       freshsnow_in(:,:)                        , & !
       t_so_in(:,:,:)                           , & !
       w_so_in(:,:,:)                           , & !
       w_i_in(:,:)                              , & !
       qv_s_in(:,:)                             , & !
       soil_in(:,:)                             , & ! GDM
       rlat_in(:,:)                             , & !
       rlon_in(:,:)                             , & !
       fr_land_in(:,:)                          , & !
       t_s_in(:,:)                              , & !
       t_m_in(:,:)                              , & !
       t_cl_in(:,:)                             , & !
       w_g1_in(:,:)                             , & !
       w_g2_in(:,:)                             , & !
       w_cl_in(:,:)                             , & !
       tch_in(:,:)                              , & !
       t_g_in(:,:)                              , & ! GDM
       tcm_in(:,:)                                  !
  LOGICAL,        ALLOCATABLE                  :: mask1(:,:)
  REAL (KIND=wp), ALLOCATABLE                  :: wmat(:,:)


  ! homogenoeus initial soil conditions
  IF (lhomoinit) THEN
     t_snow=t_soil0(1)
     t_s=t_soil0(1)
     t_g=t_soil0(1)
!    IF (lmulti_layer) THEN
        DO k=1,ke_soil
           t_so(:,:,k,:)=t_soil0(k)
        ENDDO
        t_so(:,:,0,:)=t_soil0(1)
        t_so(:,:,ke_soil+1,:)=t_cl0
!    ELSE
!       t_m=t_soil0(2)
!       t_cl=t_cl0
!    ENDIF
     w_snow=w_snow0
     rho_snow=250.0_wp
     freshsnow=0.6_wp
     w_i=w_i0
     qv_s=qvsat(t_s(1,1,nnow),REAL(pc_std_p0,wp))
     tch=1.0_wp
     tcm=1.0_wp
     IF (lrel_in) THEN
!       IF (lmulti_layer) THEN
           w_g0(1)   = w_g0(1) * czhls_const(1)* cporv(soiltyp_const)
           DO k=2,ke_soil
              w_g0(k)= w_g0(k)*(czhls_const(k)-czhls_const(k-1)) &
                       * cporv(soiltyp_const)
           ENDDO
           w_cl0     = w_cl0*(czhls_const(ke_soil+1)-czhls_const(ke_soil)) &
                       * cporv(soiltyp_const)
!       ELSE
!          w_g0(1)=w_g0(1)*cdzw12*cporv(soiltyp_const)
!          w_g0(2)=w_g0(2)*cdzw22*cporv(soiltyp_const)
!       ENDIF
     ENDIF
     IF (lvol_in) THEN
!       IF (lmulti_layer) THEN
           w_g0(1)=w_g0(1) * czhls_const(1)
           DO k=2,ke_soil
              w_g0(k)=w_g0(k)*(czhls_const(k)-czhls_const(k-1))
           ENDDO
           w_cl0=w_cl0*(czhls_const(ke_soil+1)-czhls_const(ke_soil))
!       ELSE
!          w_g0(1)=w_g0(1)*cdzw12
!          w_g0(2)=w_g0(2)*cdzw22
!       ENDIF
     ENDIF
!    IF (lmulti_layer) THEN
        DO k=1,ke_soil
           w_so(:,:,k,:)=w_g0(k)
        ENDDO
        w_so(:,:,ke_soil+1,:)=w_cl0
        w_so_ice=0.0_wp
!    ELSE
!       w_g1=w_g0(1) 
!       w_g2=w_g0(2)
!       w_g3=0.10_wp
!       w_cl=w_cl0
!    ENDIF

  ELSE ! Read gribfile including initial conditions

     filename1=TRIM(soilinitdir)//TRIM(soilinitprefix)//TRIM(ydate_ini)
     filename2=TRIM(filename1)//'00'    ! File name variant with minutes

     ! get field size and location
!! DL
     IF (yform_read == 'apix') THEN
     filename=filename1
     CALL get_grib_info(filename, &
          ie_in,je_in,rlat0_in,rlon0_in,dlat_in,dlon_in,pollat_in,pollon_in,is_found)
     IF ( .NOT. is_found ) THEN
         filename=filename2
         CALL get_grib_info(filename, &
              ie_in,je_in,rlat0_in,rlon0_in,dlat_in,dlon_in,pollat_in,pollon_in,is_found)
     ENDIF
    
     ELSE

     filename=filename1
     CALL get_gribinfo(filename, &
          ie_in,je_in,rlat0_in,rlon0_in,dlat_in,dlon_in,pollat_in,pollon_in,is_found)
     IF ( .NOT. is_found ) THEN
         filename=filename2
         CALL get_gribinfo(filename, &
              ie_in,je_in,rlat0_in,rlon0_in,dlat_in,dlon_in,pollat_in,pollon_in,is_found)
     ENDIF

     ENDIF
!! DL
     IF ( .NOT. is_found ) THEN
        WRITE(6,*) "Neither ",TRIM(filename1)," nor ", TRIM(filename2), " has been found"
        STOP
     ENDIF

     ! allocate input fields
     ALLOCATE (w_so_in(ie_in,je_in,ke_soil+1),   &
               t_so_in(ie_in,je_in,0:ke_soil+1), &
               t_g_in(ie_in,je_in),              &   ! GDM
               w_snow_in(ie_in,je_in),           &
               rho_snow_in(ie_in,je_in),         &
               freshsnow_in(ie_in,je_in),        &
               t_snow_in(ie_in,je_in),           &
               w_i_in(ie_in,je_in),              &
               qv_s_in(ie_in,je_in),             &
               soil_in(ie_in,je_in),             &   ! GDM
               rlon_in(ie_in,je_in),             &
               rlat_in(ie_in,je_in),             &
               fr_land_in(ie_in,je_in),          &
               tch_in(ie_in,je_in),              &
               tcm_in(ie_in,je_in),              &
               mask1(ie_in,je_in) )    
     w_so_in=rundef
     t_so_in=rundef
     t_g_in=rundef   ! GDM
     w_snow_in=rundef
     rho_snow_in=rundef
     freshsnow_in=rundef
     t_snow_in=rundef
     w_i_in=rundef
     qv_s_in=rundef
     tch_in=rundef
     tcm_in=rundef
     rlon_in=rundef
     rlat_in=rundef
!    IF (lmulti_layer .AND. (.NOT. lmulti_in)) THEN
!       ALLOCATE (t_s_in(ie_in,je_in),  &
!                 t_m_in(ie_in,je_in),  &
!                 t_cl_in(ie_in,je_in), &
!                 w_g1_in(ie_in,je_in), &
!                 w_g2_in(ie_in,je_in), &
!                 w_cl_in(ie_in,je_in) )
!       t_s_in=rundef
!       t_m_in=rundef
!       t_cl_in=rundef
!       w_g1_in=rundef
!       w_g2_in=rundef
!       w_cl_in=rundef
!    ENDIF

     ! read gribfile
!! DL< 
     IF (yform_read == 'apix') THEN

     CALL init_griblist(99,ie_in,je_in,1) 
     CALL add2griblist(t_g_in(:,:),      'T_G', 'surface')
     CALL add2griblist(w_i_in(:,:),      'W_I','surface')
     CALL add2griblist(qv_s_in(:,:),     'QV_S','surface')
     CALL add2griblist(w_snow_in(:,:),   'W_SNOW','surface')
     CALL add2griblist(rho_snow_in(:,:), 'RHO_SNOW','surface')
     CALL add2griblist(freshsnow_in(:,:),'FRESHSNW','surface')
     CALL add2griblist(t_snow_in(:,:),   'T_SNOW','surface')
!! DL Why FR_LAND and SOILTYP again ?? - Because grid may be different to constant one??
     CALL add2griblist(fr_land_in(:,:),  'FR_LAND','surface')
     CALL add2griblist(soil_in(:,:),     'SOILTYP','surface')

!    IF (lmulti_layer .AND. (.NOT. lmulti_in)) THEN
!       CALL add2griblist(t_s_in(:,:), 'T_S','depthBelowLand',0,0)
!       CALL add2griblist(t_m_in(:,:), 'T_M','depthBelowLand',0,9)
!       CALL add2griblist(t_cl_in(:,:),'T_CL_M','depthBelowLand',0,41)
!       CALL add2griblist(w_g1_in(:,:),'W_G1','depthBelowLandLayer',0,10)
!       CALL add2griblist(w_g2_in(:,:),'W_G2','depthBelowLandLayer',10,100)
!       CALL add2griblist(w_cl_in(:,:),'W_CL','depthBelowLandLayer',100,190)
!    ELSE
        CALL add2griblist(t_so_in(:,:,0),'T_SO','depthBelowLand',0,0)  
!!DL Be careful, because level definitions are for GRIB1 !!
        DO kso=1,ke_soil+1
           ilev=INT(czmls(kso)*100.0+0.99)
           CALL add2griblist(t_so_in(:,:,kso),'T_SO','depthBelowLand',ilev,ilev)   
           CALL add2griblist(w_so_in(:,:,kso),'W_SO','depthBelowLand',ilev,ilev)   
!!DL       NO W_SO_ICE ????
        ENDDO
!    ENDIF
     CALL add2griblist(tch_in(:,:),'TCH','surface')
     CALL add2griblist(tcm_in(:,:),'TCM','surface')

     WRITE(6,*) 'Reading ',TRIM(filename)
     CALL read_grib_eccodes(TRIM(filename))  
     CALL clear_griblist()         

     ELSE

     CALL init_griblist(99,ie_in,je_in,1)
     CALL add2griblist(t_g_in(:,:), 11,  2,1) ! GDM
     CALL add2griblist(w_i_in(:,:),200,201,1)
     CALL add2griblist(qv_s_in(:,:),51,  2,1)
     CALL add2griblist(w_snow_in(:,:),    65,  2,1)
     CALL add2griblist(rho_snow_in(:,:), 133,201,1)
     CALL add2griblist(freshsnow_in(:,:),129,201,1)
     CALL add2griblist(t_snow_in(:,:),   203,201,1)
     CALL add2griblist(fr_land_in(:,:),81,2,1)
     CALL add2griblist(soil_in(:,:), 57,202,1) ! GDM
 
!    IF (lmulti_layer .AND. (.NOT. lmulti_in)) THEN
!       CALL add2griblist(t_s_in(:,:), 85,2,111,0,0)
!       CALL add2griblist(t_m_in(:,:), 85,2,111,0,9)
!       CALL add2griblist(t_cl_in(:,:),85,2,111,0,41)
!       CALL add2griblist(w_g1_in(:,:),86,2,112,0,10)
!       CALL add2griblist(w_g2_in(:,:),86,2,112,10,100)
!       CALL add2griblist(w_cl_in(:,:),86,2,112,100,190)
!    ELSE
        CALL add2griblist(t_so_in(:,:,0),197,201,111,0,0)
        DO kso=1,ke_soil+1
           ilev=INT(czmls(kso)*100.0+0.99)
           CALL add2griblist(t_so_in(:,:,kso),197,201,111,0,ilev)
           CALL add2griblist(w_so_in(:,:,kso),198,201,111,0,ilev)
        ENDDO
!    ENDIF
     CALL add2griblist(tch_in(:,:),171,201,1)
     CALL add2griblist(tcm_in(:,:),170,201,1)
 
     WRITE(6,*) 'Reading ',TRIM(filename)
     CALL read_grib(TRIM(filename))
     CALL clear_griblist()
 
     ENDIF 
!! DL >
     ! ... check for missing fields

     IF ( ALL(rho_snow_in==rundef) ) THEN
        rho_snow_in = 250.0_wp
     ENDIF
     IF ( ALL(freshsnow_in==rundef) ) THEN
        freshsnow_in = 0.6_wp
     ENDIF
     IF ( ALL(tch_in==rundef) ) THEN
        tch_in = 1.0_wp
     ENDIF
     IF ( ALL(tcm_in==rundef) ) THEN
        tcm_in = 1.0_wp
     ENDIF
     IF ( ALL(qv_s_in==rundef) ) THEN
        qv_s_in = 0.0_wp
     ENDIF
!! DL
     ! ... Stop if other fields are missing !!
     IF ( (ALL(fr_land_in==rundef)) .OR. (ALL(soil_in==rundef))) THEN
        WRITE(6,*) 'fr_land_in and/or soil_in are missing in initial fields -- STOP!!'
        STOP 'read_initial_fields'
     ENDIF

! DL Next part has to be moved after conversion of input data befor lcrop
!    or wsnow has to be devided by 1000.
     IF ( ALL(t_g_in==rundef)) THEN
        mask1=fr_land_in>=0.5_wp
!! DL   CALL tgcom ( t_g_in, t_snow_in, t_s_in,    &
        CALL tgcom ( t_g_in(:,:), t_snow_in(:,:), t_so_in(:,:,0),    &
! DL                w_snow_in(:,:), llandmask(:,:), ie_in, je_in, cf_snow, &
                    w_snow_in(:,:)/1000.0_wp, mask1(:,:), ie_in, je_in, cf_snow, &
                    1, ie, 1, je )
        WRITE(6,*) 'DL: tgcom called to compute missing t_g'
     ENDIF

     ! YUCHKDAT
     IF (lcheck) THEN
        OPEN(UNIT=444,FILE='YUCHKDAT',FORM='FORMATTED',STATUS='OLD',POSITION='APPEND')
        WRITE(444,*)
        WRITE(444,*) 'Reading ',TRIM(filename)
        CALL stats("t_g",t_g_in,UNIT=444) ! GDM
        CALL stats("w_i",w_i_in,UNIT=444)
        CALL stats("qv_s",qv_s_in,UNIT=444)
        CALL stats("w_snow",w_snow_in,UNIT=444)
        CALL stats("rho_snow",rho_snow_in,UNIT=444)
        CALL stats("freshsnow",freshsnow_in,UNIT=444)
        CALL stats("t_snow",t_snow_in,UNIT=444)
        CALL stats("fr_land",fr_land_in,UNIT=444)
        CALL stats("soiltyp",soil_in,UNIT=444) ! GDM
!       IF (lmulti_layer.AND.(.NOT.lmulti_in)) THEN
!          CALL stats("t_s",t_s_in,UNIT=444)
!          CALL stats("t_m",t_m_in,UNIT=444)
!          CALL stats("t_cl",t_cl_in,UNIT=444)
!          CALL stats("w_g1",w_g1_in,UNIT=444)
!          CALL stats("w_g2",w_g2_in,UNIT=444)
!          CALL stats("w_cl",w_cl_in,UNIT=444)
!       ELSE
           DO kso=0,ke_soil+1
              WRITE(str,'(i2.2)') kso
              CALL stats("t_so("//str//")",t_so_in(:,:,kso),UNIT=444)
           ENDDO
           DO kso=1,ke_soil+1
              WRITE(str,'(i2.2)') kso
              CALL stats("w_so("//str//")",w_so_in(:,:,kso),UNIT=444)
           ENDDO
!       ENDIF
        CALL stats("tch",tch_in,UNIT=444)
        CALL stats("tcm",tcm_in,UNIT=444)
        CLOSE(UNIT=444)
     ENDIF

     ! convert two-layer soil mositure and soil temperature to
     !  multi-layer data by interpolation
!    IF (lmulti_layer.AND.(.NOT.lmulti_in)) THEN

!       ALLOCATE ( wmat(ke_soil+1,ke_soil_old) ) ! weight matrix for interpolation
!       ! soil moisture
!       CALL gen_interpolmat(ke_soil_old+1,zsm_level_old,         &
!            ke_soil+2,(/ 0.0_wp ,czhls_const(:) /),.TRUE.,wmat)
!       w_g1_in=w_g1_in/(zsm_level_old(2)-zsm_level_old(1))
!       w_g2_in=w_g2_in/(zsm_level_old(3)-zsm_level_old(2))
!       w_cl_in=w_cl_in/(zsm_level_old(4)-zsm_level_old(3))
!       DO idest=1,ke_soil+1
!          w_so_in(:,:,idest)=              &
!               wmat(idest,1)*w_g1_in(:,:)+ &
!               wmat(idest,2)*w_g2_in(:,:)+ &
!               wmat(idest,3)*w_cl_in(:,:)
!          IF (idest==1) THEN
!             w_so_in(:,:,idest)= w_so_in(:,:,idest)*czhls_const(1)
!          ELSE
!             w_so_in(:,:,idest)= w_so_in(:,:,idest)*      & 
!                  (czhls_const(idest)-czhls_const(idest-1))
!          ENDIF
!       ENDDO

!       ! soil temperature
!       CALL gen_interpolmat(ke_soil_old,zt_level_old,           &
!            ke_soil+2,(/ 0.0_wp, czhls_const(:) /),.FALSE.,wmat)
!       DO idest=1,ke_soil+1
!          t_so_in(:,:,idest)=             &
!               wmat(idest,1)*t_s_in(:,:)+ &
!               wmat(idest,2)*t_m_in(:,:)+ &
!               wmat(idest,3)*t_cl_in(:,:)
!       ENDDO
!       t_so_in(:,:,0)=t_s_in(:,:)

!       DEALLOCATE(wmat)
!    ENDIF

     ! conversions of input variables
     w_so_in=w_so_in/1000.0_wp
     w_snow_in=w_snow_in/1000.0_wp
     w_i_in=w_i_in/1000.0_wp

     ! tell user what we are doing with initial parameters
     lcrop=.FALSE.
     IF ( (ABS(pollat-pollat_in)<1.0E-3_wp)       &
          .AND. (ABS(pollon-pollon_in)<1.0E-3_wp) &
          .AND. (ABS(dlon_in-dlon)<1.0E-3_wp)     &
          .AND. (ABS(dlat_in-dlat)<1.0E-3_wp) )   THEN
        i0_crop=NINT((startlon-rlon0_in)/dlon+1)
        j0_crop=NINT((startlat-rlat0_in)/dlat+1)
        IF ( (i0_crop>=1)                                                &
             .AND. (j0_crop>=1)                                          &
             .AND. (ABS((i0_crop-1)*dlon-startlon+rlon0_in)<1.0E-3_wp)   &
             .AND. (ABS((j0_crop-1)*dlat-startlat+rlat0_in)<1.0E-3_wp) ) THEN
           lcrop=.TRUE.
        ENDIF
     ENDIF

     IF (lcrop) THEN
        IF (i0_crop==1 .AND. j0_crop==1 .AND. ie==ie_in .AND. je==je_in) THEN
           WRITE(6,*) "Perfect match of initial data!"
        ELSE
           WRITE(6,*) "Croping initial data!"
        ENDIF
     ELSE
        WRITE(6,*) "Interpolating initial data!"
     ENDIF

     ! Calculate coordinates of each input grid point
     IF (.NOT.lcrop) THEN 
        IF ((pollat==pollat_in).AND.(pollon==pollon_in)) THEN
           DO i=1,ie_in
              rlon_in(i,:)=rlon0_in+dlon_in*(i-1)
           ENDDO
           DO j=1,je_in  
              rlat_in(:,j)=rlat0_in+dlat_in*(j-1)
           ENDDO
        ELSE
           DO i=1,ie_in 
              rlon_rot_in=(i-1)*dlon_in+rlon0_in
              DO j=1,je_in
                 rlat_rot_in=(j-1)*dlat_in+rlat0_in
                 ! XYZ> added polgam=0 for compatability with utilities
                 rlat_geo=phirot2phi(rlat_rot_in,rlon_rot_in,pollat_in,pollon_in,0.0_wp)
                 rlon_geo=rlarot2rla(rlat_rot_in,rlon_rot_in,pollat_in,pollon_in,0.0_wp)
                 rlat_in(i,j)=phi2phirot(rlat_geo,rlon_geo,pollat,pollon)
                 rlon_in(i,j)=rla2rlarot(rlat_geo,rlon_geo,pollat,pollon,0.0_wp)
                 ! XYZ<
              ENDDO
           ENDDO
        ENDIF
     ENDIF

    ! Spatial interpolation  
     mask1=fr_land_in>=0.5_wp    !SZB
!DL
     icount=0
!DL
     IF (lcrop) THEN
        !DR HACK LANDMASK
        fr_land=fr_land_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1) 
        !END DR HACK 
        t_so(:,:,:,nnew)   = t_so_in  (i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1,:)
        w_so(:,:,:,nnew)   = w_so_in  (i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1,:)
        w_snow(:,:,nnew)   = w_snow_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)
        rho_snow(:,:,nnew) = rho_snow_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)
        freshsnow(:,:)     = freshsnow_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)
        t_snow(:,:,nnew)   = t_snow_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)
        t_g(:,:,nnew)      = t_g_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)     ! GDM
        w_i(:,:,nnew)      = w_i_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)
        qv_s(:,:,nnew)     = qv_s_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)
        WHERE ( NINT(soil_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)) >=9   &
               .AND. llandmask) ! water
           w_so(:,:,1,nnew)=rundef
           t_snow(:,:,nnew)=rundef
           t_g(:,:,nnew)=rundef  ! GDM
           w_i(:,:,nnew)=rundef
           qv_s(:,:,nnew)=rundef
           t_so(:,:,0,nnew)=rundef
           w_snow(:,:,nnew)=rundef
           rho_snow(:,:,nnew)=rundef
           freshsnow(:,:)=rundef
        END WHERE
!DL
        icount = icount + COUNT(NINT(soil_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)) >=9   &
               .AND. llandmask)
        WHERE ( NINT(soil_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)) <=2  &
                .AND. llandmask) ! ice & rock
           w_so(:,:,1,nnew)=rundef
        END WHERE
!DL
        icount = icount + COUNT(NINT(soil_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)) <=2  &
                .AND. llandmask)
        DO k=1,ke_soil+1
           WHERE (w_so(:,:,1,nnew)==rundef)
              w_so(:,:,k,nnew)=rundef
           END WHERE
           WHERE (t_so(:,:,0,nnew)==rundef)
              t_so(:,:,k,nnew)=rundef
           END WHERE
        ENDDO
        tch(:,:) = tch_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)
        tcm(:,:) = tcm_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)
     ELSE
        fac=dlat_in/dlat
        CALL general_interpol2d(rlat_in,rlon_in,t_so_in(:,:,0),      &
             startlat,startlon,dlat,dlon,t_so(:,:,0,nnew),           &
             2,2,0.51_wp*fac,1.0_wp,mask1,llandmask)
        CALL general_interpol2d(rlat_in,rlon_in,w_snow_in(:,:),      &
             startlat,startlon,dlat,dlon,w_snow(:,:,nnew),           &
             2,2,0.51_wp*fac,1.0_wp,mask1,llandmask)
        CALL general_interpol2d(rlat_in,rlon_in,rho_snow_in(:,:),    &
             startlat,startlon,dlat,dlon,rho_snow(:,:,nnew),         &
             2,2,0.51_wp*fac,1.0_wp,mask1,llandmask)
        CALL general_interpol2d(rlat_in,rlon_in,freshsnow_in(:,:),   &
             startlat,startlon,dlat,dlon,freshsnow(:,:),             &
             2,2,0.51_wp*fac,1.0_wp,mask1,llandmask)
        CALL general_interpol2d(rlat_in,rlon_in,t_snow_in(:,:),      &
             startlat,startlon,dlat,dlon,t_snow(:,:,nnew),           &
             2,2,0.51_wp*fac,1.0_wp,mask1,llandmask)
        CALL general_interpol2d(rlat_in,rlon_in,w_i_in(:,:),         &
             startlat,startlon,dlat,dlon,w_i(:,:,nnew),              &
             2,2,0.51_wp*fac,1.0_wp,mask1,llandmask)
        CALL general_interpol2d(rlat_in,rlon_in,qv_s_in(:,:),        &
             startlat,startlon,dlat,dlon,qv_s(:,:,nnew),             &
             2,2,0.51_wp*fac,1.0_wp,mask1,llandmask)
        DO kso=1,ke_soil+1
           CALL general_interpol2d(rlat_in,rlon_in,t_so_in(:,:,kso), &
                startlat,startlon,dlat,dlon,t_so(:,:,kso,nnew),      &
                2,2,0.51_wp*fac,1.0_wp,mask1,llandmask)
           CALL general_interpol2d(rlat_in,rlon_in,w_so_in(:,:,kso), &
                startlat,startlon,dlat,dlon,w_so(:,:,kso,nnew),      &
                2,2,0.51_wp*fac,1.0_wp,mask1,llandmask)
        ENDDO
        CALL general_interpol2d(rlat_in,rlon_in,tch_in(:,:),         &
             startlat,startlon,dlat,dlon,tch(:,:),                   &
             2,2,0.51_wp*fac,1.0_wp,mask1,llandmask)
        CALL general_interpol2d(rlat_in,rlon_in,tcm_in(:,:),         &
             startlat,startlon,dlat,dlon,tcm(:,:),                   &
             2,2,0.51_wp*fac,1.0_wp,mask1,llandmask)
     ENDIF

     ! Deallocate input fields that are no longer needed
     DEALLOCATE(w_so_in,t_so_in,w_snow_in,rho_snow_in,t_snow_in,freshsnow_in, &
                t_g_in,w_i_in,qv_s_in,soil_in,rlat_in,rlon_in,fr_land_in,tch_in,tcm_in,mask1)
!    IF (lmulti_layer.AND.(.NOT.lmulti_in)) THEN
!       DEALLOCATE(t_s_in,t_m_in,t_cl_in,w_g1_in,w_g2_in,w_cl_in)
!    ENDIF
     !DR HACK LANDMASK
     ! ensure consistency between fr_land (llandmask) and soiltyp
     WHERE (fr_land<0.5_wp)     !SZB
        soiltyp=9
     END WHERE
     WHERE (fr_land>=0.5_wp .AND. NINT(soiltyp)>=9)     !SZB
        fr_land=0.0_wp
     END WHERE
     ! final land-sea mask according to soil-type map;
     ! all point with llandmask=.FALSE. will not be considered for TERRA simulations 
     llandmask= fr_land>=0.5_wp   !SZB
     !END DR HACK LANDMASK

     ! BUG fix
     DO i=1,ie
        DO j=1,je
           IF ( ANY(t_so(i,j,:,nnew) < 230.0_wp) ) THEN
              t_so(i,j,:,nnew)=rundef
              w_so(i,j,:,nnew)=rundef
           ENDIF
        ENDDO
     ENDDO

!DL
     WRITE(6,*) "In read_initial_fields: Undef.points on land: ", icount
!DL
     ! Gap filling
     IF (COUNT((w_so(:,:,1,nnew)==rundef).AND.llandmask(:,:))>0) THEN
        WRITE(6,*) "Performing gap filling in read_initial_fields ..."   
        w_so(:,:,:,nnow)   = w_so(:,:,:,nnew)
        t_so(:,:,:,nnow)   = t_so(:,:,:,nnew)
        w_i(:,:,nnow)      = w_i(:,:,nnew)
        qv_s(:,:,nnow)     = qv_s(:,:,nnew)
        w_snow(:,:,nnow)   = w_snow(:,:,nnew)
        rho_snow(:,:,nnow) = rho_snow(:,:,nnew)
        t_snow(:,:,nnow)   = t_snow(:,:,nnew)

        DO i=1,ie
           DO j=1,je
              IF ((w_so(i,j,1,nnew)==rundef).AND.llandmask(i,j)) THEN
                 DO ir=1,MAX(i,j,ie-i,je-j)
                    imin=MAX(1,i-ir)
                    jmin=MAX(1,j-ir)
                    imax=MIN(ie,i+ir)
                    jmax=MIN(je,j+ir)
                    cvalid=COUNT((w_so(imin:imax,jmin:jmax,1,nnow)>rundef).AND. &
                                 (llandmask(imin:imax,jmin:jmax)))
                    IF (cvalid>0) THEN
                       DO k=0,ke_soil+1
                          IF (k>0) THEN
                               w_so(i,j,k,nnew)=                               &
                               SUM(w_so(imin:imax,jmin:jmax,k,nnow),           &
                               ((w_so(imin:imax,jmin:jmax,k,nnow)>rundef)      &
                               .AND. (llandmask(imin:imax,jmin:jmax))))/       &
                               FLOAT(cvalid)
                          ENDIF
                          t_so(i,j,k,nnew)=                                    &
                          SUM(t_so(imin:imax,jmin:jmax,k,nnow),                &
                          ((t_so(imin:imax,jmin:jmax,k,nnow)>rundef)           &
                          .AND.(llandmask(imin:imax,jmin:jmax))))/             &
                          FLOAT(COUNT((t_so(imin:imax,jmin:jmax,k,nnow)>rundef)&
                          .AND. (llandmask(imin:imax,jmin:jmax))))
                       ENDDO
                       IF (t_snow(i,j,nnew)==rundef) THEN
                            t_snow(i,j,nnew) = SUM(t_snow(imin:imax,jmin:jmax,nnow), &
                            ((t_snow(imin:imax,jmin:jmax,nnow)>rundef)               &
                            .AND. (llandmask(imin:imax,jmin:jmax)))) /               &
                            FLOAT(COUNT((t_snow(imin:imax,jmin:jmax,nnow)>rundef)    &
                            .AND. (llandmask(imin:imax,jmin:jmax))))
                       ENDIF
                       IF (w_snow(i,j,nnew)==rundef) THEN
                            w_snow(i,j,nnew) = SUM(w_snow(imin:imax,jmin:jmax,nnow), &
                            ((w_snow(imin:imax,jmin:jmax,nnow)>rundef)               &
                            .AND. (llandmask(imin:imax,jmin:jmax)))) /               &
                            FLOAT(COUNT((w_snow(imin:imax,jmin:jmax,nnow)>rundef)    &
                            .AND. (llandmask(imin:imax,jmin:jmax))))
                       ENDIF
                       IF (rho_snow(i,j,nnew)==rundef) THEN
                            rho_snow(i,j,nnew) = SUM(rho_snow(imin:imax,jmin:jmax,nnow), &
                            ((rho_snow(imin:imax,jmin:jmax,nnow)>rundef)             &
                            .AND. (llandmask(imin:imax,jmin:jmax)))) /               &
                            FLOAT(COUNT((rho_snow(imin:imax,jmin:jmax,nnow)>rundef)  &
                            .AND. (llandmask(imin:imax,jmin:jmax))))
                       ENDIF
                       IF (freshsnow(i,j)==rundef) THEN
                            freshsnow(i,j) = SUM(freshsnow(imin:imax,jmin:jmax),     &
                            ((freshsnow(imin:imax,jmin:jmax)>rundef)                 &
                            .AND. (llandmask(imin:imax,jmin:jmax)))) /               &
                            FLOAT(COUNT((freshsnow(imin:imax,jmin:jmax)>rundef)      &
                            .AND. (llandmask(imin:imax,jmin:jmax))))
                       ENDIF
                       IF (w_i(i,j,nnew)==rundef) THEN
                            w_i(i,j,nnew) = SUM(w_i(imin:imax,jmin:jmax,nnow),       &
                            ((w_i(imin:imax,jmin:jmax,nnow)>rundef)                  &
                            .AND. (llandmask(imin:imax,jmin:jmax)))) /               &
                            FLOAT(COUNT((w_i(imin:imax,jmin:jmax,nnow)>rundef)       &
                            .AND. (llandmask(imin:imax,jmin:jmax))))
                       ENDIF
                       IF (qv_s(i,j,nnew)==rundef) THEN
                            qv_s(i,j,nnew) = SUM(qv_s(imin:imax,jmin:jmax,nnow),     &
                            ((qv_s(imin:imax,jmin:jmax,nnow)>rundef)                 &
                            .AND. (llandmask(imin:imax,jmin:jmax)))) /               &
                            FLOAT(COUNT((qv_s(imin:imax,jmin:jmax,nnow)>rundef)      &
                            .AND. (llandmask(imin:imax,jmin:jmax))))
                       ENDIF
                       IF (tch(i,j)==rundef) THEN
                            tch(i,j) = SUM(tch(imin:imax,jmin:jmax),                 &
                            ((tch(imin:imax,jmin:jmax)>rundef)                       &
                            .AND. (llandmask(imin:imax,jmin:jmax)))) /               &
                            FLOAT(COUNT((tch(imin:imax,jmin:jmax)>rundef)            &
                            .AND. (llandmask(imin:imax,jmin:jmax))))
                       ENDIF
                       IF (tcm(i,j)==rundef) THEN
                            tcm(i,j) = SUM(tcm(imin:imax,jmin:jmax),                 &
                            ((tcm(imin:imax,jmin:jmax)>rundef)                       &
                            .AND. (llandmask(imin:imax,jmin:jmax)))) /               &
                            FLOAT(COUNT((tcm(imin:imax,jmin:jmax)>rundef)            &
                            .AND. (llandmask(imin:imax,jmin:jmax))))
                       ENDIF
                       EXIT
                    ENDIF  !cvalid > 0
                 ENDDO  !ir
              ENDIF  ! w_so...==rundef...
           ENDDO  ! j
        ENDDO  ! i
     ENDIF  ! count(w_so)...==rundef...

     ! Override with t_cl from external parameters if requested
     IF (lgettcl) THEN
        WRITE(6,*) "Overriding climatological deep soil temperature from external parameters"
        t_so(:,:,ke_soil+1,nnew)=t_cl(:,:)
     ENDIF

! GDM> tgcom begin 
     ir=0
     DO i=1,ie 
       DO j=1,je 
         IF (t_g(i,j,nnew)<=0 .AND. t_g(i,j,nnew)/=rundef) THEN
            ir=ir+1
         ENDIF
       ENDDO
     ENDDO
     IF (ir /= 0) THEN
       !XYZ> tgcom was commented by GDM. in order to uncomment,
       !      declare USE meteo_utilities ONLY :: tgcom in terra program
!      CALL tgcom ( t_g(:,:,nnew), t_snow(:,:,nnew), t_s(:,:,nnew),    &
!                   w_snow(:,:,nnew), llandmask(:,:), ie, je, cf_snow, &
!                   1, ie, 1, je )
       WRITE (*,*) "GDM: tgcom was not called to initialize t_g"
!!DL   WRITE (*,*) "DL : tgcom was     called to initialize t_g"
     ELSE
       IF (ntype_radinput==1) THEN
          WRITE (6,*) "GDM: tgcom is not called as t_g is available"
       ENDIF
     ENDIF 
! GDM< tgcom end

     ! Copy time levels
     w_so(:,:,:,nnow)   = w_so(:,:,:,nnew)
     t_so(:,:,:,nnow)   = t_so(:,:,:,nnew)
     w_i(:,:,nnow)      = w_i(:,:,nnew)
     qv_s(:,:,nnow)     = qv_s(:,:,nnew)
     w_snow(:,:,nnow)   = w_snow(:,:,nnew)
     rho_snow(:,:,nnow) = rho_snow(:,:,nnew)
     t_snow(:,:,nnow)   = t_snow(:,:,nnew)
     t_g(:,:,nnow)      = t_g(:,:,nnew)
     t_s(:,:,nnew)      = t_so(:,:,0,nnew)
     t_s(:,:,nnow)      = t_s(:,:,nnew)

     ! Attention: right now, no frozen soil water is initialized
     w_so_ice=0.0_wp

     ! Diagnostic Output
     IF (lcheck) THEN
        CALL stats("t_g",t_g(:,:,nnew))
        CALL stats("w_i",w_i(:,:,nnew))
        CALL stats("qv_s",qv_s(:,:,nnew))
        CALL stats("w_snow",w_snow(:,:,nnew))
        CALL stats("rho_snow",rho_snow(:,:,nnew))
        CALL stats("freshsnow",freshsnow(:,:))
        CALL stats("t_s",t_s(:,:,nnew))
        CALL stats("t_snow",t_snow(:,:,nnew))
        DO kso=0,ke_soil+1
           WRITE(str,'(i2.2)') kso
           CALL stats("t_so("//str//")",t_so(:,:,kso,nnew))
        ENDDO
        DO kso=1,ke_soil+1
           WRITE(str,'(i2.2)') kso
           CALL stats("w_so("//str//")",w_so(:,:,kso,nnew))
        ENDDO
        CALL stats("tch",tch(:,:))
        CALL stats("tcm",tcm(:,:))
     ENDIF

  ENDIF  ! lhomoinit

END SUBROUTINE read_initial_fields
!----------------------------------------------------------------------

SUBROUTINE read_initial_fields_icon
!----------------------------------------------------------------------
! Description:
!
! Read initial ICON fields and assign them to corresponding TERRA fields 
!----------------------------------------------------------------------


  IMPLICIT NONE

  ! Parameters
  INTEGER                 , PARAMETER          :: ke_soil_old=3
  REAL    (KIND=wp),        PARAMETER          :: &
       zsm_level_old(ke_soil_old+1)             = &
         (/ 0.00_wp, 0.10_wp, 1.00_wp, 1.90_wp/), &
          zt_level_old(ke_soil_old)             = &
         (/ 0.00_wp, 0.09_wp, 0.49_wp/)

  ! Local Variables
  LOGICAL                                      :: &
       is_found                                     !
  REAL (KIND=wp)                               :: &
       dt                                       , & !
       fac                                          !
  CHARACTER (LEN=100)                          :: &
       filename                                 , & ! 
       filename1                                , & ! 
       filename2                                    !
  CHARACTER (LEN=2)                            :: str
  INTEGER                                      :: &
       ie_in                                    , & !
       je_in                                    , & !
       kso                                      , & !
       ilev                                     , & !
       idest                                    , & !
       k                                        , & !
       k1                                       , & !
       i                                        , & !
       j                                        , & !
       imin                                     , & !
       imax                                     , & !
       jmin                                     , & !
       jmax                                     , & !
       cvalid                                   , & !
       ir                                       , & !
!DL for tests
       icount                                   , & !
       i2                                       , & !
       j2                                           !
!!DL new
  INTEGER                                      :: &
       ngp_in
  CHARACTER (LEN=1)                            :: &
       uuid_in(16)
  REAL (KIND=wp)                               :: &
       rlat0_in                                 , & !
       rlon0_in                                 , & !
       dlat_in                                  , & !
       dlon_in                                  , & !
       pollon_in                                , & !
       pollat_in                                , & !
       rlat_geo                                 , & !
       rlon_geo                                 , & !
       rlat_rot_in                              , & !
       rlon_rot_in                                  !
  REAL (KIND=wp), TARGET, ALLOCATABLE          :: &
       t_snow_in(:,:)                           , & !
       w_snow_in(:,:)                           , & !
       rho_snow_in(:,:)                         , & !
       freshsnow_in(:,:)                        , & !
       t_so_in(:,:,:)                           , & !
       w_so_in(:,:,:)                           , & !
       w_i_in(:,:)                              , & !
       qv_s_in(:,:)                             , & !
!!DL   soil_in(:,:)                             , & ! GDM
       rlat_in(:,:)                             , & !
       rlon_in(:,:)                             , & !
!!DL   fr_land_in(:,:)                          , & !
       t_s_in(:,:)                              , & !
       t_m_in(:,:)                              , & !
       t_cl_in(:,:)                             , & !
       w_g1_in(:,:)                             , & !
       w_g2_in(:,:)                             , & !
       w_cl_in(:,:)                             , & !
       tch_in(:,:)                              , & !
       t_g_in(:,:)                              , & ! GDM
       tcm_in(:,:)                                  !
!!DL  LOGICAL,        ALLOCATABLE                  :: mask1(:,:)
  REAL (KIND=wp), ALLOCATABLE                  :: wmat(:,:)


  ! homogenoeus initial soil conditions
  IF (lhomoinit) THEN
     t_snow=t_soil0(1)
     t_s=t_soil0(1)
     t_g=t_soil0(1)
!    IF (lmulti_layer) THEN
        DO k=1,ke_soil
           t_so(:,:,k,:)=t_soil0(k)
        ENDDO
        t_so(:,:,0,:)=t_soil0(1)
        t_so(:,:,ke_soil+1,:)=t_cl0
!    ELSE
!       t_m=t_soil0(2)
!       t_cl=t_cl0
!    ENDIF
     w_snow=w_snow0
     rho_snow=250.0_wp
     freshsnow=0.6_wp
     w_i=w_i0
     qv_s=qvsat(t_s(1,1,nnow),REAL(pc_std_p0,wp))
     tch=1.0_wp
     tcm=1.0_wp
     IF (lrel_in) THEN
!       IF (lmulti_layer) THEN
           w_g0(1)   = w_g0(1) * czhls_const(1)* cporv(soiltyp_const)
           DO k=2,ke_soil
              w_g0(k)= w_g0(k)*(czhls_const(k)-czhls_const(k-1)) &
                       * cporv(soiltyp_const)
           ENDDO
           w_cl0     = w_cl0*(czhls_const(ke_soil+1)-czhls_const(ke_soil)) &
                       * cporv(soiltyp_const)
!       ELSE
!          w_g0(1)=w_g0(1)*cdzw12*cporv(soiltyp_const)
!          w_g0(2)=w_g0(2)*cdzw22*cporv(soiltyp_const)
!       ENDIF
     ENDIF
     IF (lvol_in) THEN
!       IF (lmulti_layer) THEN
           w_g0(1)=w_g0(1) * czhls_const(1)
           DO k=2,ke_soil
              w_g0(k)=w_g0(k)*(czhls_const(k)-czhls_const(k-1))
           ENDDO
           w_cl0=w_cl0*(czhls_const(ke_soil+1)-czhls_const(ke_soil))
!       ELSE
!          w_g0(1)=w_g0(1)*cdzw12
!          w_g0(2)=w_g0(2)*cdzw22
!       ENDIF
     ENDIF
!    IF (lmulti_layer) THEN
        DO k=1,ke_soil
           w_so(:,:,k,:)=w_g0(k)
        ENDDO
        w_so(:,:,ke_soil+1,:)=w_cl0
        w_so_ice=0.0_wp
!    ELSE
!       w_g1=w_g0(1) 
!       w_g2=w_g0(2)
!       w_g3=0.10_wp
!       w_cl=w_cl0
!    ENDIF

  ELSE ! Read gribfile including initial conditions

     filename1=TRIM(soilinitdir)//TRIM(soilinitprefix)//TRIM(ydate_ini)
     filename2=TRIM(filename1)//'00'    ! File name variant with minutes

     ! get field size and location
     filename=filename1
     CALL get_grib_info_icon(filename, &
          ngrid,ngpi,ngp_in,uuid_in,is_found)
     IF ( .NOT. is_found ) THEN
         filename=filename2
         CALL get_grib_info_icon(filename, &
              ngrid,ngpi,ngp_in,uuid_in,is_found)
     ENDIF
    
     IF ( .NOT. is_found ) THEN
        WRITE(6,*) "Neither ",TRIM(filename1)," nor ", TRIM(filename2), " has been found"
        STOP
     ENDIF
     IF (ANY(uuid_in(:) /= uuid(:))) THEN
        WRITE(6,*) "Mismatch in uuid of ICON grid for initial fields!"
        STOP
     ENDIF

     ie_in = nproma
     je_in = nblock

     ! allocate input fields
     ALLOCATE (w_so_in(ie_in,je_in,ke_soil+1),   &
               t_so_in(ie_in,je_in,0:ke_soil+1), &
               t_g_in(ie_in,je_in),              &   ! GDM
               w_snow_in(ie_in,je_in),           &
               rho_snow_in(ie_in,je_in),         &
               freshsnow_in(ie_in,je_in),        &
               t_snow_in(ie_in,je_in),           &
               w_i_in(ie_in,je_in),              &
               qv_s_in(ie_in,je_in),             &
!!DL           soil_in(ie_in,je_in),             &   ! GDM
               rlon_in(ie_in,je_in),             &
               rlat_in(ie_in,je_in),             &
!!DL           fr_land_in(ie_in,je_in),          &
               tch_in(ie_in,je_in),              &
               tcm_in(ie_in,je_in)               )
     w_so_in=rundef
     t_so_in=rundef
     t_g_in=rundef   ! GDM
     w_snow_in=rundef
     rho_snow_in=rundef
     freshsnow_in=rundef
     t_snow_in=rundef
     w_i_in=rundef
     qv_s_in=rundef
     tch_in=rundef
     tcm_in=rundef
     rlon_in=rundef
     rlat_in=rundef
!    IF (lmulti_layer .AND. (.NOT. lmulti_in)) THEN
!       ALLOCATE (t_s_in(ie_in,je_in),  &
!                 t_m_in(ie_in,je_in),  &
!                 t_cl_in(ie_in,je_in), &
!                 w_g1_in(ie_in,je_in), &
!                 w_g2_in(ie_in,je_in), &
!                 w_cl_in(ie_in,je_in) )
!       t_s_in=rundef
!       t_m_in=rundef
!       t_cl_in=rundef
!       w_g1_in=rundef
!       w_g2_in=rundef
!       w_cl_in=rundef
!    ENDIF

     ! Define griblist and read gribfile

!!DL CALL init_griblist(99,ie_in,je_in,1) 
     CALL init_griblist(99,nproma,nblock,1) 
     CALL add2griblist(t_g_in(:,:),      'T_G', 'surface')
     CALL add2griblist(w_i_in(:,:),      'W_I','surface')
     CALL add2griblist(qv_s_in(:,:),     'QV_S','surface')
     CALL add2griblist(w_snow_in(:,:),   'W_SNOW','surface')
     CALL add2griblist(rho_snow_in(:,:), 'RHO_SNOW','surface')
     CALL add2griblist(freshsnow_in(:,:),'FRESHSNW','surface')
     CALL add2griblist(t_snow_in(:,:),   'T_SNOW','surface')
!! DL Why FR_LAND and SOILTYP again ?? - Because grid may be different to constant one??
!!DL CALL add2griblist(fr_land_in(:,:),  'FR_LAND','surface')
!!DL CALL add2griblist(soil_in(:,:),     'SOILTYP','surface')

!    IF (lmulti_layer .AND. (.NOT. lmulti_in)) THEN
!       CALL add2griblist(t_s_in(:,:), 'T_S','depthBelowLand',0,0)
!       CALL add2griblist(t_m_in(:,:), 'T_M','depthBelowLand',0,9)
!       CALL add2griblist(t_cl_in(:,:),'T_CL_M','depthBelowLand',0,41)
!       CALL add2griblist(w_g1_in(:,:),'W_G1','depthBelowLandLayer',0,10)
!       CALL add2griblist(w_g2_in(:,:),'W_G2','depthBelowLandLayer',10,100)
!       CALL add2griblist(w_cl_in(:,:),'W_CL','depthBelowLandLayer',100,190)
!    ELSE
        CALL add2griblist(t_so_in(:,:,0),'T_SO','depthBelowLand',0,0)  
!!DL Be careful, because level definitions are for GRIB1 !!
        DO kso=1,ke_soil+1
           ilev=INT(czmls(kso)*100.0+0.99)
           CALL add2griblist(t_so_in(:,:,kso),'T_SO','depthBelowLand',ilev,ilev)   
           CALL add2griblist(w_so_in(:,:,kso),'W_SO','depthBelowLand',ilev,ilev)   
!!DL       NO W_SO_ICE ????
        ENDDO
!    ENDIF
     CALL add2griblist(tch_in(:,:),'TCH','surface')
     CALL add2griblist(tcm_in(:,:),'TCM','surface')

     WRITE(6,*) 'Reading ',TRIM(filename)
     CALL read_grib_eccodes(TRIM(filename))  
     CALL clear_griblist()         

     ! ... check for missing fields

     IF ( ALL(rho_snow_in==rundef) ) THEN
        rho_snow_in = 250.0_wp
     ENDIF
     IF ( ALL(freshsnow_in==rundef) ) THEN
        freshsnow_in = 0.6_wp
     ENDIF
     IF ( ALL(tch_in==rundef) ) THEN
        tch_in = 1.0_wp
     ENDIF
     IF ( ALL(tcm_in==rundef) ) THEN
        tcm_in = 1.0_wp
     ENDIF
     IF ( ALL(qv_s_in==rundef) ) THEN
        qv_s_in = 0.0_wp
     ENDIF
!! DL
!!DL ! ... Stop if other fields are missing !!
!!DL IF ( (ALL(fr_land_in==rundef)) .OR. (ALL(soil_in==rundef))) THEN
!!DL    WRITE(6,*) 'fr_land_in and/or soil_in are missing in initial fields -- STOP!!'
!!DL    STOP 'read_initial_fields_icon'
!!DL ENDIF

! DL Next part has to be moved after conversion of input data befor lcrop
!    or wsnow has to be devided by 1000.
     IF ( ALL(t_g_in==rundef)) THEN
!! DL   CALL tgcom ( t_g_in, t_snow_in, t_s_in,    &
        CALL tgcom ( t_g_in(:,:), t_snow_in(:,:), t_so_in(:,:,0),    &
! DL                w_snow_in(:,:), llandmask(:,:), ie_in, je_in, cf_snow, &
                    w_snow_in(:,:)/1000.0_wp, llandmask(:,:), ie_in, je_in, cf_snow, &
                    1, ie, 1, je )
        WRITE(6,*) 'DL: tgcom called to compute missing t_g'
     ENDIF

     ! YUCHKDAT
     IF (lcheck) THEN
        OPEN(UNIT=444,FILE='YUCHKDAT',FORM='FORMATTED',STATUS='OLD',POSITION='APPEND')
        WRITE(444,*)
        WRITE(444,*) 'Reading ',TRIM(filename)
        CALL stats("t_g",t_g_in,UNIT=444) ! GDM
        CALL stats("w_i",w_i_in,UNIT=444)
        CALL stats("qv_s",qv_s_in,UNIT=444)
        CALL stats("w_snow",w_snow_in,UNIT=444)
        CALL stats("rho_snow",rho_snow_in,UNIT=444)
        CALL stats("freshsnow",freshsnow_in,UNIT=444)
        CALL stats("t_snow",t_snow_in,UNIT=444)
!!DL    CALL stats("fr_land",fr_land_in,UNIT=444)
!!DL    CALL stats("soiltyp",soil_in,UNIT=444) ! GDM
!       IF (lmulti_layer.AND.(.NOT.lmulti_in)) THEN
!          CALL stats("t_s",t_s_in,UNIT=444)
!          CALL stats("t_m",t_m_in,UNIT=444)
!          CALL stats("t_cl",t_cl_in,UNIT=444)
!          CALL stats("w_g1",w_g1_in,UNIT=444)
!          CALL stats("w_g2",w_g2_in,UNIT=444)
!          CALL stats("w_cl",w_cl_in,UNIT=444)
!       ELSE
           DO kso=0,ke_soil+1
              WRITE(str,'(i2.2)') kso
              CALL stats("t_so("//str//")",t_so_in(:,:,kso),UNIT=444)
           ENDDO
           DO kso=1,ke_soil+1
              WRITE(str,'(i2.2)') kso
              CALL stats("w_so("//str//")",w_so_in(:,:,kso),UNIT=444)
           ENDDO
!       ENDIF
        CALL stats("tch",tch_in,UNIT=444)
        CALL stats("tcm",tcm_in,UNIT=444)
        CLOSE(UNIT=444)
     ENDIF

     ! convert two-layer soil mositure and soil temperature to
     !  multi-layer data by interpolation
!    IF (lmulti_layer.AND.(.NOT.lmulti_in)) THEN

!       ALLOCATE ( wmat(ke_soil+1,ke_soil_old) ) ! weight matrix for interpolation
!       ! soil moisture
!       CALL gen_interpolmat(ke_soil_old+1,zsm_level_old,         &
!            ke_soil+2,(/ 0.0_wp ,czhls_const(:) /),.TRUE.,wmat)
!       w_g1_in=w_g1_in/(zsm_level_old(2)-zsm_level_old(1))
!       w_g2_in=w_g2_in/(zsm_level_old(3)-zsm_level_old(2))
!       w_cl_in=w_cl_in/(zsm_level_old(4)-zsm_level_old(3))
!       DO idest=1,ke_soil+1
!          w_so_in(:,:,idest)=              &
!               wmat(idest,1)*w_g1_in(:,:)+ &
!               wmat(idest,2)*w_g2_in(:,:)+ &
!               wmat(idest,3)*w_cl_in(:,:)
!          IF (idest==1) THEN
!             w_so_in(:,:,idest)= w_so_in(:,:,idest)*czhls_const(1)
!          ELSE
!             w_so_in(:,:,idest)= w_so_in(:,:,idest)*      & 
!                  (czhls_const(idest)-czhls_const(idest-1))
!          ENDIF
!       ENDDO

!       ! soil temperature
!       CALL gen_interpolmat(ke_soil_old,zt_level_old,           &
!            ke_soil+2,(/ 0.0_wp, czhls_const(:) /),.FALSE.,wmat)
!       DO idest=1,ke_soil+1
!          t_so_in(:,:,idest)=             &
!               wmat(idest,1)*t_s_in(:,:)+ &
!               wmat(idest,2)*t_m_in(:,:)+ &
!               wmat(idest,3)*t_cl_in(:,:)
!       ENDDO
!       t_so_in(:,:,0)=t_s_in(:,:)

!       DEALLOCATE(wmat)
!    ENDIF

     ! conversions of input variables
     w_so_in=w_so_in/1000.0_wp
     w_snow_in=w_snow_in/1000.0_wp
     w_i_in=w_i_in/1000.0_wp

     ! tell user what we are doing with initial parameters
        t_so(:,:,:,nnew)   = t_so_in(:,:,:)
        w_so(:,:,:,nnew)   = w_so_in(:,:,:)
        w_snow(:,:,nnew)   = w_snow_in(:,:)
        rho_snow(:,:,nnew) = rho_snow_in(:,:)
        freshsnow(:,:)     = freshsnow_in(:,:)
        t_snow(:,:,nnew)   = t_snow_in(:,:)
        t_g(:,:,nnew)      = t_g_in(:,:)     ! GDM
        w_i(:,:,nnew)      = w_i_in(:,:)
        qv_s(:,:,nnew)     = qv_s_in(:,:)
        WHERE ( NINT(soiltyp) >=9   &
               .AND. llandmask) ! water
           w_so(:,:,1,nnew)=rundef
           t_snow(:,:,nnew)=rundef
           t_g(:,:,nnew)=rundef  ! GDM
           w_i(:,:,nnew)=rundef
           qv_s(:,:,nnew)=rundef
           t_so(:,:,0,nnew)=rundef
           w_snow(:,:,nnew)=rundef
           rho_snow(:,:,nnew)=rundef
           freshsnow(:,:)=rundef
        END WHERE
        tch(:,:) = tch_in(:,:)
        tcm(:,:) = tcm_in(:,:)

     ! Deallocate input fields that are no longer needed
     DEALLOCATE(w_so_in,t_so_in,w_snow_in,rho_snow_in,t_snow_in,freshsnow_in, &
                t_g_in,w_i_in,qv_s_in,rlat_in,rlon_in,tch_in,tcm_in)
!!DL            t_g_in,w_i_in,qv_s_in,soil_in,rlat_in,rlon_in,fr_land_in,tch_in,tcm_in)
!    IF (lmulti_layer.AND.(.NOT.lmulti_in)) THEN
!       DEALLOCATE(t_s_in,t_m_in,t_cl_in,w_g1_in,w_g2_in,w_cl_in)
!    ENDIF

     ! BUG fix
!!   DO i=1,ie
!!      DO j=1,je
     DO i=1,nproma
        DO j=1,nblock
           IF ( ANY(t_so(i,j,:,nnew) < 230.0_wp) ) THEN
              t_so(i,j,:,nnew)=rundef
              w_so(i,j,:,nnew)=rundef
           ENDIF
        ENDDO
     ENDDO


     ! Override with t_cl from external parameters if requested
     IF (lgettcl) THEN
        WRITE(6,*) "Overriding climatological deep soil temperature from external parameters"
        t_so(:,:,ke_soil+1,nnew)=t_cl(:,:)
     ENDIF

! GDM> tgcom begin 
!!DL ir=0
!!DL DO i=1,ie 
!!DL   DO j=1,je 
!!DL     IF (t_g(i,j,nnew)<=0 .AND. t_g(i,j,nnew)/=rundef) THEN
!!DL        ir=ir+1
!!DL     ENDIF
!!DL   ENDDO
!!DL ENDDO
!!DL IF (ir /= 0) THEN
       !XYZ> tgcom was commented by GDM. in order to uncomment,
       !      declare USE meteo_utilities ONLY :: tgcom in terra program
!      CALL tgcom ( t_g(:,:,nnew), t_snow(:,:,nnew), t_s(:,:,nnew),    &
!                   w_snow(:,:,nnew), llandmask(:,:), ie, je, cf_snow, &
!                   1, ie, 1, je )
!!DL   WRITE (*,*) "GDM: tgcom was not called to initialize t_g"
!!DL   WRITE (*,*) "DL : tgcom was     called to initialize t_g"
!!DL ELSE
!!DL   IF (ntype_radinput==1) THEN
!!DL      WRITE (6,*) "GDM: tgcom is not called as t_g is available"
!!DL   ENDIF
!!DL ENDIF 
! GDM< tgcom end

     ! Copy time levels
     w_so(:,:,:,nnow)   = w_so(:,:,:,nnew)
     t_so(:,:,:,nnow)   = t_so(:,:,:,nnew)
     w_i(:,:,nnow)      = w_i(:,:,nnew)
     qv_s(:,:,nnow)     = qv_s(:,:,nnew)
     w_snow(:,:,nnow)   = w_snow(:,:,nnew)
     rho_snow(:,:,nnow) = rho_snow(:,:,nnew)
     t_snow(:,:,nnow)   = t_snow(:,:,nnew)
     t_g(:,:,nnow)      = t_g(:,:,nnew)
     t_s(:,:,nnew)      = t_so(:,:,0,nnew)
     t_s(:,:,nnow)      = t_s(:,:,nnew)

     ! Attention: right now, no frozen soil water is initialized
     w_so_ice=0.0_wp

     ! Diagnostic Output
     IF (lcheck) THEN
        CALL stats("t_g",t_g(:,:,nnew))
        CALL stats("w_i",w_i(:,:,nnew))
        CALL stats("qv_s",qv_s(:,:,nnew))
        CALL stats("w_snow",w_snow(:,:,nnew))
        CALL stats("rho_snow",rho_snow(:,:,nnew))
        CALL stats("freshsnow",freshsnow(:,:))
        CALL stats("t_s",t_s(:,:,nnew))
        CALL stats("t_snow",t_snow(:,:,nnew))
        DO kso=0,ke_soil+1
           WRITE(str,'(i2.2)') kso
           CALL stats("t_so("//str//")",t_so(:,:,kso,nnew))
        ENDDO
        DO kso=1,ke_soil+1
           WRITE(str,'(i2.2)') kso
           CALL stats("w_so("//str//")",w_so(:,:,kso,nnew))
        ENDDO
        CALL stats("tch",tch(:,:))
        CALL stats("tcm",tcm(:,:))
     ENDIF

  ENDIF  ! lhomoinit

END SUBROUTINE read_initial_fields_icon
!----------------------------------------------------------------------



SUBROUTINE read_lmgrib(date, nx, effective_date, nave, look_forward)
!----------------------------------------------------------------------
! Description:
!
! Read local model grib files and assign fields to corresponding TERRA fields 
!----------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
  CHARACTER (LEN=10),       INTENT(IN)  :: date
  CHARACTER (LEN=10),       INTENT(OUT) :: effective_date
  INTEGER                 , INTENT(IN)  :: nx
  INTEGER                 , INTENT(OUT) :: nave
  LOGICAL ,                 INTENT(IN)  :: look_forward

! Local Variables
  CHARACTER (LEN=100)         :: &
       filename                , & ! 
       filename1               , & ! 
       filename2                   ! 
  INTEGER                     :: &
       ie_in                   , & ! 
       je_in                   , & ! 
       i                       , & ! 
       j                       , & ! 
       i0_crop                 , & ! 
       j0_crop                 , & ! 
       ir                      , & ! 
       imin                    , & ! 
       imax                    , & ! 
       jmin                    , & ! 
       jmax                    , & ! 
       cvalid                  , & ! 
       tincr                       ! 
!DL
  INTEGER                     :: &
       imiss, imiss_tg
!DL
  INTEGER (KIND=kind_idate)   :: &
       idate                   , & ! 
       new_date                    ! 
  LOGICAL                     :: &
       lcrop                   , & ! 
       is_found                , & ! 
!DL
       lhourly_in              , & ! one hour time interval data ??
!DL
       any_undef                   ! 
  REAL (KIND=wp)              :: &
       rlat_geo                , & ! 
       rlon_geo                , & ! 
       rlat_rot_in             , & ! 
       rlon_rot_in             , & ! 
       rlat0_in                , & ! 
       rlon0_in                , & ! 
       dlat_in                 , & ! 
       dlon_in                 , & ! 
       pollon_in               , & ! 
       pollat_in                   ! 
  REAL (KIND=wp)              :: fac
  REAL (KIND=wp), ALLOCATABLE :: &
       rlon_in(:,:)            , & ! 
       rlat_in(:,:)            , & ! 
       u_in(:,:)               , & ! 
       v_in(:,:)               , & ! 
       t_in(:,:)               , & ! 
       qv_in(:,:)              , & ! 
       ps_in(:,:)              , & ! 
       rain_gsp_in(:,:)        , & ! 
       rain_con_in(:,:)        , & ! 
       snow_gsp_in(:,:)        , & ! 
       snow_con_in(:,:)        , & ! 
       sobs2_in(:,:)           , & ! 
!DL provide average fields as well asob2, athbs2, apabs
       asobs2_in(:,:)          , & ! 
       thbs2_in(:,:)           , & ! 
       athbs2_in(:,:)          , & ! 
       pabs_in(:,:)            , & ! 
       apabs_in(:,:)           , & ! 
       t_g_in(:,:)             , & ! 
       albrad_in(:,:)          , & ! 
       t_s_in(:,:)             , & ! 
!DL alternative to t_s_in
       t_so_0_in (:,:)         , & ! 
       w_snow_in(:,:)          , & ! 
       t_snow_in(:,:)          , & ! 
       gfill(:,:,:)                ! 
  LOGICAL, ALLOCATABLE        :: mask1(:,:)

!Declaration of STATEMENT-FUNCTIONS
!--------------------
  REAL (KIND=wp) :: fpvsw, zst, fqvs, zsge, zsp
  fpvsw(zst)       = b1 * EXP( b2w * (zst - b3) / (zst - b4w) )
  fqvs (zsge, zsp) = rdv * zsge / ( zsp - o_m_rdv * zsge )

!DL
  ! preset lhourly_in
  lhourly_in = lhourly_data
!DL

  ! Look for file, look forward when requested
  effective_date = date
  READ(date,*) idate
  idate = idate * 100
  tincr = 0

  DO       ! find undef if any 
     DO     ! find file
        filename1=TRIM(metfiledir)//TRIM(metfileprefix)//TRIM(effective_date)
        filename2=TRIM(filename1)//'00'    ! File name variant with minutes
        filename=filename1
        IF (yform_read == 'apix') THEN
        CALL get_grib_info(filename, ie_in,je_in,rlat0_in,rlon0_in,    &
                          dlat_in,dlon_in,pollat_in,pollon_in,is_found)
        IF ( .NOT. is_found ) THEN
          filename=filename2
          CALL get_grib_info(filename, ie_in,je_in,rlat0_in,rlon0_in,    &
                            dlat_in,dlon_in,pollat_in,pollon_in,is_found)
        ENDIF
        ELSE
        CALL get_gribinfo(filename, ie_in,je_in,rlat0_in,rlon0_in,    &
                          dlat_in,dlon_in,pollat_in,pollon_in,is_found)
        IF ( .NOT. is_found ) THEN
          filename=filename2
          CALL get_gribinfo(filename, ie_in,je_in,rlat0_in,rlon0_in,    &
                            dlat_in,dlon_in,pollat_in,pollon_in,is_found)
        ENDIF
        ENDIF

        IF ( is_found ) THEN
           EXIT
        ENDIF
        IF ( .NOT. look_forward ) THEN
           WRITE(6,*) "Neither ",TRIM(filename1)," nor ", TRIM(filename2), " have been found"
           STOP
        ELSE
print *, 'about tincr:  ', tincr, tincr_max, filename1, filename2
           IF ( tincr < tincr_max*60 ) THEN
              tincr = tincr + 60
              new_date = date_SUM(idate,tincr)
              WRITE(effective_date,'(i10.10)') new_date/100
              WRITE(6,*) "Neither ",TRIM(filename1)," nor ", TRIM(filename2), " have been found, ", &
                         "looking at ", TRIM(effective_date)
           ELSE
              WRITE(6,*) "Maximum gap reached without finding file (tincr_max= ", tincr_max, ")"
              STOP
           ENDIF
        ENDIF
     ENDDO  ! find file

     ! Allocation of fields depending on configuration defined by
     !  ntype_atminput, ntype_raininput, ntype_radinput
     IF ( ALLOCATED(rlon_in) )     DEALLOCATE(rlon_in)
     IF ( ALLOCATED(rlat_in) )     DEALLOCATE(rlat_in)
     IF ( ALLOCATED(mask1) )       DEALLOCATE(mask1)
     IF ( ALLOCATED(u_in) )        DEALLOCATE(u_in)
     IF ( ALLOCATED(v_in) )        DEALLOCATE(v_in)
     IF ( ALLOCATED(t_in) )        DEALLOCATE(t_in)
     IF ( ALLOCATED(qv_in) )       DEALLOCATE(qv_in)
     IF ( ALLOCATED(ps_in) )       DEALLOCATE(ps_in)
     IF ( ALLOCATED(rain_gsp_in) ) DEALLOCATE(rain_gsp_in)
     IF ( ALLOCATED(rain_con_in) ) DEALLOCATE(rain_con_in)
     IF ( ALLOCATED(snow_gsp_in) ) DEALLOCATE(snow_gsp_in)
     IF ( ALLOCATED(snow_con_in) ) DEALLOCATE(snow_con_in)
!DL also for averaged fields a...
     IF ( ALLOCATED(asobs2_in) )    DEALLOCATE(asobs2_in)
     IF ( ALLOCATED(athbs2_in) )    DEALLOCATE(athbs2_in)
     IF ( ALLOCATED(apabs_in) )     DEALLOCATE(apabs_in)
!DL
     IF ( ALLOCATED(sobs2_in) )    DEALLOCATE(sobs2_in)
     IF ( ALLOCATED(thbs2_in) )    DEALLOCATE(thbs2_in)
     IF ( ALLOCATED(pabs_in) )     DEALLOCATE(pabs_in)
     IF ( ALLOCATED(t_g_in) )      DEALLOCATE(t_g_in)
     IF ( ALLOCATED(albrad_in) )   DEALLOCATE(albrad_in)
     IF ( ALLOCATED(w_snow_in) )   DEALLOCATE(w_snow_in)
     IF ( ALLOCATED(t_snow_in) )   DEALLOCATE(t_snow_in)
     IF ( ALLOCATED(t_s_in) )      DEALLOCATE(t_s_in)
!DL
     IF ( ALLOCATED(t_so_0_in) )   DEALLOCATE(t_so_0_in)
     ALLOCATE (                       &
               rlon_in(ie_in,je_in),  & 
               rlat_in(ie_in,je_in),  &
               mask1(ie_in,je_in) )
     rlon_in=rundef
     rlat_in=rundef
     IF (ntype_atminput<=3) THEN
        ALLOCATE (               &
             u_in(ie_in,je_in),  &
             v_in(ie_in,je_in),  &
             t_in(ie_in,je_in),  &
             qv_in(ie_in,je_in), &
             ps_in(ie_in,je_in) )
        u_in=rundef
        v_in=rundef
        t_in=rundef
        qv_in=rundef
        ps_in=rundef
     ENDIF
     IF (ntype_raininput==1) THEN 
        ALLOCATE (                     &
             rain_gsp_in(ie_in,je_in), &
             rain_con_in(ie_in,je_in), &
             snow_gsp_in(ie_in,je_in), &
             snow_con_in(ie_in,je_in) )
        rain_gsp_in=rundef
        rain_con_in=rundef
        snow_gsp_in=rundef
        snow_con_in=rundef
     ENDIF
     IF (ntype_radinput<=2) THEN 
        ALLOCATE (                  &
             sobs2_in(ie_in,je_in), &
             thbs2_in(ie_in,je_in), &
!DL also for averaged fields
             asobs2_in(ie_in,je_in), &
             athbs2_in(ie_in,je_in), &
             apabs_in(ie_in,je_in),  &
!DL
             pabs_in(ie_in,je_in),  &
             t_g_in(ie_in,je_in),   &
             albrad_in(ie_in,je_in) )
        sobs2_in=rundef
        thbs2_in=rundef
        asobs2_in=rundef
        athbs2_in=rundef
        pabs_in=rundef
        apabs_in=rundef
        t_g_in=rundef
        albrad_in=rundef
     ENDIF
     IF (ntype_radinput==2) THEN 
        ALLOCATE (                   &
             w_snow_in(ie_in,je_in), &
             t_snow_in(ie_in,je_in), &
             t_s_in(ie_in,je_in)   , &
!DL
             t_so_0_in(ie_in,je_in) )
        w_snow_in=rundef
        t_snow_in=rundef
        t_s_in=rundef
        t_so_0_in=rundef
     ENDIF

     ! Reading gribfile
!! DL
     IF (yform_read == 'apix') THEN

     CALL init_griblist(99,ie_in,je_in,1)
     SELECT CASE (ntype_atminput)
     CASE (1)
        IF (ke_model>0) THEN
           CALL add2griblist(u_in,     'U' ,'generalVerticalLayer',ke_model,ke_model+1)   ! zonal wind
           CALL add2griblist(v_in,     'V' ,'generalVerticalLayer',ke_model,ke_model+1)   ! meridional wind
           CALL add2griblist(t_in,     'T' ,'generalVerticalLayer',ke_model,ke_model+1)   ! temperature
           CALL add2griblist(qv_in,    'QV','generalVerticalLayer',ke_model,ke_model+1)   ! specific humidity
        ELSE
           CALL add2griblist(u_in,     'U' ,'generalVerticalLayer')   ! zonal wind
           CALL add2griblist(v_in,     'V' ,'generalVerticalLayer')   ! meridional wind
           CALL add2griblist(t_in,     'T' ,'generalVerticalLayer')   ! temperature
           CALL add2griblist(qv_in,    'QV','generalVerticalLayer')   ! specific humidity
        ENDIF
        CALL add2griblist(ps_in,       'PS','surface')   ! pressure at the surface
     CASE (2)
        CALL add2griblist(u_in,     'SP_10M','heightAboveGround',0,10,'nudging')! wind speed
        CALL add2griblist(t_in,     'T_2M','heightAboveGround',0,2,'nudging') ! temperature
        CALL add2griblist(qv_in,    'TD_2M','heightAboveGround',0,2,'nudging') ! dew point
        CALL add2griblist(ps_in,    'PS','surface')   ! pressure at the surface
     CASE (3)
        CALL add2griblist(u_in,     'SP_10M','heightAboveGround',0,10,'analysis') ! wind speed
        CALL add2griblist(t_in,     'T_2M','heightAboveGround',0,2,'analysis')  ! temperature
        CALL add2griblist(qv_in,    'RELHUM_2M','heightAboveGround',0,2,'analysis')  ! relative humidity (!)
        CALL add2griblist(ps_in,    'PS','surface')    ! pressure at the surface
     END SELECT

     IF (ntype_raininput==1) THEN
        CALL add2griblist(rain_gsp_in,'RAIN_GSP','surface',-1,-1,'accum') ! stratiform rain
        CALL add2griblist(rain_con_in,'RAIN_CON','surface',-1,-1,'accum') ! convective rain
        CALL add2griblist(snow_gsp_in,'SNOW_GSP','surface',-1,-1,'accum') ! stratiform snow
        CALL add2griblist(snow_con_in,'SNOW_CON','surface',-1,-1,'accum') ! convective snow
     ENDIF

     SELECT CASE (ntype_radinput)
     CASE (1)
        CALL add2griblist(sobs2_in,   'SOBS_RAD','surface') ! net shortwave radiation
        CALL add2griblist(thbs2_in,   'THBS_RAD','surface') ! net longwave radiation
!DL also averaged fields
        CALL add2griblist(asobs2_in,   'ASOB_S','surface',-1,-1,'avg') ! net shortwave radiation
        CALL add2griblist(athbs2_in,   'ATHB_S','surface',-1,-1,'avg') ! net longwave radiation
!DL avg
        CALL add2griblist(albrad_in,  'ALB_RAD','surface') ! albeDO
!!DL OR CALL add2griblist(albrad_in,  'ALBEDO_B','surface')
        CALL add2griblist(t_g_in,     'T_G','surface') ! surface temperature
     CASE (2)
        CALL add2griblist(sobs2_in,   'SOBS_RAD','surface') ! net shortwave radiation
        CALL add2griblist(thbs2_in,   'THBS_RAD','surface') ! net longwave radiation
!DL also averaged fields
        CALL add2griblist(asobs2_in,   'ASOB_S','surface',-1,-1,'avg') ! net shortwave radiation
        CALL add2griblist(athbs2_in,   'ATHB_S','surface',-1,-1,'avg') ! net longwave radiation
!DL avg
        CALL add2griblist(albrad_in,  'ALB_RAD','surface') ! albeDO
!!DL OR CALL add2griblist(albrad_in,  'ALBEDO_B','surface')
        CALL add2griblist(t_s_in,     'T_S','depthBelowLand',0,0,'analysis')   ! surface temperature
!DL also trying T_SO(0)
        CALL add2griblist(t_so_0_in,  'T_SO','depthBelowLand',0,0,'analysis')
        CALL add2griblist(t_snow_in,  'T_SNOW','surface') ! snow temperature
        CALL add2griblist(w_snow_in,  'W_SNOW','surface') ! water content of snow storage
     END SELECT

     IF (lpar) THEN
        CALL add2griblist(pabs_in,    'PABS_RAD','surface') ! PAR
!DL also averaged field
        CALL add2griblist(apabs_in,   'APAB_S','surface',-1,-1,'avg') ! averaged PAR
!DL avg
     ENDIF

     ! Try to read data from file
     WRITE(6,*) 'Reading ',TRIM(filename)
     CALL read_grib_eccodes (TRIM(filename),ivv=nave,anyundef=any_undef,lhourly=lhourly_in)
     CALL clear_griblist()

     ELSE
 
     CALL init_griblist(99,ie_in,je_in,1)  
     SELECT CASE (ntype_atminput) 
     CASE (1)
        IF (ke_model>0) THEN
           CALL add2griblist(u_in,     33,  2,110,ke_model,ke_model+1)   ! zonal wind
           CALL add2griblist(v_in,     34,  2,110,ke_model,ke_model+1)   ! meridional wind
           CALL add2griblist(t_in,     11,  2,110,ke_model,ke_model+1)   ! temperature
           CALL add2griblist(qv_in,    51,  2,110,ke_model,ke_model+1)   ! specific humidity
        ELSE
           CALL add2griblist(u_in,     33,  2,110)   ! zonal wind
           CALL add2griblist(v_in,     34,  2,110)   ! meridional wind
           CALL add2griblist(t_in,     11,  2,110)   ! temperature
           CALL add2griblist(qv_in,    51,  2,110)   ! specific humidity
        ENDIF
        CALL add2griblist(ps_in,     1,  2,  1)   ! pressure at the surface
     CASE (2)
        CALL add2griblist(u_in,     32,  2,105,0,10,13)! wind speed
        CALL add2griblist(t_in,     11,  2,105,0,2,13) ! temperature
        CALL add2griblist(qv_in,    17,  2,105,0,2,13) ! dew point
        CALL add2griblist(ps_in,     1,  2,  1)   ! pressure at the surface
     CASE (3)
        CALL add2griblist(u_in,     32,  2,105,0,10,0) ! wind speed
        CALL add2griblist(t_in,     11,  2,105,0,2,0)  ! temperature
        CALL add2griblist(qv_in,    52,  2,105,0,2,0)  ! relative humidity (!)
        CALL add2griblist(ps_in,     1,  2,  1)    ! pressure at the surface
     END SELECT
 
     IF (ntype_raininput==1) THEN 
        CALL add2griblist(rain_gsp_in,102,201,  1) ! stratiform rain
        CALL add2griblist(rain_con_in,113,201,  1) ! convective rain
        CALL add2griblist(snow_gsp_in, 79,  2,  1) ! stratiform snow
        CALL add2griblist(snow_con_in, 78,  2,  1) ! convective snow
     ENDIF
 
     SELECT CASE (ntype_radinput)  
     CASE (1)
        CALL add2griblist(sobs2_in,   111,  2,  1) ! net shortwave radiation
        CALL add2griblist(thbs2_in,   112,  2,  1) ! net longwave radiation
        CALL add2griblist(albrad_in,   84,  2,  1) ! albeDO
        CALL add2griblist(t_g_in,      11,  2,  1) ! surface temperature
     CASE (2)
        CALL add2griblist(sobs2_in,   111,  2,  1) ! net shortwave radiation
        CALL add2griblist(thbs2_in,   112,  2,  1) ! net longwave radiation
        CALL add2griblist(albrad_in,   84,  2,  1) ! albeDO
        CALL add2griblist(t_s_in,      85,  2,111,0,0,0)   ! surface temperature
        CALL add2griblist(t_snow_in,  203,201,  1) ! snow temperature
        CALL add2griblist(w_snow_in,   65,  2,  1) ! water content of snow storage
     END SELECT
 
     IF (lpar) THEN
        CALL add2griblist(pabs_in,      5,201,  1) ! PAR
     ENDIF
 
     ! Try to read data from file
     WRITE(6,*) 'Reading ',TRIM(filename)
     CALL read_grib(TRIM(filename),ipds18=nave,anyundef=any_undef)  
     CALL clear_griblist()
 
     ENDIF ! yform_read

!DL####### NEW ##### check if all required fields are available, otherwise STOP #####
     imiss=0
     imiss_tg=0
     IF (ALL(u_in==rundef)) THEN
       imiss = imiss + 1
       WRITE(6,*) "U is missing!"
!#   ELSEIF ((ntype_atminput==1) .AND. (ALL(v_in==rundef))) THEN
!#     imiss = imiss + 1
!#     WRITE(6,*) "V is missing!"
     ELSEIF (ALL(t_in==rundef)) THEN
       imiss = imiss + 1
       WRITE(6,*) "T is missing!"
     ELSEIF (ALL(qv_in==rundef)) THEN
       imiss = imiss + 1
       WRITE(6,*) "QV is missing!"
     ELSEIF (ALL(ps_in==rundef)) THEN
       imiss = imiss + 1
       WRITE(6,*) "PS is missing!"
     ENDIF
     IF (ntype_atminput==1) THEN
       IF (ALL(v_in==rundef)) THEN
         imiss = imiss + 1
         WRITE(6,*) "V is missing!"
       ENDIF
     ENDIF
     IF (ntype_radinput==1) THEN
       IF (ALL(t_g_in==rundef)) THEN
         imiss_tg = 1
         WRITE(6,*) "T_G is missing!"
       ENDIF
     ENDIF
     IF (ntype_radinput==2) THEN
       IF (ALL(t_s_in==rundef)) THEN
         IF (ALL(t_so_0_in==rundef)) THEN
           imiss_tg = imiss_tg + 1
           WRITE(6,*) "T_S is missing!"
         ELSE
           t_s_in=t_so_0_in
           WRITE(6,*) "T_SO(0) taken as T_S!"
         ENDIF
       ELSEIF (ALL(t_snow_in==rundef)) THEN
         imiss_tg = imiss_tg + 1
         WRITE(6,*) "T_SNOW is missing!"
       ELSEIF (ALL(w_snow_in==rundef)) THEN
         imiss_tg = imiss_tg + 1
         WRITE(6,*) "W_SNOW is missing!"
       ENDIF
     ENDIF
     If (imiss /= 0) THEN
       WRITE(6,*) " !!! There are ",imiss," required fields missing in read_lmgrib (met. forcing) !!!"
       STOP 'read_lmgrib'
     ENDIF
     IF ((imiss_tg==1) .AND. (ntype_radinput==1)) THEN
       WRITE(6,*) " !!! STOP because T_G is missing in read_lmgrib !!!"
       STOP 'read_lmgrib'
     ENDIF
     IF ((imiss_tg/=0) .AND. (ntype_radinput==2)) THEN
       WRITE(6,*) " !!! There are ",imiss_tg," required fields missing to compute T_G in read_lmgrib !!!"
       STOP 'read_lmgrib'
     ENDIF
!DL####### NEW ##### check if all required fields are available, otherwise STOP #####

     ! YUCHKDAT
     IF (lcheck) THEN
        OPEN(UNIT=444,FILE='YUCHKDAT',FORM='FORMATTED',STATUS='OLD',POSITION='APPEND')
        WRITE(444,*)
        WRITE(444,*) 'Reading ',TRIM(filename)
        CALL stats("u",u_in,llandmask,UNIT=444)
        IF (ntype_atminput==1) THEN
           CALL stats("v",v_in,llandmask,UNIT=444)
        ENDIF
        CALL stats("t",t_in,llandmask,UNIT=444)
        CALL stats("qv",qv_in,llandmask,UNIT=444)
        CALL stats("ps",ps_in,llandmask,UNIT=444)
        IF (ntype_raininput==1) THEN 
           CALL stats("rain_gsp",rain_gsp_in,llandmask,UNIT=444)
           CALL stats("rain_con",rain_con_in,llandmask,UNIT=444)
           CALL stats("snow_gsp",snow_gsp_in,llandmask,UNIT=444)
           CALL stats("snow_con",snow_con_in,llandmask,UNIT=444)
        ENDIF
        CALL stats("sobs2",sobs2_in,llandmask,UNIT=444)
        CALL stats("thbs2",thbs2_in,llandmask,UNIT=444)
!DL
        CALL stats("asobs2",asobs2_in,llandmask,UNIT=444)
        CALL stats("athbs2",athbs2_in,llandmask,UNIT=444)
!DL
        IF (lpar) THEN
           CALL stats("pabs",pabs_in,llandmask,UNIT=444)
!DL
           CALL stats("apabs",apabs_in,llandmask,UNIT=444)
!DL
        ENDIF
        CALL stats("albrad",albrad_in,llandmask,UNIT=444)

        SELECT CASE (ntype_radinput)
        CASE (1)
           CALL stats("t_g",t_g_in,llandmask,UNIT=444)
        CASE (2)
           CALL stats("t_s",t_s_in,llandmask,UNIT=444)
           CALL stats("t_so_0",t_so_0_in,llandmask,UNIT=444)   !DL
           CALL stats("t_snow",t_snow_in,llandmask,UNIT=444)
           CALL stats("w_snow",w_snow_in,llandmask,UNIT=444)
        END SELECT
        CLOSE(UNIT=444)
     ENDIF

     ! check if this file is ok
     IF ( any_undef ) THEN
        IF ( .NOT. look_forward ) THEN
           WRITE(6,*) "File ",TRIM(filename)," contains undef values"
           STOP
        ENDIF
        IF ( tincr < tincr_max*60 ) THEN
           tincr = tincr + 60
           new_date = date_SUM(idate,tincr)
           WRITE(effective_date,'(i10.10)') new_date/100
           WRITE(6,*) "File ",TRIM(filename)," contains undef values, looking at ", TRIM(effective_date)
        ELSE
           WRITE(6,*) "Maximum gap reached without finding file (tincr_max= ", tincr_max, ")"
           STOP
        ENDIF
     ELSE
        EXIT
     ENDIF
  ENDDO  ! find undef if any

  ! conversions of input variables
  SELECT CASE (ntype_atminput) 
  CASE (1) ! destagger velocities, if requested
     IF (ldestaggeruv) THEN
        DO i=ie_in,1,-1
           DO j=je_in,1,-1
              u_in(i,j)=0.5_wp*(u_in(MAX(1,i-1),j)+u_in(i,j))
              v_in(i,j)=0.5_wp*(v_in(i,MAX(1,j-1))+v_in(i,j))
           ENDDO
        ENDDO
     ENDIF
     u_in(:,:)=SQRT(u_in(:,:)*u_in(:,:)+v_in(:,:)*v_in(:,:)) 
  CASE (2) ! conversion from dew point to specific humidity 
     DO i=1,ie_in
        DO j=1,je_in
           qv_in(i,j)= fqvs ( fpvsw ( qv_in(i,j) ) , ps_in (i,j) )
        ENDDO
     ENDDO
  CASE (3) ! conversion from relative to specific humidity 
     DO i=1,ie_in
        DO j=1,je_in
           qv_in(i,j)= MIN( 1.0_wp, qv_in(i,j) ) * &
                fqvs ( fpvsw ( t_in(i,j) ) , ps_in (i,j) )
        ENDDO
     ENDDO
  END SELECT

  IF (ntype_radinput==2) THEN ! determine t_g
     t_g_in= t_snow_in                                  &
             + (1.0_wp - MIN(1.0_wp,w_snow_in/cf_snow)) &
             * (t_s_in - t_snow_in)
  ENDIF

  IF (ntype_radinput<=2) THEN ! convert net radiation to downwelling radiation
!DL First, check missing fields !!
    IF ( ALL(albrad_in==rundef) ) THEN
      albrad_in = csalb_p
      WRITE(6,*) "ALB_RAD is missing, therefore set to csalb_p!!"
    ENDIF

!!DL    IF ( ALL(sobs2_in==rundef) ) THEN
!!DL       IF ( ALL(asobs2_in==rundef) ) THEN
!!DL         sobs2_in = 0.0
!!DL         WRITE(6,*) "sobs2_in/asobs2_in is missing, therefore set to zero!!"
!!DL       ELSE
!!DL         sobs2_in = asobs2_in
!!DL         WRITE(6,*) "asob_s is taken as sobs_rad!"
!!DL       ENDIF
!!DL    ENDIF
!!DL    sobs2_in(:,:) = sobs2_in(:,:) /(1.0_wp-albrad_in(:,:)/100.0_wp)

!DL Take ASOB_S (averaged field), if missing take SOBS_RAD (instant)
    IF ( ANY(asobs2_in/=rundef) ) THEN
      sobs2_in = asobs2_in
    ELSEIF ( ALL(asobs2_in==rundef) ) THEN
      IF ( ALL(sobs2_in==rundef) ) THEN
        sobs2_in = 0.0
        WRITE(6,*) "SOBS_RAD/ASOB_S is missing, therefore set to zero!!"
      ELSE
!DL     sob2_in already contains sobs_rad
        WRITE(6,*) "SOBS_RAD is taken as ASOB_S!"
      ENDIF
    ENDIF
    sobs2_in(:,:) = sobs2_in(:,:) /(1.0_wp-albrad_in(:,:)/100.0_wp)
    

!!DL    IF ( ALL(thbs2_in==rundef) ) THEN
!!DL      IF ( ALL(athbs2_in==rundef) ) THEN
!!DL        thbs2_in = 0.0
!!DL        WRITE(6,*) "thbs2_in is missing, therefore set to zero!!"
!!DL      ELSE
!!DL        thbs2_in = athbs2_in
!!DL        WRITE(6,*) "athb_s is taken as thbs_rad!"
!!DL        thbs2_in(:,:) = (thbs2_in(:,:)+(1.0_wp-ctalb)*sigma*t_g_in(:,:)**4)/(1.0_wp-ctalb)
!!DL      ENDIF
!!DL    ELSE
!!DL       thbs2_in(:,:) = (thbs2_in(:,:)+(1.0_wp-ctalb)*sigma*t_g_in(:,:)**4)/(1.0_wp-ctalb)
!!DL    ENDIF

!DL Take ATHB_S (averaged field), if missing take THBS_RAD (instant)
    IF ( ANY(athbs2_in/=rundef) ) THEN
      thbs2_in = athbs2_in
    ELSEIF ( ALL(athbs2_in==rundef) ) THEN
      IF ( ALL(thbs2_in==rundef) ) THEN
        thbs2_in = 0.0
        WRITE(6,*) "THBS_RAD/ATHB_S is missing, therefore set to zero!!"
      ELSE
!DL     thbs2_in already contains thbs_rad
        WRITE(6,*) "THBS_RAD taken as ATHB_S!"
      ENDIF
    ENDIF
    thbs2_in(:,:) = (thbs2_in(:,:)+(1.0_wp-ctalb)*sigma*t_g_in(:,:)**4)/(1.0_wp-ctalb)

!!DL sobs2_in(:,:) = sobs2_in(:,:) /(1.0_wp-albrad_in(:,:)/100.0_wp)
!!DL thbs2_in(:,:) = (thbs2_in(:,:)+(1.0_wp-ctalb)*sigma*t_g_in(:,:)**4)/(1.0_wp-ctalb)

  ENDIF  ! ntype_radinput<=2

!!DL also for missing pabs_in
!DL Take APAB_S (averaged field), if missing take PABS_RAD (instant)
!DL if missing or .NOT.lpar take 0.5*sobs2_in
    IF ( ANY(apabs_in/=rundef) ) THEN
      pabs_in = apabs_in
    ELSEIF ( ALL(apabs_in==rundef) ) THEN
      IF ( ALL(pabs_in==rundef) ) THEN
        pabs_in = 0.5_wp*sobs2_in
        WRITE(6,*) "PABS_RAD/APAB_S is missing, therefore set to 0.5*sob2_in!!"
      ELSE
!DL     pabs_in already contains pabs_rad
        WRITE(6,*) "PABS_RAD taken as APAB_S!"
      ENDIF
    ENDIF

  IF (ntype_raininput==1) THEN
!DL First, check missing fields !!
    IF ( ALL(rain_gsp_in==rundef) ) THEN 
      rain_gsp_in = 0.0
      WRITE(6,*) "RAIN_GSP is missing, set to zero !"
    ENDIF
    IF ( ALL(rain_con_in==rundef) ) THEN 
      WRITE(6,*) "RAIN_CON is missing, set to zero !"
      rain_con_in = 0.0
    ENDIF
    IF ( ALL(snow_gsp_in==rundef) ) THEN 
      WRITE(6,*) "SNOW_GSP is missing, set to zero !"
      snow_gsp_in = 0.0
    ENDIF
    IF ( ALL(snow_con_in==rundef) ) THEN 
      WRITE(6,*) "SNOW_CON is missing, set to zero !"
      snow_con_in = 0.0
    ENDIF
    rain_gsp_in=rain_gsp_in/3600.0
    rain_con_in=rain_con_in/3600.0
    snow_gsp_in=snow_gsp_in/3600.0
    snow_con_in=snow_con_in/3600.0
  ENDIF

!DL
  ! check if precipitation / radiation are hourly data
  IF (lhourly_in .NEQV. lhourly_data) THEN
    lhourly_data = lhourly_in
    WRITE(6,*) "lhourly_data changed to ",lhourly_in
  ENDIF
!DL

  ! check domains and if simple croping is possible
  lcrop=.FALSE.
  IF ((ABS(pollat-pollat_in)<1.0E-3_wp)      &
     .AND. (ABS(pollon-pollon_in)<1.0E-3_wp) &
     .AND. (ABS(dlon_in-dlon)<1.0E-3_wp)     &
     .AND. (ABS(dlat_in-dlat)<1.0E-3_wp))    THEN
     i0_crop=NINT((startlon-rlon0_in)/dlon+1)
     j0_crop=NINT((startlat-rlat0_in)/dlat+1)
     IF ((i0_crop>=1)                                              &
        .AND. (j0_crop>=1)                                         &
        .AND. (ABS((i0_crop-1)*dlon-startlon+rlon0_in)<1.0E-3_wp)  &
        .AND. (ABS((j0_crop-1)*dlat-startlat+rlat0_in)<1.0E-3_wp)) THEN
        lcrop=.TRUE.
     ENDIF
  ENDIF

  ! tell user what we are doing with forcing data
  IF (lcrop) THEN
     IF (i0_crop==1 .AND. j0_crop==1 .AND. ie==ie_in .AND. je==je_in) THEN
        WRITE(6,*) "Perfect match of forcing data!"
     ELSE
        WRITE(6,*) "Croping forcing data!"
     ENDIF
  ELSE
     WRITE(6,*) "Interpolating forcing data!"
  ENDIF

  ! Calculate coordinates of each input grid point
  IF (.NOT.lcrop) THEN 
     IF ((pollat==pollat_in).AND.(pollon==pollon_in)) THEN
        DO i=1,ie_in
           rlon_in(i,:)=rlon0_in+dlon_in*(i-1)
        ENDDO
        DO j=1,je_in  
           rlat_in(:,j)=rlat0_in+dlat_in*(j-1)
        ENDDO
     ELSE
        DO i=1,ie_in 
           rlon_rot_in=(i-1)*dlon_in+rlon0_in
           DO j=1,je_in
              rlat_rot_in=(j-1)*dlat_in+rlat0_in
              ! XYZ> added polgam=0 for compatability with utilities
              rlat_geo=phirot2phi(rlat_rot_in,rlon_rot_in,pollat_in,pollon_in,0.0_wp)
              rlon_geo=rlarot2rla(rlat_rot_in,rlon_rot_in,pollat_in,pollon_in,0.0_wp)
              rlat_in(i,j)=phi2phirot(rlat_geo,rlon_geo,pollat,pollon)
              rlon_in(i,j)=rla2rlarot(rlat_geo,rlon_geo,pollat,pollon,0.0_wp)
           ENDDO
        ENDDO
     ENDIF
  ENDIF

  ! Spatial Interpolation
  mask1=.TRUE.
  fac=dlat_in/dlat 
  IF (ntype_atminput<=3) THEN
     IF (lcrop) THEN
        u_bd(:,:,ke,nx) =u_in (i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)
        t_bd(:,:,ke,nx) =t_in (i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)
        qv_bd(:,:,ke,nx)=qv_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)
!       qv_bd(:,:,ke)=qv_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)    !XYZ> in TSA, qv is 4 dimensional
        ps_bd(:,:,nx)   =ps_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)
     ELSE
        CALL general_interpol2d(rlat_in,rlon_in,u_in,      &
             startlat,startlon,dlat,dlon,u_bd(:,:,ke,nx),  &
             2,2,0.51_wp*fac,1.0_wp)                          
        CALL general_interpol2d(rlat_in,rlon_in,t_in,      &
             startlat,startlon,dlat,dlon,t_bd(:,:,ke,nx),  &
             2,2,0.51_wp*fac,1.0_wp,mask1,llandmask)
        CALL general_interpol2d(rlat_in,rlon_in,qv_in,     &
             startlat,startlon,dlat,dlon,qv_bd(:,:,ke,nx), &
             2,2,0.51_wp*fac,1.0_wp,mask1,llandmask)
!       CALL general_interpol2d(rlat_in,rlon_in,qv_in,     &    !XYZ> in TSA, qv is 4 dimensional
!            startlat,startlon,dlat,dlon,qv_bd(:,:,ke),    &
!            2,2,0.51_wp*fac,1.0_wp,mask1,llandmask)
        CALL general_interpol2d(rlat_in,rlon_in,ps_in,     &
             startlat,startlon,dlat,dlon,ps_bd(:,:,nx),    &
             2,2,0.51_wp*fac,1.0_wp,mask1,llandmask)
     ENDIF
  ENDIF
  IF (ntype_raininput==1) THEN
     IF (lcrop) THEN
        prr_gsp_bd(:,:,nx) =(rain_gsp_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)+ &
                             rain_con_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1))
        prs_gsp_bd(:,:,nx) =(snow_gsp_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)+ &
                             snow_con_in(i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1))
     ELSE
        CALL general_interpol2d(rlat_in,rlon_in,             &
             (rain_gsp_in+rain_con_in),                      &
             startlat,startlon,dlat,dlon,prr_gsp_bd(:,:,nx), &
             2,2,0.51_wp*fac,1.0_wp,mask1,llandmask)
        CALL general_interpol2d(rlat_in,rlon_in,             &
             (snow_gsp_in+snow_con_in),                      &
             startlat,startlon,dlat,dlon,prs_gsp_bd(:,:,nx), &
             2,2,0.51_wp*fac,1.0_wp,mask1,llandmask)
     ENDIF
  ENDIF
  IF (ntype_radinput<=2) THEN
     IF (lcrop) THEN
        so_down_bd(:,:,nx) =sobs2_in (i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)
        th_down_bd(:,:,nx) =thbs2_in (i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)
        pabs_bd   (:,:,nx) = pabs_in (i0_crop:i0_crop+ie-1,j0_crop:j0_crop+je-1)
     ELSE
        CALL general_interpol2d(rlat_in,rlon_in,sobs2_in,    &
             startlat,startlon,dlat,dlon,so_down_bd(:,:,nx), &
             2,2,0.51_wp*fac,1.0_wp,mask1,llandmask)
        CALL general_interpol2d(rlat_in,rlon_in,thbs2_in,    &
             startlat,startlon,dlat,dlon,th_down_bd(:,:,nx), &
             2,2,0.51_wp*fac,1.0_wp,mask1,llandmask)
        CALL general_interpol2d(rlat_in,rlon_in,pabs_in,     &
             startlat,startlon,dlat,dlon,pabs_bd(:,:,nx),    &
             2,2,0.51_wp*fac,1.0_wp,mask1,llandmask)      
     ENDIF
  ENDIF

  ! Clean up
  DEALLOCATE(rlon_in,rlat_in,mask1)
  IF (ntype_atminput<=3) THEN
     DEALLOCATE(u_in, v_in,t_in,qv_in,ps_in )
  ENDIF
  IF (ntype_raininput==1) THEN 
     DEALLOCATE(rain_gsp_in,rain_con_in, snow_gsp_in, snow_con_in)
  ENDIF
  IF (ntype_radinput<=2) THEN 
!DL  DEALLOCATE(sobs2_in, thbs2_in, pabs_in, t_g_in,albrad_in)
     DEALLOCATE(sobs2_in, asobs2_in, thbs2_in, athbs2_in, pabs_in, apabs_in, t_g_in,albrad_in)
  ENDIF
  IF (ntype_radinput==2) THEN
!DL  DEALLOCATE(w_snow_in,t_snow_in,t_s_in)
     DEALLOCATE(w_snow_in,t_snow_in,t_s_in,t_so_0_in)
  ENDIF

  ! Gap filling
  IF ((ntype_atminput<=3).AND.(ntype_raininput==1).AND.(ntype_radinput<=2)) THEN
     IF (COUNT((t_bd(:,:,ke,nx)==rundef).AND.llandmask(:,:))>0) THEN
        WRITE(6,*) "Performing gap filling in read_lmgrib ..."   
        ALLOCATE( gfill(ie,je,9) )
        gfill(:,:,1)=u_bd(:,:,ke,nx)
        gfill(:,:,2)=t_bd(:,:,ke,nx)
        gfill(:,:,3)=qv_bd(:,:,ke,nx)
!       gfill(:,:,3)=qv_bd(:,:,ke)   !XYZ> in TSA, qv is 4 dimensional
        gfill(:,:,4)=ps_bd(:,:,nx)
        gfill(:,:,5)=prr_gsp_bd(:,:,nx)
        gfill(:,:,6)=prs_gsp_bd(:,:,nx)
        gfill(:,:,7)=so_down_bd(:,:,nx)
        gfill(:,:,8)=th_down_bd(:,:,nx)
        gfill(:,:,9)=pabs_bd(:,:,nx)

        DO i=1,ie
           DO j=1,je
              IF ((t_bd(i,j,ke,nx)==rundef).AND.llandmask(i,j)) THEN
                 DO ir=1,MAX(i,j,ie-i,je-j)
                    imin=MAX(1,i-ir)
                    jmin=MAX(1,j-ir)
                    imax=MIN(ie,i+ir)
                    jmax=MIN(je,j+ir)
                    cvalid=COUNT((gfill(imin:imax,jmin:jmax,2)>rundef).AND.   &
                         (llandmask(imin:imax,jmin:jmax)))
                    IF (cvalid>0) THEN
                       u_bd(i,j,ke,nx) = SUM(gfill(imin:imax,jmin:jmax,1),    &
                            ((gfill(imin:imax,jmin:jmax,1)>rundef)            &
                            .AND. (llandmask(imin:imax,jmin:jmax)))) /        &
                            FLOAT(COUNT((gfill(imin:imax,jmin:jmax,1)>rundef) &
                            .AND. (llandmask(imin:imax,jmin:jmax))))
                       t_bd(i,j,ke,nx) = SUM(gfill(imin:imax,jmin:jmax,2),    &
                            ((gfill(imin:imax,jmin:jmax,2)>rundef)            &
                            .AND. (llandmask(imin:imax,jmin:jmax)))) /        &
                            FLOAT(COUNT((gfill(imin:imax,jmin:jmax,2)>rundef) &
                            .AND. (llandmask(imin:imax,jmin:jmax))))
                       qv_bd(i,j,ke,nx) = SUM(gfill(imin:imax,jmin:jmax,3),   &
                            ((gfill(imin:imax,jmin:jmax,3)>rundef)            &
                            .AND. (llandmask(imin:imax,jmin:jmax)))) /        &
                            FLOAT(COUNT((gfill(imin:imax,jmin:jmax,3)>rundef) &
                            .AND. (llandmask(imin:imax,jmin:jmax))))
!XYZ> in TSA, qv is 4 dimensional
!                      qv_bd(i,j,ke) = SUM(gfill(imin:imax,jmin:jmax,3),      &
!                           ((gfill(imin:imax,jmin:jmax,3)>rundef)            &
!                           .AND. (llandmask(imin:imax,jmin:jmax)))) /        & 
!                          FLOAT(COUNT((gfill(imin:imax,jmin:jmax,3)>rundef)  &
!                           .AND. (llandmask(imin:imax,jmin:jmax))))
!XYZ<
                       ps_bd(i,j,nx) = SUM(gfill(imin:imax,jmin:jmax,4),      &
                            ((gfill(imin:imax,jmin:jmax,4)>rundef)            &
                            .AND. (llandmask(imin:imax,jmin:jmax)))) /        &
                            FLOAT(COUNT((gfill(imin:imax,jmin:jmax,4)>rundef) &
                            .AND. (llandmask(imin:imax,jmin:jmax))))
                       prr_gsp_bd(i,j,nx) = SUM(gfill(imin:imax,jmin:jmax,5), &
                            ((gfill(imin:imax,jmin:jmax,5)>rundef)            &
                            .AND. (llandmask(imin:imax,jmin:jmax)))) /        &
                            FLOAT(COUNT((gfill(imin:imax,jmin:jmax,5)>rundef) &
                            .AND. (llandmask(imin:imax,jmin:jmax))))
                       prs_gsp_bd(i,j,nx) = SUM(gfill(imin:imax,jmin:jmax,6), &
                            ((gfill(imin:imax,jmin:jmax,6)>rundef)            &
                            .AND. (llandmask(imin:imax,jmin:jmax)))) /        &
                            FLOAT(COUNT((gfill(imin:imax,jmin:jmax,6)>rundef) &
                            .AND. (llandmask(imin:imax,jmin:jmax))))
                       so_down_bd(i,j,nx) = SUM(gfill(imin:imax,jmin:jmax,7), &
                            ((gfill(imin:imax,jmin:jmax,7)>rundef)            &
                            .AND. (llandmask(imin:imax,jmin:jmax)))) /        &
                            FLOAT(COUNT((gfill(imin:imax,jmin:jmax,7)>rundef) &
                            .AND. (llandmask(imin:imax,jmin:jmax))))
                       th_down_bd(i,j,nx) = SUM(gfill(imin:imax,jmin:jmax,8), &
                            ((gfill(imin:imax,jmin:jmax,8)>rundef)            &
                            .AND. (llandmask(imin:imax,jmin:jmax)))) /        &
                            FLOAT(COUNT((gfill(imin:imax,jmin:jmax,8)>rundef) &
                            .AND. (llandmask(imin:imax,jmin:jmax))))
                       pabs_bd(i,j,nx) = SUM(gfill(imin:imax,jmin:jmax,9),    &
                            ((gfill(imin:imax,jmin:jmax,9)>rundef)            &
                            .AND. (llandmask(imin:imax,jmin:jmax)))) /        &
                            FLOAT(COUNT((gfill(imin:imax,jmin:jmax,9)>rundef) &
                            .AND. (llandmask(imin:imax,jmin:jmax))))
                       EXIT
                    ENDIF  ! cvalid >0
                 ENDDO  ! ir
              ENDIF  ! t_bd==rundef...
           ENDDO  ! j
        ENDDO  ! i
        DEALLOCATE(gfill)
     ENDIF  ! count t_bd...
  ENDIF  ! ntype

  ! Control output
  IF (lcheck) THEN
     IF (ntype_atminput<=3) THEN
        CALL stats("u_bd",u_bd(:,:,ke,nx),llandmask)
        CALL stats("t_bd",t_bd(:,:,ke,nx),llandmask)
        CALL stats("qv_bd",qv_bd(:,:,ke,nx),llandmask)
!       CALL stats("qv_bd",qv_bd(:,:,ke),llandmask)   !XYZ> in TSA, qv is 4 dimensional
        CALL stats("ps_bd",ps_bd(:,:,nx),llandmask)
     ENDIF
     IF (ntype_raininput==1) THEN 
        CALL stats("prr_gsp_bd",prr_gsp_bd(:,:,nx),llandmask)
        CALL stats("prs_gsp_bd",prs_gsp_bd(:,:,nx),llandmask)
     ENDIF
     IF (ntype_radinput<=2) THEN 
        CALL stats("so_down_bd",so_down_bd(:,:,nx),llandmask)
        CALL stats("th_down_bd",th_down_bd(:,:,nx),llandmask)
        CALL stats("pabs_bd",pabs_bd(:,:,nx),llandmask)
     ENDIF
  ENDIF

END SUBROUTINE read_lmgrib
!----------------------------------------------------------------------


SUBROUTINE read_icongrib(date, nx, effective_date, nave, look_forward)
!----------------------------------------------------------------------
! Description:
!
! Read icon model grib files and assign fields to corresponding TERRA fields 
!----------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
  CHARACTER (LEN=10),       INTENT(IN)  :: date
  CHARACTER (LEN=10),       INTENT(OUT) :: effective_date
  INTEGER ,                 INTENT(IN)  :: nx
  INTEGER ,                 INTENT(OUT) :: nave
  LOGICAL ,                 INTENT(IN)  :: look_forward

! Local Variables
  CHARACTER (LEN=1)           :: &
       uuid_in(16)
  CHARACTER (LEN=100)         :: &
       filename                , & ! 
       filename1               , & ! 
       filename2                   ! 
  INTEGER                     :: &
       ie_in                   , & ! 
       je_in                   , & ! 
       i                       , & ! 
       j                       , & ! 
       tincr                       ! 
!DL
  INTEGER                     :: ngp
  INTEGER                     :: &
       imiss, imiss_tg
!DL
  INTEGER (KIND=kind_idate)   :: &
       idate                   , & ! 
       new_date                    ! 
  LOGICAL                     :: &
       is_found                , & ! 
!DL tot_prec available
       ltot_prec               , & !
       any_undef                   ! 
!!REAL (KIND=wp)              :: &
!!     rlat_geo                , & ! 
!!     rlon_geo                , & ! 
!!     rlat_rot_in             , & ! 
!!     rlon_rot_in             , & ! 
!!     rlat0_in                , & ! 
!!     rlon0_in                , & ! 
!!     dlat_in                 , & ! 
!!     dlon_in                 , & ! 
!!     pollon_in               , & ! 
!!     pollat_in                   ! 
  REAL (KIND=wp), ALLOCATABLE :: &
       u_in(:,:)               , & ! 
       v_in(:,:)               , & ! 
       t_in(:,:)               , & ! 
       qv_in(:,:)              , & ! 
       ps_in(:,:)              , & ! 
       rain_gsp_in(:,:)        , & ! 
       rain_con_in(:,:)        , & ! 
       snow_gsp_in(:,:)        , & ! 
       snow_con_in(:,:)        , & ! 
!!DL use also tot_prec
       tot_prec_in(:,:)        , & ! 
       sobs2_in(:,:)           , & ! 
!DL provide average fields as well asob2, athbs2, apabs
       asobs2_in(:,:)          , & ! 
       thbs2_in(:,:)           , & ! 
       athbs2_in(:,:)          , & ! 
       apabs_in(:,:)           , & ! 
       pabs_in(:,:)            , & ! 
       t_g_in(:,:)             , & ! 
       albrad_in(:,:)          , & ! 
       t_s_in(:,:)             , & ! 
!DL alternative to t_s_in
       t_so_0_in (:,:)         , & ! 
       w_snow_in(:,:)          , & ! 
       t_snow_in(:,:)          , & ! 
       gfill(:,:,:)                ! 

!Declaration of STATEMENT-FUNCTIONS
!--------------------
  REAL (KIND=wp) :: fpvsw, zst, fqvs, zsge, zsp
  fpvsw(zst)       = b1 * EXP( b2w * (zst - b3) / (zst - b4w) )
  fqvs (zsge, zsp) = rdv * zsge / ( zsp - o_m_rdv * zsge )


  ! Look for file, look forward when requested
  effective_date = date
  READ(date,*) idate
  idate = idate * 100
  tincr = 0

  DO       ! find undef if any 
     DO     ! find file
        filename1=TRIM(metfiledir)//TRIM(metfileprefix)//TRIM(effective_date)
        filename2=TRIM(filename1)//'00'    ! File name variant with minutes
        filename=filename1
        CALL get_grib_info_icon(filename,    &
                          ngrid,ngpi,ngp,uuid_in,is_found)
        IF ( .NOT. is_found ) THEN
          filename=filename2
          CALL get_grib_info_icon(filename,  &
                            ngrid,ngpi,ngp,uuid_in,is_found)
        ENDIF

        IF (ANY(uuid_in(:) /= uuid(:))) THEN
          WRITE(6,*) "uuid_in in read_icongrib does not match"
          STOP 'WRONG UUID'
        ENDIF

        IF ( is_found ) THEN
           EXIT
        ENDIF
        IF ( .NOT. look_forward ) THEN
           WRITE(6,*) "Neither ",TRIM(filename1)," nor ", TRIM(filename2), " have been found"
           STOP
        ELSE
print *, 'about tincr:  ', tincr, tincr_max, filename1, filename2
           IF ( tincr < tincr_max*60 ) THEN
              tincr = tincr + 60
              new_date = date_SUM(idate,tincr)
              WRITE(effective_date,'(i10.10)') new_date/100
              WRITE(6,*) "Neither ",TRIM(filename1)," nor ", TRIM(filename2), " have been found, ", &
                         "looking at ", TRIM(effective_date)
           ELSE
              WRITE(6,*) "Maximum gap reached without finding file (tincr_max= ", tincr_max, ")"
              STOP
           ENDIF
        ENDIF
     ENDDO  ! find file

     ! Allocation of fields depending on configuration defined by
     !  ntype_atminput, ntype_raininput, ntype_radinput
     IF ( ALLOCATED(u_in) )        DEALLOCATE(u_in)
     IF ( ALLOCATED(v_in) )        DEALLOCATE(v_in)
     IF ( ALLOCATED(t_in) )        DEALLOCATE(t_in)
     IF ( ALLOCATED(qv_in) )       DEALLOCATE(qv_in)
     IF ( ALLOCATED(ps_in) )       DEALLOCATE(ps_in)
     IF ( ALLOCATED(rain_gsp_in) ) DEALLOCATE(rain_gsp_in)
     IF ( ALLOCATED(rain_con_in) ) DEALLOCATE(rain_con_in)
     IF ( ALLOCATED(snow_gsp_in) ) DEALLOCATE(snow_gsp_in)
     IF ( ALLOCATED(snow_con_in) ) DEALLOCATE(snow_con_in)
!DL tot_prec as alternative
     IF ( ALLOCATED(tot_prec_in) ) DEALLOCATE(tot_prec_in)
!DL also for averaged fields a...
     IF ( ALLOCATED(asobs2_in) )    DEALLOCATE(asobs2_in)
     IF ( ALLOCATED(athbs2_in) )    DEALLOCATE(athbs2_in)
     IF ( ALLOCATED(apabs_in)  )    DEALLOCATE(athbs2_in)
!DL
     IF ( ALLOCATED(sobs2_in) )    DEALLOCATE(sobs2_in)
     IF ( ALLOCATED(thbs2_in) )    DEALLOCATE(thbs2_in)
     IF ( ALLOCATED(pabs_in) )     DEALLOCATE(pabs_in)
     IF ( ALLOCATED(t_g_in) )      DEALLOCATE(t_g_in)
     IF ( ALLOCATED(albrad_in) )   DEALLOCATE(albrad_in)
     IF ( ALLOCATED(w_snow_in) )   DEALLOCATE(w_snow_in)
     IF ( ALLOCATED(t_snow_in) )   DEALLOCATE(t_snow_in)
     IF ( ALLOCATED(t_s_in) )      DEALLOCATE(t_s_in)
!DL
     IF ( ALLOCATED(t_so_0_in) )   DEALLOCATE(t_so_0_in)

     ie_in=nproma
     je_in=nblock 

     IF (ntype_atminput<=3) THEN
        ALLOCATE (               &
             u_in(ie_in,je_in),  &
             v_in(ie_in,je_in),  &
             t_in(ie_in,je_in),  &
             qv_in(ie_in,je_in), &
             ps_in(ie_in,je_in) )
        u_in=rundef
        v_in=rundef
        t_in=rundef
        qv_in=rundef
        ps_in=rundef
     ENDIF
     IF (ntype_raininput==1) THEN 
        ALLOCATE (                     &
             rain_gsp_in(ie_in,je_in), &
             rain_con_in(ie_in,je_in), &
             snow_gsp_in(ie_in,je_in), &
             snow_con_in(ie_in,je_in), &
             tot_prec_in(ie_in,je_in)  )
        rain_gsp_in=rundef
        rain_con_in=rundef
        snow_gsp_in=rundef
        snow_con_in=rundef
!DL total precipitation
        tot_prec_in=rundef
     ENDIF
     IF (ntype_radinput<=2) THEN 
        ALLOCATE (                  &
             sobs2_in(ie_in,je_in), &
             thbs2_in(ie_in,je_in), &
!DL also for averaged fields
             asobs2_in(ie_in,je_in), &
             athbs2_in(ie_in,je_in), &
             apabs_in(ie_in,je_in), &
!DL
             pabs_in(ie_in,je_in),  &
             t_g_in(ie_in,je_in),   &
             albrad_in(ie_in,je_in) )
        sobs2_in=rundef
        thbs2_in=rundef
        asobs2_in=rundef
        athbs2_in=rundef
        pabs_in=rundef
        apabs_in=rundef
        t_g_in=rundef
        albrad_in=rundef
     ENDIF
     IF (ntype_radinput==2) THEN 
        ALLOCATE (                   &
             w_snow_in(ie_in,je_in), &
             t_snow_in(ie_in,je_in), &
             t_s_in(ie_in,je_in)   , &
!DL
             t_so_0_in(ie_in,je_in) )
        w_snow_in=rundef
        t_snow_in=rundef
        t_s_in=rundef
        t_so_0_in=rundef
     ENDIF

     ! Reading gribfile

     CALL init_griblist(99,ie_in,je_in,1)
     SELECT CASE (ntype_atminput)
     CASE (1)
        IF (ke_model>0) THEN
           CALL add2griblist(u_in,     'U' ,'generalVerticalLayer',ke_model,ke_model+1)   ! zonal wind
           CALL add2griblist(v_in,     'V' ,'generalVerticalLayer',ke_model,ke_model+1)   ! meridional wind
           CALL add2griblist(t_in,     'T' ,'generalVerticalLayer',ke_model,ke_model+1)   ! temperature
           CALL add2griblist(qv_in,    'QV','generalVerticalLayer',ke_model,ke_model+1)   ! specific humidity
        ELSE
           CALL add2griblist(u_in,     'U' ,'generalVerticalLayer')   ! zonal wind
           CALL add2griblist(v_in,     'V' ,'generalVerticalLayer')   ! meridional wind
           CALL add2griblist(t_in,     'T' ,'generalVerticalLayer')   ! temperature
           CALL add2griblist(qv_in,    'QV','generalVerticalLayer')   ! specific humidity
        ENDIF
        CALL add2griblist(ps_in,       'PS','surface')   ! pressure at the surface
     CASE (2)
        CALL add2griblist(u_in,     'SP_10M','heightAboveGround',0,10,'nudging')! wind speed
        CALL add2griblist(t_in,     'T_2M','heightAboveGround',0,2,'nudging') ! temperature
        CALL add2griblist(qv_in,    'TD_2M','heightAboveGround',0,2,'nudging') ! dew point
        CALL add2griblist(ps_in,    'PS','surface')   ! pressure at the surface
     CASE (3)
        CALL add2griblist(u_in,     'SP_10M','heightAboveGround',0,10,'analysis') ! wind speed
        CALL add2griblist(t_in,     'T_2M','heightAboveGround',0,2,'analysis')  ! temperature
        CALL add2griblist(qv_in,    'RELHUM_2M','heightAboveGround',0,2,'analysis')  ! relative humidity (!)
        CALL add2griblist(ps_in,    'PS','surface')    ! pressure at the surface
     END SELECT

     IF (ntype_raininput==1) THEN
        CALL add2griblist(rain_gsp_in,'RAIN_GSP','surface',-1,-1,'accum') ! stratiform rain
        CALL add2griblist(rain_con_in,'RAIN_CON','surface',-1,-1,'accum') ! convective rain
        CALL add2griblist(snow_gsp_in,'SNOW_GSP','surface',-1,-1,'accum') ! stratiform snow
        CALL add2griblist(snow_con_in,'SNOW_CON','surface',-1,-1,'accum') ! convective snow
!DL
        CALL add2griblist(tot_prec_in,'TOT_PREC','surface',-1,-1,'accum') ! total precipitation
     ENDIF

     SELECT CASE (ntype_radinput)
     CASE (1)
        CALL add2griblist(sobs2_in,   'SOBS_RAD','surface') ! net shortwave radiation
        CALL add2griblist(thbs2_in,   'THBS_RAD','surface') ! net longwave radiation
!DL also averaged fields
        CALL add2griblist(asobs2_in,   'ASOB_S','surface',-1,-1,'avg') ! net shortwave radiation
        CALL add2griblist(athbs2_in,   'ATHB_S','surface',-1,-1,'avg') ! net longwave radiation
!DL avg
        CALL add2griblist(albrad_in,  'ALB_RAD','surface') ! albeDO
!!DL OR CALL add2griblist(albrad_in,  'ALBEDO_B','surface')
        CALL add2griblist(t_g_in,     'T_G','surface') ! surface temperature
     CASE (2)
        CALL add2griblist(sobs2_in,   'SOBS_RAD','surface') ! net shortwave radiation
        CALL add2griblist(thbs2_in,   'THBS_RAD','surface') ! net longwave radiation
!DL also averaged fields
        CALL add2griblist(asobs2_in,   'ASOB_S','surface',-1,-1,'avg') ! net shortwave radiation
        CALL add2griblist(athbs2_in,   'ATHB_S','surface',-1,-1,'avg') ! net longwave radiation
!DL avg
        CALL add2griblist(albrad_in,  'ALB_RAD','surface') ! albeDO
!!DL OR CALL add2griblist(albrad_in,  'ALBEDO_B','surface')
        CALL add2griblist(t_s_in,     'T_S','depthBelowLand',0,0,'analysis')   ! surface temperature
!DL also trying T_SO(0)
        CALL add2griblist(t_so_0_in,  'T_SO','depthBelowLand',0,0,'analysis')
        CALL add2griblist(t_snow_in,  'T_SNOW','surface') ! snow temperature
        CALL add2griblist(w_snow_in,  'W_SNOW','surface') ! water content of snow storage
     END SELECT

     IF (lpar) THEN
        CALL add2griblist(pabs_in,     'PABS_RAD','surface') ! PAR
        CALL add2griblist(apabs_in,    'APAB_S','surface',-1,-1,'avg') ! PAR - average
     ENDIF

     ! Try to read data from file
     WRITE(6,*) 'Reading ',TRIM(filename)
     CALL read_grib_eccodes (TRIM(filename),ivv=nave,anyundef=any_undef)
     CALL clear_griblist()


!DL####### NEW ##### check if all required fields are available, otherwise STOP #####
     imiss=0
     imiss_tg=0
     IF (ALL(u_in==rundef)) THEN
       imiss = imiss + 1
       WRITE(6,*) "U is missing!"
!#   ELSEIF ((ntype_atminput==1) .AND. (ALL(v_in==rundef))) THEN
!#     imiss = imiss + 1
!#     WRITE(6,*) "V is missing!"
     ELSEIF (ALL(t_in==rundef)) THEN
       imiss = imiss + 1
       WRITE(6,*) "T is missing!"
     ELSEIF (ALL(qv_in==rundef)) THEN
       imiss = imiss + 1
       WRITE(6,*) "QV is missing!"
     ELSEIF (ALL(ps_in==rundef)) THEN
       imiss = imiss + 1
       WRITE(6,*) "PS is missing!"
     ENDIF
     IF (ntype_atminput==1) THEN
       IF (ALL(v_in==rundef)) THEN
         imiss = imiss + 1
         WRITE(6,*) "V is missing!"
       ENDIF
     ENDIF
     IF (ntype_radinput==1) THEN
       IF (ALL(t_g_in==rundef)) THEN
         imiss_tg = 1
         WRITE(6,*) "T_G is missing!"
       ENDIF
     ENDIF
     IF (ntype_radinput==2) THEN
       IF (ALL(t_s_in==rundef)) THEN
         IF (ALL(t_so_0_in==rundef)) THEN
           imiss_tg = imiss_tg + 1
           WRITE(6,*) "T_S is missing!"
         ELSE
           t_s_in=t_so_0_in
           WRITE(6,*) "T_SO(0) taken as T_S!"
         ENDIF
       ELSEIF (ALL(t_snow_in==rundef)) THEN
         imiss_tg = imiss_tg + 1
         WRITE(6,*) "T_SNOW is missing!"
       ELSEIF (ALL(w_snow_in==rundef)) THEN
         imiss_tg = imiss_tg + 1
         WRITE(6,*) "W_SNOW is missing!"
       ENDIF
     ENDIF
     If (imiss /= 0) THEN
       WRITE(6,*) " !!! There are ",imiss," required fields missing in read_icongrib (met. forcing) !!!"
       STOP 'read_icongrib'
     ENDIF
     IF ((imiss_tg==1) .AND. (ntype_radinput==1)) THEN
       WRITE(6,*) " !!! STOP because T_G is missing in read_icongrib !!!"
       STOP 'read_icongrib'
     ENDIF
     IF ((imiss_tg/=0) .AND. (ntype_radinput==2)) THEN
       WRITE(6,*) " !!! There are ",imiss_tg," required fields missing to compute T_G in read_icongrib !!!"
       STOP 'read_icongrib'
     ENDIF
!DL####### NEW ##### check if all required fields are available, otherwise STOP #####

     ! YUCHKDAT
     IF (lcheck) THEN
        OPEN(UNIT=444,FILE='YUCHKDAT',FORM='FORMATTED',STATUS='OLD',POSITION='APPEND')
        WRITE(444,*)
        WRITE(444,*) 'Reading ',TRIM(filename)
        CALL stats("u",u_in,llandmask,UNIT=444)
        IF (ntype_atminput==1) THEN
           CALL stats("v",v_in,llandmask,UNIT=444)
        ENDIF
        CALL stats("t",t_in,llandmask,UNIT=444)
        CALL stats("qv",qv_in,llandmask,UNIT=444)
        CALL stats("ps",ps_in,llandmask,UNIT=444)
        IF (ntype_raininput==1) THEN 
           CALL stats("rain_gsp",rain_gsp_in,llandmask,UNIT=444)
           CALL stats("rain_con",rain_con_in,llandmask,UNIT=444)
           CALL stats("snow_gsp",snow_gsp_in,llandmask,UNIT=444)
           CALL stats("snow_con",snow_con_in,llandmask,UNIT=444)
           CALL stats("TOT_PREC",tot_prec_in,llandmask,UNIT=444)
        ENDIF
        CALL stats("sobs2",sobs2_in,llandmask,UNIT=444)
        CALL stats("thbs2",thbs2_in,llandmask,UNIT=444)
!DL
        CALL stats("asobs2",asobs2_in,llandmask,UNIT=444)
        CALL stats("athbs2",athbs2_in,llandmask,UNIT=444)
!DL
        IF (lpar) THEN
           CALL stats("pabs",pabs_in,llandmask,UNIT=444)
           CALL stats("apabs",apabs_in,llandmask,UNIT=444)
        ENDIF
        CALL stats("albrad",albrad_in,llandmask,UNIT=444)

        SELECT CASE (ntype_radinput)
        CASE (1)
           CALL stats("t_g",t_g_in,llandmask,UNIT=444)
        CASE (2)
           CALL stats("t_s",t_s_in,llandmask,UNIT=444)
           CALL stats("t_so_0",t_so_0_in,llandmask,UNIT=444)   !DL
           CALL stats("t_snow",t_snow_in,llandmask,UNIT=444)
           CALL stats("w_snow",w_snow_in,llandmask,UNIT=444)
        END SELECT
        CLOSE(UNIT=444)
     ENDIF

     ! check if this file is ok
     IF ( any_undef ) THEN
        IF ( .NOT. look_forward ) THEN
           WRITE(6,*) "File ",TRIM(filename)," contains undef values"
           STOP
        ENDIF
        IF ( tincr < tincr_max*60 ) THEN
           tincr = tincr + 60
           new_date = date_SUM(idate,tincr)
           WRITE(effective_date,'(i10.10)') new_date/100
           WRITE(6,*) "File ",TRIM(filename)," contains undef values, looking at ", TRIM(effective_date)
        ELSE
           WRITE(6,*) "Maximum gap reached without finding file (tincr_max= ", tincr_max, ")"
           STOP
        ENDIF
     ELSE
        EXIT
     ENDIF
  ENDDO  ! find undef if any

  ! conversions of input variables
  SELECT CASE (ntype_atminput) 
  CASE (1) ! destagger velocities, if requested
!!DL IF (ldestaggeruv) THEN
!!DL    DO i=ie_in,1,-1
!!DL       DO j=je_in,1,-1
!!DL          u_in(i,j)=0.5_wp*(u_in(MAX(1,i-1),j)+u_in(i,j))
!!DL          v_in(i,j)=0.5_wp*(v_in(i,MAX(1,j-1))+v_in(i,j))
!!DL       ENDDO
!!DL    ENDDO
!!DL ENDIF
     u_in(:,:)=SQRT(u_in(:,:)*u_in(:,:)+v_in(:,:)*v_in(:,:)) 
  CASE (2) ! conversion from dew point to specific humidity 
     DO i=1,ie_in
        DO j=1,je_in
           qv_in(i,j)= fqvs ( fpvsw ( qv_in(i,j) ) , ps_in (i,j) )
        ENDDO
     ENDDO
  CASE (3) ! conversion from relative to specific humidity 
     DO i=1,ie_in
        DO j=1,je_in
           qv_in(i,j)= MIN( 1.0_wp, qv_in(i,j) ) * &
                fqvs ( fpvsw ( t_in(i,j) ) , ps_in (i,j) )
        ENDDO
     ENDDO
  END SELECT

  IF (ntype_radinput==2) THEN ! determine t_g
     t_g_in= t_snow_in                                  &
             + (1.0_wp - MIN(1.0_wp,w_snow_in/cf_snow)) &
             * (t_s_in - t_snow_in)
  ENDIF

  IF (ntype_radinput<=2) THEN ! convert net radiation to downwelling radiation
!DL First, check missing fields !!
    IF ( ALL(albrad_in==rundef) ) THEN
      albrad_in = csalb_p
      WRITE(6,*) "ALB_RAD is missing, therefore set to csalb_p!!"
    ENDIF

!!DL    IF ( ALL(sobs2_in==rundef) ) THEN
!!DL       IF ( ALL(asobs2_in==rundef) ) THEN
!!DL         sobs2_in = 0.0
!!DL         WRITE(6,*) "sobs2_in/asobs2_in is missing, therefore set to zero!!"
!!DL       ELSE
!!DL         sobs2_in = asobs2_in
!!DL         WRITE(6,*) "asob_s is taken as sobs_rad!"
!!DL       ENDIF
!!DL    ENDIF
!!DL    sobs2_in(:,:) = sobs2_in(:,:) /(1.0_wp-albrad_in(:,:)/100.0_wp)

!DL Take ASOB_S (averaged field), if missing take SOBS_RAD (instant)
    IF ( ANY(asobs2_in/=rundef) ) THEN
      sobs2_in = asobs2_in
    ELSEIF ( ALL(asobs2_in==rundef) ) THEN
      IF ( ALL(sobs2_in==rundef) ) THEN
        sobs2_in = 0.0_wp
        WRITE(6,*) "SOBS_RAD/ASOB_S is missing, therefore set to zero!!"
      ELSE
!DL     sob2_in already contains sobs_rad
        WRITE(6,*) "SOBS_RAD is taken as ASOB_S!"
      ENDIF
    ENDIF
    sobs2_in(:,:) = sobs2_in(:,:) /(1.0_wp-albrad_in(:,:)/100.0_wp)
    

!!DL    IF ( ALL(thbs2_in==rundef) ) THEN
!!DL      IF ( ALL(athbs2_in==rundef) ) THEN
!!DL        thbs2_in = 0.0
!!DL        WRITE(6,*) "thbs2_in is missing, therefore set to zero!!"
!!DL      ELSE
!!DL        thbs2_in = athbs2_in
!!DL        WRITE(6,*) "athb_s is taken as thbs_rad!"
!!DL        thbs2_in(:,:) = (thbs2_in(:,:)+(1.0_wp-ctalb)*sigma*t_g_in(:,:)**4)/(1.0_wp-ctalb)
!!DL      ENDIF
!!DL    ELSE
!!DL       thbs2_in(:,:) = (thbs2_in(:,:)+(1.0_wp-ctalb)*sigma*t_g_in(:,:)**4)/(1.0_wp-ctalb)
!!DL    ENDIF

!DL Take ATHB_S (averaged field), if missing take THBS_RAD (instant)
    IF ( ANY(athbs2_in/=rundef) ) THEN
      thbs2_in = athbs2_in
    ELSEIF ( ALL(athbs2_in==rundef) ) THEN
      IF ( ALL(thbs2_in==rundef) ) THEN
        thbs2_in = 0.0_wp
        WRITE(6,*) "THBS_RAD/ATHB_S is missing, therefore set to zero!!"
      ELSE
!DL     thbs2_in already contains thbs_rad
        WRITE(6,*) "THBS_RAD taken as ATHB_S!"
      ENDIF
    ENDIF
    thbs2_in(:,:) = (thbs2_in(:,:)+(1.0_wp-ctalb)*sigma*t_g_in(:,:)**4)/(1.0_wp-ctalb)

!!DL sobs2_in(:,:) = sobs2_in(:,:) /(1.0_wp-albrad_in(:,:)/100.0_wp)
!!DL thbs2_in(:,:) = (thbs2_in(:,:)+(1.0_wp-ctalb)*sigma*t_g_in(:,:)**4)/(1.0_wp-ctalb)

  ENDIF  ! ntype_radinput<=2

!!DL also for missing pabs_in
!DL Take APAB_S (averaged field), if missing take PABS_RAD (instant)
!DL if missing or .NOT.lpar take 0.5*sobs2_in
    IF ( ANY(apabs_in/=rundef) ) THEN
      pabs_in = apabs_in
    ELSEIF ( ALL(apabs_in==rundef) ) THEN
      IF ( ALL(pabs_in==rundef) ) THEN
        pabs_in = 0.5_wp*sobs2_in
        WRITE(6,*) "PABS_RAD/APAB_S is missing, therefore set to 0.5*sob2_in!!"
      ELSE
!DL     pabs_in already contains pabs_rad
        WRITE(6,*) "PABS_RAD taken as APAB_S!"
      ENDIF
    ENDIF
 
  IF (ntype_raininput==1) THEN
    ltot_prec=.TRUE.
!DL First, check missing fields !!
    IF ( ALL(rain_gsp_in==rundef) ) THEN 
      rain_gsp_in = 0.0_wp
      WRITE(6,*) "RAIN_GSP is missing, set to zero !"
    ENDIF
    IF ( ALL(rain_con_in==rundef) ) THEN 
      WRITE(6,*) "RAIN_CON is missing, set to zero !"
      rain_con_in = 0.0_wp
    ENDIF
    IF ( ALL(snow_gsp_in==rundef) ) THEN 
      WRITE(6,*) "SNOW_GSP is missing, set to zero !"
      snow_gsp_in = 0.0_wp
    ENDIF
    IF ( ALL(snow_con_in==rundef) ) THEN 
      WRITE(6,*) "SNOW_CON is missing, set to zero !"
      snow_con_in = 0.0_wp
    ENDIF
    IF ( ALL(tot_prec_in==rundef) ) THEN
      WRITE(6,*) "TOT_PREC is missing, set to zero !"
      ltot_prec=.FALSE.
      tot_prec_in = 0.0_wp
    ENDIF
    rain_gsp_in=rain_gsp_in/3600.0_wp
    rain_con_in=rain_con_in/3600.0_wp
    snow_gsp_in=snow_gsp_in/3600.0_wp
    snow_con_in=snow_con_in/3600.0_wp
    tot_prec_in=tot_prec_in/3600.0_wp
  ENDIF


  ! tell user what we are doing with forcing data

  IF (ntype_atminput<=3) THEN
        u_bd(:,:,ke,nx) = u_in 
        t_bd(:,:,ke,nx) = t_in
        qv_bd(:,:,ke,nx)= qv_in
        ps_bd(:,:,nx)   = ps_in
  ENDIF
  IF (ntype_raininput==1) THEN
!.....tot_prec.....
    IF (ltot_prec) THEN
        prr_gsp_bd(:,:,nx) = tot_prec_in
        prs_gsp_bd(:,:,nx) = 0.0_wp
    ELSE
        prr_gsp_bd(:,:,nx) = (rain_gsp_in + rain_con_in)
        prs_gsp_bd(:,:,nx) = (snow_gsp_in + snow_con_in)
    ENDIF

  ENDIF
  IF (ntype_radinput<=2) THEN
        so_down_bd(:,:,nx) = sobs2_in 
        th_down_bd(:,:,nx) = thbs2_in
        pabs_bd   (:,:,nx) = pabs_in
  ENDIF

  ! Clean up
  IF (ntype_atminput<=3) THEN
     DEALLOCATE(u_in, v_in,t_in,qv_in,ps_in )
  ENDIF
  IF (ntype_raininput==1) THEN 
     DEALLOCATE(rain_gsp_in,rain_con_in, snow_gsp_in, snow_con_in, tot_prec_in)
  ENDIF
  IF (ntype_radinput<=2) THEN 
!DL  DEALLOCATE(sobs2_in, thbs2_in, pabs_in, t_g_in,albrad_in)
     DEALLOCATE(sobs2_in, asobs2_in, thbs2_in, athbs2_in, pabs_in, apabs_in, t_g_in,albrad_in)
  ENDIF
  IF (ntype_radinput==2) THEN
!DL  DEALLOCATE(w_snow_in,t_snow_in,t_s_in)
     DEALLOCATE(w_snow_in,t_snow_in,t_s_in,t_so_0_in)
  ENDIF


  ! Control output
  IF (lcheck) THEN
     IF (ntype_atminput<=3) THEN
        CALL stats("u_bd",u_bd(:,:,ke,nx),llandmask)
        CALL stats("t_bd",t_bd(:,:,ke,nx),llandmask)
        CALL stats("qv_bd",qv_bd(:,:,ke,nx),llandmask)
!       CALL stats("qv_bd",qv_bd(:,:,ke),llandmask)   !XYZ> in TSA, qv is 4 dimensional
        CALL stats("ps_bd",ps_bd(:,:,nx),llandmask)
     ENDIF
     IF (ntype_raininput==1) THEN 
        CALL stats("prr_gsp_bd",prr_gsp_bd(:,:,nx),llandmask)
        CALL stats("prs_gsp_bd",prs_gsp_bd(:,:,nx),llandmask)
     ENDIF
     IF (ntype_radinput<=2) THEN 
        CALL stats("so_down_bd",so_down_bd(:,:,nx),llandmask)
        CALL stats("th_down_bd",th_down_bd(:,:,nx),llandmask)
        CALL stats("pabs_bd",pabs_bd(:,:,nx),llandmask)
     ENDIF
  ENDIF

END SUBROUTINE read_icongrib
!----------------------------------------------------------------------



SUBROUTINE read_radolan(date,nx)
!----------------------------------------------------------------------
! Description:
!
! Reads radiation files - RADOLAN 
!----------------------------------------------------------------------


  IMPLICIT NONE 

  ! Subroutine Arguments
  CHARACTER (LEN=*),INTENT(IN)          :: date
  INTEGER                 , INTENT(IN)  :: nx

  ! Parameters
  INTEGER                 , PARAMETER   :: &
       ie_in=900                         , & !
       je_in=900                             !

  ! Local variables
  INTEGER                               :: &
       i,j,k,l,n                         , & ! do loops index variables
       ilaenge                           , & !
       ierr                                  !
  CHARACTER (LEN=150)                   :: head
  REAL (KIND=wp)                        :: obsfield(ie_in,je_in)  
  INTEGER                 , ALLOCATABLE :: w_buffer(:)
  REAL (KIND=wp)                        :: &
       rlat_rad_in                       , & ! 
       rlon_rad_in                       , & ! 
       rlon_geo                          , & ! 
       rlat_geo                          , & ! 
       rlat0_in                          , & ! 
       rlon0_in                          , & ! 
       dlat_in                           , & ! 
       dlon_in                           , & ! 
       lon,pi                                ! 
  REAL (KIND=wp)                        :: fac
  REAL (KIND=wp), ALLOCATABLE, SAVE     :: &
       rlon_in(:,:)                      , & ! 
       rlat_in(:,:)                          ! 
  LOGICAL, ALLOCATABLE                  :: mask1(:,:)
  LOGICAL, SAVE                         :: lfirst=.TRUE.


  ! Initialization
  ALLOCATE ( mask1(ie_in,je_in) )
  IF (lfirst)  THEN
     ALLOCATE ( rlon_in(ie_in,je_in),  &
                rlat_in(ie_in,je_in) )
  ENDIF

  ! Reading Radolan file
  OPEN(10,FILE =TRIM(radofiledir)//'/radorw'//TRIM(date),        &
       FORM='unformatted', ACCESS='direct', RECL=150, STATUS='old')
  READ(10,REC=1,IOSTAT=ierr) head
  CLOSE(10)

  IF (ierr/=0) THEN
     prr_gsp_bd(:,:,nx)=0.0_wp
     WRITE(6,*) "Radolanfile empty: "//TRIM(date)
  ELSE
     READ(head(20:26),'(i7)') ilaenge
     CLOSE(10)

     OPEN(10,FILE =TRIM(radofiledir)//'/radorw'//TRIM(date),&
          FORM='unformatted', ACCESS='direct', RECL=ilaenge, STATUS='old')
     ALLOCATE(w_buffer(ilaenge/2))
     READ(10,REC=1) w_buffer
     CLOSE(10)
     
     j=(ilaenge-2*ie_in*je_in+2)/2   ! XYZ> Replaced 900 with ie_in, je_in
     IF (j<1) THEN
        prr_gsp_bd(:,:,nx)=0.0_wp
        WRITE(6,*) "Radolanfile errorneous: "//TRIM(date)
     ELSE
        fac=0.1_wp/3600.0_wp
        DO l=1,ie_in                 ! XYZ> Replaced 900 with ie_in
           DO k=1,je_in              ! XYZ> Replaced 900 with je_in
              IF ( w_buffer(j) == 2500        &
                   .OR. w_buffer(j) == 8192   &
                   .OR. w_buffer(j) == 10692  &
                   .OR.w_buffer(j) < 0)       THEN
                 obsfield(k,l) = -9999
              ELSE
                 obsfield(k,l) = w_buffer(j)*fac    
              ENDIF
              j=j+1
           ENDDO
        ENDDO
        DEALLOCATE(w_buffer)

        IF (lfirst) THEN
           ! Calulation of coordinates
           rlon0_in=-521.952_wp
           rlat0_in=-4657.868_wp
           dlon_in=1.0_wp
           dlat_in=1.0_wp
!           pi=4.0*atan(1.0)
           DO i=1,ie_in
              rlon_rad_in=(i-1)*dlon_in+rlon0_in  
              DO j=1,je_in     
                 rlat_rad_in=(j-1)*dlat_in+rlat0_in    
                 rlon_geo=10.0_wp+((ATAN((rlon_rad_in)/(-rlat_rad_in)))/pi*180.0_wp)
                 rlat_geo=ASIN((((6370.04_wp**2)*((1.0_wp+sin(60.0_wp*pi/180.0_wp))**2)) &
                      - ((rlon_rad_in**2)+(rlat_rad_in**2)))                             &
                      / (((6370.04_wp**2)*((1+sin(60.0_wp*pi/180.0_wp))**2))             &   
                      + ((rlon_rad_in**2)+(rlat_rad_in**2))))*180.0_wp/pi
                 rlat_in(i,j)=phi2phirot(rlat_geo,rlon_geo,pollat,pollon)
                 rlon_in(i,j)=rla2rlarot(rlat_geo,rlon_geo,pollat,pollon,0.0_wp)  ! XYZ> added polgam=0 for compatability with utilities
              ENDDO
           ENDDO
           lfirst=.FALSE.
        ENDIF

        ! Spatial Interpolations
        mask1=.TRUE.
        WHERE (obsfield(:,:)<0) 
           mask1(:,:)=.FALSE.
        END WHERE
        CALL general_interpol2d(rlat_in,rlon_in,             &
             REAL(obsfield,wp),                              &                  
             startlat,startlon,dlat,dlon,prr_gsp_bd(:,:,nx), &
             2,2,0.51_wp,1.0_wp,mask1,llandmask)
        WHERE (prr_gsp_bd(:,:,nx)<0) 
           prr_gsp_bd(:,:,nx)=0.0_wp
        END WHERE

        ! Clean up 
        DEALLOCATE(mask1)
     ENDIF  ! j<1
  ENDIF  ! iserr/=0

  ! Control output
  IF (lcheck) THEN
     WRITE(6,*) "read_radolan: prr_gsp",&
          MINVAL(prr_gsp(:,:),llandmask), &
          MAXVAL(prr_gsp(:,:),llandmask)
  ENDIF

END SUBROUTINE read_radolan

!==============================================================================

END MODULE tsa_input
