! =====================================================================
! Subroutines for Output of Terra Stand Alone (TSA)
! =====================================================================

MODULE  tsa_output

!------------------------------------------------------------------------------
!
! Description:
!  This file contains subroutines for I/O operation
!   of the stand Alone version of the soil module Terra.
!  This file is included in terra when compiling
!
!  Currently included: Subroutines for output:
!
!    - grbout
!      Writes output as GRIB file using libDWD
!
!    - grbout_eccodes
!      Writes output as GRIB file using eccodes
!
!    - binout
!      Writes output as BIN file
!
!    - ascout
!      Writes output as ASCII file
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
!  Added eccodes for writing grib files
!  General refactoring: tsa_output now contains only routines for output
!    from former terra_io.f90
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
    so_down_bd, th_down_bd, pabs_bd, ps_bd, uuid, astore, store,            &
    runoff, rlat_geo, rlon_geo, ydate_ini, ydate_end,                       &
    nout_interval, outdir, outprefix, ntype_output, lconstout,              &
    radofiledir, metfiledir, metfileprefix, &
    ntype_atminput, ntype_radinput, ntype_raininput, rundef,                &
    w_g0, w_i0, t_soil0, w_snow0, t_cl0, w_cl0,                             &
    soilinitdir, soilinitprefix, soiltyp_const,                             &
    ke_model, ngp, ngpi, ngrid, ymodel, yform_write, yform_read,            &
    ntstep_max, tincr_max, constfilename, nsub_x, nsub_y, which_subreg,     &
    lrel_in        , lvol_in        , lvegadapt      , lhomosoil      ,     &
    lhomoinit      , lconstout      , lrootadapt     , lcheck         ,     &
    lmulti_in      ,                  lz0local       , lconstvegalb   ,     &
    lext_monthly   , linfil_revised , lhourly_data   , ldestaggeruv   ,     &
    lpar           , lgettcl        , lcalc          ,                      &
    lai_const      , lai_min_const  , lai_max_const  , plcov_const    ,     &
    plcov_min_const, plcov_max_const, rootdp_const   , plcov_mn       ,     &
    plcov_mx       , rootdp0        , lai_mn         , lai_mx         ,     &
    z0             , rstom_mx       , rstom_mn       , vegalb         ,     &
    z0_const       , rstom_mn_const , rstom_mx_const , vegalb_const   ,     &
    dz             , dz_u           , nrecmax        , hsurface       ,     &
    rain_fac       , lat

USE data_io,               ONLY:                                            &
    intgribf, intgribc, irealgrib, iblock, lfd, lds, iwlength,              &
    iwlength

USE data_fields,           ONLY:                                            &
    u, t, ps, p0, prr_gsp, prs_gsp, u_10m, v_10m, t_2m, t_s, t_so, w_so,    &
    w_snow, rho_snow, freshsnow, t_snow, t_g, w_i, qv_s, tch, tcm,          &
    pabs, sobs, thbs, llandmask, lai, plcov, rootdp, soiltyp, hsurf,        &
    fr_land, rlon, rlat, t_cl, w_cl, w_so_ice, rsmin2d,                     &
    u_bd, t_bd , h_snow

USE data_modelconfig,      ONLY:                                            &
    pollon, pollat, polgam, dlon, dlat, startlon, startlat, ke_soil, czmls, &
    ie, je, ke, dt
    
USE data_runcontrol,       ONLY:                                            &
    nnew, nnow, nblock, nproma, nlastproma, lstomata, itype_canopy,         &
    itype_evsl, itype_heatcond, itype_hydbound, itype_mire, itype_root,     &
    lmelt, lmelt_var, lmulti_snow, ntstep, itype_calendar

USE tsa_lmparam,            ONLY:                                           &
    stats

USE utilities,              ONLY:                                           &
    get_utc_date

!==============================================================================

IMPLICIT NONE

!==============================================================================

CONTAINS

!======================================================================
! OUTPUT
!======================================================================


SUBROUTINE grbout(nt)
!----------------------------------------------------------------------
! Description:
!
! Writes output as GRIB file
!----------------------------------------------------------------------


  IMPLICIT NONE

  ! Subroutine Arguments
  INTEGER                 , INTENT(IN) :: nt

  ! Local variables
  INTEGER                              :: &
       idayofyear                       , & ! 
       minute                           , & ! 
       nvar                             , & ! 
       ii                                   ! 
  REAL (KIND=wp)                       :: acthour
  CHARACTER (LEN=14)                   :: yactdate1
  CHARACTER (LEN=11)                   :: str          ! GDM
  CHARACTER (LEN=28)                   :: yactdate2
  CHARACTER (LEN=14), SAVE             :: yactdate1_old
  LOGICAL                              :: lconst
  CHARACTER (LEN=100)                  :: filename

  ! Grib variables
  INTEGER (KIND=intgribf),  PARAMETER  :: &
       npds   =    321_intgribf         , & ! Dimension for product definition section (pds)
       ngds   =    626_intgribf         , & ! Dimension for grid description section (gds)
       nbms   =      3_intgribf         , & ! Dimension for bit map section (bms)
       nbds   =     11_intgribf         , & ! Dimension for binary data section
       ndsup  =     73_intgribf         , & ! Dimension for dsup
       ndims  =     20_intgribf             ! Dimension for idims (contains all dimensions)

  REAL    (KIND=irealgrib), PARAMETER  :: &
       undefgrib =  -1.0E7_irealgrib        ! value for "undefined" in the grib routines  
  INTEGER (KIND=intgribf)              :: &
       iednr     =   1                      ! grib edition number 
  INTEGER (KIND=intgribf)              :: &
       iz_ps=1                          , & ! 
       ierrf                                ! 

  ! The following dimensions are set during program execution
  INTEGER (KIND=intgribf)              :: &
       lfd                              , & ! Dimension for iblock
       lbm                              , & ! Dimension for bitmap: this value has to be at
                                            !   least the same as lbmax in subroutine grbin1
       lds                                  ! Dimension for unpacked data

  INTEGER (KIND=intgribf)              :: &
       idims (ndims)                    , & ! array for all dimensions (output)
       ipds  (npds)                     , & ! product description section
       igds  (ngds)                     , & ! grid description section
       ibms  (nbms)                     , & ! bit map section
       ibds  (nbds)                         ! binary data section

  ! Arrays that are allocated during program execution
  INTEGER (KIND=intgribf), ALLOCATABLE :: &
       iblock    (:)                    , & ! array for gribed data
       ibmap     (:)                        ! array for bit map

  REAL   (KIND=irealgrib), ALLOCATABLE :: &
       dsup (:)                         , & ! array for special data
       ds   (:)                         , & ! array for unpacked data
       ds2(:,:)                             !

  INTEGER (KIND=intgribf)              :: &
       ilen                             , & ! 
       ierr                             , & ! 
       istat                                ! 
  INTEGER  (KIND=intgribc)             :: & ! corresponding variables for C-routines
       nudatc                           , & ! 
       maxlenc                          , & ! 
       ilenc                            , & ! 
       ierrc                                !     
  REAL     (KIND=wp)                   :: pollon_grib
  INTEGER (KIND=intgribf), EXTERNAL    :: IREFTS
  REAL   (KIND=irealgrib)              :: rmin
  ! Looping indices
  INTEGER                              :: k,i,j

! Body of subroutine:
  ! time calculations and handling
  IF (nt==1) THEN
     yactdate1_old=""
  ENDIF
  CALL get_utc_date (nt, ydate_ini, dt, itype_calendar,        & 
       yactdate1, yactdate2, idayofyear,acthour)
  lconst=((yactdate1(1:8) /= yactdate1_old(1:8)).AND.(lconstout))
  yactdate1_old=yactdate1
  minute=MODULO(NINT(acthour*60.0_wp),60)
  IF (minute<10) THEN
     WRITE(filename,'(a,a,a,a,a,I1)') &
          TRIM(outdir),'/',TRIM(outprefix),yactdate1(1:10),'0',minute
  ELSE
     WRITE(filename,'(a,a,a,a,I2)')   &
          TRIM(outdir),'/',TRIM(outprefix),yactdate1(1:10),minute
  ENDIF

  ! tell user what we are doing
  WRITE(6,*) "Writing Gribfile ",TRIM(filename)

  ! Grib Initializations
  lds = ie * je
  lbm = 1875
  lfd = lds * 2 / iwlength + 2000   ! the "2" means 2 bytes per word
                                    !   the "2000" just is safety

  ALLOCATE (iblock(lfd), ibmap(lbm), STAT=istat)
  ALLOCATE (ds(lds), dsup(ndsup),ds2(ie,je), STAT=istat)

  ! dimensions for grib routines
  idims( 1) = npds
  idims( 2) = ngds
  idims( 3) = nbms
  idims( 4) = nbds
  idims( 5) = lbm
  idims( 6) = ndsup
  idims( 7) = lds
  idims( 8) = lfd

  !  REAL dimensions
  idims(11) = 47
  idims(12) = 25+1+1+4
  idims(13) = 3
  idims(14) = 5
  idims(15) = ie*je
  idims(16) = 0
  idims(17) = ie*je

  ! gridpoints, simple packing, floating point data
  ibds(2)   = 0
  ! nrbit, number of bits
  ibds(5)   = 16  ! GDM
  ! no bitmap
  ibms(3)   = -2

  nudatc=251
  CALL copen (nudatc,TRIM(filename),'w  ',ierrc)

  ! Grid Description Sector
  igds(:) = -999999
  ! length of igds in bytes
  igds(1) = 42 + (4 + 1+1) * 4
  ! number of the vertical coordinate parameters
  !     igds(2) = ke_tot+1 + 4
  igds(2) = 0
  ! location of the list of vertical coordinate parameters in bytes
  igds(3) = 43
  ! data representation type
  igds(4) = 10

  ! calculation of the left bottom corner and
  ! left bottom and right upper corner in millidegrees:
  ! This depends on the variable (U- and V-related variables are shifted in
  ! the Arakawa C-grid) and is determined during the output step.
  ! number of gridpoints
  igds( 5) = ie
  igds( 6) = je
  igds( 9) = 0

  ! increments
  igds(12) = 0
  igds(13) = 0

  igds(14) = 64
  igds(15:19) = 0

  ! coordinates of the pole: in the grib code the southern pole has to be
  ! specified: for the latitude this is the negative value, for the
  ! longitude it is the value + 180.0 and then limiting the value again
  ! to -180.0 ... 180.0
  pollon_grib = pollon + 180.0_wp
  IF (pollon_grib > 180.0_wp) THEN
     pollon_grib = pollon_grib - 360.0_wp
  ENDIF

  igds(20) = NINT(-pollat      * 1000.0_wp)
  igds(21) = NINT( pollon_grib * 1000.0_wp)
  igds(22) = IREFTS(0.0_irealgrib)

  ! vertical coordinate parameters
  igds(26) = irefts( REAL (0,   irealgrib) )
  igds(27) = irefts( REAL (0,   irealgrib) )
  igds(28) = irefts( REAL (0,  irealgrib) )
  igds(29) = irefts( REAL (0, irealgrib) )

  DO k = 1, 1+1
     igds(29 + k) = irefts( REAL (0, irealgrib) )
  ENDDO

  igds(7)=NINT(startlat*1000.0_wp)
  igds(8)=NINT(startlon*1000.0_wp)
  igds(10)=NINT((startlat+(je-1)*dlat)*1000.0_wp)
  igds(11)=NINT((startlon+(ie-1)*dlon)*1000.0_wp)

  ! Product Definitions Sector
  ipds(1)=54
  ipds(3)=78
  ipds(4)=131
  ipds(5)=255
  ipds(6)=128
  ipds(37)=254

  READ(yactdate1(1:10),"(5(I2))") ipds(22), &
       ipds(11),ipds(12),ipds(13),    &
       ipds(14)
  ipds(15) = minute
  ipds(22) = ipds(22)+1
  ipds(19) = 0
  ipds(17) = 0
  ipds(18) = 0
  ipds(16) = 1
  ipds(20) = 0
  ipds(21) = 0
  ipds(24) = 0
  ipds(38:46) = 0
  ipds(47) = 1

  ! Output of fields

  ! convert sums into averages
  astore(:,:,1)=astore(:,:,1)/FLOAT(nout_interval)  ! lhfl_s
  astore(:,:,2)=astore(:,:,2)/FLOAT(nout_interval)  ! shfl_s
  astore(:,:,5)=astore(:,:,5)/FLOAT(nout_interval)  ! sobs
  astore(:,:,6)=astore(:,:,6)/FLOAT(nout_interval)  ! thbs
  astore(:,:,7)=astore(:,:,7)/FLOAT(nout_interval)  ! runoff_s
  astore(:,:,8)=astore(:,:,8)/FLOAT(nout_interval)  ! runoff_g
   
  ! convert rates into totals
  astore(:,:,3)=astore(:,:,3)*dt                    ! prr_gsp
  astore(:,:,4)=astore(:,:,4)*dt                    ! prs_gsp
  
  ! GDM> produce YUCHKDAT with stats
  IF (lcheck) THEN
    OPEN(UNIT=444,FILE='YUCHKDAT',FORM='FORMATTED',STATUS='OLD',POSITION='APPEND')
    WRITE(444,*)
    WRITE(444,*) 'Writing ',TRIM(filename)
  ENDIF

  SELECT CASE (ntype_output)
  CASE (1)
     nvar=ke_soil+1
  CASE (2)
     nvar=2*(ke_soil+1)+4
  CASE (3)
     nvar=3*(ke_soil+1)+28
  END SELECT

  DO i=1,nvar

     IF (i<=ke_soil+1) THEN
        ipds(7)=198
        ipds(2)=201
        ipds(8)=111
        ipds(9)=0
        ipds(10)=INT(czmls(i)*100.0+0.99)
        ds2=REAL(w_so(:,:,i,nnew)*1000.0,irealgrib)         
        WRITE(str,'(A5,I4)') "w_so_",INT(czmls(i)*100.0+0.99)
        IF (lcheck) THEN
           CALL stats(str,w_so(:,:,i,nnew)*1000.0,llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==(ke_soil+1)+1) THEN
        ipds(7)=197
        ipds(2)=201
        ipds(8)=111
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(t_so(:,:,0,nnew),irealgrib)   
        IF (lcheck) THEN
           CALL stats("t_so_0",t_so(:,:,0,nnew),llandmask,UNIT=444)
        ENDIF
      ELSE IF (i<=2*(ke_soil+1)+1) THEN
        ii=i-(ke_soil+2)
        ipds(7)=197
        ipds(2)=201
        ipds(8)=111
        ipds(9)=0
        ipds(10)=INT(czmls(ii)*100.0+0.99)
        ds2=REAL(t_so(:,:,ii,nnew),irealgrib)   
        WRITE(str,'(A5,I4)') "t_so_",INT(czmls(ii)*100.0+0.99)
        IF (lcheck) THEN
           CALL stats(str,t_so(:,:,ii,nnew),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+2) THEN
        ipds(7)=203
        ipds(2)=201
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(t_snow(:,:,nnew),irealgrib)   
        IF (lcheck) THEN
           CALL stats("t_snow",t_snow(:,:,nnew),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+3) THEN
        ipds(7)=65
        ipds(2)=2
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(w_snow(:,:,nnew)*1000.0,irealgrib)   
        IF (lcheck) THEN
           CALL stats("w_snow",w_snow(:,:,nnew)*1000.0,llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+4) THEN
        ipds(7)=200
        ipds(2)=201
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(w_i(:,:,nnew)*1000.0,irealgrib)   
        IF (lcheck) THEN
           CALL stats("w_i",w_i(:,:,nnew)*1000.0,llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+5) THEN
! GDM> ke instead of 30m
        ipds(7)=11
        ipds(2)=2
        ipds(8)=110
        ipds(9)=ke_model
        ipds(10)=ke_model+1
        ds2=REAL(t(:,:,ke,nnew),irealgrib)   
        WRITE(str,'(A2,I2)') "t_k",ke_model
        IF (lcheck) THEN
           CALL stats(str,t(:,:,ke,nnew),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+6) THEN
        ipds(7)=51
        ipds(2)=2
        ipds(8)=110
        ipds(9)=ke_model
        ipds(10)=ke_model+1
        ds2=REAL(qv(:,:,ke,nnew),irealgrib)
!       ds2=REAL(qv(:,:,ke),irealgrib)   
     WRITE(str,'(A3,I2)') "qv_k",ke_model
        IF (lcheck) THEN
           CALL stats(str,qv(:,:,ke,nnew),llandmask,UNIT=444)
        ENDIF
!XYZ> qv in TSA is 4 dimensional
!       IF (lcheck) THEN
!          CALL stats(str,qv(:,:,ke),llandmask,UNIT=444)
!       ENDIF
!XYZ<
     ELSE IF (i==2*(ke_soil+1)+7) THEN
        ipds(7)=32
        ipds(2)=2
        ipds(8)=110
        ipds(9)=ke_model
        ipds(10)=ke_model+1
        ds2=REAL(u(:,:,ke,nnew),irealgrib)   
        WRITE(str,'(A2,I2)') "u_k",ke_model
        IF (lcheck) THEN
           CALL stats(str,u(:,:,ke,nnew),llandmask,UNIT=444)
        ENDIF
!GDM<
     ELSE IF (i==2*(ke_soil+1)+8) THEN
        ipds(7)=1
        ipds(2)=2
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(ps(:,:,nnew),irealgrib)   
        IF (lcheck) THEN
           CALL stats("ps",ps(:,:,nnew),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+9) THEN
        ipds(7)=102
        ipds(2)=201
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(astore(:,:,3),irealgrib)         ! prr/rain_gsp
        IF (lcheck) THEN
           CALL stats("prr_gsp",astore(:,:,3),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+10) THEN
        ipds(7)=79
        ipds(2)=2
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(astore(:,:,4),irealgrib)         ! prs/snow_gsp
        IF (lcheck) THEN
           CALL stats("prs_gsp",astore(:,:,4),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+11) THEN
        ipds(7)=111
        ipds(2)=2
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(astore(:,:,5),irealgrib)         ! sobs_s
        IF (lcheck) THEN
           CALL stats("sobs_s",astore(:,:,5),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+12) THEN
        ipds(7)=112
        ipds(2)=2
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(astore(:,:,6),irealgrib)         ! thbs_s
        IF (lcheck) THEN
            CALL stats("thbs_s",astore(:,:,6),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+13) THEN
        ipds(7)=121
        ipds(2)=2
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(astore(:,:,1),irealgrib)         ! lhfl_s
        IF (lcheck) THEN
           CALL stats("lhfl_s",astore(:,:,1),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+14) THEN
        ipds(7)=122
        ipds(2)=2
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(astore(:,:,2),irealgrib)         ! shfl_s
        IF (lcheck) THEN
            CALL stats("shfl_s",astore(:,:,2),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+15) THEN
        ipds(7)=90
        ipds(2)=2
        ipds(8)=112
        ipds(9)=0
        ipds(10)=10
        ds2=REAL(astore(:,:,7),irealgrib)         ! runoff_s
        IF (lcheck) THEN
           CALL stats("runoff_s",astore(:,:,7),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+16) THEN
        ipds(7)=90
        ipds(2)=2
        ipds(8)=112
        ipds(9)=10
        ipds(10)=100
        ds2=REAL(astore(:,:,8),irealgrib)         ! runoff_g
        IF (lcheck) THEN
           CALL stats("runoff_g",astore(:,:,8),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+17) THEN
        ipds(7)=171
        ipds(2)=201
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(tch(:,:),irealgrib)             ! tch
        IF (lcheck) THEN
           CALL stats("tch",tch(:,:),llandmask,UNIT=444)
        ENDIF
!GDM>  store  
     ELSE IF (i==2*(ke_soil+1)+18) THEN
        ipds(7)=21
        ipds(2)=201
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(store(:,:,1),irealgrib)         ! rstom
        IF (lcheck) THEN
           CALL stats("rstom",store(:,:,1),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+19) THEN
        ipds(7)=251
        ipds(2)=2
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(store(:,:,2),irealgrib)         ! zf_wat
        IF (lcheck) THEN
           CALL stats("zf_wat",store(:,:,2),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+20) THEN
        ipds(7)=252
        ipds(2)=2
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(store(:,:,3),irealgrib)         ! zf_tem
        IF (lcheck) THEN
           CALL stats("zf_tem",store(:,:,3),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+21) THEN
        ipds(7)=253
        ipds(2)=2
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(store(:,:,4),irealgrib)         ! zf_rad
        IF (lcheck) THEN
           CALL stats("zf_rad",store(:,:,4),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+22) THEN
        ipds(7)=254
        ipds(2)=2
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(store(:,:,5),irealgrib)         ! zep_s 
        IF (lcheck) THEN
           CALL stats("zep_s",store(:,:,5),llandmask,UNIT=444)
        ENDIF
! GDM< store
     ELSE IF (i<=3*(ke_soil+1)+22) THEN
        ii=i-(2*(ke_soil+1)+22)
        ipds(7)=255
        ipds(2)=201
        ipds(8)=111
        ipds(9)=0
        ipds(10)=INT(czmls(ii)*100.0+0.99)
        ds2=REAL(runoff(:,:,ii),irealgrib)         
        WRITE(str,'(A7,I4)') "runoff_",INT(czmls(ii)*100.0+0.99)
        IF (lcheck) THEN
           CALL stats(str,runoff(:,:,ii),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==3*(ke_soil+1)+23) THEN
        ipds(7)=11
        ipds(2)=2
        ipds(8)=105
        ipds(9)=0
        ipds(10)=2
        ds2=REAL(t_2m,irealgrib)         
        IF (lcheck) THEN
           CALL stats("t_2m",t_2m,llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==3*(ke_soil+1)+24) THEN
        ipds(7)=133
        ipds(2)=201
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(rho_snow(:,:,nnew),irealgrib)   
        IF (lcheck) THEN
           CALL stats("rho_snow",rho_snow(:,:,nnew),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==3*(ke_soil+1)+25) THEN
        ipds(7)=129
        ipds(2)=201
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(freshsnow(:,:),irealgrib)   
        IF (lcheck) THEN
           CALL stats("freshsnow",freshsnow(:,:),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==3*(ke_soil+1)+26) THEN
        ipds(7)=170
        ipds(2)=201
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(tcm(:,:),irealgrib)   
        IF (lcheck) THEN
          CALL stats("tcm",tcm(:,:),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==3*(ke_soil+1)+27) THEN
        ipds(7)=51
        ipds(2)=2
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(qv_s(:,:,nnew),irealgrib)   
        IF (lcheck) THEN
          CALL stats("qv_s",qv_s(:,:,nnew),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==3*(ke_soil+1)+28) THEN
        ipds(7)=11
        ipds(2)=2
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(t_g(:,:,nnew),irealgrib)   
        IF (lcheck) THEN
          CALL stats("t_g",t_g(:,:,nnew),llandmask,UNIT=444)
        ENDIF
     ENDIF

     ! undef handling for non-landpoints
!DL uncomment for testing
     !rmin=MINVAL(ds2,llandmask)
     rmin=MINVAL(ds2,llandmask)
     WHERE (.NOT.llandmask)
!DL for testing ....
!DL     ds2(:,:)=undefgrib !rmin
        ds2(:,:)=rmin
     END WHERE

     ! write data
     ds=RESHAPE(ds2,(/ ie*je /))
     CALL grbex1(iednr, iz_ps, undefgrib, ndims, idims, &                  ! Coding to GRIB
          ipds, igds, ibms, ibds, ibmap, dsup , ds, iblock, ierrf)
     ilenc=idims(18)*iwlength
     CALL cuegex (nudatc, iblock, ilenc, ierrc)                            ! Write to file

  ENDDO

  ! write constant fields
  nvar=0
  IF (lconst) THEN
     nvar=nvar+9
  ENDIF
  DO i=1,nvar

     IF (i==1) THEN
        ipds(7)=61
        ipds(2)=202
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(lai(:,:),irealgrib)   
        IF (lcheck) THEN
          CALL stats("lai_out",lai,llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2) THEN
        ipds(7)=87
        ipds(2)=2
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(plcov(:,:),irealgrib)   
        IF (lcheck) THEN
          CALL stats("plcov_out",plcov,llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==3) THEN
        ipds(7)=62
        ipds(2)=202
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(rootdp(:,:),irealgrib)   
        IF (lcheck) THEN
          CALL stats("rootdp_out",rootdp,llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==4) THEN
        ipds(7)=57
        ipds(2)=202
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(soiltyp(:,:),irealgrib)   
        IF (lcheck) THEN
          CALL stats("soiltyp_out",soiltyp,llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==5) THEN
        ipds(7)=8
        ipds(2)=2
        ipds(8)=109
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(hsurf(:,:),irealgrib)   
        IF (lcheck) THEN
          CALL stats("hsurf_out",hsurf,llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==6) THEN
        ipds(7)=83
        ipds(2)=2
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(z0(:,:),irealgrib)   
        IF (lcheck) THEN
          CALL stats("z0_out",z0,llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==7) THEN
        ipds(7)=81
        ipds(2)=2
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(fr_land(:,:),irealgrib)   
        IF (lcheck) THEN
          CALL stats("fr_land_out",fr_land,llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==8) THEN
        ipds(7)=115
        ipds(2)=202
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(rlon(:,:),irealgrib)   
        IF (lcheck) THEN
          CALL stats("rlon_out",rlon,llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==9) THEN
        ipds(7)=114
        ipds(2)=202
        ipds(8)=1
        ipds(9)=0
        ipds(10)=0
        ds2=REAL(rlat(:,:),irealgrib)   
        IF (lcheck) THEN
          CALL stats("rlat_out",rlat,llandmask,UNIT=444)
        ENDIF
     ENDIF

     ! undef handling for non-land points
     IF (i<=3) THEN
        !rmin=MINVAL(ds2,llandmask)
        WHERE (.NOT.llandmask)
           ds2(:,:)=undefgrib !rmin
        END WHERE
     ENDIF

     ! write data
     ds=RESHAPE(ds2,(/ ie*je /))
     CALL grbex1(iednr, iz_ps, undefgrib, ndims, idims, &                  ! Coding to GRIB
          ipds, igds, ibms, ibds, ibmap, dsup , ds, iblock, ierrf)
     ilenc=idims(18)*iwlength
     CALL cuegex (nudatc, iblock, ilenc, ierrc)                            ! Write to file

  ENDDO

  ! reset summation storage
  astore=0.0

  ! YUCHKDAT
  IF (lcheck) THEN
    CLOSE(UNIT=444)
  ENDIF

  CALL cclose(nudatc,'exi',ierrc)
  DEALLOCATE (iblock, ibmap)
  DEALLOCATE (ds, dsup)

END SUBROUTINE grbout
!----------------------------------------------------------------------


SUBROUTINE grbout_eccodes(nt)
!----------------------------------------------------------------------
! Description:
!
! Writes output as GRIB file using eccodes (GRIB1 AND GRIB2)
!----------------------------------------------------------------------

  USE grib_api
  IMPLICIT NONE

  ! Subroutine Arguments
  INTEGER                 , INTENT(IN) :: nt

  ! Local variables
  INTEGER                              :: &
       idayofyear                       , & ! 
       nzminute, nzdate, nzhour         , & ! 
       nvar                             , & ! 
       idim, jdim, ii                                   ! 

  REAL (KIND=wp)                       :: acthour
  CHARACTER (LEN=14)                   :: yactdate1
  CHARACTER (LEN=11)                   :: str          ! GDM
  CHARACTER (LEN=28)                   :: yactdate2
  CHARACTER (LEN=14), SAVE             :: yactdate1_old
  LOGICAL                              :: lconst, lbitmap
  CHARACTER (LEN=100)                  :: filename

  INTEGER                              :: &
       igribo, igrib1_id, igrib2_id         ! grib handles

  REAL    (KIND=irealgrib), PARAMETER  :: &
       undefgrib =  -1.0E7_irealgrib        ! value for "undefined" in the grib routines  
  INTEGER (KIND=intgribf)              :: &
       iednr                                ! grib edition number 

  ! The following dimensions are set during program execution
  INTEGER (KIND=intgribf)              :: &
       lds                                  ! Dimension for unpacked data

  REAL   (KIND=irealgrib), ALLOCATABLE :: &
       ds   (:)                         , & ! array for unpacked data
       ds2(:,:)                             !

  INTEGER (KIND=intgribf)              :: &
       ilen                             , & ! 
       ierr                             , & ! 
       ifac                             , & !
       istat                                ! 
     INTEGER  (KIND=intgribc)             :: & ! corresponding variables for C-routines
          nudatc                           , & ! 
          ierrc                                !     

  REAL     (KIND=wp)                   :: pollon_grib, zfac

!DL  INTEGER (KIND=intgribf), EXTERNAL    :: IREFTS
  REAL   (KIND=irealgrib)              :: rmin
  ! Looping indices
  INTEGER                              :: k,i,j

! Body of subroutine:
  ! time calculations and handling
  IF (nt==1) THEN
     yactdate1_old=""
  ENDIF
  CALL get_utc_date (nt, ydate_ini, dt, itype_calendar,       & 
       yactdate1, yactdate2, idayofyear,acthour)
  lconst=((yactdate1(1:8) /= yactdate1_old(1:8)).AND.(lconstout))
  yactdate1_old=yactdate1
  nzminute=MODULO(NINT(acthour*60.0_wp),60)
  IF (nzminute<10) THEN
     WRITE(filename,'(a,a,a,a,a,I1)') &
          TRIM(outdir),'/',TRIM(outprefix),yactdate1(1:10),'0',nzminute
  ELSE
     WRITE(filename,'(a,a,a,a,I2)')   &
          TRIM(outdir),'/',TRIM(outprefix),yactdate1(1:10),nzminute
  ENDIF

  ! tell user what we are doing
  WRITE(6,*) "Writing Gribfile ",TRIM(filename)

  ! Grib Initializations
!DL NEW grib edition comes with NL
  IF (yform_write == 'api2') THEN
    iednr = 2
  ELSE
    iednr = 1
  ENDIF

!DL Distinguish between COSMO and ICON
  IF (ymodel == 'COSMO') THEN
    idim = ie
    jdim = je
    lds  = ie * je
  ELSE
    idim = nproma
    jdim = nblock
    lds  = ngp
  ENDIF
  ALLOCATE (ds(lds), ds2(idim,jdim), STAT=istat)

! Open file with grib_open_file for writing output foelds
     CALL grib_open_file (nudatc, TRIM(filename),'w  ',ierrc)

! Get the grib handle from sample and clone it - dependend on edition number iednr and ymodel                
  ierr = 0
  IF (iednr == 1) THEN
    CALL grib_new_from_samples (igrib1_id, 'DWD_rotated_ll_7km_G_grib1', ierr)
    IF (ierr /= GRIB_SUCCESS) THEN
      PRINT *,   ' *** Error in grib_new_from_sample (DWD_rotated_ll_7km_G_grib1): ', ierr
    ENDIF
    CALL grib_clone(igrib1_id,igribo, ierr)
    IF (ierr /= GRIB_SUCCESS) THEN
      PRINT *,   ' *** Error in grib_clone (GRIB1): ', ierr
    ENDIF
  ! Factor for startstep
    zfac=1./3600.
  ELSE
    IF (ymodel == 'COSMO') THEN
      CALL grib_new_from_samples (igrib2_id, 'DWD_rotated_ll_7km_G_grib2', ierr)
      IF (ierr /= GRIB_SUCCESS) THEN
        PRINT *,   ' *** Error in grib_new_from_sample (DWD_rotated_ll_7km_G_grib2): ', ierr
      ENDIF
    ELSE
      CALL grib_new_from_samples (igrib2_id, 'DWD_ICON_0026_R03B07_G_sfc', ierr)
      IF (ierr /= GRIB_SUCCESS) THEN
        PRINT *,   ' *** Error in grib_new_from_sample (DWD_ICON_0026_R03B07_G_sfc): ', ierr
      ENDIF
    ENDIF
    CALL grib_clone(igrib2_id,igribo, ierr)
    IF (ierr /= GRIB_SUCCESS) THEN
      PRINT *,   ' *** Error in grib_clone (GRIB2): ', ierr
    ENDIF
  ! Factor for startstep in hours- could be negative in GRIB2
    zfac=-1./3600.
  ENDIF

  IF (ierr /= 0) THEN
    WRITE(6,*) "Problems with GRIB sample or clone - STOP!!"
    STOP 
  ENDIF

!DL.............................................................................

! Define the grid - here rotated_ll or unstructured_grid (see sample above)

!..................COSMO.......................................................
IF (ymodel == 'COSMO') THEN

  ! no vertical coordinate parameters
  CALL grib_set (igribo,'deletePV', 1)

  CALL grib_set (igribo, 'Ni', ie)
  CALL grib_set (igribo, 'Nj', je)

  IF (iednr == 1) THEN
    CALL grib_set (igribo, 'resolutionAndComponentFlags', 8)
    CALL grib_set (igribo, 'Di', 0)
    CALL grib_set (igribo, 'Dj', 0)
  ELSE
    CALL grib_set (igribo, 'ijDirectionIncrementGiven',       1)
    CALL grib_set (igribo, 'uvRelativeToGrid',                1)
    CALL grib_set (igribo, 'iDirectionIncrementInDegrees', dlon)
    CALL grib_set (igribo, 'jDirectionIncrementInDegrees', dlat)
  ENDIF

  ! scanning mode
  CALL grib_set (igribo, 'scanningMode', 64)

  ! coordinates of the pole: in the grib code the southern pole has to be
  ! specified: for the latitude this is the negative value, for the
  ! longitude it is the value + 180.0 and then limiting the value again
  ! to -180.0 ... 180.0
  pollon_grib = pollon + 180.0_wp
  IF (pollon_grib > 180.0_wp) THEN
     pollon_grib = pollon_grib - 360.0_wp
  ENDIF

  CALL grib_set (igribo, 'latitudeOfSouthernPoleInDegrees', -pollat)
  CALL grib_set (igribo, 'longitudeOfSouthernPoleInDegrees', pollon_grib)
  CALL grib_set (igribo, 'angleOfRotationInDegrees',         polgam)

  ! start and end longitudes and latitudes
  CALL grib_set (igribo, 'latitudeOfFirstGridPointInDegrees',  startlat)
  CALL grib_set (igribo, 'longitudeOfFirstGridPointInDegrees', startlon)
  CALL grib_set (igribo, 'latitudeOfLastGridPointInDegrees',   startlat+(je-1)*dlat)
  CALL grib_set (igribo, 'longitudeOfLastGridPointInDegrees',  startlon+(ie-1)*dlon)
                            
!....................................ICON......................................
ELSE
  CALL grib_set (igribo, 'numberOfGridUsed',        ngrid)
  CALL grib_set (igribo, 'numberOfGridInReference', 1)
  CALL grib_set (igribo, 'uuidOfHGrid',             uuid)
ENDIF

!DL.............................................................................

  ! Product Definitions Sector
     CALL grib_set (igribo, 'centre','edzw')
     
     IF (ymodel == 'COSMO') THEN
       CALL grib_set (igribo, 'generatingProcessIdentifier', 131)
     ELSE
       CALL grib_set (igribo, 'generatingProcessIdentifier', ngpi)
     ENDIF 

!DL Reference time of data
  READ(yactdate1( 1: 8),'(I8)') nzdate
  READ(yactdate1( 9:10),'(I2)') nzhour
  CALL grib_set (igribo, 'dataDate',        nzdate)   ! yyyymmdd
  CALL grib_set (igribo, 'hour',            nzhour)   ! hh
  CALL grib_set (igribo, 'minute',        nzminute)   ! mm
  CALL grib_set (igribo, 'stepUnits',          'h')

!DL Additional meta information for GRIB2
  IF (iednr == 2) THEN
    CALL grib_set (igribo, 'significanceOfReferenceTime',     0)   ! Analysis
    CALL grib_set (igribo, 'productionStatusOfProcessedData', 2)   ! Research
    CALL grib_set (igribo, 'typeOfProcessedData',             2)   ! Analysis and forecast product
  ENDIF

!DL.............................................................................
  ! Some other GRIB metadata coming NOT with the sample
  CALL grib_set (igribo, 'bitsPerValue',16)
!DL.............................................................................
!DL.............................................................................

!------------------------------------------------------------------------------------
  ! Output of fields

  ! convert sums into averages
  astore(:,:,1)=astore(:,:,1)/FLOAT(nout_interval)  ! lhfl_s
  astore(:,:,2)=astore(:,:,2)/FLOAT(nout_interval)  ! shfl_s
  astore(:,:,5)=astore(:,:,5)/FLOAT(nout_interval)  ! sobs
  astore(:,:,6)=astore(:,:,6)/FLOAT(nout_interval)  ! thbs
   
  ! convert rates into totals
  astore(:,:,3)=astore(:,:,3)*dt                    ! prr_gsp
  astore(:,:,4)=astore(:,:,4)*dt                    ! prs_gsp
!DL ??
  astore(:,:,7)=astore(:,:,7)*dt                    ! runoff_s
  astore(:,:,8)=astore(:,:,8)*dt                    ! runoff_g
!DL ??
  
  ! GDM> produce YUCHKDAT with stats
  IF (lcheck) THEN
    OPEN(UNIT=444,FILE='YUCHKDAT',FORM='FORMATTED',STATUS='OLD',POSITION='APPEND')
    WRITE(444,*)
    WRITE(444,*) 'Writing ',TRIM(filename)
  ENDIF

  SELECT CASE (ntype_output)
  CASE (1)
     nvar=ke_soil+1
  CASE (2)
     nvar=2*(ke_soil+1)+4
  CASE (3)
     nvar=3*(ke_soil+1)+28
  END SELECT

  DO i=1,nvar

  !  reset some Keys
     CALL grib_set (igribo, 'stepType','instant')
     IF (iednr == 2) THEN
       CALL grib_set (igribo, 'shortName', 'T')
       CALL grib_set (igribo, 'productDefinitionTemplateNumber',0)
     ENDIF
     CALL grib_set (igribo, 'startStep',              0)   ! vv=0
     CALL grib_set (igribo, 'endStep',                0)   ! vv=0
     CALL grib_set (igribo, 'typeOfLevel',    'surface')
     CALL grib_set (igribo, 'level',                  0)
     CALL grib_set (igribo, 'topLevel',               0)
     CALL grib_set (igribo, 'bottomLevel',            0)
     IF (iednr == 2) THEN
       CALL grib_set_missing (igribo, 'scaleFactorOfFirstFixedSurface')
       CALL grib_set_missing (igribo, 'scaledValueOfFirstFixedSurface')
       CALL grib_set_missing (igribo, 'scaleFactorOfSecondFixedSurface')
       CALL grib_set_missing (igribo, 'scaledValueOfSecondFixedSurface')
     ENDIF

     IF (i<=ke_soil+1) THEN

        CALL grib_set (igribo, 'shortName', 'W_SO')
        IF (iednr == 1) THEN
           CALL grib_set (igribo, 'typeOfLevel', 'depthBelowLand')
           CALL grib_set (igribo, 'level', INT(czmls(i)*100.0+0.99))
        ELSE
           CALL grib_set (igribo, 'typeOfLevel', 'depthBelowLandLayer')
           IF (i==1) THEN
             CALL grib_set (igribo, 'topLevel',0)
             CALL grib_set (igribo, 'bottomLevel',czhls_const(i))
           ELSE
             CALL grib_set (igribo, 'topLevel',czhls_const(i-1))
             CALL grib_set (igribo, 'bottomLevel',czhls_const(i))
           ENDIF
        ENDIF
        ds2=REAL(w_so(:,:,i,nnew)*1000.0,irealgrib)         
        WRITE(str,'(A5,I4)') "w_so_",INT(czmls(i)*100.0+0.99)
        IF (lcheck) THEN
           CALL stats(str,w_so(:,:,i,nnew)*1000.0,llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==(ke_soil+1)+1) THEN

        CALL grib_set (igribo, 'shortName', 'T_SO')
        CALL grib_set (igribo, 'typeOfLevel', 'depthBelowLand')
        CALL grib_set (igribo, 'level', 0)
        ds2=REAL(t_so(:,:,0,nnew),irealgrib)   
        IF (lcheck) THEN
           CALL stats("t_so_0",t_so(:,:,0,nnew),llandmask,UNIT=444)
        ENDIF
      ELSE IF (i<=2*(ke_soil+1)+1) THEN
        ii=i-(ke_soil+2)
        CALL grib_set (igribo, 'shortName', 'T_SO')
        CALL grib_set (igribo, 'typeOfLevel', 'depthBelowLand')
        IF (iednr == 1) THEN
          CALL grib_set (igribo, 'level', INT(czmls(ii)*100.0+0.99))
        ELSE
          IF(ii==1) THEN
            CALL grib_set (igribo, 'scaleFactorOfFirstFixedSurface',3)
            CALL grib_set (igribo, 'scaledValueOfFirstFixedSurface',5)
          ELSE
            CALL grib_set (igribo, 'scaleFactorOfFirstFixedSurface',2)
            CALL grib_set (igribo, 'scaledValueOfFirstFixedSurface',INT(czmls(ii)*100.0+0.99))
          ENDIF
        ENDIF
        ds2=REAL(t_so(:,:,ii,nnew),irealgrib)   
        WRITE(str,'(A5,I4)') "t_so_",INT(czmls(ii)*100.0+0.99)
        IF (lcheck) THEN
           CALL stats(str,t_so(:,:,ii,nnew),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+2) THEN
        CALL grib_set (igribo, 'shortName', 'T_SNOW')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        ds2=REAL(t_snow(:,:,nnew),irealgrib)   
        IF (lcheck) THEN
           CALL stats("t_snow",t_snow(:,:,nnew),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+3) THEN
        CALL grib_set (igribo, 'shortName', 'W_SNOW')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        ds2=REAL(w_snow(:,:,nnew)*1000.0,irealgrib)   
        IF (lcheck) THEN
           CALL stats("w_snow",w_snow(:,:,nnew)*1000.0,llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+4) THEN
        CALL grib_set (igribo, 'shortName', 'W_I')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        ds2=REAL(w_i(:,:,nnew)*1000.0,irealgrib)   
        IF (lcheck) THEN
           CALL stats("w_i",w_i(:,:,nnew)*1000.0,llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+5) THEN
! GDM> ke instead of 30m
        CALL grib_set (igribo, 'shortName', 'T')
        IF (iednr == 1) THEN
          CALL grib_set (igribo, 'typeOfLevel', 'hybridLayer')
        ELSE
           CALL grib_set (igribo, 'typeOfLevel', 'generalVerticalLayer')
        ENDIF
        CALL grib_set (igribo, 'topLevel',ke_model)
        CALL grib_set (igribo, 'bottomLevel',ke_model+1)
        ds2=REAL(t(:,:,ke,nnew),irealgrib)   
        WRITE(str,'(A2,I2)') "t_k",ke_model
        IF (lcheck) THEN
           CALL stats(str,t(:,:,ke,nnew),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+6) THEN
        CALL grib_set (igribo, 'shortName', 'QV')
        IF (iednr == 1) THEN
          CALL grib_set (igribo, 'typeOfLevel', 'hybridLayer')
        ELSE
           CALL grib_set (igribo, 'typeOfLevel', 'generalVerticalLayer')
        ENDIF
        CALL grib_set (igribo, 'topLevel',ke_model)
        CALL grib_set (igribo, 'bottomLevel',ke_model+1)
        ds2=REAL(qv(:,:,ke,nnew),irealgrib)
!       ds2=REAL(qv(:,:,ke),irealgrib)   
     WRITE(str,'(A3,I2)') "qv_k",ke_model
        IF (lcheck) THEN
           CALL stats(str,qv(:,:,ke,nnew),llandmask,UNIT=444)
        ENDIF
!XYZ> qv in TSA is 4 dimensional
!       IF (lcheck) THEN
!          CALL stats(str,qv(:,:,ke),llandmask,UNIT=444)
!       ENDIF
!XYZ<
     ELSE IF (i==2*(ke_soil+1)+7) THEN
        CALL grib_set (igribo, 'shortName', 'SP')
        IF (iednr == 1) THEN
          CALL grib_set (igribo, 'typeOfLevel', 'hybridLayer')
        ELSE
           CALL grib_set (igribo, 'typeOfLevel', 'generalVerticalLayer')
        ENDIF
        CALL grib_set (igribo, 'topLevel',ke_model)
        CALL grib_set (igribo, 'bottomLevel',ke_model+1)
        ds2=REAL(u(:,:,ke,nnew),irealgrib)   
        WRITE(str,'(A2,I2)') "u_k",ke_model
        IF (lcheck) THEN
           CALL stats(str,u(:,:,ke,nnew),llandmask,UNIT=444)
        ENDIF
!GDM<
     ELSE IF (i==2*(ke_soil+1)+8) THEN
        CALL grib_set (igribo, 'shortName', 'PS')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        ds2=REAL(ps(:,:,nnew),irealgrib)   
        IF (lcheck) THEN
           CALL stats("ps",ps(:,:,nnew),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+9) THEN
        IF (iednr == 2) CALL grib_set (igribo, 'productDefinitionTemplateNumber',8)
        CALL grib_set (igribo, 'shortName', 'TOT_RAIN')
        IF (iednr == 1) CALL grib_set (igribo, 'timeRangeIndicator',4)
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        CALL grib_set (igribo, 'endStep',0)
        CALL grib_set (igribo, 'startStep', int(zfac * dt * real(nout_interval,wp)))
        ds2=REAL(astore(:,:,3),irealgrib)         ! prr/rain_gsp
        IF (lcheck) THEN
           CALL stats("prr_gsp",astore(:,:,3),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+10) THEN
        IF (iednr == 2) CALL grib_set (igribo, 'productDefinitionTemplateNumber',8)
        CALL grib_set (igribo, 'shortName', 'TOT_SNOW')
        IF (iednr == 1) CALL grib_set (igribo, 'timeRangeIndicator',4)
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        CALL grib_set (igribo, 'endStep',0)
        CALL grib_set (igribo, 'startStep',  int(zfac * dt * real(nout_interval,wp)))
        ds2=REAL(astore(:,:,4),irealgrib)         ! prs/snow_gsp
        IF (lcheck) THEN
           CALL stats("prs_gsp",astore(:,:,4),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+11) THEN
        IF (iednr == 2) CALL grib_set (igribo, 'productDefinitionTemplateNumber',8)
        CALL grib_set (igribo, 'shortName', 'ASOB_S')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        CALL grib_set (igribo, 'endStep',0)
        CALL grib_set (igribo, 'startStep',  int(zfac * dt * real(nout_interval,wp)))
        ds2=REAL(astore(:,:,5),irealgrib)         ! sobs_s
        IF (lcheck) THEN
           CALL stats("sobs_s",astore(:,:,5),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+12) THEN
        IF (iednr == 2) CALL grib_set (igribo, 'productDefinitionTemplateNumber',8)
        CALL grib_set (igribo, 'shortName', 'ATHB_S')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        CALL grib_set (igribo, 'endStep',0)
        CALL grib_set (igribo, 'startStep',  int(zfac * dt * real(nout_interval,wp)))
        ds2=REAL(astore(:,:,6),irealgrib)         ! thbs_s
        IF (lcheck) THEN
            CALL stats("thbs_s",astore(:,:,6),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+13) THEN
        IF (iednr == 2) CALL grib_set (igribo, 'productDefinitionTemplateNumber',8)
        CALL grib_set (igribo, 'shortName', 'ALHFL_S')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        CALL grib_set (igribo, 'endStep',0)
        CALL grib_set (igribo, 'startStep',  int(zfac * dt * real(nout_interval,wp)))
        ds2=REAL(astore(:,:,1),irealgrib)         ! lhfl_s
        IF (lcheck) THEN
           CALL stats("lhfl_s",astore(:,:,1),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+14) THEN
        IF (iednr == 2) CALL grib_set (igribo, 'productDefinitionTemplateNumber',8)
        CALL grib_set (igribo, 'shortName', 'ASHFL_S')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        CALL grib_set (igribo, 'endStep',0)
        CALL grib_set (igribo, 'startStep',  int(zfac * dt * real(nout_interval,wp)))
        ds2=REAL(astore(:,:,2),irealgrib)         ! shfl_s
        IF (lcheck) THEN
            CALL stats("shfl_s",astore(:,:,2),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+15) THEN
!DL level setting has to be checked for GRIB1/2 
        IF (iednr == 2) CALL grib_set (igribo, 'productDefinitionTemplateNumber',8)
        CALL grib_set (igribo, 'shortName', 'RUNOFF_S')
        IF (iednr == 1) CALL grib_set (igribo, 'timeRangeIndicator',4)
        CALL grib_set (igribo, 'typeOfLevel', 'depthBelowLandLayer')
        CALL grib_set (igribo, 'topLevel',0.)
        CALL grib_set (igribo, 'bottomLevel', 0.1)
        CALL grib_set (igribo, 'endStep',0)
        CALL grib_set (igribo, 'startStep',  int(zfac * dt * real(nout_interval,wp)))
        ds2=REAL(astore(:,:,7),irealgrib)         ! runoff_s
        IF (lcheck) THEN
           CALL stats("runoff_s",astore(:,:,7),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+16) THEN
        IF (iednr == 2) CALL grib_set (igribo, 'productDefinitionTemplateNumber',8)
        CALL grib_set (igribo, 'shortName', 'RUNOFF_G')
        IF (iednr == 1) CALL grib_set (igribo, 'timeRangeIndicator',4)
        CALL grib_set (igribo, 'typeOfLevel', 'depthBelowLandLayer')
        CALL grib_set (igribo, 'topLevel',0.1)
        CALL grib_set (igribo, 'bottomLevel', 1.0)
        CALL grib_set (igribo, 'endStep',0)
        CALL grib_set (igribo, 'startStep',  int(zfac * dt * real(nout_interval,wp)))
        ds2=REAL(astore(:,:,8),irealgrib)         ! runoff_g
        IF (lcheck) THEN
           CALL stats("runoff_g",astore(:,:,8),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+17) THEN
        CALL grib_set (igribo, 'shortName', 'TCH')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        ds2=REAL(tch(:,:),irealgrib)             ! tch
        IF (lcheck) THEN
           CALL stats("tch",tch(:,:),llandmask,UNIT=444)
        ENDIF
!GDM>  store  
     ELSE IF (i==2*(ke_soil+1)+18) THEN
        CALL grib_set (igribo, 'shortName', 'RSTOM')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        ds2=REAL(store(:,:,1),irealgrib)         ! rstom
        IF (lcheck) THEN
           CALL stats("rstom",store(:,:,1),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+19) THEN
        CALL grib_set (igribo, 'shortName', 'DUMMY_251')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        ds2=REAL(store(:,:,2),irealgrib)         ! zf_wat
        IF (lcheck) THEN
           CALL stats("zf_wat",store(:,:,2),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+20) THEN
        CALL grib_set (igribo, 'shortName', 'DUMMY_252')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        ds2=REAL(store(:,:,3),irealgrib)         ! zf_tem
        IF (lcheck) THEN
           CALL stats("zf_tem",store(:,:,3),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+21) THEN
        CALL grib_set (igribo, 'shortName', 'DUMMY_253')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        ds2=REAL(store(:,:,4),irealgrib)         ! zf_rad
        IF (lcheck) THEN
           CALL stats("zf_rad",store(:,:,4),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2*(ke_soil+1)+22) THEN
        CALL grib_set (igribo, 'shortName', 'DUMMY_254')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        ds2=REAL(store(:,:,5),irealgrib)         ! zep_s 
        IF (lcheck) THEN
           CALL stats("zep_s",store(:,:,5),llandmask,UNIT=444)
        ENDIF
! GDM< store
     ELSE IF (i<=3*(ke_soil+1)+22) THEN
        ii=i-(2*(ke_soil+1)+22)
        IF (iednr == 1) THEN
           CALL grib_set (igribo, 'shortName', 'RUNOFF')
           CALL grib_set (igribo, 'typeOfLevel', 'depthBelowLand')
           CALL grib_set (igribo, 'level', INT(czmls(ii)*100.0+0.99))
        ELSE
!DL as long as <ro> is not available (workstation)
           CALL grib_set (igribo, 'discipline',2)
           CALL grib_set (igribo, 'parameterCategory',0)
           CALL grib_set (igribo, 'parameterNumber',5)
           CALL grib_set (igribo, 'typeOfLevel', 'depthBelowLandLayer')
           IF (ii==1) THEN
             CALL grib_set (igribo, 'topLevel',0)
             CALL grib_set (igribo, 'bottomLevel',czhls_const(ii))
           ELSE
             CALL grib_set (igribo, 'topLevel',czhls_const(ii-1))
             CALL grib_set (igribo, 'bottomLevel',czhls_const(ii))
           ENDIF
        ENDIF
        ds2=REAL(runoff(:,:,ii),irealgrib)         
        WRITE(str,'(A7,I4)') "runoff_",INT(czmls(ii)*100.0+0.99)
        IF (lcheck) THEN
           CALL stats(str,runoff(:,:,ii),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==3*(ke_soil+1)+23) THEN
        CALL grib_set (igribo, 'shortName', 'T_2M')
        ds2=REAL(t_2m,irealgrib)         
        IF (lcheck) THEN
           CALL stats("t_2m",t_2m,llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==3*(ke_soil+1)+24) THEN
        CALL grib_set (igribo, 'shortName', 'RHO_SNOW')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        ds2=REAL(rho_snow(:,:,nnew),irealgrib)   
        IF (lcheck) THEN
           CALL stats("rho_snow",rho_snow(:,:,nnew),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==3*(ke_soil+1)+25) THEN
        CALL grib_set (igribo, 'shortName', 'FRESHSNW')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        ds2=REAL(freshsnow(:,:),irealgrib)   
        IF (lcheck) THEN
           CALL stats("freshsnow",freshsnow(:,:),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==3*(ke_soil+1)+26) THEN
        CALL grib_set (igribo, 'shortName', 'TCM')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        ds2=REAL(tcm(:,:),irealgrib)   
        IF (lcheck) THEN
          CALL stats("tcm",tcm(:,:),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==3*(ke_soil+1)+27) THEN
        CALL grib_set (igribo, 'shortName', 'QV_S')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        ds2=REAL(qv_s(:,:,nnew),irealgrib)   
        IF (lcheck) THEN
          CALL stats("qv_s",qv_s(:,:,nnew),llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==3*(ke_soil+1)+28) THEN
        CALL grib_set (igribo, 'shortName', 'T_G')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        ds2=REAL(t_g(:,:,nnew),irealgrib)   
        IF (lcheck) THEN
          CALL stats("t_g",t_g(:,:,nnew),llandmask,UNIT=444)
        ENDIF
     ENDIF

!DL NEW
     ! undef handling for non-land points
     lbitmap=.FALSE.
     !  rmin=MINVAL(ds2,llandmask)
     WHERE (.NOT.llandmask)
       ds2(:,:)=undefgrib !rmin
     END WHERE
     IF (ANY (ds2 == undefgrib)) lbitmap=.TRUE.
!DL NEW create bitmap
     IF (lbitmap) THEN
       CALL grib_set (igribo, 'missingValue',undefgrib)
       CALL grib_set (igribo, 'bitmapPresent',1)
     ENDIF
!DL NEW

     ! write data
     ds=RESHAPE(ds2,(/ idim*jdim /))

     CALL grib_set (igribo, 'values', ds(1:lds), ierr)
     IF (ierr /= GRIB_SUCCESS) THEN
       WRITE(6,*) "Error in grib_set (values) for field ",i," of processed fields"
     ENDIF

     ! Write the encoded fields to file
     CALL grib_write (igribo, nudatc)

  ENDDO
    

  ! write constant fields
  nvar=0
  IF (lconst) THEN
     nvar=nvar+9
  ENDIF
  DO i=1,nvar

     IF (i==1) THEN
        CALL grib_set (igribo, 'shortName', 'LAI')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        ds2=REAL(lai(:,:),irealgrib)   
        IF (lcheck) THEN
          CALL stats("lai_out",lai,llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==2) THEN
        CALL grib_set (igribo, 'shortName', 'PLCOV')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        ds2=REAL(plcov(:,:),irealgrib)   
        IF (lcheck) THEN
          CALL stats("plcov_out",plcov,llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==3) THEN
        CALL grib_set (igribo, 'shortName', 'ROOTDP')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        ds2=REAL(rootdp(:,:),irealgrib)   
        IF (lcheck) THEN
          CALL stats("rootdp_out",rootdp,llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==4) THEN
        CALL grib_set (igribo, 'shortName', 'SOILTYP')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        ds2=REAL(soiltyp(:,:),irealgrib)   
        IF (lcheck) THEN
          CALL stats("soiltyp_out",soiltyp,llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==5) THEN
        CALL grib_set (igribo, 'shortName', 'HSURF')
!DL all the level coding is defined in shortName --> no call for typeOfLevel
        ds2=REAL(hsurf(:,:),irealgrib)   
        IF (lcheck) THEN
          CALL stats("hsurf_out",hsurf,llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==6) THEN
        CALL grib_set (igribo, 'shortName', 'Z0')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        ds2=REAL(z0(:,:),irealgrib)   
        IF (lcheck) THEN
          CALL stats("z0_out",z0,llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==7) THEN
        CALL grib_set (igribo, 'shortName', 'FR_LAND')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        ds2=REAL(fr_land(:,:),irealgrib)   
        IF (lcheck) THEN
          CALL stats("fr_land_out",fr_land,llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==8) THEN
        CALL grib_set (igribo, 'shortName', 'RLON')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        ds2=REAL(rlon(:,:),irealgrib)   
        IF (lcheck) THEN
          CALL stats("rlon_out",rlon,llandmask,UNIT=444)
        ENDIF
     ELSE IF (i==9) THEN
        CALL grib_set (igribo, 'shortName', 'RLAT')
        CALL grib_set (igribo, 'typeOfLevel', 'surface')
        ds2=REAL(rlat(:,:),irealgrib)   
        IF (lcheck) THEN
          CALL stats("rlat_out",rlat,llandmask,UNIT=444)
        ENDIF
     ENDIF

     ! undef handling for non-land points
     lbitmap=.FALSE.
     IF (i<=3) THEN
     !  rmin=MINVAL(ds2,llandmask)
        WHERE (.NOT.llandmask)
          ds2(:,:)=undefgrib !rmin
        END WHERE
        IF (ANY (ds2 == undefgrib)) lbitmap=.TRUE.
     ENDIF

!DL NEW create bitmap
     IF (lbitmap) THEN
       CALL grib_set (igribo, 'missingValue',undefgrib)
       CALL grib_set (igribo, 'bitmapPresent',1)
     ENDIF
!DL
     ! write data
     ds=RESHAPE(ds2,(/ idim*jdim /))

     CALL grib_set (igribo, 'values', ds(1:lds), ierr)
     IF (ierr /= GRIB_SUCCESS) THEN
       WRITE(6,*) "Error in grib_set (values) for field ",i," of constant fields"
     ENDIF
     ! Write the encoded fields to file
     CALL grib_write (igribo, nudatc)

  ENDDO

!DL .........................................................................................
  ! Release the grib handle
  CALL grib_release (igribo)
  CALL grib_release (igrib1_id)
  CALL grib_release (igrib2_id)

  ! Close the grib file
  CALL grib_close_file (nudatc)
!DL .........................................................................................

  ! reset summation storage
  astore=0.0

  ! YUCHKDAT
  IF (lcheck) THEN
    CLOSE(UNIT=444)
  ENDIF

  DEALLOCATE (ds, ds2)

END SUBROUTINE grbout_eccodes
!----------------------------------------------------------------------

SUBROUTINE binout(nt)
!----------------------------------------------------------------------
! Description:
!
! Writes output as BIN file
!----------------------------------------------------------------------

  IMPLICIT NONE
  ! Subroutine Arguments:
  INTEGER                 , INTENT(IN) :: nt
  ! Local Variables
  INTEGER                              :: &
      noff                              , & !
      nvar                              , & !
      i,j                                   ! Looping indices

  nvar = 3*(ke_soil+1)+22
  IF (nt>=0) THEN
    noff = nt/nout_interval*nvar
    OPEN(91,FILE=TRIM(outdir)//'/'//TRIM(outprefix),ACCESS='direct',RECL=ie*je*4)
  ELSE
    OPEN(91,FILE=TRIM(outdir)//'/'//TRIM(outprefix)//'_init',ACCESS='direct',RECL=ie*je*4)
  ENDIF

  DO i=1,ke_soil+1
     WRITE(91,REC=noff+i) REAL(w_so(:,:,i,nnew),4)         
  ENDDO

  DO i=1,ke_soil+1
     WRITE(91,REC=noff+i+ke_soil+1) REAL(t_so(:,:,i,nnew),4)         
  ENDDO

  WRITE(91,REC=noff+2*(ke_soil+1)+1) REAL(t_snow(:,:,nnew),4)   
  WRITE(91,REC=noff+2*(ke_soil+1)+2) REAL(w_snow(:,:,nnew),4)   
  WRITE(91,REC=noff+2*(ke_soil+1)+3) REAL(w_i(:,:,nnew),4)   
  WRITE(91,REC=noff+2*(ke_soil+1)+4) REAL(t(:,:,ke,nnew),4)   
  WRITE(91,REC=noff+2*(ke_soil+1)+5) REAL(qv(:,:,ke,nnew),4)
! WRITE(91,REC=noff+2*(ke_soil+1)+5) REAL(qv(:,:,ke),4)       !XYZ> qv in TSA is 4 dimensional
  WRITE(91,REC=noff+2*(ke_soil+1)+6) REAL(u(:,:,ke,nnew),4)   
! WRITE(91,REC=noff+2*(ke_soil+1)+6) REAL(u_10m(:,:),4)   
  WRITE(91,REC=noff+2*(ke_soil+1)+7) REAL(ps(:,:,nnew),4)

  astore=astore/FLOAT(nout_interval)   
  DO i=1,8
     WRITE(91,REC=noff+2*(ke_soil+1)+7+i) REAL(astore(:,:,i),4)   
  ENDDO

  WRITE(91,REC=noff+2*(ke_soil+1)+16) REAL(tch(:,:),4) 

  DO i=1,5
     WRITE(91,REC=noff+2*(ke_soil+1)+16+i) REAL(store(:,:,i),4)  
  ENDDO

  DO i=1,ke_soil+1
     WRITE(91,REC=noff+2*(ke_soil+1)+21+i) REAL(runoff(:,:,i),4)         
  ENDDO

  WRITE(91,REC=noff+3*(ke_soil+1)+22) REAL(t_2m(:,:),4)         
  CLOSE(91)
  ! reset summation storage
  astore=0.0
END SUBROUTINE binout
!----------------------------------------------------------------------


SUBROUTINE ascout(nt)
!----------------------------------------------------------------------
! Description:
!
! Writes output as ASCII file
!----------------------------------------------------------------------

  IMPLICIT NONE
  ! Subroutine Arguments:
  INTEGER                 , INTENT(IN) :: nt
  ! Local Variables
  INTEGER                              :: i

  IF (nt<0) THEN
     OPEN(91,FILE=TRIM(outdir)//'/'//TRIM(outprefix)//'_init',STATUS='replace')
  ELSEIF (nt==0) THEN
     OPEN(91,FILE=TRIM(outdir)//'/'//TRIM(outprefix),STATUS='replace')
  ELSE
     OPEN(91,FILE=TRIM(outdir)//'/'//TRIM(outprefix),POSITION='append')
  ENDIF

  WRITE(91,'(I6,1X)',ADVANCE='no') nt
  WRITE(91,'(F9.4,1X)',ADVANCE='no') w_so(:,:,1,nnew)*100 / REAL(czhls_const(1))
     
  DO i=2,ke_soil+1
     WRITE(91,'(F9.4,1X)',ADVANCE='no') &
          w_so(:,:,i,nnew)*100 / REAL(czhls_const(i)- czhls_const(i-1))        
  ENDDO

  DO i=1,ke_soil+1
     WRITE(91,'(F9.4,1X)',ADVANCE='no') t_so(:,:,i,nnew)
  ENDDO


  write(91,'(F9.4,1X)',ADVANCE='no') rho_snow(:,:,nnew)
  write(91,'(F9.4,1X)',ADVANCE='no') w_snow(:,:,nnew)
  write(91,'(F9.4,1X)',ADVANCE='no') t_snow(:,:,nnew)
  write(91,'(F9.4,1X)',ADVANCE='no') h_snow(:,:,nnew)


  WRITE(91,*)
  CLOSE(91)

END SUBROUTINE ascout
!----------------------------------------------------------------------

!----------------------------------------------------------------------

END MODULE tsa_output
