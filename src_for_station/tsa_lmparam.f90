!+ Subroutines and functions which correspond to the COSMO-Parameterisations
!==============================================================================

MODULE tsa_lmparam

!==============================================================================
!
! Description:
!  This file contains subroutines that correspond to COSMO-Parameterisations
!  and allow to run as Stand Alone version.
!  This file is included in terra when compiling
!
!  Currently included:
!
!    - vegadapt:
!      vegetation parameters calculations for TSA
!
!    - qvsat (Function):
!      calculates saturation vapor quantity
!
!    - parturs:
!      calculates transfer coefficients (tcm,tch)
!      according to louis scheme.
!      substitutes for TURB_DIFF in COSMO
!
!    - calc_albedo:
!      calculates solar and thermal albedos for the surface
!
!    - gen_area_indices:
!      generates area indices - sai,eai & tai according to lai
!
!    - near_surface (from near_surface of COSMO)
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
!  Revised code to correspond to COSMO version 5.01
!  Arranged code to adhere to coding standards.
! V5.3       2015-04-01 Yiftach Ziv, IMS (XYZ)
!  Added subroutine parturs_new: a more stable version of partur.
!   this subroutine was developed and tested by Julian Todter (JT)
!   and colleagues in GUF.
! V5.7_SNOW 2021-01-08 Varun Sharma
!   albedo forced to be fixed @ 0.9
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
    ymodel, rundef, z0, dz, dz_u, lhomosoil, lrootadapt, qv, vegalb

USE data_constants,        ONLY:                                            &
    g, pi, cpdr, rvd_m_o, cp_d, t0_melt

USE data_fields,           ONLY:                                            &
    tcm, tch, t_g, llandmask, u_10m, v_10m, w_so, w_snow, t_snow, plcov,    &
    lai, eai, sai, tai, u, v, t, p0, pp, t_2m, qv_s, soiltyp, hsurf,        &
    rlat

USE data_modelconfig,      ONLY:                                            &
    ie, je, ke, istartpar, iendpar, jstartpar, jendpar, czmls

USE data_runcontrol,       ONLY:                                            &
    nblock, nproma, nlastproma, nnow

USE sfc_terra_data,        ONLY:                                            &
    cf_snow, csalb_snow, ctalb, csalb, csalbw

USE turb_data,             ONLY:                                            &
    akt

!==============================================================================

IMPLICIT NONE

!==============================================================================

CONTAINS

!==============================================================================

SUBROUTINE vegadapt(jahrestag,hsurface,lat,lai_min,lai_max,         &
                    plcov_min,plcov_max,rootdp_in,lai,plcov,rootdp)


!-------------------------------------------------------------------------------
!
! Description:
!   calculates Leaf Area Index (lai), PLant COVer (plcov) and ROOT DePth (rootdp)
!   based on annual cycle for the standalone version
!-------------------------------------------------------------------------------

! Subroutine arguments:
  REAL (KIND=wp), INTENT(IN)            :: &
       hsurface                          , & ! Surface height
       lat                               , & ! Latitude
       lai_min(:,:)                      , & ! minimal LAI
       lai_max(:,:)                      , & ! maximal LAI
       plcov_min(:,:)                    , & ! minimal plant cover
       plcov_max(:,:)                    , & ! maximal plant cover
       rootdp_in(:,:)
  REAL (KIND=wp), INTENT(OUT)           :: &
       lai(:,:)                          , & ! actual calculated LAI
       plcov(:,:)                        , & ! actual calculated plant cover
       rootdp(:,:)                           ! actual calculated root depth
  INTEGER, INTENT(IN)                   :: &
       jahrestag                             ! 

  REAL (KIND=wp)                        :: &
       zhred                             , & !
       zvegfac                           , & !
       zbvp                              , & !
       zdvp                                  !

  INTEGER                               :: & !
       i, j

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE vegadapt
!------------------------------------------------------------------------------

!DL Take care for ICON array structure
  IF (ymodel /= 'COSMO') THEN
    ie = nproma
    je = nblock
  ENDIF

  DO i=1,ie
     DO j=1,je
        IF (lhomosoil) THEN
           IF ((i==1).AND.(j==1)) THEN
              zhred=EXP(-5.0e-9_wp* (hsurface*g)**2)
              zbvp=MAX(1.0_wp,3.0*(ABS(lat)-20.0)) 
              zdvp=MIN(365.0_wp,345.0_wp-4.5*(ABS(lat)-20.0))
           ENDIF
        ELSE
           zhred=EXP(-5.0e-9_wp* (hsurf(i,j)*g)**2)
           zbvp=MAX(1.0_wp,3.0*(ABS(rlat(i,j))-20.0)) 
           zdvp=MIN(365.0_wp,345.0_wp-4.5*(ABS(rlat(i,j))-20.0))
        ENDIF
        IF ((jahrestag<zbvp).OR.(jahrestag>zbvp+zdvp)) THEN
           zvegfac=0.0
        ELSE
           zvegfac=MAX(0.0_wp,MIN(1.0_wp, &
                   1.12_wp*SIN(pi*MAX(0.0_wp,(FLOAT(jahrestag)-zbvp)/zdvp))*zhred))
!                   1.12_wp*SIN(3.14159*MAX(0.0_wp,(FLOAT(jahrestag)-zbvp)/zdvp))*zhred))
        ENDIF
        lai(i,j)   =lai_min(i,j) + zvegfac * (lai_max(i,j) - lai_min(i,j))
        plcov(i,j) =plcov_min(i,j) +  zvegfac * (plcov_max(i,j) - plcov_min(i,j))
        IF (lrootadapt) THEN
           rootdp(i,j)=MIN(rootdp_in(i,j) , 0.12_wp+zvegfac**2*0.58_wp)
        ENDIF
     ENDDO
  ENDDO
END SUBROUTINE vegadapt

! -------------------------------------------------------------------


REAL (KIND=wp) FUNCTION qvsat(t,p)
!-------------------------------------------------------------------------------
!
! Description:
!   calculates saturation vapor quantity from temperature and pressure
!-------------------------------------------------------------------------------

! Function arguments:
  REAL (KIND=wp), INTENT(IN) :: &
       t                      , & ! air temperature
       p                          ! air pressure

! Local variable
  REAL (KIND=wp)             :: esat ! saturation vapor pressure

! Begin Function qvsat
  esat=610.78*EXP((17.1*(t-273.15))/(235.0+t-273.15))
  qvsat=0.622*esat/(p-0.378*esat)

END FUNCTION qvsat

!--------------------------------------------------------


SUBROUTINE  parturs

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure calculates the transfer coefficients for momentum
!   (TCM) and heat and moisture (TCH) at the surface as well as the roughness
!   length over sea (GZ0).
!
! Method:
!   Dyer-Businger equations as modified by Louis.
!
!------------------------------------------------------------------------------

! Subroutine arguments: None
! --------------------
!
! Local parameters:
! ----------------
  REAL    (KIND=wp   ), PARAMETER ::  &
    ! basic constants of the parameterization scheme and meteorological and
    ! numerical parameters
    zah     = 5.300_wp,    & !
    zgz0hh  = 0.980_wp,    & ! upper limit for roughness length for heat
    zalphaf = 1.000_wp,    & !
    zalpha0 = 0.0150_wp,   & ! Charnock constant for momentum
    zalphah = 0.60_wp  ,   & !
    zviscos = 1.5E-5   ,   & ! kinematic viscosity constant (m**2/s)
    zbeta10 = 0.0420_wp,   & !
    z10     = 10.0_wp  ,   & !
    zvmin   = 0.01_wp  ,   & ! Minimum wind velocity
    ztmmin  = 0.140_wp ,   & ! Minimum value for TCM
    zthmin  = 0.010_wp ,   & ! Minimum value for TCH
    zed3    = 1.0_wp/3.0_wp,  & !
    zd3     = 2.0_wp/3.0_wp,  & !
    zctke   = 0.516_wp

! Local scalars:
! -------------
  INTEGER          ::  &
    i     ,            & ! loop index in x-direction              
    j     ,            & ! loop index in y-direction              
    nx    ,            & ! time-level for computation
    izg   ,            & !
    jzg   ,            & !
    mzg   ,            & !
    nzpa  ,            & !
!FEA
    im1   ,            & !
    jm1                  !
!FEE

  REAL (KIND=wp)   ::  &
    zum   ,            & ! 
    zvm   ,            & ! 
    ztvg  ,            & ! 
    ztvke ,            & ! 
    zxi   ,            & ! 
    zxih  ,            & ! 
    zy    ,            & !
    zgz0d ,            & !
    zgz0dd,            & !
    zustar,            & !
    zdz   ,            & !
!FEA
    psi_m_dz,          & !
    psi_m_dzu,         & !
    x     ,            & !
    x0    ,            & !
    u_temp,            & !
    dtemp ,            & !
    l_monin              !
!FEE

! Local (automatic) arrays:
! -------------------------
  REAL    (KIND=wp   ) ::  &
    zvpb    (ie,je),   & !  
    zx      (ie,je),   & ! 
    ztcm    (ie,je),   & ! 
    ztch    (ie,je),   & ! 
    zdfip   (ie,je),   & ! 
    zris    (ie,je),   & !
    zgz0m   (ie,je),   & ! Roughness length for momentum
    zgz0h   (ie,je)      ! Roughness length for scalars (heat, moisture)

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine parturs             
!------------------------------------------------------------------------------
 
 nx = nnow

!------------------------------------------------------------------------------
! Section 1: Averaging the horizontal velocities (C-grid) at lowest layer
!------------------------------------------------------------------------------

  ! In the interior of the subdomain and the northern and eastern boundary of 
  ! the total domain
  DO j = jstartpar , jendpar
     jm1  = MAX( 1, j-1 )
   DO i = istartpar , iendpar
      im1  = MAX( 1, i-1 )
!FEA
      IF (llandmask(i,j)) THEN

      ztvg       = t_g(i,j,   nx)*(1.0 + rvd_m_o*qv_s(i,j,   nx))
      ztvke      = t  (i,j,ke,nx)*(1.0 + rvd_m_o*qv  (i,j,ke,nx))
!     ztvke      = t  (i,j,ke,nx)*(1.0 + rvd_m_o*qv  (i,j,ke))

!FEA wind adaption
      im1  = MAX( 1, i-1 )
      zum  = 0.5*(  u(i,j  ,ke,nx) + u(im1,j  ,ke,nx) )
      zvm  = 0.5*(  v(i,j  ,ke,nx) + v(i  ,jm1,ke,nx) )
      zvpb(i,j) =  SQRT(zum**2 + zvm**2)
      IF (dz_u/=dz) THEN
         dtemp=ztvke - ztvg + cpdr*zdfip(i,j) 
         u_temp  =  zvpb(i,j)*LOG(dz/z0(i,j))/LOG(dz_u/z0(i,j))
         IF (( zvpb(i,j)==0).OR.( dtemp==0)) THEN 
            zvpb(i,j)  =  MAX( u_temp, zvmin ) 
         ELSE
            l_monin=tcm(i,j)/tch(i,j) * &
                ztvg/(akt*akt*g*u_temp*dtemp)
            IF (l_monin>0) THEN
               psi_m_dz  = 4.7*(dz/l_monin-z0(i,j)/l_monin)
               psi_m_dzu = 4.7*(dz_u/l_monin-z0(i,j)/l_monin)
            ELSE
               x=EXP(0.25*LOG(1-15*dz/l_monin))
               x0=EXP(0.25*LOG(1-15*z0(i,j)/l_monin))
               psi_m_dz  = -2*LOG((1.0+x)/(1.0+x0))-LOG((1.0+x*x)/(1.0+x0*x0)) + &
                            2*ATAN(x)-2*ATAN(x0)
               x= EXP(0.25*LOG(1-15*dz_u/l_monin))    
               psi_m_dzu = -2*LOG((1.0+x)/(1.0+x0))-LOG((1.0+x*x)/(1.0+x0*x0)) + &
                            2*ATAN(x)-2*ATAN(x0)              
            ENDIF
            zvpb(i,j)  =  MAX(zvmin, &
                          zvpb(i,j)*(LOG(dz/z0(i,j))+psi_m_dz)/(LOG(dz_u/z0(i,j))+psi_m_dzu))
         ENDIF
      ENDIF
!FEE

!FEA
!      zdfip(i,j) = 0.5*g*( hhl(i,j,ke) - hhl(i,j,ke1) ) 
      zdfip(i,j) = g*dz
!FEE
      zx   (i,j) = ( ztvke - ztvg + cpdr*zdfip(i,j) )*zdfip(i,j)/t_g(i,j,nx)
      zris (i,j) = zx(i,j)/zvpb(i,j)**2
      zx   (i,j) = ABS(zx(i,j))
   
      ! Limit the roughness length for momentum to a value smaller than
      ! half of the lower layer thickness

!FEA
!     zgz0m(i,j) = MIN( gz0(i,j), zdfip(i,j) - g )
      zgz0m(i,j) = MIN( g*z0(i,j), zdfip(i,j) - g )
!FEE

      ! Derive roughness length for heat over open sea using the friction
      ! velocity ustar derived from z0 for momentum (Charnock formula)
      ! Over land an sea ice, z0h = z0m.

!!$      IF (  fr_land(i,j) < 0.5 .AND.  t_g(i,j,nx) >=  t0_melt-1.7 ) THEN
!!$        zustar = SQRT ( zgz0m(i,j)/zalpha0 )
!!$        zgz0h(i,j) = g*zviscos*zalphah / MAX( 1.0E-8_wp, zustar )
!!$      ELSE
        zgz0h(i,j) = zgz0m(i,j)
!!$      ENDIF

      ! Limit the roughness length for momentum to a fixed value

      zgz0h(i,j)  = MIN ( zgz0h(i,j), zgz0hh)

      zxi         = zdfip(i,j)/ zgz0m(i,j)
      zxih        = zdfip(i,j)/ zgz0h(i,j)
      zy          = (akt/LOG(zxi))**2
   
      ! Stable case for land and water 
      ! ------------------------------
      IF ( zris(i,j) >= 0.0 ) THEN
  
        ztcm(i,j) = zy*zvpb(i,j)*MAX ( ztmmin, 1.0/  &
                   (1.0 + 10.0*zris(i,j)/SQRT(1.0 + 5.0*zris(i,j)) ) ) 
        ztch(i,j) = akt**2/(LOG(zxi)*LOG(zxih))*zvpb(i,j)*    &
                     MAX ( zthmin, 1.0/(1.0 + 15.0*zris(i,j)*  &
                     SQRT( 1.0 + 5.0*zris(i,j) ) ) ) 

      ! Unstable case
      ELSE
   
        ! Land
!!$        IF ( fr_land(i,j) >= 0.5 ) THEN
          ztcm(i,j) = zy*zvpb(i,j)*(1.0 - 10.0*zris(i,j)/   &
                      (1.0 + 75.0*zy*(zxi**zed3-1.0)**1.5*SQRT(-zris(i,j)) ))
          ztch(i,j) = akt**2/(LOG(zxi)*LOG(zxih))*zvpb(i,j)*   &
                      (1.0-15.0*zris(i,j)/(1.0 + 75.0*SQRT(zy)*akt/LOG(zxih)* &
                      (zxih**zed3-1.0)**1.5*SQRT(-zris(i,j)) ))
      ENDIF
   
      ! Store TCM, TCH, GZ0 and the skin layer values of TKVH, TKVM and TKE
   
      tcm  (i,j) = ztcm(i,j)/zvpb(i,j)
      tch  (i,j) = ztch(i,j)/zvpb(i,j)

!FEA
      ENDIF ! ( IF (llandmask(i,j) )
!FEE

    ENDDO
  ENDDO         

!------------------------------------------------------------------------------
! End of module procedure parturs
!------------------------------------------------------------------------------

END SUBROUTINE parturs
!--------------------------------------------------------


SUBROUTINE calc_albedo(im,jm,nzx,also,alth)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine returns the solar (shortwave) and thermal (longwave)
!    albedo of the surface - bare soil, snow and vegetation included
!
! Method:
!
!
!------------------------------------------------------------------------------

! Subroutine arguments: 
! --------------------
  INTEGER                 , INTENT(IN)  :: &
          im                             , & ! 
          jm                             , & !
          nzx
  REAL    (KIND=wp)       , INTENT(OUT) :: &
          also                           , & ! total solar albedo
          alth                               ! total thermal albedo

! Local variables: 
! --------------------
  INTEGER                               :: ist

  REAL (KIND=wp)                        :: &
       zdzwb                             , & !
       zsnow                             , & ! partial cover of snow
       zvege                                 ! partial cover of vegetation

  INTEGER                  ::  zdalbs(10)    ! Surface albedo for diffenrent soil types (maximum 10)

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine calc_albedo             
!------------------------------------------------------------------------------

  ! Albedo of soil type (without vegetation) as function of soil water content
  ! (in mH2O) and depth of the upper soil layer       

!! IF (lmulti_layer) THEN
     zdzwb = 1.0_wp / (2*czmls(1))
!!ELSE
!!   IF ( nlgw == 2 ) THEN
!!      zdzwb = 1.0_wp / cdzw12
!!   ELSE
!!      zdzwb = 1.0_wp / cdzw13
!!   ENDIF
!! ENDIF
  DO ist = 1, 10
     zdalbs(ist) = csalbw(ist) * zdzwb
  ENDDO

  ist = 10
  IF ( llandmask(im,jm) .OR. (t_snow(im,jm,nzx) >= t0_melt - 1.7 )) THEN 
     ist = NINT(soiltyp(im,jm))
  ENDIF
! IF (lmulti_layer) THEN
     also = csalb(ist) - zdalbs(ist)*w_so(im,jm,1,nzx)
! ELSE
!    also = csalb(ist) - zdalbs(ist)*w_g1(im,jm,nzx)
! ENDIF

  ! Snow cover and vegetation
  zsnow= 0.0
  IF ( llandmask(im,jm) .AND. w_snow(im,jm,nzx) > 0.0 ) THEN
!!!  replaced by (14.01.02)
     zsnow = MIN(1.0_wp,w_snow(im,jm,nzx)/cf_snow)
  ENDIF
  zvege    = plcov(im,jm)
!FEA
!  also    =  zsnow * csalb_snow                                           &
!          + (1.-zsnow) * (zvege * csalb_p + (1.-zvege) * also)
  also    =  zsnow * csalb_snow                                           &
          + (1.-zsnow) * (zvege * vegalb(im,jm) + (1.-zvege) * also)
!FEE
  alth = ctalb

  also = 0.9 
END SUBROUTINE calc_albedo
!--------------------------------------------------------


SUBROUTINE gen_area_indices
!------------------------------------------------------------------------------
!
! Description:
!   This subroutine generates area indices for the Stand Alone version:
!   
!------------------------------------------------------------------------------

! Begin Subroutine gen_area_indices
  sai=MAX(1.0E-6_wp,lai)
  eai=(1.0_wp-plcov)*sai
  tai=MAX(1.0E-6_wp,lai)
  tai=plcov*tai
  !              tai(i,j)=MAX(1.0E-6_wp,lai(i,j))
  !              tai(i,j)=plcov(i,j)*tai(i,j)  ! transpiration area index
  !              eai(i,j)=c_soil               ! evaporation area index
  !              sai(i,j)=c_lnd+tai(i,j)       ! surface area index
END SUBROUTINE gen_area_indices
!--------------------------------------------------------


SUBROUTINE near_surface (nx,lt2m,lu10m)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine is adapted from COSMO "near_surface" for TSE requirements.
!   This routine calculates values near the surface.
!   - winds in 10 m
!   - temperature, dew-point temperature and specific water vapor content in 2m
!   - accumulation of precipitation rates
!   - minimal and maximal temperature in 2m
!   - maximal expected squall
!   - mean values over forecast for solar and thermal heating and radiation
!   - in case of llm, some variables are set to the lowest level or to some
!     unrealistic values
!
! Method:
!   See Comments in the Sections.
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8236 1493
!  email:  uschaettler@dwd.d400.de
!
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.29       1999/05/11 Ulrich Schaettler
!  Initial release
! 1.30       1999/06/24 Matthias Raschendofer
!  Alternative call of the diagnosis of 2m- and 10m-values: call synop_diag.
! 1.32       1999/08/24 Guenther Doms
!  Reset of top_con and bas_con at every hour; 
!  direct inclusion of former routines 'low_winds_temp' and 'synop_diag'.
! 1.33       1999/10/14 Matthias Raschendorfer
!  The former routine 'synop_diag' is now included in routine 'turbtran' in
!  a more sophisticated way. Here, the previous code is left only for the
!  case of not running the turbulence routines.
! 2.8        2001/07/06 Ulrich Schaettler
!  Added summarization of surface fluxes
! 2.18       2002/07/16 Ulrich Schaettler
!  Changed treatment of dew-point temperature (by Lucio Torrisi, UGM Rome)
! 2.19       2002/10/24 Ulrich Schaettler
!  Moved re-initialization of vmax_10m to lmorg (initialize_loop)
! 3.5        2003/09/02 Ulrich Schaettler
!  Compute surface pressure here (instead of in the relaxation);
!  Thus the communication for ps can be avoided
! 3.13       2004/12/03 Thorsten Reinhardt
!  Changes for the new graupel scheme
! 3.15       2005/03/03 Jan-Peter Schulz
!  Adapt calculation of vmax_10m to the GME 40km/40L formulation
! 3.18       2006/03/03 Matthias Raschendorfer, Ulrich Schaettler
!  Introduction of rh_2m; Additional fields for Climate-LM Version.
!  Near surace levels also possible between the lowest two model levels
!  in the case of 'itype_synd==1'.
! 3.21       2006/12/04 Burkhardt Rockel, Ulrich Schaettler
!  Renamed sunshhrs, sodwdir to dursun, sodwddm
!  Do not calculate summations and meanvalues for ntstep==0
! V3_23        2007/03/30 Matthias Raschendorfer
!  Moved 'akt' to MOULE data_turbulence.
!  Accumulation of fields for topographic radiation correction (Matteo Buzzi)
! V3_24        2007/04/26 Ulrich Schaettler
!  Eliminated nincmxt and introduced control as for other increments
! V4_1         2007/12/04 Jan-Peter Schulz
!  Re-tuning of gusts by using wind speed at 10 m instead of 30 m
! V4_4         2008/07/16 Jan-Peter Schulz
!  Accumulation of fields for sub-grid scale orography scheme
! V4_8         2009/02/16 Ulrich Schaettler
!  Compute additional averaged values for radiation
!  Introduced several options for wind gusts
! V4_9         2009/07/16 Ulrich Schaettler, Christian Bollmann
!  Inserted a NEC compiler option directive
! V4_11        2009/11/30 Lucio Torrisi
!  Computation of averages for additional fields
! V4_12        2010/05/11 Oliver Fuhrer
!  Additional computations for sunshine duration
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_17        2011/02/24 Ulrich Blahak
!  Fixed index bounds (istartpar/jstartpar instead of istart/jstart) in computation
!  of zum1, zum2, ...
! V4_18        2011/05/26 Ulrich Schaettler
!  After adjusting all loop boundaries, the additional loops for setting some
!  fields at the boundary zone can be eliminated
!  Bug fix for itype_diag_gusts=2: use values in 10m instead of 30m (Oli Fuhrer)
! V4_21        2011/12/06 Jan-Peter Schulz
!  Introduce gust option itype_diag_gusts = 4. Here the gust factor weakly
!  depends on the mean wind speed at 10 m.
! V4_23        2012/05/10 Burkhardt Rockel (CLM)
!  Introduction of new diagnostic variable for maximum wind speed in 10m height
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer, Ulrich Blahak
!  Replaced qx-variables by using them from the tracer module
!  For the 2-moment microphysics scheme, added hail_gsp and prh_gsp (UB).
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  Bug Fix: summation of fluxes in multi-layer soil must only be done,
!    if lmulti_layer is TRUE
!  Introduced MESSy interface
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
! V_TSA        2015-12-01 Yiftach Ziv (XYZ)
!  updated near_surface subroutine according to COSMO "near_surface" V5_1

! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with INTENT(IN):
  INTEGER, INTENT (IN)       :: nx
  LOGICAL, INTENT (IN)       :: lt2m,lu10m 
!------------------------------------------------------------------------------
!
! Local scalars:
  REAL (KIND=wp)             :: &
    zvmin

  REAL (KIND=wp)             :: &
    z10g, z2g, zpke, zcm, zch , &
    zh, zhs, zsqcm, zchdcm    , &
    zh1, zlnz1, zdz, zh05m    , &
    zh2m
!   zh2

  INTEGER                    :: &
    i, j                      , &
    iend                      , & ! local end index for ICON block structure
    im1, jm1                      ! ii, jj in COSMO "near_surface". notation here adheres better to src_soil_multilay

! Local arrays:
!DL Take care for ICON block structure
  REAL (KIND=wp), ALLOCATABLE:: &
    zum1 ( :, :)              , & ! u-wind on the lowest level at the mass position
    zvm1 ( :, :)                  ! v-wind on the lowest level at the mass position
!!DLzum1 (ie,je)              , & ! u-wind on the lowest level at the mass position
!!DLzvm1 (ie,je)                  ! v-wind on the lowest level at the mass position

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE near_surface
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!  Section 1: Initializations
!------------------------------------------------------------------------------

  ! Allocate zum1 and zvm1
  IF (ymodel == 'COSMO') THEN
    ALLOCATE (zum1 (ie,je))
    ALLOCATE (zvm1 (ie,je))
  ELSE
    ALLOCATE (zum1 (nproma,nblock))
    ALLOCATE (zvm1 (nproma,nblock))
  ENDIF

  ! minimal wind velocity
  zvmin = 0.01_wp

  ! Velocity of the wind in the Prandtl layer (k=ke) and in the layer
  ! above (k=ke-1) for time step nnow.
  ! Interpolate wind speed in 30 m above ground (for gust computation).
  ! And interpolation of u and v
  iend = iendpar
  DO j = jstartpar,jendpar
    IF ((ymodel /= 'COSMO') .AND. (j==nblock)) iend = nlastproma
    DO i = istartpar,iend
       im1=MAX(i-1,1)
       jm1=MAX(j-1,1)

       zum1(i,j) = 0.5_wp * ( u(i,j,ke,nx) + u(im1,j,ke,nx) )
       zvm1(i,j) = 0.5_wp * ( v(i,j,ke,nx) + v(i,jm1,ke,nx) )
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
!  Section 3a: Compute wind, temperature and humidity at screen levels 
!              (10m and 2m) using a default scheme (itype_synopd = 1) based
!              on similarity theory
!------------------------------------------------------------------------------

  z10g = 10.0_wp * g
  z2g  =  2.0_wp * g

    iend = iendpar
    DO j = jstartpar, jendpar
      IF ((ymodel /= 'COSMO') .AND. (j==nblock)) iend = nlastproma
      DO i = istartpar, iendpar
!FEA>
      IF (.NOT.llandmask(i,j)) THEN
        u_10m(i,j) = rundef
        v_10m(i,j) = rundef
        t_2m(i,j)  = rundef
        CYCLE
      END IF
!FEE<
      ! some initializations
      ! check tcm and tch
        zcm       = MAX ( tcm(i,j) , 5.0E-4_wp )
        zch       = MAX ( tch(i,j) , 4.0E-5_wp )
      ! other variables
      zpke   = p0(i,j,ke) + pp(i,j,ke,nx)
      zh     = cp_d* t(i,j,ke,nx) + g* dz ! zhll(i,j)
      zhs    = cp_d* t_g(i,j,nx)  !  + g*hsurf(i,j)
      zsqcm  = SQRT (zcm)
      zchdcm = zch / ( akt * zsqcm )
      zh1    = g* dz ! ( zhll(i,j) - hsurf(i,j) )

      IF (zhs <= zh) THEN
        ! stable case
        zlnz1      = LOG ( (z10g + g*z0(i,j)) / g*z0(i,j) )  &        ! g*z0=gz0 in COSMO
                   - z10g / zh1 * LOG ( (zh1+g*z0(i,j)) / g*z0(i,j) )
        zlnz1      =(z10g/zh1 + zsqcm/akt * zlnz1)
        zh05m      = zhs + 0.25_wp * (zh - zhs)
        zh2m       = zh05m + 1.5_wp * g * (zh - zh05m) / (zh1 - 0.5_wp*g)
      ELSE
        ! unstable case
        zsqcm     = MAX ( zsqcm  , 0.01_wp )
        zchdcm    = MAX ( zchdcm , 0.01_wp )

        zdz       = MAX(zh1-z10g, 0.0_wp)
        zlnz1     = 1.0_wp - zsqcm/akt * LOG ( 1.0 + (EXP ( akt/zsqcm ) - 1.0_wp) &
                    *g*z0(i,j) * zdz / (zh1*(z10g + g*z0(i,j))))
        zdz       = MAX(zh1 - z2g, 0.0_wp)
        zh2m      = zhs + (zh - zhs) * (1.0 - zchdcm                              &
                    * LOG ( 1.0_wp + (EXP ( 1.0_wp/zchdcm) - 1.0_wp) * (g*z0(i,j) &
                        * zdz / (zh1 * (z2g + g*z0(i,j))) ) ))
      ENDIF

      ! wind in 10 m
      IF (lu10m) THEN
         u_10m(i,j) = zum1(i,j) * zlnz1
         v_10m(i,j) = zvm1(i,j) * zlnz1
      ENDIF
      ! temperature in 2 m
      IF (lt2m) THEN
         t_2m(i,j)  = ( zh2m - z2g )*cpdr ! - g*hsurf(i,j)) * cpdr
      ENDIF

    ENDDO
  ENDDO

!------------------------------------------------------------------------------
!- End Subroutine near_surface
!------------------------------------------------------------------------------
END SUBROUTINE near_surface

!>XYZ following JT of GUF
SUBROUTINE parturs_new (zh)
!------------------------------------------------------------------------------
!
! Description:
!   This module procedure calculates the transfer coefficients for momentum
!   (TCM) and heat and moisture (TCH) at the surface as well as the roughness
!   length over sea (GZ0).
!   This is a revised version that requires height of measurements- zh (usually dz_u=10m here).
!   temperature, humidity and wind are assumed to be on the same height.
!   This version is supposed to prevent computational instabilities over
!   high roughness terrain.
!
! Method:
!   Dyer-Businger equations as modified by Louis.
!
!------------------------------------------------------------------------------

! Subroutine arguments: None
! --------------------
!
! Local parameters:
! --------------------
  REAL (KIND=wp), PARAMETER :: &
    zed3    = 1.0_wp/3.0_wp,   & ! = 1/3 (needed as exponent) [--]
    z0hmax  = 0.1_wp,          & ! upper limit for roughness length of heat/moisture [m]
    zvmin   = 0.5_wp,          & ! Minimum wind velocity for transfer scheme [m/s]
    tcm_min = 5.0E-4_wp,       & ! Minimum value for tcm [--]
    tcm_max = 1.0_wp,          & ! Maximum value for tcm [--]
    tch_min = 5.0E-5_wp,       & ! Minimum value for tch [--]
    tch_max = 0.1_wp               ! Maximum value for tch [--]
  
! Local scalars:
! --------------------
  INTEGER                  ::  &
    i ,                        & ! loop index in x-direction              
    j ,                        & ! loop index in y-direction              
    nx                           ! time-level for computation
 
  REAL (KIND=wp), INTENT(IN) :: zh ! height to where to compute the transfer [m]

  REAL (KIND=wp)             :: &
    tv_g  ,                     & ! virtual temperature of ground [K]
    tv_ke ,                     & ! virtual temperature of air [K]
    zwind ,                     & ! horizontal wind speed [m/s]
    fm    ,                     & ! stability function value to modify neutral momentum transfer coefficient [--]
    fh    ,                     & ! stability function value to modify neutral heat transfer coefficient [--]
    tcm_n ,                     & ! momentum transfer coefficient in neutral case [--]
    tch_n ,                     & ! heat transfer coefficient in neutral case [--]
    z0m   ,                     & ! roughness length for momentum [m]
    z0h   ,                     & ! roughness length for heat/moisture [m]
    brn   ,                     & ! bulk richardson number [--]
    abrn  ,                     & ! absolute bulk richardson number [--]
    z0mmax                        ! upper limit for roughness length of momentum [m]
  

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine parturs_new
!------------------------------------------------------------------------------

! Fixed settings
 nx = nnow   !current time level
 z0mmax = zh/2.0_wp-1.0_wp   !limit z0 for momentum to less than 1/2 of atmospheric level height
                                    !z0hmin is already specified as parameter in header
                       
! Start the loop over all landpoints
 DO i = 1 , ie
   DO j = 1 , je 
     IF (llandmask(i,j)) THEN

     
     ! (1) PREPARATIONS
     ! (this section could also be externalized as constant fields, as values remains same all the time at each gridpoint)
     ! Setting of roughness length for heat,
     ! Ensure limitations of roughness lengths for momentum and heat 
   
       z0m = min( z0(i,j), z0mmax )   !roughness length for momentum
       z0h = min( z0m    , z0hmax )   !heat roughness length is same as for momentum but further limited to z0hmax=0.1m

     ! Compute transfer coefficients for the neutral case
       tcm_n = ( akt / log(zh/z0m) )**2
       tch_n = akt**2 / ( log(zh/z0m) * log(zh/z0h) )


     ! (2) BULK RICHARDSON NUMBER
     ! Compute wind speed from components, minimum value of 0.01
       zwind = max( sqrt(u(i,j,ke,nx)**2 + v(i,j,ke,nx)**2) , zvmin)
  
     ! Transform temperatures to virtual temperatures:
     ! last term corrects to potential temperature in dz=2m height
       tv_g  = t_g(i,j   ,nx) * (1.0_wp + rvd_m_o*qv_s(i,j   ,nx)) + cpdr*g*dz
       tv_ke =   t(i,j,ke,nx) * (1.0_wp + rvd_m_o*  qv(i,j,ke,nx))

     ! Compute bulk Richardson number (BRN) 
       brn = ( g*(tv_ke - tv_g)*(zh-z0m) ) / ( t_g(i,j,nx) * zwind**2 )
     ! brn = ( tv_ke - tv_g + cpdr*g*zh )*g*zh/ ( t_g(i,j,nx) *zwind**2 )  #DWD version of subroutine parturs

  
     ! (3) STABILITY FUNCTIONS
     ! Evaluate stability functions fm and fh depending on sign of BRN
       IF ( brn >= 0 ) THEN   !STABLE CASE: reduction of transfer coefficients
         fm = 1.0_wp/( 1.0_wp + 10.0_wp*brn/sqrt(1.0_wp + 5.0_wp*brn) )
         fh = 1.0_wp/( 1.0_wp + 15.0_wp*brn/sqrt(1.0_wp + 5.0_wp*brn) )  
     
       ELSE  ! UNSTABLE CASE: increase of transfer coefficients
         abrn = abs(brn)  !take absolute value 
         fm = 1.0_wp + 10.0_wp*abrn / ( 1.0_wp + 75.0_wp*tcm_n*((zh/z0m)**zed3 -1.0_wp)**1.5_wp *sqrt(abrn) ) 
         fh = 1.0_wp + 15.0_wp*abrn / ( 1.0_wp + 75.0_wp*tch_n*((zh/z0h)**zed3 -1.0_wp)**1.5_wp *sqrt(abrn) )     
  
       ENDIF  ! sign of brn


     ! (4) FINAL TRANSFER COEFFICIENTS
     ! Neutral coefficients modified by stability function values
       tcm(i,j) = tcm_n * fm
       tch(i,j) = tch_n * fh
  
     ! Final checks for upper / lower bounds
       tcm(i,j) = max( min(tcm(i,j),tcm_max) , tcm_min)
       tch(i,j) = max( min(tch(i,j),tch_max) , tch_min)

     ENDIF
   ENDDO
 ENDDO !loop over land points

!------------------------------------------------------------------------------
! End of module procedure parturs_new
!------------------------------------------------------------------------------

END SUBROUTINE parturs_new
!--------------------------------------------------------
!<XYZ following JT of GUF

!#################################### NEW #################################################
!>XYZ following JT of GUF
SUBROUTINE parturs_newblock (zh, iblock, nvec, fr_land, z0, u, v, t, qv, qv_s, t_g, tcm, tch)
!------------------------------------------------------------------------------
!
! Description:
!   This module procedure calculates the transfer coefficients for momentum
!   (TCM) and heat and moisture (TCH) at the surface as well as the roughness
!   length over sea (GZ0).
!   This is a revised version that requires height of measurements- zh (usually dz_u=10m here).
!   temperature, humidity and wind are assumed to be on the same height.
!   This version is supposed to prevent computational instabilities over
!   high roughness terrain.
!
! Method:
!   Dyer-Businger equations as modified by Louis.
!
!------------------------------------------------------------------------------

! Subroutine arguments: 
! --------------------
  INTEGER, INTENT(IN)        :: iblock, nvec
  REAL (KIND=wp), INTENT(IN) :: zh ! height to where to compute the transfer [m]
  REAL (KIND=wp), DIMENSION(nvec), INTENT(IN) :: &
             fr_land          , & ! land fraction
             z0               , & ! vegetation roughness length                    ( m )
             u                , & ! zonal wind speed                              ( m/s )
             v                , & ! meridional wind speed                         ( m/s )
             t                , & ! temperature                                   (  k  )
             qv               , & ! specific water vapor content                  (kg/kg)
             qv_s             , & ! specific water vapor content (surface)        (kg/kg)
             t_g                  ! weighted surface temperature                  (  K  )
  REAL (KIND = wp), DIMENSION(nvec), INTENT(INOUT) :: &
                  tch         , & ! turbulent transfer coefficient for heat       ( -- )
                  tcm             ! turbulent transfer coefficient for momentum   ( -- )
             

!
! Local parameters:
! --------------------
  REAL (KIND=wp), PARAMETER :: &
    zed3    = 1.0_wp/3.0_wp,   & ! = 1/3 (needed as exponent) [--]
    z0hmax  = 0.1_wp,          & ! upper limit for roughness length of heat/moisture [m]
    zvmin   = 0.5_wp,          & ! Minimum wind velocity for transfer scheme [m/s]
    tcm_min = 5.0E-4_wp,       & ! Minimum value for tcm [--]
    tcm_max = 1.0_wp,          & ! Maximum value for tcm [--]
    tch_min = 5.0E-5_wp,       & ! Minimum value for tch [--]
    tch_max = 0.1_wp             ! Maximum value for tch [--]
  
! Local scalars:
! --------------------
  INTEGER                  ::  &
    i                            ! loop index
 
  REAL (KIND=wp)             :: &
    tv_g  ,                     & ! virtual temperature of ground [K]
    tv_ke ,                     & ! virtual temperature of air [K]
    zwind ,                     & ! horizontal wind speed [m/s]
    fm    ,                     & ! stability function value to modify neutral momentum transfer coefficient [--]
    fh    ,                     & ! stability function value to modify neutral heat transfer coefficient [--]
    tcm_n ,                     & ! momentum transfer coefficient in neutral case [--]
    tch_n ,                     & ! heat transfer coefficient in neutral case [--]
    z0m   ,                     & ! roughness length for momentum [m]
    z0h   ,                     & ! roughness length for heat/moisture [m]
    brn   ,                     & ! bulk richardson number [--]
    abrn  ,                     & ! absolute bulk richardson number [--]
    z0mmax                        ! upper limit for roughness length of momentum [m]
  
  LOGICAL, DIMENSION(nvec)   :: lsm ! land sea mask

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine parturs_newblock
!------------------------------------------------------------------------------

! Fixed settings
 z0mmax = zh/2.0_wp-1.0_wp   !limit z0 for momentum to less than 1/2 of atmospheric level height
                                    !z0hmin is already specified as parameter in header
                       
! Start the loop over all landpoints
   lsm = fr_land >= 0.5_wp
   DO i = 1 , nvec
     IF (lsm(i)) THEN

     
     ! (1) PREPARATIONS
     ! (this section could also be externalized as constant fields, as values remains same all the time at each gridpoint)
     ! Setting of roughness length for heat,
     ! Ensure limitations of roughness lengths for momentum and heat 
   
       z0m = min( z0(i), z0mmax )   !roughness length for momentum
       z0h = min( z0m    , z0hmax )   !heat roughness length is same as for momentum but further limited to z0hmax=0.1m

     ! Compute transfer coefficients for the neutral case
       tcm_n = ( akt / log(zh/z0m) )**2
       tch_n = akt**2 / ( log(zh/z0m) * log(zh/z0h) )


     ! (2) BULK RICHARDSON NUMBER
     ! Compute wind speed from components, minimum value of 0.01
       zwind = max( sqrt(u(i)**2 + v(i)**2) , zvmin)
  
     ! Transform temperatures to virtual temperatures:
     ! last term corrects to potential temperature in dz=2m height
       tv_g  = t_g(i) * (1.0_wp + rvd_m_o*qv_s(i)) + cpdr*g*dz
       tv_ke =   t(i) * (1.0_wp + rvd_m_o*  qv(i))

     ! Compute bulk Richardson number (BRN) 
       brn = ( g*(tv_ke - tv_g)*(zh-z0m) ) / ( t_g(i) * zwind**2 )
     ! brn = ( tv_ke - tv_g + cpdr*g*zh )*g*zh/ ( t_g(i) *zwind**2 )  #DWD version of subroutine parturs

  
     ! (3) STABILITY FUNCTIONS
     ! Evaluate stability functions fm and fh depending on sign of BRN
       IF ( brn >= 0 ) THEN   !STABLE CASE: reduction of transfer coefficients
         fm = 1.0_wp/( 1.0_wp + 10.0_wp*brn/sqrt(1.0_wp + 5.0_wp*brn) )
         fh = 1.0_wp/( 1.0_wp + 15.0_wp*brn/sqrt(1.0_wp + 5.0_wp*brn) )  
     
       ELSE  ! UNSTABLE CASE: increase of transfer coefficients
         abrn = abs(brn)  !take absolute value 
         fm = 1.0_wp + 10.0_wp*abrn / ( 1.0_wp + 75.0_wp*tcm_n*((zh/z0m)**zed3 -1.0_wp)**1.5_wp *sqrt(abrn) ) 
         fh = 1.0_wp + 15.0_wp*abrn / ( 1.0_wp + 75.0_wp*tch_n*((zh/z0h)**zed3 -1.0_wp)**1.5_wp *sqrt(abrn) )     
  
       ENDIF  ! sign of brn


     ! (4) FINAL TRANSFER COEFFICIENTS
     ! Neutral coefficients modified by stability function values
       tcm(i) = tcm_n * fm
       tch(i) = tch_n * fh
  
     ! Final checks for upper / lower bounds
       tcm(i) = max( min(tcm(i),tcm_max) , tcm_min)
       tch(i) = max( min(tch(i),tch_max) , tch_min)

     ENDIF
 ENDDO !loop over land points

!------------------------------------------------------------------------------
! End of module procedure parturs_newblock
!------------------------------------------------------------------------------

END SUBROUTINE parturs_newblock
!--------------------------------------------------------
!<XYZ following JT of GUF

!------------------------------------------------------------------------------
! Utility routine for checking fields
!------------------------------------------------------------------------------
SUBROUTINE stats(msg,field,mask,UNIT)
!----------------------------------------------------------------------
! Description:
!
! Writes Status of variables 
!----------------------------------------------------------------------

  IMPLICIT NONE

  ! Subroutine Arguments
  CHARACTER (LEN=*),        INTENT(IN)           :: msg
  REAL (KIND=wp),           INTENT(IN)           :: field(:,:)
  LOGICAL,                  INTENT(IN), OPTIONAL :: mask(:,:)
  INTEGER,                  INTENT(IN), OPTIONAL :: UNIT

  ! parameters
  REAL (KIND=wp), PARAMETER :: undefgrib=-1.0E7_wp

  ! local variables
  INTEGER                   :: &
     lb(2)                   , & ! 
     ub(2)                   , & ! 
     lbm(2)                  , & ! 
     ubm(2)                      ! 
  LOGICAL, ALLOCATABLE      :: &
     m(:,:)                  , & ! 
     mm(:,:)                     ! 
  ! local scalars:
  INTEGER                   :: &
     locmin(2)                  , & ! location of min in the 2D field
     locmax(2)                      ! location of max in the 2D field
  REAL (KIND=wp)            :: &
     rmin                    , & ! 
     rmax                    , & ! 
     rmean                       ! 
  INTEGER                   :: &
     lun                     , & ! 
     nn                          ! 

  locmin = (1,1)
  locmax = (1,1)

  ! get bounds
  lb=LBOUND(field)
  ub=UBOUND(field)

!!!
!!!  WRITE(6,*) 'size of ub: ',size(ub)
!!!  WRITE(6,*) 'size of field: ',size(field)
!!!

  ! set UNIT
  IF (PRESENT(UNIT)) THEN
     lun=UNIT
  ELSE
     lun=6
  ENDIF

  ! make mask
  ALLOCATE(m(lb(1):ub(1),lb(2):ub(2)))
  IF (PRESENT(mask)) THEN
     lbm=LBOUND(mask)
     ubm=UBOUND(mask)
     IF (lb(1) /= lbm(1) .OR. lb(2) /= lbm(2)     &
        .OR.ub(1) /= ubm(1) .OR. ub(2) /= ubm(2)) THEN
!!!        WRITE(5,*) 'ERROR: In subroutine stats, bound of mask do not match field'
        WRITE(6,*) 'ERROR: In subroutine stats, bound of mask do not match field'
!!!        WRITE(6,*) 'lb(1): ',lb(1)
!!!        WRITE(6,*) 'lbm(1): ',lbm(1)
!!!        WRITE(6,*) 'lb(2): ',lb(2)
!!!        WRITE(6,*) 'lbm(2): ',lbm(2)
!!!        WRITE(6,*) 'ub(1): ',ub(1)
!!!        WRITE(6,*) 'ubm(1): ',ubm(1)
!!!        WRITE(6,*) 'ub(2): ',ub(2)
!!!        WRITE(6,*) 'ubm(2): ',ubm(2)
!!!
        RETURN
     ENDIF
     m(:,:)=mask(:,:)
  ELSE
     m(:,:)=.TRUE.
  ENDIF

  ! make undef and NaN mask
  ALLOCATE(mm(lb(1):ub(1),lb(2):ub(2)))
  mm=(m .AND. (field/=rundef) .AND. (field/=undefgrib) .AND. (field==field))

  ! compute stats
  nn=COUNT(mm)
  IF (lb(1)==ub(1) .AND. lb(2)==ub(2)) THEN
    nn=1
!   mm=rundef
!US error with gnu: Can't convert REAL(8) to LOGICAL(4)
!   do you mean:
    mm = (field/=rundef)
  ENDIF
  IF (nn>0) THEN
     rmin=MINVAL(field,mm)
     locmin=MINLOC(field,mm)
     rmax=MAXVAL(field,mm)
     locmax=MAXLOC(field,mm)
     rmean=SUM(field,mm)/REAL(nn,wp)
  ELSE
     rmin=rundef
     locmin(1)=-1
     locmin(2)=-1
     rmax=rundef
     locmax(1)=-1
     locmax(2)=-1
     rmean=rundef
  ENDIF

  ! write message
  WRITE(lun,'(a15,1x,g14.5,2(g14.5,2i4),5i7)')                                      &
       TRIM(msg),rmean,rmin,locmin(1),locmin(2),rmax,locmax(1),locmax(2),           &
       COUNT(m), nn, COUNT(m .AND. field==rundef), COUNT(m .AND. field==undefgrib), &
       COUNT(m .AND. field/=field)

  ! cleanup
  DEALLOCATE(m,mm)

END SUBROUTINE stats

!==============================================================================

SUBROUTINE model_abort (my_id,ierrorcode,yerrorstring,yroutine, implerrorcode)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine stops the program in case of errors. Model_abort prints 
!   a self-defined error message, the errorcode and the calling routine.
!   If the error occurs in the routine of the message passing library used, 
!   the optional parameter implerrorcode must be present. In this case
!   also an appropriate error message is printed. This is determined via 
!   MPI_ERROR_STRING and implerrorcode.
!   If the error occurs in several nodes, this message will also be printed
!   several times.
!
! Method:
!   Printing an error message and then MPI_ABORT.
!
!------------------------------------------------------------------------------

! Parameter list:

INTEGER,           INTENT(IN) ::     &
  my_id,        & ! id of this processor
  ierrorcode      ! self-defined integer code of the error detected

CHARACTER (LEN=*), INTENT(IN) ::     &
   yerrorstring   ! self-defined error message
CHARACTER (LEN=*), INTENT(IN) ::     &
   yroutine       ! calling routine

INTEGER, OPTIONAL, INTENT(IN) ::     &
  implerrorcode   ! error-code of the message passing library

! Local variables:
INTEGER                    :: i, nzerrcode, nzlen
LOGICAL                    :: lzopen
CHARACTER (LEN=100) ymplmsg   ! for MPI error message

!------------------------------------------------------------------------------

! here we could even implement a parallel abort
  IF (my_id == 0) THEN

    ! print the error message
!   IF (PRESENT (implerrorcode)) THEN
!     ! this is parallel mode
!     CALL MPI_ERROR_STRING (implerrorcode, ymplmsg, nzlen, nzerrcode)

!     PRINT *,'*------------------------------------------------------------*'
!     PRINT *,'*    PROGRAM TERMINATED BECAUSE OF MPI ERRORS DETECTED'
!     PRINT *,'*              IN ROUTINE:   ',yroutine
!     PRINT *,'*'
!     PRINT *,'*    ERROR CODE is ',ierrorcode,': '
!     PRINT *,'*    MPI ROUTINE:  ',yerrorstring(1:LEN_TRIM(yerrorstring))
!     PRINT *,'*    MPI ERROR CODE is ', implerrorcode,':  ',ymplmsg
!     PRINT *,'*------------------------------------------------------------*'
!   ELSE
      PRINT *,'*------------------------------------------------------------*'
      PRINT *,'*    PROGRAM TERMINATED BECAUSE OF ERRORS DETECTED'
      PRINT *,'*              IN ROUTINE:   ',yroutine
      PRINT *,'*'
      PRINT *,'*    ERROR CODE is ',ierrorcode
      PRINT *,'*    ', yerrorstring(1:LEN_TRIM(yerrorstring))
      PRINT *,'*------------------------------------------------------------*'
!   ENDIF

!   ! Check, whether there are open files and close them
!   DO i = iunit_start, iunit_end
!     IF (iunit_table(i) == 1) THEN
!       ! inquire, whether unit i is open
!       lzopen = .FALSE.
!       INQUIRE (UNIT=i, OPENED=lzopen)
!       IF (lzopen) THEN
!         CLOSE (i)
!       ENDIF
!     ENDIF
!   ENDDO

! ELSE
  ENDIF

  CALL EXIT (ierrorcode)

END SUBROUTINE model_abort


!==============================================================================
END MODULE tsa_lmparam
