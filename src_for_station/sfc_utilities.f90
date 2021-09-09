!! Utility routines related to the surface schemes
!!  contains (up to now)
!!   - diag_snowfrac_tg
!!       diagnosis of snow fraction and computation of t_g
!!   - diag_lh_surfflux
!!       calculation of latent heat surface flux (lhfl_s) from 
!!       water vapour surface flux (qvfl_s)
!!
!! Note : All OpenACC kernels in this module are asynchronous and use the
!!        default async queue. When modifying or adding kernels, we have to
!!        make sure they are asynchronous and that the correct queue is used.
!
! Current Code Owner: DWD, Juergen Helmert
!  phone:  +49  69  8062 2704
!  fax:    +49  69  8062 3721
!  email:  juergen.helmert@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! V5_4e        2017-03-23 Ulrich Schaettler
!  Initial release for COSMO
! V5_4h        2017-12-15 Xavier Lapillonne
!  Port subroutine diag_lh_surfflux to GPU
! V5_5         2018-02-23 Ulrich Schaettler
!  Modifications to run the full block of parameterizations on GPU
!  OpenACC port of SR diag_snowfrac_tg
! V5_7         2020-02-21 Stefan Ruedisuehli, Yannick Boetzel
!  Use kind parameter vpp for variable precision physics
!  Made OpenACC kernels asynchronous (YB)
!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!==============================================================================

MODULE  sfc_utilities

!==============================================================================
!
! Description:
! 
!   Routines (module procedures) currently contained:

#ifdef __COSMO__
USE kind_parameters, ONLY : &
  wp       , & ! KIND-type parameters for real variables (working precision)
  vpp          ! KIND-type parameters for real variables (variable precision physics)

USE data_constants,  ONLY : dbl_eps => rprecision,                            &
                            t0_melt,    &   ! absolute zero for temperature
                            tf_salt,    &   ! salt water freezing point
                            lh_v, lh_s, lh_f
USE data_runcontrol, ONLY : lmulti_snow
USE sfc_terra_data,  ONLY : cf_snow_vpp          => cf_snow, &
                            tune_minsnowfrac_vpp => tune_minsnowfrac, &
                            idiag_snowfrac
#endif

#ifdef __ICON__
USE mo_kind,                ONLY: wp, vpp   ! KIND-type parameter for real variables
USE mo_physical_constants,  ONLY: t0_melt=>tmelt, &
                                  lh_v=>alv     , &
                                  lh_s=>als     , &
                                  lh_f=>alf
USE mo_phyparam_soil,       ONLY: cf_snow      ! soil and vegetation parameters for TILES
USE mo_math_constants,      ONLY: dbl_eps
USE mo_lnd_nwp_config,      ONLY: lmulti_snow, idiag_snowfrac
USE mo_nwp_tuning_config,   ONLY: tune_minsnowfrac
#endif

! do not import these modules to TSA
!USE sfc_flake_data, ONLY : h_Ice_min_flk_vpp => h_Ice_min_flk
!USE sfc_seaice,     ONLY : hice_min

#ifdef ALLOC_WKARR
USE sfc_terra_data, ONLY : lzurban
#endif

!==============================================================================

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: diag_snowfrac_tg, diag_lh_surfflux

!==============================================================================

! define the variables from sfc_flake_data, sfc_seaice here for TSA:
REAL (KIND=vpp) ::    &
  h_Ice_min_flk_vpp  = 1.0E-9_vpp            ! Minimum ice thickness [m]

REAL (KIND=wp)  ::    &
  hice_min           = 0.05_wp               ! minimum sea-ice thickness [m]

!==============================================================================

CONTAINS

!==============================================================================

!-------------------------------------------------------------------------

SUBROUTINE diag_snowfrac_tg (ivstart, ivend, lc_class, i_lc_urban, t_snow,      &
                             t_soiltop, w_snow, rho_snow, freshsnow, sso_sigma, &
                             tai, snowfrac, t_g, meltrate)

!-------------------------------------------------------------------------

  INTEGER, INTENT (IN) :: ivstart, ivend ! start and end-indices of the computation

  REAL(vpp), DIMENSION(:), INTENT(IN) :: t_snow, t_soiltop, w_snow, rho_snow, &
                                        freshsnow, sso_sigma, tai

  INTEGER,  INTENT(IN), OPTIONAL :: lc_class(:)    ! list of land-cover classes
  INTEGER,  INTENT(IN), OPTIONAL :: i_lc_urban     ! land-cover class index for urban / artificial surface

  REAL(vpp), INTENT(IN), OPTIONAL :: meltrate(:)    ! snow melting rate in kg/(m**2*s)

  REAL(vpp), INTENT(OUT) :: snowfrac(:), t_g(:)

!-------------------------------------------------------------------------

  INTEGER  :: ic
  REAL(vpp) :: h_snow, snowdepth_fac, sso_fac, lc_fac, lc_limit, tai_mod(ivend)

#ifndef ALLOC_WKARR
  LOGICAL  :: lzurban(ivend)
#endif

!-------------------------------------------------------------------------

  ! Subroutine parameters IN
  !$acc enter data async create(tai_mod)

  IF (PRESENT(lc_class) .AND. PRESENT(i_lc_urban)) THEN
    !$acc parallel async default(present)
    !$acc loop gang vector
    DO ic = ivstart, ivend
      lzurban(ic) = (lc_class(ic) == i_lc_urban)
    ENDDO
    !$acc end parallel
  ELSE
    !$acc parallel async default(present)
    !$acc loop gang vector
    DO ic = ivstart, ivend
      lzurban(ic) = .FALSE.
    ENDDO
    !$acc end parallel
  ENDIF

  ! Modified tai that is increased in urban areas
  !$acc parallel async default(present)
  !$acc loop gang vector
  DO ic = ivstart, ivend
    IF (lzurban(ic)) THEN
      tai_mod(ic) = MAX(3.0_vpp,tai(ic))
    ELSE
      tai_mod(ic) = tai(ic)
    ENDIF
  ENDDO
  !$acc end parallel

  SELECT CASE (idiag_snowfrac)
  CASE (1) ! old parameterization depending on SWE only
    !$acc parallel async default(present)
    !$acc loop gang vector
    DO ic = ivstart, ivend
      snowfrac(ic) = MIN(1.0_vpp, w_snow(ic)/cf_snow_vpp)
      t_g(ic) = t_snow(ic) + (1.0_vpp - snowfrac(ic))*(t_soiltop(ic) - t_snow(ic))
    ENDDO
    !$acc end parallel
  CASE (2, 20) ! more advanced parameterization depending on snow depth, accounts also for vegetation and SSO
    !$acc parallel async default(present)
    !$acc loop gang vector
    DO ic = ivstart, ivend
      IF (w_snow(ic) <= 1.e-6_vpp) THEN
        snowfrac(ic) = 0.0_vpp
      ELSE
        h_snow = 1000.0_vpp*w_snow(ic)/rho_snow(ic)  ! snow depth in m
        sso_fac = SQRT(0.025_vpp*MAX(25.0_vpp,sso_sigma(ic)*(1.0_vpp-freshsnow(ic))))
        snowdepth_fac = h_snow*(17.5_vpp*freshsnow(ic)+5.0_vpp+5.0_vpp/sso_fac*(1.0_vpp-freshsnow(ic)))
        lc_fac   = MAX(1.0_vpp,SQRT(2.5_vpp*tai_mod(ic)))
        IF (lzurban(ic)) THEN
          lc_limit = 0.875_vpp ! this accounts for the effect of human activities on snow cover
        ELSE
          lc_limit = 1.0_vpp
        ENDIF
        snowfrac(ic) = MIN(lc_limit,snowdepth_fac/lc_fac)
      ENDIF
      t_g(ic) = t_snow(ic) + (1.0_vpp - snowfrac(ic))*(t_soiltop(ic) - t_snow(ic))
    ENDDO
    !$acc end parallel
  CASE (3, 30)  ! similar to option 2, but somewhat less snow cover and limit over high vegetation
    !$acc parallel async default(present)
    !$acc loop gang vector
    DO ic = ivstart, ivend
      IF (w_snow(ic) <= 1.e-6_vpp) THEN
        snowfrac(ic) = 0.0_vpp
      ELSE
        h_snow = 1000.0_vpp*w_snow(ic)/rho_snow(ic)  ! snow depth in m
        sso_fac = SQRT(0.025_vpp*MAX(25.0_vpp,sso_sigma(ic)*(1.0_vpp-freshsnow(ic))))
        snowdepth_fac = h_snow*(17.5_vpp*freshsnow(ic)+5.0_vpp+5.0_vpp/sso_fac*(1.0_vpp-freshsnow(ic)))
        lc_fac   = MAX(1.0_vpp,SQRT(5.0_vpp*tai_mod(ic)))
        IF (lzurban(ic)) THEN
          lc_limit = 0.8_vpp ! this accounts for the effect of human activities on snow cover
        ELSE
          lc_limit = MAX(0.925_vpp,MIN(1.0_vpp,1.0_vpp/MAX(0.1_vpp,2.5_vpp*tai_mod(ic))**0.125_vpp))
        ENDIF
        snowfrac(ic) = MIN(lc_limit,snowdepth_fac/lc_fac)
      ENDIF
      t_g(ic) = t_snow(ic) + (1.0_vpp - snowfrac(ic))*(t_soiltop(ic) - t_snow(ic))
    ENDDO
    !$acc end parallel
  CASE (4, 40)  ! same as option 3, but even more restrictive snow cover limit over high vegetation
    !$acc parallel async default(present)
    !$acc loop gang vector
    DO ic = ivstart, ivend
      IF (w_snow(ic) <= 1.e-6_vpp) THEN
        snowfrac(ic) = 0.0_vpp
      ELSE
        h_snow = 1000.0_vpp*w_snow(ic)/rho_snow(ic)  ! snow depth in m
        sso_fac = SQRT(0.025_vpp*MAX(25.0_vpp,sso_sigma(ic)*(1.0_vpp-freshsnow(ic))))
        snowdepth_fac = h_snow*(17.5_vpp*freshsnow(ic)+5.0_vpp+5.0_vpp/sso_fac*(1.0_vpp-freshsnow(ic)))
        lc_fac   = MAX(1.0_vpp,SQRT(5.0_vpp*tai_mod(ic)))
        IF (lzurban(ic)) THEN
          lc_limit = 0.8_vpp ! this accounts for the effect of human activities on snow cover
        ELSE
          lc_limit = MAX(0.85_vpp,MIN(1.0_vpp,1.0_vpp/MAX(0.1_vpp,2.5_vpp*tai_mod(ic))**0.125_vpp))
        ENDIF
        snowfrac(ic) = MIN(lc_limit,snowdepth_fac/lc_fac)
      ENDIF
      t_g(ic) = t_snow(ic) + (1.0_vpp - snowfrac(ic))*(t_soiltop(ic) - t_snow(ic))
    ENDDO
    !$acc end parallel
  END SELECT

  ! For the single-layer scheme, t_soiltop represents the weighted average between the snow-covered
  ! part and the snow-free part in the case of melting snow concurrent with above-freezing soil temperatures
  IF (.NOT. lmulti_snow) THEN
    !$acc parallel async default(present)
    !$acc loop gang vector
    DO ic = ivstart, ivend
      IF (t_soiltop(ic) > t0_melt .AND. t_snow(ic) >= t0_melt - dbl_eps) t_g(ic) = t_soiltop(ic)
    ENDDO
    !$acc end parallel
  ENDIF

  SELECT CASE (idiag_snowfrac)
  CASE (20, 30, 40)
    ! Artificially reduce snow-cover fraction in case of melting snow in order to reduce the ubiquitous
    ! cold bias in such situations
    IF (PRESENT(meltrate)) THEN
      !$acc parallel async default(present)
      !$acc loop gang vector
      DO ic = ivstart, ivend
        snowfrac(ic) = MIN(snowfrac(ic),MAX(tune_minsnowfrac_vpp,1.0_vpp-9000.0_vpp*tai_mod(ic)*meltrate(ic)))
      ENDDO
      !$acc end parallel
    ENDIF
  END SELECT

  !$acc exit data async delete(tai_mod)

END SUBROUTINE diag_snowfrac_tg

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE diag_lh_surfflux (ivstart, ivend, llake, lseaice,                &
                             qvfl_s, fr_land, depth_lk, t_g, h_ice, w_snow, &
                             lhfl_s)

!------------------------------------------------------------------------------
!
! Description:
!  This subroutine calculates the latent heat flux at the surface from
!  the water vapour surface flux. Different values for latent heat of 
!     - vapourization:    lh_v    (liquid to gas)
!     - sublimation:      lh_s    (solid  to gas)
!     - fusion:           lh_f    (solid  to liquid)
!  are taken into account. 
!
!  Depending on whether FLake or sea-ice scheme are running, there are 
!  different methods to check, whether ice is present on a water point or not.
!
!  There is a difficulty to identify lake grid points vs. sea grid points, 
!  if FLake is not running. Because then we do not have the depth_lk parameter
!  available. The following routine tries to set the correct factor for 
!  computing lhfl_s in all cases.
!
!------------------------------------------------------------------------------

  INTEGER, INTENT (IN) :: ivstart, ivend ! start and end-indices of the computation

  REAL(wp), DIMENSION(:), INTENT(IN)  :: qvfl_s, fr_land, depth_lk,        &
                                         t_g, h_ice, w_snow
  REAL(wp), DIMENSION(:), INTENT(OUT) :: lhfl_s

  LOGICAL,  INTENT(IN) :: llake, lseaice

!------------------------------------------------------------------------------

  INTEGER  :: iv
  REAL(wp) :: zlh

  !$acc parallel async default(present)
  !$acc loop gang vector
  DO iv = ivstart, ivend
    zlh = lh_v

    IF (fr_land(iv) < 0.5_wp) THEN   ! Wasserpunkte

      ! in this block, also lake points are treated. 
      ! then these points are set, even if Flake is not running
      IF (.NOT. lseaice) THEN
        ! now all water points (lake and sea) are treated
        IF (t_g(iv) < tf_salt) THEN
          zlh = lh_s
        ENDIF
      ELSE
        ! here, also lake points are treated
        IF (h_ice(iv) >= hice_min) THEN
          zlh = lh_s
        ENDIF
      ENDIF

      IF (llake) THEN
        ! now we have to treat only flake points
        IF (depth_lk(iv) > 0.0_wp) THEN
          ! reset zlh, which might have been set above
          zlh = lh_v
          IF (h_ice(iv) >= REAL(h_Ice_min_flk_vpp,wp)) THEN
            zlh = lh_s
          ENDIF
        ENDIF
      ENDIF
    ELSE          ! Landpunkte
      ! This setting is only in effect if TERRA is not running
      ! TERRA recomputes this flux (with due regard for partial snow cover)
      IF (w_snow(iv) > 0.0_wp) THEN
        zlh = lh_s
      ENDIF
    ENDIF
    lhfl_s(iv) = qvfl_s(iv) * zlh
  ENDDO
  !$acc end parallel

END SUBROUTINE diag_lh_surfflux

!==============================================================================

END MODULE sfc_utilities
