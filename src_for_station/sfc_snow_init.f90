!------------------------------------------------------------------------------
!+ Source module  "sfc_snow_init"
!------------------------------------------------------------------------------

MODULE sfc_snow_init

!------------------------------------------------------------------------------
!
! Description:
!   This module contains the initialisation procedure for the 
!   multi-layer snow cover scheme (SNOWPOLINO).!
!
! Note:
!
!
! Current Code Owner: MeteoSwiss, Sascha Bellaire
!   phone:
!   fax:
!   email:  sascha.bellaire@meteoswiss.com
!
!
! History:
!   Version    Date       Name
!   ---------- ---------- ----
!   1.0       yyy/mm/dd Sascha Bellaire
!   Initial release, i.e. first stable version
!
!
! Code Description:
!   Language: Fortran 90.
!   Software Standards: "European Standards for Writing and
!   Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! =============================================================================
! + Begin loading modules
! =============================================================================

! General modules
#ifdef __COSMO__
USE kind_parameters, ONLY :   &
#elif __ICON__
USE mo_kind, ONLY:     &
#endif
    vpp           ! KIND-type parameter for real variables (variable precision physics)

!
USE data_runcontrol, ONLY:  &
    ntstep

! Additional snow scheme modules
USE sfc_snow_data           ! contains all relevant data and setting
USE sfc_snow_utilities      ! contains all required subroutines (see list in module header)


! Physical and other constants
USE data_constants, ONLY: &
    rho_w_wp => rho_w



! =============================================================================
! - End loading modules
! =============================================================================

  IMPLICIT NONE

  PUBLIC           ! All constants and variables in this module are public

  CONTAINS

! =============================================================================
!  + Begin subroutine  - snow_init
! =============================================================================

  SUBROUTINE snow_init (                   &
                       
                        ! general
                        nvec             , & ! array dimensions
                        ivstart          , & ! start index
                        ivend            , & ! end index
                        ke_soil          , & ! without lowermost (climat.) soil layer      US ??
                        dt               , & ! integration timestep

                        ! snow
                        t_sn_now         , & ! snow temperature (main level)             (  K  )
                        t_sn_new         , & ! snow temperature (main level)             (  K  )

                        theta_i_now      , & ! volumetric ice content                    (  -  )
                        theta_i_new      , & ! volumetric ice content                    (  -  )

                        theta_w_now      , & ! volumetric water content                  (  -  )
                        theta_w_new      , & ! volumetric water content                  (  -  )

                        theta_a_now      , & ! volumetric air content                    (  -  )
                        theta_a_new      , & ! volumetric air content                    (  -  )

                        dzm_sn_now       , & ! snow layer thickness                      (  -  )
                        dzm_sn_new       , & ! snow layer thickness                      (  -  )

                        hn_sn_now        , & ! new snow amount (storage)                 (  m  )
                        hn_sn_new        , & ! new snow amount (storage)                 (  m  )

                        top_sn_now       , & ! index of the first (top) snow layer       (  -  )
                        top_sn_new       , & ! index of the first (top) snow layer       (  -  )

                        h_snow           , & ! snow height

                        ! soil
                        t_so_now         , & ! soil temperature (main level)             (  K  )
                        t_so_new           ) ! soil temperature (main level)             (  K  ) 


                     
! =============================================================================
!
! Description:
!
! Method:
!
! Dependencies:
!
! Refernces:
!
! Notes:
!
! =============================================================================

! -----------------------------------------------------------------------------
! + Begin declarations
! -----------------------------------------------------------------------------

! ------------------------
! + Global
! -----------------------

  ! generell
  INTEGER, INTENT(IN)  ::  &
                  nvec             , & ! array dimensions
                  ivstart          , & ! start index
                  ivend            , & ! end index
                  ke_soil              ! number of soil layers

  REAL    (KIND = vpp), INTENT(IN)  ::  &
                  dt                   ! time step

  ! snow
   REAL    (KIND = vpp)    , DIMENSION(nvec), INTENT(INOUT)  :: &
                  h_snow                  ! snow height


   REAL    (KIND = vpp)    , DIMENSION(nvec,1:n_layers), INTENT(INOUT) :: &
                  t_sn_now          , &   ! snow temperature (main level)               (K)
                  theta_i_now       , &   ! volumetric ice content                      (-)
                  theta_w_now       , &   ! water ice content                           (-)
                  theta_a_now       , &   ! air ice content                             (-)
                  dzm_sn_now              ! layer thickness between main levels         (m)

   REAL    (KIND = vpp)    , DIMENSION(nvec,1:n_layers), INTENT(OUT) :: &
                  t_sn_new          , &   ! snow temperature (main level)               (K)
                  theta_i_new       , &   ! volumetric ice content                      (-)
                  theta_w_new       , &   ! water ice content                           (-)
                  theta_a_new       , &   ! air ice content                             (-)
                  dzm_sn_new              ! layer thickness between main levels         (m)

   REAL    (KIND = vpp)    , DIMENSION(nvec), INTENT(INOUT) :: &
                  hn_sn_now               ! new snow amounts (storage)                  (m)
   REAL    (KIND = vpp)    , DIMENSION(nvec), INTENT(OUT)   :: &
                  hn_sn_new

   REAL    (KIND = vpp)    , DIMENSION(nvec), INTENT(INOUT) :: &
                  top_sn_now              ! index of the first (top) layer index        (-)
   REAL    (KIND = vpp)    , DIMENSION(nvec), INTENT(OUT)   :: &
                  top_sn_new

   ! soil
   REAL    (KIND = vpp), DIMENSION(nvec,0:ke_soil+1), INTENT(INOUT) :: &
                  t_so_now             ! soil temperature (main level)                  (  K  )
   REAL    (KIND = vpp), DIMENSION(nvec,0:ke_soil+1), INTENT(OUT) :: &
                  t_so_new             ! soil temperature (main level)                  (  K  )


! ------------------------
! + Local
! -----------------------

  INTEGER      :: &

    ! Indices
    i               , & ! loop index in x-drection
    ksn                 ! loop index for snow layers

  INTEGER, DIMENSION(nvec) ::  &
    top                 ! index of first (top) snow layer

  REAL    (KIND=vpp) ::             &

    hm_sn      (nvec,1:n_layers)  , &   ! height (from bottom) of snow layers (main levels)      (m)
    zm_sn      (nvec,1:n_layers)  , &   ! depth  (from top) of snow layers (main levels)         (m)
    rho_sn     (nvec,1:n_layers)  , &   ! density of snow layers                                 (kg/m**3)
    m_sn       (nvec,1:n_layers)  , &   ! snow layer mass                                        (kg)
    hcap_sn    (nvec,1:n_layers)  , &   ! snow layer heat capacity
    hcon_sn    (nvec,1:n_layers)        ! snow layer heat conductivity


   ! physical and other  constants
   REAL    (KIND = vpp) :: &
                  rho_w                ! density of water                                        (kg / m**3)

   INTEGER (kind = vpp) :: l_top
! ------------------------------------------------------------------------------
! - End declarations
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + Begin subroutine arguments
! ------------------------------------------------------------------------------

  ! local arrays
  !$acc enter data async                                                       &
  !$acc create(top, hm_sn, zm_sn, rho_sn, m_sn, hcap_sn, hcon_sn)
  !$acc data                                                                   &
  !$acc present(top, hm_sn, zm_sn, rho_sn, m_sn, hcap_sn, hcon_sn)

  ! subroutine arguments
  !$acc data                                                                   &
  !$acc present(h_snow, t_sn_now, theta_i_now, theta_w_now, theta_a_now)       &
  !$acc present(dzm_sn_now, t_sn_new, theta_i_new, theta_w_new, theta_a_new)   &
  !$acc present(dzm_sn_new, hn_sn_now, hn_sn_new, top_sn_now, top_sn_new)      &
  !$acc present(t_so_now, t_so_new)


  ! Double precision equivalents of module variables
  rho_w   = REAL(rho_w_wp  , vpp)

  hn_sn_now = REAL(0.0 , vpp)
! ------------------------------------------------------------------------------
! + Initiate fields
! ------------------------------------------------------------------------------
dzm_sn_now = 0.0_vpp
t_sn_now = 0.0_vpp
theta_i_now = 0.0_vpp
theta_w_now = 0.0_vpp
theta_a_now = 0.9_vpp

   ! ------------------------
   ! Intiate local fields to zero
   ! -----------------------

   !$acc parallel async default(none)
   !$acc loop gang vector
   DO ksn = 1, n_layers
     DO i = ivstart, ivend

       hm_sn   (i,ksn) = 0.0_vpp
       zm_sn   (i,ksn) = 0.0_vpp
       rho_sn  (i,ksn) = 0.0_vpp
       m_sn    (i,ksn) = 0.0_vpp
       hcon_sn (i,ksn) = 0.0_vpp
       hcap_sn (i,ksn) = 0.0_vpp
       !hdif_sn (i,ksn) = 0.0_vpp
     ENDDO
   ENDDO
   !$acc end parallel

  ! --------------------------
  ! Index of top snow layer
  ! -------------------------

   !$acc parallel async default(none)
   !$acc loop gang vector
   DO i = ivstart, ivend
     ksn_loop : do ksn = 1, n_layers
        if ( dzm_sn_now(i,ksn) .le. 0.000  ) then
          exit ksn_loop
        endif
     enddo ksn_loop
     top(i) = ksn-1 !NINT( top_sn_now(i) )
     top_sn_now(i) = top(i)


   ENDDO
   !$acc end parallel

  ! --------------------------
  ! Height of snow layer (main) levels
  ! -------------------------
 
   !$acc parallel async default(none)
   !$acc loop gang vector private(l_top)
   DO i = ivstart, ivend
      l_top = top(i)
      DO ksn = 1, l_top

        IF(ksn .EQ. 1) THEN
       
          hm_sn(i,ksn) = dzm_sn_now(i,ksn)
      
        ELSE

          hm_sn(i,ksn) = hm_sn(i,ksn-1) + dzm_sn_now(i,ksn)

        ENDIF

      ENDDO
   ENDDO 
   !$acc end parallel


  ! --------------------------
  ! Depth of snow layer (main) levels
  ! -------------------------

   !$acc parallel async default(none)
   !$acc loop gang vector private(l_top)
   DO i = ivstart, ivend
      l_top = top(i)
      DO ksn = l_top , 1 , -1
       zm_sn(i, (l_top+1) - ksn) = hm_sn(i,ksn) !invert height vector

     ENDDO
   ENDDO
   !$acc end parallel


  ! --------------------------
  ! Snow layer density 
  ! -------------------------

   !$acc parallel async default(none)
   !$acc loop gang vector private(l_top)
   DO i = ivstart, ivend
      l_top = top(i)
      DO ksn = 1, l_top, 1

        IF(theta_i_now(i,ksn) .EQ. 0.0_vpp) THEN

          rho_sn(i,ksn) = 0.0_vpp                              

        ELSE

          rho_sn(i,ksn) = theta_i_now(i,ksn)*rho_i + theta_w_now(i,ksn)*rho_w

        ENDIF

     ENDDO
   ENDDO
   !$acc end parallel


  ! --------------------------
  ! Heat capacity
  ! --------------------------

   !$acc parallel async default(none)
   !$acc loop gang vector private(l_top)
   DO i = ivstart, ivend
     l_top = top(i)
     DO ksn = 1, l_top, 1
       IF(rho_sn(i,ksn) .LT. eps) THEN

         hcap_sn(i,ksn) = 0.0_vpp

       ELSE

         hcap_sn(i,ksn) = (  rho_a   * theta_a_now(i,ksn) * specific_heat_air     &
                           + rho_i   * theta_i_now(i,ksn) * specific_heat_ice     &
                           + rho_w   * theta_w_now(i,ksn) * specific_heat_water ) &
                           / rho_sn(i,ksn)

       ENDIF

     ENDDO
   ENDDO
   !$acc end parallel


  ! --------------------------
  ! Heat conductivity
  ! --------------------------

   !$acc parallel async default(none)
   !$acc loop gang vector private(l_top)
   DO i = ivstart, ivend
    l_top = top(i)
    do ksn = 1,l_top,1
      if(rho_sn(i,ksn) .lt. eps) then
           hcon_sn(i,ksn) = 0.0_vpp
      else
           hcon_sn(i,ksn) = 2.22_vpp * EXP(1.88_vpp * LOG(rho_sn(i,ksn)/rho_i))
      endif
     ENDDO
   ENDDO
   !$acc end parallel

   
  ! --------------------------
  ! Snow layer mass
  ! --------------------------

   !$acc parallel async default(none)
   !$acc loop gang vector private(l_top)
   DO i = ivstart, ivend
    l_top = top(i)
    do ksn = 1,l_top,1

     IF(ksn .EQ. 1) THEN

        m_sn(i,ksn) = zm_sn(i,ksn) * ((theta_i_now(i,ksn) * rho_i) + (theta_w_now(i,ksn) * rho_w))

      ELSE

        m_sn(i,ksn) = ABS(zm_sn(i,ksn) - zm_sn(i,ksn-1))  * ((theta_i_now(i,ksn) * rho_i) + (theta_w_now(i,ksn) * rho_w))

      ENDIF

     ENDDO
   ENDDO
   !$acc end parallel


! ------------------------------------------------------------------------------
! - End subroutine arguments
! ------------------------------------------------------------------------------

  !$acc end data
  !$acc end data

  !$acc exit data async                                                        &
  !$acc delete(top, hm_sn, zm_sn, rho_sn, m_sn, hcap_sn, hcon_sn)


  END SUBROUTINE snow_init

! =============================================================================
!  - End subroutine  - snow_init
! =============================================================================





! DONE, DONE!!!
END MODULE sfc_snow_init

!==============================================================================
! + End of module sfc_snow_init
!==============================================================================











