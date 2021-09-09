!+ Source module "gribio_eccodes" - Subroutines for input and output of Grib Code
!+  (libDWD(GRIB1) and ecCodes-Library (GRIB1 and GRIB2), ECMWF)
!------------------------------------------------------------------------------

MODULE tsa_gribio

!------------------------------------------------------------------------------
!
!! Description:
! contains subroutines for processing of Grib:
! initialyzes, resets, adds fields and clears "griblist"-
! a data structure which points to all data fields
! and information on the PDS for GRIB1 (libDWD) or sname etc
! for eccodes GRIB1 and GRIB2.
! read_grib (read_grib_eccodes) reads a GRIB file according to "griblist"
!
!    Currently included:
!    - init_griblist:
!      initialyzes a griblist
!    - clear_griblist:
!      clears (deletes) a griblist
!    - reset_griblist:
!      resets (empties) a griblist
!    - add2griblist2D/add2griblist2D_grb1:
!      adds field to a 2 Dimensional griblist
!    - add2griblist3D/add2griblist3D_grb1:
!      adds field to a 3 Dimensional griblist
!
!    - get_gribinfo --> now get_grib_info!
!      reads GRIB info
!    - get_grib_info (NEW)
!      reads GRIB info with eccodes library from ECMWF
!    - get_grib_info_icon (NEW)
!      reads GRIB info of ICON fields with eccodes library from ECMWF
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
! V1.0       2006/10/01 Felix Ament, MeteoSwiss
!  Initial Version
! V_TSA      2015/12/01 Yiftach Ziv, IMS
!  Using KIND definitions from data_parameters
!  Changed code to adhere to coding standards
! 5.07       2020/02/21     Doerte Liermann, Ulrich Schaettler
!  Replaced data_parameters by kind_parameters and removed iintegers
! V5.8       2021-02-01     Varun Sharma
!  Implementation of snowpolino in TSA
!  bug fix : added a missing grib_release in read_grib_eccodes ( line 1078 )
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

USE kind_parameters, ONLY :   &
    wp           ! KIND-type parameter for real variables

USE data_io, ONLY :   &
    intgribc                , & ! Kind type for C integer variables used in the GRIB library
    intgribf                , & ! Kind type for Fortran integer variables used in the GRIB librar
    irealgrib               , & ! Kind type for Fortran real variables used in the GRIB library
    iwlength                    ! length of integers used in the griblib in byte
                                !   4: for dwdlib on SGI systems and griblib on all systems
!DL
USE tsa_data,              ONLY : &
    rundef                      ! Use rundef to preset the grib fields
!DL

!==============================================================================

USE grib_api
IMPLICIT NONE

!==============================================================================

!DL griblist_type defines all the things needed to deal with GRIB libDWD and eccodes

TYPE griblist_type
    INTEGER(KIND=intgribf)          :: pds(20)
!!  INTEGER(KIND=intgribf)          :: gds(20)
  CHARACTER (LEN=30)              :: sname
  CHARACTER (LEN=30)              :: levtype
  INTEGER                         :: toplev
  INTEGER                         :: botlev
  CHARACTER (LEN=30)              :: stattype
  REAL   (KIND=wp), POINTER       :: p3(:,:,:)
  REAL   (KIND=wp), POINTER       :: p2(:,:)
END TYPE griblist_type


TYPE (griblist_type), ALLOCATABLE :: griblist(:)
  INTEGER                         :: &
      nvarmax                      , & ! Maximum number of Grib-fields
      ie_grib                      , & ! Number of gridpoints in x-direction
      je_grib                      , & ! Number of gridpoints in y-direction
      ke_grib                          ! Number of model vertical layers
  REAL (KIND=wp)    , ALLOCATABLE :: vcoord(:)
  REAL (KIND=wp)                  :: &
      p0sl                         , & !
      t0sl                         , & !
      dt0lp                        , & !
      vcflat                           !
  INTEGER                         :: ivctype  ! Type of vertical coordinates


!===============================================================
! Interface Blocks
!DL NEW: subroutines for libDWD and eccodes!!
INTERFACE add2griblist
   MODULE PROCEDURE        &
     add2griblist2d,       &
     add2griblist2d_grb1,  &
     add2griblist3d,       &
     add2griblist3d_grb1
END INTERFACE

!========================================================================

CONTAINS

!========================================================================
SUBROUTINE init_griblist(nvarmax2,ie2,je2,ke2)
!
!   Initializes griblist for further read/WRITE
!
   IMPLICIT NONE
!========================================================================
   INTEGER                  :: &
       nvarmax2              , & ! Maximum number of Grib-fields
       ie2                   , & ! Number of gridpoints in x-direction
       je2                   , & ! Number of gridpoints in y-direction
       ke2                       ! Number of vertical model layers

   nvarmax=nvarmax2
   ie_grib=ie2
   je_grib=je2
   ke_grib=ke2
   ALLOCATE(griblist(nvarmax))
   ALLOCATE(vcoord(ke_grib+1))
   CALL reset_griblist()
END SUBROUTINE init_griblist

!========================================================================
SUBROUTINE clear_griblist()
!
!   clears griblist
!
   IMPLICIT NONE

!========================================================================
   DEALLOCATE(griblist,vcoord)
END SUBROUTINE clear_griblist


!========================================================================
SUBROUTINE reset_griblist()
!
!   resets griblist
!
   IMPLICIT NONE
!========================================================================
   INTEGER                  :: i
   DO i=1,nvarmax
   griblist(i)%pds=-1
   griblist(i)%sname='unknown'
   griblist(i)%levtype='unknown'
   griblist(i)%toplev=-1
   griblist(i)%botlev=-1
   griblist(i)%stattype='unknown'
   NULLIFY(griblist(i)%p2,griblist(i)%p3)
   ENDDO
END SUBROUTINE reset_griblist


!========================================================================
SUBROUTINE add2griblist2d_grb1(field,ee,tab,lty,lvt,lvb,ntri)
!
!   adds fields to 2 dimensional griblist - if GRIB1 libDWD
!
   IMPLICIT NONE

!========================================================================

   INTEGER, INTENT(IN) :: &
      ee,  & ! Element Identifier (GRIB)
      tab, & ! Table Number       (GRIB)
      lty    ! Leveltype          (GRIB)
   INTEGER, INTENT(IN), OPTIONAL :: &
      lvt, & ! Top Level (GRIB)
      lvb, & ! Bottom Level (GRIB)
      ntri
   REAL (KIND=wp), TARGET ::    field(:,:) ! datafield to be processed

   INTEGER :: i
   LOGICAL :: ladded

   IF (.NOT.ALLOCATED(griblist)) THEN
       WRITE(6,*) "Error: Griblist is not initialized!"
      STOP
   ENDIF

   ladded=.FALSE.
   DO i=1,nvarmax
      IF (.NOT.(ASSOCIATED(griblist(i)%p2).OR.ASSOCIATED(griblist(i)%p3))) THEN
         griblist(i)%p2=>field
         griblist(i)%pds(7)=ee
         griblist(i)%pds(2)=tab
         griblist(i)%pds(8)=lty
         IF (PRESENT(lvt))  THEN
            IF (lvt>=0) griblist(i)%pds(9)=lvt
         ENDIF
         IF (PRESENT(lvb)) THEN
            IF (lvb>=0) griblist(i)%pds(10)=lvb
         ENDIF
         IF (PRESENT(ntri)) THEN
            griblist(i)%pds(19)=ntri
         ENDIF
         ladded=.TRUE.
         EXIT
      ENDIF
   ENDDO
   IF (.NOT.ladded) THEN
      WRITE(6,*) "Error: nvarmax to small!"
      STOP
   ENDIF
END SUBROUTINE add2griblist2d_grb1


!========================================================================
SUBROUTINE add2griblist3d_grb1(field,ee,tab,lty,lvt)
!
!   adds fields to 3 dimensional griblist - if GRIB1 libDWD
!
   IMPLICIT NONE

!========================================================================
  INTEGER, INTENT(IN) :: &
      ee,  & ! Element Identifier (GRIB)
      tab, & ! Table Number       (GRIB)
      lty    ! Leveltype          (GRIB)
   INTEGER, INTENT(IN), OPTIONAL :: &
      lvt    ! Top Level (GRIB)
   REAL (KIND=wp), TARGET ::    field(:,:,:) ! datafield to be processed

   INTEGER :: j
   LOGICAL :: ladded

   IF (.NOT.ALLOCATED(griblist)) THEN
      WRITE(6,*) "Error: Griblist is not initialized!"
      STOP
   ENDIF

   ladded=.FALSE.
   DO j=1,nvarmax
      IF (.NOT.(ASSOCIATED(griblist(j)%p2).OR.ASSOCIATED(griblist(j)%p3))) THEN
         griblist(j)%p3=>field
         griblist(j)%pds(7)=ee
         griblist(j)%pds(2)=tab
         griblist(j)%pds(8)=lty
         IF (PRESENT(lvt)) THEN
            IF (lvt>=0) THEN
               griblist(j)%pds(9)=lvt
            ENDIF
         ENDIF
         ladded=.TRUE.
         EXIT
      ENDIF
   ENDDO
   IF (.NOT.ladded) THEN
      WRITE(6,*) "Error: nvarmax to small!"
      STOP
   ENDIF

END SUBROUTINE add2griblist3d_grb1

!========================================================================

SUBROUTINE read_grib(dateiname,i02,j02,di2,dj2,ipds18,anyundef)
!
!   Reads GRIB according to griblist - GRIB1 libDWD
!
   IMPLICIT NONE


!========================================================================

   CHARACTER (LEN=*), INTENT(IN)                   :: dateiname
   INTEGER, INTENT(IN) , OPTIONAL  ::  &
       i02,j02,di2,dj2
   INTEGER, INTENT(OUT), OPTIONAL :: ipds18
   LOGICAL, INTENT(OUT), OPTIONAL                  :: anyundef


   INTEGER (KIND=intgribf),  PARAMETER :: &
     npds   =    321_intgribf, & ! Dimension for product definition section (pds)
     ngds   =    626_intgribf, & ! Dimension for grid description section (gds)
     nbms   =      3_intgribf, & ! Dimension for bit map section (bms)
     nbds   =     11_intgribf, & ! Dimension for binary data section
     ndsup  =     73_intgribf, & ! Dimension for dsup
     ndims  =     20_intgribf    ! Dimension for idims (contains all dimensions)

   REAL    (KIND=irealgrib), PARAMETER :: &
     undefgrib =  -1e7_irealgrib ! value for "undefined" in the grib routines
   INTEGER (KIND=intgribf)             :: &
     iednr     =   1             ! grib edition number 

! The following dimensions are set during program execution
   INTEGER (KIND=intgribf)             :: &
     lfd,                      & ! Dimension for iblock
     lbm,                      & ! Dimension for bitmap: this value has to be at
                              ! least the same as lbmax in SUBROUTINE grbin1
     lds                         ! Dimension for unpacked data

! Global Arrays:
   INTEGER (KIND=intgribf)              :: &
     idims     (ndims),  & ! array for all dimensions
     ipds      (npds),   & ! product definition section
     igds      (ngds),   & ! grid description section 
     ibms      (nbms),   & ! bit map section
     ibds      (nbds)      ! binary data section

! Arrays that are ALLOCATED during program execution
   INTEGER (KIND=intgribf), ALLOCATABLE :: &
     iblock    (:),      & ! array for gribed data
     ibmap     (:)         ! array for bit map

   REAL   (KIND=irealgrib), ALLOCATABLE :: &
     dsup (:),       & ! array for special data
     ds   (:),       & ! array for unpacked data
     dstemp(:,:)



   INTEGER                              :: &
!   INTEGER              :: &
      nelement, nlevtyp, iegb, jegb, kegb, &
      nlevel, nyear, nmonth, nday, ntime,  &
      ntflag, nvv, nvv1, nvv2, iz_countl,  &
      ip,nv,k, i0,j0,di,dj,i
   REAL (KIND=wp)                       :: &
      philugb, rlalugb, dphigb,            &
      dlamgb, refstf
   INTEGER (KIND=intgribc)              :: &
      nudatc=100,                          &
      maxlenc,lenc,ierrc
   INTEGER (KIND=intgribf)              :: ierrf
   LOGICAL                              :: leof,found

   LOGICAL, ALLOCATABLE                 :: lvarfound(:)

! Beginning of subroutine
!--------------------------

!!!>dbg XYZ
   WRITE(6,*) "gribFile: ",TRIM(dateiname)
!!!<dbg XYZ

   IF (.NOT.ALLOCATED(griblist)) THEN
      WRITE(6,*) "Error: Griblist is not initialized!"
      STOP
   ENDIF

   IF (PRESENT(i02)) THEN
      i0=i02
      j0=j02
      di=di2
      dj=dj2
   ELSE
      i0=1
      j0=1
      di=ie_grib
      dj=je_grib
   ENDIF
   IF (PRESENT(ipds18)) THEN
      ipds18=0
   ENDIF

   IF ( PRESENT(anyundef) ) anyundef = .FALSE.    ! Keep track of undef values

   lds = ie_grib * je_grib
   lbm = 1875
   lfd = lds * 2 / iwlength + 2000   ! the "2" means 2 bytes per word
                                     ! the "2000" is just for safety
   ALLOCATE (iblock(lfd), ibmap(lbm) )
   ALLOCATE (ds(ie_grib*je_grib), dsup(ndsup)  )
   ALLOCATE (lvarfound(nvarmax));
   DO nv=1,nvarmax
      lvarfound(nv)=.NOT.(ASSOCIATED(griblist(nv)%p2).OR.ASSOCIATED(griblist(nv)%p3))
   ENDDO

! moving arraydimensions into idims
   idims( 1)   = npds
   idims( 2)   = ngds
   idims( 3)   = nbms
   idims( 4)   = nbds
   idims( 5)   = lbm
   idims( 6)   = ndsup
   idims( 7)   = lds
   idims( 8)   = lfd
   idims(9:20) = 0

!open the grib file with copen
   CALL copen(nudatc,TRIM(dateiname),'r  ',ierrc)
   IF (ierrc/=0) THEN
    WRITE(6,*) "Can't open file ",TRIM(dateiname),"! Program aborted!"
    STOP
   ENDIF

! initial levelcounter and EOF
   iz_countl = 0
   leof      = .FALSE.

! reading the data records
! endless loop until EOF is reached

   DO
      ! Read one grib record with cuegin ( in DWD-libgrib1)
      ! Input  to cuegin: internal file descriptor, max. length of data array, array to be filled  
      ! Output of cuegin: length of data actually read, error flag
      maxlenc=lfd*iwlength
      ierrc=0
      CALL cuegin(nudatc,maxlenc,iblock,lenc,ierrc)
      iz_countl = iz_countl + 1
      IF (lenc == 0) THEN
         leof = .TRUE.   ! EOF reached
      ENDIF
      ! decrypt read bitstream, if we got any data
      IF (ierrc/=0) THEN
         WRITE(6,*) "Error at cuegin. Program aborted! Error code: ",ierrc
         STOP
      ENDIF
!XYZ>   substituting GOTO with EXIT   
      IF (leof) THEN!GOTO 999
         EXIT
      ENDIF
!XYZ<
      ds = 0.0_irealgrib

      CALL grbin1(iednr, undefgrib, ndims, idims, iblock, ibmap,         &
           ipds, igds, ibms, ibds, dsup, ds, ierrf)
      IF (ierrf /= 0) THEN
         STOP
      ENDIF


      !     getting sigma from gdb
      IF ((ke_grib>1).AND.(iz_countl == 1).AND.(igds(26)/=-999999)) THEN
         p0sl = refstf(igds(26))
         t0sl = refstf(igds(27))
         dt0lp = refstf(igds(28))
         vcflat = refstf(igds(29))
         DO k = 1,ke_grib+1
           vcoord(k) = refstf(igds(29 + k))
         ENDDO

      ! Check for the type and consistency of the vertical coordinate parameters
         IF ( vcoord(2) > vcoord(1) ) THEN
            ivctype = 1
         ELSEIF ( vcoord(2) < vcoord(1) ) THEN
            ivctype = 2
         ELSE
            WRITE(6,*) ' ERROR *** Type ivctype of vertical coordinate not available*** '
            STOP
          ENDIF
      ENDIF

      !     Reading the meshsizes (DPHIGB, DLAMGB) and the lower left
      !     corner out of the Grib-Description-Block (IGDB)
      PHILUGB  =  IGDS( 7)*0.001
      RLALUGB  =  IGDS( 8)*0.001
      DPHIGB   =  IGDS(13)*0.001
      DLAMGB   =  IGDS(12)*0.001

      !     Reading the number of gridpoints (IE, JE) out of the IGDB
      IEGB     =  IGDS( 5)
      JEGB     =  IGDS( 6)

      !     Getting the number of layers (KE)
      !     Only for data on model levels, not for pressure level data
      NLEVTYP  = IPDS( 8)
      IF (NLEVTYP.NE.100 .OR. NLEVTYP.NE.102) THEN
         KEGB  = (IGDS(1) - 42)/4 - 4
      ELSE
         KEGB  = ke_grib
      ENDIF

      !     Reading the element, the level typ, level, date and time
      !     out of the Product-Definition-Block (IPDS)
      NELEMENT = IPDS( 7)

      NLEVTYP  = IPDS( 8)
      IF (NLEVTYP.NE.100 .OR. NLEVTYP.NE.102) THEN
         NLEVEL= IPDS( 9)
      ELSE
         NLEVEL= IPDS(10)
      ENDIF

      IF (NLEVTYP.EQ.109) THEN
         NLEVEL= IPDS(10)
      ENDIF

      NYEAR    = IPDS(11)
      NMONTH   = IPDS(12)
      NDAY     = IPDS(13)
      NTIME    = IPDS(14)
      NTFLAG   = IPDS(19)
      IF ( NTFLAG.EQ.10 ) THEN
         NVV    = IPDS(18)
      ELSE IF ( NTFLAG.EQ.3 .OR. NTFLAG.EQ.4 ) THEN
         NVV1   = IPDS(17)
         NVV2   = IPDS(18)
      ENDIF

      IF (PRESENT(ipds18)) THEN
         IF (ipds18<ipds(18)) THEN
            ipds18=ipds(18)
         ENDIF
      ENDIF

      DO nv=1,nvarmax
         IF (ASSOCIATED(griblist(nv)%p2).OR.ASSOCIATED(griblist(nv)%p3)) THEN
           found=.TRUE.
            DO ip=1,20
               IF (griblist(nv)%pds(ip)>=0) THEN
                  IF (griblist(nv)%pds(ip)/=ipds(ip)) THEN
                     found=.FALSE.
                  ENDIF
               ENDIF
            ENDDO
            IF (found) THEN
               lvarfound(nv)=.TRUE.
               IF ( ANY(ds(:ie_grib*je_grib) == undefgrib).AND.PRESENT(anyundef) ) THEN
                  anyundef = .TRUE.
               ENDIF
               IF (PRESENT(i02)) THEN
                  ALLOCATE (dstemp(ie_grib,je_grib))
                  dstemp=RESHAPE(ds(1:ie_grib*je_grib),(/ ie_grib, je_grib /))
                  IF (ASSOCIATED(griblist(nv)%p2)) THEN
                     griblist(nv)%p2(:,:)=dstemp(i0:i0+di-1,j0:j0+dj-1)
                  ELSE
                     IF (ipds(8)==109) THEN
                        griblist(nv)%p3(:,:,ipds(10))=dstemp(i0:i0+di-1,j0:j0+dj-1)
                     ELSE
                        griblist(nv)%p3(:,:,ipds(9))=dstemp(i0:i0+di-1,j0:j0+dj-1)
                     ENDIF
                  ENDIF
                  DEALLOCATE(dstemp)
               ELSE
                  IF (ASSOCIATED(griblist(nv)%p2)) THEN
                     griblist(nv)%p2(:,:)=RESHAPE(ds(1:ie_grib*je_grib),(/ ie_grib, je_grib /))
                  ELSE
                     IF (ipds(8)==109) THEN
                        griblist(nv)%p3(:,:,ipds(10))=RESHAPE(ds(1:ie_grib*je_grib),(/ ie_grib, je_grib /))
                     ELSE
                        griblist(nv)%p3(:,:,ipds(9))=RESHAPE(ds(1:ie_grib*je_grib),(/ ie_grib, je_grib /))
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO
   ENDDO !read_loop

!   999 CONTINUE  XYZ> cancellation of GOTO 999

   IF (COUNT(lvarfound)<nvarmax) THEN
      WRITE(6,*) "WARNING: The following fields of griblist are not available:"
       DO nv=1,nvarmax
          IF (.NOT.lvarfound(nv)) &
               WRITE(6,*) "!!! ee:",griblist(nv)%pds(7),"lty",griblist(nv)%pds(8),"tab:",griblist(nv)%pds(2)
       ENDDO
   ENDIF

   CALL cclose(nudatc,'exi',ierrc)
   IF (ierrc /= 0) THEN
      WRITE(6,*) 'Error in cclose!'
      STOP
   ENDIF

   DEALLOCATE (iblock,ibmap,lvarfound )
   DEALLOCATE (dsup,ds )

END SUBROUTINE read_grib

!==============================================================================
!========================================================================
SUBROUTINE add2griblist2d(field,sname,lty,lvt,lvb,stattype)
!
!   adds fields to 2 dimensional griblist
!
   IMPLICIT NONE

!========================================================================
!! INTEGER, INTENT(IN) :: &
!!    ee,  & ! Element Identifier (GRIB)
!!    tab, & ! Table Number       (GRIB)
!!    lty    ! Leveltype          (GRIB)
   CHARACTER (LEN=*),       INTENT(IN) :: &
      sname, & ! shortName (GRIB)
      lty      ! Leveltype (GRIB)

   INTEGER, INTENT(IN), OPTIONAL :: &
      lvt, &   ! Top Level (GRIB)
      lvb      ! Bottom Level (GRIB)
!!    ntri 
   CHARACTER (LEN=*),       INTENT(IN), OPTIONAL :: &
      stattype ! stepType (GRIB)

   REAL (KIND=wp), TARGET ::    field(:,:) ! datafield to be processed

   INTEGER :: i
   LOGICAL :: ladded

   IF (.NOT.ALLOCATED(griblist)) THEN
       WRITE(6,*) "Error: Griblist is not initialized!"
      STOP
   ENDIF

   ladded=.FALSE.
   DO i=1,nvarmax
      IF (.NOT.(ASSOCIATED(griblist(i)%p2).OR.ASSOCIATED(griblist(i)%p3))) THEN
         griblist(i)%p2=>field
         griblist(i)%sname=sname
         griblist(i)%levtype=lty
         IF (PRESENT(lvt))  THEN 
            IF (lvt>=0) griblist(i)%toplev=lvt
         ENDIF
         IF (PRESENT(lvb)) THEN
            IF (lvb>=0) griblist(i)%botlev=lvb
         ENDIF
         IF (PRESENT(stattype)) THEN
            griblist(i)%stattype=stattype
         ENDIF
         ladded=.TRUE.
         EXIT
      ENDIF
   ENDDO
   IF (.NOT.ladded) THEN
      WRITE(6,*) "Error: nvarmax to small!"
      STOP
   ENDIF
END SUBROUTINE add2griblist2d


!========================================================================
SUBROUTINE add2griblist3d(field,sname,lty,lvt)
!
!   adds fields to 3 dimensional griblist
!
   IMPLICIT NONE

!========================================================================
  CHARACTER (LEN=30),       INTENT(IN) :: &
      sname, & ! shortName          (GRIB)
      lty      ! Leveltype          (GRIB)
!!INTEGER, INTENT(IN) :: &
!!    ee,  & ! Element Identifier (GRIB)
!!    tab, & ! Table Number       (GRIB)
!!    lty    ! Leveltype          (GRIB)
   INTEGER, INTENT(IN), OPTIONAL :: &
      lvt    ! Top Level (GRIB)
   REAL (KIND=wp), TARGET ::    field(:,:,:) ! datafield to be processed

   INTEGER :: j
   LOGICAL :: ladded

   IF (.NOT.ALLOCATED(griblist)) THEN
      WRITE(6,*) "Error: Griblist is not initialized!"
      STOP
   ENDIF

   ladded=.FALSE.
   DO j=1,nvarmax
      IF (.NOT.(ASSOCIATED(griblist(j)%p2).OR.ASSOCIATED(griblist(j)%p3))) THEN
         griblist(j)%p3=>field
         griblist(j)%sname=sname
         griblist(j)%levtype=lty
         IF (PRESENT(lvt)) THEN
            IF (lvt>=0) THEN
               griblist(j)%toplev=lvt
            ENDIF
         ENDIF
         ladded=.TRUE.
         EXIT
      ENDIF
   ENDDO
   IF (.NOT.ladded) THEN
      WRITE(6,*) "Error: nvarmax to small!"
      STOP
   ENDIF
   
END SUBROUTINE add2griblist3d


!========================================================================

!========================================================================
SUBROUTINE read_grib_eccodes(dateiname,i02,j02,di2,dj2,ivv,anyundef,lhourly)
!
!   Reads GRIB according to griblist with eccodes from ECMWF
!
   USE grib_api

   IMPLICIT NONE


!========================================================================

   CHARACTER (LEN=*), INTENT(IN)                   :: dateiname 
   INTEGER, INTENT(IN), OPTIONAL  ::  &
       i02,j02,di2,dj2
   INTEGER, INTENT(OUT), OPTIONAL :: ivv
   LOGICAL, INTENT(OUT), OPTIONAL                  :: anyundef, lhourly

   REAL    (KIND=irealgrib), PARAMETER :: &
     undefgrib =  -1e7_irealgrib ! value for "undefined" in the grib routines
   INTEGER (KIND=intgribf)             :: &
     iednr                       ! grib edition number 

   INTEGER (KIND=intgribf)             :: &
     lds                         ! Dimension for unpacked data

! Local variables for grib_get calls
   CHARACTER (LEN=30)                  :: sname, levtype, stype, gridtype
   INTEGER                             :: &
     botlev, toplev, startstep, endstep,  &
     nkey
   REAL    (KIND=irealgrib)            :: &
     zbotlev, ztoplev

! Arrays that are ALLOCATED during program execution
!! INTEGER (KIND=intgribf), ALLOCATABLE :: &
!!   iblock    (:)       & ! array for gribbed data

   REAL   (KIND=irealgrib), ALLOCATABLE :: &
     ds   (:),       & ! array for unpacked data
     dstemp(:,:)


      
   INTEGER                              :: &
!   INTEGER              :: &
      nelement, nlevtyp, iegb, jegb, kegb, &
      nlevel, nyear, nmonth, nday, ntime,  &
      ntflag, nvv, nvv1, nvv2, iz_countl,  &
      ip,nv,k, i0,j0,di,dj,i

!  Handling of missing fields
   INTEGER                              :: &
    imiss_soil, imiss_tg, imiss_tsnow, imiss_wsnow, imiss_ts

!  Check if lhourly_data is true
   INTEGER                             :: &
     ntime_int(20), nno

   REAL (KIND=wp)                       :: &
      philugb, rlalugb, dphigb,            &
      dlamgb, refstf
   INTEGER                               :: &
      nudat,                              &
      igrib,ierr
   INTEGER (KIND=intgribf)              :: ierrf
   LOGICAL                              :: leof,found

   LOGICAL, ALLOCATABLE                 :: lvarfound(:)

! Beginning of subroutine
!--------------------------

!!!>dbg XYZ
   WRITE(6,*) "gribFile: ",TRIM(dateiname)
!!!<dbg XYZ

   IF (.NOT.ALLOCATED(griblist)) THEN
      WRITE(6,*) "Error: Griblist is not initialized!"
      STOP
   ENDIF

   IF (PRESENT(i02)) THEN
      i0=i02
      j0=j02
      di=di2
      dj=dj2
   ELSE
      i0=1
      j0=1
      di=ie_grib
      dj=je_grib
   ENDIF
   IF (PRESENT(ivv)) THEN
      ivv=0
   ENDIF

   IF ( PRESENT(anyundef) ) anyundef = .FALSE.    ! Keep track of undef values
   IF ( PRESENT(lhourly)  ) lhourly  = .FALSE.    ! No hourly data of precipitation or radiation

!!  ALLOCATE (iblock(lfd), ibmap(lbm) )
   ALLOCATE (ds(ie_grib*je_grib))


   ALLOCATE (lvarfound(nvarmax));
   DO nv=1,nvarmax
      lvarfound(nv)=.NOT.(ASSOCIATED(griblist(nv)%p2).OR.ASSOCIATED(griblist(nv)%p3))
   ENDDO 


! open the grib file with grib_open_file
   CALL grib_open_file(nudat,TRIM(dateiname),'r  ',ierr)
   IF (ierr/=0) THEN
    WRITE(6,*) "Can't open file ",TRIM(dateiname),"! Program aborted!"
    STOP
   ENDIF

! initial levelcounter and EOF
   iz_countl = 0
   leof      = .FALSE.

! presetings for checking hourly data
  ntime_int(:) = -99
  nno          = 0

! reading the data records
! endless loop until EOF is reached

   DO
      ! Read one grib record with  grib_new_from_file (eccodes)
      ierr=0
      CALL grib_new_from_file(nudat,igrib,ierr)
      iz_countl = iz_countl + 1 
      IF (ierr == GRIB_END_OF_FILE) THEN 
         leof = .TRUE.   ! EOF reached
      ENDIF
      IF (.NOT. leof .AND. (ierr /= GRIB_SUCCESS)) THEN
        WRITE(6,*) "read_grib_eccodes: Error with grib_new_from_file. Program aborted! Error code: ",ierr
        STOP
      ENDIF
!XYZ>   substituting GOTO with EXIT   
      IF (leof) THEN!GOTO 999
         EXIT
      ENDIF
!XYZ<
      
      CALL grib_get (igrib, 'edition',     iednr)
      CALL grib_get (igrib, 'shortName',   sname)
      CALL grib_get (igrib, 'typeOfLevel', levtype)
      CALL grib_get (igrib, 'bottomLevel', botlev)
      CALL grib_get (igrib, 'topLevel',    toplev)
      IF (iednr == 2 .AND. sname == 'W_SO') THEN
!!DL     CALL grib_get (igrib, 'bottomLevel', zbotlev)
!!DL     botlev = NINT (zbotlev * 100._wp)
         CALL grib_get (igrib, 'topLevel',    ztoplev)
         toplev = NINT (ztoplev * 100._wp)
      ENDIF
      IF (iednr == 2 .AND. sname == 'T_SO') THEN
         CALL grib_get (igrib, 'level', ztoplev)
         toplev = NINT ((ztoplev + 0.0005_wp) * 100._wp)
         botlev = toplev
      ENDIF
      CALL grib_get (igrib, 'stepType',    stype)
      CALL grib_get (igrib, 'startStep',   startstep)
!!DL
      CALL grib_set (igrib, 'stepUnits',   'h')
!!DL
      CALL grib_get (igrib, 'endStep',     endstep)
      CALL grib_get (igrib, 'gridType',    gridtype)


! Take care for special cases in GRIB1 / GRIB2, so that griblist will match
! -------------------------------------------------------------------------
! Identify, if stype in griblist is 'analysis' or 'nudging' and NOT statistical process
      IF (stype == 'instant') THEN
         SELECT CASE (iednr)
         CASE (1)
           CALL grib_get (igrib, 'timeRangeIndicator',nkey)
           IF (nkey == 13) stype = 'nudging'
           IF (nkey ==  0) stype = 'analysis'
         CASE (2)
           CALL grib_get (igrib, 'typeOfGeneratingProcess',nkey)
           IF (nkey == 202 .OR. nkey == 203) stype = 'nudging'     ! could be nudging or nudgecast
           IF (nkey ==  0)  stype = 'analysis'
           IF (nkey ==  2)  stype = 'forecast'                     ! ??? forecast used for ICON 
         END SELECT
      ENDIF
! Set correct level type for model levels - GRIB1 coding (hybrid/hybridlayer) required in griblist
!!    IF (levtype == 'generalVertical')      levtype='hybrid'
!!    IF (levtype == 'generalVerticalLayer') levtype='hybridLayer'
! Set correct level type for model levels - GRIB2 coding (generalVertical, generalVerticalLayer) required in griblist
      IF (levtype == 'hybrid')      levtype='generalVertical'
      IF (levtype == 'hybridLayer') levtype='generalVerticalLayer'
! Take care for W_SO in GRIB1 (cm, level) and GRIB2 (m, layer), griblist definition is in GRIB1
      IF (sname == 'W_SO' .AND. iednr == 2) THEN
         SELECT CASE (toplev)
         CASE (0)
           toplev = 1
         CASE (1)
           toplev = 2
         CASE (3)
           toplev = 6
         CASE (9)
           toplev = 18
         CASE (27)
           toplev = 54
         CASE (81)
           toplev = 162
         CASE (243)
           toplev = 486
         CASE (729)
           toplev = 1458
         CASE DEFAULT
           WRITE(6,*) "read_grib_eccodes: Unknown W_SO level in GRIB2: ",toplev
           STOP
         END SELECT
         botlev = toplev
         levtype= 'depthBelowLand'
      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    part of old subroutine saved via comments - for later use ??
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !     getting sigma from gdb
!!    IF ((ke_grib>1).AND.(iz_countl == 1).AND.(igds(26)/=-999999)) THEN
!!       p0sl = refstf(igds(26))
!!       t0sl = refstf(igds(27))
!!       dt0lp = refstf(igds(28))
!!       vcflat = refstf(igds(29))      
!!       DO k = 1,ke_grib+1
!!         vcoord(k) = refstf(igds(29 + k))
!!       ENDDO
     
      ! Check for the type and consistency of the vertical coordinate parameters
!!       IF ( vcoord(2) > vcoord(1) ) THEN 
!!          ivctype = 1
!!       ELSEIF ( vcoord(2) < vcoord(1) ) THEN
!!          ivctype = 2
!!       ELSE
!!          WRITE(6,*) ' ERROR *** Type ivctype of vertical coordinate not available*** '
!!          STOP
!!        ENDIF
!!    ENDIF
   
      !     Reading the meshsizes (DPHIGB, DLAMGB) and the lower left
      !     corner out of the Grib-Description-Block (IGDB)
!!    PHILUGB  =  IGDS( 7)*0.001
!!    RLALUGB  =  IGDS( 8)*0.001
!!    DPHIGB   =  IGDS(13)*0.001
!!    DLAMGB   =  IGDS(12)*0.001
 
      !     Reading the number of gridpoints (IE, JE) out of the IGDB
!!    IEGB     =  IGDS( 5)
!!    JEGB     =  IGDS( 6)
 
      !     Getting the number of layers (KE)
      !     Only for data on model levels, not for pressure level data
!!    NLEVTYP  = IPDS( 8)
!!    IF (NLEVTYP.NE.100 .OR. NLEVTYP.NE.102) THEN
!!       KEGB  = (IGDS(1) - 42)/4 - 4
!!    ELSE
!!       KEGB  = ke_grib
!!    ENDIF

      !     Reading the element, the level typ, level, date and time
      !     out of the Product-Definition-Block (IPDS)
!!    NELEMENT = IPDS( 7)
!! 
!!    NLEVTYP  = IPDS( 8)
!!    IF (NLEVTYP.NE.100 .OR. NLEVTYP.NE.102) THEN
!!       NLEVEL= IPDS( 9)
!!    ELSE
!!       NLEVEL= IPDS(10)
!!    ENDIF
   
!!    IF (NLEVTYP.EQ.109) THEN
!!       NLEVEL= IPDS(10)
!!    ENDIF

!!    NYEAR    = IPDS(11)
!!    NMONTH   = IPDS(12)
!!    NDAY     = IPDS(13)
!!    NTIME    = IPDS(14)
!!    NTFLAG   = IPDS(19)
!!    IF ( NTFLAG.EQ.10 ) THEN
!!       NVV    = IPDS(18)
!!    ELSE IF ( NTFLAG.EQ.3 .OR. NTFLAG.EQ.4 ) THEN
!!       NVV1   = IPDS(17)
!!       NVV2   = IPDS(18)
!!    ENDIF
   
!! should ivv IN/OUT instead of OUT ?????????
      IF (PRESENT(ivv)) THEN
         IF (ivv<endStep) THEN
             ivv=endStep    
         ENDIF
      ENDIF

!..............check variable against griblist, is variable required?.............................
      DO nv=1,nvarmax
         IF (ASSOCIATED(griblist(nv)%p2).OR.ASSOCIATED(griblist(nv)%p3)) THEN
            found=.TRUE.
            IF (TRIM(griblist(nv)%sname) /= TRIM(sname)) found=.FALSE.
            IF (TRIM(griblist(nv)%levtype) /= TRIM(levtype)) found=.FALSE.
            IF (griblist(nv)%toplev /= -1 .AND. griblist(nv)%toplev /= toplev) found=.FALSE.

            IF (ASSOCIATED(griblist(nv)%p2)) THEN
              IF (griblist(nv)%botlev /= -1 .AND. griblist(nv)%botlev /= botlev) found=.FALSE.
              IF (TRIM(griblist(nv)%stattype) /= 'unknown' .AND. TRIM(griblist(nv)%stattype) /= TRIM(stype)) found=.FALSE.
            ENDIF
           
            IF (found) THEN
               lvarfound(nv)=.TRUE.

            ! Check if radiation and precipitation are hourly data (lhourly_data=.true.)
            ! --------------------------------------------------------------------------
              IF (stype == 'accum' .OR. stype == 'avg') THEN
                IF (endstep /= startstep) THEN
                  nno = nno + 1
                  ntime_int(nno) = ABS (endstep - startstep)
                  WRITE (6,*) " Time interval for ",sname," is ",ntime_int(nno)
                ENDIF
              ENDIF

               CALL grib_get_size (igrib, 'values', lds)
!!DL ie/je_grib are ie/je for COSMO and nblock/nproma for ICON
               IF (lds > ie_grib*je_grib) THEN
                 WRITE(6,*) "read_grib_eccodes: ERROR: size of message is too big for allocated field: ", &
                             lds, ie_grib*je_grib
                 STOP
               ENDIF
!!DL because of icon block structure preset field elements with rundef (like in reading routines)
!!DL           ds = 0.0_irealgrib
               ds = rundef
               ds(1:lds) = 0.0_irealgrib
               CALL grib_set (igrib, 'missingValue', undefgrib)
               CALL grib_get (igrib, 'values', ds, ierr)
               IF (ierr /= 0) THEN
                  WRITE(6,*) 'read_grib_eccodes: Error in grib_get for values!'
                  STOP
               ENDIF

!!DL           IF ( ANY(ds(:ie_grib*je_grib) == undefgrib).AND.PRESENT(anyundef) ) THEN
!!DL - to be used with ICON
               IF ( ANY(ds(1:lds) == undefgrib).AND.PRESENT(anyundef) ) THEN
                  anyundef = .TRUE.
               ENDIF
               IF (PRESENT(i02)) THEN
                  ALLOCATE (dstemp(ie_grib,je_grib))
                  dstemp=RESHAPE(ds(1:ie_grib*je_grib),(/ ie_grib, je_grib /))
                  IF (ASSOCIATED(griblist(nv)%p2)) THEN
                     griblist(nv)%p2(:,:)=dstemp(i0:i0+di-1,j0:j0+dj-1)
                  ELSE
                     griblist(nv)%p3(:,:,botlev)=dstemp(i0:i0+di-1,j0:j0+dj-1)
                  ENDIF
                  DEALLOCATE(dstemp)
               ELSE
                  IF (ASSOCIATED(griblist(nv)%p2)) THEN
                     griblist(nv)%p2(:,:)=RESHAPE(ds(1:ie_grib*je_grib),(/ ie_grib, je_grib /))
                  ELSE
PRINT * ,"IN read_grib_eccodes: ie_grib=",ie_grib," je_grib=",je_grib
                     griblist(nv)%p3(:,:,botlev)=RESHAPE(ds(1:ie_grib*je_grib),(/ ie_grib, je_grib /))
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO ! nv

     call grib_release(igrib)

   ENDDO !read_loop

! End of checking hourly data
   IF (nno > 0) THEN
     IF (ALL (ntime_int(1:nno) == 1)) THEN
       WRITE (6,*) " lhourly_data should be true!"
       IF(PRESENT(lhourly)) lhourly = .TRUE.
     ELSE
       WRITE (6,*) " lhourly_data should be false!"
     ENDIF
   ENDIF


!   999 CONTINUE  XYZ> cancellation of GOTO 999

   IF (COUNT(lvarfound)<nvarmax) THEN
      imiss_soil=0
      imiss_tg=0
      imiss_tsnow=0
      imiss_wsnow=0
      imiss_ts=0
      WRITE(6,*) "WARNING: The following fields of griblist are not available:"
       DO nv=1,nvarmax
          IF (.NOT.lvarfound(nv)) THEN
               WRITE(6,*) "!!! shortName:",griblist(nv)%sname,"lty:",griblist(nv)%levtype,"levtop, levbot:",griblist(nv)%toplev,griblist(nv)%botlev,&
                          "stepType:",griblist(nv)%stattype
!DL NEW : Do some checks for initial data
             IF ((griblist(nv)%sname=='T_SO') .OR. (griblist(nv)%sname=='W_SO')) imiss_soil=imiss_soil+1
             IF ( griblist(nv)%sname=='T_G'                                    ) imiss_tg=1
             IF ( griblist(nv)%sname=='T_S'                                    ) imiss_ts=1
             IF ( griblist(nv)%sname=='T_SNOW'                                 ) imiss_tsnow=1
             IF ( griblist(nv)%sname=='W_SNOW'                                 ) imiss_wsnow=1
          ENDIF
       ENDDO
               WRITE(6,*) "!!! imiss_tg, imiss_ts, imiss_tsnow, imiss_wsnow, imiss_soil: ", imiss_tg, imiss_ts, imiss_tsnow, imiss_wsnow, imiss_soil 
       IF (imiss_soil >= 1) THEN
          WRITE(6,*) "STOP: ",imiss_soil, " soil variables are missing!!, STOP"
          STOP
       ELSEIF ((imiss_tg == 1) .AND. ((imiss_ts+imiss_tsnow+imiss_wsnow) >= 1)) THEN
          WRITE(6,*) "STOP: T_G can not be computed because of missing T_SNOW or W_SNOW or T_S or T_SO(0) !!, STOP"
          STOP
       ENDIF
   ENDIF

   CALL grib_close_file (nudat,ierr)
   IF (ierr /= 0) THEN
      WRITE(6,*) 'read_grib_eccodes: Error in grib_close_file!'
      STOP
   ENDIF

   DEALLOCATE (lvarfound)
   DEALLOCATE (ds)

END SUBROUTINE read_grib_eccodes

!==============================================================================
SUBROUTINE get_gribinfo(filename,ie,je,rlat0,rlon0,dlat,dlon,pollat,pollon,isfound)
!----------------------------------------------------------------------
! Description:
!
! Extracts GRIB info for later GRIB read
!----------------------------------------------------------------------

  IMPLICIT NONE

  ! Subroutine Arguments:
  CHARACTER (LEN=*),        INTENT(IN) :: filename
  INTEGER                              :: &
       ie                               , & !
       je                                   ! 
  REAL (KIND=wp),          INTENT(OUT) :: &
       rlat0                            , & ! 
       rlon0                            , & ! 
       dlat                             , & ! 
       dlon                             , & ! 
       pollat                           , & ! 
       pollon                               ! 
  LOGICAL,                 INTENT(OUT) :: isfound

  ! Local variables
  INTEGER                              :: &
       i                                , & ! 
       n                                , & ! 
       b                                , & ! 
       pos                              , & ! 
       griblen                              ! 
  CHARACTER (LEN=1)                    :: c
  CHARACTER (LEN=4), PARAMETER         :: grib="GRIB"
  LOGICAL                              :: ex


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
       iednr  =  1                          ! grib edition number 
  INTEGER (KIND=intgribf)              :: &
       lfd                              , & ! Dimension for iblock
       lbm                              , & ! Dimension for bitmap: this value has to be at
                                            !   least the same as lbmax in SUBROUTINE grbin1
       lds                                  ! Dimension for unpacked data

  INTEGER (KIND=intgribf)              :: &
       idims_in  (ndims)                , & ! array for all dimensions (input)
       ipds      (npds)                 , & ! product definition section
       igds_in   (ngds)                 , & ! grid description section for output
       ibms      (nbms)                 , & ! bit map section
       ibds      (nbds)                     ! binary data section


  ! Arrays that are allocated during program execution
  INTEGER (KIND=intgribf), ALLOCATABLE :: &
       iblock    (:)                    , & ! array for gribed data
       ibmap     (:)                        ! array for bit map

  REAL   (KIND=irealgrib), ALLOCATABLE :: &
       dsup (:)                         , & ! array for special data
       ds   (:)                             ! array for unpacked data


  INTEGER (KIND=intgribf)              :: &
       ilen                             , & ! 
       ierr                             , & ! 
       istat                                ! 
  INTEGER  (KIND=intgribc)             :: & ! corresponding variables for C-routines
       nudatc                           , & ! 
       maxlenc                          , & ! 
       ilenc                            , & ! 
       ierrc                                !     

! Body of subroutine:
  ! Determine length of first grib record
  isfound=.TRUE.
  pos=1
  i=1
  griblen=0
  INQUIRE(FILE=TRIM(filename),EXIST=ex)
  IF (.NOT.ex) THEN
     isfound=.FALSE.
     RETURN
  ENDIF

  OPEN(99,FILE=filename,ACCESS="direct",STATUS='old',RECL=1)

  DO WHILE ((pos<=4).AND.(i<30))
     READ (99,REC=i) c
     b=ICHAR(c)
     IF (grib(pos:pos)==c) THEN
        pos=pos+1
     ELSE
        pos=1
     ENDIF
     i=i+1
  ENDDO

  DO n=i,i+2
     READ (99,REC=n) c
     b=ICHAR(c)
     griblen=griblen*256+b
  ENDDO
  CLOSE(99)

  ! Read first grib record
  lfd=griblen/intgribf+100
  lds=lfd * iwlength / 2
  lbm = 1875

  ALLOCATE (iblock(lfd), ibmap(lbm), STAT=istat)
  ALLOCATE (ds(lds), dsup(ndsup), STAT=istat)

  ! dimensions for grib routines
  idims_in( 1) = npds
  idims_in( 2) = ngds
  idims_in( 3) = nbms
  idims_in( 4) = nbds
  idims_in( 5) = lbm
  idims_in( 6) = ndsup
  idims_in( 7) = -1             ! Do not unpack data set
  idims_in( 8) = lfd
  idims_in(9:20) = 0

  nudatc=250
  CALL copen (nudatc,TRIM(filename),'r  ',ierrc)

  maxlenc=lfd*iwlength
  CALL cuegin(nudatc,maxlenc, iblock, ilenc, ierrc)
  ilen=ilenc
  IF (ilen == 0) THEN
     WRITE(6,*) "Error in get_gribinfo!"
     STOP
  ENDIF
  ds = 0.0_irealgrib
  CALL grbin1(iednr, undefgrib, ndims, idims_in, iblock, ibmap,         &
       ipds, igds_in, ibms, ibds, dsup, ds, ierr)

  ! Interpretation of GDS
  rlat0=igds_in(7)/1000.0_wp
  rlon0=igds_in(8)/1000.0_wp
  ie=igds_in(5)
  je=igds_in(6)
  dlat=(igds_in(10)-igds_in(7))/1000.0_wp/FLOAT(je-1)
  dlon=(igds_in(11)-igds_in(8))/1000.0_wp/FLOAT(ie-1)

  pollat=-igds_in(20)/1000.0_wp
  pollon=igds_in(21)/1000.0_wp-180.0_wp
  IF (pollon<-180.0_wp) pollon=pollon+360.0_wp

  ! Clean up
  CALL cclose(nudatc,'exi',ierrc)
  DEALLOCATE (iblock, ibmap)
  DEALLOCATE (ds, dsup)

END SUBROUTINE get_gribinfo
!------------------------------------------------------------------------------

SUBROUTINE get_grib_info(filename,ie,je,rlat0,rlon0,dlat,dlon,pollat,pollon,isfound)
!----------------------------------------------------------------------
! Description:
!
! Extracts GRIB info for later GRIB read - using eccodes library from ECMWF
!----------------------------------------------------------------------
  USE grib_api

  IMPLICIT NONE

  ! Subroutine Arguments:
  CHARACTER (LEN=*),        INTENT(IN) :: filename
  INTEGER                              :: &
       ie                               , & !
       je                                   ! 
  REAL (KIND=wp),          INTENT(OUT) :: &
       rlat0                            , & ! 
       rlon0                            , & ! 
       dlat                             , & ! 
       dlon                             , & ! 
       pollat                           , & ! 
       pollon                               ! 
  LOGICAL,                 INTENT(OUT) :: isfound

  ! Local variables
  LOGICAL                              :: ex

  ! Grib variables
  REAL    (KIND=irealgrib), PARAMETER  :: &
       undefgrib =  -1.0E7_irealgrib        ! value for "undefined" in the grib routines
  REAL (KIND=wp)                       :: &
       rlatl, rlonl                         ! last latitude and longitude of grid 
  INTEGER                              :: &
       igrib                            , & ! grib handle
       iednr                                ! grib edition 

  INTEGER (KIND=intgribf)              :: &
       nudatc                           , & ! file unit
       ierr                                 ! error status
  CHARACTER (LEN=30)                    :: gridtype

! Body of subroutine:
! ******************
  ! Does the file exists?                     
  isfound=.TRUE.
  INQUIRE(FILE=TRIM(filename),EXIST=ex)
  IF (.NOT.ex) THEN
     isfound=.FALSE.
     RETURN
  ENDIF

  ! Read grid information from first grib record

  ierr=0

  CALL grib_open_file (nudatc,TRIM(filename),'r  ',ierr)
  CALL grib_new_from_file (nudatc, igrib, ierr)

  IF (ierr /= 0) THEN
     WRITE(6,*) "Error in get_grib_info with grib_open_file or grib_new_from_file!"
     STOP
  ENDIF

  CALL grib_get(igrib, 'edition'                           , iednr , ierr)
  CALL grib_get(igrib, 'gridType'                          , gridtype, ierr)
 
  IF (gridtype == 'rotated_ll') THEN
    CALL grib_get(igrib, 'latitudeOfFirstGridPointInDegrees' , rlat0 , ierr)  
    CALL grib_get(igrib, 'longitudeOfFirstGridPointInDegrees', rlon0 , ierr)
    CALL grib_get(igrib, 'Ni'                                , ie    , ierr)
    CALL grib_get(igrib, 'Nj'                                , je    , ierr)
    CALL grib_get(igrib, 'iDirectionIncrementInDegrees'      , dlat  , ierr)
    CALL grib_get(igrib, 'jDirectionIncrementInDegrees'      , dlon  , ierr)
    CALL grib_get(igrib, 'latitudeOfSouthernPoleInDegrees'   , pollat, ierr)
    CALL grib_get(igrib, 'longitudeOfSouthernPoleInDegrees'  , pollon, ierr)
      CALL grib_get(igrib, 'latitudeOfLastGridPointInDegrees'  , rlatl , ierr)
      CALL grib_get(igrib, 'longitudeOfLastGridPointInDegrees' , rlonl , ierr)
    IF (iednr == 1) THEN
      dlat=(rlatl - rlat0)/FLOAT(je-1)
      dlon=(rlonl - rlon0)/FLOAT(ie-1)
    ENDIF
  ELSE
    WRITE(6,*) "get_grib_info: Implementation only for rotated_ll, ICON grid missing!"
    STOP
  ENDIF

!!DL for GRIB2 input
  IF (rlon0 > 180.0_wp) rlon0 = rlon0 - 360.0_wp
!!DL for GRIB2 input
  pollat= -pollat
  pollon= pollon - 180.0_wp
  IF (pollon<-180.0_wp) pollon=pollon+360.0_wp

!!DL
  WRITE(6,*) "get_grib_info: rlat0,rlon0,ie,je,dlat,dlon,pollat,pollon,rlatl,rlonl", &
                             rlat0,rlon0,ie,je,dlat,dlon,pollat,pollon,rlatl,rlonl
!!DL

  IF (ierr /= 0) THEN
     WRITE(6,*) "Error in get_grib_info with grib_get!"
     STOP
  ENDIF

  ! Clean up
  CALL grib_close_file (nudatc)

END SUBROUTINE get_grib_info

!------------------------------------------------------------------------------

SUBROUTINE get_grib_info_icon(filename,ngrid,ngpi,ngp,uuid_hor_grid,isfound)
!----------------------------------------------------------------------
! Description:
!
! Extracts GRIB info of ICON fields for later GRIB read - using eccodes library from ECMWF
!----------------------------------------------------------------------
  USE grib_api

  IMPLICIT NONE

  ! Subroutine Arguments:
  CHARACTER (LEN=*),        INTENT(IN) :: filename
  INTEGER                 , INTENT(IN) :: &
       ngrid                            , & ! numberOfGridUsed - number of horizontal grid
       ngpi                                 ! generatingProcessIdentifier - to identify domain
  INTEGER                 ,INTENT(OUT) :: &
       ngp                                  ! number of grid points = numberOfValues
  CHARACTER (LEN=1),       INTENT(OUT) :: &
       uuid_hor_grid(16)                        ! identifier for horizontal grid for later comparision
  LOGICAL,                 INTENT(OUT) :: isfound

  ! Local variables
  LOGICAL                              :: ex

  ! Grib variables
!!  REAL    (KIND=irealgrib), PARAMETER  :: &
!!       undefgrib =  -1.0E7_irealgrib        ! value for "undefined" in the grib routines
  INTEGER                              :: &
       igrib                                ! grib handle

  INTEGER (KIND=intgribf)              :: &
       nudatc                           , & ! file unit
       ierr                                 ! error status

  INTEGER                              :: &
       ngrid_in, ngpi_in                    ! 

  CHARACTER (LEN=30)                    :: gridtype

! Body of subroutine:
! ******************
  ! Does the file exists?                     
  isfound=.TRUE.
  INQUIRE(FILE=TRIM(filename),EXIST=ex)
  IF (.NOT.ex) THEN
     isfound=.FALSE.
     RETURN
  ENDIF

  ! Read grid information from first grib record

  ierr=0

  CALL grib_open_file (nudatc,TRIM(filename),'r  ',ierr)
  CALL grib_new_from_file (nudatc, igrib, ierr)

  IF (ierr /= 0) THEN
     WRITE(6,*) "Error in get_grib_info_icon with grib_open_file or grib_new_from_file!"
     STOP
  ENDIF

  CALL grib_get(igrib, 'gridType'                     , gridtype      , ierr)
 
  IF (gridtype == 'unstructured_grid') THEN
    CALL grib_get(igrib, 'numberOfGridUsed'           , ngrid_in      , ierr)  
    CALL grib_get(igrib, 'generatingProcessIdentifier', ngpi_in       , ierr)
    CALL grib_get(igrib, 'uuidOfHGrid'                , uuid_hor_grid , ierr)
    CALL grib_get(igrib, 'numberOfDataPoints'         , ngp           , ierr)
  ELSE
    WRITE(6,*) "get_grib_info_icon: Implementation for ICON grid!"
    STOP
  ENDIF

  IF (ierr /= 0) THEN
     WRITE(6,*) "Error in get_grib_info_icon with grib_get!"
     STOP
  ENDIF

  ! Check values with NL input
  IF (ngrid /= ngrid_in) THEN
    WRITE (6,*) "get_grid_info_icon: ngrid .ne. ngrid_in: ", ngrid, ngrid_in
    STOP 'WRONG ICON grid!'
  ENDIF
  IF (ngpi_in /= ngpi) THEN
    WRITE (6,*) "get_grid_info_icon: ngpi .ne. ngpi_in: ", ngpi, ngpi_in
    STOP 'WRONG ICON domain!'
  ENDIF 

  ! Clean up
  CALL grib_close_file (nudatc)

END SUBROUTINE get_grib_info_icon

!------------------------------------------------------------------------------

!==============================================================================

END MODULE tsa_gribio
