!+ module for interpolation of input data, 
!------------------------------------------------------------------------------

MODULE tsa_interpol

!------------------------------------------------------------------------------
!
! Description:
!  This module contains subroutines for interpolation of input data 
!  according to user selected interpolation method.
!  Details of methods inside.
!
!    Currently included:
!    - general_interpol2d:
!      general interpolation of input data
!    - polyfit_solve (Function):
!      calculates linear or polynomial fit
!    - polyfit_fill:
!      ---
!    - gaussj:
!      calculates the gaussian for polyfit_solve
!
!
! Current Code Owner: DWD, ???
!  phone:  +49  69  
!  fax:    +49  69  
!  email:  
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
! V5.01      2015/12/01 Yiftach Ziv, IMS (XYZ)
!  Arranged code to adhere to coding standards.
! V5.07      2020/02/21     Ulrich Schaettler
!  Replaced data_parameters by kind_parameters and removed iintegers
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
USE kind_parameters , ONLY :   &
  wp           ! KIND-type parameters for real variables

!=======================================================================

IMPLICIT NONE

!=======================================================================
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
! Global (i.e. public) Declarations:

INTEGER,  PARAMETER :: nhistmax=50

REAL (KIND=wp), PRIVATE,   PARAMETER :: rundef=-9999.0

CONTAINS
!==============================================================================
!+ Subroutine for general interpolation of input data
!------------------------------------------------------------------------------

SUBROUTINE general_interpol2d(rlat1,rlon1,field1,                   &
                              startlat_tot,startlon_tot,dlat,dlon,  &
                              field2,type_neigh,type_int,rinf,rexp, &
                              lmask1,lmask2)


!------------------------------------------------------------------------------
!
! Description:
!   This routine performs interpolation of input data according 
!   to selected type arguments
!
! Method:
!   Arithmetical statements
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  REAL (KIND=wp), INTENT(IN)            ::   &
       field1(:,:)                         , & ! Values of the output field
       rlat1(:,:)                          , & ! Coordinates of the output field
       rlon1(:,:)                              ! ---       "       ---
  REAL (KIND=wp), INTENT(OUT)           ::   &
       field2(:,:)                             ! Interpolation Field 
  INTEGER,        INTENT(IN)            ::   &
       type_neigh                          , & ! Neighborhood Type
       type_int                                ! Interpolation Type
  REAL (KIND=wp), INTENT(IN)            ::   &
       startlat_tot                        , & ! Definition of output grid
       startlon_tot                        , & ! ---       "       ---
       dlat                                , & ! ---       "       ---
       dlon                                    ! ---       "       ---
  REAL (KIND=wp), INTENT(IN)            ::   &
       rinf                                , & ! Radius of influence 
                                               ! (IN gridpoints des output grid)
       rexp                                    ! parameter for exponential weighting

  LOGICAL, OPTIONAL ,INTENT(IN)         ::   &
       lmask1(:,:)                         , & ! Masking of points which are
       lmask2(:,:)                             ! not used for interpolation

!!! == Type Parameter Explanation: =======================================
!!! type_int: 
!!! 1 = next neighbour (irregular, 2d only)
!!! 2 = average (irregular, 2d only)
!!! 3 = exponential weighting  (rexp, irregular, 2d only)          
!!! 4 = dominant type  (integer values only, irregular, 2d only)
!!! 5 = logarhitmic average (irregular, 2d only)
!!! 6 = crop (fill part of total field) (irregular, 2d only)
!!! 7 = polynomial fitting (linear)
!!! 8 = polynomial fitting (quadratic)
!!!
!!! type_neigh:
!!! 1 = NEAREST: each original grid point is assigned 
!!!              only to the nearest interpolation point
!!! 2 = SQUARE: all point within a square of size rinf are used
!!! 3 = CIRCLE: all point within a sphere of radius rinf are used
!!!======================================================================


! Local variable
! --------------------
  INTEGER                                ::  &
     i,j,imin,imax,jmin,jmax,iin,jin       , & !
     ie1,je1,ie2,je2                       , & !
     a(1),nval,ncount                      , & !
     norder,nsys
  REAL (KIND=wp)                         ::  &
       ri,rj,dist2,rinf2,weight                !

  REAL (KIND=wp), ALLOCATABLE            ::  &
        rsum(:,:)                          , & !
        rcount(:,:)                        , & !
        ramat(:,:,:,:)                     , & !
        rbvec(:,:,:)                           !
  INTEGER, ALLOCATABLE                   ::  &
        rhist(:,:,:)                           !
  LOGICAL                                ::  &
        lok                                    !

! Body of Subroutine
! --------------------

! Initialization

  ie1=size(field1,1)
  je1=size(field1,2)
  ie2=size(field2,1)
  je2=size(field2,2)


  ALLOCATE(rcount(ie2,je2))
  rcount=0
  SELECT CASE (type_int)
    CASE (1,2,3,5) 
      ALLOCATE(rsum(ie2,je2))
      rsum=0.0
      IF (type_int==1) THEN
         rcount=-1.0
      ENDIF
    CASE (4,9)
      ALLOCATE(rhist(ie2,je2,nhistmax))
      rhist=0.0
    CASE (7,8)
      IF (type_int==7) THEN
         norder=1
      ELSE
         norder=2
      ENDIF
      nsys=NINT((norder+1)*(norder+2)*0.5)
      ALLOCATE(ramat(nsys,nsys,ie2,je2), &
               rbvec(nsys,ie2,je2))
      ramat=0.0
      rbvec=0.0
  END SELECT

  field2=rundef

! Calculations

  IF (type_int/=6) THEN

      rinf2=rinf*rinf
      ncount=0
      DO iin=1,ie1
         DO jin=1,je1
            
            lok=(field1(iin,jin)/=rundef)
            IF (PRESENT(lmask1)) THEN
               lok=lok.AND.lmask1(iin,jin)
            ENDIF

            IF (lok) THEN
            ri=(rlon1(iin,jin)-startlon_tot)/dlon+1.0
            rj=(rlat1(iin,jin)-startlat_tot)/dlat+1.0
            
               IF (type_neigh==1) THEN
                  i=NINT(ri)
                  j=NINT(rj)
                  IF ((i>=1).AND.(j>=1).AND.(i<=ie2).AND.(j<=je2)) THEN
                     imin=i;imax=i
                     jmin=j;jmax=j
                  ELSE
                     imin=0;imax=-1
                     jmin=0;jmax=-1
                  ENDIF
               ELSE
                  imin=MAX(1,FLOOR(ri-rinf)+1)
                  imax=MIN(ie2,FLOOR(ri+rinf))
                  jmin=MAX(1,FLOOR(rj-rinf)+1)
                  jmax=MIN(je2,FLOOR(rj+rinf))
               ENDIF
            
!              write(6,*) iin,jin,imin,imax,jmin,jmax
               DO i=imin,imax
                  DO j=jmin,jmax

                     IF (PRESENT(lmask2)) THEN
                        IF (.NOT.lmask2(i,j)) CYCLE
                     ENDIF
                     IF ((type_neigh>1).OR.(type_int==1).OR.(type_int==3)) THEN
                        dist2=(ri-REAL(i))**2+(rj-REAL(j))**2
                     ELSE
                        dist2=1.0
                     ENDIF

                     IF ((type_neigh/=3).OR.(dist2<rinf2)) THEN
                        SELECT CASE (type_int)
                           CASE (1)
                              IF ((rcount(i,j)<0).OR.(rcount(i,j)>dist2)) THEN
                                 rsum(i,j)=field1(iin,jin)
                                 rcount(i,j)=dist2
                              ENDIF
                           CASE (2)
                              rsum(i,j)  = rsum(i,j)+field1(iin,jin)
                              rcount(i,j)= rcount(i,j)+1
                           CASE (3)
                              weight=exp(-sqrt(dist2)/rexp)
                              rsum(i,j)  = rsum(i,j)+&
                                           weight*field1(iin,jin)
                              rcount(i,j)= rcount(i,j)+weight
                           CASE (4,9)
                              nval=NINT(field1(iin,jin))
                              IF ((nval>0).AND.(nval<=nhistmax)) THEN
                                 rhist(i,j,nval)=rhist(i,j,nval)+1
                                 rcount(i,j)= rcount(i,j)+1
                              ENDIF
                           CASE (5)
                              IF (field1(iin,jin)<1e-4) THEN 
                                 rsum(i,j)  = rsum(i,j) + log(1e-4)
                              ELSE
                                 rsum(i,j)  = rsum(i,j)+log(field1(iin,jin))
                              ENDIF
                              rcount(i,j)= rcount(i,j)+1
                           CASE (7,8)
                              CALL polyfit_fill(ramat(:,:,i,j),rbvec(:,i,j), &
                                                ri-i,rj-j,field1(iin,jin),norder)
                              rcount(i,j)= rcount(i,j)+1
                        END SELECT
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO
   
!!! Must be revised !!!
!   ELSE ! type_int==6 (partly crop)

!!$      IF ((ds(n)%iin_crop<1).OR.(ds(n)%jin_crop<1)) THEN
!!$         IF ((nlam==0).OR.(nphi==0)) THEN
!!$            write(6,*) "Coordinates missing! partly crop 1"
!!$            stop
!!$         ENDIF
!!$         IF (((.NOT.associated(ds(n)%field(nlam)%r)).OR. &
!!$              (.NOT.associated(ds(n)%field(nphi)%r)))) THEN
!!$            write(6,*) "Coordinates missing! partly crop 2"
!!$            stop
!!$         ENDIF
!!$         ds(n)%iin_crop=NINT((ds(n)%field(nlam)%r(1,1,1)-startlon_tot)/dlon)+1.0
!!$         ds(n)%jin_crop=NINT((ds(n)%field(nphi)%r(1,1,1)-startlat_tot)/dlat)+1.0        
!!$         write(6,*) "New iin_crop jin_crop: ", ds(n)%iin_crop,ds(n)%jin_crop
!!$         ds(n)%ie_crop=ds(n)%ie; ds(n)%je_crop=ds(n)%je
!!$      ENDIF
!!$      IF ((ds(n)%i0_crop<1).OR.(ds(n)%j0_crop<1)) THEN
!!$         ifac=ie_tot/ds(n)%ie
!!$         jfac=je_tot/ds(n)%je
!!$         IF (type_neigh/=6) THEN
!!$            DO i=1,ds(n)%ie
!!$               imin=int((i-1)*ifac)+1
!!$               imax=int(i*ifac)
!!$               DO j=1,ds(n)%je
!!$                  jmin=int((j-1)*jfac)+1
!!$                  jmax=int(j*jfac)
!!$                  rfield(imin:imax,jmin:jmax) = ds(n)%field(nvar)%r(i,j,klev)
!!$               ENDDO
!!$            ENDDO
!!$         ELSE
!!$            DO i=1,ds(n)%ie
!!$               imin=int((i-1)*ifac)+1
!!$               imax=int(i*ifac)
!!$               DO j=1,ds(n)%je
!!$                  jmin=int((j-1)*jfac)+1
!!$                  jmax=int(j*jfac)
!!$                  DO i2=imin,imax
!!$                     DO j2=jmin,jmax
!!$                        IF ((ds(n)%field(nlAND)%r(i,j,1)<=0.5).eqv. &
!!$                             (fr_lAND(i2,j2)<=0.5)) THEN
!!$                           rfield(i2,j2) = ds(n)%field(nvar)%r(i,j,klev)
!!$                        ENDIF
!!$                     ENDDO
!!$                  ENDDO
!!$               ENDDO
!!$            ENDDO
!!$         ENDIF
!!$       ELSE
!!$         rfield(ds(n)%iin_crop:ds(n)%iin_crop+ds(n)%ie_crop-1, &
!!$              ds(n)%jin_crop:ds(n)%jin_crop+ds(n)%je_crop-1)= &
!!$              ds(n)%field(nvar)%r(ds(n)%i0_crop:ds(n)%i0_crop+ds(n)%ie_crop-1, &
!!$              ds(n)%j0_crop:ds(n)%j0_crop+ds(n)%je_crop-1,klev)
!!$      ENDIF
!!$      fntile=REAL(ntile)
!!$      protect(NINT((ds(n)%iin_crop-1)/fntile+1): &
!!$              NINT((ds(n)%iin_crop-1)/fntile+1+ds(n)%ie_crop/fntile)-1, &
!!$              NINT((ds(n)%jin_crop-1)/fntile+1): &
!!$              NINT((ds(n)%jin_crop-1)/fntile+1+ds(n)%je_crop/fntile)-1)=n
!!$   ENDIF

  ENDIF
   

  IF (type_int/=6)  THEN
      DO i=1,ie2
         DO j=1,je2
            SELECT CASE (type_int)
               CASE (1)
                  IF (rcount(i,j)>0) THEN
                     field2(i,j)=rsum(i,j)
                  ELSE
                     field2(i,j)=rundef
                  ENDIF
               CASE (2,3)
                  IF ( rcount(i,j)>0) THEN
                     field2(i,j)  = rsum(i,j)/ rcount(i,j)
                  ELSE
                     field2(i,j)=rundef
                  ENDIF
               CASE (4)
                  IF (SUM(rhist(i,j,:))>0) THEN
                     a=MAXLOC(rhist(i,j,:))
                     field2(i,j)=a(1)
                  ELSE
                     field2(i,j)=rundef
                  ENDIF
               CASE (5)
                  IF ( rcount(i,j)>0) THEN
                     field2(i,j)  = EXP(rsum(i,j)/ rcount(i,j))
                  ELSE
                     field2(i,j)=rundef
                  ENDIF
               CASE (7,8)
                  SELECT CASE (type_int)
                     CASE (7)
                        norder=1
                     CASE (8)
                        norder=2
                  END SELECT

                  IF (rcount(i,j)<(norder+1)*(norder+2)*0.5) THEN
                     IF (rcount(i,j)>0) THEN
                        field2(i,j)=rbvec(1,i,j)/rcount(i,j)
                     ELSE
                        field2(i,j)=rundef
                     ENDIF
                  ELSE
                     field2(i,j)=polyfit_solve(ramat(:,:,i,j),rbvec(:,i,j), &
                                 NINT(rcount(i,j)),norder,.TRUE.)
                  ENDIF
            END SELECT
         ENDDO
      ENDDO
   ENDIF

! Finalizing
! --------------------

   SELECT CASE (type_int)
      CASE (1,2,3,5)
         DEALLOCATE(rsum)
      CASE (4)
         DEALLOCATE(rhist)
      CASE (7,8)
         DEALLOCATE(ramat,rbvec)
   END SELECT
   DEALLOCATE(rcount)

!   write(6,*) "Min. and Max.: ",minval(field2),maxval(field2)

END SUBROUTINE general_interpol2d


!------------------------------------------------------------------------------

REAL FUNCTION polyfit_solve(a,b,npoint,norder,lmean)

!------------------------------------------------------------------------------
!
! Description:
!   This function calculates linear or polynomial fit and returns 
!   the fit
!
!------------------------------------------------------------------------------

! Function arguments:
! --------------------

  INTEGER, INTENT(IN)                  ::  &
          norder,npoint
  REAL    (KIND=wp), INTENT(INOUT)     ::  &
          a(:,:),b(:)
  LOGICAL, INTENT(IN)                  ::  &
          lmean

! Local variable
! --------------------
  INTEGER                              ::  &
          nsys,istat
  REAL    (KIND=wp)                    ::  &
          mean

! Body of Function
! --------------------

  nsys=NINT((norder+1)*(norder+2)*0.5)

  a=a/REAL(npoint)
  b=b/REAL(npoint)
  mean=b(1)

  CALL gaussj(a,nsys,nsys,b,1,1,istat)

  IF (istat==0) THEN
     polyfit_solve=b(1)
  ELSE
     IF (lmean) THEN
        polyfit_solve=mean
     ELSE
        polyfit_solve=rundef
     ENDIF
  ENDIF

END FUNCTION polyfit_solve


!------------------------------------------------------------------------------

SUBROUTINE polyfit_fill(a,b,xkor,ykor,val,norder)

!------------------------------------------------------------------------------
!
! Description:
!   
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  REAL    (KIND=wp)       , INTENT(IN)    ::  &
          xkor,ykor
  REAL    (KIND=wp)       , INTENT(IN)    ::  &
          val
  INTEGER                 , INTENT(IN)    ::  &
          norder
  REAL    (KIND=wp)       , INTENT(INOUT) ::  &
          a(:,:),b(:)

! Local variable
! --------------------

  INTEGER                                 ::  &
          nsx,nsy,nox,noy,nxx,nxy,nyx,nyy,i
  REAL    (KIND=wp)                       ::  &
          temp

! Body of Function
! --------------------

   DO nox=0,norder
      DO nxx=0,nox
         nyx=nox-nxx
         nsx=NINT((nox+1)*nox*0.5)+nxx+1
            
         temp=val
         DO i=1,nxx
            temp = temp*xkor
         ENDDO
         DO i=1,nyx
            temp = temp*ykor
         ENDDO
         b(nsx)=b(nsx)+temp
      
         DO noy=0,nox
            DO nxy=0,noy
               nyy=noy-nxy
               nsy=NINT((noy+1)*noy*0.5)+nxy+1

               IF (nsy<=nsx) THEN
                  temp=1.0
                  DO i=1,nxx+nxy
                     temp = temp*xkor
                  ENDDO
                  DO i=1,nyx+nyy
                     temp = temp*ykor
                  ENDDO
                  a(nsx,nsy)=a(nsx,nsy)+temp
                  a(nsy,nsx)=a(nsx,nsy)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO

END SUBROUTINE polyfit_fill


!------------------------------------------------------------------------------

SUBROUTINE gaussj(a,n,np,b,m,mp,istat)

!------------------------------------------------------------------------------
!
! Description:
!   The routine calculates the gaussian for polyfit_solve
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER                 , INTENT(IN)    ::  &
          n,np,m,mp
  INTEGER                 , INTENT(INOUT) ::  &
          istat
  REAL    (KIND=wp)       , INTENT(INOUT) ::  &
          a(np,np),b(np,mp)
  INTEGER                 , PARAMETER     :: NMAX=50


! Local variable
! --------------------

  INTEGER                                 ::  &
          i,icol,irow,j,k,l,ll             ,  &
          indxc(NMAX),indxr(NMAX),ipiv(NMAX)
  REAL    (KIND=wp)                       ::  &
          big,dum,pivinv
 
! Body of Function
! --------------------
 
  istat=0
  DO j=1,n
     ipiv(j)=0
  ENDDO
  DO i=1,n
     big=0.
     DO j=1,n
        IF (ipiv(j).NE.1) THEN
           DO k=1,n
              IF (ipiv(k).EQ.0) THEN
                 IF (ABS(a(j,k)).GE.big) THEN
                    big=ABS(a(j,k))
                    irow=j
                    icol=k
                 ENDIF
              ELSE IF (ipiv(k).GT.1) THEN
!             write(6,*) 'singular matrix IN gaussj 1'
              istat=1
              EXIT
              ENDIF
           ENDDO
        ENDIF
     ENDDO
     ipiv(icol)=ipiv(icol)+1
     IF (irow.NE.icol) THEN
        DO l=1,n
           dum=a(irow,l)
           a(irow,l)=a(icol,l)
           a(icol,l)=dum
        ENDDO
        DO l=1,m
           dum=b(irow,l)
           b(irow,l)=b(icol,l)
           b(icol,l)=dum
        ENDDO
     ENDIF
     indxr(i)=irow
     indxc(i)=icol
!FEA
!    IF (a(icol,icol)).EQ.0.) THEN
     IF (ABS(a(icol,icol)).LT.1e-6) THEN
!FEE
!        write(6,*) 'singular matrix IN gaussj 2'
         istat=2
         EXIT
     ENDIF
     pivinv=1./a(icol,icol)
     a(icol,icol)=1.
     DO l=1,n
        a(icol,l)=a(icol,l)*pivinv
     ENDDO
     DO l=1,m
        b(icol,l)=b(icol,l)*pivinv
     ENDDO
     DO ll=1,n
        IF (ll.NE.icol) THEN
           dum=a(ll,icol)
           a(ll,icol)=0.
           DO l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
           ENDDO
           DO l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
           ENDDO
        ENDIF
     ENDDO
  ENDDO
     DO l=n,1,-1
        IF (indxr(l).NE.indxc(l)) THEN
           DO k=1,n
              dum=a(k,indxr(l))
              a(k,indxr(l))=a(k,indxc(l))
              a(k,indxc(l))=dum
           ENDDO
        ENDIF
     ENDDO
     RETURN
END SUBROUTINE gaussj

END MODULE tsa_interpol
