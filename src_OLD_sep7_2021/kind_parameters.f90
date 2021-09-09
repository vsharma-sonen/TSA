! Module to determine kinds for different precisions
!------------------------------------------------------------------------------

MODULE kind_parameters

!==============================================================================
!
! Description for TSA version
!  This is kind_parameters from COSMO version 5.07 without modifications
!
!==============================================================================
!>
!!  Module determines kinds for different precisions.
!!  Number model from which the SELECTED_*\\_KIND are requested: <br>
!! @f{tabular}{{r@{\hspace*{3em}}c@{\hspace*{3em}}c}
!!                     &4 byte REAL     &8 byte REAL        \\\
!!        IEEE:        &precision = 6   &precision =   15   \\\
!!                     &exponent  = 37  &exponent  =  307
!! @f}
!! \\medskip
!!
!!
!! Current Code Owner
!! @author  Ulrich Schaettler, DWD
!!  email:  ulrich.schaettler@dwd.de
!!
!! @par Revision History
!
! History:
! Version      Date       Name
! ----------   ---------- ----
! V5_1         2014-11-28 Ulrich Schaettler
!  Initial Release
! V5_7         2020-02-21 Stefan Ruedisuehli
!  Introduced kind parameter vpp for variable precision physics
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

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: wp, sp, dp, vp, vpp, wi, i1, i2, i4, i8

!==============================================================================

  ! Floating point section
  ! ----------------------

  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND( 6, 37) !< single precision
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307) !< double precision

#ifdef SINGLEPRECISION
  INTEGER, PARAMETER :: wp = sp          !< working precision is single precision

  ! define variable physics precision, which might be different
#ifdef VAR_PHYSICS_PREC
  INTEGER, PARAMETER :: vpp = dp         !< variable precision physics is double precision
#else 
  INTEGER, PARAMETER :: vpp = sp         !< variable precision physics is working precision
#endif

#else
  ! not SINGLEPRECISION
  INTEGER, PARAMETER :: wp  = dp         !< working precision is double precision
  INTEGER, PARAMETER :: vpp = dp         !< variable precision physics is double precision
#endif

#ifdef __MIXED_PRECISION
  INTEGER, PARAMETER :: vp = sp          !< variable precision is single precision
#else
  INTEGER, PARAMETER :: vp = wp          !< variable precision is working precision
#endif


  ! Integer section
  ! ---------------

  INTEGER, PARAMETER :: i1 = SELECTED_INT_KIND(  2)   !< at least 1 byte integer
  INTEGER, PARAMETER :: i2 = SELECTED_INT_KIND(  4)   !< at least 2 byte integer
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(  9)   !< at least 4 byte integer
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND( 18)   !< at least 8 byte integer

  INTEGER, PARAMETER :: wi = i4                       !< selected working precision

!==============================================================================

END MODULE kind_parameters
