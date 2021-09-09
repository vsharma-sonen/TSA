!+ Generic support: procedures to manipulate date and time
!
!+****************************************************************************
MODULE support_datetime
!=============================================================================
!
! Collection of procedures to manipulate data and time information:
!  > transformation between calendar, julian date, day of year
!  > name of months
!  > sum and difference of dates
!  > manipulate date string
!
!
! Current Code Owner
! ------------------
! Jean-Marie Bettems
! MeteoSwiss, Zurich, Switzerland
! jean-marie.bettems@meteoswiss.ch
!
!-----------------------------------------------------------------------------
USE fxtr_definition, ONLY: &
    kind_idate,            &
    iundef



IMPLICIT NONE
PRIVATE



! Public entities
!================
! transformation between calendar, julian date, day of year
PUBLIC :: calendar, julian, day_of_year
! name of months
PUBLIC :: month_name3
! sum and difference of dates
PUBLIC :: time_sum, date_sum, date_delta
! manipulate date string
PUBLIC :: check_date, split_date, format_date, iso_date



CONTAINS

!+****************************************************************************
FUNCTION split_date(date_and_time)
!=============================================================================
!
! Split integer date_and_time, returns (YYYY,MM,DD,hh,mm)
!   YYYY: year, MM: month, DD: day, hh: hour, mm: minute
!
! Format of date_and_time may be either YYYYMMDDhhmm or YYYYMMDDhh.
! If date_and_time is undefined (<0), it is first reset to 0.
!
!-----------------------------------------------------------------------------
  ! Dummy arguments
  INTEGER(KIND=kind_idate),INTENT(IN)  :: date_and_time
  INTEGER, DIMENSION(5)                :: split_date

  ! Local variables
  INTEGER(KIND=kind_idate)      :: n_date_and_time
  INTEGER                       :: remain


  ! Normalize date_and_time value
  IF      ( date_and_time < 0_kind_idate ) THEN
    ! ... take care of undefined values
    n_date_and_time = 0_kind_idate
  ELSE IF ( date_and_time < 10000000000_kind_idate ) THEN
    ! ... take care of missing minutes
    n_date_and_time = date_and_time * 100
  ELSE 
    n_date_and_time = date_and_time
  ENDIF

  ! Split date
  split_date(1) = n_date_and_time / 100000000_kind_idate
  remain = MOD(n_date_and_time, 100000000_kind_idate)
  split_date(2) = remain / 1000000  ;  remain = MOD(remain, 1000000)
  split_date(3) = remain / 10000    ;  remain = MOD(remain, 10000)
  split_date(4) = remain / 100    ;  split_date(5) = MOD(remain, 100)

END FUNCTION split_date



!+****************************************************************************
FUNCTION check_date(date_and_time)
!=============================================================================
!
! Check that date_and_time is compatible with YYYYMMDDHHmm,
!   YYYY: year, MM: month, DD: day, hh: hour, mm: minute
! return .false. otherwise
!
!-----------------------------------------------------------------------------
  ! Dummy arguments
  INTEGER(KIND=kind_idate),INTENT(IN)  :: date_and_time
  LOGICAL                              :: check_date

  ! Local variables
  INTEGER, DIMENSION(5)                 :: ddtt


  check_date = .TRUE.

  IF ( date_and_time <= 0 .OR. date_and_time > 999999999999_kind_idate ) THEN
    check_date = .FALSE.
  ELSE
    ddtt(:) = split_date(date_and_time)
    IF ( ddtt(1) < 1900 .OR.                                       &
         (ddtt(2) < 1 .OR. ddtt(2) > 12) .OR.                      &
         (ddtt(3) < 1 .OR. ddtt(3) > 31) .OR.                      &
         ddtt(4) > 23                    .OR.                      &
         ddtt(5) > 59                          ) check_date = .FALSE.
  ENDIF

END FUNCTION check_date



!+****************************************************************************
FUNCTION iso_date(date_and_time)
!=============================================================================
!
! Return date_and_time formatted according to ISO specification
!
!-----------------------------------------------------------------------------
  ! Dummy arguments
  INTEGER(KIND=kind_idate), INTENT(IN)  :: date_and_time
  CHARACTER(LEN=23)                     :: iso_date

  ! Local parameters
  ! ISO format for date and time
  CHARACTER(LEN=*), PARAMETER            :: iso_date_fmt =          &
               '(I4,2("-",I2.2),1X,I2,":",I2.2,1X,SP,I3,":",SS,I2.2)'
  ! Local variables
  INTEGER, DIMENSION(5)                   :: ddtt


  ddtt(:) = split_date(date_and_time)
  WRITE(iso_date,iso_date_fmt) ddtt(1), ddtt(2), ddtt(3), ddtt(4), ddtt(5), 0, 0


END FUNCTION iso_date



!+****************************************************************************
FUNCTION format_date(date_and_time, format)
!=============================================================================
!
! Return formatted date_and_time
!
! Supported formats:
!   format = "std8" :   date_and_time  ->  YYMMDDhh
!   format = "std10":   date_and_time  ->  YYYYMMDDhh
!   format = "std12":   date_and_time  ->  YYYYMMDDhhmm
!   format = "dat8" :   date_and_time  ->  YYYYMMDD
!   format = "mm10" :   date_and_time  ->  YYMMDDhhmm
!   format = "mm2"  :   date_and_time  ->  mm
!
! If date_and_time is undefined (<0), it is first reset to 0
!
!-----------------------------------------------------------------------------
  ! Dummy arguments
  INTEGER(KIND=kind_idate),INTENT(IN)  :: date_and_time
  CHARACTER(LEN=*),INTENT(IN)          :: format
  INTEGER(KIND=kind_idate)             :: format_date

  ! Local variables
  INTEGER(KIND=kind_idate)             :: n_date_and_time


  ! Take care of undefined values
  n_date_and_time = MAX( date_and_time , 0_kind_idate )

  ! Format date
  SELECT CASE ( format )
  CASE ("std8")
  ! YYMMDDhh
    format_date = MOD( n_date_and_time/100 , 100000000_kind_idate )

  CASE ("std10")
  ! YYYYMMDDhh
    format_date = n_date_and_time/100

  CASE ("std12")
  ! YYYYMMDDhhmm
    format_date = n_date_and_time

  CASE ("dat8")
  ! YYYYMMDD
    format_date = n_date_and_time/10000

  CASE ("mm")
  ! mm
    format_date = MOD( n_date_and_time , 100_kind_idate )

  CASE ("mm10")
  ! YYMMDDhhmm
    format_date = MOD( n_date_and_time/100 , 10000000000_kind_idate )

  CASE DEFAULT
    format_date = n_date_and_time

  END SELECT

END FUNCTION format_date



!+****************************************************************************
FUNCTION month_name3(mm)
!=============================================================================
!
! Returns 3-characters month name
!
!-----------------------------------------------------------------------------
  ! Dummy arguments
  INTEGER,INTENT(IN)               :: mm
  CHARACTER(LEN=3)                 :: month_name3

  ! Local parameters
  CHARACTER(LEN=*),PARAMETER       :: month_names(12) =  (/                   &
    "Jan","Feb","Mar","Apr","May","Jun", "Jul","Aug","Sep","Oct","Nov","Dec" /)

  IF ( mm > 0 .AND. mm < 13 ) THEN
    month_name3 = month_names(mm)
  ELSE
    month_name3 = "???"
  ENDIF

END FUNCTION month_name3



!+****************************************************************************
SUBROUTINE calendar(julian, yyyy, mm, dd)
!=============================================================================
!
! Converts Julian day to calendar date. Algorithm according to
! Numerical Recipes.
! Sequence of arguments adapted to ISO 8601 standard date format (YYYY-MM-DD).
!
! See also: julian
!
! ____________________________________________________________________________
! I/O   Name            Type    Description
! ____________________________________________________________________________
! I     julian          I       Julian day
! O     yyyy            I       Year
! O     mm              I       Month
! O     dd              I       Day
!
! Version history:
! 1.0   Basic version
!-----------------------------------------------------------------------------
  ! Dummy arguments
  INTEGER,INTENT(IN)  :: julian
  INTEGER,INTENT(OUT) :: yyyy, mm, dd
  
  ! Local parameters
  INTEGER,PARAMETER :: igreg = 2299161

  ! Local variables
  INTEGER :: ja, jalpha, jb, jc, jd, je
  
  IF ( julian >= igreg ) THEN
    jalpha = INT(((julian - 1867216) - 0.25) / 36524.25)
    ja = julian + 1 + jalpha - INT(0.25 * jalpha)
  ELSE
    ja = julian
  ENDIF
  jb = ja + 1524
  jc = INT(6680. + ((jb-2439870)-122.1)/365.25)
  jd = 365 * jc + INT(0.25 * jc)
  je = INT((jb - jd) / 30.6001)
  dd = jb - jd - INT(30.6001 * je)
  mm = je - 1
  IF ( mm > 12 ) mm = mm - 12
  yyyy = jc - 4715
  IF ( mm > 2 ) yyyy = yyyy - 1
  IF ( yyyy <= 0 ) yyyy = yyyy - 1

END SUBROUTINE calendar



!+****************************************************************************
FUNCTION day_of_year(yyyy, mm, dd)
!=============================================================================
!
! Compute day of current year
! Counting starts with 0:   day_of_year(yyyy,01,01) := 0
!
!-----------------------------------------------------------------------------
  ! Dummy arguments
  INTEGER                     :: day_of_year
  INTEGER, INTENT(IN)         :: yyyy, mm, dd

  ! Local variables
  LOGICAL                     :: is_leap

  ! Local data statements
  INTEGER, DIMENSION(12)      :: month, month_l
  DATA month    / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /
  DATA month_l  / 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /



  ! Leap year ?
  is_leap = ( MOD(yyyy,4) == 0 )
  IF ( MOD(yyyy,100) == 0 ) THEN
    IF ( MOD(yyyy,400) == 0 ) THEN
      is_leap = .TRUE.
    ELSE
      is_leap = .FALSE.
    ENDIF
  ENDIF

  ! Compute number of days since January first
  IF     ( mm /= 1 .AND. .NOT. is_leap ) THEN
    day_of_year = SUM(month(1:(mm-1))) + dd - 1
  ELSEIF ( mm /= 1 .AND. is_leap ) THEN
    day_of_year = SUM(month_l(1:(mm-1))) + dd - 1
  ELSE
    day_of_year = dd - 1
  ENDIF


END FUNCTION day_of_year



!+****************************************************************************
FUNCTION julian(yyyy, mm, dd)
!=============================================================================
!
! Converts date to Julian day. Algorithm according to Numerical Recipes.
! Sequence of arguments adapted to ISO 8601 standard date format (YYYY-MM-DD).
!
! See also: calendar
!
! ____________________________________________________________________________
! I/O   Name            Type    Description
! ____________________________________________________________________________
! I     yyyy            I       Year
! I     mm              I       Month
! I     dd              I       Day
! O     julian          I       Julian day
!
! Version history:
! 1.0   Basic version
!-----------------------------------------------------------------------------
  ! Dummy arguments
  INTEGER             :: julian
  INTEGER,INTENT(IN)  :: yyyy, mm, dd
  
  ! Local parameters
  INTEGER, PARAMETER  :: igreg = 15 + 31 * (10 + 12 * 1582)

  ! Local variables
  INTEGER             :: ja, jm, jy
  
  jy = yyyy
  IF ( jy < 0 ) jy = jy + 1
  IF ( mm > 2 ) THEN
    jm = mm + 1
  ELSE
    jy = jy - 1
    jm = mm + 13
  ENDIF
  julian = INT(365.25 * jy) + INT(30.6001 * jm) + dd + 1720995
  IF ( dd + 31 * (mm + 12 * yyyy) >= igreg ) THEN
    ja = INT(0.01 * jy)
    julian = julian + 2 - ja + INT(0.25 * ja)
  ENDIF

END FUNCTION julian



!+****************************************************************************
FUNCTION date_delta(date1,date2)
!=============================================================================
!
! Returns the difference in minutes between date1 and date2
!   date1 and date2 are represented as YYYYMMDDhhmm
!     YYYY: 4-digits year, MM: 2-digits month, DD: 2-digits day
!     hh: 2-digits hour, mm: 2-digits minutes
!
!-----------------------------------------------------------------------------
  ! Dummy arguments
  INTEGER                               :: date_delta
  INTEGER(KIND=kind_idate), INTENT(IN)  :: date1, date2

  ! Local variables
  INTEGER, DIMENSION(5)                 :: ddtt1, ddtt2


  date_delta = iundef
  IF ( date1 /= iundef .AND. date2 /= iundef ) THEN
    ddtt1(:) = split_date(date1)
    ddtt2(:) = split_date(date2)
    date_delta = ( (julian(ddtt1(1),ddtt1(2),ddtt1(3))*24 + ddtt1(4))*60 + ddtt1(5) ) - &
                 ( (julian(ddtt2(1),ddtt2(2),ddtt2(3))*24 + ddtt2(4))*60 + ddtt2(5) )  
  ENDIF

END FUNCTION date_delta



!+****************************************************************************
FUNCTION date_sum(date,delta_minute)
!=============================================================================
!
! Returns the sum of date and delta_minute.
!   date and date_sum are represented as YYYYMMDDhhmm
!     YYYY: 4-digits year, MM: 2-digits month, DD: 2-digits day
!     hh: 2-digits hour, mm: 2-digits minutes
!   delta_minute is expressed in minutes
!
!-----------------------------------------------------------------------------
  ! Dummy arguments
  INTEGER(KIND=kind_idate)              :: date_sum
  INTEGER(KIND=kind_idate), INTENT(IN)  :: date
  INTEGER, INTENT(IN)                   :: delta_minute

  ! Local variables
  INTEGER, DIMENSION(8)                 :: t1, t2, t3
  INTEGER                               :: year, month, day, hour, minute, remain


  date_sum = iundef
  IF ( date /= iundef .AND. delta_minute /= iundef ) THEN
    year  = date / 100000000_kind_idate  
    remain = MOD(date, 100000000_kind_idate)
    month = remain / 1000000  ; remain = MOD(remain, 1000000)
    day   = remain / 10000    ; remain = MOD(remain, 10000)
    hour  = remain / 100      ; minute = MOD(remain, 100)

    t1 = (/ year, month, day, 0, hour, minute, 0, 0 /)
    t2 = (/ 0, 0, 0, 0, 0, delta_minute, 0, 0 /)
    t3 = time_sum(t1, t2)
    date_sum = t3(1)*100000000_kind_idate + t3(2)*1000000_kind_idate + &
               t3(3)*10000_kind_idate + t3(5)*100_kind_idate + t3(6)
  ENDIF

END FUNCTION date_sum



!+****************************************************************************
FUNCTION time_sum(time1,time2)
!=============================================================================
!
! Adds two times.
!
! See also: julian, calendar
!
! ____________________________________________________________________________
! I/O   Name            Type    Description
! ____________________________________________________________________________
! I     time1           I(:)    Date-time array
! I     time2           I(:)    Date-time array
! O     time_sum        I(8)    Date-time array, sum of time1 and time2
!
! Element ordering in date-time array according to DATE_AND_TIME intrinsic:
! Year, month, day, minutes difference form UTC, hour, minutes, seconds,
! milliseconds.
! One of time1 and time2 must be a delta time, 
! i.e. month == 0 and minutes difference from UTC == 0
!
! CAUTION:
! Some years have an additional second at the end. This is NOT considered here.
!
! Author: Pirmin Kaufmann, SMI
!
! Version history:
! 1.0   Basic version
!-----------------------------------------------------------------------------
  ! Dummy arguments
  INTEGER, DIMENSION(8)            :: time_sum
  INTEGER, DIMENSION(:),INTENT(IN) :: time1, time2
  
  ! Local variables
  INTEGER, DIMENSION(8)            :: t1, t2, tsum
  INTEGER                          :: jday
  LOGICAL, DIMENSION(2)            :: abstime
  


  t1 = 0
  t2 = 0
  t1(1:size(time1)) = time1
  t2(1:size(time2)) = time2
  
  abstime(1) = .NOT. (t1(1) == 0 .AND. t1(2) == 0 .AND. t1(4) == 0)
  abstime(2) = .NOT. (t2(1) == 0 .AND. t2(2) == 0 .AND. t2(4) == 0)

  ! One variable must be delta-time, 
  ! i.e. year, month and difference from UTC == 0
  IF ( ALL(abstime) ) THEN
    time_sum = -HUGE(0)
    RETURN
  ENDIF
  IF ( abstime(1) ) THEN
    t1(3) = julian(t1(1),t1(2),t1(3))
    t1(1:2) = 0
  ENDIF
  IF ( abstime(2) ) THEN
    t2(3) = julian(t2(1),t2(2),t2(3))
    t2(1:2) = 0
  ENDIF
  
  tsum = t1 + t2
  tsum(7) = tsum(7) + FLOOR((tsum(8)+0.5)/1000.)
  tsum(8) = MODULO(tsum(8), 1000)
  tsum(6) = tsum(6) + FLOOR((tsum(7)+0.5)/60.)
  tsum(7) = MODULO(tsum(7), 60)
  tsum(5) = tsum(5) + FLOOR((tsum(6)+0.5)/60.)
  tsum(6) = MODULO(tsum(6), 60)
  tsum(3) = tsum(3) + FLOOR((tsum(5)+0.5)/24.)
  tsum(5) = MODULO(tsum(5), 24)
  
  IF ( abstime(1) .NEQV. abstime(2) ) THEN
    jday = tsum(3)
    CALL calendar(jday,tsum(1),tsum(2),tsum(3))
  ENDIF
  
  time_sum = tsum

END FUNCTION time_sum


END MODULE support_datetime
