MODULE constants_mod
    
  implicit none
  INTEGER,PARAMETER :: DOF       = 3    ! Degree of Freedoms within a 1D element

  REAL,PARAMETER    :: gravity   = 9.80616
  REAL,PARAMETER    :: pi        = 2.*asin(1.)
  
  REAL,PARAMETER    :: radius    = 6371229.
  REAL,PARAMETER    :: D2R       = PI/180.    ! convert degree into radian
  REAL,PARAMETER    :: R2D       = 180./PI    ! convert radian into degree
  REAL,PARAMETER    :: Omega     = 7.292E-5

  REAL,PARAMETER    :: FillValue = -9999999999999999.  
  
  integer  ,parameter :: i2  = 2
  integer  ,parameter :: i4  = 4
  integer  ,parameter :: i8  = 8
  integer  ,parameter :: i16 = 16
  integer  ,parameter :: r2  = 2
  integer  ,parameter :: r4  = 4
  integer  ,parameter :: r8  = 8
  integer  ,parameter :: r16 = 16
  
  real(r8),parameter :: piq       = 0.25_r16*pi
  real(r8),parameter :: pih       = 0.5_r16 *pi
  real(r8),parameter :: pi2       = 2._r16  *pi
  real(r8),parameter :: Inf       = huge(Inf)
  real(r8),parameter :: tolerance = 1.E-15         ! tolerant parameter
  
END MODULE constants_mod
