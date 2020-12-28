MODULE constants_mod
    
  implicit none
  integer  ,parameter :: nVar = 3
  
  integer  ,parameter :: i2  = 2
  integer  ,parameter :: i4  = 4
  integer  ,parameter :: i8  = 8
  integer  ,parameter :: i16 = 16
  integer  ,parameter :: r2  = 2
  integer  ,parameter :: r4  = 4
  integer  ,parameter :: r8  = 8
  integer  ,parameter :: r16 = 16
  
  integer  ,parameter :: i_kind = i4
  integer  ,parameter :: r_kind = r8

  real(r_kind),parameter :: pi        = 2.*asin(1.)
  real(r_kind),parameter :: D2R       = PI/180.    ! convert degree into radian
  real(r_kind),parameter :: R2D       = 180./PI    ! convert radian into degree
  real(r_kind),parameter :: FillValue = -999999999999999.  
  
  real(r_kind),parameter :: piq       = 0.25_r16*pi
  real(r_kind),parameter :: pih       = 0.5_r16 *pi
  real(r_kind),parameter :: pi2       = 2._r16  *pi
  real(r_kind),parameter :: Inf       = huge(Inf)
  real(r_kind),parameter :: tolerance = 10.**(-r_kind*2+1) ! Tolerant parameter, choose from tiny(1._r_kind), 10.**(-r_kind*2+1)
  
  real(r_kind),parameter :: gravity   = 9.80616
  
  !! thermal dynamic constants
  !real(r_kind),parameter :: Rd        = 287.0583
  !real(r_kind),parameter :: Rv        = 461.5233
  !real(r_kind),parameter :: p0        = 100000. ! base pressure in Pa
  !real(r_kind),parameter :: Md        = 0.0289644
  !real(r_kind),parameter :: Mv        = 0.01801528
  !real(r_kind),parameter :: eq        = Md / Mv - 1.
  !real(r_kind),parameter :: Cpd       = 1005.
  !real(r_kind),parameter :: Cpv       = 1846.
  !real(r_kind),parameter :: Cvd       = 718.
  !real(r_kind),parameter :: Cvv       = 1385.
  
  ! thermal dynamic constants
  real(r_kind),parameter :: Rd        = 287.
  real(r_kind),parameter :: Rv        = 461.5233
  real(r_kind),parameter :: p0        = 100000. ! base pressure in Pa
  real(r_kind),parameter :: Md        = 0.0289644
  real(r_kind),parameter :: Mv        = 0.01801528
  real(r_kind),parameter :: eq        = Md / Mv - 1.
  real(r_kind),parameter :: Cpd       = 1004.5
  real(r_kind),parameter :: Cpv       = 1846.
  real(r_kind),parameter :: Cvd       = 717.5
  real(r_kind),parameter :: Cvv       = 1385.
END MODULE constants_mod
