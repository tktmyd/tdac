!! ------------------------------------------------------------------------- !!
!>
!! Common parameters used in llw2d and tdac
!<
!! --
module m_params

  implicit none
  public

  integer, parameter, public :: DP = selected_real_kind(13)
  integer, parameter, public :: STDERR  = 0
  integer, parameter, public :: STDOUT  = 6
  real(DP), parameter :: g0 = 9.80665_DP   !<  Gravity Constant

  !! finite differnce method parametrs
  integer,  parameter :: nx = 200         !<  grid number (NS)
  integer,  parameter :: ny = 200         !<  grid number (EW)
  integer,  parameter :: no = 36           !<  Number of stations
  real(DP), parameter :: dx = 2000._DP     !<  grid width of x-direction (m)
  real(DP), parameter :: dy = 2000._DP     !<  grid width of y-direction (m)
  real(DP), parameter :: dt = 1.0_DP       !<  time step width (sec)

end module m_params
!! ------------------------------------------------------------------------- !!
