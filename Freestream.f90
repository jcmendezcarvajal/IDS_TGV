!********* freesream conditions and the constants **************
Subroutine Freestream
  use Variables
  implicit none

  gamma   = 1.4
  R_gas	= 287.0 !AIR
  Cp = gamma*R_gas/(gamma - 1)
  r_inf   = 0.5323
  Pr_inf = 0.711

  ! cM_inf  = 0.5 to be defined in WAG
  ! T_inf   = 362 to be defined in WAG

   Re_L = 500 !to be defined in WAG

  ! P_inf	= r_inf*R_gas*T_inf       !AIR
  ! a_inf   = sqrt(gamma*R_gas*T_inf) to be defined in WAG
  ! u_inf   = a_inf*cM_inf to be defined in WAG


  x_actual    = 2.0*Pi()!  in meters

  write(*,*)  Re_L,  'Re_L Number'
  write(*,*)  P_inf,  'P_inf'
  write(*,*)  u_inf,  'u_inf'
  write(*,*)  x_actual,  'x_actual'



  !	write(*,*) dx, ddy(3), dx/dy, 'dx, dy, aspect ratio ...... GOT IT?'

End	Subroutine Freestream

!     ****************** Generating the Grids ********************
