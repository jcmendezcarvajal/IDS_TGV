!     ****************** Generating the Grids ********************
Subroutine GridNonDimensional
  use Variables
  implicit none

  yL		            = 1.0
  xL		            = 1.0

  ! Evaluate the incremental change in x and y direction
  dx = xL/(Imax-1)
  dy = yL/(Jmax-1)

  ! The numerical domain information

  !Define the Grid points
  do i = 1, Imax
    xx(i) = 0.0 + (i-1)*dx       !Regular
  enddo

  do j = 1, Jmax
    yy(j) = 0.0 + (j-1)*dy
  enddo

  !Transfer to Old System
  do j = 1,Jmax
    do i = 1,Imax
      y(i,j) = yy(j)
      x(i,j) = xx(i)
    enddo
  enddo

  !Evaluate the dxi and dyj only
  ddx(1:Imax-1)   = xx(2:Imax) - xx(1:Imax-1)
  ddx(Imax)   = ddx(Imax-1)

  ddy(1:Jmax-1)   = yy(2:Jmax) - yy(1:Jmax-1)
  ddy(Jmax)   = ddy(Jmax-1)

  write(*,*) ddx(1), ddy(1)," aspect ratio ...... GOT IT?"
  
End Subroutine GridNonDimensional

Subroutine Grid
  !This takes care of the units for the dimanional part of the harmonic functions
  use Variables
  implicit none

  yL		            = x_actual
  xL		            = x_actual

  ! Evaluate the incremental change in x and y direction
  dx = xL/(Imax-1)
  dy = yL/(Jmax-1)

  ! The numerical domain information

  !Define the Grid points
  do i = 1, Imax
    xx(i) = 0.0 + (i-1)*dx       !Regular
  enddo

  do j = 1, Jmax
    yy(j) = 0.0 + (j-1)*dy
  enddo

  !Transfer to Old System
  do j = 1,Jmax
    do i = 1,Imax
      y(i,j) = yy(j)
      x(i,j) = xx(i)
    enddo
  enddo

  !Evaluate the dxi and dyj only
  ddx(1:Imax-1)   = xx(2:Imax) - xx(1:Imax-1)
  ddx(Imax)   = ddx(Imax-1)

  ddy(1:Jmax-1)   = yy(2:Jmax) - yy(1:Jmax-1)
  ddy(Jmax)   = ddy(Jmax-1)

  write(*,*) ddx(1), ddy(1)," aspect ratio ...... GOT IT?"
End Subroutine Grid
