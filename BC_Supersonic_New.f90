!	*******************  Boundary Conditions **************
Subroutine BC_Supersonic_New
  use Variables
  implicit none


  !set the Periodic BC for the vertical direction
  do j = 1,Jmax
    r_new(1,j) = r_new(Imax-1,j)
    u_new(1,j) = u_new(Imax-1,j)
    v_new(1,j) = v_new(Imax-1,j)
    T_new(1,j) = T_new(Imax-1,j)

    r_new(Imax,j) = r_new(2,j)
    u_new(Imax,j) = u_new(2,j)
    v_new(Imax,j) = v_new(2,j)
    T_new(Imax,j) = T_new(2,j)
  enddo



  !set the Periodic BC for the horizontal direction
  do i = 1,Imax
    r_new(i,jmax) = r_new(i,2)
    u_new(i,jmax) = u_new(i,2)
    v_new(i,jmax) = v_new(i,2)
    T_new(i,jmax) = T_new(i,2)

    r_new(i,1) = r_new(i,jmax-1)
    u_new(i,1) = u_new(i,jmax-1)
    v_new(i,1) = v_new(i,jmax-1)
    T_new(i,1) = T_new(i,jmax-1)
  enddo


End Subroutine BC_Supersonic_New
