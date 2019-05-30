!	*******************  Boundary Conditions **************
Subroutine BC_Supersonic_Temp
  use Variables
  implicit none


  !set the Periodic BC for the vertical direction
  do j = 1,Jmax
    r_temp(1,j) = r_temp(Imax-1,j)
    u_temp(1,j) = u_temp(Imax-1,j)
    v_temp(1,j) = v_temp(Imax-1,j)
    T_temp(1,j) = T_temp(Imax-1,j)

    r_temp(Imax,j) = r_temp(2,j)
    u_temp(Imax,j) = u_temp(2,j)
    v_temp(Imax,j) = v_temp(2,j)
    T_temp(Imax,j) = T_temp(2,j)
  enddo



  !set the Periodic BC for the horizontal direction
  do i = 1,Imax
    r_temp(i,jmax) = r_temp(i,2)
    u_temp(i,jmax) = u_temp(i,2)
    v_temp(i,jmax) = v_temp(i,2)
    T_temp(i,jmax) = T_temp(i,2)

    r_temp(i,1) = r_temp(i,jmax-1)
    u_temp(i,1) = u_temp(i,jmax-1)
    v_temp(i,1) = v_temp(i,jmax-1)
    T_temp(i,1) = T_temp(i,jmax-1)
  enddo

End Subroutine BC_Supersonic_Temp
