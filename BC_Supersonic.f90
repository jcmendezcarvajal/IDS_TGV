!	*******************  Boundary Conditions **************
Subroutine BC_Supersonic
  use Variables
  implicit none


  !set the Periodic BC for the vertical direction
  do j = 1,Jmax
    r_old(1,j) = r_old(Imax-1,j)
    u_old(1,j) = u_old(Imax-1,j)
    v_old(1,j) = v_old(Imax-1,j)
    T_old(1,j) = T_old(Imax-1,j)

    r_old(Imax,j) = r_old(2,j)
    u_old(Imax,j) = u_old(2,j)
    v_old(Imax,j) = v_old(2,j)
    T_old(Imax,j) = T_old(2,j)
  enddo



  !set the Periodic BC for the horizontal direction
  do i = 1,Imax
    r_old(i,jmax) = r_old(i,2)
    u_old(i,jmax) = u_old(i,2)
    v_old(i,jmax) = v_old(i,2)
    T_old(i,jmax) = T_old(i,2)

    r_old(i,1) = r_old(i,jmax-1)
    u_old(i,1) = u_old(i,jmax-1)
    v_old(i,1) = v_old(i,jmax-1)
    T_old(i,1) = T_old(i,jmax-1)
  enddo

End Subroutine BC_Supersonic
