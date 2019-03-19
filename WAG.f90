!	******************* WAG  **************

! Setup an initial guess for internal points

Subroutine WAG
use Variables
implicit none

! These are dimensional quantities
	do j = 2,Jmax-1
	do i = 2,Imax-1
    r_old(i,j) = 1.0
		u_old(i,j) = SIN(x(i,j))*COS(y(i,j))
		v_old(i,j) = -COS(x(i,j))*SIN(y(i,j))
		T_old(i,j) = (COS(2*x(i,j)) + COS(2*y(i,j)))*r_inf/(4.0) ! Visualize this as Pressure
	enddo
	enddo

	do j = 2,Jmax-1
	do i = 2,Imax-1
		T_old(i,j) = T_old(i,j)/(0.287*r_inf) ! This is Temperature (Dimensional)
	enddo
	enddo

  u_inf = MAXVAL(u_old(2:Imax-1,2:jmax-1)) !Max velocity
  T_inf = MAXVAL(T_old(2:Imax-1,2:jmax-1)) !Max Temperature
  print*,T_inf
  ! This are non dimanional quantities

	do j = 2,Jmax-1
	do i = 2,Imax-1
    r_old(i,j) = 1.0
		u_old(i,j) =  u_old(i,j)/u_inf
		v_old(i,j) = 	v_old(i,j)/u_inf
		T_old(i,j) =  T_old(i,j)/T_inf
	enddo
	enddo

  a_inf   = sqrt(gamma*R_gas*T_inf)
	cM_inf = u_inf/a_inf
	Re_L    = (r_inf*u_inf*x_actual)/vMu_inf        ! Reynolds number


  Call BC_Supersonic

End Subroutine WAG
