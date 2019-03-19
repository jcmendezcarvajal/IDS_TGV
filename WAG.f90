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
		Pres(i,j) = 101325 + ((COS(2*x(i,j)) + COS(2*y(i,j)))*r_inf/(4.0))*9.80
	enddo
	enddo

  u_inf = MAXVAL(u_old(2:Imax-1,2:jmax-1)) !Max velocity
  P_inf = MAXVAL(Pres(2:Imax-1,2:jmax-1)) !Max Pressure
	P_inf = P_inf/1000
	T_inf = P_inf/(0.287*r_inf)

	do j = 2,Jmax-1
	do i = 2,Imax-1
		T_old(i,j) = Pres(i,j)/(0.287*r_inf) ! This is Temperature (Dimensional)
	enddo
	enddo

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
  vMu_inf = (r_inf*u_inf*x_actual)/Re_L       ! Reynolds number


  Call BC_Supersonic

End Subroutine WAG
