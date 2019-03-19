!      ******************  Calculating the Time Step **************************
Subroutine Time_Step
	use Variables
	use omp_lib
	implicit none

	!$OMP PARALLEL PRIVATE(i,j)
	!$OMP DO
	do j = 2,Jmax-1
		do i = 2,Imax-1
			t1(i,j) = abs(u_old(i,j))/ddx(i)+abs(v_old(i,j))/ddy(j)
			t2(i,j) = a_old(i,j)*sqrt(1.0/ddx(i)**2+1.0/ddy(j)**2)
			tNu_max(i,j) =max((4.0/3.0)*vMu(i,j)/(r_old(i,j)*Re_L),&
			gamma*vMu(i,j)/(r_old(i,j)*Pr_inf*Re_L))
		enddo
	enddo
	!$OMP END DO

	!$OMP DO
	do j = 2,Jmax-1
		do i = 2,Imax-1
			small_time(i,j) = (t1(i,j)+t2(i,j)+2.0*tNu_max(i,j)*&
			(1.0/(ddx(i)**2)+1.0/(ddy(j)**2)))**(-1)
		enddo
	enddo
	!$OMP END PARALLEL

	delta_t = small_time(2,2)

	do j = 3,Jmax-1
		do i = 3,Imax-1
			delta_t = min(delta_t,small_time(i,j))
		enddo
	enddo

	delta_t = 0.7*delta_t

End Subroutine Time_Step
