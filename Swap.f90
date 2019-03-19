!  Switch the old values to the new values
Subroutine Swap
use Variables
use omp_lib
implicit none

!$OMP PARALLEL DO PRIVATE(i,j)
do j = 1,Jmax
	do i = 1,Imax
		r_old(i,j) = r_temp(i,j)
		u_old(i,j) = u_temp(i,j)
		v_old(i,j) = v_temp(i,j)
		T_old(i,j) = T_temp(i,j)
	enddo
	enddo
!$OMP END PARALLEL DO
End Subroutine Swap
