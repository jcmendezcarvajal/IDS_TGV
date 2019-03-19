!	******************* The solution flux (U)  ***************

Subroutine Flux_U
use Variables
use omp_lib
implicit none

	alpha = gamma*(gamma-1.0)*cM_inf**2

	!$OMP PARALLEL PRIVATE (i,j)
	!$OMP DO
	do j = 1,Jmax
	do i = 1,Imax
	H_old(i,j) = r_old(i,j)*(T_old(i,j)/alpha +0.5*(u_old(i,j)**2+&
                  v_old(i,j)**2))
	enddo
	enddo
	!$OMP END DO

	!$OMP DO
	do j = 1,Jmax
	do i = 1,Imax
	U1_old(i,j) = r_old(i,j)
	U2_old(i,j) = r_old(i,j)*u_old(i,j)
	U3_old(i,j) = r_old(i,j)*v_old(i,j)
	U4_old(i,j) = H_old(i,j)
  enddo
  enddo
	!$OMP END DO
	!$OMP END PARALLEL

End Subroutine Flux_U
