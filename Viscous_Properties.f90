!      *****************  Evaluate the  viscous properties ****************
Subroutine Viscous_Properties
use Variables
use omp_lib
implicit none

  !$OMP PARALLEL PRIVATE (i,j)
  !$OMP DO
	do j = 1,Jmax
	do i = 1,Imax
	vMu(i,j)	= ((T_old(i,j))**1.5)*((1.0+110.0/T_inf)/&
  (T_old(i,j)+110.0/T_inf))
  enddo
  enddo
  !$OMP END DO

  !$OMP DO
  do j = 1,Jmax
	do i = 1,Imax
  tK(i,j)		= vMu(i,j)
  vNu(i,j)	= vMu(i,j)/r_old(i,j)
  a_old(i,j)  = sqrt(T_old(i,j))/cM_inf
  enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  End Subroutine Viscous_Properties
