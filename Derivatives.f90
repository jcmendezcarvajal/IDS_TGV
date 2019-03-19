! Calculate the rate of change of solution flux (Um) of (a,b,c,d) cells
Subroutine Derivatives
use Variables
use omp_lib
implicit none

!$OMP PARALLEL PRIVATE (i,j)
!$OMP DO
do j =2,Jmax-1
do i =2,Imax-1
	! cell(a)

    dU1dt_a(i,j) =  -(E1_inv_a_right(i,j)-E1_inv_a_left(i,j))/ddx(i-1)      &
                    -(F1_inv_a_upper(i,j)-F1_inv_a_lower(i,j))/ddy(j-1)

    dU2dt_a(i,j) =  -(E2_inv_a_right(i,j)-E2_inv_a_left(i,j))/ddx(i-1)      &
                     -(F2_inv_a_upper(i,j)-F2_inv_a_lower(i,j))/ddy(j-1)


    dU3dt_a(i,j) =  -(E3_inv_a_right(i,j)-E3_inv_a_left(i,j))/ddx(i-1)      &
                     -(F3_inv_a_upper(i,j)-F3_inv_a_lower(i,j))/ddy(j-1)


	  dU4dt_a(i,j) =  -(E4_inv_a_right(i,j)-E4_inv_a_left(i,j))/ddx(i-1)      &
                     -(F4_inv_a_upper(i,j)-F4_inv_a_lower(i,j))/ddy(j-1)


	! cell(b)

    dU1dt_b(i,j) =  -(E1_inv_b_right(i,j) - E1_inv_b_left(i,j)) /ddx(i)       &
                    -(F1_inv_b_upper(i,j) - F1_inv_b_lower(i,j))/ddy(j-1)

    dU2dt_b(i,j) =  -(E2_inv_b_right(i,j) - E2_inv_b_left(i,j)) /ddx(i)       &
                    -(F2_inv_b_upper(i,j) - F2_inv_b_lower(i,j))/ddy(j-1)


	dU3dt_b(i,j) =  -(E3_inv_b_right(i,j)-E3_inv_b_left(i,j))/ddx(i)        &
                     -(F3_inv_b_upper(i,j)-F3_inv_b_lower(i,j))/ddy(j-1)

	dU4dt_b(i,j) =  -(E4_inv_b_right(i,j)-E4_inv_b_left(i,j))/ddx(i)        &
                     -(F4_inv_b_upper(i,j)-F4_inv_b_lower(i,j))/ddy(j-1)

	! cell(c)

	dU1dt_c(i,j) =  -(E1_inv_c_right(i,j)-E1_inv_c_left(i,j))/ddx(i) &
                     -(F1_inv_c_upper(i,j)-F1_inv_c_lower(i,j))/ddy(j)

      dU2dt_c(i,j) =  -(E2_inv_c_right(i,j)-E2_inv_c_left(i,j))/ddx(i) &
                     -(F2_inv_c_upper(i,j)-F2_inv_c_lower(i,j))/ddy(j)

	dU3dt_c(i,j) =  -(E3_inv_c_right(i,j)-E3_inv_c_left(i,j))/ddx(i) &
                     -(F3_inv_c_upper(i,j)-F3_inv_c_lower(i,j))/ddy(i)

	dU4dt_c(i,j) =  -(E4_inv_c_right(i,j)-E4_inv_c_left(i,j))/ddx(i) &
                     -(F4_inv_c_upper(i,j)-F4_inv_c_lower(i,j))/ddy(j)

	! cell(d)

	dU1dt_d(i,j) =  -(E1_inv_d_right(i,j)-E1_inv_d_left(i,j))/ddx(i-1) &
                     -(F1_inv_d_upper(i,j)-F1_inv_d_lower(i,j))/ddy(j)

      dU2dt_d(i,j) = -(E2_inv_d_right(i,j)-E2_inv_d_left(i,j))/ddx(i-1) &
                    -(F2_inv_d_upper(i,j)-F2_inv_d_lower(i,j))/ddy(j)

	dU3dt_d(i,j) =  -(E3_inv_d_right(i,j)-E3_inv_d_left(i,j))/ddx(i-1) &
                     -(F3_inv_d_upper(i,j)-F3_inv_d_lower(i,j))/ddy(j)

	dU4dt_d(i,j) =  -(E4_inv_d_right(i,j)-E4_inv_d_left(i,j))/ddx(i-1) &
                     -(F4_inv_d_upper(i,j)-F4_inv_d_lower(i,j))/ddy(j) 
Enddo
Enddo
!$OMP END DO

 ! Calculate the average derivative of solution flux (Um)
 !$OMP DO
 do j =2,Jmax-1
 do i =2,Imax-1
 dU1dt(i,j)= 0.25*(dU1dt_a(i,j)+dU1dt_b(i,j)+ &
                       dU1dt_c(i,j)+dU1dt_d(i,j))

 dU2dt(i,j)= 0.25*(dU2dt_a(i,j)+dU2dt_b(i,j)+ &
                      dU2dt_c(i,j)+dU2dt_d(i,j))

 dU3dt(i,j)= 0.25*(dU3dt_a(i,j)+dU3dt_b(i,j)+ &
                      dU3dt_c(i,j)+dU3dt_d(i,j))

 dU4dt(i,j)= 0.25*(dU4dt_a(i,j)+dU4dt_b(i,j)+ &
                       dU4dt_c(i,j)+dU4dt_d(i,j))
Enddo
Enddo
!$OMP END DO
!$OMP END PARALLEL
End Subroutine Derivatives
