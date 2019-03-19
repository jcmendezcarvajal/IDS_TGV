!Evaluate the invicid spatial fluxes on the left and right side of (a,b,c,d) cells
Subroutine Invicid_X_fluxes
	use Variables
	use omp_lib
	implicit none

	!$OMP PARALLEL PRIVATE (i,j)
	!$OMP DO
	do j =2,Jmax-1
		do i =2,Imax-1

			

			! cell(a) left                                              ! average mass into cell (a)
			E1_inv_a_left(i,j)=  0.5*(r_old(i-1,j)*u_old(i-1,j) &     ! mass  @ j = j
			+r_old(i-1,j-1)*u_old(i-1,j-1))     ! mass  @ j = j-1

			! average x_mom & pressure into cell (a)
			E2_inv_a_left(i,j)=  0.5* (r_old(i-1,j)*u_old(i-1,j)**2 & ! x_mom  @ j = j
			+ r_old(i-1,j)*T_old(i-1,j)/gamma/cM_inf**2 &    ! pressure  @ j = j

			+ r_old(i-1,j-1)*u_old(i-1,j-1)**2 &     ! x_mom  @ j = j-1
			+ r_old(i-1,j-1)*T_old(i-1,j-1)/gamma/cM_inf**2) ! pressure  @ j = j-1

			! average y_mom  into cell (a)
			E3_inv_a_left(i,j)=  0.5* (r_old(i-1,j)*u_old(i-1,j)*v_old(i-1,j) & !y_mom @j
			+ r_old(i-1,j-1)*u_old(i-1,j-1)*v_old(i-1,j-1))          !y_mom @j-1


			! average Energy & pressure into cell (a)
			E4_inv_a_left(i,j)=  0.5* (H_old(i-1,j)*u_old(i-1,j) & ! energy @ j = j
			+ r_old(i-1,j)*T_old(i-1,j)*u_old(i-1,j)/gamma/cM_inf**2 & ! press @ j = j

			+ H_old(i-1,j-1)*u_old(i-1,j-1)  &  ! energy @ j = j-1
			+ r_old(i-1,j-1)*T_old(i-1,j-1)*u_old(i-1,j-1)/gamma/cM_inf**2)!pres @ j-1

			! cell(a) right
			! average mass outof cell (a)
			E1_inv_a_right(i,j)=  0.5*(r_old(i,j)*u_old(i,j)  &   ! mass  @ j = j
			+r_old(i,j-1)*u_old(i,j-1))     ! mass  @ j = j-1

			! average x_mom & pressure outof cell (a)
			E2_inv_a_right(i,j)=  0.5* (r_old(i,j)*u_old(i,j)**2 & ! x_mom  @ j = j
			+ r_old(i,j)*T_old(i,j)/gamma/cM_inf**2  &   ! pressure  @ j = j

			+ r_old(i,j-1)*u_old(i,j-1)**2  &    ! x_mom  @ j = j-1
			+ r_old(i,j-1)*T_old(i,j-1)/gamma/cM_inf**2) ! pressure  @ j = j-1

			! average y_mom  into cell (a)
			E3_inv_a_right(i,j)=  0.5* (r_old(i,j)*u_old(i,j)*v_old(i,j) & !y_mom @j
			+ r_old(i,j-1)*u_old(i,j-1)*v_old(i,j-1))          !y_mom @j-1


			! average Energy & pressure outof cell (a)
			E4_inv_a_right(i,j)=  0.5* (H_old(i,j)*u_old(i,j) & ! energy @ j = j
			+ r_old(i,j)*T_old(i,j)*u_old(i,j)/gamma/cM_inf**2 & ! press @ j = j

			+ H_old(i,j-1)*u_old(i,j-1)  &  ! energy @ j = j-1
			+ r_old(i,j-1)*T_old(i,j-1)*u_old(i,j-1)/gamma/cM_inf**2)!pres @ j-1
			! *****************************************
			! cell(b) left                                              ! average mass into cell (b)
			E1_inv_b_left(i,j)=  0.5*(r_old(i,j)*u_old(i,j)  &   ! mass  @ j = j
			+r_old(i,j-1)*u_old(i,j-1))     ! mass  @ j = j-1

			! average x_mom & pressure into cell (b)
			E2_inv_b_left(i,j)=  0.5* (r_old(i,j)*u_old(i,j)**2 & ! x_mom  @ j = j
			+ r_old(i,j)*T_old(i,j)/gamma/cM_inf**2   &  ! pressure  @ j = j

			+ r_old(i,j-1)*u_old(i,j-1)**2  &    ! x_mom  @ j = j-1
			+ r_old(i,j-1)*T_old(i,j-1)/gamma/cM_inf**2) ! pressure  @ j = j-1

			! average y_mom  into cell (b)
			E3_inv_b_left(i,j)=  0.5* (r_old(i,j)*u_old(i,j)*v_old(i,j) & !y_mom @j
			+ r_old(i,j-1)*u_old(i,j-1)*v_old(i,j-1))          !y_mom @j-1


			! average Energy & pressure into cell (b)
			E4_inv_b_left(i,j)=  0.5* (H_old(i,j)*u_old(i,j) & ! energy @ j = j
			+ r_old(i,j)*T_old(i,j)*u_old(i,j)/gamma/cM_inf**2 & ! press @ j = j

			+ H_old(i,j-1)*u_old(i,j-1)  &  ! energy @ j = j-1
			+ r_old(i,j-1)*T_old(i,j-1)*u_old(i,j-1)/gamma/cM_inf**2)!pres @ j-1

			! cell(b) right
			! average mass outof cell (b)
			E1_inv_b_right(i,j)=  0.5*(r_old(i+1,j)*u_old(i+1,j) &    ! mass  @ j = j
			+r_old(i+1,j-1)*u_old(i+1,j-1))     ! mass  @ j = j-1

			! average x_mom & pressure outof cell (b)
			E2_inv_b_right(i,j)=  0.5* (r_old(i+1,j)*u_old(i+1,j)**2 & ! x_mom  @ j = j
			+ r_old(i+1,j)*T_old(i+1,j)/gamma/cM_inf**2 &    ! pressure  @ j = j

			+ r_old(i+1,j-1)*u_old(i+1,j-1)**2  &    ! x_mom  @ j = j-1
			+ r_old(i+1,j-1)*T_old(i+1,j-1)/gamma/cM_inf**2) ! pressure  @ j = j-1

			! average y_mom  outof cell (b)
			E3_inv_b_right(i,j)= 0.5* (r_old(i+1,j)*u_old(i+1,j)*v_old(i+1,j) & !y_mom @j
			+ r_old(i+1,j-1)*u_old(i+1,j-1)*v_old(i+1,j-1))          !y_mom @j-1


			! average Energy & pressure outof cell (b)
			E4_inv_b_right(i,j)=  0.5* (H_old(i+1,j)*u_old(i+1,j) & ! energy @ j = j
			+ r_old(i+1,j)*T_old(i+1,j)*u_old(i+1,j)/gamma/cM_inf**2 & ! press @ j = j

			+ H_old(i+1,j-1)*u_old(i+1,j-1)  &  ! energy @ j = j-1
			+r_old(i+1,j-1)*T_old(i+1,j-1)*u_old(i+1,j-1)/gamma/cM_inf**2)!pres @ j-1

		Enddo
	Enddo
	!$OMP END DO


	! *****************************************
	!$OMP DO
	do j =2,Jmax-1
		do i =2,Imax-1
			! cell(c) left                                              ! average mass into cell (c)
			E1_inv_c_left(i,j)=  0.5*(r_old(i,j)*u_old(i,j)  &   ! mass  @ j = j
			+r_old(i,j+1)*u_old(i,j+1))     ! mass  @ j = j+1

			! average x_mom & pressure into cell (c)
			E2_inv_c_left(i,j)=  0.5* (r_old(i,j)*u_old(i,j)**2 & ! x_mom  @ j = j
			+ r_old(i,j)*T_old(i,j)/gamma/cM_inf**2  &   ! pressure  @ j = j

			+ r_old(i,j+1)*u_old(i,j+1)**2  &    ! x_mom  @ j = j+1
			+ r_old(i,j+1)*T_old(i,j+1)/gamma/cM_inf**2) ! pressure  @ j = j+1

			! average y_mom  into cell (c)
			E3_inv_c_left(i,j)=  0.5* (r_old(i,j)*u_old(i,j)*v_old(i,j) & !y_mom @j
			+ r_old(i,j+1)*u_old(i,j+1)*v_old(i,j+1))          !y_mom @j+1


			! average Energy & pressure into cell (c)
			E4_inv_c_left(i,j)=  0.5* (H_old(i,j)*u_old(i,j) & ! energy @ j = j
			+ r_old(i,j)*T_old(i,j)*u_old(i,j)/gamma/cM_inf**2 & ! press @ j = j

			+ H_old(i,j+1)*u_old(i,j+1) &   ! energy @ j = j+1
			+ r_old(i,j+1)*T_old(i,j+1)*u_old(i,j+1)/gamma/cM_inf**2)!pres @ j+1

			! cell(c) right
			! average mass outof cell (c)
			E1_inv_c_right(i,j)=  0.5*(r_old(i+1,j)*u_old(i+1,j) &    ! mass  @ j = j
			+r_old(i+1,j+1)*u_old(i+1,j+1))     ! mass  @ j = j+1

			! average x_mom & pressure outof cell (c)
			E2_inv_c_right(i,j)=  0.5* (r_old(i+1,j)*u_old(i+1,j)**2 & ! x_mom  @ j = j
			+ r_old(i+1,j)*T_old(i+1,j)/gamma/cM_inf**2  &   ! pressure  @ j = j

			+ r_old(i+1,j+1)*u_old(i+1,j+1)**2  &    ! x_mom  @ j = j+1
			+ r_old(i+1,j+1)*T_old(i+1,j+1)/gamma/cM_inf**2) ! pressure  @ j = j+1

			! average y_mom  outof cell (c)
			E3_inv_c_right(i,j)= 0.5* (r_old(i+1,j)*u_old(i+1,j)*v_old(i+1,j) & !y_mom @j
			+ r_old(i+1,j+1)*u_old(i+1,j+1)*v_old(i+1,j+1))          !y_mom @j+1


			! average Energy & pressure outof cell (c)
			E4_inv_c_right(i,j)=  0.5* (H_old(i+1,j)*u_old(i+1,j) & ! energy @ j = j
			+ r_old(i+1,j)*T_old(i+1,j)*u_old(i+1,j)/gamma/cM_inf**2 & ! press @ j = j

			+ H_old(i+1,j+1)*u_old(i+1,j+1)   & ! energy @ j = j+1
			+r_old(i+1,j+1)*T_old(i+1,j+1)*u_old(i+1,j+1)/gamma/cM_inf**2)!pres @ j+1

			! *****************************************
			! cell(d) left                                              ! average mass into cell (d)
			E1_inv_d_left(i,j)=  0.5*(r_old(i-1,j)*u_old(i-1,j) &    ! mass  @ j = j
			+r_old(i-1,j+1)*u_old(i-1,j+1))     ! mass  @ j = j+1

			! average x_mom & pressure into cell (d)
			E2_inv_d_left(i,j)=  0.5* (r_old(i-1,j)*u_old(i-1,j)**2 & ! x_mom  @ j = j
			+ r_old(i-1,j)*T_old(i-1,j)/gamma/cM_inf**2  &   ! pressure  @ j = j

			+ r_old(i-1,j+1)*u_old(i-1,j+1)**2  &    ! x_mom  @ j = j+1
			+ r_old(i-1,j+1)*T_old(i-1,j+1)/gamma/cM_inf**2) ! pressure  @ j = j+1

			! average y_mom  into cell (d)
			E3_inv_d_left(i,j)=  0.5* (r_old(i-1,j)*u_old(i-1,j)*v_old(i-1,j) & !y_mom @j
			+ r_old(i-1,j+1)*u_old(i-1,j+1)*v_old(i-1,j+1))          !y_mom @j+1


			! average Energy & pressure into cell (d)
			E4_inv_d_left(i,j)=  0.5* (H_old(i-1,j)*u_old(i-1,j) & ! energy @ j = j
			+ r_old(i-1,j)*T_old(i-1,j)*u_old(i-1,j)/gamma/cM_inf**2 & ! press @ j = j

			+ H_old(i-1,j+1)*u_old(i-1,j+1)  &  ! energy @ j = j+1
			+ r_old(i-1,j+1)*T_old(i-1,j+1)*u_old(i-1,j+1)/gamma/cM_inf**2)!pres @ j+1

			! cell(d) right
			! average mass outof cell (d)
			E1_inv_d_right(i,j)=  0.5*(r_old(i,j)*u_old(i,j)  &   ! mass  @ j = j
			+r_old(i,j+1)*u_old(i,j+1))     ! mass  @ j = j+1

			! average x_mom & pressure outof cell (d)
			E2_inv_d_right(i,j)=  0.5* (r_old(i,j)*u_old(i,j)**2 & ! x_mom  @ j = j
			+ r_old(i,j)*T_old(i,j)/gamma/cM_inf**2 &   ! pressure  @ j = j

			+ r_old(i,j+1)*u_old(i,j+1)**2  &    ! x_mom  @ j = j+1
			+ r_old(i,j+1)*T_old(i,j+1)/gamma/cM_inf**2) ! pressure  @ j = j+1

			! average y_mom  outof cell (d)
			E3_inv_d_right(i,j)= 0.5* (r_old(i,j)*u_old(i,j)*v_old(i,j) & !y_mom @j
			+ r_old(i,j+1)*u_old(i,j+1)*v_old(i,j+1))          !y_mom @j+1


			! average Energy & pressure outof cell (d)
			E4_inv_d_right(i,j)=  0.5* (H_old(i,j)*u_old(i,j) & ! energy @ j = j
			+ r_old(i,j)*T_old(i,j)*u_old(i,j)/gamma/cM_inf**2 & ! press @ j = j

			+ H_old(i,j+1)*u_old(i,j+1)  &  ! energy @ j = j+1
			+r_old(i,j+1)*T_old(i,j+1)*u_old(i,j+1)/gamma/cM_inf**2)!pres @ j+1
		Enddo
	Enddo
	!$OMP END DO
	!$OMP END PARALLEL
End Subroutine Invicid_X_fluxes
