!Evaluate the invicid spatial fluxes on the lower and upper side of (a,b,c,d) cells
Subroutine Invicid_Y_fluxes
	use Variables
	use omp_lib
	implicit none

	!$OMP PARALLEL PRIVATE (i,j)
	!$OMP DO
	do j =2,Jmax-1
		do i =2,Imax-1

			! cell(a) lower  ***********                               ! average mass into cell (a)
			F1_inv_a_lower(i,j)=  0.5*(r_old(i,j-1)*v_old(i,j-1) &     ! mass  @ i = i
			+r_old(i-1,j-1)*v_old(i-1,j-1))     ! mass  @ i = i-1


			! average x_mom  into cell (a)
			F2_inv_a_lower(i,j)= 0.5* (r_old(i,j-1)*u_old(i,j-1)*v_old(i,j-1) & !x_mom @i
			+ r_old(i-1,j-1)*u_old(i-1,j-1)*v_old(i-1,j-1))          !x_mom @i-1

			! average y_mom & pressure into cell (a)
			F3_inv_a_lower(i,j)=  0.5* (r_old(i,j-1)*v_old(i,j-1)**2 & ! y_mom  @ i = i
			+ r_old(i,j-1)*T_old(i,j-1)/gamma/cM_inf**2  &  ! pressure  @ i = i

			+ r_old(i-1,j-1)*v_old(i-1,j-1)**2  &    ! y_mom  @ i = i-1
			+ r_old(i-1,j-1)*T_old(i-1,j-1)/gamma/cM_inf**2) ! pressure  @ i = i-1


			! average Energy & pressure into cell (a)
			F4_inv_a_lower(i,j)=  0.5* (H_old(i,j-1)*v_old(i,j-1) & ! energy @ i = i
			+ r_old(i,j-1)*T_old(i,j-1)*v_old(i,j-1)/gamma/cM_inf**2 & ! press @ i = i

			+ H_old(i-1,j-1)*v_old(i-1,j-1) &   ! energy @ i = i-1
			+ r_old(i-1,j-1)*T_old(i-1,j-1)*v_old(i-1,j-1)/gamma/cM_inf**2)!pres @ i-1

			! cell(a) upper
			! average mass outof cell (a)
			F1_inv_a_upper(i,j)=  0.5*(r_old(i,j)*v_old(i,j)  &   ! mass  @ i = i
			+r_old(i-1,j)*v_old(i-1,j))     ! mass  @ i = i-1

			! average x_mom  into cell (a)
			F2_inv_a_upper(i,j)=  0.5* (r_old(i,j)*u_old(i,j)*v_old(i,j) & !x_mom @i
			+ r_old(i-1,j)*u_old(i-1,j)*v_old(i-1,j))          !x_mom @i-1

			! average x_mom & pressure outof cell (a)
			F3_inv_a_upper(i,j)=  0.5* (r_old(i,j)*v_old(i,j)**2 & ! x_mom  @ i = i
			+ r_old(i,j)*T_old(i,j)/gamma/cM_inf**2  &   ! pressure  @ i = i

			+ r_old(i-1,j)*v_old(i-1,j)**2  &    ! x_mom  @ i = i-1
			+ r_old(i-1,j)*T_old(i-1,j)/gamma/cM_inf**2) ! pressure  @ i = i-1


			! average Energy & pressure outof cell (a)
			F4_inv_a_upper(i,j)=  0.5* (H_old(i,j)*v_old(i,j) & ! energy @ i = i
			+ r_old(i,j)*T_old(i,j)*v_old(i,j)/gamma/cM_inf**2 & ! press @ i = i

			+ H_old(i-1,j)*v_old(i-1,j)  &  ! energy @ i = i-1
			+ r_old(i-1,j)*T_old(i-1,j)*v_old(i-1,j)/gamma/cM_inf**2)!pres @ i-1

			! cell(b) lower ****************                              ! average mass into cell (b)
			F1_inv_b_lower(i,j)=  0.5*(r_old(i,j-1)*v_old(i,j-1)  &   ! mass  @ i = i
			+r_old(i+1,j-1)*v_old(i+1,j-1))     ! mass  @ i = i+1


			! average x_mom  into cell (b)
			F2_inv_b_lower(i,j)= 0.5* (r_old(i,j-1)*u_old(i,j-1)*v_old(i,j-1) & !x_mom @i
			+ r_old(i+1,j-1)*u_old(i+1,j-1)*v_old(i+1,j-1))          !x_mom @i+1

			! average y_mom & pressure into cell (b)
			F3_inv_b_lower(i,j)=  0.5* (r_old(i,j-1)*v_old(i,j-1)**2 & ! y_mom  @ i = i
			+ r_old(i,j-1)*T_old(i,j-1)/gamma/cM_inf**2   & ! pressure  @ i = i

			+ r_old(i+1,j-1)*v_old(i+1,j-1)**2   &   ! y_mom  @ i = i+1
			+ r_old(i+1,j-1)*T_old(i+1,j-1)/gamma/cM_inf**2) ! pressure  @ i = i+1


			! average Energy & pressure into cell (b)
			F4_inv_b_lower(i,j)=  0.5* (H_old(i,j-1)*v_old(i,j-1) & ! energy @ i = i
			+ r_old(i,j-1)*T_old(i,j-1)*v_old(i,j-1)/gamma/cM_inf**2 & ! press @ i = i

			+ H_old(i+1,j-1)*v_old(i+1,j-1)  &  ! energy @ i = i+1
			+r_old(i+1,j-1)*T_old(i+1,j-1)*v_old(i+1,j-1)/gamma/cM_inf**2)!pres @ i+1

			! cell(b) upper
			! average mass outof cell (b)
			F1_inv_b_upper(i,j)=  0.5*(r_old(i,j)*v_old(i,j)  &   ! mass  @ i = i
			+r_old(i+1,j)*v_old(i+1,j))     ! mass  @ i = i+1

			! average x_mom  into cell (b)
			F2_inv_b_upper(i,j)=  0.5* (r_old(i,j)*u_old(i,j)*v_old(i,j) & !x_mom @i
			+ r_old(i+1,j)*u_old(i+1,j)*v_old(i+1,j))          !x_mom @i+1

			! average x_mom & pressure outof cell (b)
			F3_inv_b_upper(i,j)=  0.5* (r_old(i,j)*v_old(i,j)**2 & ! x_mom  @ i = i
			+ r_old(i,j)*T_old(i,j)/gamma/cM_inf**2  &   ! pressure  @ i = i

			+ r_old(i+1,j)*v_old(i+1,j)**2   &   ! x_mom  @ i = i+1
			+ r_old(i+1,j)*T_old(i+1,j)/gamma/cM_inf**2) ! pressure  @ i = i+1



			! average Energy & pressure outof cell (b)
			F4_inv_b_upper(i,j)=  0.5* (H_old(i,j)*v_old(i,j) & ! energy @ i = i
			+ r_old(i,j)*T_old(i,j)*v_old(i,j)/gamma/cM_inf**2 & ! press @ i = i

			+ H_old(i+1,j)*v_old(i+1,j) &   ! energy @ i = i+1
			+ r_old(i+1,j)*T_old(i+1,j)*v_old(i+1,j)/gamma/cM_inf**2)!pres @ i+1

		Enddo
	Enddo
	!$OMP END DO



	!$OMP DO
	do j =2,Jmax-1
		do i =2,Imax-1
			! cell(c) lower ****************                              ! average mass into cell (c)
			F1_inv_c_lower(i,j)=  0.5*(r_old(i,j)*v_old(i,j) &    ! mass  @ i = i
			+r_old(i+1,j)*v_old(i+1,j))     ! mass  @ i = i+1


			! average x_mom  into cell (c)
			F2_inv_c_lower(i,j)= 0.5* (r_old(i,j)*u_old(i,j)*v_old(i,j) & !x_mom @i
			+ r_old(i+1,j)*u_old(i+1,j)*v_old(i+1,j))          !x_mom @i+1

			! average y_mom & pressure into cell (c)
			F3_inv_c_lower(i,j)=  0.5* (r_old(i,j)*v_old(i,j)**2  & ! y_mom  @ i = i
			+ r_old(i,j)*T_old(i,j)/gamma/cM_inf**2  &  ! pressure  @ i = i

			+ r_old(i+1,j)*v_old(i+1,j)**2   &   ! y_mom  @ i = i+1
			+ r_old(i+1,j)*T_old(i+1,j)/gamma/cM_inf**2) ! pressure  @ i = i+1


			! average Energy & pressure into cell (c)
			F4_inv_c_lower(i,j)=  0.5* (H_old(i,j)*v_old(i,j) & ! energy @ i = i
			+ r_old(i,j)*T_old(i,j)*v_old(i,j)/gamma/cM_inf**2 & ! press @ i = i

			+ H_old(i+1,j)*v_old(i+1,j)  &  ! energy @ i = i+1
			+r_old(i+1,j)*T_old(i+1,j)*v_old(i+1,j)/gamma/cM_inf**2)!pres @ i+1

			! cell(c) upper
			! average mass outof cell (c)
			F1_inv_c_upper(i,j)=  0.5*(r_old(i,j+1)*v_old(i,j+1) &    ! mass  @ i = i
			+r_old(i+1,j+1)*v_old(i+1,j+1))     ! mass  @ i = i+1

			! average x_mom  into cell (c)
			F2_inv_c_upper(i,j)= 0.5* (r_old(i,j+1)*u_old(i,j+1)*v_old(i,j+1) & !x_mom @i
			+ r_old(i+1,j+1)*u_old(i+1,j+1)*v_old(i+1,j+1))          !x_mom @i+1

			! average x_mom & pressure outof cell (c)
			F3_inv_c_upper(i,j)=  0.5* (r_old(i,j+1)*v_old(i,j+1)**2 & ! x_mom  @ i = i
			+ r_old(i,j+1)*T_old(i,j+1)/gamma/cM_inf**2 &    ! pressure  @ i = i

			+ r_old(i+1,j+1)*v_old(i+1,j+1)**2  &    ! x_mom  @ i = i+1
			+ r_old(i+1,j+1)*T_old(i+1,j+1)/gamma/cM_inf**2) ! pressure  @ i = i+1



			! average Energy & pressure outof cell (c)
			F4_inv_c_upper(i,j)=  0.5* (H_old(i,j+1)*v_old(i,j+1) & ! energy @ i = i
			+ r_old(i,j+1)*T_old(i,j+1)*v_old(i,j+1)/gamma/cM_inf**2 & ! press @ i = i

			+ H_old(i+1,j+1)*v_old(i+1,j+1) &   ! energy @ i = i+1
			+ r_old(i+1,j+1)*T_old(i+1,j+1)*v_old(i+1,j+1)/gamma/cM_inf**2)!pres @ i+1

			! cell(d) lower ****************                              ! average mass into cell (d)
			F1_inv_d_lower(i,j)=  0.5*(r_old(i,j)*v_old(i,j)  &   ! mass  @ i = i
			+r_old(i-1,j)*v_old(i-1,j))     ! mass  @ i = i-1


			! average x_mom  into cell (d)
			F2_inv_d_lower(i,j)= 0.5* (r_old(i,j)*u_old(i,j)*v_old(i,j) & !x_mom @i
			+ r_old(i-1,j)*u_old(i-1,j)*v_old(i-1,j))          !x_mom @i-1

			! average y_mom & pressure into cell (d)
			F3_inv_d_lower(i,j)=  0.5* (r_old(i,j)*v_old(i,j)**2 & ! y_mom  @ i = i
			+ r_old(i,j)*T_old(i,j)/gamma/cM_inf**2  &  ! pressure  @ i = i

			+ r_old(i-1,j)*v_old(i-1,j)**2 &     ! y_mom  @ i = i-1
			+ r_old(i-1,j)*T_old(i-1,j)/gamma/cM_inf**2) ! pressure  @ i = i-1


			! average Energy & pressure into cell (d)
			F4_inv_d_lower(i,j)=  0.5* (H_old(i,j)*v_old(i,j) & ! energy @ i = i
			+ r_old(i,j)*T_old(i,j)*v_old(i,j)/gamma/cM_inf**2 & ! press @ i = i

			+ H_old(i-1,j)*v_old(i-1,j)  &  ! energy @ i = i-1
			+r_old(i-1,j)*T_old(i-1,j)*v_old(i-1,j)/gamma/cM_inf**2)!pres @ i-1

			! cell(d) upper
			! average mass outof cell (d)
			F1_inv_d_upper(i,j)=  0.5*(r_old(i,j+1)*v_old(i,j+1)  &   ! mass  @ i = i
			+r_old(i-1,j+1)*v_old(i-1,j+1))     ! mass  @ i = i-1

			! average x_mom  into cell (d)
			F2_inv_d_upper(i,j)= 0.5* (r_old(i,j+1)*u_old(i,j+1)*v_old(i,j+1) & !x_mom @i
			+ r_old(i-1,j+1)*u_old(i-1,j+1)*v_old(i-1,j+1))          !x_mom @i-1

			! average x_mom & pressure outof cell (d)
			F3_inv_d_upper(i,j)=  0.5* (r_old(i,j+1)*v_old(i,j+1)**2  & ! x_mom  @ i = i
			+ r_old(i,j+1)*T_old(i,j+1)/gamma/cM_inf**2  &   ! pressure  @ i = i

			+ r_old(i-1,j+1)*v_old(i-1,j+1)**2 &     ! x_mom  @ i = i-1
			+ r_old(i-1,j+1)*T_old(i-1,j+1)/gamma/cM_inf**2) ! pressure  @ i = i-1



			! average Energy & pressure outof cell (d)
			F4_inv_d_upper(i,j)=  0.5* (H_old(i,j+1)*v_old(i,j+1) & ! energy @ i = i
			+ r_old(i,j+1)*T_old(i,j+1)*v_old(i,j+1)/gamma/cM_inf**2 & ! press @ i = i

			+ H_old(i-1,j+1)*v_old(i-1,j+1) &   ! energy @ i = i-1
			+ r_old(i-1,j+1)*T_old(i-1,j+1)*v_old(i-1,j+1)/gamma/cM_inf**2)!pres @ i-1
		Enddo
	Enddo
	!$OMP END DO
	!$OMP END PARALLEL
End Subroutine Invicid_Y_fluxes
