!Evaluate the viscous spatial fluxes on the lower and upper side of (a,b,c,d) cells
Subroutine Viscous_Y_fluxes
use Variables
use omp_lib
implicit none

!$OMP PARALLEL PRIVATE (i,j)
!$OMP DO
do j =2,Jmax-1
do i =2,Imax-1
	! cell(a) lower
	F1_vis_a_lower(i,j)= 0.0

	F2_vis_a_lower(i,j)=-Tauxy_a_lower(i,j)/Re_L

	F3_vis_a_lower(i,j)= -Tauyy_a_lower(i,j)/Re_L

	F4_vis_a_lower(i,j)= - (Tauyy_a_lower(i,j)*0.5*(v_old(i,j-1)+ &
     v_old(i-1,j-1))/Re_L+Tauxy_a_lower(i,j)*0.5*(u_old(i,j-1)+ &
     u_old(i-1,j-1))/Re_L - qyy_a_lower(i,j)/Re_L/pr_inf/cM_inf**2 &
          /(gamma-1))
       ! cell(a) upper

	F1_vis_a_upper(i,j)=  0.0

	F2_vis_a_upper(i,j)=  Tauxy_a_upper(i,j)/Re_L

	F3_vis_a_upper(i,j)=  Tauyy_a_upper(i,j)/Re_L

	F4_vis_a_upper(i,j)=   Tauyy_a_upper(i,j)*0.5*(v_old(i,j)+ &
     v_old(i-1,j))/Re_L+Tauxy_a_upper(i,j)*0.5*(u_old(i,j)+ &
     u_old(i-1,j))/Re_L - qyy_a_upper(i,j)/Re_L/pr_inf/cM_inf**2 &
          /(gamma-1)

	! cell(b) lower
	F1_vis_b_lower(i,j)= 0.0

	F2_vis_b_lower(i,j)=-Tauxy_b_lower(i,j)/Re_L

	F3_vis_b_lower(i,j)= -Tauyy_b_lower(i,j)/Re_L

	F4_vis_b_lower(i,j)= - (Tauyy_b_lower(i,j)*0.5*(v_old(i,j-1)+ &
     v_old(i+1,j-1))/Re_L+Tauxy_b_lower(i,j)*0.5*(u_old(i,j-1)+ &
     u_old(i+1,j-1))/Re_L - qyy_b_lower(i,j)/Re_L/pr_inf/cM_inf**2 &
          /(gamma-1))
       ! cell(b) upper

	F1_vis_b_upper(i,j)=  0.0

	F2_vis_b_upper(i,j)=  Tauxy_b_upper(i,j)/Re_L

	F3_vis_b_upper(i,j)=  Tauyy_b_upper(i,j)/Re_L

	F4_vis_b_upper(i,j)=   Tauyy_b_upper(i,j)*0.5*(v_old(i,j)+ &
     v_old(i+1,j))/Re_L+Tauxy_b_upper(i,j)*0.5*(u_old(i,j)+ &
     u_old(i+1,j))/Re_L - qyy_b_upper(i,j)/Re_L/pr_inf/cM_inf**2 &
          /(gamma-1)

 enddo
 enddo
 !$OMP END DO

 !$OMP DO
 do j =2,Jmax-1
 do i =2,Imax-1

	! cell(c) lower
	F1_vis_c_lower(i,j)= 0.0

	F2_vis_c_lower(i,j)=-Tauxy_c_lower(i,j)/Re_L

	F3_vis_c_lower(i,j)= -Tauyy_c_lower(i,j)/Re_L

	F4_vis_c_lower(i,j)= - (Tauyy_c_lower(i,j)*0.5*(v_old(i,j)+ &
     v_old(i+1,j))/Re_L+Tauxy_c_lower(i,j)*0.5*(u_old(i,j)+ &
     u_old(i+1,j))/Re_L - qyy_c_lower(i,j)/Re_L/pr_inf/cM_inf**2 &
          /(gamma-1))
       ! cell(c) upper

	F1_vis_c_upper(i,j)=  0.0

	F2_vis_c_upper(i,j)=  Tauxy_c_upper(i,j)/Re_L

	F3_vis_c_upper(i,j)=  Tauyy_c_upper(i,j)/Re_L

	F4_vis_c_upper(i,j)=   Tauyy_c_upper(i,j)*0.5*(v_old(i,j+1)+ &
     v_old(i+1,j+1))/Re_L+Tauxy_c_upper(i,j)*0.5*(u_old(i,j+1)+ &
     u_old(i+1,j+1))/Re_L - qyy_c_upper(i,j)/Re_L/pr_inf/cM_inf**2 &
          /(gamma-1)

	! cell(d) lower
	F1_vis_d_lower(i,j)= 0.0

	F2_vis_d_lower(i,j)=-Tauxy_d_lower(i,j)/Re_L

	F3_vis_d_lower(i,j)= -Tauyy_d_lower(i,j)/Re_L

	F4_vis_d_lower(i,j)= - (Tauyy_d_lower(i,j)*0.5*(v_old(i,j)+ &
     v_old(i-1,j))/Re_L+Tauxy_d_lower(i,j)*0.5*(u_old(i,j)+ &
     u_old(i-1,j))/Re_L - qyy_d_lower(i,j)/Re_L/pr_inf/cM_inf**2 &
          /(gamma-1))
       ! cell(d) upper

	F1_vis_d_upper(i,j)=  0.0

	F2_vis_d_upper(i,j)=  Tauxy_d_upper(i,j)/Re_L

	F3_vis_d_upper(i,j)=  Tauyy_d_upper(i,j)/Re_L

	F4_vis_d_upper(i,j)=   Tauyy_d_upper(i,j)*0.5*(v_old(i,j+1)+ &
     v_old(i-1,j+1))/Re_L+Tauxy_d_upper(i,j)*0.5*(u_old(i,j+1)+ &
     u_old(i-1,j+1))/Re_L - qyy_d_upper(i,j)/Re_L/pr_inf/cM_inf**2 &
          /(gamma-1)
Enddo
Enddo
!$OMP END DO
!$OMP END PARALLEL 
End Subroutine Viscous_Y_fluxes
