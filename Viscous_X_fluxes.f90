!Evaluate the viscous spatial fluxes on the left and right side of (a,b,c,d) cells
Subroutine Viscous_X_fluxes
use Variables
use omp_lib
implicit none

!$OMP PARALLEL PRIVATE (i,j)
!$OMP DO
do j =2,Jmax-1
do i =2,Imax-1

	! cell(a) left
	E1_vis_a_left(i,j)= 0.0

	E2_vis_a_left(i,j)=-Tauxx_a_left(i,j)/Re_L

	E3_vis_a_left(i,j)= -Tauxy_a_left(i,j)/Re_L

	E4_vis_a_left(i,j)= - (Tauxx_a_left(i,j)*0.5*(u_old(i-1,j)+ &
     u_old(i-1,j-1))/Re_L+Tauxy_a_left(i,j)*0.5*(v_old(i-1,j)+ &
     v_old(i-1,j-1))/Re_L - qxx_a_left(i,j)/Re_L/pr_inf/cM_inf**2 &
          /(gamma-1))

       ! cell(a) right

	E1_vis_a_right(i,j)=  0.0

	E2_vis_a_right(i,j)=  Tauxx_a_right(i,j)/Re_L

	E3_vis_a_right(i,j)=  Tauxy_a_right(i,j)/Re_L

	E4_vis_a_right(i,j)=   Tauxx_a_right(i,j)*0.5*(u_old(i,j)+ &
     u_old(i,j-1))/Re_L+Tauxy_a_right(i,j)*0.5*(v_old(i,j)+ &
     v_old(i,j-1))/Re_L - qxx_a_right(i,j)/Re_L/pr_inf/cM_inf**2 &
          /(gamma-1)


	! cell(b) left
	E1_vis_b_left(i,j)= 0.0

	E2_vis_b_left(i,j)=-Tauxx_b_left(i,j)/Re_L

	E3_vis_b_left(i,j)= -Tauxy_b_left(i,j)/Re_L

	E4_vis_b_left(i,j)= - (Tauxx_b_left(i,j)*0.5*(u_old(i,j)+ &
     u_old(i,j-1))/Re_L+Tauxy_b_left(i,j)*0.5*(v_old(i,j)+ &
     v_old(i,j-1))/Re_L - qxx_b_left(i,j)/Re_L/pr_inf/cM_inf**2 &
          /(gamma-1))
       ! cell(b) right

	E1_vis_b_right(i,j)=  0.0

	E2_vis_b_right(i,j)=  Tauxx_b_right(i,j)/Re_L

	E3_vis_b_right(i,j)=  Tauxy_b_right(i,j)/Re_L

	E4_vis_b_right(i,j)=   Tauxx_b_right(i,j)*0.5*(u_old(i+1,j)+ &
     u_old(i+1,j-1))/Re_L+Tauxy_b_right(i,j)*0.5*(v_old(i+1,j)+ &
     v_old(i+1,j-1))/Re_L - qxx_b_right(i,j)/Re_L/pr_inf/cM_inf**2 &
          /(gamma-1)

 enddo
 enddo
 !$OMP END DO

 !$OMP DO
 do j =2,Jmax-1
 do i =2,Imax-1
	! cell(c) left
	E1_vis_c_left(i,j)= 0.0

	E2_vis_c_left(i,j)=-Tauxx_c_left(i,j)/Re_L

	E3_vis_c_left(i,j)= -Tauxy_c_left(i,j)/Re_L

	E4_vis_c_left(i,j)= - (Tauxx_c_left(i,j)*0.5*(u_old(i,j)+ &
     u_old(i,j+1))/Re_L+Tauxy_c_left(i,j)*0.5*(v_old(i,j)+ &
     v_old(i,j+1))/Re_L - qxx_c_left(i,j)/Re_L/pr_inf/cM_inf**2 &
          /(gamma-1))
       ! cell(c) right

	E1_vis_c_right(i,j)=  0.0

	E2_vis_c_right(i,j)=  Tauxx_c_right(i,j)/Re_L

	E3_vis_c_right(i,j)=  Tauxy_c_right(i,j)/Re_L

	E4_vis_c_right(i,j)=   Tauxx_c_right(i,j)*0.5*(u_old(i+1,j)+ &
     u_old(i+1,j+1))/Re_L+Tauxy_c_right(i,j)*0.5*(v_old(i+1,j)+ &
     v_old(i+1,j+1))/Re_L - qxx_c_right(i,j)/Re_L/pr_inf/cM_inf**2 &
          /(gamma-1)
	! cell(d) left
	E1_vis_d_left(i,j)= 0.0

	E2_vis_d_left(i,j)=-Tauxx_d_left(i,j)/Re_L

	E3_vis_d_left(i,j)= -Tauxy_d_left(i,j)/Re_L

	E4_vis_d_left(i,j)= - (Tauxx_d_left(i,j)*0.5*(u_old(i-1,j)+ &
     u_old(i-1,j+1))/Re_L+Tauxy_d_left(i,j)*0.5*(v_old(i-1,j)+ &
     v_old(i-1,j+1))/Re_L - qxx_d_left(i,j)/Re_L/pr_inf/cM_inf**2 &
          /(gamma-1))
       ! cell(d) right

	E1_vis_d_right(i,j)=  0.0

	E2_vis_d_right(i,j)=  Tauxx_d_right(i,j)/Re_L

	E3_vis_d_right(i,j)=  Tauxy_d_right(i,j)/Re_L

	E4_vis_d_right(i,j)=   Tauxx_d_right(i,j)*0.5*(u_old(i,j)+ &
     u_old(i,j+1))/Re_L+Tauxy_d_right(i,j)*0.5*(v_old(i,j)+ &
     v_old(i,j+1))/Re_L - qxx_d_right(i,j)/Re_L/pr_inf/cM_inf**2 &
          /(gamma-1)
Enddo
Enddo
!$OMP END DO
!$OMP END PARALLEL
End Subroutine Viscous_X_fluxes
