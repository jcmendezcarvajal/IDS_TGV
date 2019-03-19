! Calculate the stresses and heat transfer on the right side of (a,b,c,d) boundary cells
Subroutine R_vis_bc
use Variables
use omp_lib
implicit none
	do j =2,Jmax-1

	! cell (a)                                                 !Mu_avg @ Imax-1
	Tauxy_a_right(Imax-1,j)= 0.5*(vMu(Imax-1,j)+vMu(Imax-1,j-1)) &
															  !Avg forw v_Grad
         *(0.5 * (v_old(Imax,j)-v_old(Imax-1,j)  &              !forw_Grad @ j
                +v_old(Imax,j-1)-v_old(Imax-1,j-1)  )/ddx(Imax-1) &      !Forw_Grad @ j-1
                                               !u_Grad along Seg @ fixed i= Imax-1
     	+  (u_old(Imax-1,j)-u_old(Imax-1,j-1))/ddy(j-1))
                                                               !Mu_avg @ i = Imax-1
	Tauxx_a_right(Imax-1,j)= (1.0/3.0)*(vMu(Imax-1,j)+vMu(Imax-1,j-1))&
                                                              !Avg Forw u_Grad
         * ((u_old(Imax,j)-u_old(Imax-1,j)&                   !Forw_Grad @ j
               +u_old(Imax,j-1)-u_old(Imax-1,j-1)  )/ddx(Imax-1) &     !Forw_Grad @ j-1
                                              !v_Grad along Seg @ fixed i= Imax-1
     	-  (v_old(Imax-1,j)-v_old(Imax-1,j-1))/ddy(j-1))
                                                            !tK_avg @ i = Imax-1
	qxx_a_right(Imax-1,j)= -0.25*(tK(Imax-1,j)+tK(Imax-1,j-1)) &
                      * ((T_old(Imax,j)-T_old(Imax-1,j) &       !T_forw_Grad @ j
                     + T_old(Imax,j-1)-T_old(Imax-1,j-1))/ddx(Imax-1))!T_forw_Grad @ j-1

	! cell (b)                                                !Mu_avg @ i= Imax
	Tauxy_b_right(Imax-1,j)= 0.5*(vMu(Imax,j)+vMu(Imax,j-1)  ) &
													      !Avg back_v_Grad
         *(0.5 * (v_old(Imax,j)-v_old(Imax-1,j) &           !backw_Grad @ j
                +v_old(Imax,j-1)-v_old(Imax-1,j-1)  )/ddx(Imax-1) &  !backw_Grad @ j-1
                                               !u_Grad along Seg @ fixed i= Imax
     	+  (u_old(Imax,j)-u_old(Imax,j-1))/ddy(i-1))
                                                               !Mu_avg @ i= Imax
	Tauxx_b_right(Imax-1,j)= (1.0/3.0)*(vMu(Imax,j)+vMu(Imax,j-1))&
                                                              !Avg Forw u_Grad
         * ((u_old(Imax,j)-u_old(Imax-1,j)&                    !backw_Grad @ j
               +u_old(Imax,j-1)-u_old(Imax-1,j-1)  )/ddx(Imax-1) &      !backw_Grad @ j-1
                                              !v_Grad along Seg @ fixed i= i+1
     	-  (v_old(Imax,j)-v_old(Imax,j-1))/ddy(j-1))

	                                                          !tK_avg @ i= Imax
	qxx_b_right(Imax-1,j)= -0.25*(tK(Imax,j)+tK(Imax,j-1)) &
                      * ((T_old(Imax,j)-T_old(Imax-1,j) &     !T_backw_Grad @ j
                     +T_old(Imax,j-1)-T_old(Imax-1,j-1))/ddx(Imax-1))!T_backw_Grad @ j-1


	! cell (c)                                               !Mu_avg@ i= Imax

	Tauxy_c_right(Imax-1,j)= 0.5*(vMu(Imax,j)+vMu(Imax,j+1)  ) &
															!Avg Forw v_Grad
         *(0.5 * (v_old(Imax,j)-v_old(Imax-1,j) &             !backw_Grad @ j
                +v_old(Imax,j+1)-v_old(Imax-1,j+1)  )/ddx(Imax-1) &   !backw_Grad @ j+1
                                               !u_Grad along Seg @ fixed i= Imax
     	+  (u_old(Imax,j+1)-u_old(Imax,j))/ddy(j))
                                                                !Mu_avg@ i= Imax
	Tauxx_c_right(Imax-1,j)= (1.0/3.0)*(vMu(Imax,j)+vMu(Imax,j+1)) &
                                                              !Avg Forw u_Grad
         * ((u_old(Imax,j)-u_old(Imax-1,j)  &                 !back_Grad @ j
               +u_old(Imax,j+1)-u_old(Imax-1,j+1)  )/ddx(Imax-1)  &    !back_Grad @ j+1
                                              !v_Grad along Seg @ fixed i= Imax
     	-  (v_old(Imax,j+1)-v_old(Imax,j))/ddy(j))
                                                              !tK_avg @ i= Imax
	qxx_c_right(Imax-1,j)= -0.25*(tK(Imax,j)+tK(Imax,j+1)) &
                      * ((T_old(Imax,j)-T_old(Imax-1,j) &       !T_back_Grad @ j
                    +T_old(Imax,j+1)-T_old(Imax-1,j+1))/ddx(Imax-1))!T_back_Grad @ j+1


	! cell (d)                                              !Mu_avg @ Imax-1
	Tauxy_d_right(Imax-1,j)= 0.5*(vMu(Imax-1,j)+vMu(Imax-1,j+1)  ) &
															!Avg Forw v_Grad
         *(0.5 * (v_old(Imax,j)-v_old(Imax-1,j) &              !Forw_Grad @ j
                +v_old(Imax,j+1)-v_old(Imax-1,j+1)  )/ddx(Imax-1)  &     !Forw_Grad @ j+1
                                               !u_Grad along Seg @ fixed i= Imax-1
     	+  (u_old(Imax-1,j+1)-u_old(Imax-1,j))/ddy(j))
                                                               !Mu_avg @ Imax-1
	Tauxx_d_right(Imax-1,j)= (1.0/3.0)*(vMu(Imax-1,j)+vMu(Imax-1,j+1)) &
                                                              !Avg Forw u_Grad
         * ((u_old(Imax,j)-u_old(Imax-1,j)  &                 !Forw_Grad @ j=j
               +u_old(Imax,j+1)-u_old(Imax-1,j+1)  )/ddx(Imax-1) &     !Forw_Grad @ j=j+1
                                              !v_Grad along Seg @ fixed i= Imax-1
     	-  (v_old(Imax-1,j+1)-v_old(Imax-1,j))/ddy(j))
                                                                !tK_avg @ i= Imax-1
	qxx_d_right(Imax-1,j)= -0.25*(tK(Imax-1,j)+tK(Imax-1,j+1)) &
                      * ((T_old(Imax,j)-T_old(Imax-1,j) &     !T_forw_Grad @ j
                    + T_old(Imax,j+1)-T_old(Imax-1,j+1))/ddx(Imax-1))!T_forw_Grad @ j+1
	enddo
End Subroutine R_vis_bc
