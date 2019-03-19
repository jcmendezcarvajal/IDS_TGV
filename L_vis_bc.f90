! Calculate the stresses and heat transfer on the left side of (a,b,c,d) boundary cells

Subroutine L_vis_bc
use Variables
use omp_lib
implicit none

  do j =2,Jmax-1

! cell (a)
	Tauxy_a_left(2,j)   = 0.5*( vMu(1,j) + vMu(1,j-1) )&  !Mu_avg @ i = 1
															 !Avg Back v_Grad
         *(0.5 * (v_old(2,j)-v_old(1,j)&                   !Back_Grad @ j
     			+v_old(2,j-1)-v_old(1,j-1) )/ddx(1) &         !Back_Grad @ j-1
                                               !u_Grad along Seg @ fixed i= 1
     	+ (u_old(1,j)-u_old(1,j-1))/ddy(j-1))


	Tauxx_a_left(2,j)= (1.0/3.0)*(vMu(1,j)+vMu(1,j-1)) & !Mu_avg @  i=1
                                                              !Avg Back u_Grad
         * ((u_old(2,j)-u_old(1,j) &                      !Back_Grad @ j
               +u_old(2,j-1)-u_old(1,j-1)  )/ddx(1)   &       !Back_Grad @ j-1
	                                         !v_Grad along Seg @ fixed i= 1
     	-  (v_old(1,j)-v_old(1,j-1))/ddy(j-1))


	qxx_a_left(2,j)= -0.25*(tK(1,j)+tK(1,j-1)) &         !tK_avg @ i = 1
                      * ((T_old(2,j)-T_old(1,j)  &        !T_back_Grad @ j
                      +  T_old(2,j-1)-T_old(1,j-1))/ddx(1))  !T_back_Grad @ j-1

	! cell (b)
	Tauxy_b_left(2,j)   = 0.5*( vMu(2,j) + vMu(2,j-1) ) &    !Mu_avg @ i = 2
															!Avg Back v_Grad
         *(0.5 * (v_old(2,j)-v_old(1,j) &                   !Back_Grad @ j
     			+v_old(2,j-1)-v_old(1,j-1) )/ddx(2)  &         !Back_Grad @ j-1
                                               !u_Grad along Seg @ fixed i= 2
     	+ (u_old(2,j)-u_old(2,j-1))/ddy(j))

	Tauxx_b_left(2,j)= (1.0/3.0)*(vMu(2,j)+vMu(2,j)) &     !Mu_avg @  i=2
                                                              !Avg Back u_Grad
         * ((u_old(2,j)-u_old(1,j)  &                       !Back_Grad @ j
               +u_old(2,j-1)-u_old(1,j-1)  )/ddx(2)  &          !Back_Grad @ j-1
	                                         !v_Grad along Seg @ fixed i= 2
     	-  (v_old(2,j)-v_old(2,j-1))/ddy(j))

	qxx_b_left(2,j)= -0.25*(tK(2,j)+tK(2,j-1))  &            !tK_avg @ i = 2
                      * ((T_old(2,j)-T_old(1,j)   &         !T_back_Grad @ j
                      +  T_old(2,j-1)-T_old(1,j-1))/ddx(2))    !T_back_Grad @ j-1

	! cell (c)

	Tauxy_c_left(2,j)   = 0.5*( vMu(2,j+1) + vMu(2,j) ) &    !Mu_avg @ i = 2
															!Avg Back v_Grad
         *(0.5 * (v_old(2,j)-v_old(1,j)  &                  !Back_Grad @ j
     			+v_old(2,j+1)-v_old(1,j+1) )/ddx(2)  &         !Back_Grad @ j+1
                                               !u_Grad along Seg @ fixed i= 2
     	+ (u_old(2,j+1)-u_old(2,j))/ddy(j))

	Tauxx_c_left(2,j)= (1.0/3.0)*(vMu(2,j)+vMu(2,j+1)) &     !Mu_avg @  i=2
                                                              !Avg Back u_Grad
         * ((u_old(2,j)-u_old(1,j) &                        !Back_Grad @ j
               +u_old(2,j+1)-u_old(1,j+1)  )/ddx(2)   &         !Back_Grad @ j+1
	                                         !v_Grad along Seg @ fixed i= 2
     	-  (v_old(2,j+1)-v_old(2,j))/ddy(j))

	qxx_c_left(2,j)= -0.25*(tK(2,j)+tK(2,j+1))  &          !tK_avg @ i = 2
                      * ((T_old(2,j)-T_old(1,j)  &        !T_back_Grad @ j=j
                      +  T_old(2,j+1)-T_old(1,j+1))/ddx(2))  !T_back_Grad @ j=j+1

	! cell (d)

	Tauxy_d_left(2,j)   = 0.5*( vMu(1,j+1) + vMu(1,j) ) & !Mu_avg @ i = 1
															!Avg Back v_Grad
         *(0.5 * (v_old(2,j+1)-v_old(1,j+1) &             !Back_Grad @ j=j+1
     			+v_old(2,j)-v_old(1,j) )/ddx(1)   &          !Back_Grad @ j=j
                                               !u_Grad along Seg @ fixed i= 1
     	+ (u_old(1,j+1)-u_old(1,j))/ddy(j))

	Tauxx_d_left(2,j)= (1.0/3.0)*(vMu(1,j)+vMu(1,j+1)) & !Mu_avg @  i=1
                                                              !Avg Back u_Grad
         * ((u_old(2,j)-u_old(1,j) &                      !Back_Grad @ j
               +u_old(2,j+1)-u_old(1,j+1)  )/ddx(1)   &       !Back_Grad @ j+1
	                                         !v_Grad along Seg @ fixed i= 1
     	-  (v_old(1,j+1)-v_old(1,j))/ddy(j))

	qxx_d_left(2,j)= -0.25*(tK(1,j)+tK(1,j+1))  &      !tK_avg @ i = i-1
                      * ((T_old(2,j)-T_old(1,j) &       !T_back_Grad @ j=j
                      +  T_old(2,j+1)-T_old(1,j+1))/ddx(1))!T_back_Grad @ j=j+1
	enddo
End Subroutine L_vis_bc
