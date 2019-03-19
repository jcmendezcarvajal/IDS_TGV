! Calculate the stresses and heat transfer on the lower side of (a,b,c,d) boundary cells

Subroutine Lower_vis_bc
use Variables
use omp_lib
implicit none
	do i =2,Imax-1

	! cell (a)
	Tauxy_a_lower(i,2)	= 0.5*( vMu(i-1,1) + vMu(i,1) ) &	!Mu_avg @ j = 1
															!Avg Back U_Grad
     	*( 0.5 * (u_old(i-1,2)-u_old(i-1,1)	&			!Back_Grad @ i-1
                  +u_old(i,2)-u_old(i,1)   )/ddy(1) &			!Back_Grad @ i
											!V_Grad along Seg @ fixed j= 1
     	+  (v_old(i,1)-v_old(i-1,1))/ddx(i-1))

	Tauyy_a_lower(i,2)= (1.0/3.0)*(vMu(i-1,1)+vMu(i,1)) & !Mu_avg @  j=1
                                                              !Avg Back v_Grad
         * ((v_old(i-1,2)-v_old(i-1,1)   &                !Back_Grad @ i-1
               +v_old(i,2)-v_old(i,1))/ddy(1)   &             !Back_Grad @ i
                                              !u_Grad along Seg @ fixed j= 1
     	-  (u_old(i,1)-u_old(i-1,1))/ddx(i-1))

	qyy_a_lower(i,2)= -0.25*(tK(i-1,1)+tK(i,1))  &       !tK_avg @ j = 1
                      * ((T_old(i-1,2)-T_old(i-1,1) &     !T_back_Grad @ i-1
                      +  T_old(i,2)-T_old(i,1))/ddy(1))      !T_back_Grad @ i

	!cell (b)

	Tauxy_b_lower(i,2)	= 0.5*( vMu(i+1,1) + vMu(i,1) )	& !Mu_avg @ j = 1
															!Avg Back U_Grad
     	*( 0.5 * (u_old(i+1,2)-u_old(i+1,1)	&			!Back_Grad @ i+1
                  +u_old(i,2)-u_old(i,1)   )/ddy(1) &			!Back_Grad @ i
											!V_Grad along Seg @ fixed j= 1
     	+  (v_old(i+1,1)-v_old(i,1))/ddx(i)     )

	Tauyy_b_lower(i,2)= (1.0/3.0)*(vMu(i,1)+vMu(i+1,1)) & !Mu_avg @  j=1
                                                              !Avg Back v_Grad
         * ((v_old(i,2)-v_old(i,1)  &                     !Back_Grad @ i
               +v_old(i+1,2)-v_old(i+1,1))/ddy(1)  &          !Back_Grad @ i+1
                                              !u_Grad along Seg @ fixed j= 1
     	-  (u_old(i+1,1)-u_old(i,1))/ddx(i))

	qyy_b_lower(i,2)= -0.25*(tK(i,1)+tK(i+1,1))   &      !tK_avg @ j = 1
                      * ((T_old(i,2)-T_old(i,1)  &        !T_back_Grad @ i=i
                      +  T_old(i+1,2)-T_old(i+1,1))/ddy(1))  !T_back_Grad @ i+1

	!cell (c)

	Tauxy_c_lower(i,2)	= 0.5*( vMu(i,2) + vMu(i+1,2) )	& !Mu_avg @ j = 2
															!Avg Back U_Grad
     	*( 0.5 * (u_old(i+1,2)-u_old(i+1,1)	&			!Back_Grad @ i+1
                  +u_old(i,2)-u_old(i,1)   )/ddy(1) &			!Back_Grad @ i
											!V_Grad along Seg @ fixed j= 2
     	+  (v_old(i+1,2)-v_old(i,2))/ddx(i)     )

	Tauyy_c_lower(i,2)= (1.0/3.0)*(vMu(i,2)+vMu(i+1,2)) & !Mu_avg @  j=2
                                                          !Avg Back v_Grad
         * ((v_old(i,2)-v_old(i,1)     &                !Back_Grad @ i=i
               +v_old(i+1,2)-v_old(i+1,1))/ddy(1) &         !Back_Grad @ i+1
                                              !u_Grad along Seg @ fixed j= 2
     	-  (u_old(i+1,2)-u_old(i,2))/ddx(i))

	qyy_c_lower(i,2)= -0.25*(tK(i,2)+tK(i+1,2)) &          !tK_avg @ j = 2
                      * ((T_old(i,2)-T_old(i,1)  &        !T_back_Grad @ i=i
                      +  T_old(i+1,2)-T_old(i+1,1))/ddy(1))  !T_back_Grad @ i=i+1

	!cell(d)

	Tauxy_d_lower(i,2)	= 0.5*( vMu(i,2) + vMu(i-1,2) )	 &   !Mu_avg @ j = 2
															!Avg Back U_Grad
     	*( 0.5 * (u_old(i-1,2)-u_old(i-1,1)	&			!Back_Grad @ i-1
                  +u_old(i,2)-u_old(i,1)   )/ddy(1) &			!Back_Grad @ i=i
											!V_Grad along Seg @ fixed j= 2
     	+  (v_old(i,2)-v_old(i-1,2))/ddx(i-1)     )

      Tauyy_d_lower(i,2)= (1.0/3.0)*(vMu(i,2)+vMu(i-1,2))  &   !Mu_avg @  j=2
                                                              !Avg Back v_Grad
         * ((v_old(i,2)-v_old(i,1)   &                      !Back_Grad @ i=i
               +v_old(i-1,2)-v_old(i-1,1))/ddy(1)  &            !Back_Grad @ i=i-1
                                              !u_Grad along Seg @ fixed j= 2
     	-  (u_old(i,2)-u_old(i-1,2))/ddx(i-1))

	qyy_d_lower(i,2)= -0.25*(tK(i,2)+tK(i-1,2)) &          !tK_avg @ j = 2
                      * ((T_old(i,2)-T_old(i,1) &         !T_back_Grad @ i=i
                      +  T_old(i-1,2)-T_old(i-1,1))/ddy(1))  !T_back_Grad @ i=i-1
      enddo
End Subroutine Lower_vis_bc
