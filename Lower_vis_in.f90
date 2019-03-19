! Calculate the stresses and heat transfer on the lower side of (a,b,c,d) cells
Subroutine Lower_vis_in
use Variables
use omp_lib
implicit none

	!$OMP PARALLEL PRIVATE (i,j)
	!$OMP DO
	do j =3,Jmax-1
	do i =2,Imax-1

	! cell (a)
	Tauxy_a_lower(i,j)	= 0.5*( vMu(i-1,j-1) + vMu(i,j-1) ) &	!Mu_avg @ j = j-1
															!Avg Back U_Grad
     	*( 0.5 * (u_old(i-1,j-1)-u_old(i-1,j-2)&				!Back_Grad @ i-1
                  +u_old(i,j-1)-u_old(i,j-2)   )/ddy(j-2) &			!Back_Grad @ i
											!V_Grad along Seg @ fixed j= j-1
     	+  (v_old(i,j-1)-v_old(i-1,j-1))/ddx(i-1))

	Tauyy_a_lower(i,j)= (1.0/3.0)*(vMu(i-1,j-1)+vMu(i,j-1)) & !Mu_avg @  j-1
                                                              !Avg Back v_Grad
         * ((v_old(i-1,j-1)-v_old(i-1,j-2) &                   !Back_Grad @ i-1
               +v_old(i,j-1)-v_old(i,j-2))/ddy(j-2) &               !Back_Grad @ i
                                              !u_Grad along Seg @ fixed j= j-1
     	-  (u_old(i,j-1)-u_old(i-1,j-1))/ddx(i-1))

	qyy_a_lower(i,j)= -0.25*(tK(i-1,j-1)+tK(i,j-1)) &        !tK_avg @ j = j-1
                      * ((T_old(i-1,j-1)-T_old(i-1,j-2) &     !T_back_Grad @ i-1
                      +  T_old(i,j-1)-T_old(i,j-2))/ddy(j-2))      !T_back_Grad @ i

	!cell (b)

	Tauxy_b_lower(i,j)	= 0.5*( vMu(i+1,j-1) + vMu(i,j-1) )&	!Mu_avg @ j = j-1
															!Avg Back U_Grad
     	*( 0.5 * (u_old(i+1,j-1)-u_old(i+1,j-2)&				!Back_Grad @ i+1
                  +u_old(i,j-1)-u_old(i,j-2)   )/ddy(j-2) &			!Back_Grad @ i
											!V_Grad along Seg @ fixed j= j-1
     	+  (v_old(i+1,j-1)-v_old(i,j-1))/ddx(i)     )

	Tauyy_b_lower(i,j)= (1.0/3.0)*(vMu(i,j-1)+vMu(i+1,j-1))& !Mu_avg @  j-1
                                                              !Avg Back v_Grad
         * ((v_old(i,j-1)-v_old(i,j-2)&                       !Back_Grad @ i
               +v_old(i+1,j-1)-v_old(i+1,j-2))/ddy(j-2) &            !Back_Grad @ i+1
                                              !u_Grad along Seg @ fixed j= j-1
     	-  (u_old(i+1,j-1)-u_old(i,j-1))/ddx(i))

	qyy_b_lower(i,j)= -0.25*(tK(i,j-1)+tK(i+1,j-1)) &        !tK_avg @ j = j-1
                      * ((T_old(i,j-1)-T_old(i,j-2) &         !T_back_Grad @ i=i
                      +  T_old(i+1,j-1)-T_old(i+1,j-2))/ddy(j-2))  !T_back_Grad @ i+1

 enddo
 enddo
 !$OMP END DO

 !$OMP DO
 do j =3,Jmax-1
 do i =2,Imax-1

 !cell (c)

	Tauxy_c_lower(i,j)	= 0.5*( vMu(i,j) + vMu(i+1,j) )     &	!Mu_avg @ j = j
															!Avg Back U_Grad
     	*( 0.5 * (u_old(i+1,j)-u_old(i+1,j-1)&				!Back_Grad @ i+1
                  +u_old(i,j)-u_old(i,j-1)   )/ddy(j-1)     &			!Back_Grad @ i
											!V_Grad along Seg @ fixed j= j
     	+  (v_old(i+1,j)-v_old(i,j))/ddx(i)     )

	Tauyy_c_lower(i,j)= (1.0/3.0)*(vMu(i,j)+vMu(i+1,j))& !Mu_avg @  j=j
                                                          !Avg Back v_Grad
         * ((v_old(i,j)-v_old(i,j-1)&                     !Back_Grad @ i=i
               +v_old(i+1,j)-v_old(i+1,j-1))/ddy(j-1) &         !Back_Grad @ i+1
                                              !u_Grad along Seg @ fixed j= j
     	-  (u_old(i+1,j)-u_old(i,j))/ddx(i))

	qyy_c_lower(i,j)= -0.25*(tK(i,j)+tK(i+1,j)) &          !tK_avg @ j = j
                      * ((T_old(i,j)-T_old(i,j-1) &         !T_back_Grad @ i=i
                      +  T_old(i+1,j)-T_old(i+1,j-1))/ddy(j-1))  !T_back_Grad @ i=i+1

	!cell(d)

	Tauxy_d_lower(i,j)	= 0.5*( vMu(i,j) + vMu(i-1,j) )	&    !Mu_avg @ j = j
															!Avg Back U_Grad
     	*( 0.5 * (u_old(i-1,j)-u_old(i-1,j-1)&				!Back_Grad @ i-1
                  +u_old(i,j)-u_old(i,j-1)   )/ddy(j-1) &			!Back_Grad @ i=i
											!V_Grad along Seg @ fixed j= j
     	+  (v_old(i,j)-v_old(i-1,j))/ddx(i-1)     )

      Tauyy_d_lower(i,j)= (1.0/3.0)*(vMu(i,j)+vMu(i-1,j)) &    !Mu_avg @  j=j
                                                              !Avg Back v_Grad
         * ((v_old(i,j)-v_old(i,j-1) &                        !Back_Grad @ i=i
               +v_old(i-1,j)-v_old(i-1,j-1))/ddy(j-1) &             !Back_Grad @ i=i-1
                                              !u_Grad along Seg @ fixed j= j
     	-  (u_old(i,j)-u_old(i-1,j))/ddx(i-1))

	qyy_d_lower(i,j)= -0.25*(tK(i,j)+tK(i-1,j)) &          !tK_avg @ j = j
                      * ((T_old(i,j)-T_old(i,j-1) &         !T_back_Grad @ i=i
                      +  T_old(i-1,j)-T_old(i-1,j-1))/ddy(j-1))  !T_back_Grad @ i=i-1
  enddo
	enddo
	!$OMP END DO
	!$OMP END PARALLEL
End Subroutine Lower_vis_in
