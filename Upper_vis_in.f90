! Calculate the stresses and heat transfer on the upper side of (a,b,c,d) cells
Subroutine Upper_vis_in
use Variables
use omp_lib
implicit none

  !$OMP PARALLEL PRIVATE (i,j)
  !$OMP DO
	do j =2,Jmax-2
	do i =2,Imax-1
	!cell (a)
	Tauxy_a_upper(i,j)	= 0.5*(vMu(i-1,j)+vMu(i,j))&			!Mu_avg @ j = j
															!Avg Forw U_Grad
         *( 0.5 * (u_old(i-1,j+1)-u_old(i-1,j)&				!Forw_Grad @ i-1
                  +u_old(i,j+1)-u_old(i,j)   )/ddy(j)   &			!Forw_Grad @ i
											!V_Grad along Seg @ fixed j= j
     	+  (v_old(i,j)-v_old(i-1,j))/ddx(i-1))

	Tauyy_a_upper(i,j)= (1.0/3.0)*(vMu(i-1,j)+vMu(i,j)) &    !Mu_avg @ j = j
                                                              !Avg Forw v_Grad
         * ((v_old(i-1,j+1)-v_old(i-1,j)  &                   !Forw_Grad @ i-1
               +v_old(i,j+1)-v_old(i,j))/ddy(j)  &                !Forw_Grad @ i-1
                                               !u_Grad along Seg @ fixed j= j
     	-  (u_old(i,j)-u_old(i-1,j))/ddx(i-1)  )


	qyy_a_upper(i,j)= -0.25*(tK(i-1,j)+tK(i,j))&             !tK_avg @ j = j
                      * ((T_old(i,j+1)-T_old(i,j)&            !T_forw_Grad @ i
                      +  T_old(i-1,j+1)-T_old(i-1,j))/ddy(j))    !T_forw_Grad @ i-1

	! cell (b)

	Tauxy_b_upper(i,j)	= 0.5*(vMu(i+1,j)+vMu(i,j))&			!Mu_avg @ j = j
															!Avg Forw U_Grad
         *( 0.5 * (u_old(i+1,j+1)-u_old(i+1,j)&				!Forw_Grad @ i+1
                  +u_old(i,j+1)-u_old(i,j)   )/ddy(j)   &			!Forw_Grad @ i
											!V_Grad along Seg @ fixed j= j
     	+  (v_old(i+1,j)-v_old(i,j))/ddx(i))


	Tauyy_b_upper(i,j)= (1.0/3.0)*(vMu(i,j)+vMu(i+1,j))&     !Mu_avg @ j = j
                                                              !Avg Forw v_Grad
         * ((v_old(i,j+1)-v_old(i,j) &                        !Forw_Grad @ i=i
               +v_old(i+1,j+1)-v_old(i+1,j))/ddy(j) &             !Forw_Grad @ i+1
                                             !u_Grad along Seg @ fixed j= j
     	-  (u_old(i+1,j)-u_old(i,j))/ddx(i)  )

	qyy_b_upper(i,j)= -0.25*(tK(i,j)+tK(i+1,j))&            !tK_avg @ j = j
                      * ((T_old(i,j+1)-T_old(i,j)&            !T_forw_Grad @ i
                      +  T_old(i+1,j+1)-T_old(i+1,j))/ddy(j))    !T_forw_Grad @ i+1

  enddo
  enddo
  !$OMP END DO

  !$OMP DO
	do j =2,Jmax-2
	do i =2,Imax-1

  !cell (c)
	Tauxy_c_upper(i,j)	= 0.5*(vMu(i,j+1)+vMu(i+1,j+1))&		!Mu_avg @ j = j+1
															!Avg Forw U_Grad
         *( 0.5 * (u_old(i+1,j+2)-u_old(i+1,j+1)&				!Forw_Grad @ i+1
                  +u_old(i,j+2)-u_old(i,j+1)   )/ddy(j+1)&			!Forw_Grad @ i
											!V_Grad along Seg @ fixed j= j+1
     	+  (v_old(i+1,j+1)-v_old(i,j+1))/ddx(i))



	Tauyy_c_upper(i,j)= (1.0/3.0)*(vMu(i,j+1)+vMu(i+1,j+1))& !Mu_avg @ j = j+1
                                                             !Avg Forw v_Grad
         * ((v_old(i,j+2)-v_old(i,j+1) &                     !Forw_Grad @ i=i
               +v_old(i+1,j+2)-v_old(i+1,j+1))/ddy(j+1) &          !Forw_Grad @ i+1
                                             !u_Grad along Seg @ fixed j= j+1
     	-  (u_old(i+1,j+1)-u_old(i,j+1))/ddx(i)  )

	qyy_c_upper(i,j)= -0.25*(tK(i,j+1)+tK(i+1,j+1)) &      !tK_avg @ j = j+1
                      * ((T_old(i,j+2)-T_old(i,j+1) &       !T_forw_Grad @ i=i
                      +  T_old(i+1,j+2)-T_old(i+1,j+1))/ddy(j+1)) !T_forw_Grad @ i+1

	! cell (d)

	Tauxy_d_upper(i,j)	= 0.5*(vMu(i,j+1)+vMu(i-1,j+1))	&	!Mu_avg @ j = j+1
															!Avg Forw U_Grad
         *( 0.5 * (u_old(i-1,j+2)-u_old(i-1,j+1) &				!Forw_Grad @ i-1
                  +u_old(i,j+2)-u_old(i,j+1)   )/ddy(j+1) &			!Forw_Grad @ i
											!V_Grad along Seg @ fixed j= j+1
     	+  (v_old(i,j+1)-v_old(i-1,j+1))/ddx(i-1))

	Tauyy_d_upper(i,j)= (1.0/3.0)*(vMu(i,j+1)+vMu(i-1,j+1)) & !Mu_avg @ j = j+1
                                                              !Avg Forw v_Grad
         * ((v_old(i,j+2)-v_old(i,j+1) &                      !Forw_Grad @ i=i
               +v_old(i-1,j+2)-v_old(i-1,j+1))/ddy(j+1) &           !Forw_Grad @ i-1
                                             !u_Grad along Seg @ fixed j= j+1
     	-  (u_old(i,j+1)-u_old(i-1,j+1))/ddx(i-1)  )

	qyy_d_upper(i,j)= -0.25*(tK(i,j+1)+tK(i-1,j+1))  &     !tK_avg @ j = j+1
                      * ((T_old(i,j+2)-T_old(i,j+1) &       !T_forw_Grad @ i=i
                      +  T_old(i-1,j+2)-T_old(i-1,j+1))/ddy(j+1))!T_forw_Grad @ i-1
	enddo
	enddo
	!$OMP END DO
  !$OMP END PARALLEL
End Subroutine Upper_vis_in
