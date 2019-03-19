! Calculate the stresses and heat transfer on the upper side of (a,b,c,d) boundary cells
Subroutine Upper_vis_bc
use Variables
use omp_lib
implicit none
	do i =2,Imax-1
	!cell (a)                                                    !Mu_avg @ j =Jmax-1
	Tauxy_a_upper(i,Jmax-1)	= 0.5*(vMu(i-1,Jmax-1)+vMu(i,Jmax-1)) &
															!Avg Forw U_Grad
         *( 0.5 * (u_old(i-1,Jmax)-u_old(i-1,Jmax-1) &			!Forw_Grad @ i-1
                  +u_old(i,Jmax)-u_old(i,Jmax-1)   )/ddy(Jmax-1) &		!Forw_Grad @ i
											!V_Grad along Seg @ fixed j= Jmax-1
     	+  (v_old(i,Jmax-1)-v_old(i-1,Jmax-1))/ddx(i-1))
                                                                      !Mu_avg @ j =Jmax-1
	Tauyy_a_upper(i,Jmax-1)= (1.0/3.0)*(vMu(i-1,Jmax-1)+vMu(i,Jmax-1)) &                                                                  !Avg Forw v_Grad
         * ((v_old(i-1,Jmax)-v_old(i-1,Jmax-1)    &             !Forw_Grad @ i-1
               +v_old(i,Jmax)-v_old(i,Jmax-1))/ddy(Jmax-1)  &            !Forw_Grad @ i-1
                                               !u_Grad along Seg @ fixed j= Jmax-1
     	-  (u_old(i,Jmax-1)-u_old(i-1,Jmax-1))/ddx(i-1)  )

	                                                          !tK_avg @ j = Jmax-1
	qyy_a_upper(i,Jmax-1)= -0.25*(tK(i-1,Jmax-1)+tK(i,Jmax-1)) &
                      * ((T_old(i,Jmax)-T_old(i,Jmax-1)  &     !T_forw_Grad @ i
                    +T_old(i-1,Jmax)-T_old(i-1,Jmax-1))/ddy(Jmax-1))  !T_forw_Grad @ i-1

	! cell (b)
	                                                               !Mu_avg @ j = Jmax-1
	Tauxy_b_upper(i,Jmax-1)	= 0.5*(vMu(i+1,Jmax-1)+vMu(i,Jmax-1)) &
															!Avg Forw U_Grad
         *( 0.5 * (u_old(i+1,Jmax)-u_old(i+1,Jmax-1) &			!Forw_Grad @ i+1
                  +u_old(i,Jmax)-u_old(i,Jmax-1)   )/ddy(Jmax-1) &	    !Forw_Grad @ i
											!V_Grad along Seg @ fixed j= Jmax-1
     	+  (v_old(i+1,Jmax-1)-v_old(i,Jmax-1))/ddx(i))

                                                                        !Mu_avg @ j = Jmax-1
	Tauyy_b_upper(i,Jmax-1)= (1.0/3.0)*(vMu(i+1,Jmax-1)+vMu(i,Jmax-1)) &
                                                              !Avg Forw v_Grad
         * ((v_old(i,Jmax)-v_old(i,Jmax-1)    &                !Forw_Grad @ i=i
               +v_old(i+1,Jmax)-v_old(i+1,Jmax-1))/ddy(Jmax-1)  &       !Forw_Grad @ i+1
                                             !u_Grad along Seg @ fixed j= Jmax-1
     	-  (u_old(i+1,Jmax-1)-u_old(i,Jmax-1))/ddx(i)  )
                                                                 !tK_avg @ j = Jmax-1
	qyy_b_upper(i,Jmax-1)= -0.25*(tK(i,Jmax-1)+tK(i+1,Jmax-1))  &
                      * ((T_old(i,Jmax)-T_old(i,Jmax-1)   &   !T_forw_Grad @ i
                     + T_old(i+1,Jmax)-T_old(i+1,Jmax-1))/ddy(Jmax-1))!T_forw_Grad @ i+1

	!cell (c)	                                               !Mu_avg @ j = Jmax
	Tauxy_c_upper(i,Jmax-1)	= 0.5*(vMu(i,Jmax)+vMu(i+1,Jmax)) &
															!Avg Forw U_Grad
         *( 0.5 * (u_old(i+1,Jmax)-u_old(i+1,Jmax-1)	&		!Forw_Grad @ i+1
                  +u_old(i,Jmax)-u_old(i,Jmax-1)   )/ddy(Jmax-1) &		!Forw_Grad @ i
											!V_Grad along Seg @ fixed j= Jmax
     	+  (v_old(i+1,Jmax)-v_old(i,Jmax))/ddx(i))

	                                                              !Mu_avg @ j = Jmax
	Tauyy_c_upper(i,Jmax-1)= (1.0/3.0)*(vMu(i,Jmax)+vMu(i+1,Jmax)) &
                                                             !Avg Forw v_Grad
         * ((v_old(i,Jmax)-v_old(i,Jmax-1)  &                  !Forw_Grad @ i=i
               +v_old(i+1,Jmax)-v_old(i+1,Jmax-1))/ddy(Jmax-1)   &      !Forw_Grad @ i+1
                                             !u_Grad along Seg @ fixed j= Jmax
     	-  (u_old(i+1,Jmax)-u_old(i,Jmax))/ddx(i)  )
                                                            !tK_avg @ j = Jmax
	qyy_c_upper(i,Jmax-1)= -0.25*(tK(i,Jmax)+tK(i+1,Jmax)) &
                      * ((T_old(i,Jmax)-T_old(i,Jmax-1) &   !T_forw_Grad @ i=i
                   + T_old(i+1,Jmax)-T_old(i+1,Jmax-1))/ddy(Jmax-1))!T_forw_Grad @ i+1

	! cell (d)
	                                                         !Mu_avg @ j = Jmax
	Tauxy_d_upper(i,Jmax-1)	= 0.5*(vMu(i,Jmax)+vMu(i-1,Jmax)) &
															!Avg Forw U_Grad
         *( 0.5 * (u_old(i-1,Jmax)-u_old(i-1,Jmax-1) &			!Forw_Grad @ i-1
                  +u_old(i,Jmax)-u_old(i,Jmax-1)   )/ddy(Jmax-1) &		!Forw_Grad @ i
											!V_Grad along Seg @ fixed j= Jmax
     	+  (v_old(i,Jmax)-v_old(i-1,Jmax))/ddx(i-1))
                                                                    !Mu_avg @ j = Jmax
	Tauyy_d_upper(i,Jmax-1)= (1.0/3.0)*(vMu(i,Jmax)+vMu(i-1,Jmax)) &
                                                              !Avg Forw v_Grad
         * ((v_old(i,Jmax)-v_old(i,Jmax-1)  &                 !Forw_Grad @ i=i
               +v_old(i-1,Jmax)-v_old(i-1,Jmax-1))/ddy(Jmax-1) &       !Forw_Grad @ i-1
                                             !u_Grad along Seg @ fixed j= Jmax
     	-  (u_old(i,Jmax)-u_old(i-1,Jmax))/ddx(i-1)  )

	qyy_d_upper(i,Jmax-1)= -0.25*(tK(i,Jmax)+tK(i-1,Jmax)) & !tK_avg @ j = Jmax
                      * ((T_old(i,Jmax)-T_old(i,Jmax-1) &     !T_forw_Grad @ i=i
                    + T_old(i-1,Jmax)-T_old(i-1,Jmax-1))/ddy(Jmax-1)) !T_forw_Grad @ i-1
	enddo
End Subroutine Upper_vis_bc
