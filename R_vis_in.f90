! Calculate the stresses and heat transfer on the right side of (a,b,c,d) cells
Subroutine R_vis_in
use Variables
use omp_lib
implicit none

 !$OMP PARALLEL PRIVATE (i,j)
 !$OMP DO
	do j =2,Jmax-1
	do i =2,Imax-2

	! cell (a)
	Tauxy_a_right(i,j)= 0.5*(vMu(i,j)+vMu(i,j-1)  )&           !Mu_avg @ i = i
															  !Avg Forw v_Grad
         *(0.5 * (v_old(i+1,j)-v_old(i,j)&                      !Forw_Grad @ j
                +v_old(i+1,j-1)-v_old(i,j-1)  )/ddx(i) &            !Forw_Grad @ j-1
                                               !u_Grad along Seg @ fixed i= i
     	+  (u_old(i,j)-u_old(i,j-1))/ddy(j-1))

	Tauxx_a_right(i,j)= (1.0/3.0)*(vMu(i,j)+vMu(i,j-1))&     !Mu_avg @ i = i
                                                              !Avg Forw u_Grad
         * ((u_old(i+1,j)-u_old(i,j)&                         !Forw_Grad @ j
               +u_old(i+1,j-1)-u_old(i,j-1)  )/ddx(i) &            !Forw_Grad @ j-1
                                              !v_Grad along Seg @ fixed i= i
     	-  (v_old(i,j)-v_old(i,j-1))/ddy(j-1))

	qxx_a_right(i,j)= -0.25*(tK(i,j)+tK(i,j-1)) &            !tK_avg @ i = i
                      * ((T_old(i+1,j)-T_old(i,j) &           !T_forw_Grad @ j
                      +  T_old(i+1,j-1)-T_old(i,j-1))/ddx(i))    !T_forw_Grad @ j-1

	! cell (b)
	Tauxy_b_right(i,j)= 0.5*(vMu(i+1,j)+vMu(i+1,j-1)  ) &    !Mu_avg @ i= i+1
															!Avg Forw v_Grad
         *(0.5 * (v_old(i+2,j)-v_old(i+1,j)  &                !Forw_Grad @ j
                +v_old(i+2,j-1)-v_old(i+1,j-1)  )/ddx(i+1) &        !Forw_Grad @ j-1
                                               !u_Grad along Seg @ fixed i= i+1
     	+  (u_old(i+1,j)-u_old(i+1,j-1))/ddy(j-1))

	Tauxx_b_right(i,j)= (1.0/3.0)*(vMu(i+1,j)+vMu(i+1,j-1)) & !Mu_avg @ i = i+1
                                                              !Avg Forw u_Grad
         * ((u_old(i+2,j)-u_old(i+1,j) &                      !Forw_Grad @ j
               +u_old(i+2,j-1)-u_old(i+1,j-1)  )/ddx(i+1)  &        !Forw_Grad @ j-1
                                              !v_Grad along Seg @ fixed i= i+1
     	-  (v_old(i+1,j)-v_old(i+1,j-1))/ddy(j-1))


	qxx_b_right(i,j)= -0.25*(tK(i+1,j)+tK(i+1,j-1)) &        !tK_avg @ i = i+1
                      * ((T_old(i+2,j)-T_old(i+1,j) &         !T_forw_Grad @ j
                      +  T_old(i+2,j-1)-T_old(i+1,j-1))/ddx(i+1))  !T_forw_Grad @ j-1

	enddo
	enddo
	!$OMP END DO

	!$OMP DO
 	do j =2,Jmax-1
 	do i =2,Imax-2
	! cell (c)

	Tauxy_c_right(i,j)= 0.5*(vMu(i+1,j)+vMu(i+1,j+1)  ) &    !Mu_avg @ i= i+1
															!Avg Forw v_Grad
         *(0.5 * (v_old(i+2,j)-v_old(i+1,j) &                !Forw_Grad @ j
                +v_old(i+2,j+1)-v_old(i+1,j+1)  )/ddx(i+1)  &       !Forw_Grad @ j+1
                                               !u_Grad along Seg @ fixed i= i+1
     	+  (u_old(i+1,j+1)-u_old(i+1,j))/ddy(j))

	Tauxx_c_right(i,j)= (1.0/3.0)*(vMu(i+1,j)+vMu(i+1,j+1)) & !Mu_avg @ i = i+1
                                                              !Avg Forw u_Grad
         * ((u_old(i+2,j)-u_old(i+1,j) &                      !Forw_Grad @ j
               +u_old(i+2,j+1)-u_old(i+1,j+1)  )/ddx(i+1) &         !Forw_Grad @ j+1
                                              !v_Grad along Seg @ fixed i= i+1
     	-  (v_old(i+1,j+1)-v_old(i+1,j))/ddy(j))

	qxx_c_right(i,j)= -0.25*(tK(i+1,j)+tK(i+1,j+1)) &       !tK_avg @ i = i+1
                      * ((T_old(i+2,j)-T_old(i+1,j)  &      !T_forw_Grad @ j
                      +  T_old(i+2,j+1)-T_old(i+1,j+1))/ddx(i+1)) !T_forw_Grad @ j+1


	! cell (d)
	Tauxy_d_right(i,j)= 0.5*(vMu(i,j)+vMu(i,j+1)  ) &        !Mu_avg @ i= i
															!Avg Forw v_Grad
        *(0.5 * (v_old(i+1,j)-v_old(i,j)  &                  !Forw_Grad @ j
                +v_old(i+1,j+1)-v_old(i,j+1)  )/ddx(i)  &         !Forw_Grad @ j+1
                                               !u_Grad along Seg @ fixed i= i
     	+  (u_old(i,j+1)-u_old(i,j))/ddy(j))

	Tauxx_d_right(i,j)= (1.0/3.0)*(vMu(i,j)+vMu(i,j+1)) &    !Mu_avg @ i = i
                                                              !Avg Forw u_Grad
         * ((u_old(i+1,j)-u_old(i,j)  &                       !Forw_Grad @ j=j
               +u_old(i+1,j+1)-u_old(i,j+1)  )/ddx(i)  &          !Forw_Grad @ j=j+1
                                              !v_Grad along Seg @ fixed i= i
     	-  (v_old(i,j+1)-v_old(i,j))/ddy(j))

	qxx_d_right(i,j)= -0.25*(tK(i,j)+tK(i,j+1)) &          !tK_avg @ i = i
                      * ((T_old(i+1,j)-T_old(i,j) &         !T_forw_Grad @ j
                      +  T_old(i+1,j+1)-T_old(i,j+1))/ddx(i))  !T_forw_Grad @ j+1
	enddo
	enddo
	!$OMP END DO
	!$OMP END PARALLEL 
End Subroutine R_vis_in
