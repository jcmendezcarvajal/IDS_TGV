!	  stresses and heat transfer on the left side of internal (a,b,c,d) cells
Subroutine L_vis_in
use Variables
use omp_lib
implicit none

  !$OMP PARALLEL PRIVATE (i,j)
  !$OMP DO
	do j =2,Jmax-1
	do i =3,Imax-1
! Calculate the stresses and heat transfer on the left side of (a,b,c,d) cells

! cell (a)
    Tauxy_a_left(i,j)   = 0.5*( vMu(i-1,j) + vMu(i-1,j-1) )     &   !Mu_avg @ i = i-1, !Avg Back v_Grad
         *(0.5 * (v_old(i-1,  j) - v_old(i-2,  j)               &
     			+ v_old(i-1,j-1) - v_old(i-2,j-1) )/ddx(i-2)    &   !Back_Grad @ j, !Back_Grad @ j-1
     	+ (u_old(i-1,j)-u_old(i-1,j-1))/ddy(j-1))               !u_Grad along Seg @ fixed i= i-1


	Tauxx_a_left(i,j)= (1.0/3.0)*(vMu(i-1,j)+vMu(i-1,j-1))&  !Mu_avg @  i-1
                                                              !Avg Back u_Grad
         * ((u_old(i-1,j)-u_old(i-2,j)&                       !Back_Grad @ j
               +u_old(i-1,j-1)-u_old(i-2,j-1)  )/ddx(i-2) &          !Back_Grad @ j-1
	                                         !v_Grad along Seg @ fixed i= i-1
     	-  (v_old(i-1,j)-v_old(i-1,j-1))/ddy(j-1))


	qxx_a_left(i,j)= -0.25*(tK(i-1,j)+tK(i-1,j-1))&          !tK_avg @ i = i-1
                      * ((T_old(i-1,j)-T_old(i-2,j)&          !T_back_Grad @ j
                      +  T_old(i-1,j-1)-T_old(i-2,j-1))/ddx(i-2))  !T_back_Grad @ j-1

	! cell (b)
	Tauxy_b_left(i,j)   = 0.5*( vMu(i,j) + vMu(i,j-1) )&     !Mu_avg @ i = i
															!Avg Back v_Grad
         *(0.5 * (v_old(i,j)-v_old(i-1,j)&                    !Back_Grad @ j
     			+v_old(i,j-1)-v_old(i-1,j-1) )/ddx(i-1) &          !Back_Grad @ j-1
                                               !u_Grad along Seg @ fixed i= i
     	+ (u_old(i,j)-u_old(i,j-1))/ddy(j-1))

	Tauxx_b_left(i,j)= (1.0/3.0)*(vMu(i,j)+vMu(i,j-1))&      !Mu_avg @  i=i
                                                              !Avg Back u_Grad
         * ((u_old(i,j)-u_old(i-1,j)&                         !Back_Grad @ j
               +u_old(i,j-1)-u_old(i-1,j-1)  )/ddx(i-1) &            !Back_Grad @ j-1
	                                         !v_Grad along Seg @ fixed i= i
     	-  (v_old(i,j)-v_old(i,j-1))/ddy(j-1))

	qxx_b_left(i,j)= -0.25*(tK(i,j)+tK(i,j-1))&              !tK_avg @ i = i
                      * ((T_old(i,j)-T_old(i-1,j)&            !T_back_Grad @ j
                      +  T_old(i,j-1)-T_old(i-1,j-1))/ddx(i-1))    !T_back_Grad @ j-1

  enddo
  enddo
  !$OMP END DO

  !$OMP DO
	do j =2,Jmax-1
	do i =3,Imax-1
  ! cell (c)

	Tauxy_c_left(i,j)   = 0.5*( vMu(i,j+1) + vMu(i,j) )&     !Mu_avg @ i = i
															!Avg Back v_Grad
         *(0.5 * (v_old(i,j)-v_old(i-1,j)&                    !Back_Grad @ j
     			+v_old(i,j+1)-v_old(i-1,j+1) )/ddx(i-1) &          !Back_Grad @ j+1
                                               !u_Grad along Seg @ fixed i= i
     	+ (u_old(i,j+1)-u_old(i,j))/ddy(j))

	Tauxx_c_left(i,j)= (1.0/3.0)*(vMu(i,j)+vMu(i,j+1))&      !Mu_avg @  i=i
                                                              !Avg Back u_Grad
         * ((u_old(i,j)-u_old(i-1,j)&                         !Back_Grad @ j
               +u_old(i,j+1)-u_old(i-1,j+1)  )/ddx(i-1) &           !Back_Grad @ j+1
	                                         !v_Grad along Seg @ fixed i= i
     	-  (v_old(i,j+1)-v_old(i,j))/ddy(j))
	qxx_c_left(i,j)= -0.25*(tK(i,j)+tK(i,j+1))&            !tK_avg @ i = i
                      * ((T_old(i,j)-T_old(i-1,j)&          !T_back_Grad @ j=j
                      +  T_old(i,j+1)-T_old(i-1,j+1))/ddx(i-1))  !T_back_Grad @ j=j+1

	! cell (d)

	Tauxy_d_left(i,j)   = 0.5*( vMu(i-1,j+1) + vMu(i-1,j) )& !Mu_avg @ i = i-1
															!Avg Back v_Grad
         *(0.5 * (v_old(i-1,j+1)-v_old(i-2,j+1)&              !Back_Grad @ j=j+1
     			+v_old(i-1,j)-v_old(i-2,j) )/ddx(i-2)&             !Back_Grad @ j=j
                                               !u_Grad along Seg @ fixed i= i-1
     	+ (u_old(i-1,j+1)-u_old(i-1,j))/ddy(j))

	Tauxx_d_left(i,j)= (1.0/3.0)*(vMu(i-1,j)+vMu(i-1,j+1))&  !Mu_avg @  i=i-1
                                                              !Avg Back u_Grad
         * ((u_old(i-1,j)-u_old(i-2,j)&                       !Back_Grad @ j
               +u_old(i-1,j+1)-u_old(i-2,j+1)  )/ddx(i-2) &         !Back_Grad @ j+1
	                                         !v_Grad along Seg @ fixed i= i-1
     	-  (v_old(i-1,j+1)-v_old(i-1,j))/ddy(j))

	qxx_d_left(i,j)= -0.25*(tK(i-1,j)+tK(i-1,j+1))&        !tK_avg @ i = i-1
                      * ((T_old(i-1,j)-T_old(i-2,j)&        !T_back_Grad @ j=j
                      +  T_old(i-1,j+1)-T_old(i-2,j+1))/ddx(i-2))!T_back_Grad @ j=j+1
	enddo
	enddo
	!$OMP END DO
  !$OMP END PARALLEL
End Subroutine L_vis_in
