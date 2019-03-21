!This subroutine Allocates the arrays required for the code
Subroutine InitializeArrays
use Variables
Implicit none

! Variable for geometric references (LOCAL)
ALLOCATE(x(Imax,Jmax),y(Imax,Jmax),xDim(Imax,Jmax),yDim(Imax,Jmax))
ALLOCATE(xx(Imax), ddx(Imax))
ALLOCATE(yy(Jmax), ddy(Jmax))
! Set an initial guess for the unknown primitive variables
ALLOCATE(r_old(Imax,Jmax),u_old(Imax,Jmax),v_old(Imax,Jmax),T_old(Imax,Jmax))
ALLOCATE(a_old(Imax,Jmax),Pres(Imax,Jmax))
ALLOCATE(r_temp(Imax,Jmax),u_temp(Imax,Jmax),v_temp(Imax,Jmax),T_temp(Imax,Jmax))
ALLOCATE(r_new(Imax,Jmax),u_new(Imax,Jmax),v_new(Imax,Jmax),T_new(Imax,Jmax))

!Analytical part
ALLOCATE(vorticity_exact(Imax,Jmax),psi_exact(Imax,Jmax))
ALLOCATE(u_exact(Imax,Jmax),v_exact(Imax,Jmax))


! The Time-flux Vector_U_old
ALLOCATE(U1_old(Imax,Jmax),U2_old(Imax,Jmax),U3_old(Imax,Jmax),U4_old(Imax,Jmax))
ALLOCATE(H_old(Imax,Jmax),H_temp(Imax,Jmax),Vorticity(Imax,Jmax),Psi(Imax,Jmax))
! Calculate the viscous properties
ALLOCATE(vMu(Imax,Jmax),tK(Imax,Jmax),vNu(Imax,Jmax))
! Calculate the time_step
ALLOCATE(small_time(Imax,Jmax),tNu_max(Imax,Jmax),t1(Imax,Jmax),t2(Imax,Jmax))
! Shear and Normal stresses(Tauxy),(Tauxx),(Tauyy)
ALLOCATE(Tauxy_a_lower(Imax,Jmax),Tauxy_a_upper(Imax,Jmax))
ALLOCATE(Tauxy_a_left(Imax,Jmax),Tauxy_a_right(Imax,Jmax))
ALLOCATE(Tauxx_a_left(Imax,Jmax),Tauxx_a_right(Imax,Jmax))
ALLOCATE(Tauyy_a_lower(Imax,Jmax),Tauyy_a_upper(Imax,Jmax))
ALLOCATE(Tauxy_b_lower(Imax,Jmax),Tauxy_b_upper(Imax,Jmax))
ALLOCATE(Tauxy_b_left(Imax,Jmax),Tauxy_b_right(Imax,Jmax))
ALLOCATE(Tauxx_b_left(Imax,Jmax),Tauxx_b_right(Imax,Jmax))
ALLOCATE(Tauyy_b_lower(Imax,Jmax),Tauyy_b_upper(Imax,Jmax))
ALLOCATE(Tauxy_c_lower(Imax,Jmax),Tauxy_c_upper(Imax,Jmax))
ALLOCATE(Tauxy_c_left(Imax,Jmax),Tauxy_c_right(Imax,Jmax))
ALLOCATE(Tauxx_c_left(Imax,Jmax),Tauxx_c_right(Imax,Jmax))
ALLOCATE(Tauyy_c_lower(Imax,Jmax),Tauyy_c_upper(Imax,Jmax))
ALLOCATE(Tauxy_d_lower(Imax,Jmax),Tauxy_d_upper(Imax,Jmax))
ALLOCATE(Tauxy_d_left(Imax,Jmax),Tauxy_d_right(Imax,Jmax))
ALLOCATE(Tauxx_d_left(Imax,Jmax),Tauxx_d_right(Imax,Jmax))
ALLOCATE(Tauyy_d_lower(Imax,Jmax),Tauyy_d_upper(Imax,Jmax))
! Heat transfer
ALLOCATE(qxx_a_left(Imax,Jmax),qxx_a_right(Imax,Jmax))
ALLOCATE(qyy_a_lower(Imax,Jmax),qyy_a_upper(Imax,Jmax))
ALLOCATE(qxx_b_left(Imax,Jmax),qxx_b_right(Imax,Jmax))
ALLOCATE(qyy_b_lower(Imax,Jmax),qyy_b_upper(Imax,Jmax))
ALLOCATE(qxx_c_left(Imax,Jmax),qxx_c_right(Imax,Jmax))
ALLOCATE(qyy_c_lower(Imax,Jmax),qyy_c_upper(Imax,Jmax))
ALLOCATE(qxx_d_left(Imax,Jmax),qxx_d_right(Imax,Jmax))
ALLOCATE(qyy_d_lower(Imax,Jmax),qyy_d_upper(Imax,Jmax))
! Calculate the x_direction invicid fluxes (E_inv)
ALLOCATE(E1_inv_a_left(Imax,Jmax),E2_inv_a_left(Imax,Jmax))
ALLOCATE(E3_inv_a_left(Imax,Jmax),E4_inv_a_left(Imax,Jmax))
ALLOCATE(E1_inv_a_right(Imax,Jmax),E2_inv_a_right(Imax,Jmax))
ALLOCATE(E3_inv_a_right(Imax,Jmax),E4_inv_a_right(Imax,Jmax))
ALLOCATE(E1_inv_b_left(Imax,Jmax),E2_inv_b_left(Imax,Jmax))
ALLOCATE(E3_inv_b_left(Imax,Jmax),E4_inv_b_left(Imax,Jmax))
ALLOCATE(E1_inv_b_right(Imax,Jmax),E2_inv_b_right(Imax,Jmax))
ALLOCATE(E3_inv_b_right(Imax,Jmax),E4_inv_b_right(Imax,Jmax))
ALLOCATE(E1_inv_c_left(Imax,Jmax),E2_inv_c_left(Imax,Jmax))
ALLOCATE(E3_inv_c_left(Imax,Jmax),E4_inv_c_left(Imax,Jmax))
ALLOCATE(E1_inv_c_right(Imax,Jmax),E2_inv_c_right(Imax,Jmax))
ALLOCATE(E3_inv_c_right(Imax,Jmax),E4_inv_c_right(Imax,Jmax))
ALLOCATE(E1_inv_d_left(Imax,Jmax),E2_inv_d_left(Imax,Jmax))
ALLOCATE(E3_inv_d_left(Imax,Jmax),E4_inv_d_left(Imax,Jmax))
ALLOCATE(E1_inv_d_right(Imax,Jmax),E2_inv_d_right(Imax,Jmax))
ALLOCATE(E3_inv_d_right(Imax,Jmax),E4_inv_d_right(Imax,Jmax))
! Calculate the y_direction invicid fluxes (F_inv)
ALLOCATE(F1_inv_a_lower(Imax,Jmax),F2_inv_a_lower(Imax,Jmax))
ALLOCATE(F3_inv_a_lower(Imax,Jmax),F4_inv_a_lower(Imax,Jmax))
ALLOCATE(F1_inv_a_upper(Imax,Jmax),F2_inv_a_upper(Imax,Jmax))
ALLOCATE(F3_inv_a_upper(Imax,Jmax),F4_inv_a_upper(Imax,Jmax))
ALLOCATE(F1_inv_b_lower(Imax,Jmax),F2_inv_b_lower(Imax,Jmax))
ALLOCATE(F3_inv_b_lower(Imax,Jmax),F4_inv_b_lower(Imax,Jmax))
ALLOCATE(F1_inv_b_upper(Imax,Jmax),F2_inv_b_upper(Imax,Jmax))
ALLOCATE(F3_inv_b_upper(Imax,Jmax),F4_inv_b_upper(Imax,Jmax))
ALLOCATE(F1_inv_c_lower(Imax,Jmax),F2_inv_c_lower(Imax,Jmax))
ALLOCATE(F3_inv_c_lower(Imax,Jmax),F4_inv_c_lower(Imax,Jmax))
ALLOCATE(F1_inv_c_upper(Imax,Jmax),F2_inv_c_upper(Imax,Jmax))
ALLOCATE(F3_inv_c_upper(Imax,Jmax),F4_inv_c_upper(Imax,Jmax))
ALLOCATE(F1_inv_d_lower(Imax,Jmax),F2_inv_d_lower(Imax,Jmax))
ALLOCATE(F3_inv_d_lower(Imax,Jmax),F4_inv_d_lower(Imax,Jmax))
ALLOCATE(F1_inv_d_upper(Imax,Jmax),F2_inv_d_upper(Imax,Jmax))
ALLOCATE(F3_inv_d_upper(Imax,Jmax),F4_inv_d_upper(Imax,Jmax))
! Calculate the x_direction viscous fluxes (E_vis)
ALLOCATE(E1_vis_a_left(Imax,Jmax),E2_vis_a_left(Imax,Jmax))
ALLOCATE(E3_vis_a_left(Imax,Jmax),E4_vis_a_left(Imax,Jmax))
ALLOCATE(E1_vis_a_right(Imax,Jmax),E2_vis_a_right(Imax,Jmax))
ALLOCATE(E3_vis_a_right(Imax,Jmax),E4_vis_a_right(Imax,Jmax))
ALLOCATE(E1_vis_b_left(Imax,Jmax),E2_vis_b_left(Imax,Jmax))
ALLOCATE(E3_vis_b_left(Imax,Jmax),E4_vis_b_left(Imax,Jmax))
ALLOCATE(E1_vis_b_right(Imax,Jmax),E2_vis_b_right(Imax,Jmax))
ALLOCATE(E3_vis_b_right(Imax,Jmax),E4_vis_b_right(Imax,Jmax))
ALLOCATE(E1_vis_c_left(Imax,Jmax),E2_vis_c_left(Imax,Jmax))
ALLOCATE(E3_vis_c_left(Imax,Jmax),E4_vis_c_left(Imax,Jmax))
ALLOCATE(E1_vis_c_right(Imax,Jmax),E2_vis_c_right(Imax,Jmax))
ALLOCATE(E3_vis_c_right(Imax,Jmax),E4_vis_c_right(Imax,Jmax))
ALLOCATE(E1_vis_d_left(Imax,Jmax),E2_vis_d_left(Imax,Jmax))
ALLOCATE(E3_vis_d_left(Imax,Jmax),E4_vis_d_left(Imax,Jmax))
ALLOCATE(E1_vis_d_right(Imax,Jmax),E2_vis_d_right(Imax,Jmax))
ALLOCATE(E3_vis_d_right(Imax,Jmax),E4_vis_d_right(Imax,Jmax))
! Calculate the y_direction viscous fluxes (F_vis)
ALLOCATE(F1_vis_a_lower(Imax,Jmax),F2_vis_a_lower(Imax,Jmax))
ALLOCATE(F3_vis_a_lower(Imax,Jmax),F4_vis_a_lower(Imax,Jmax))
ALLOCATE(F1_vis_a_upper(Imax,Jmax),F2_vis_a_upper(Imax,Jmax))
ALLOCATE(F3_vis_a_upper(Imax,Jmax),F4_vis_a_upper(Imax,Jmax))
ALLOCATE(F1_vis_b_lower(Imax,Jmax),F2_vis_b_lower(Imax,Jmax))
ALLOCATE(F3_vis_b_lower(Imax,Jmax),F4_vis_b_lower(Imax,Jmax))
ALLOCATE(F1_vis_b_upper(Imax,Jmax),F2_vis_b_upper(Imax,Jmax))
ALLOCATE(F3_vis_b_upper(Imax,Jmax),F4_vis_b_upper(Imax,Jmax))
ALLOCATE(F1_vis_c_lower(Imax,Jmax),F2_vis_c_lower(Imax,Jmax))
ALLOCATE(F3_vis_c_lower(Imax,Jmax),F4_vis_c_lower(Imax,Jmax))
ALLOCATE(F1_vis_c_upper(Imax,Jmax),F2_vis_c_upper(Imax,Jmax))
ALLOCATE(F3_vis_c_upper(Imax,Jmax),F4_vis_c_upper(Imax,Jmax))
ALLOCATE(F1_vis_d_lower(Imax,Jmax),F2_vis_d_lower(Imax,Jmax))
ALLOCATE(F3_vis_d_lower(Imax,Jmax),F4_vis_d_lower(Imax,Jmax))
ALLOCATE(F1_vis_d_upper(Imax,Jmax),F2_vis_d_upper(Imax,Jmax))
ALLOCATE(F3_vis_d_upper(Imax,Jmax),F4_vis_d_upper(Imax,Jmax))
! The fluxes derivatives
ALLOCATE(dU1dt_a(Imax,Jmax),dU2dt_a(Imax,Jmax),dU3dt_a(Imax,Jmax),dU4dt_a(Imax,Jmax))
ALLOCATE(dU1dt_b(Imax,Jmax),dU2dt_b(Imax,Jmax),dU3dt_b(Imax,Jmax),dU4dt_b(Imax,Jmax))
ALLOCATE(dU1dt_c(Imax,Jmax),dU2dt_c(Imax,Jmax),dU3dt_c(Imax,Jmax),dU4dt_c(Imax,Jmax))
ALLOCATE(dU1dt_d(Imax,Jmax),dU2dt_d(Imax,Jmax),dU3dt_d(Imax,Jmax),dU4dt_d(Imax,Jmax))
ALLOCATE(dU1dt(Imax,Jmax),dU2dt(Imax,Jmax),dU3dt(Imax,Jmax),dU4dt(Imax,Jmax))
! The Time-flux Vector_U_new
ALLOCATE(U1_new(Imax,Jmax),U2_new(Imax,Jmax),U3_new(Imax,Jmax),U4_new(Imax,Jmax))
! Stopping criteria
ALLOCATE(err1(Imax,Jmax),err2(Imax,Jmax),err3(Imax,Jmax),err4(Imax,Jmax))
ALLOCATE(Cf(Imax,Jmax),Cf_res(Imax,Jmax),Ch(Imax,Jmax),Ch_res(Imax,Jmax))

End Subroutine InitializeArrays

!This subroutine Deallocates the arrays required for the code
Subroutine KillArrays
use Variables
Implicit none

DEALLOCATE(r_temp,u_temp,v_temp,T_temp)
DEALLOCATE(r_new,u_new,v_new,T_new)
DEALLOCATE(U1_old,U2_old,U3_old,U4_old,H_old,H_temp)
DEALLOCATE(small_time,tNu_max,t1,t2)
DEALLOCATE(Tauxy_a_lower,Tauxy_a_upper,Tauxy_a_left,Tauxy_a_right)
DEALLOCATE(Tauxx_a_left,Tauxx_a_right,Tauyy_a_lower,Tauyy_a_upper)
DEALLOCATE(Tauxy_b_lower,Tauxy_b_upper,Tauxy_b_left,Tauxy_b_right)
DEALLOCATE(Tauxx_b_left,Tauxx_b_right,Tauyy_b_lower,Tauyy_b_upper)
DEALLOCATE(Tauxy_c_lower,Tauxy_c_upper,Tauxy_c_left,Tauxy_c_right)
DEALLOCATE(Tauxx_c_left,Tauxx_c_right,Tauyy_c_lower,Tauyy_c_upper)
DEALLOCATE(Tauxy_d_lower,Tauxy_d_upper,Tauxy_d_left,Tauxy_d_right)
DEALLOCATE(Tauxx_d_left,Tauxx_d_right,Tauyy_d_lower,Tauyy_d_upper)
DEALLOCATE(qxx_a_left,qxx_a_right,qyy_a_lower,qyy_a_upper)
DEALLOCATE(qxx_b_left,qxx_b_right,qyy_b_lower,qyy_b_upper)
DEALLOCATE(qxx_c_left,qxx_c_right,qyy_c_lower,qyy_c_upper)
DEALLOCATE(qxx_d_left,qxx_d_right,qyy_d_lower,qyy_d_upper)
DEALLOCATE(E1_inv_a_left,E2_inv_a_left,E3_inv_a_left,E4_inv_a_left)
DEALLOCATE(E1_inv_a_right,E2_inv_a_right,E3_inv_a_right,E4_inv_a_right)
DEALLOCATE(E1_inv_b_left,E2_inv_b_left,E3_inv_b_left,E4_inv_b_left)
DEALLOCATE(E1_inv_b_right,E2_inv_b_right,E3_inv_b_right,E4_inv_b_right)
DEALLOCATE(E1_inv_c_left,E2_inv_c_left,E3_inv_c_left,E4_inv_c_left)
DEALLOCATE(E1_inv_c_right,E2_inv_c_right,E3_inv_c_right,E4_inv_c_right)
DEALLOCATE(E1_inv_d_left,E2_inv_d_left,E3_inv_d_left,E4_inv_d_left)
DEALLOCATE(E1_inv_d_right,E2_inv_d_right,E3_inv_d_right,E4_inv_d_right)
DEALLOCATE(F1_inv_a_lower,F2_inv_a_lower,F3_inv_a_lower,F4_inv_a_lower)
DEALLOCATE(F1_inv_a_upper,F2_inv_a_upper,F3_inv_a_upper,F4_inv_a_upper)
DEALLOCATE(F1_inv_b_lower,F2_inv_b_lower,F3_inv_b_lower,F4_inv_b_lower)
DEALLOCATE(F1_inv_b_upper,F2_inv_b_upper,F3_inv_b_upper,F4_inv_b_upper)
DEALLOCATE(F1_inv_c_lower,F2_inv_c_lower,F3_inv_c_lower,F4_inv_c_lower)
DEALLOCATE(F1_inv_c_upper,F2_inv_c_upper,F3_inv_c_upper,F4_inv_c_upper)
DEALLOCATE(F1_inv_d_lower,F2_inv_d_lower,F3_inv_d_lower,F4_inv_d_lower)
DEALLOCATE(F1_inv_d_upper,F2_inv_d_upper,F3_inv_d_upper,F4_inv_d_upper)
DEALLOCATE(E1_vis_a_left,E2_vis_a_left,E3_vis_a_left,E4_vis_a_left)
DEALLOCATE(E1_vis_a_right,E2_vis_a_right,E3_vis_a_right,E4_vis_a_right)
DEALLOCATE(E1_vis_b_left,E2_vis_b_left,E3_vis_b_left,E4_vis_b_left)
DEALLOCATE(E1_vis_b_right,E2_vis_b_right,E3_vis_b_right,E4_vis_b_right)
DEALLOCATE(E1_vis_c_left,E2_vis_c_left,E3_vis_c_left,E4_vis_c_left)
DEALLOCATE(E1_vis_c_right,E2_vis_c_right,E3_vis_c_right,E4_vis_c_right)
DEALLOCATE(E1_vis_d_left,E2_vis_d_left,E3_vis_d_left,E4_vis_d_left)
DEALLOCATE(E1_vis_d_right,E2_vis_d_right,E3_vis_d_right,E4_vis_d_right)
DEALLOCATE(F1_vis_a_lower,F2_vis_a_lower,F3_vis_a_lower,F4_vis_a_lower)
DEALLOCATE(F1_vis_a_upper,F2_vis_a_upper,F3_vis_a_upper,F4_vis_a_upper)
DEALLOCATE(F1_vis_b_lower,F2_vis_b_lower,F3_vis_b_lower,F4_vis_b_lower)
DEALLOCATE(F1_vis_b_upper,F2_vis_b_upper,F3_vis_b_upper,F4_vis_b_upper)
DEALLOCATE(F1_vis_c_lower,F2_vis_c_lower,F3_vis_c_lower,F4_vis_c_lower)
DEALLOCATE(F1_vis_c_upper,F2_vis_c_upper,F3_vis_c_upper,F4_vis_c_upper)
DEALLOCATE(F1_vis_d_lower,F2_vis_d_lower,F3_vis_d_lower,F4_vis_d_lower)
DEALLOCATE(F1_vis_d_upper,F2_vis_d_upper,F3_vis_d_upper,F4_vis_d_upper)
DEALLOCATE(dU1dt_a,dU2dt_a,dU3dt_a,dU4dt_a)
DEALLOCATE(dU1dt_b,dU2dt_b,dU3dt_b,dU4dt_b)
DEALLOCATE(dU1dt_c,dU2dt_c,dU3dt_c,dU4dt_c)
DEALLOCATE(dU1dt_d,dU2dt_d,dU3dt_d,dU4dt_d)
DEALLOCATE(U1_new,U2_new,U3_new,U4_new)
DEALLOCATE(err1,err2,err3,err4)
DEALLOCATE(Cf,Cf_res,Ch,Ch_res)

End Subroutine KillArrays
