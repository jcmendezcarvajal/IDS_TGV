! Update the solution flux (Um)
Subroutine Update
  use Variables
  use omp_lib
  implicit none

  !$OMP PARALLEL PRIVATE(i,j)
  !$OMP DO
  do j =2,Jmax-1
    do i =2,Imax-1
      U1_new(i,j) = U1_old(i,j)+dU1dt(i,j)*delta_t
      U2_new(i,j) = U2_old(i,j)+dU2dt(i,j)*delta_t
      U3_new(i,j) = U3_old(i,j)+dU3dt(i,j)*delta_t
      U4_new(i,j) = 1.0 !U4_old(i,j)+dU4dt(i,j)*delta_t
    enddo
  enddo
  !$OMP END DO

  !calculate the new flow variables
  !$OMP DO
  do j =2,Jmax-1
    do i =2,Imax-1
      r_new(i,j) = U1_new(i,j)
      u_new(i,j) = U2_new(i,j)/U1_new(i,j)
      v_new(i,j) = U3_new(i,j)/U1_new(i,j)
      T_new(i,j) = 1.0 !((U4_new(i,j)/U1_new(i,j))-0.5*(u_new(i,j)**2+ &
      !v_new(i,j)**2))*(gamma*(gamma-1.0)*cM_inf**2)
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  Call BC_Supersonic_New


  ! Adjust the primative variables at the inner nodes
  !$OMP PARALLEL PRIVATE(i,j)
  !$OMP DO
  do j = 2,Jmax-1
    do i = 2,Imax-1
      r_temp(i,j)	=  1.0*r_new(i,j)/4.0 &
      +  1.0*(r_new(i,j-1)   + r_new(i,j+1) &
      +  r_new(i-1,j)   + r_new(i+1,j))/8.0 &
      +  1.0*(r_new(i-1,j-1) + r_new(i-1,j+1) &
      +  r_new(i+1,j-1) + r_new(i+1,j+1))/16.0

      u_temp(i,j)=    1.0*u_new(i,j)/4.0 &
      +  1.0*(u_new(i,j-1)   + u_new(i,j+1) &
      +  u_new(i-1,j)   + u_new(i+1,j))/8.0 &
      +  1.0*(u_new(i-1,j-1) + u_new(i-1,j+1) &
      +  u_new(i+1,j-1) + u_new(i+1,j+1))/16.0

      v_temp(i,j)=    1.0*v_new(i,j)/4.0 &
      +  1.0*(v_new(i,j-1)   + v_new(i,j+1) &
      +  v_new(i-1,j)   + v_new(i+1,j))/8.0 &
      +  1.0*(v_new(i-1,j-1) + v_new(i-1,j+1) &
      +  v_new(i+1,j-1) + v_new(i+1,j+1))/16.0

      T_temp(i,j)=    1.0*T_new(i,j)/4.0 &
      +  1.0*(T_new(i,j-1)   + T_new(i,j+1) &
      +  T_new(i-1,j)   + T_new(i+1,j))/8.0 &
      +  1.0*(T_new(i-1,j-1) + T_new(i-1,j+1) &
      +  T_new(i+1,j-1) + T_new(i+1,j+1))/16.0
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  ! Adjust the predicted boundary conditions
  Call BC_Supersonic_Temp



  ! Calculate the maximum error in the fluxes

  !$OMP PARALLEL PRIVATE(i,j)
  !$OMP DO
  do j = 2,Jmax-1
    do i = 2,Imax-1
      H_temp(i,j) = 1.0 !r_temp(i,j)*(T_temp(i,j)/alpha +0.5*(u_temp(i,j)**2+ &
      !v_temp(i,j)**2))
    enddo
  enddo
  !$OMP END DO

  !This section computes the RMS of the residuals

  !$OMP DO
  do j = 2,Jmax-1
    do i = 2,Imax-1
      err1(i,j) = (r_temp(i,j)-r_old(i,j))/delta_t
      err2(i,j) = (r_temp(i,j)*u_temp(i,j)-r_old(i,j)*u_old(i,j))/delta_t
      err3(i,j) = (r_temp(i,j)*v_temp(i,j)-r_old(i,j)*v_old(i,j))/delta_t
      err4(i,j) = (H_temp(i,j)-H_old(i,j))/delta_t
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  NormU1 = NORM2(err1(2:Imax-1,2:Jmax-1))/sqrt(real(Imax-2)*real(Jmax-2))
  NormU2 = NORM2(err2(2:Imax-1,2:Jmax-1))/sqrt(real(Imax-2)*real(Jmax-2))
  NormU3 = NORM2(err3(2:Imax-1,2:Jmax-1))/sqrt(real(Imax-2)*real(Jmax-2))
  NormU4 = NORM2(err4(2:Imax-1,2:Jmax-1))/sqrt(real(Imax-2)*real(Jmax-2))

  eps=max(NormU1, NormU2, NormU3, NormU4)


End Subroutine Update
