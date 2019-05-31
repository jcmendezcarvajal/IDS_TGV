!This subroutine writes transient files
Subroutine Enstrophy_Computation
use Variables
Implicit none

!Computing Vorticity
!$OMP PARALLEL PRIVATE (i,j)
!$OMP DO
! Let's try with the actual value and not the magnitude of the vector
do j = 2,Jmax-1
    do i = 2,Imax-1
        Vorticity(i,j)	= ((v_old(i+1,j) - v_old(i-1,j))/(ddx(i)+ddx(i-1)))		&
        - ((u_old(i,j+1) - u_old(i,j-1))/(ddy(j)+ddy(j-1)))
    enddo
enddo
!$OMP END DO
!$OMP END PARALLEL

Max_Vorticity_numerical = MAXVAL(Vorticity(2:Imax-1,2:jmax-1))

Enstrophy = 0.0
!$OMP DO REDUCTION(+:Enstrophy)
do j = 2,Jmax-1
do i = 2,Imax-1
  Enstrophy = Enstrophy + &
  (1.0/2.0)*(Vorticity(i,j)*Vorticity(i,j))
enddo
enddo
!$OMP END DO

Enstrophy = Enstrophy/((Imax-2)*(Jmax-2))
Enstrophy = Enstrophy/Enstrophy_Initial

write(33,*) AccumulatedTime, Enstrophy

End Subroutine Enstrophy_Computation

Subroutine Initial_Enstrophy
use Variables
Implicit none
! Computing the initial Enstrophy
Enstrophy_Initial = 0.0

! Let's try with the actual value and not the magnitude of the vector
do j = 2,Jmax-1
    do i = 2,Imax-1
        Vorticity(i,j)	= ((v_old(i+1,j) - v_old(i-1,j))/(ddx(i)+ddx(i-1)))		&
        - ((u_old(i,j+1) - u_old(i,j-1))/(ddy(j)+ddy(j-1)))
    enddo
enddo

do j = 2,Jmax-1
do i = 2,Imax-1
  Enstrophy_Initial = Enstrophy_Initial + &
  (1.0/2.0)*(Vorticity(i,j)*Vorticity(i,j))
enddo
enddo

Enstrophy_Initial = Enstrophy_Initial/((Imax-2)*(Jmax-2))
write(33,*) 0.0, Enstrophy_Initial/Enstrophy_Initial

End Subroutine Initial_Enstrophy

Subroutine InitialCondition
use Variables
Implicit none
open (Unit = 131,file='Initial.dat')
write(131,*)'  zone T = "zone", I = ',Imax,' J= ',Jmax,' F = point'

do j = 1, Jmax
  do i = 1, Imax
      write(131,*) x(i,j),y(i,j),u_old(i,j),v_old(i,j),T_old(i,j),r_old(i,j)
  enddo
enddo
close(131)

End Subroutine InitialCondition


Subroutine Transient_Primitive
use Variables
Implicit none

WRITE(fileout,'(f10.5)') (PrintFrecuency)
NAME = trim('Transient_Primitive')//'_'//trim(fileout)

open(62,file = trim(NAME)//'.dat')

write(62,*) 'TITLE = "Transient Output"'
write(62,*) 'VARIABLES = "X"'
write(62,*) '"Y"'
write(62,*) '"U Velocity"'
write(62,*) '"V Velocity"'
write(62,*) '"Temperature"'
write(62,*) '"Density"'
write(62,*) '" Vorticity "'
write(62,*) '" Stream Function "'
write(62,*)'  zone T = "zone", I = ',Imax,' J= ',Jmax,' F = point'

!Computing the stream function variable Psi
Psi(1,1)    = 0.0
do j = 2,Jmax
    Psi(1,j)    = Psi(1,j-1) - u_old(1,j)*ddy(j)
enddo
do j = 2, Jmax
    do i = 2, Imax
        Psi(i,j)    = Psi(i-1,j) + v_old(i,j)*ddx(i)
    enddo
enddo

!Computing Vorticity
!$OMP PARALLEL PRIVATE (i,j)
!$OMP DO
! Let's try with the actual value and not the magnitude of the vector
do j = 2,Jmax-1
    do i = 2,Imax-1
        Vorticity(i,j)	= ((v_old(i+1,j) - v_old(i-1,j))/(ddx(i)+ddx(i-1)))		&
        - ((u_old(i,j+1) - u_old(i,j-1))/(ddy(j)+ddy(j-1)))
    enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

Max_Vorticity_numerical = MAXVAL(Vorticity(2:Imax-1,2:jmax-1))

Vorticity(1,   2:Jmax-1)	= Vorticity(Imax-1,     2:Jmax-1)
Vorticity(Imax,2:Jmax-1)	= Vorticity(2,2:Jmax-1)
Vorticity(1:Imax, 1    )	= Vorticity(1:Imax, Jmax-1     )
Vorticity(1:Imax, Jmax )	= Vorticity(1:Imax,  2)


do j = 1, Jmax
  do i = 1, Imax
    write(62,*) x(i,j),y(i,j),u_old(i,j),v_old(i,j),T_old(i,j),r_old(i,j),Vorticity(i,j)/Max_Vorticity_numerical,Psi(i,j)
  enddo
enddo

close(62)
PrintFrecuency =  TimeToPrint + PrintFrecuency

! Closing and opening the files again, to have access to such files more often.

close(33)
close(34)
close(35)
close(36)
close(37)
!Opening the files again

open(33, file = 'Enstrophy.dat',Access = 'append',Status='old')
open(34, file = 'VorticityRMS.dat',Access = 'append',Status='old')
open(35, file = 'U_RMS.dat',Access = 'append',Status='old')
open(36, file = 'V_RMS.dat',Access = 'append',Status='old')
open(37, file = 'EnstrophyTheoretical.dat',Access = 'append',Status='old')

End Subroutine Transient_Primitive

subroutine Analytical_Solution
Use variables
Implicit None

!Dimensional
F = -2.0*((AccumulatedTime*(x_actual/u_inf))*vNu_inf)
!$OMP PARALLEL PRIVATE (i,j)
!$OMP DO
do j = 2,Jmax-1
    do i = 2,Imax-1
        vorticity_exact(i,j) = (2.0*SIN(xDim(i,j))*SIN(yDim(i,j)))&
        *EXP(F*vNu(i,j))

        psi_exact(i,j)       = (SIN(xDim(i,j))*SIN(yDim(i,j)))&
        *EXP(F*vNu(i,j))

        u_exact(i,j)         = (SIN(xDim(i,j))*COS(yDim(i,j)))&
        *EXP(F*vNu(i,j))

        v_exact(i,j)         = (-COS(xDim(i,j))*SIN(yDim(i,j)))&
        *EXP(F*vNu(i,j))
    end do
end do
!$OMP END DO
!$OMP END PARALLEL

Enstrophy = 0.0
!$OMP DO REDUCTION(+:Enstrophy)
do j = 2,Jmax-1
do i = 2,Imax-1
  Enstrophy = Enstrophy + &
  (1.0/2.0)*(vorticity_exact(i,j)*vorticity_exact(i,j))
enddo
enddo
!$OMP END DO

Enstrophy = Enstrophy/((Imax-2)*(Jmax-2))
Enstrophy = Enstrophy/Enstrophy_Initial

!Writting the Theoretical Enstrophy
write(37,*) AccumulatedTime, Enstrophy


uexact_inf =  MAXVAL(u_exact(2:imax-1,2:Jmax-1))
vorticity_inf = MAXVAL(vorticity_exact(2:Imax-1,2:jmax-1))
psi_inf = MAXVAL(psi_exact(2:Imax-1,2:jmax-1))

!Non_Dimensional
!$OMP PARALLEL PRIVATE (i,j)
!$OMP DO
	do j = 2,Jmax-1
	do i = 2,Imax-1
		u_exact(i,j) =  u_exact(i,j)/uexact_inf
		v_exact(i,j) = 	v_exact(i,j)/uexact_inf
	  vorticity_exact(i,j) = vorticity_exact(i,j)/vorticity_inf
    psi_exact(i,j) = psi_exact(i,j)/psi_inf
	enddo
	enddo
  !$OMP END DO
  !$OMP END PARALLEL

!Vorticity_exact Boundary
vorticity_exact(1,   2:Jmax-1)	= vorticity_exact(Imax-1,     2:Jmax-1) !Left wall
vorticity_exact(Imax,2:Jmax-1)	= vorticity_exact(2,2:Jmax-1)           !Right wall
vorticity_exact(1:Imax, 1    )	= vorticity_exact(1:Imax, Jmax-1     )  !Bottom Wall
vorticity_exact(1:Imax, Jmax )	= vorticity_exact(1:Imax,  2)           !Top Wall
!Psi_exact Boundary
psi_exact(1,   2:Jmax-1)	= psi_exact(Imax-1,     2:Jmax-1)
psi_exact(Imax,2:Jmax-1)	= psi_exact(2,2:Jmax-1)
psi_exact(1:Imax, 1    )	= psi_exact(1:Imax, Jmax-1     )
psi_exact(1:Imax, Jmax )	= psi_exact(1:Imax,  2)
!u_exact boundary
u_exact(1,   2:Jmax-1)	= u_exact(Imax-1,     2:Jmax-1)
u_exact(Imax,2:Jmax-1)	= u_exact(2,2:Jmax-1)
u_exact(1:Imax, 1    )	= u_exact(1:Imax, Jmax-1     )
u_exact(1:Imax, Jmax )	= u_exact(1:Imax,  2)
!v_exact boundary
v_exact(1,   2:Jmax-1)	= v_exact(Imax-1,     2:Jmax-1)
v_exact(Imax,2:Jmax-1)	= v_exact(2,2:Jmax-1)
v_exact(1:Imax, 1    )	= v_exact(1:Imax, Jmax-1     )
v_exact(1:Imax, Jmax )	= v_exact(1:Imax,  2)


if(((PrintFrecuency-(AccumulatedTime))/PrintFrecuency).LT.1.0*10E-2) then

  WRITE(fileout,'(f10.5)') (PrintFrecuency)
  NAME = trim('Transient_Properties')//'_'//trim(fileout)

  open(64,file = trim(NAME)//'.dat')

  write(64,*) 'TITLE = "Transient Output"'
  write(64,*) 'VARIABLES = "X"'
  write(64,*) '"Y"'
  write(64,*) '"Vorticity Exact"'
  write(64,*) '"Psi exact"'
  write(64,*) '"U exact"'
  write(64,*) '"V exact "'
  write(64,*)'  zone T = "zone", I = ',Imax,' J= ',Jmax,' F = point'

  do j = 1, Jmax
    do i = 1, Imax
      write(64,*) x(i,j),y(i,j),vorticity_exact(i,j),psi_exact(i,j),u_exact(i,j),v_exact(i,j)
    enddo
  enddo

  close(64)
end if

end subroutine Analytical_Solution

Subroutine Error

Use variables
Implicit None
!
!$OMP PARALLEL PRIVATE (i,j)
!$OMP DO
do j = 2, Jmax-1
    do i= 2, Imax-1
        Vor_err(i,j) = ABS(vorticity_exact(i,j)-Vorticity(i,j)/Max_Vorticity_numerical)
        u_err(i,j)   = ABS(u_exact(i,j)-u_old(i,j))
        v_err(i,j)   = ABS(v_exact(i,j)-v_old(i,j))
    end do
end do
!$OMP END DO
!$OMP END PARALLEL

    Vor_erms = SQRT(SUM(Vor_err(2:imax-1,2:Jmax-1)**2)/(Imax-2)**2)
    u_erms   = SQRT(SUM(u_err(2:imax-1,2:Jmax-1)**2)/(Imax-2)**2)
    v_erms   = SQRT(SUM(v_err(2:imax-1,2:Jmax-1)**2)/(Imax-2)**2)

    !!!!Storing Result

    write(34,*) AccumulatedTime, Vor_erms
    write(35,*) AccumulatedTime, u_erms
    write(36,*) AccumulatedTime, v_erms

end subroutine Error
