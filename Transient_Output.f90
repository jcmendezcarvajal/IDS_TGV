!This subroutine writes transient files
Subroutine KineticEnergy_Computation
use Variables
Implicit none

KineticEnergy = 0.0
do j = 2,Jmax-1
do i = 2,Imax-1
  KineticEnergy = KineticEnergy + &
  (1.0/2.0)*(u_old(i,j)**2.0 + v_old(i,j)**2.0)
enddo
enddo

KineticEnergy = KineticEnergy/((Imax-2)*(Jmax-2))
KineticEnergy = KineticEnergy/KineticEnergy_Initial

write(33,*) (kk*delta_t), KineticEnergy

End Subroutine KineticEnergy_Computation

Subroutine Initial_KineticEnergy
use Variables
Implicit none
! Computing the initial Kinetic Energy
KineticEnergy_Initial = 0.0

do j = 2,Jmax-1
do i = 2,Imax-1
  KineticEnergy_Initial = KineticEnergy_Initial + &
  (1.0/2.0)*(u_old(i,j)**2.0 + v_old(i,j)**2.0)
enddo
enddo

KineticEnergy_Initial = KineticEnergy_Initial/((Imax-2)*(Jmax-2))
write(33,*) 0.0, KineticEnergy_Initial/KineticEnergy_Initial

End Subroutine Initial_KineticEnergy

Subroutine InitialCondition
use Variables
Implicit none
open (Unit = 131,file='Initial.dat')
write(131,*)'  zone T = "zone", I = ',Imax,' J= ',Jmax,' F = point'

do j = 1, Jmax
  do i = 1, Imax
      write(131,*) x(i,j),y(i,j),u_old(i,j),v_old(i,j),T_old(i,j)
  enddo
enddo
close(131)

End Subroutine InitialCondition


Subroutine Transient_Primitive
use Variables
Implicit none

WRITE(fileout,'(f7.5)') (PrintFrecuency)
NAME = trim('Transient_Primitive')//'_'//trim(fileout)

open(62,file = trim(NAME)//'.dat')

write(62,*) 'TITLE = "Transient Output"'
write(62,*) 'VARIABLES = "X"'
write(62,*) '"Y"'
write(62,*) '"U Velocity"'
write(62,*) '"V Velocity"'
write(62,*) '"Temperature"'
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
do j = 2,Jmax-1
    do i = 2,Imax-1
        Vorticity(i,j)	= sqrt(((v_old(i+1,j) - v_old(i-1,j))/(ddx(i)+ddx(i-1)))**2		&
        + ((u_old(i,j+1) - u_old(i,j-1))/(ddy(j)+ddy(j-1)))**2 )
    enddo
enddo

Vorticity(1,   2:Jmax-1)	= Vorticity(Imax-1,     2:Jmax-1)
Vorticity(Imax,2:Jmax-1)	= Vorticity(2,2:Jmax-1)
Vorticity(1:Imax, 1    )	= Vorticity(1:Imax, Jmax-1     )
Vorticity(1:Imax, Jmax )	= Vorticity(1:Imax,  2)


do j = 1, Jmax
  do i = 1, Imax
    write(62,*) x(i,j),y(i,j),u_old(i,j),v_old(i,j),T_old(i,j),Vorticity(i,j),Psi(i,j)
  enddo
enddo

close(62)
PrintFrecuency =  0.001 + PrintFrecuency

End Subroutine Transient_Primitive




subroutine Analytical_Solution
Use variables
Implicit None





!Dimensional
do j = 2,Jmax-1
    do i = 2,Imax-1
    
        vorticity_exact(i,j) = (2.0*SIN(xDim(i,j))*SIN(yDim(i,j))*EXP(-2.0*(kk*delta_t)*vNu(i,j)/vNu_inf)) !/vorticity_inf
        psi_exact(i,j)       = (SIN(xDim(i,j))*SIN(yDim(i,j))    *EXP(-2.0*(kk*delta_t)*vNu(i,j)/vNu_inf)) !/psi_inf
        u_exact(i,j)         = (SIN(xDim(i,j))*COS(yDim(i,j))    *EXP(-2.0*(kk*delta_t)*vNu(i,j)/vNu_inf)) !/u_inf 
        v_exact(i,j)         = (-COS(xDim(i,j))*SIN(yDim(i,j))   *EXP(-2.0*(kk*delta_t)*vNu(i,j)/vNu_inf)) !/u_inf


    end do
end do

uexact_inf =  MAXVAL(u_exact(2:imax-1,2:Jmax-1))
!vorticity_inf = MAXVAL(vorticity_exact(2:Imax-1,2:jmax-1))
!psi_inf = MAXVAL(psi_exact(2:Imax-1,2:jmax-1))

!Non_Dimensional
	do j = 2,Jmax-1
	do i = 2,Imax-1

		u_exact(i,j) =  u_exact(i,j)/uexact_inf
		v_exact(i,j) = 	v_exact(i,j)/uexact_inf
	    vorticity_exact(i,j) = vorticity_exact(i,j)/vorticity_inf
        psi_exact(i,j) = psi_exact(i,j)/psi_inf
        
	enddo
	enddo




    
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



end subroutine Analytical_Solution

subroutine Error

Use variables
Implicit None
!
do j = 1, Jmax
    do i= 1, Imax
    
        Vor_err(i,j) = ABS(vorticity_exact(i,j)-Vorticity(i,j))
        u_err(i,j)   = ABS(u_exact(i,j)-u_old(i,j))
        v_err(i,j)   = ABS(v_exact(i,j)-v_old(i,j))
        
        
    end do
end do

    Vor_erms = SQRT(SUM(Vor_err**2)/Imax**2)
    u_erms   = SQRT(SUM(u_err**2)/Imax**2)
    v_erms   = SQRT(SUM(v_err**2)/Imax**2)

!!!!Storing Result

    
    write(991,*) (kk*delta_t), Vor_erms
    write(992,*) (kk*delta_t), u_erms
    write(993,*) (kk*delta_t), v_erms
    



end subroutine Error

