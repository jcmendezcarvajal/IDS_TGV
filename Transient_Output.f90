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
write(62,*)'  zone T = "zone", I = ',Imax,' J= ',Jmax,' F = point'
write(62,*) 'TITLE = "Transient Output"'
write(62,*) 'VARIABLES = "X"'
write(62,*) '"Y"'
write(62,*) '"U Velocity"'
write(62,*) '"V Velocity"'
write(62,*) '"Temperature"'
write(62,*) '" Vorticity "'


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
    write(62,*) x(i,j),y(i,j),u_old(i,j),v_old(i,j),T_old(i,j),Vorticity(i,j)
  enddo
enddo

close(62)
PrintFrecuency =  0.001 + PrintFrecuency

End Subroutine Transient_Primitive
