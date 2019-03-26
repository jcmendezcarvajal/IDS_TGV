! Generate the output
Subroutine Output
    use Variables
    implicit none

    ! Formating section

    398   format(6(e22.15,2x))
    399   format(7(e22.15,2x))
    39    format(5(e22.15,2x))

    !!! Computing Cp for postprocessing purposes
    Cp = gamma*R_gas/(gamma-1)


    !! ****************************************************************************** !!
    !Printing Analytical Solution

    open (113, file = 'Analytical_Solutions_Final.dat')
    write(113,*) 'TITLE = "Analytical_Solutions"'
    write(113,*) 'VARIABLES = "X"'
    write(113,*) '"Y"'
    write(113,*) '"U Velocity Exact"'
    write(113,*) '"V Velocity Exact"'
    write(113,*) '"Vorticity Exact"'
    write(113,*) '"Stream Function Exact"'

    write(113,*)'  zone T = "zone1", I = ',Imax,' J= ',Jmax,' F = point'

    do j = 1,Jmax
        do i = 1,Imax
            write(113,398) x(i,j),y(i,j),u_exact(i,j),v_exact(i,j),vorticity_exact(i,j),psi_exact(i,j)
        enddo
    enddo
    close(113)




    !! ****************************************************************************** !!
    ! Writing the Primitive Variables

    open(330, file = 'Primitive_Variables.dat')
    write(330,*) 'TITLE = "Primitive Variables Contours"'
    write(330,*) 'VARIABLES = "X"'
    write(330,*) '"Y"'
    write(330,*) '"U Velocity"'
    write(330,*) '"V Velocity"'
    write(330,*) '"Density"'
    write(330,*) '"Temperature"'

    write(330,*)'  zone T = "zone1", I = ',Imax,' J= ',Jmax,' F = point'

    do j = 1,Jmax
        do i = 1,Imax
            write(330,398) x(i,j),y(i,j),u_old(i,j),v_old(i,j),r_old(i,j),T_old(i,j)
        enddo
    enddo
    close(330)

    !! ****************************************************************************** !!
    !Writing the Compound variables; time fluxes.

    open(331, file = 'Time_Fluxes.dat')
    write(331,*) 'TITLE = "Time Fluxes"'
    write(331,*) 'VARIABLES = "X"'
    write(331,*) '"Y"'
    write(331,*) '"dU1dt"'
    write(331,*) '"dU2dt"'
    write(331,*) '"dU3dt"'
    write(331,*) '"dU4dt"'

    write(331,*)'  zone T = "zone1", I = ',Imax,' J= ',Jmax,' F = point'

    do j = 1,Jmax
        do i = 1,Imax
            write(331,398) x(i,j),y(i,j),dU1dt(i,j),dU2dt(i,j),dU3dt(i,j),dU4dt(i,j)
        enddo
    enddo
    close(331)

    !! ****************************************************************************** !!
    !Writing more Compound variables: Stream function, Vorticity, Mach (resultant) and Mach vector
    allocate (MachRe(Imax,Jmax))
    allocate (MachHo(Imax,Jmax))
    allocate (MachVe(Imax,Jmax))
  !  allocate (Psi(Imax,Jmax))
  !  allocate (Vorticity(Imax,Jmax))

    !Computing the Mach Number, based on the resultant velocity
    MachRe(:,:)= sqrt((u_old(:,:)*u_inf)**2 + (v_old(:,:)*u_inf)**2)/sqrt(gamma*R_gas*T_old(:,:)*T_inf)

    !Computing the Mach Number in the horizontal and vertical direction, thus it is a vector.
    MachHo(:,:) = u_old(:,:)*u_inf/sqrt(gamma*R_gas*T_old(:,:)*T_inf)
    MachVe(:,:) = v_old(:,:)*u_inf/sqrt(gamma*R_gas*T_old(:,:)*T_inf)

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

    Max_Vorticity = MAXVAL(Vorticity(1:Imax,1:jmax))

    open(332, file = 'CompoundData_I.dat')
    write(332,*) 'TITLE = "Compound Physics Quantities"'
    write(332,*) 'VARIABLES = "X"'
    write(332,*) '"Y"'
    write(332,*) '"Stream Function"'
    write(332,*) '"Vorticity"'
    write(332,*) '"Mach (resultant)"'
    write(332,*) '"Mach Horizontal (i)"'
    write(332,*) '"Mach Vertical   (j)"'

    write(332,*)'  zone T = "zone1", I = ',Imax,' J= ',Jmax,' F = point'

    do j = 1,Jmax
        do i = 1,Imax
            write(332,399) x(i,j),y(i,j),Psi(i,j),Vorticity(i,j)/Max_Vorticity,MachRe(i,j),MachHo(i,j),MachVe(i,j)
        enddo
    enddo
    close(332)
    deallocate (Psi)

    !! ****************************************************************************** !!
    !Writing more Compound variables: Pressure, Total properties, Delta Entropy and Delta Enthalpy
    ! delta entropy is normalized by the maximum entropy values whereas enthalpy is nondimensionalized by
    ! dividing  delta enthalpy by the enthalpy of the freestream.


    !allocate (Pres(Imax,Jmax))
    allocate (P_total(Imax,Jmax))
    allocate (T_total(Imax,Jmax))
    allocate (DeltaS(Imax,Jmax))
    allocate (DeltaH(Imax,Jmax))

    Pres(:,:) = r_old(:,:)*T_old(:,:)
    T_total(:,:) = (1 + ((gamma-1.0)/2.0)*MachRe(:,:)**2.0)
    P_total(:,:) = (1 + ((gamma-1.0)/2.0)*MachRe(:,:)**2.0)**(gamma/(gamma-1))
    DeltaS(:,:) = Cp*LOG((T_old(:,:)*T_inf)/T_inf) - R_gas*LOG((Pres(:,:)*P_inf)/P_inf)
    DeltaH(:,:) = (Cp*(T_old(:,:)*T_inf - T_inf))/(Cp*T_inf)
    MaxDeltaS = abs(MAXVAL(DeltaS))

    open(333, file = 'CompoundData_II.dat')
    write(333,*) 'TITLE = "Compound Physics Quantities"'
    write(333,*) 'VARIABLES = "X"'
    write(333,*) '"Y"'
    write(333,*) '"Pressure"'
    write(333,*) '"Total Pressure"'
    write(333,*) '"Total Temperature"'
    write(333,*) '"Delta Entropy"'
    write(333,*) '"Delta Enthalpy"'

    write(333,*)'  zone T = "zone1", I = ',Imax,' J= ',Jmax,' F = point'

    do j = 1,Jmax
        do i = 1,Imax
            write(333,399) x(i,j),y(i,j),Pres(i,j),P_total(i,j),T_total(i,j),DeltaS(i,j)/MaxDeltaS,DeltaH(i,j)
        enddo
    enddo
    close(333)

    P_TotalFree = (1 + ((gamma-1.0)/2.0)*cM_inf**2.0)**(gamma/(gamma-1))
    DeltaS(:,:) = -1.0*R_gas*LOG((P_total(:,:))/P_TotalFree)

    open(333, file = 'CompoundData_III.dat')
    write(333,*) 'TITLE = "Compound Physics Quantities"'
    write(333,*) 'VARIABLES = "X"'
    write(333,*) '"Y"'
    write(333,*) '"Efficiency"'
    write(333,*) '"Delta Entropy"'

    write(333,*)'  zone T = "zone1", I = ',Imax,' J= ',Jmax,' F = point'

    do j = 1,Jmax
        do i = 1,Imax
            write(333,39) x(i,j),y(i,j),P_total(i,j)/P_TotalFree,DeltaS(i,j)
        enddo
    enddo
    close(333)

    deallocate (P_total, T_total, DeltaS, DeltaH)

    !! ****************************************************************************** !!!!
    !Computing flow feature extraction functions.

    allocate (GradRho(Imax,Jmax))
    allocate (NGradRho(Imax,Jmax))
    allocate (threshold(Imax,Jmax))
    allocate (MagPresGrad(Imax,Jmax))
    allocate (MachNormal(Imax,Jmax))
    allocate (Enthalpy_Stag(Imax,Jmax))
    allocate (GradEntroI(Imax,Jmax))
    allocate (GradEntroJ(Imax,Jmax))
    allocate (MagGradEntropy(Imax,Jmax))
    allocate (Q_crit(Imax,Jmax))

    !Computing the gradient of density; The printed variable is the magnitude of the gradient.
    do j = 2,Jmax-1
        do i = 2,Imax-1
            GradRho(i,j)	= sqrt( ((r_old(i+1,j) - r_old(i-1,j))/(ddx(i)+ddx(i-1)))**2 &
            +((r_old(i,j+1) - r_old(i,j-1))/(ddy(j)+ddy(j-1)))**2)
        enddo
    enddo
    GradRho(1,   2:Jmax-1)	= GradRho(2,     2:Jmax-1)
    GradRho(Imax,2:Jmax-1)	= GradRho(Imax-1,2:Jmax-1)
    GradRho(1:Imax, 1    )	= GradRho(1:Imax, 2      )
    GradRho(1:Imax, Jmax )	= GradRho(1:Imax,  Jmax-1)

    Max_GradRho = GradRho(5,5)
    do j = 1,Jmax
        do i = 1,Imax
            Max_GradRho = max(Max_GradRho, GradRho(i,j))
        enddo
    enddo

    !Computing the normal gradient of density. This computes the gradient of density in the direction
    ! of the velocity. This was proposed by Pagendarm and Seits 1993.
    do j = 2,Jmax-1
        do i = 2,Imax-1
            NGradRho(i,j) = (((r_old(i+1,j) - r_old(i-1,j))/(ddx(i)+ddx(i-1)))        &
            *   (u_old(i,j)/sqrt((u_old(i,j))**2+(v_old(i,j))**2))      &
            +  ((r_old(i,j+1) - r_old(i,j-1))/(ddy(j)+ddy(j-1)))        &
            *   (v_old(i,j)/sqrt((u_old(i,j))**2+(v_old(i,j))**2)) )
        enddo
    enddo
    NGradRho(1,   2:Jmax-1)	= NGradRho(2,     2:Jmax-1)
    NGradRho(Imax,2:Jmax-1)	= NGradRho(Imax-1,2:Jmax-1)
    NGradRho(1:Imax, 1    )	= NGradRho(1:Imax, 2      )
    NGradRho(1:Imax, Jmax )	= NGradRho(1:Imax,  Jmax-1)

    !Computing the normal Mach number. It computes the Mach number normal to the shock, whose value is 1
    !This technique was proposed by Lovely and Haimes
    do j = 2,Jmax-1
        do i = 2,Imax-1
            threshold(i,j) = 0.007*((((Pres(i,j)*P_inf) + P_inf)/P_inf)/ &
            (sqrt((ddx(i)+ddx(i-1))**2 + (ddy(j)+ddy(j-1))**2)))
        enddo
    enddo

    do j = 2,Jmax-1
        do i = 2,Imax-1
            MagPresGrad(i,j) = sqrt((((Pres(i+1,j) - Pres(i-1,j))/(ddx(i)+ddx(i-1))))**2 &
            + (((Pres(i,j+1) - Pres(i,j-1))/(ddy(j)+ddy(j-1))))**2)
        enddo
    enddo

    do j = 2,Jmax-1
        do i = 2,Imax-1
            MachNormal(i,j) = cM_inf*(((MachHo(i,j)*((Pres(i+1,j) - Pres(i-1,j))/(ddx(i)+ddx(i-1)))) &
            + (MachVe(i,j)*((Pres(i,j+1) - Pres(i,j-1))/(ddy(j)+ddy(j-1))))) &
            /max(MagPresGrad(i,j),threshold(i,j)))
        enddo
    enddo

    MachNormal(1,   2:Jmax-1)	= MachNormal(2,     2:Jmax-1)
    MachNormal(Imax,2:Jmax-1)	= MachNormal(Imax-1,2:Jmax-1)
    MachNormal(1:Imax, 1    )	= MachNormal(1:Imax, 2      )
    MachNormal(1:Imax, Jmax )	= MachNormal(1:Imax,  Jmax-1)

    ! This section computes the magnitude of the entropy gradient.
    !Computing Enthalpy (cp= cte; calorically perfect gas)
    Enthalpy_Stag(:,:) = Cp*T_old(:,:)+((sqrt((u_old(:,:))**2+(v_old(:,:))**2))**2)/2.0

    !Computing the gradient of entropy based on Croccos Theorem

    do j = 2,Jmax-1
        do i = 2,Imax-1
            ! I direction
            GradEntroI(i,j) = ((v_old(i,j)*Vorticity(i,j))/T_old(i,j)) - &
            ((Enthalpy_Stag(i+1,j) - Enthalpy_Stag(i-1,j)) &
            /((ddx(i)+ddx(i-1))*T_old(i,j)))
            ! J direction
            GradEntroJ(i,j) = ((u_old(i,j)*Vorticity(i,j))/T_old(i,j)) - &
            ((Enthalpy_Stag(i,j+1) - Enthalpy_Stag(i,j-1)) &
            /((ddy(j)+ddy(j-1))*T_old(i,j)))
        enddo
    enddo

    GradEntroI(1,   2:Jmax-1)	= GradEntroI(2,     2:Jmax-1)
    GradEntroI(Imax,2:Jmax-1)	= GradEntroI(Imax-1,2:Jmax-1)
    GradEntroI(1:Imax, 1    )	= GradEntroI(1:Imax, 2      )
    GradEntroI(1:Imax, Jmax )	= GradEntroI(1:Imax,  Jmax-1)

    GradEntroJ(1,   2:Jmax-1)	= GradEntroJ(2,     2:Jmax-1)
    GradEntroJ(Imax,2:Jmax-1)	= GradEntroJ(Imax-1,2:Jmax-1)
    GradEntroJ(1:Imax, 1    )	= GradEntroJ(1:Imax, 2      )
    GradEntroJ(1:Imax, Jmax )	= GradEntroJ(1:Imax,  Jmax-1)


    !Computing the magnitude of the entropy gradient

    MagGradEntropy(:,:) = sqrt((GradEntroI(:,:))**2 + (GradEntroJ(:,:))**2 )

    ! This section computes the Q criterion
    do j = 2,Jmax-1
        do i = 2,Imax-1
            S11 = 0.5*(((u_old(i+1,j) - u_old(i-1,j))/(ddx(i)+ddx(i-1))) + &
            ((u_old(i+1,j) - u_old(i-1,j))/(ddx(i)+ddx(i-1))))
            S12 = 0.5*(((u_old(i,j+1) - u_old(i,j-1))/(ddy(j)+ddy(j-1))) + &
            ((v_old(i+1,j) - v_old(i-1,j))/(ddx(i)+ddx(i-1))))
            S22 = 0.5*(((v_old(i,j+1) - v_old(i,j-1))/(ddy(j)+ddy(j-1))) + &
            ((v_old(i,j+1) - v_old(i,j-1))/(ddy(j)+ddy(j-1))))

            W11= 0.0
            W12 = 0.5*(((u_old(i,j+1) - u_old(i,j-1))/(ddy(j)+ddy(j-1))) - &
            ((v_old(i+1,j) - v_old(i-1,j))/(ddx(i)+ddx(i-1))))
            W22= 0.0

            Q_crit(i,j) = 0.5*((W11*W11+W12*W12+W12*W12+W22*W22) - &
            (S11*S11+S12*S12+S12*S12+S22*S22))
        enddo
    enddo
    Q_crit(1,   2:Jmax-1)	= Q_crit(2,     2:Jmax-1)
    Q_crit(Imax,2:Jmax-1)	= Q_crit(Imax-1,2:Jmax-1)
    Q_crit(1:Imax, 1    )	= Q_crit(1:Imax, 2      )
    Q_crit(1:Imax, Jmax )	= Q_crit(1:Imax,  Jmax-1)

    open(334, file = 'FlowFeaturesExtraction_Variables.dat')
    write(334,*) 'TITLE = "Flow Feature Extraction Variables"'
    write(334,*) 'VARIABLES = "X"'
    write(334,*) '"Y"'
    write(334,*) '"Density Gradient"'
    write(334,*) '"Normal Density Gradient"'
    write(334,*) '"Magnitude of Entropy Gradient"'
    write(334,*) '"Normal Mach number"'
    write(334,*) '"Q criterion"'


    write(334,*)'  zone T = "zone1", I = ',Imax,' J= ',Jmax,' F = point'

    do j = 1,Jmax
        do i = 1,Imax
            write(334,399) x(i,j),y(i,j),GradRho(i,j),NGradRho(i,j),MagGradEntropy(i,j),MachNormal(i,j),Q_crit(i,j)
        enddo
    enddo
    close(334)



    deallocate (GradRho,NGradRho,threshold,MagPresGrad,MachNormal,Enthalpy_Stag,GradEntroI,GradEntroJ)
    deallocate (MagGradEntropy, Q_crit)



End Subroutine Output
