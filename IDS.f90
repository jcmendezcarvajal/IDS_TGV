! A program to the Navier_Stokes equations flate plate with ports at quarter markers of the plate
!
Program IDS
  use Variables
  use omp_lib
  implicit none
  open(32, file = 'error.dat')
  open(33, file = 'Enstrophy.dat')
  open(34, file = 'VorticityRMS.dat')
  open(35, file = 'U_RMS.dat')
  open(36, file = 'V_RMS.dat')
  open(37, file = 'EnstrophyTheoretical.dat')


  ! Initializing the Arrays
  Call InitializeArrays

  ! Defining the number of threads for OpenMP
  CALL omp_set_num_threads(NumThreads)

  Call Freestream
  iRestart    = 0
  AccumulatedTime = 0.0 ! This will be used to accumulate the time elapsed since dt is not constant

  !Settting up the frecuency for the temporal plot
  TimeToPrint = 0.1
  PrintFrecuency =  TimeToPrint  !Every nondimensional time will be printed the sol.

  if(iRestart.eq.0)then
    Call Grid
    Call WAG
    Call GridNonDimensional
  else
    Call Grid
    Call Read_SolutionRestart
  endif

  Call InitialCondition

  kk = 0
  call Initial_Enstrophy

  DO
    kk= kk + 1

    Call Flux_U
    Call Viscous_Properties
    Call Time_Step
    Call L_vis_in
    Call R_vis_in
    Call Lower_vis_in
    Call Upper_vis_in
    Call L_vis_bc
    Call R_vis_bc
    Call Lower_vis_bc
    Call Upper_vis_bc
    Call Invicid_X_fluxes
    Call Invicid_Y_fluxes
    Call Viscous_X_fluxes
    Call Viscous_Y_fluxes
    Call Derivatives
    Call Update

    If (mod(kk,1000)==0) then
      write(*,*) kk,eps,(kk*delta_t)
    endif
    Call Swap

    ! Calculating the accumulated time.
    AccumulatedTime = AccumulatedTime + delta_t

    Call Enstrophy_Computation
    CALL Analytical_Solution
    CALL Error


    ! Time check for temporal plot
    if(((PrintFrecuency-(AccumulatedTime))/PrintFrecuency).LT.1.0*10E-2) Call Transient_Primitive

  ! Checking Convergence or computational time.
    if ((AccumulatedTime).GE.1.0) then  !This represents the nondimensional time
      exit
    endif
  END DO

  Call KillArrays   ! Deallocating the arrays not needed for writting output
  Call Output
  Call Write_SolutionRestart

End Program IDS
