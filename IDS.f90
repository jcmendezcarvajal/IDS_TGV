! A program to the Navier_Stokes equations flate plate with ports at quarter markers of the plate
!
Program IDS
  use Variables
  use omp_lib
  implicit none
  open(32, file = 'error.dat')
  open(33, file = 'Total_KineticEnergy.dat')


  ! Initializing the Arrays
  Call InitializeArrays

  ! Defining the number of threads for OpenMP
  CALL omp_set_num_threads(NumThreads)

  Call Freestream
  iRestart    = 0

  !Settting up the frecuency for the temporal plot
  PrintFrecuency =  0.01  !Every nondimensional time will be printed the sol.

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
  call Initial_KineticEnergy
!  DO
!    kk= kk + 1
   do KK = 1,5000
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
    Call KineticEnergy_Computation

    ! Time check for temporal plot
!    if(((PrintFrecuency-(kk*delta_t))/PrintFrecuency).LT.1.0*10E-2) Call Transient_Primitive

  ! Checking Convergence or computational time.
    if ((kk*delta_t).GE.3) then  !This represents the nondimensional time
      exit
    endif
  END DO

  Call KillArrays   ! Deallocating the arrays not needed for writting output
  Call Output
  Call Write_SolutionRestart

End Program IDS
