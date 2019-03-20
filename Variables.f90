Module variables
    Implicit none
    ! Start the Initialization Block
    integer, parameter :: dp = 8
    INTEGER,PARAMETER::Imax=2001,Jmax=2001
    INTEGER::I,J,kk, iRestart
    INTEGER:: I_0Q, I_1Q,I_2Q, I_3Q, I_4Q

    !OpenMP Variables
    INTEGER,PARAMETER :: NumThreads = 32

    ! Freestream Input and calculated variables
    real(kind=dp)::r_inf,u_inf,T_inf,P_inf,vMu_inf,TwTinf,Pr_inf,cM_inf,a_inf
    real(kind=dp)::Re_L,gamma,alpha,T_ref,C_ref,Ce,eps,epsGlobal,R_gas,Cp
    real(kind=dp)::dx,dy, BetaG, Eta
    real(kind=dp)::delta_t,delta_tGlobal,xL,yL, x_actual, Leading_Edge_Gap ,delta, x_min, x_max
    ! Local Mesh Variables
    real(kind=dp),dimension(:,:),allocatable ::x,y,Re_cell
    real(kind=dp),dimension(:), allocatable ::xx, ddx
    real(kind=dp),dimension(:), allocatable ::yy, ddy
    ! Global Mesh Variables
    real(kind=dp),dimension(:,:),allocatable ::xGlobal,yGlobal
    real(kind=dp),dimension(:), allocatable ::xxGlobal, ddxGlobal
    real(kind=dp),dimension(:), allocatable ::yyGlobal, ddyGlobal
    ! Set an initial guess for the unknown primitive variables
    real(kind=dp),dimension(:,:),allocatable::r_old,u_old,v_old,T_old,a_old
    real(kind=dp),dimension(:,:),allocatable::r_temp,u_temp,v_temp,T_temp
    real(kind=dp),dimension(:,:),allocatable::r_new,u_new,v_new,T_new
    ! Set the Global Variables for the Boundary conditions
    real(kind=dp),dimension(:),allocatable::xxGlobalBC
    ! The Time-flux Vector_U_old
    real(kind=dp),dimension(:,:),allocatable::U1_old,U2_old,U3_old,U4_old,H_old,H_temp
    ! Calculate the viscous properties
    real(kind=dp),dimension(:,:),allocatable::vMu,tK,vNu
    ! Calculate the time_step
    real(kind=dp),dimension(:,:),allocatable::small_time,tNu_max,t1,t2
    ! Shear and Normal stresses(Tauxy),(Tauxx),(Tauyy)
    real(kind=dp),dimension(:,:),allocatable::Tauxy_a_lower,Tauxy_a_upper,Tauxy_a_left,Tauxy_a_right
    real(kind=dp),dimension(:,:),allocatable::Tauxx_a_left,Tauxx_a_right,Tauyy_a_lower,Tauyy_a_upper
    real(kind=dp),dimension(:,:),allocatable::Tauxy_b_lower,Tauxy_b_upper,Tauxy_b_left,Tauxy_b_right
    real(kind=dp),dimension(:,:),allocatable::Tauxx_b_left,Tauxx_b_right,Tauyy_b_lower,Tauyy_b_upper
    real(kind=dp),dimension(:,:),allocatable::Tauxy_c_lower,Tauxy_c_upper,Tauxy_c_left,Tauxy_c_right
    real(kind=dp),dimension(:,:),allocatable::Tauxx_c_left,Tauxx_c_right,Tauyy_c_lower,Tauyy_c_upper
    real(kind=dp),dimension(:,:),allocatable::Tauxy_d_lower,Tauxy_d_upper,Tauxy_d_left,Tauxy_d_right
    real(kind=dp),dimension(:,:),allocatable::Tauxx_d_left,Tauxx_d_right,Tauyy_d_lower,Tauyy_d_upper
    ! Heat transfer
    real(kind=dp),dimension(:,:),allocatable::qxx_a_left,qxx_a_right,qyy_a_lower,qyy_a_upper
    real(kind=dp),dimension(:,:),allocatable::qxx_b_left,qxx_b_right,qyy_b_lower,qyy_b_upper
    real(kind=dp),dimension(:,:),allocatable::qxx_c_left,qxx_c_right,qyy_c_lower,qyy_c_upper
    real(kind=dp),dimension(:,:),allocatable::qxx_d_left,qxx_d_right,qyy_d_lower,qyy_d_upper
    ! Calculate the x_direction invicid fluxes (E_inv)
    real(kind=dp),dimension(:,:),allocatable::E1_inv_a_left,E2_inv_a_left,E3_inv_a_left,E4_inv_a_left
    real(kind=dp),dimension(:,:),allocatable::E1_inv_a_right,E2_inv_a_right,E3_inv_a_right,E4_inv_a_right
    real(kind=dp),dimension(:,:),allocatable::E1_inv_b_left,E2_inv_b_left,E3_inv_b_left,E4_inv_b_left
    real(kind=dp),dimension(:,:),allocatable::E1_inv_b_right,E2_inv_b_right,E3_inv_b_right,E4_inv_b_right
    real(kind=dp),dimension(:,:),allocatable::E1_inv_c_left,E2_inv_c_left,E3_inv_c_left,E4_inv_c_left
    real(kind=dp),dimension(:,:),allocatable::E1_inv_c_right,E2_inv_c_right,E3_inv_c_right,E4_inv_c_right
    real(kind=dp),dimension(:,:),allocatable::E1_inv_d_left,E2_inv_d_left,E3_inv_d_left,E4_inv_d_left
    real(kind=dp),dimension(:,:),allocatable::E1_inv_d_right,E2_inv_d_right,E3_inv_d_right,E4_inv_d_right
    ! Calculate the y_direction invicid fluxes (F_inv)
    real(kind=dp),dimension(:,:),allocatable::F1_inv_a_lower,F2_inv_a_lower,F3_inv_a_lower,F4_inv_a_lower
    real(kind=dp),dimension(:,:),allocatable::F1_inv_a_upper,F2_inv_a_upper,F3_inv_a_upper,F4_inv_a_upper
    real(kind=dp),dimension(:,:),allocatable::F1_inv_b_lower,F2_inv_b_lower,F3_inv_b_lower,F4_inv_b_lower
    real(kind=dp),dimension(:,:),allocatable::F1_inv_b_upper,F2_inv_b_upper,F3_inv_b_upper,F4_inv_b_upper
    real(kind=dp),dimension(:,:),allocatable::F1_inv_c_lower,F2_inv_c_lower,F3_inv_c_lower,F4_inv_c_lower
    real(kind=dp),dimension(:,:),allocatable::F1_inv_c_upper,F2_inv_c_upper,F3_inv_c_upper,F4_inv_c_upper
    real(kind=dp),dimension(:,:),allocatable::F1_inv_d_lower,F2_inv_d_lower,F3_inv_d_lower,F4_inv_d_lower
    real(kind=dp),dimension(:,:),allocatable::F1_inv_d_upper,F2_inv_d_upper,F3_inv_d_upper,F4_inv_d_upper
    ! Calculate the x_direction viscous fluxes (E_vis)
    real(kind=dp),dimension(:,:),allocatable::E1_vis_a_left,E2_vis_a_left,E3_vis_a_left,E4_vis_a_left
    real(kind=dp),dimension(:,:),allocatable::E1_vis_a_right,E2_vis_a_right,E3_vis_a_right,E4_vis_a_right
    real(kind=dp),dimension(:,:),allocatable::E1_vis_b_left,E2_vis_b_left,E3_vis_b_left,E4_vis_b_left
    real(kind=dp),dimension(:,:),allocatable::E1_vis_b_right,E2_vis_b_right,E3_vis_b_right,E4_vis_b_right
    real(kind=dp),dimension(:,:),allocatable::E1_vis_c_left,E2_vis_c_left,E3_vis_c_left,E4_vis_c_left
    real(kind=dp),dimension(:,:),allocatable::E1_vis_c_right,E2_vis_c_right,E3_vis_c_right,E4_vis_c_right
    real(kind=dp),dimension(:,:),allocatable::E1_vis_d_left,E2_vis_d_left,E3_vis_d_left,E4_vis_d_left
    real(kind=dp),dimension(:,:),allocatable::E1_vis_d_right,E2_vis_d_right,E3_vis_d_right,E4_vis_d_right
    ! Calculate the y_direction viscous fluxes (F_vis)
    real(kind=dp),dimension(:,:),allocatable::F1_vis_a_lower,F2_vis_a_lower,F3_vis_a_lower,F4_vis_a_lower
    real(kind=dp),dimension(:,:),allocatable::F1_vis_a_upper,F2_vis_a_upper,F3_vis_a_upper,F4_vis_a_upper
    real(kind=dp),dimension(:,:),allocatable::F1_vis_b_lower,F2_vis_b_lower,F3_vis_b_lower,F4_vis_b_lower
    real(kind=dp),dimension(:,:),allocatable::F1_vis_b_upper,F2_vis_b_upper,F3_vis_b_upper,F4_vis_b_upper
    real(kind=dp),dimension(:,:),allocatable::F1_vis_c_lower,F2_vis_c_lower,F3_vis_c_lower,F4_vis_c_lower
    real(kind=dp),dimension(:,:),allocatable::F1_vis_c_upper,F2_vis_c_upper,F3_vis_c_upper,F4_vis_c_upper
    real(kind=dp),dimension(:,:),allocatable::F1_vis_d_lower,F2_vis_d_lower,F3_vis_d_lower,F4_vis_d_lower
    real(kind=dp),dimension(:,:),allocatable::F1_vis_d_upper,F2_vis_d_upper,F3_vis_d_upper,F4_vis_d_upper
    ! The fluxes derivatives
    real(kind=dp),dimension(:,:),allocatable::dU1dt_a,dU2dt_a,dU3dt_a,dU4dt_a
    real(kind=dp),dimension(:,:),allocatable::dU1dt_b,dU2dt_b,dU3dt_b,dU4dt_b
    real(kind=dp),dimension(:,:),allocatable::dU1dt_c,dU2dt_c,dU3dt_c,dU4dt_c
    real(kind=dp),dimension(:,:),allocatable::dU1dt_d,dU2dt_d,dU3dt_d,dU4dt_d
    real(kind=dp),dimension(:,:),allocatable::dU1dt,dU2dt,dU3dt,dU4dt
    ! The Time-flux Vector_U_new
    real(kind=dp),dimension(:,:),allocatable::U1_new,U2_new,U3_new,U4_new
    ! Stopping criteria
    real(kind=dp),dimension(:,:),allocatable::err1,err2,err3,err4
    real(kind=dp),dimension(:,:),allocatable::Cf,Cf_res,Ch,Ch_res
    real(kind=dp) :: NormU1, NormU2, NormU3, NormU4


    !Variables for post processing purposes
    real(kind=dp),dimension(:,:),allocatable :: u_oldGlobal,v_oldGlobal,r_oldGlobal,T_oldGlobal
    real(kind=dp),dimension(:,:),allocatable :: dU1dtGlobal,dU2dtGlobal,dU3dtGlobal,dU4dtGlobal
    real(kind=dp),dimension(:,:),allocatable :: Vorticity,GradRho, NGradRho, Psi, MachRe, MachHo,MachVe
    real(kind=dp),dimension(:,:),allocatable :: Pres,MachNormal,MagGradEntropy,GradEntroI,GradEntroJ,Enthalpy_Stag
    real(kind=dp),dimension(:,:),allocatable :: threshold, MagPresGrad, Q_crit, P_total, T_total,DeltaS,DeltaH
    real(kind=dp),dimension(:,:), allocatable :: qflux_wall, Tau, SkinFriction, Ustar, yplus
    real(kind=dp),dimension(:,:), allocatable :: qflux_wallGlobal, TauGlobal, SkinFrictionGlobal, yplusGlobal
    real(kind=dp)::Max_Vorticity, Max_GradRho,S11,S12,S22,W11,W12,W22,MaxDeltaS, P_TotalFree
    real(kind=dp):: KineticEnergy_Initial, KineticEnergy, PrintFrecuency
    CHARACTER (len=10)  :: fileout
    CHARACTER (len=30)  :: name



    contains

    Function Pi()
        real(kind=dp)::Pi
        Pi = 4.0*atan(1.0)
    End Function

End Module variables
