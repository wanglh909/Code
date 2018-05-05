program main
  use kind
  use omp_lib
  use data
  use NOP_mod
  use basis_f
  use front_mod, only: debug_NAN, load_balance_flag, SINGLE, DUAL, DOMAINS, &
       THREADS_FRONT, multifront, determine_offsets, SWAP_LOCAL_IJ!, seed
  implicit none

  real(kind=rk):: t1, t_program

  t1 = REAL(omp_get_wtime(),rk)

  !initial droplet shape
  angle_c_degree = 50.0_rk
  angle_c = angle_c_degree /180.0_rk*pi
  !ths = 2
  no_vapor = 1       !no_vapor = 1: do not solve for vapor phase, impose function flux
  uniflux = 0  !determine if the imposed flux is uniform. Notice: If flux is uniform, still use divergent heat flux.
  call data_folder  !determine xe

  !set mesh parameters
  !2 6 4 2 2 simple mesh
  NEL = 8!10    !input   
  NEM = 200  !input     !decide in data_folder
  NEV = 50 !200!     1000!  !input
  NES = 5   !10!  30  !input
  NEM_alge = 90!450!90   !  135!NEM/3*2    !decide in data_folder
  
!----------------------------------front settings-------------------------------------
  SWAP_LOCAL_IJ = .FALSE.             !This determines if your local is transposed
  debug_NAN = .FALSE.                !True for NaN debugging (HUGE PERFORMANCE PENALTY)
  load_balance_flag = .TRUE.        !True for load balancing, good to have on generally
  ths = 4
  THREADS_FRONT = ths                       !THREADS for front, must be <= ths
  SOLVER_MODE = DOMAINS
  !Call init_front( .... ) with argument of SINGLE, DUAL, or DOMAINS (Those are module vars) (ex: call init_front(DOMAINS) )
  !call before first use, when changing solvers, or when switch which vars are solved
  !Only call multifront(L2res,time) now, it will use the initialized solver

  !Basically always use DOMAINS, load_balance_flag = .TRUE.
  !Use single or dual for debugging, same with debug_NAN
!-----------------------------------done settings------------------------------------
  
  !set terms option
  NStrans = 1
  Inert = 1
  Capil = 1
  Viscous = 1
  GravI = 1
  if(Inert.eq.0 .and. Viscous.eq.0) then
     write(*,*) 'Inertia and Viscous cannot be 0 at the same time'
     stop
  end if
  Ttime = 1
  Tconv = 1
  Tdiff = 1
  TtimeS = 1
  TdiffS = 1
  !set flow parameters
  call parameter_values

  alge_corner = 1


  graph_step = 10  !graph every 'graph_step' steps
  dt = 1.0e-5_rk!0.01_rk   !dt in first 5 steps
  FTS = 5 !fixed timesteps

  !debug flag
  simple_mesh = 1     ! !1: use simple mesh for quicker calculation
  graph_mode = 0    !1: graph each step; 0: graph each timestep
  check_0_in_Jac = 0   !1: put 'sj's together as Jac, check Jac
  if( simple_mesh.eq.1 ) then
     angle_c_degree = 10.0_rk!90.0_rk
     angle_c = angle_c_degree /180.0_rk*pi
     NEL = 4
     NEM = 5
     NEV = 3!6
     NES = 2
     NEM_alge = NEM/3*2
     !--------solver setting-------------
     ths = 1
     THREADS_FRONT = ths                      !THREADS for front, must be <= ths
     SOLVER_MODE = single
     !-------end solver------------------
     dt = 1.0e-5_rk!0.01_rk   !dt in first 5 steps
  end if
  if( no_vapor.eq.0 ) NEV = 3 !?
  if( ( outer.eq.1.0_rk .or. substrate.eq.0.0_rk ) .and. no_vapor.eq.1 ) NEV = 0
  
  if(substrate.eq.0.0_rk) NES = 0

!-----------------------------------calculation start--------------------------------------

  call vapor_mesh_size
  call NOP   ! NNR1, NNR1p, NNX, NTN, NTE, globalNM, rowNM(eleN), columnNM, regN(eleN), WFLAG(eleN)
  call variableN    !NOPP(NTN), MDF(NTN), NVar, rNOP(NTN,4,2), iBW
  call basis_function
  call initialization   !allocate, make everything 0, call multifront preparation(init_front, assembler, associater, excluder, custom_order, var_finder, seed)
  call initial_condition
! call split_sol

!  seed = NTE  !??see if can be put into initialization
  
  !switch to mesh-only mode
  MDF(:) = 2
  s_mode = 1
  !call NOPP_define
  call determine_offsets()  !in multifront module, call whenever switching between mesh and full dynamics solving
  call graph  !graph starting mesh

  write(*,*) 'enter loops'
  do    !loop for time step
     
     write(*,*)'------------------------next timestep-----------------------------'

     !start solving for solutions for each timestep
     step = 0
     call flag_mesh   !use size_function (viz. fsize)
     !flag adjustment for mesh establish(mesh_size_change & initial_vapor_solv), use size_function

     if(graph_mode.eq.1 .or. initial_vapor_solving.eq.1)  call graph
     if(initial_vapor_solved.eq.1) call prediction  !timestep+1, time+dt

     call newton_raphson  !calculating part, use multifront, jac_check_0, l2_error, split_sol, graph

     if(graph_mode.eq.1 .and. diverge.eq.0) then
        if(s_mode.eq.0) then
           write(*,*) 'pause for every timestep if solve for dynamics and graph each step'
           pause
        end if
     else if( timestep.le.20 .or. mod(timestep,graph_step).eq.0 .or. diverge.eq.1 ) then  
        call graph          !graph every 'graph_step' timesteps or right before divergence
        if(diverge.eq.1)  graph_mode = 1
     end if

     call flag_mesh   !flag for mesh establish(mesh_size_change & initial_vapor_solv)
      
     !calculate contact angle & flux
     if((initial_vapor_solving.eq.1 .or. initial_vapor_solved.eq.1).and. diverge.eq.0)  call variable_cal

     t_program = REAL(omp_get_wtime(),rk) - t1
     open(unit = 20, file = trim(folder)//'cal_time.dat', status = 'old', access = 'append')
     write(20, '(A)') ' '
     write(20,'(A,es13.6,A,es13.6,A,es13.6,A,es13.6,A)') 'total:', t_program,'s, ', &
          t_program/60.0,'min, ', t_program/60.0/60.0,'hr, ', t_program/60.0/60.0/24.0,'d'
     close(20)

     !***********************************conditions to stop time loop*********************************
     !if(timestep.eq.20) stop
     if(angle_c.le.0.0_rk) stop

     !****************************************************************************************

  end do

end program main
