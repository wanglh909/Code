subroutine newton_raphson
  use kind
  use data
  use NOP_mod, only: damfac
  use front_mod, only: multifront, load_dum, debug_NAN

  implicit none
  real(kind=rk):: cal_time
  real(kind=rk), parameter:: TOL = 5.0e-5_rk
  integer(kind=ik):: i


  do       !newtons iteration loop

     step = step + 1

     !define soldot
     if (timestep.le.FTS) then
        soldot = ( sol - solp )/dt
     else
        soldot = 2.0_rk*( sol - solp )/dt - soldotp
     end if
    
     call multifront(error1, cal_time) !using subroutine assemble
     !print *, "Solver time:", cal_time
     
! if(initial_vapor_solved.eq.1) pause

     if(debug_NAN) pause 
     
     !for debug, put 'sj's together as Jac, check Jac
     if( check_0_in_Jac.eq.1 .and. s_mode.eq.0) then!.and. initial_vapor_solved.eq.1  .and. pack_start.eq.1 
        call jac_check_0
        Jac(:,:) = 0.0_rk
     end if
     ! if(check_0_in_Jac.eq.1) 
     

     call L2_error(cal_time)
     !check Jacobian whole if error NaN shows
     if(error1.ne.error1 .or. error2.ne.error2) then
        print *, 'error NaN'
        check_0_in_Jac = 1
        if(diverge.eq.0) then
           diverge = 1
           allocate( Jac(NVar,NVar), Res(NVar) )
           Jac = 0.0_rk
           Res = 0.0_rk
           exit
        else !diverge.eq.1
           stop
        end if
     end if

     !update sol
     if( s_mode.eq.1 ) then
        sol = sol + damfac(error1) * dsol
     else
        sol = sol + dsol
     end if
     call split_sol        !made u,v,r,z at boundary 0       !?do before L2 error
     !need to do for every iteration to prepare for (r,z,u,v,p)local
     if( graph_mode.eq.1 .or. initial_vapor_solving.eq.1) call graph


     !pause
     
     !bug exsistance, (stop program) or ( diverge=1 and redo this timestep to record )
     !if ( ( step.gt.20 .and. (timestep.ne.0 .and. timestep.ne.1) ) .or. error2.gt.1.0e8_rk  )  then 
     if ( ( step.gt.10 .and. (timestep.ne.0) ) .or. error2.gt.1.0e8_rk  )  then
        print *, 'did not converge, diverge number', diverge, 's_mode', s_mode
        if(diverge.eq.1 .or. s_mode.eq.1) then
           print *, 'data wrote in', folder
           ! call system('preplot '//trim(folder)//'dynamics.dat')
           ! call system('preplot '//trim(folder)//'mesh.dat')
           ! call system('preplot '//trim(folder)//'divergence.dat')
           ! call system('preplot '//trim(folder)//'temperature.dat')
           ! call system('preplot '//trim(folder)//'vapor_concentration.dat')
           stop
        else
           sol = solp
           call split_sol
           timestep = timestep - 1
           timestep_stable = timestep_stable - 1  !?useful?
           time = time - dt
           diverge = 1
! print *, 'diverge exit'
           exit
        end if 
     end if

     !convergence reached
     if ( (error1.lt.TOL.and.error2.lt.TOL) .or. &
          ( s_mode.eq.1 .and. final_size.eq.0 .and. (error1.lt.2.0_rk .or. error2.lt.2.0_rk) ) ) then
        print *, 'this timestep converged'
        if(diverge.eq.1) diverge=0
        !update soldot
        if (timestep.le.FTS) then
           soldot = ( sol - solp )/dt
        else
           soldot = 2.0_rk*( sol - solp )/dt - soldotp
        end if
! print *, final_size, 'converge exit'
        exit
     end if
     ! if(final_size.eq.0 .and. ( error1.lt.1.0e-2_rk.and.error2.lt.1.0e-2_rk )) exit
  end do

  
  return
end subroutine newton_raphson
