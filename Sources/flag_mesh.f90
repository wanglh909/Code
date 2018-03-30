subroutine flag_mesh
  use kind
  use data, only: step, final_size, size_function_change, s_mode, &
       initial_vapor_solved, initial_vapor_solving, MDF, MDFd, diverge, angle_c, pi, no_vapor
  use front_mod, only: determine_offsets
  use NOP_mod, only: NOPP_define
  implicit none

  if(step.eq.0) then

     !mesh size change
     if(final_size.eq.0) then                
        call size_function!(size_function_change)  
        write(*,*) 'size_function_change=', size_function_change, ', contact angle =', angle_c/pi*180.0_rk
        size_function_change = size_function_change + 1
     end if
     !adjust vapor concentration
     if(no_vapor.eq.0) then
        if(s_mode.eq.0 .and. initial_vapor_solved.eq.0) then
           write(*,*) 'solving for initial vapor concentration'
           initial_vapor_solving = 1!
        end if
     else ! no_vapor=1
        if(s_mode.eq.0) initial_vapor_solved = 1
     end if

     if(diverge.eq.1) write(*,*) 'repeat diverged timestep'


  else

     !mesh size change finished
     if (s_mode.eq.1 .and. final_size.eq.1) then
        MDF(:) = MDFd(:)
        s_mode = 0!
        !call NOPP_define
        call determine_offsets()
     end if
     !adjust vapor concentration finished
     if(no_vapor.eq.0) then
        if(initial_vapor_solving.eq.1 .and. initial_vapor_solved.eq.0) then
           initial_vapor_solved = 1!
           initial_vapor_solving = 0!
        end if
     end if


  end if


  return
end subroutine flag_mesh
