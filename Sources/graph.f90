


subroutine graph
  use kind
  use data!, only: timestep, time, dt, NTN, NTE, rcoordinate, zcoordinate, usol, vsol, psol, globalNM, step, size_function_change, final_size


  implicit none

  integer(kind=ik):: i, step_indicator
  character(LEN=4):: fileN
  character(LEN=40):: file_mesh, temp


  if(s_mode.eq.1) then

     if(graph_mode.eq.0) then
        step_indicator = size_function_change
     else
        step_indicator = step
     end if

     if (step_indicator.eq.0) then
        open(unit = 10, file = trim(folder)//'mesh.dat', status = 'replace')
     else
        open(unit = 10, file = trim(folder)//'mesh.dat', status = 'old', access = 'append')
     end if
     call graph_mesh(step_indicator)
write(*,*) 'write mesh'
    
     if( final_size.eq.1) then

        if(read_coordinate_value.eq.0) then
           open(unit = 20, file = 'Sources/elliptic_mesh.dat', status = 'replace')
           write(20,'(4i4,2es13.6)') NEL, NEM, NEM_alge, NEV, outer, substrate
           do i = 1, NTN
              write(20,'(2es15.7)') rcoordinate(i), zcoordinate(i)
           end do
           close(20)
           print *, 'write elliptic mesh result'
        end if

        if( no_vapor.eq.1 ) then
           open(unit = 11, file = trim(folder)//'dynamics.dat', status = 'replace')
           zone = 1
           open(unit = 12, file = trim(folder)//'temperature.dat', status = 'replace')
           if(no_vapor.eq.0) &
                open(unit = 13, file = trim(folder)//'vapor_concentration.dat', status = 'replace' )
           call graph_dynamics(0)
write(*,*) 'dynamics'
        end if
     end if

  else   !s_mode = 0

     ! if(initial_vapor_solving.eq.1) then
     !    if(step.eq.0) then
     !       open(unit = 11, file = trim(folder)//'vapor_initial.dat', status = 'replace')
     !    else
     !       open(unit = 11, file = trim(folder)//'vapor_initial.dat', status = 'old', access = 'append')
     !    end if
     !    call graph_dynamics(step)
     ! else

        if(graph_mode.eq.0) then
           step_indicator = timestep
        else
           step_indicator = step
        end if

        if(diverge.eq.1) then
           if(step_indicator.eq.1) then
              open(unit = 11, file = trim(folder)//'divergence.dat', status = 'replace')
           else if(graph_mode.eq.1) then
              open(unit = 11, file = trim(folder)//'divergence.dat', status = 'old', access = 'append')
           end if
        else  !not diverge
           if(no_vapor.eq.0 .and. step_indicator.eq.0) then
              open(unit = 11, file = trim(folder)//'dynamics.dat', status = 'replace')
              zone = 1
              open(unit = 12, file = trim(folder)//'temperature.dat', status = 'replace')
              if(no_vapor.eq.0) &
                   open(unit = 13, file = trim(folder)//'vapor_concentration.dat', status = 'replace')
           else
              open(unit = 11, file = trim(folder)//'dynamics.dat', status = 'old', access = 'append')
              zone = zone + 1
              open(unit = 12, file = trim(folder)//'temperature.dat', status = 'old', access = 'append')
              if(no_vapor.eq.0) &
                   open(unit = 13, file = trim(folder)//'vapor_concentration.dat', status = 'old', access = 'append')
           end if
        end if

        ! if( step_indicator.eq.1 ) then
        !    if(diverge.eq.1 ) then
        !       open(unit = 11, file = trim(folder)//'divergence.dat', status = 'replace')
        !    else
        !       open(unit = 11, file = trim(folder)//'dynamics.dat', status = 'replace')
        !       zone = 1
        !       open(unit = 12, file = trim(folder)//'temperature.dat', status = 'replace')
        !       if(no_vapor.eq.0) &
        !            open(unit = 13, file = trim(folder)//'vapor_concentration.dat', status = 'replace')
        !    end if
        ! else
        !    if(diverge.eq.1 .and. graph_mode.eq.1) then
        !       open(unit = 11, file = trim(folder)//'divergence.dat', status = 'old', access = 'append')
        !    else
        !       open(unit = 11, file = trim(folder)//'dynamics.dat', status = 'old', access = 'append')
        !       zone = zone + 1
        !       open(unit = 12, file = trim(folder)//'temperature.dat', status = 'old', access = 'append')
        !       if(no_vapor.eq.0) &
        !            open(unit = 13, file = trim(folder)//'vapor_concentration.dat', status = 'old', access = 'append')
        !    end if
        ! end if
write(*,*) 'write dnamics'
        call graph_dynamics(step_indicator)

     ! end if

  end if


  return



contains

  subroutine graph_mesh(step_indicator)
    use data, only: graph_mode

    integer(kind=ik), intent(in)::step_indicator 
    integer(kind=ik):: Nregion, Nnode, Nele


    do Nregion = 1, 3
       if(substrate.eq.0.0_rk .and. Nregion.eq.3) cycle 
       if(no_vapor.eq.1 .and. Nregion.eq.2) cycle

       select case(Nregion)
       case(1)
          Nnode = NND
          Nele = 4*ED
       case(2)
          Nnode = NNV
          Nele = 4*EV
       case(3)
          Nnode = NNS
          Nele = 4*ES
       end select


       write(10, '(A)') 'variables = "r", "z"'
       write(10,200) 'Zone T = "step:', step_indicator, ', RGN: ', Nregion, &
            '", STRANDID = 1, SOLUTIONTIME =', step_indicator, &
            ', Datapacking = Point, Zonetype = FEQuadrilateral, N =', Nnode, ', E =', Nele, &
            ', DT = (double,double)'
200    format(A,i8,A,i8,A,i8,A,i8,A,i8,A)

       do i = 1, NTN, 1
          if (Nregion.eq.1 .and. (VN(i).eq.1 .or. VN(i).eq.5) ) cycle
          if (Nregion.eq.2 .and. (VN(i).eq.0 .or. VN(i).eq.5) ) cycle
          if (Nregion.eq.3 .and. (.not.( VN(i).eq.5 .or. BCflagN(i,2).ne.0 )) ) cycle          
          write(10,'(2es15.7)') rcoordinate(i), zcoordinate(i)
       end do

       !form the quadrilateral with nodes
       nodes: do i = 1, NTE, 1
          if (Nregion.eq.1 .and. VE(i).ne.0 ) cycle
          if (Nregion.eq.2 .and. VE(i).ne.1 ) cycle
          if (Nregion.eq.3 .and. VE(i).ne.5 ) cycle
          ! write(10,'(4i8)') globalNM(i,1), globalNM(i,3), globalNM(i,9), globalNM(i,7)
          write(10,'(4i8)') NOPDV( globalNM(i,1), Nregion ), NOPDV( globalNM(i,2), Nregion ), &
               NOPDV( globalNM(i,5), Nregion ), NOPDV( globalNM(i,4), Nregion )
          write(10,'(4i8)') NOPDV( globalNM(i,2), Nregion ), NOPDV( globalNM(i,3), Nregion ), &
               NOPDV( globalNM(i,6), Nregion ), NOPDV( globalNM(i,5), Nregion )
          write(10,'(4i8)') NOPDV( globalNM(i,4), Nregion ), NOPDV( globalNM(i,5), Nregion ), &
               NOPDV( globalNM(i,8), Nregion ), NOPDV( globalNM(i,7), Nregion )
          write(10,'(4i8)') NOPDV( globalNM(i,5), Nregion ), NOPDV( globalNM(i,6), Nregion ), &
               NOPDV( globalNM(i,9), Nregion ), NOPDV( globalNM(i,8), Nregion )
       end do nodes

    end do

    close(10)


    return
  end subroutine graph_mesh




  subroutine graph_dynamics(step_indicator)
    use data, only: graph_mode

    integer(kind=ik), intent(in):: step_indicator
    integer(kind=ik):: Nregion, Nnode, Nele
    real(kind=rk):: time_indicator

    if(graph_mode.eq.0 .and. ( initial_vapor_solved.eq.1 .or. no_vapor.eq.1 ) ) then
       ! step_indicator = timestep
       time_indicator = time
    else
       ! step_indicator = step
       time_indicator = real(step,rk)
    end if
    
    if(diverge.eq.0) then
       if(Nregion.eq.1) fileN = '11'
       if(Nregion.eq.1 .or. Nregion.eq.3) fileN = '12'
       if(Nregion.eq.2) fileN = '13'
    else
       fileN = '11'
    end if
          
    do Nregion = 1, 3
       if(substrate.eq.0.0_rk .and. Nregion.eq.3) cycle 
       if(no_vapor.eq.1 .and. Nregion.eq.2) cycle 

       select case(Nregion)
       case(1)
          Nnode = NND
          Nele = ED
       case(2)
          Nnode = NNV
          Nele = EV
       case(3)
          Nnode = NNS
          Nele = ES
       end select

       !write headlines
       if(diverge.eq.0) then
          if(Nregion.eq.1) then
             write(11, '(A)') 'variables = "r", "z", "u", "v", "T", "p"'
             ! write(11, '(A,es13.6,A)') 'DATASETAUXDATA Umax = "', umax, '"'
             ! write(11, '(A,es13.6,A)') 'DATASETAUXDATA Vmax = "', vmax, '"'

             write(11,202) 'Zone T = "step:', step_indicator, &!', RGN: ', Nregion, &
                  '", STRANDID = 1, SOLUTIONTIME =', time_indicator, &
                  ', Datapacking = Point, Zonetype = FEQuadrilateral, N =', Nnode, ', E =', Nele, &
                  ', DT = (double,double,double,double,double,double), &
                  AUXDATA angle = "',angle_c_degree, '", AUXDATA UMAX = "', umax, '", &
                  AUXDATA VMAX = "', vmax, '", AUXDATA ZTOP = "', ztop, '"'
202    format(A,i8,A,es14.7,A,i8,A,i8,A, f7.3, A,es13.6,A,es13.6,A,ES13.6,A)

          end if
          if(Nregion.eq.1 .or. Nregion.eq.3) then
             write(12, '(A)') 'variables = "r", "z", "T"'
             write(12,203) 'Zone T = "step:', step_indicator, ', RGN: ', Nregion, &
                  '", STRANDID = 1, SOLUTIONTIME =', time_indicator, &
                  ', Datapacking = Point, Zonetype = FEQuadrilateral, N =', Nnode, ', E =', Nele, &
                  ', DT = (double,double,double), &
                  AUXDATA angle = "',angle_c_degree, '"'
203    format(A,i8,A,i8,A,es14.7,A,i8,A,i8,A, f7.3, A)
          end if
          if(Nregion.eq.2) then
             ! pause
             ! write(*,*) 'headline'
             write(13, '(A)') 'variables = "r", "z", "c"'
             write(13,201) 'Zone T = "step:', step_indicator, &!', RGN: ', Nregion, &
                  '", STRANDID = 1, SOLUTIONTIME =', time_indicator, &
                  ', Datapacking = Point, Zonetype = FEQuadrilateral, N =', Nnode, ', E =', Nele, &
                  ', DT = (double,double,double)'
          end if
       else
          write(11, '(A)') 'variables = "r", "z", "u", "v", "T", "p", "c"'
          write(11,201) 'Zone T = "step:', step_indicator, &!', RGN: ', Nregion, &
               '", STRANDID = 1, SOLUTIONTIME =', time_indicator, &
               ', Datapacking = Point, Zonetype = FEQuadrilateral, N =', Nnode, ', E =', Nele, &
               ', DT = (double,double,double,double,double,double)'
       end if
       ! write(11,201) 'Zone T = "step:', step_indicator, ', RGN: ', Nregion, &
       !      '", STRANDID = 1, SOLUTIONTIME =', time_indicator, &
       !      ', Datapacking = Point, Zonetype = FEQuadrilateral, N =', Nnode, ', E =', Nele, &
       !      ', DT = (double,double,double,double,double,double)'
201    format(A,i8,A,es14.7,A,i8,A,i8,A)

       do i = 1, NTN, 1
          if (Nregion.eq.1 .and. (VN(i).eq.1 .or. VN(i).eq.5) ) cycle
          if (Nregion.eq.2 .and. (VN(i).eq.0 .or. VN(i).eq.5) ) cycle
          if (Nregion.eq.3 .and. (.not.( VN(i).eq.5 .or. BCflagN(i,2).ne.0 )) ) cycle

          if(diverge.eq.0) then
             if(Nregion.eq.1) then
                write(11,'(7es15.7)') rcoordinate(i), zcoordinate(i), usol(i), vsol(i), Tsol(i), psol(i)
             end if
             if(Nregion.eq.1 .or. Nregion.eq.3) then
                write(12,'(7es15.7)') rcoordinate(i), zcoordinate(i), Tsol(i)
             end if
             if(Nregion.eq.2) then
                write(13,'(7es15.7)') rcoordinate(i), zcoordinate(i), csol(i)
             end if
          else
             write(11,'(7es15.7)') rcoordinate(i), zcoordinate(i), usol(i), vsol(i), Tsol(i), psol(i), csol(i)
          end if
          
       end do

       !form the quadrilateral with nodes
       do i = 1, NTE, 1
          if (Nregion.eq.1 .and. VE(i).ne.0 ) cycle
          if (Nregion.eq.2 .and. VE(i).ne.1 ) cycle
          if (Nregion.eq.3 .and. VE(i).ne.5 ) cycle

          if(diverge.eq.0) then
             if(Nregion.eq.1) then
                write(11,'(4i8)') NOPDV(globalNM(i,1),Nregion), NOPDV(globalNM(i,3),Nregion), &
                     NOPDV(globalNM(i,9),Nregion), NOPDV(globalNM(i,7),Nregion)
             end if
             if(Nregion.eq.1 .or. Nregion.eq.3) then
                write(12,'(4i8)') NOPDV(globalNM(i,1),Nregion), NOPDV(globalNM(i,3),Nregion), &
                     NOPDV(globalNM(i,9),Nregion), NOPDV(globalNM(i,7),Nregion)
             end if
             if(Nregion.eq.2) then
                write(13,'(4i8)') NOPDV(globalNM(i,1),Nregion), NOPDV(globalNM(i,3),Nregion), &
                     NOPDV(globalNM(i,9),Nregion), NOPDV(globalNM(i,7),Nregion)
             end if
          else
             write(11,'(4i8)') NOPDV(globalNM(i,1),Nregion), NOPDV(globalNM(i,3),Nregion), &
                  NOPDV(globalNM(i,9),Nregion), NOPDV(globalNM(i,7),Nregion)
          end if
          ! write(fileN,'(4i8)') NOPDV(globalNM(i,1),Nregion), NOPDV(globalNM(i,3),Nregion), &
          !      NOPDV(globalNM(i,9),Nregion), NOPDV(globalNM(i,7),Nregion)
       end do

    end do

    write(11,'(A,f7.3,A,i7,A,i7)') 'Text X=40, Y=90, F=Times, T= "contact angle =', angle_c_degree, '", ZN= ', zone

    close(11)
    close(12)
    if(no_vapor.eq.0) close(13)

    return
  end subroutine graph_dynamics




end subroutine graph
