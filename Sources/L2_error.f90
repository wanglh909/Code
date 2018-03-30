subroutine L2_error(cal_time)
  use kind
  use data
  use front_mod, only: load_dum

  implicit none
  real(kind=rk), intent(in):: cal_time
  integer(kind=ik):: i, imaxr, imaxz, imaxu, imaxv, imaxT, imaxp, imaxc, write_node
  integer(kind=ik):: imaxRsi, imaxReta, imaxRu, imaxRv, imaxRt, imaxRp, imaxRc
  real(kind=rk):: dr(NTN), dz(NTN), du(NTN), dv(NTN), dTemp(NTN), dp(NTN), dc(NTN), &   
       drmax, dzmax, dumax, dvmax, dTmax, dpmax, dcmax      !dTemp dstinct from dt in data
  real(kind=rk):: Rsi(NTN), Reta(NTN), Ru(NTN), Rv(NTN), Rt(NTN), Rp(NTN), Rc(NTN), &
       Rsimax, Retamax, Rumax, Rvmax, Rtmax, Rpmax, Rcmax
  
  write_node = 1  !1: write which node gives the greatest error
  
  !find maximum value of each variable, used to normalize error
  rmax = 1.0_rk!0.0_rk
  zmax = 1.0_rk!0.0_rk
  umax = 0.0_rk
  vmax = 0.0_rk
  Tmax = 0.0_rk
  pmax = 0.0_rk
  cmax = 0.0_rk
  do i = 1, NTN
     if(s_mode.eq.0) then
        if( abs(rcoordinate(i)) .gt. rmax ) rmax = abs(rcoordinate(i))
        if( abs(zcoordinate(i)) .gt. zmax ) zmax = abs(zcoordinate(i))
     end if
     if( abs(usol(i)) .gt. umax ) umax = abs(usol(i))
     if( abs(vsol(i)) .gt. vmax ) vmax = abs(vsol(i))
     if( abs(Tsol(i)) .gt. Tmax ) Tmax = abs(Tsol(i))
     if( abs(psol(i)) .gt. pmax ) pmax = abs(psol(i))
     if( no_vapor.eq.0 .and. abs(csol(i)) .gt. cmax ) cmax = abs(csol(i))
  end do

  if (timestep.eq.0 .and. step.eq.1 .and. ( size_function_change.eq.1 .or. read_coordinate_value.eq.1 ) ) then
     open(unit = 11, file = trim(folder)//'max_value.dat', status = 'replace')
     write(11, '(A)') 'variables = "rmax", "zmax", "umax", "vmax", "Tmax", "pmax", "cmax" '
  else
     open(unit = 11, file = trim(folder)//'max_value.dat', status = 'old', access = 'append')
  end if
  if(step.eq.1) then
     write(11, '(A)') ' '
     write(11, '(A, 2es15.7, i8, A)') '----------------', time, dt, timestep, '---------------------'
     write(11, '(A)') ' '
  end if
  write(11,'(7es14.6)') rmax, zmax, umax, vmax, Tmax, pmax, cmax
  close(11)
  ! write(*,*) rmax, zmax, umax, vmax, Tmax, pmax, cmax  


  
  !calculate solution error error2, L2 norm
  if(Tmax.eq.0.0_rk) Tmax = 1.0_rk
  if(umax.eq.0.0_rk) umax = 1.0_rk
  if(vmax.eq.0.0_rk) vmax = 1.0_rk
  if(pmax.eq.0.0_rk) pmax = 1.0_rk
  if(cmax.eq.0.0_rk) cmax = 1.0_rk
  
  error2 = 0.0_rk
  do i = 1, NTN   !?change split_sol accordingly
     error2 = error2 + ( dsol(NOPP(i)+Nr)/rmax )**2
     if(NTs.eq.2 .or. VN(i).ne.5)  error2 = error2 + ( dsol(NOPP(i)+Nz)/zmax )**2
     if(s_mode.eq.0) then
        if(initial_vapor_solved.eq.1) then
           if( VN(i).eq.0 .or. VN(i).eq.2 .or. BCflagN(i,2).eq.1) then
              error2 = error2 + ( dsol(NOPP(i)+NT)/Tmax )**2
              if( VN(i).eq.0 .or. VN(i).eq.2 ) then
                 error2 = error2 + ( dsol(NOPP(i)+Nu)/umax )**2
                 error2 = error2 + ( dsol(NOPP(i)+Nv)/vmax )**2
                 if(PN(i).eq.1) error2 = error2 + ( dsol(NOPP(i)+Np)/pmax )**2
              end if
           end if
           if(VN(i).eq.5) error2 = error2 + ( dsol(NOPP(i)+NTs)/Tmax )**2
        end if
        if(no_vapor.eq.0 .and. ( VN(i).eq.1 .or. VN(i).eq.2 ) ) &
             error2 = error2 + ( dsol(NOPP(i)+MDF(i)-1)/cmax )**2
     end if
  end do
  ! do i = 1, NVar, 1      
  !    error2 = error2 + ( dsol(i)/sol(i) )**2
  ! end do
  error2 = sqrt(error2)             !absolute value of dsol

  !write Res & Sol error in file error.dat
  write(*,*) 'Res:', error1, ' Sol:', error2, ' step:', step
  if (timestep.eq.0 .and. step.eq.1 .and. ( size_function_change.eq.1 .or. read_coordinate_value.eq.1 ) ) then
     open(unit = 11, file = trim(folder)//'error.dat', status = 'replace')
     write(11, '(A)') 'variables = "ResError", "SolError", "step" '
  else
     open(unit = 11, file = trim(folder)//'error.dat', status = 'old', access = 'append')
  end if
  if(step.eq.1) then
     write(11, '(A)') ' '
     write(11, '(A, 2es15.7, i8, A)') '----------------', time, dt, timestep, '---------------------'
     write(11, '(A)') ' '
  end if
  write(11,200) 'Res. Error =', error1, '    Sol. Error =', error2, '    step =', step
  close(11)
200 format (A,es11.5,A,es11.5,A,i4)

  !calculate Sol error seperately, Linf norm
  dr(:) = 0.0_rk         
  dz(:) = 0.0_rk
  du(:) = 0.0_rk
  dv(:) = 0.0_rk
  dTemp(:) = 0.0_rk
  dp(:) = 0.0_rk
  dc(:) = 0.0_rk
  do i = 1, NTN, 1
     dr(i) = abs( dsol( NOPP(i) + Nr ) / rmax )
     if( NTs.eq.2 .or. VN(i).ne.5) dz(i) = abs( dsol( NOPP(i) + Nz ) / zmax )
     if(s_mode.eq.0) then
        if(VN(i).eq.0 .or. VN(i).eq.2) then
           du(i) = abs( dsol( NOPP(i) + Nu ) / umax )
           dv(i) = abs( dsol( NOPP(i) + Nv ) / vmax )
           dTemp(i) = abs( dsol( NOPP(i) + NT ) / Tmax )
           if ( PN(i).eq.1 ) dp(i) = abs( dsol( NOPP(i) + Np ) / pmax )
        end if
        if( (VN(i).eq.1 .or. VN(i).eq.2) ) then
           if(no_vapor.eq.0) dc(i) = abs( dsol( NOPP(i) + MDF(i)-1 ) / cmax )
           if( BCflagN(i,2).eq.1 )  dTemp(i) = abs( dsol( NOPP(i) + NT ) / Tmax )
        end if
        if (VN(i).eq.5) then
           dTemp(i) = abs( dsol( NOPP(i) + NTs ) / Tmax )
        end if
     end if
  end do

  !?calculate error2 here, when du, dv... have been splited

  drmax = 0.0_rk
  dzmax = 0.0_rk
  dumax = 0.0_rk
  dvmax = 0.0_rk
  dTmax = 0.0_rk
  dpmax = 0.0_rk
  dcmax = 0.0_rk
  imaxr = 0
  imaxz = 0
  imaxu = 0
  imaxv = 0
  imaxT = 0
  imaxp = 0
  imaxc = 0

  do i = 1, NTN
     if( dr(i) .gt. drmax ) then
        drmax = dr(i)
        imaxr = i
     end if
     if( dz(i) .gt. dzmax ) then
        dzmax = dz(i)
        imaxz = i
     end if
     if( du(i) .gt. dumax ) then
        ! if(BCflagN(i,3).eq.1) cycle   !???
        dumax = du(i)
        imaxu = i
     end if
     if( dv(i) .gt. dvmax ) then
        ! if(BCflagN(i,3).eq.1) cycle
        dvmax = dv(i)
        imaxv = i
     end if
     if( dTemp(i) .gt. dTmax ) then
        dTmax = dTemp(i)
        imaxT = i
     end if
     if( dp(i) .gt. dpmax ) then
        dpmax = dp(i)
        imaxp = i
     end if
     if( no_vapor.eq.0 .and. dc(i) .gt. dcmax ) then
        dcmax = dc(i)
        imaxc = i
     end if

  end do

  if (timestep.eq.0 .and. step.eq.1 .and. ( size_function_change.eq.1 .or. read_coordinate_value.eq.1 ) ) then
     open(unit = 10, file = trim(folder)//'error_Sol.dat', status = 'replace')
     write(10, '(A)') 'variables = "dr", "dz", "du", "dv", "dT", "dp", "dc" '
  else
     open(unit = 10, file = trim(folder)//'error_Sol.dat', status = 'old', access = 'append')
  end if
  if(step.eq.1) then
     write(10, '(A)') ' '
     write(10, '(A, 2es15.7, i8, A)') '----------------', time, dt, timestep, '---------------------'
     write(10, '(A)') ' '
  end if
  write(10,'(7es12.5)') drmax, dzmax, dumax, dvmax, dTmax, dpmax, dcmax
  if(write_node.eq.1) write(10,'(7i12)') imaxr, imaxz, imaxu, imaxv, imaxT, imaxp, imaxc
  write(10,'(es14.7)') error2
  write(10,'(A)') ' '
  close(10)


  !calculate Res error seperately, Linf norm
  Rsi(:) = 0.0_rk       
  Reta(:) = 0.0_rk
  Ru(:) = 0.0_rk
  Rv(:) = 0.0_rk
  Rt(:) = 0.0_rk
  Rp(:) = 0.0_rk
  Rc(:) = 0.0_rk
  do i = 1, NTN, 1
     Rsi(i) = abs( load_dum( NOPP(i) + Nr ) )
     Reta(i) = abs( load_dum( NOPP(i) + Nz ) )
     if(VN(i).eq.0 .or. VN(i).eq.2) then
        Ru(i) = abs( load_dum( NOPP(i) + Nu ) )
        Rv(i) = abs( load_dum( NOPP(i) + Nv ) )
        Rt(i) = abs( load_dum( NOPP(i) + NT ) )
        if ( MDF(i) .gt. Np )     Rp(i) = abs( load_dum( NOPP(i) + Np ) )
     end if
     if( (VN(i).eq.1 .or. VN(i).eq.2) .and. s_mode.eq.0) then
        if(no_vapor.eq.0) Rc(i) = abs( load_dum( NOPP(i) + MDF(i)-1 ) )
        if( BCflagN(i,2).eq.1 )  Rt(i) = abs( load_dum( NOPP(i) + NT ) )
     end if
     if (VN(i).eq.5 .and. s_mode.eq.0) then
        Rt(i) = abs( load_dum( NOPP(i) + NTs ) )
     end if
  end do
  Rsimax = 0.0_rk
  Retamax = 0.0_rk
  Rumax = 0.0_rk
  Rvmax = 0.0_rk
  Rtmax = 0.0_rk
  Rpmax = 0.0_rk
  Rcmax = 0.0_rk
  imaxRsi = 0
  imaxReta = 0
  imaxRu = 0
  imaxRv = 0
  imaxRt = 0
  imaxRp = 0
  imaxRc = 0
  do i = 1, NTN
     if( Rsi(i) .gt. Rsimax ) then
        Rsimax = Rsi(i)
        imaxRsi = i
     end if
     if( Reta(i) .gt. Retamax ) then
        Retamax = Reta(i)
        imaxReta = i
     end if
     if( Ru(i) .gt. Rumax ) then
        if(BCflagN(i,3).eq.1) cycle
        Rumax = Ru(i)
        imaxRu = i
     end if
     if( Rv(i) .gt. Rvmax ) then
        if(BCflagN(i,3).eq.1) cycle
        Rvmax = Rv(i)
        imaxRv = i
     end if
     if( Rt(i) .gt. Rtmax ) then
        Rtmax = Rt(i)
        imaxRt = i
     end if
     if( Rp(i) .gt. Rpmax ) then
        Rpmax = Rp(i)
        imaxRp = i
     end if
     if( no_vapor.eq.0 .and. Rc(i) .gt. Rcmax ) then
        Rcmax = Rc(i)
        imaxRc = i
     end if
  end do

  if (timestep.eq.0 .and. step.eq.1 .and. ( size_function_change.eq.1 .or. read_coordinate_value.eq.1 ) ) then
     open(unit = 12, file = trim(folder)//'error_Res.dat', status = 'replace')
     write(12, '(A)') 'variables = "Rsi", "Reta", "Ru", "Rv", "Rt", "Rp", "Rc" '
  else
     open(unit = 12, file = trim(folder)//'error_Res.dat', status = 'old', access = 'append')
  end if
  if(step.eq.1) then
     write(12, '(A)') ' '
     write(12, '(A, 2es15.7, i8, A)') '----------------', time, dt, timestep, '---------------------'
     write(12, '(A)') ' '
  end if
  write(12,'(7es12.5)') Rsimax, Retamax, Rumax, Rvmax, Rtmax, Rpmax, Rcmax
  if(write_node.eq.1) write(12,'(7i12)') imaxRsi, imaxReta, imaxRu, imaxRv, imaxRt, imaxRp, imaxRc
  write(12,'(es14.7)') error1
  write(12,'(A)') ' '
  close(12)
  
  !write calculation time
  if (timestep.eq.0 .and. step.eq.1 .and. ( size_function_change.eq.1 .or. read_coordinate_value.eq.1 ) ) then
     open(unit = 20, file = trim(folder)//'cal_time.dat', status = 'replace')
     write(20, '(A)') 'variables = "step", "cal_time" '
  else
     open(unit = 20, file = trim(folder)//'cal_time.dat', status = 'old', access = 'append')
  end if
  if(step.eq.1) then
     write(20, '(A)') ' '
     write(20, '(A, 2es15.7, i8, A)') '----------------', time, dt, timestep, '---------------------'
     write(20, '(A)') ' '
  end if
  write(20,'(i4,f7.3)') step, cal_time
  close(20)


  return
end subroutine L2_error
