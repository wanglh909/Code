
  
subroutine prediction
  use kind
  use data!, only: timestep, time, sol, solp, dt, dtp, soldot, soldotp, soldotpp, solpred, CTJ, eps, trunerr, NVar, NTN, rcoordinate, zcoordinate, step, globalNM, rowNM, columnNM, RegN, Nr, Nz, Nu, Nv, NEM, NTE, FTS


  implicit none

  integer(kind=ik):: i, j, imax, jmax, m
  character(LEN=2) :: var
  ! real(kind=rk):: rmax, zmax, umax, vmax, Tmax, pmax, cmax

if(diverge.eq.0) then

  !prepare for the next time step
  if(timestep.ge.FTS) then
     
     ! !find maximum value of each variable, used to normalize error
     ! rmax = 0.0_rk
     ! zmax = 0.0_rk
     ! umax = 0.0_rk
     ! vmax = 0.0_rk
     ! Tmax = 0.0_rk
     ! pmax = 0.0_rk
     ! cmax = 0.0_rk
     ! do i = 1, NTN
     !    if( abs(rcoordinate(i)) .gt. rmax ) rmax = abs(rcoordinate(i))
     !    if( abs(zcoordinate(i)) .gt. zmax ) zmax = abs(zcoordinate(i))
     !    if( abs(usol(i)) .gt. umax ) umax = abs(usol(i))
     !    if( abs(vsol(i)) .gt. vmax ) vmax = abs(vsol(i))
     !    if( abs(Tsol(i)) .gt. Tmax ) Tmax = abs(Tsol(i))
     !    if( abs(psol(i)) .gt. pmax ) pmax = abs(psol(i))
     !    if( abs(csol(i)) .gt. cmax ) cmax = abs(csol(i))
     ! end do

  
     !calculate and write truncate error
     trunerr = 0.0_rk     !Linf norm
     do i = 1, NTN
        if (ABS( sol(NOPP(i)+Nr) - solpred(NOPP(i)+Nr) )/rmax .gt.trunerr) then
           trunerr = ABS( sol(NOPP(i)+Nr) - solpred(NOPP(i)+Nr) ) /rmax
           imax = i
           var = 'r'
        end if
        if(NTs.eq.2 .or. VN(i).ne.5)  then
           if (ABS( sol(NOPP(i)+Nz) - solpred(NOPP(i)+Nz) )/zmax .gt.trunerr) then
              trunerr = ABS( sol(NOPP(i)+Nz) - solpred(NOPP(i)+Nz) ) /zmax
              imax = i
              var = 'z'
           end if
        end if
        if(s_mode.eq.0) then
           if( VN(i).eq.0 .or. VN(i).eq.2 .or. BCflagN(i,2).eq.1) then
              if (ABS( sol(NOPP(i)+NT) - solpred(NOPP(i)+NT) )/Tmax .gt.trunerr) then
                 trunerr = ABS( sol(NOPP(i)+NT) - solpred(NOPP(i)+NT) ) /Tmax
                 imax = i
                 var = 'T'
              end if
              if( VN(i).eq.0 .or. VN(i).eq.2 ) then
                 if (ABS( sol(NOPP(i)+Nu) - solpred(NOPP(i)+Nu) )/umax .gt.trunerr) then
                    trunerr = ABS( sol(NOPP(i)+Nu) - solpred(NOPP(i)+Nu) ) /umax
                    imax = i
                    var = 'u'
                 end if
                 if (ABS( sol(NOPP(i)+Nv) - solpred(NOPP(i)+Nv) )/vmax .gt.trunerr) then
                    trunerr = ABS( sol(NOPP(i)+Nv) - solpred(NOPP(i)+Nv) ) /vmax
                    imax = i
                    var = 'v'
                 end if
                 ! if(PN(i).eq.1) then    ! skip calculating P change
                 !    if (ABS( sol(NOPP(i)+Np) - solpred(NOPP(i)+Np) )/pmax .gt.trunerr) then
                 !       trunerr = ABS( sol(NOPP(i)+Np) - solpred(NOPP(i)+Np) ) /pmax
                 !       imax = i
                 !       var = 'p'
                 !    end if
                 ! end if
              end if
           end if
           if(VN(i).eq.5) then
              if (ABS( sol(NOPP(i)+NTs) - solpred(NOPP(i)+NTs) )/Tmax .gt.trunerr) then
                 trunerr = ABS( sol(NOPP(i)+NTs) - solpred(NOPP(i)+NTs) ) /Tmax
                 imax = i
                 var = 'T'
              end if
           end if
           if(VN(i).eq.1 .or. VN(i).eq.2) then
              if (ABS( sol(NOPP(i)+MDF(i)-1) - solpred(NOPP(i)+MDF(i)-1) )/cmax .gt.trunerr) then
                 trunerr = ABS( sol(NOPP(i)+MDF(i)-1) - solpred(NOPP(i)+MDF(i)-1) ) /cmax
                 imax = i
                 var = 'c'
              end if
           end if
        end if
        ! do j = 1, MDF(i)
        !    if( PN(i).eq.1 .and. j.eq.Np+1 ) cycle   ! skip calculating P change
        !    if (ABS( sol( NOPP(i)+j-1 )  - solpred( NOPP(i)+j-1 ) ) .gt.trunerr) then
        !       trunerr = ABS( sol( NOPP(i)+j-1 ) - solpred( NOPP(i)+j-1 ) )
        !       imax = i
        !       jmax = j
        !    end if
        ! end do
     end do
     ! if(jmax.eq.Nr+1) var = 'r'
     ! if(jmax.eq.Nz+1) var = 'z'
     ! if( VN(imax).ne.1 .and. jmax.eq.Nu+1 ) var = 'u'
     ! if( VN(imax).ne.1 .and. jmax.eq.Nv+1 ) var = 'v'
     ! if( VN(imax).ne.1 .and. jmax.eq.NT+1 ) var = 'T'
     ! if( VN(imax).ne.0 .and. jmax.eq.MDF(imax) ) var = 'c'
     ! ! do i = 1, NVar, 1
     ! !    if (ABS( sol(i) - solpred(i) ).gt.trunerr) trunerr = ABS( sol(i) - solpred(i) )
     ! ! end do
     trunerr = trunerr/3.0_rk/( 1.0_rk + dtp/dt )

     !wirte trunerr in file
     if (timestep.eq.FTS .and. step.eq.0) then
        open(unit = 10, file = trim(folder)//'trun_error.dat', status = 'replace')
        write(10, '(A)') 'variables = "dmax", "variable", "node", "rcoordinate", "zcoordinate" '
     else
        open(unit = 10, file = trim(folder)//'trun_error.dat', status = 'old', access = 'append')
     end if
     if(step.eq.0) then
        write(10, '(A)') ' '
        write(10, '(A, 2es15.7, i8, A)') '----------------', time, dt, timestep, '---------------------'
        write(10, '(A)') ' '
     end if
     write(10,'(es15.7, A, A, i8, 2es15.7)') trunerr, '   ', var, imax, rcoordinate(imax), zcoordinate(imax)
     close(10)



     !define soldot
     ! because after the last loop of newton's method for one time step, soldot wasn't calculated
     ! in order to use 'soldotp = soldot' in the next time step, soldotp need to be calculated one more time. 
     if (timestep.le.FTS) then
        soldot = (sol - solp)/dt
     else
        soldot = 2.0_rk*( sol - solp )/dt - soldotp
     end if

  end if



  !-------------------------------------next timestep----------------------------------------------------

  timestep = timestep + 1
  solp = sol

  !define dt
  dtp = dt
  !for first 5 steps, dt doesn't need to be redefined; after 5 steps, do the following
  if (timestep.gt.FTS) then
     dt = dtp*( eps/trunerr )**(1.0_rk/3.0_rk)
  end if
  time = time + dt
  write(*,*) 'time:', time, '    dt:', dt, '    timestep:', timestep

  soldotpp = soldotp
  soldotp = soldot

  !define solpred, viz. the initial guess for each time step
  if (timestep.le.FTS) then
     solpred = solp
     !solpred = solp
  else
     solpred = solp + 0.5_rk*dt*( ( 2.0_rk + dt/dtp )*soldotp - dt/dtp*soldotpp )
     !solpred = solp + dt*soldotp
  end if

  sol = solpred
  call split_sol
  !u v r z at boundary is already 0

  !define CTJ: the coefficient of time term in Jac. CTJ will be used when difining Jac
  if (timestep.le.FTS) then
     CTJ = 1.0_rk
  else
     CTJ = 2.0_rk
  end if

else  !diverge=1

  timestep = timestep + 1
  time = time + dt
  sol = solpred
  call split_sol

end if

  return
end subroutine prediction
