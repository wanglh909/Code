
  
subroutine prediction
  use kind
  use data!, only: timestep, time, sol, solp, dt, dtp, soldot, soldotp, soldotpp, solpred, CTJ, eps, trunerr, NVar, NTN, rcoordinate, zcoordinate, step, globalNM, rowNM, columnNM, RegN, Nr, Nz, Nu, Nv, NEM, NTE, FTS


  implicit none

  integer(kind=ik):: i, j, imax, jmax, m
  character(LEN=2) :: var
  real(kind=rk):: dtmin=1.0e-8_rk, dtmax=1.0e-1_rk
  ! real(kind=rk):: rmax, zmax, umax, vmax, Tmax, pmax, cmax

if(diverge.eq.0) then
   
   !define trunerr & soldot (prepare for the next time step)
   if(timestep.ge.FTS) then

      !------------------calculate and write truncate error----------------------------
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
               if(solve_T.eq.1) then
                  if (ABS( sol(NOPP(i)+NT) - solpred(NOPP(i)+NT) )/Tmax .gt.trunerr) then
                     trunerr = ABS( sol(NOPP(i)+NT) - solpred(NOPP(i)+NT) ) /Tmax
                     imax = i
                     var = 'T'
                  end if
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
                  !skip calculating cp change?
                  if(solve_cp.eq.1) then
                     if (ABS( sol(NOPP(i)+Ncp) - solpred(NOPP(i)+Ncp) )/cpmax .gt.trunerr) then
                        trunerr = ABS( sol(NOPP(i)+Ncp) - solpred(NOPP(i)+Ncp) ) /cpmax
                        imax = i
                        var = 'cp'
                     end if
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
            if(VN(i).eq.5 .and. solve_T.eq.1) then
               if (ABS( sol(NOPP(i)+NTs) - solpred(NOPP(i)+NTs) )/Tmax .gt.trunerr) then
                  trunerr = ABS( sol(NOPP(i)+NTs) - solpred(NOPP(i)+NTs) ) /Tmax
                  imax = i
                  var = 'T'
               end if
            end if
            if( no_vapor.eq.0 .and. ( VN(i).eq.1 .or. VN(i).eq.2 ) ) then
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
      ! if( VN(imax).ne.1 .and. jmax.eq.NT+1 .and. solve_T.eq.1) var = 'T'
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

   end if   !timestep.gt.FTS
   
end if   !diverge=0
!--------------------------------------------------------------------------------


!for the first timestep, solp was not defined other than initialization
if(timestep.eq.0) solp = sol

!define soldot
! because after the last loop of newton's method for one time step, soldot wasn't calculated
! in order to use 'soldotp = soldot' in the next time step, soldotp need to be calculated one more time.
if(diverge.eq.0) then
   if (timestep.le.FTS) then
      soldot = (sol - solp)/dt
   else
      soldot = 2.0_rk*( sol - solp )/dt - soldotp
   end if
end if   !diverge=0


  !-------------------------------------next timestep----------------------------------------------------

  timestep = timestep + 1

  if(timestep_stable.eq.1 .and. timestep_stable.le.FTS) then
     timestep_stable = timestep_stable + 1
  end if

  !--------------------------------define dt & solpred-------------------------------
  if(diverge.eq.0) then
     dtp = dt
     solp = sol
     !dt. for first 5 steps, dt doesn't need to be redefined; after 5 steps:
     if (timestep.gt.FTS .and. (timestep_stable.eq.0 .or. timestep_stable.gt.FTS)) then
        dt = dtp*( eps/trunerr )**(1.0_rk/3.0_rk)
        if(dt.lt.dtmin) dt = dtmin
        !if(dt.gt.dtmax) dt = dtmax
     end if
     ! !wirte dt in file
     ! if (timestep.eq.FTS .and. step.eq.0) then
     !    open(unit = 10, file = trim(folder)//'dt.dat', status = 'replace')
     !    write(10, '(A)') 'variables = "time", "dt", "timestep" '
     ! else
     !    open(unit = 10, file = trim(folder)//'dt.dat', status = 'old', access = 'append')
     ! end if
     !    write(10, '(2es15.7, i8)') time, dt, timestep
     ! close(10)


     !define solpred, viz. the initial guess for each time step
     soldotpp = soldotp
     soldotp = soldot
     if (timestep.le.FTS .or. (timestep_stable.gt.0 .and. timestep_stable.le.FTS) ) then
        solpred = solp
        !solpred = solp
     else
        solpred = solp + 0.5_rk*dt*( ( 2.0_rk + dt/dtp )*soldotp - dt/dtp*soldotpp )
        !solpred = solp + dt*soldotp
     end if

  end if!diverge=0
  !-----------------------------------------------------------------------------------

  time = time + dt
  if(diverge.eq.0) then
     write(*,*) 'time:', time, '    dt:', dt, '    timestep:', timestep
  else
     write(*,*) 'time:', time, '    dt:', dt, '    timestep:', timestep, '(repeat)'
  end if
  
  sol = solpred
  call split_sol
  !u v r z at boundary is already 0

  !define CTJ: the coefficient of time term in Jac. CTJ will be used when difining Jac
  if (timestep.le.FTS .or. (timestep_stable.gt.0 .and. timestep_stable.le.FTS) ) then
     CTJ = 1.0_rk
  else
     CTJ = 2.0_rk
  end if

  

  return
end subroutine prediction
