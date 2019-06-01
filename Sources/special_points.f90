module special_points
  use kind
  use data
  use basis_f
contains

  subroutine turning_points(cut, delta)
    implicit none

    integer(kind=ik):: i,j
    real(kind=rk):: r_change, dudn(3), eta1,eta2,eta3
    real(kind=rk), intent(in):: cut, delta


    if(timestep.gt.0) then
       if(timestep.eq.1) then
          open(unit = 30, file = trim(folder)//'turning.dat', status = 'replace')
          write(30, '(A)') 'variables = "contact angle", "r", "time", "element", "eta" '
       else
          open(unit = 30, file = trim(folder)//'turning.dat', status = 'old', access = 'append')
       end if


if(angle_c_degree.le.40.0_rk) then
       r_change = 0.0_rk
       do i = 1, NTE

          !the element where direction changes
          if( ( BCflagE(i,2).eq.1 .or. BCflagE(i,2).eq.2 ) .and. VE(i).eq.0 ) then !interface element
             !xi or eta = 0
             if(BCflagE(i,2).eq.1) then  !region 1 interface element
                dudn(1) = -3.0_rk*usol(globalNM(i,1)) + 4.0_rk*usol(globalNM(i,4)) - usol(globalNM(i,7))  !at xi=0
             else  !region 3 interface element
                dudn(1) = -3.0_rk*usol(globalNM(i,1)) + 4.0_rk*usol(globalNM(i,2)) - usol(globalNM(i,3))  !at eta=0
             end if

             !xi or eta = 1
             if(BCflagE(i,2).eq.1) then  !region 1 interface element
                dudn(2) = -3.0_rk*usol(globalNM(i,3)) + 4.0_rk*usol(globalNM(i,6)) - usol(globalNM(i,9))  !at xi=0
             else  !region 3 interface element
                dudn(2) = -3.0_rk*usol(globalNM(i,7)) + 4.0_rk*usol(globalNM(i,8)) - usol(globalNM(i,9))  !at eta=0
             end if

             ! !this will be useful if actually do du/dz
             !            !at eta = 0
             !            retap(1) = -3.0_rk*rcoordinate(globalNM(i,3)) + 4.0_rk*rcoordinate(globalNM(i,6)) - rcoordinate(globalNM(i,9))
             !            zetap(1) = -3.0_rk*zcoordinate(globalNM(i,3)) + 4.0_rk*zcoordinate(globalNM(i,6)) - zcoordinate(globalNM(i,9))

             !            !at eta = 1
             !            retap(2) = rcoordinate(globalNM(i,3)) - 4.0_rk*rcoordinate(globalNM(i,6)) + 3.0_rk*rcoordinate(globalNM(i,9))
             !            zetap(2) = zcoordinate(globalNM(i,3)) - 4.0_rk*zcoordinate(globalNM(i,6)) + 3.0_rk*zcoordinate(globalNM(i,9))

             !dudn direction change
             if( dudn(1) * dudn(2) .lt. 0.0_rk .and. i.ne.CL_element ) then
                ! write(*,*)  dudn(1), dudn(2)
                eta1 = 0.0_rk !this stands for xi in region 1
                eta2 = 1.0_rk
                do while ( abs(eta1-eta2).gt.delta )
                   eta3 = (eta1+eta2)/2.0_rk
                   dudn(3) = 0.0_rk
                   do j = 1, 9
                      if(BCflagE(i,2).eq.1) then  !region 1 interface element
                         dudn(3) = dudn(3) + usol(globalNM(i,j))*phiieta(eta3,0.0_rk,j)
                         !eta3 stands for xi in region1
                      else  !region 3 interface element
                         dudn(3) = dudn(3) + usol(globalNM(i,j))*phiisi(0.0_rk,eta3,j)
                      end if
                   end do
! !debug
! ! if(angle_c_degree.lt.12.2_rk .and. angle_c_degree.gt.12.1_rk) then
!    print *, eta1, eta2, eta3
!    print *, dudn(1), dudn(2), dudn(3)
!    pause
! ! end if
                   if( dudn(3) * dudn(1) .lt. 0.0_rk ) then
                      eta2 = eta3
                   else  !dudn(3) * dudn(2) .le. 0.0_rk
                      eta1 = eta3
                      dudn(1) = dudn(3)   !update dudn(1) to be with the new eta1, because the criterion is to compare with dudn(1).
                   end if                   
                end do !to reach eta1 eta2 very close
                r_change = 0.0_rk
                do j = 1, 3
                   if(BCflagE(i,2).eq.1) then  !region 1 interface element
                      r_change= r_change + rcoordinate(globalNM(i,j))*phii_1d(eta1,j)
                   else  !region 3 interface element
                      r_change= r_change + rcoordinate(globalNM(i,3*j-2))*phii_1d(eta1,j)
                   end if
                end do

                if( (1.0_rk-r_change)*tan(angle_c) .gt. cut ) &
                     write(30, '(f9.3,2es15.7,i6,es15.7)') angle_c_degree, r_change, time, i, eta1
                print *, 'r_change for dudn =', r_change, 'element:', i, 'eta:', eta1
             end if   !dudn change element

          end if   !interface element
       end do  !element loop


end if  !angle < 40.0_rk

       close(30)

    end if !timestep>0

    return
  end subroutine turning_points



  




subroutine stagnation_and_extremum(cut, delta)
    implicit none

    integer(kind=ik):: i,j, flag
    real(kind=rk):: r_change, eta1, eta2, eta3, usolp, vsolp, zsolp
    real(kind=rk):: Teta(3), gradT(3), v_surf_p(3), dPdr(3), retap(3),zetap(3), v_surf(3)
    real(kind=rk), intent(in):: cut, delta

  !--------velocity direction change location on free surface & grad(T) direction------
  !we only care about stagnation points and extremum points:
  !1. after initial chaos are passed. viz. angle<40
  !2. not too close to CL. viz. 1.0-r>cut
  if(timestep.gt.0) then
   
     if(solve_T.eq.1) then
        if(timestep.eq.1) then
           open(unit = 18, file = trim(folder)//'extremum.dat', status = 'replace')
           write(18, '(A)') 'variables = "contact angle", "r", "time", "element", "eta" '
        else
           open(unit = 18, file = trim(folder)//'extremum.dat', status = 'old', access = 'append')
        end if
     end if
     
     ! if(timestep.eq.1) then
     !    open(unit = 30, file = trim(folder)//'v0_lubrication.dat', status = 'replace')
     !    write(30, '(A)') 'variables = "contact angle", "r", "time", "element" '
     ! else
     !    open(unit = 30, file = trim(folder)//'v0_lubrication.dat', status = 'old', access = 'append')
     ! end if


     if(timestep.eq.1) then
        open(unit = 14, file = trim(folder)//'stagnation.dat', status = 'replace')
        write(14, '(A)') 'variables = "contact angle", "r", "time", "element", "eta" '
     else
        open(unit = 14, file = trim(folder)//'stagnation.dat', status = 'old', access = 'append')
     end if

     
if(angle_c_degree.le.40.0_rk) then
     r_change = 0.0_rk
     do i = 1, NTE
        if(BCflagE(i,3).eq.1) then  !surface element
           !the element where direction changes

           !at eta = 0
           retap(1) = -3.0_rk*rcoordinate(globalNM(i,3)) + 4.0_rk*rcoordinate(globalNM(i,6)) - rcoordinate(globalNM(i,9))
           zetap(1) = -3.0_rk*zcoordinate(globalNM(i,3)) + 4.0_rk*zcoordinate(globalNM(i,6)) - zcoordinate(globalNM(i,9))
           v_surf(1) = ( usol(globalNM(i,3))*retap(1) + vsol(globalNM(i,3))*zetap(1) ) /sqrt( retap(1)**2+zetap(1)**2 )
           if(solve_T.eq.1) then
              Teta(1) = -3.0_rk*Tsol(globalNM(i,3)) + 4.0_rk*Tsol(globalNM(i,6)) - Tsol(globalNM(i,9))
              gradT(1) = -Teta(1)/sqrt( retap(1)**2+zetap(1)**2 )
           end if
           dPdr(1) = ( psol(globalNM(i,9)) - psol(globalNM(i,3)) ) /retap(1)
           v_surf_p(1) = -0.5_rk*dPdr(1)*zcoordinate(globalNM(i,3)) + beta*gradT(1)

           !at eta = 1
           retap(2) = rcoordinate(globalNM(i,3)) - 4.0_rk*rcoordinate(globalNM(i,6)) + 3.0_rk*rcoordinate(globalNM(i,9))
           zetap(2) = zcoordinate(globalNM(i,3)) - 4.0_rk*zcoordinate(globalNM(i,6)) + 3.0_rk*zcoordinate(globalNM(i,9))
           v_surf(2) = ( usol(globalNM(i,9))*retap(2) + vsol(globalNM(i,9))*zetap(2) ) /sqrt( retap(2)**2+zetap(2)**2 )
           if(solve_T.eq.1) then
              Teta(2) = Tsol(globalNM(i,3)) - 4.0_rk*Tsol(globalNM(i,6)) + 3.0_rk*Tsol(globalNM(i,9))
              gradT(2) = -Teta(2)/sqrt( retap(2)**2+zetap(2)**2 )
           end if
           dPdr(2) = ( psol(globalNM(i,9)) - psol(globalNM(i,3)) ) /retap(2)
           v_surf_p(2) = -0.5_rk*dPdr(2)*zcoordinate(globalNM(i,9)) + beta*gradT(2)

           ! !lubrication velocity direction change
           ! if( angle_c_degree.lt.15.0_rk .and. &
           !      v_surf_p(1) * v_surf_p(2) .lt. 0.0_rk .and. &
           !      (i.ne.top_element .and. i.ne.CL_element) ) then 
 
           !    eta1 = 0.0_rk
           !    eta2 = 1.0_rk
           !    do while ( abs(eta1-eta2).gt.0.5e-1_rk )
           !       eta3 = (eta1+eta2)/2.0_rk
           !       retap(3) = 0.0_rk
           !       zetap(3) = 0.0_rk
           !       if(solve_T.eq.1) Teta(3) = 0.0_rk
           !       zsolp = 0.0_rk
           !       do j = 1, 3
           !          retap(3) = retap(3) + rcoordinate(globalNM(i,3*j))*phiix_1d(eta3,j)
           !          zetap(3) = zetap(3) + zcoordinate(globalNM(i,3*j))*phiix_1d(eta3,j)
           !          if(solve_T.eq.1) Teta(3) = Teta(3) + Tsol(globalNM(i,3*j))*phiix_1d(eta3,j)
           !          zsolp = zsolp + zcoordinate(globalNM(i,3*j))*phii_1d(eta3,j)
           !       end do
           !       if(solve_T.eq.1) gradT(3) = -Teta(3)/sqrt( retap(3)**2+zetap(3)**2 )
           !       dPdr(3) = ( psol(globalNM(i,9)) - psol(globalNM(i,3)) ) /retap(3)
           !       v_surf_p(3) = -0.5_rk*dPdr(3)*zsolp + beta*gradT(3)
           !       if( v_surf_p(3) * v_surf_p(1) .lt. 0.0_rk ) then
           !          eta2 = eta3
           !       else  !v_surf(3) * v_surf(2) .le. 0.0_rk
           !          eta1 = eta3
           !       end if
           !    end do
           !    r_change = 0.0_rk
           !    do j = 1, 3
           !       r_change = r_change + rcoordinate(globalNM(i,3*j))*phii_1d(eta1,j)
           !    end do

           !    write(30, '(f9.3,2es15.7, i6)') angle_c_degree, r_change, time, i
           !    !write(*,*) 'r_change for lubrication velocity =', r_change, 'element:', i
           !    ! exit   !?not strict

           ! end if   !u change element



           !surface flow direction change
           if( v_surf(1) * v_surf(2) .lt. 0.0_rk .and. (i.ne.top_element .and. i.ne.CL_element) ) then
! write(*,*)  v_surf(1), v_surf(2)
              eta1 = 0.0_rk
              eta2 = 1.0_rk
              do while ( abs(eta1-eta2).gt.delta )
                 eta3 = (eta1+eta2)/2.0_rk
                 retap(3) = 0.0_rk
                 zetap(3) = 0.0_rk
                 usolp = 0.0_rk
                 vsolp = 0.0_rk
                 do j = 1, 3
                    retap(3) = retap(3) + rcoordinate(globalNM(i,3*j))*phiix_1d(eta3,j)
                    zetap(3) = zetap(3) + zcoordinate(globalNM(i,3*j))*phiix_1d(eta3,j)
                    usolp = usolp + usol(globalNM(i,3*j))*phii_1d(eta3,j)
                    vsolp = vsolp + vsol(globalNM(i,3*j))*phii_1d(eta3,j)
                 end do
                 v_surf(3) = ( usolp*retap(3) + vsolp*zetap(3) ) /sqrt( retap(3)**2+zetap(3)**2 )
                 if( v_surf(3) * v_surf(1) .lt. 0.0_rk ) then
                    eta2 = eta3
                 else  !v_surf(3) * v_surf(2) .le. 0.0_rk
                    eta1 = eta3
                    v_surf(1) = v_surf(3)   !update v_surf(1) to be with the new eta1, because the criterion is to compare with v_surf(1).
                 end if
              end do
              r_change = 0.0_rk
              do j = 1, 3
                 r_change = r_change + rcoordinate(globalNM(i,3*j))*phii_1d(eta1,j)
              end do

              if( (1.0_rk-r_change)*tan(angle_c) .gt. cut ) &
              write(14, '(f9.3,2es15.7,i6,es15.7)') angle_c_degree, r_change, time, i, eta1
              !write(*,*) 'r_change for velocity =', r_change, 'element:', i
              ! exit   !?not strict

           end if   !u change element


           !grad(T) direction change
           if(solve_T.eq.1) then
              if( gradT(1) * gradT(2) .lt. 0.0_rk .and. (i.ne.top_element .and. i.ne.CL_element) ) then  
                 eta1 = 0.0_rk
                 eta2 = 1.0_rk
                 do while ( abs(eta1-eta2).gt.delta )
                    eta3 = (eta1+eta2)/2.0_rk
                    retap(3) = 0.0_rk
                    zetap(3) = 0.0_rk
                    Teta(3) = 0.0_rk
                    do j = 1, 3
                       retap(3) = retap(3) + rcoordinate(globalNM(i,3*j))*phiix_1d(eta3,j)
                       zetap(3) = zetap(3) + zcoordinate(globalNM(i,3*j))*phiix_1d(eta3,j)
                       Teta(3) = Teta(3) + Tsol(globalNM(i,3*j))*phiix_1d(eta3,j)
                    end do
                    gradT(3) = -Teta(3)/sqrt( retap(2)**2+zetap(2)**2 )
                    if( gradT(3) * gradT(1) .lt. 0.0_rk ) then
                       eta2 = eta3
                    else  !v_surf(3) * v_surf(2) .le. 0.0_rk
                       eta1 = eta3
                       gradT(1) = gradT(3)   !update gradT(1) to be with the new eta1, because the criterion is to compare with gradT(1).
                    end if
                 end do
                 r_change = 0.0_rk
                 do j = 1, 3
                    r_change = r_change + rcoordinate(globalNM(i,3*j))*phii_1d(eta1,j)
                 end do

                 
                 if( (1.0_rk-r_change)*tan(angle_c) .gt. cut ) &
                 write(18, '(f9.3,2es15.7,i4,es15.7)') angle_c_degree, r_change, time, i, eta1
                 !write(*,*) 'r_change for gradT =', r_change, 'element:', i
                 ! exit   !?not strict

              end if   !gradT change element
           end if



        end if   !surface element
     end do  !element loop

     ! do i = 1, NTN
     !    if(BCflagN(i,3).ne.1) cycle
     !    do j = i+1, NTN
     !       if(BCflagN(j,3).eq.1) exit
     !    end do
     !    if( usol(i)*usol(j).lt.0.0_rk ) then
     !       r_change = ( rcoordinate(i) + rcoordinate(j) )/2.0_rk 
     !       write(14, '(3es15.7)') time, angle_c_degree, r_change
     !       write(*,*) 'r_change =', r_change
     !       exit   !?not strict
     !    end if
     ! end do
     
end if  !angle < 40.0_rk

     !check initial stability for init_stability, could be integrated with direction change calculation
     if(init_stability.eq.0) then
        flag = 0
        do i = 1, NTN
           if(BCflagN(i,3).ne.1) cycle
           do j = i+1, NTN
              if(BCflagN(j,3).eq.1) exit
           end do
           if( usol(i)*usol(j).lt.0.0_rk .and. &
                (1.0_rk-rcoordinate(i))*(1.0_rk-rcoordinate(j)).gt.1.0e-6_rk ) then
              write(*,*) 'stag exist', rcoordinate(i), rcoordinate(j)
              flag = 1
              exit   !?not strict
           end if
        end do
        if(flag.eq.0) then
           init_stability = 1
           dt = 1.0e-7_rk
           timestep_stable = 1
           !graph_mode = 1
        end if
     end if
     
     ! close(30)
     close(14)
     close(18)
  
  end if !timestep>0


  return
end subroutine stagnation_and_extremum




  

end module special_points
