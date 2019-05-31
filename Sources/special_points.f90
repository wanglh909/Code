module special_points
  use kind
  use data
  use basis_f
contains

  subroutine turning_points(cut)
    implicit none

    integer(kind=ik):: i,j
    real(kind=rk):: r_change, dudn(3), eta1,eta2,eta3
    real(kind=rk), intent(in):: cut


    if(timestep.gt.0) then
       if(timestep.eq.1) then
          open(unit = 30, file = trim(folder)//'turning.dat', status = 'replace')
          write(30, '(A)') 'variables = "contact angle", "r", "time", "element" '
       else
          open(unit = 30, file = trim(folder)//'turning.dat', status = 'old', access = 'append')
       end if


! if(angle_c_degree.le.40.0_rk) then
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
                do while ( abs(eta1-eta2).gt.0.5e-1_rk )
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
                   if( dudn(3) * dudn(1) .lt. 0.0_rk ) then
                      eta2 = eta3
                   else  !dudn(3) * dudn(2) .le. 0.0_rk
                      eta1 = eta3
                   end if
                end do
                r_change = 0.0_rk
                do j = 1, 3
                   r_change = r_change + rcoordinate(globalNM(i,3*j))*phii_1d(eta1,j)
                end do

                ! if( (1.0_rk-r_change)*tan(angle_c) .gt. cut ) &
                     write(30, '(f9.3,2es15.7,i6)') angle_c_degree, r_change, time, i
                print *, 'r_change for dudn =', r_change, 'element:', i
             end if   !dudn change element

          end if   !interface element
       end do  !element loop


! end if  !angle < 40.0_rk

       close(30)

    end if !timestep>0

    return
  end subroutine turning_points



end module special_points
