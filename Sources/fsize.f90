subroutine vapor_mesh_size
  use kind
  use data, only: vlayer, NEV1, R11, R, NEV, no_vapor
  
  vlayer = 3
  allocate( NEV1(vlayer-1), R11(vlayer-1) )
  NEV1 = (/ NEV/6+1, 2*(NEV/6)+2 /)
  if(no_vapor.eq.1) then
     R11 = (/ 1.05_rk, 1.1_rk /)*R  
  else
     R11 = (/ 1.05_rk, 2.0_rk /)*R
  end if

  return
end subroutine vapor_mesh_size


subroutine size_function!(size_function_change)  !? not working with vapor phase
  use kind
  use data, only: fsi_size, geta_size, final_size, size_function_change, &
       rowNM, columnNM, RegN, NTE, VE, NEV, NEV1, NEL, layer, angle_c, angle_cd, pi, &
       k_alge, NEM, NES, NEM_alge, simple_mesh, alge_corner, folder, no_vapor, read_coordinate_value

  implicit none
  integer(kind=ik)::i,j, group
  real(kind=rk):: x1, x2, x3, x4, x5, x6, x7, x8, x9, y, x_change  !the value of fsi_size & geta_size at the most dense area
  real(kind=rk):: a  !temperary value

  !integer(kind=ik), intent(in):: size_function_change
  

!   if(read_coordinate_value.eq.1) then
!      size_function_change = 3
!      angle_cd = angle_c
!   end if

!   y = 1.0_rk
!   x3 = 0.3_rk
!   x4 = 0.4_rk
!   x5 = 0.5_rk
!   x6 = 0.6_rk
!   x7 = 0.7_rk
!   x8 = 0.8_rk
!   x9 = 0.9_rk

!   if(size_function_change.eq.0) then
!      angle_cd = angle_c
!      fsi_size = 0.0_rk
!      geta_size = 0.0_rk
!      x_change = y
!      angle_c = pi/2.0_rk
!      final_size = 0
!      if(simple_mesh.eq.1) final_size = 1
!   else if(size_function_change.eq.1) then
!      final_size = 0
!      x_change = x7
!      ! k_alge = 3.0_rk
!      ! if(angle_cd.lt.pi/8.0_rk) then   !22.5 degree
!      !    angle_c = pi/8.0_rk
!      ! else
!      !    angle_c = angle_cd
!      ! end if
!   else if(size_function_change.eq.2) then
!      final_size = 0
!      x_change = x7
!      ! k_alge = 3.0_rk
!      if(angle_cd.lt.pi/8.0_rk) then   !22.5 degree
!         angle_c = pi/8.0_rk
!      else
!         angle_c = angle_cd
!      end if
!   else if( size_function_change.eq.3 ) then
!      x_change = x7
!      angle_c = angle_cd
!      k_alge = 3.0_rk
!      ! if(angle_cd.lt.pi/8.0_rk) then   !22.5 degree
!      !    angle_c = pi/8.0_rk
!      ! else
!      !    angle_c = angle_cd
!      ! end if
!      if(angle_c.eq.angle_cd) final_size = 1
! print *, 'final size'
!   ! else if(size_function_change.eq.3) then
!   !    x_change = 0.5_rk
!   !    if(angle_cd.ge.pi/12.0_rk) then
!   !       angle_c = angle_cd
!   !       final_size = 1
!   !    else
!   !       angle_c = pi/12.0_rk
!   !       ! k_alge = 5.0_rk
!   !    end if
!   ! else if(size_function_change.eq.4) then
!   !    x_change = 0.5_rk
!   !    if(angle_cd.ge.pi/18.0_rk) then
!   !       angle_c = angle_cd
!   !       final_size = 1
!   !    else
!   !       angle_c = pi/18.0_rk
!   !       k_alge = 5.0_rk
!   !    end if
!   ! else if(size_function_change.eq.5) then
!   !    x_change = 0.5_rk
!   !    angle_c = angle_cd
!   !    !    final_size = 0
!   !    ! else
!   !    !    write(*,*) 'change6'
!   !    !    x = 0.1_rk
!   !    !    y = 0.5_rk
!   !    !    x_change = 1.0_rk
!   !    final_size = 1
!   end if


!   do i = 1, NTE

!      if(VE(i).eq.0) then

!         if(RegN(i).eq.1) then
!            fsi_size(i) = (y-x6)*real(rowNM(i),rk)/real(NEL-1,rk) + x6 - (y-x6)/real(NEL-1,rk)
!            geta_size(i) = (y-x6)*real(columnNM(i),rk)/real(NEL-1,rk) + x6 - (y-x6)/real(NEL-1,rk)
!         else   !region 3&4
!            if(alge_corner.eq.1) then
!               if(rowNM(i).gt.NEM_alge) then
!                  geta_size(i) = (x_change-y)*real(rowNM(i)-NEM,rk)/real(NEM_alge+1-NEM,rk) + y 
!               else
!                  geta_size(i) = 1.0_rk
!               end if
!            else
!               geta_size(i) = 1.0_rk
!            end if
!            if (RegN(i).eq.3) then
!               fsi_size(i) = (y-x8)*real(columnNM(i),rk)/real(NEL-1,rk) + x8 - (y-x8)/real(NEL-1,rk)
!            else    !RegN = 4
!               fsi_size(i) = (y-x4)*real(columnNM(i),rk)/real(1-NEL,rk) + x4 - 2.0_rk*(y-x4)*real(NEL,rk)/real(1-NEL,rk)
!            end if
!         end if

!      else if(VE(i).eq.1) then  !VE = 1
!         if( RegN(i).eq.2 .and. rowNM(i).le.NEV1(1) ) then!??
!            a = (y-x6)/real(NEV1(1)-1,rk)*real(rowNM(i)-1,rk) + x6
!            geta_size(i) = (y-a)*real(columnNM(i),rk)/real(NEL-1,rk) + a - (y-a)/real(NEL-1,rk)
!         else if(RegN(i).eq.5 .and. rowNM(i).gt.NEM_alge .and. alge_corner.eq.1) then
!            ! a = (y-x_change)/real(NEV1(1)-1,rk)*real(columnNM(i)-1,rk) + x_change
!            geta_size(i) = (x_change-y)/real(NEM_alge+1-NEM,rk)*real(rowNM(i)-NEM,rk) + y
!         else
!            geta_size(i) = 1.0_rk
!         end if
        
!         fsi_size(i) = 1.0_rk
!         ! if(RegN(i).eq.2) then
!         !    if( rowNM(i).le.NEV-NEV1(2) ) then
!         !       fsi_size(i) = (y-x5)/real(1-NEV+NEV1(2),rk)*real(rowNM(i)-1,rk) + y 
!         !    else if( rowNM(i).le.NEV-NEV1(1) ) then
!         !       fsi_size(i) = (y-x5)/real(NEV1(1)-NEV1(2)+1,rk)*real(rowNM(i)-(NEV-NEV1(2)+1),rk) + y 
!         !    else
!         !       fsi_size(i) = 1.0_rk
!         !    end if      
!         ! else
!         !    if( rowNM(i).gt.NEV1(2) ) then
!         !       fsi_size(i) = (y-x5)/real(NEV1(2)-NEV1(1)-1,rk)*real(rowNM(i)-NEV1(1)-1,rk) + y 
!         !    else if( rowNM(i).gt.NEV1(1) ) then
!         !       fsi_size(i) = (y-x5)/real(NEV-NEV1(2)-1,rk)*real(rowNM(i)-(NEV-NEV1(2)-1),rk) + y 
!         !    else
!         !       fsi_size(i) = 1.0_rk
!         !    end if

!         ! end if
        
!         if(NEV1(1).eq.1) then
!            geta_size(i) = 1.0_rk
!            fsi_size(i) = 1.0_rk
!         end if

!      else   !VE = 5
!         if(RegN(i).eq.6) then
! if(alge_corner.eq.1) then
!            if(rowNM(i).le.NEM-NEM_alge) then
!               geta_size(i) = (x_change-y)*real(rowNM(i)-1,rk)/real(NEM-NEM_alge-1,rk) + y 
!            else
!               geta_size(i) = 1.0_rk
!            end if
! else
!    geta_size(i) = 1.0_rk
! end if
!            ! geta_size(i) = 1.0_rk
!            fsi_size(i) = (x5-y)/real(NES-1,rk)*real(columnNM(i)-1,rk) + y
!         else  !RegN = 7,8
!            geta_size(i) = (x5-y)/real(NES-1,rk)*real(columnNM(i)-1,rk) + y
!            if(RegN(i).eq.7) then
!               fsi_size(i) = (x6-y)/real(NEL-1,rk)*real(rowNM(i)-NEM-1,rk) + y
!               ! fsi_size(i) = 1.0_rk
!            else !reg 8
!               if(no_vapor.eq.0)  then
!                  fsi_size(i) = 1.0_rk
!               else  !no_vapor = 1
!                  fsi_size(i) = (x3-y)/real(1-NEV,rk)*real(rowNM(i)-(NEM+NEL+NEV),rk) + y
!               end if
!            end if

!         end if

!      end if




!      ! if(RegN(i).eq.1) then

!      !    do group = 1, NEL
!      !       if( i.gt.(group-1)**2 .and. i.le.group**2 ) then
!      !          fsi_size(i) = (y-x6)*real(group,rk)/real(NEL-1,rk) + x6 - (y-x6)/real(group,rk)
!      !          geta_size(i) = fsi_size(i)
!      !       end if
!      !    end do

!      ! else
!      !    geta_size(i) = 1.0_rk!(x_change-y)*real(rowNM(i),rk)/real(NEM-1,rk) + y - (x_change-y)/real(NEM-1,rk)

!      !    if (RegN(i).eq.2) then
!      !       fsi_size(i) = y!(y-x6)*real(columnNM(i),rk)/real(NEL-1,rk) + x6 - (y-x6)/real(NEL-1,rk)
!      !    else
!      !       fsi_size(i) = y!(y-x6)*real(columnNM(i),rk)/real(1-NEL,rk) + x6 - 2.0_rk*(y-x6)*real(NEL,rk)/real(1-NEL,rk)
!      !    end if

!      ! end if


!   end do

!   ! open(unit=10, file=trim(folder)//'size.dat', status = 'replace')
!   ! do i = 1, NTE
!   !    write(10, '(i8, 2es15.7)') i, fsi_size(i), geta_size(i)
!   ! end do
!   ! close(10)
!   ! ! pause


  ! allocate( fsi_size(NTE), geta_size(NTE) )
  ! fsi_size = 0.0_rk
  ! geta_size = 0.0_rk
  fsi_size(:) = 1.0_rk
  geta_size(:) = 1.0_rk
  final_size = 1


  return
end subroutine size_function

