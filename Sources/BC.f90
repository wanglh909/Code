


subroutine Dirichlet_BC(m, locJac, locRes, LNVar, LNOPP)
  use kind
  use data
  !use front_mod, only: determine_offsets

  implicit none
  integer(kind=ik):: i, j, ipp, n

  integer(kind=ik), intent(in):: m, LNVar, LNOPP(bas)
  real(kind=rk):: locJac(LNVar, LNVar), locRes(LNVar)
  real(kind=rk):: Rp, zp, Rp1(vlayer-1), zp1(vlayer-1), angle_temp(vlayer-1), angle_tempp

  real(kind=rk):: Rpp, mpp, npp, dRp, d, ap, bp, cp, zpp

  angle_temp = (/ angle_c, pi/3.0_rk /)!(/ pi/10.0_rk, pi/3.0_rk /)

  Rp = R/sin(angle_c)
  zp = R/tan(angle_c)
  
  if(no_vapor.eq.0) then
     if( angle_c.gt.angle_temp(2)) then
        angle_temp(:) = angle_c
     end if
  end if

  ! do n = 1, vlayer-1
  !    if(angle_c.gt.angle_temp(2)) then
  !       Rp1(n) = R11(n)/sin(angle_c)
  !       zp1(n) = R11(n)/tan(angle_c)
  !    else
  !       Rp1(n) = R11(n)/sin(angle_temp(2))
  !       zp1(n) = R11(n)/tan(angle_temp(2))
  !       if(n.eq.1 ) then! .and. angle_c.lt.angle_temp(1)) then
  !          Rp1(n) = R11(n)/sin(angle_temp(1))
  !          zp1(n) = R11(n)/tan(angle_temp(1))
  !       end if
  !    end if
  ! end do

  ! if(step.eq.1.and. m.eq.50) then
  ! write(*,*) angle_c
  ! end if


  do i = 1, 9

     !axis
     if( BCflagN( globalNM(m,i), 1 ).eq.1 ) then

        j = LNOPP(i) + Nz
        locJac(j,:) = 0.0_rk
        locJac(j,j-1) = 1.0_rk        !dReta(i)/dr(i)
        locRes(j) = 0.0_rk

        if(s_mode.eq.0) then
           if(VE(m).eq.0) then
              j = LNOPP(i) + Nu
              locJac(j,:) = 0.0_rk
              locJac(j,j) = 1.0_rk
              locRes(j) = 0.0_rk
           end if
        end if

     end if

     !base
     if( BCflagN( globalNM(m,i), 2 ).ne.0 ) then

        if( BCflagN( globalNM(m,i), 2 ).eq.1 .or. BCflagN( globalNM(m,i), 2 ).eq.3 ) then
           j = LNOPP(i) + Nz      !the location of Reta(i) in locRes
           locJac(j,:) = 0.0_rk
           locJac(j,j) = 1.0_rk      !dReta(i)/dz(i)
           locRes(j) = 0.0_rk
        else   !base = 2
           j = LNOPP(i) + Nr
           locJac(j,:) = 0.0_rk
           locJac(j,j+1) = 1.0_rk    !dRsi(i)/dz(i)
           locRes(j) = 0.0_rk
        end if

        if(s_mode.eq.0) then
           if( VN( globalNM(m,i) ).ne.1 ) then   !base nodes in drop
              do j = LNOPP(i) + Nu, LNOPP(i) + Nv  !u,v,T
                 locJac(j,:) = 0.0_rk
                 locJac(j,j) = 1.0_rk
                 locRes(j) = 0.0_rk
              end do
           end if
           if( substrate.eq.0.0_rk ) then
              j = LNOPP(i) + NT 
              locJac(j,:) = 0.0_rk
              locJac(j,j) = 1.0_rk
              locRes(j) = 0.0_rk
           end if
        end if

     end if

     !pinned contact line      !when s_mode=1, this also applies
     if( BCflagN( globalNM(m,i), 3) .eq. 3 ) then
        j = LNOPP(i) + Nr
        locJac(j,:) = 0.0_rk
        locJac(j,j) = 1.0_rk
        locRes(j) = 0.0_rk
     end if

if(alge_corner.eq.1) then
     !corner algebraic mesh     
     if( algeN(globalNM(m,i)).eq.1 ) then
        j = LNOPP(i) + Nz      !The location of Reta(i) in locRes
        locJac(j,:) = 0.0_rk
        locJac(j,j-1) = k_alge    !dReta/dri
        locJac(j,j) = -1.0_rk    !dReta/dzi
        locRes(j) = k_alge *rcoordinate( globalNM(m,i) ) - zcoordinate( globalNM(m,i) ) - k_alge*x_alge
     end if
     !substrate algebraic mesh
     if( algeS(globalNM(m,i)).eq.2 ) then
        j = LNOPP(i) + Nz      !The location of Reta(i) in locRes
        locJac(j,:) = 0.0_rk
        locJac(j,j-1) = 1.0_rk    !dReta/dri
        locRes(j) = 0.0_rk
     end if
end if
     !substrate algebraic mesh at contact line and vapor region boundary
     if( algeS(globalNM(m,i)).eq.3 .or. algeS(globalNM(m,i)).eq.1 ) then
        j = LNOPP(i) + Nr      !The location of Rsi(i) in locRes
        locJac(j,:) = 0.0_rk
        locJac(j,j) = 1.0_rk    !dRsi/dri
        locRes(j) = 0.0_rk
     end if

     !outer circle
     if( no_vapor.eq.0 .and. (.not.uniflux) .and. &
          BCflagN( globalNM(m,i), 4 ).eq.1 .or. BCflagN( globalNM(m,i), 4 ).eq.3 ) then !vapor outer
        j = LNOPP(i) + Nr      !The location of Rsi(i) in locRes
        locJac(j,:) = 0.0_rk
        locJac(j,j) = 2.0_rk*rcoordinate( globalNM(m,i) )    !dRsi/dri
        locJac(j,j+1) = 2.0_rk*zcoordinate( globalNM(m,i) )    !dRsi/dzi
        locRes(j) = rcoordinate( globalNM(m,i) )**2 + zcoordinate( globalNM(m,i) )**2 - (outer*R)**2

        if(s_mode.eq.0) then
           j = LNOPP(i) + MDF(globalNM(m,i))-1      !The location of Rc(i) in locRes
           locJac(j,:) = 0.0_rk
           locJac(j,j) = 1.0_rk    !dRci/dci
           locRes(j) = 0.0_rk
        end if
     end if
     if( BCflagN( globalNM(m,i), 4 ).eq.2 .or. BCflagN( globalNM(m,i), 4 ).eq.3 ) then    !outer BC for substrate
        j = LNOPP(i) + Nr      !The location of Rsi(i) in locRes
        locJac(j,:) = 0.0_rk
        locJac(j,j) = 1.0_rk    !dRsi/dri
        locRes(j) = 0.0_rk

        ! !if use no flux BC at outer substrate surface, then comment out
        ! if(s_mode.eq.0) then
        !    j = LNOPP(i) + NTs
        !    locJac(j,:) = 0.0_rk
        !    locJac(j,j) = 1.0_rk    !dRt/dT
        !    locRes(j) = 0.0_rk
        ! end if
           

     end if

if(no_vapor.eq.0) then
     !close to free surface
     do n = 1, vlayer-1
        if( layer( globalNM(m,i) ).eq.2*NEV1(n)+1 ) then
           !find equations for the spherical cap
           dRp = R11(n) - R
           zpp = zcoordinate(top_node) + dRp
           d = sqrt( R11(n)**2 + zpp**2 )
           Rpp = d/2.0_rk / sin(angle_temp(n)/2.0_rk)
           ap = ( zpp/R11(n) )**2 + 1.0_rk
           bp = zpp*( R11(n)**2 - zpp**2 ) / R11(n)**2 - 2.0_rk*zpp
           cp = ( (R11(n)**2 - zpp**2 ) / (2.0_rk*R11(n)) )**2 + zpp**2 - Rpp**2
           npp = ( -bp - sqrt( bp**2 - 4.0_rk*ap*cp ) ) / (2.0_rk*ap)
           mpp = ( R11(n)**2 - zpp**2 + 2.0_rk*npp*zpp ) / (2.0_rk*R11(n))

           j = LNOPP(i) + Nr      !The location of Rsi(i) in locRes
           locJac(j,:) = 0.0_rk
           locJac(j,j) = 2.0_rk*( rcoordinate( globalNM(m,i) ) - mpp )   !dRsi/dri
           locJac(j,j+1) = 2.0_rk*( zcoordinate( globalNM(m,i) ) - npp )    !dRsi/dzi
           locRes(j) = ( rcoordinate( globalNM(m,i) ) - mpp )**2 + &
                ( zcoordinate( globalNM(m,i) ) - npp )**2 - Rpp**2
        end if
     end do
end if

     !free   ?can be neglected or not
     if( ( BCflagN( globalNM(m,i), 3 ).eq.1 .or. BCflagN( globalNM(m,i), 3 ).eq.3) ) then
        if(s_mode.eq.0) then
           if( VE(m).eq.1 ) then
              j = LNOPP(i) + MDF(globalNM(m,i))-1      !The location of Rc(i) in locRes
              locJac(j,:) = 0.0_rk
              locJac(j,j) = 1.0_rk    !dRci/dci
              locRes(j) = 0.0_rk
           end if
        else !s_mode=1
           j = LNOPP(i) + Nr      !The location of Rsi(i) in locRes
           locJac(j,:) = 0.0_rk
           locJac(j,j) = 2.0_rk*rcoordinate( globalNM(m,i) )    !dRsi/dri
           locJac(j,j+1) = 2.0_rk*( zcoordinate( globalNM(m,i) ) + zp )    !dRsi/dzi
           locRes(j) = rcoordinate( globalNM(m,i) )**2 + ( zcoordinate( globalNM(m,i) ) + zp )**2 - Rp**2
        end if
     end if

     ! !should be of no use when adding thermal eqns for substrate
     ! if(s_mode.eq.0) then
     !    if(BCflagN(globalNM(m,i),2).ne.0 .and. VN(globalNM(m,i)).eq.1 ) then
     !       j = LNOPP(i) + NT      !The location of Rt(i) in locRes
     !       locJac(j,:) = 0.0_rk
     !       locJac(j,j) = 1.0_rk    !dRti/dTi
     !       locRes(j) = 0.0_rk
     !    end if
     ! end if

     !substrate base
     if( BCflagN( globalNM(m,i), 5 ).ne.0 ) then
        
        if( BCflagN( globalNM(m,i), 5 ).eq.1 .or. BCflagN( globalNM(m,i), 5 ).eq.3 ) then
           !the node shared by Region1 and Region3 ( BCflagN=3 ), should use Reta
           j = LNOPP(i) + Nz      !the location of Reta(i) in locRes
           locJac(j,:) = 0.0_rk
           locJac(j,j) = 1.0_rk      !dReta(i)/dz(i)
           locRes(j) = 0.0_rk  
        else   !end if if( BCflagN( globalNM(m,i), 5 ).eq.2 ) then
           j = LNOPP(i) + Nr
           locJac(j,:) = 0.0_rk
           locJac(j,j+1) = 1.0_rk      !dRsi(i)/dz(i)
           locRes(j) = 0.0_rk
        end if

        if(s_mode.eq.0) then
           j = LNOPP(i) + NTs
           locJac(j,:) = 0.0_rk
           locJac(j,j) = 1.0_rk      !dRt(i)/dT(i)
           locRes(j) = 0.0_rk
        end if

     end if


  end do !for i

  if(initial_vapor_solved.eq.0 .and. s_mode.eq.0) then
     !make variable other than r,z,c not change
     if(VE(m).eq.0) then
        do i = 1, 9    
           do j = LNOPP(i) + NT, LNOPP(i) + Nv     !The location of R(u,v,T)(i) in locRes
              locJac(j,:) = 0.0_rk
              locJac(j,j) = 1.0_rk               !dR(u,v,T)/d(u,v,T)i
              locRes(j) = 0.0_rk
           end do
           if(i.eq.1 .or. i.eq.3 .or. i.eq.7 .or. i.eq.9) then
              j = LNOPP(i) + Np     !The location of Rc(i) in locRes
              locJac(j,:) = 0.0_rk
              locJac(j,j) = 1.0_rk               !dRp/dpi
              locRes(j) = 0.0_rk
           end if
        end do
     end if
     !?necessary to fix from vapor elements?
     if(VE(m).eq.1 .and. BCflagE(m,3).eq.2) then
        do i = 1,7,3
           do j = LNOPP(i) + NT, LNOPP(i) + Nv     !The location of R(u,v,T)(i) in locRes
              locJac(j,:) = 0.0_rk
              locJac(j,j) = 1.0_rk               !dR(u,v,T)/d(u,v,T)i
              locRes(j) = 0.0_rk
           end do
           if(i.eq.1 .or. i.eq.3 .or. i.eq.7 .or. i.eq.9) then
              j = LNOPP(i) + Np     !The location of Rc(i) in locRes
              locJac(j,:) = 0.0_rk
              locJac(j,j) = 1.0_rk               !dRp/dpi
              locRes(j) = 0.0_rk
           end if
        end do
     end if

     !fix T in substrate
     if(VE(m).eq.5) then
        do i = 1, 9
           j = LNOPP(i) + NT
           locJac(j,:) = 0.0_rk
           locJac(j,j) = 1.0_rk               !dRt/dTi
           locRes(j) = 0.0_rk
        end do
     end if

     !fix r,z for all nodes
     do i = 1, 9    
        do j = LNOPP(i) + Nr, LNOPP(i) + Nz     !r,z for all nodes do not need to be solved
           locJac(j,:) = 0.0_rk
           locJac(j,j) = 1.0_rk               !dR(r,z)/d(r,z)i
           locRes(j) = 0.0_rk
        end do
     end do

  end if


  ! !fix nodes in substrate
  ! do i = 1, 9
  !    if(VN(globalNM(m,i)).eq.5) then
  !       do j = LNOPP(i) +Nr, LNOPP(i) + Nz
  !          locJac(j,:) = 0.0_rk
  !          locJac(j,j) = 1.0_rk      !dReta(i)/dz(i)
  !          locRes(j) = 0.0_rk
  !       end do
  !    end if
  ! end do



  if(s_mode.eq.0) then

     if(Capil.eq.0) then
        if( VE(m).eq.0 ) then
           do i = 1, 9   
              if(i.eq.1 .or. i.eq.3 .or. i.eq.7 .or. i.eq.9) then
                 j = LNOPP(i) + Np     !The location of Rc(i) in locRes
                 locJac(j,:) = 0.0_rk
                 locJac(j,j) = 1.0_rk               !dRc/dci
                 locRes(j) = 0.0_rk
              end if
           end do
        end if
     end if
        

     ! !make variable other than r, z, c not change(in case a whole row of Jac for c is 0)
     ! if( VE(m).eq.0 ) then
     !    do i = 1, 9    
     !       do j = LNOPP(i) + Nu, LNOPP(i) + Nv     !The location of Rc(i) in locRes
     !          locJac(j,:) = 0.0_rk
     !          locJac(j,j) = 1.0_rk               !dRc/dci
     !          locRes(j) = 0.0_rk
     !       end do
     !       if(i.eq.1 .or. i.eq.3 .or. i.eq.7 .or. i.eq.9) then
     !          j = LNOPP(i) + Np     !The location of Rc(i) in locRes
     !          locJac(j,:) = 0.0_rk
     !          locJac(j,j) = 1.0_rk               !dRc/dci
     !          locRes(j) = 0.0_rk
     !       end if
     !    end do
     ! end if


     ! !make variable c not change(in case a whole row of Jac for c is 0)
     ! if( VE(m).eq.1 ) then
     !    do i = 1, 9    
     !       j = LNOPP(i) + MDF( globalNM(m,i) ) - 1     !The location of Rc(i) in locRes
     !       locJac(j,:) = 0.0_rk
     !       locJac(j,j) = 1.0_rk               !dRc/dci
     !       locRes(j) = 0.0_rk
     !    end do
     ! end if



     ! !make variable T not change
     ! if( VE(m).eq.0 ) then
     !    do i = 1, 9    
     !       j = LNOPP(i) + NT     !The location of Rt(i) in locRes
     !       locJac(j,:) = 0.0_rk
     !       locJac(j,j) = 1.0_rk               !dRti/dTi
     !       locRes(j) = 0.0_rk
     !    end do
     ! end if
     ! if(VE(m).eq.1 .and. BCflagE(m,3).eq.2) then
     !    do i = 1,7,3
     !       j = LNOPP(i) + NT     !The location of Rt(i) in locRes
     !       locJac(j,:) = 0.0_rk
     !       locJac(j,j) = 1.0_rk               !dRti/dTi
     !       locRes(j) = 0.0_rk
     !    end do
     ! end if

  end if  !s_mode=0

  return

end subroutine Dirichlet_BC

