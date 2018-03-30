
subroutine split_sol
  use kind
  use data!, only: NTN, NNX, NTE, rcoordinate, zcoordinate, usol, vsol, psol, sol, Nr, Nz, Nu, Nv, Np, NOPP, MDF, globalNM, RegN, rowNM, columnNM, NEM, s_mode

  implicit none
  integer(kind=ik):: i, j


  !make u, v, r, z at boundary 0 in sol
  do i = 1, NTN
     
     !axis
     if( BCflagN(i,1).eq.1 ) then
        sol( NOPP(i) + Nr ) = 0.0_rk
        if(s_mode.eq.0) then
           if(VN(i).ne.1 .and. VN(i).ne.5) then
              sol( NOPP(i) + Nu ) = 0.0_rk
           end if
        end if
     end if

     !base
     if( BCflagN(i,2).ne.0 ) then
        sol( NOPP(i) + Nz ) = 0.0_rk    !z
        if(s_mode.eq.0 .and. VN(i).ne.1) then
           sol( NOPP(i) + Nu ) = 0.0_rk
           sol( NOPP(i) + Nv ) = 0.0_rk
           if(substrate.eq.0.0_rk) sol( NOPP(i) + NT ) = 0.0_rk 
        end if
     end if
     
     !the contact line
     if( BCflagN(i,3).eq.3 )  sol( NOPP(i) + Nr ) = R!1.0_rk
     
     !substrate base
     if( BCflagN(i,5).ne.0 ) then
        sol( NOPP(i) + Nz ) = -substrate    !z
        sol( NOPP(i) + NT ) = 0.0_rk 
     end if


     if(s_mode.eq.0) then  !vapor concentration
        
        !free surface    ?can be neglected or not
        if( VN(i).eq.2 ) then
           sol( NOPP(i) + MDF(i)-1 ) = 1.0_rk
        end if

        !outer circle
        if( ( BCflagN(i,4).eq.1 .or. BCflagN(i,4).eq.3 ) .and. no_vapor.eq.0 ) then
           sol( NOPP(i) + MDF(i)-1 ) = Hum
        end if

     end if

  end do


  !extract r, z, u, v, p from sol
  do i = 1, NTN, 1
     rcoordinate(i) = sol( NOPP(i) + Nr ) 
     zcoordinate(i) = sol( NOPP(i) + Nz ) 
     if(s_mode.eq.0) then
        if(VN(i).eq.0 .or. VN(i).eq.2 .or. BCflagN(i,2).eq.1) then
           Tsol(i) = sol( NOPP(i) + NT ) 
           if( VN(i).eq.0 .or. VN(i).eq.2 ) then           
              usol(i) = sol( NOPP(i) + Nu ) 
              vsol(i) = sol( NOPP(i) + Nv ) 
              if ( PN(i).eq.1 )     psol(i) = sol( NOPP(i) + Np ) 
           end if
        end if
        if(VN(i).eq.5) Tsol(i) = sol( NOPP(i) + NTs )
        if(VN(i).eq.1 .or. VN(i).eq.2) csol(i) = sol(NOPP(i) + MDF(i) -1 )
     end if
  end do

  !interpolate p(:)
  if(s_mode.eq.0) then
     do i = 1, NTE
        if(VE(i).eq.0) then
           do j = 2, 8, 6      !j=2 & j=8
              psol( globalNM(i,j) ) = 0.5_rk*( psol( globalNM(i,j-1) ) + psol( globalNM(i,j+1) ) )
           end do
           do j = 4, 6      !j = 4,5,6
              psol( globalNM(i,j) ) = 0.5_rk*( psol( globalNM(i,j-3) ) + psol( globalNM(i,j+3) ) )
           end do
        end if
     end do
  end if
  

  ztop = zcoordinate(top_node)

  return
end subroutine split_sol

