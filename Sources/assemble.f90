
  !call assemble_local(ele,local,loc,NOPPl,NB,id,dum)
  !ele: element to be assembled (DONT CHANGE)
  !local(NB,NB): local matrix (Warning assemble as (j,i), we normally assemble i,j as i = row, j = col)
  !NB: Integer, defines sizes (# number of variables in element (DONT CHANGE)
  !loc(NB): local rhs
  !NOPPl(bas): Local NOPPl, must be filled and returned
  !id: Thread id, integer (DONT CHANGE)
  !dum: integer no purpose (DONT CHANGE)

subroutine assemble_local(m, locJacTr, locRHS, LNOPP, LNVar, id, dum)

  use kind
  use data

  implicit none

  integer(kind=ik), intent(in):: m, LNVar, id, dum
  real(kind=rk), intent(out):: locJacTr(LNVar, LNVar), locRHS(LNVar)

  integer(kind=ik):: LNOPP(bas)
  real(kind=rk):: locJac(LNVar, LNVar), locRes(LNVar)

  integer(kind=ik):: i,j, k,l, p,q
  integer(kind=ik):: jac_check, LNVart  !true local number of variables

  jac_check = 0
  if( m.eq.0 .and. initial_vapor_solved.eq.1) then  ! m.eq.28 .and. 
     jac_check = 1  
  end if
  !the element to check, set m to 0 if not check

  LNOPP = 0
  do i = 1, bas
     do j = 1, i-1
        LNOPP(i) = LNOPP(i) + MDF( globalNM(m,j) )
     end do
     LNOPP(i) = LNOPP(i) + 1
  end do
  LNVart = LNOPP(bas) + MDF(globalNM(m,bas)) - 1

  locJac = 0.0_rk
  locJacTr = 0.0_rk
  locRes = 0.0_rk
  locRHS = 0.0_rk  
  
  call values_in_an_element(m,id)

  do i = 1, bas, 1            !i is the i in the residual Res(i) & Jac(ij)

     !define locRes: R(r,z,c)i --> Res
     call define_sf(m,i, locRes, LNVar, LNOPP,id)

     !define locJac: dR(r,z)i/d(u,v,p)j --> Jac
     do j = 1, bas, 1            !stand for (u,v,p)j, j is the j in the notation of uj, vj, pj
        call values_in_sj(m,i,j,id)
        call VI_in_sj(m,i,j, locJac, LNVar, LNOPP,id)
        call SI_in_sj(m,i,j, locJac, LNVar, LNOPP,id)
        call reverse_sj(m,i,j, locJac, LNVar, LNOPP,id)
     end do         !end loop for j

     if( jac_check.eq.1 )  call jacobian_check(m, i, locJac, locRes, LNOPP, LNVar, id)
     ! if( s_mode.eq.0)  call jacobian_check(m, i, locJac, locRes, LNOPP, LNVar, id)

  end do             !end loop for i, Res(i) & Jac(i,j)

  call Dirichlet_BC(m, locJac, locRes, LNVar, LNOPP)
  
  do i = 1, LNVar
     do j = 1, LNVar
        locJacTr(i,j) = locJac(j,i)
     end do
  end do

  locRHS = -locRes  !locRes



!------------------------------------debug------------------------------------------------

!for debug, pause for jacobian check
 if( jac_check.eq.1 ) pause
  ! if(  m.eq.ele_c ) pause

!for debug, put 'sj's together as Jac
if( check_0_in_Jac.eq.1 ) then!.and. initial_vapor_solved.eq.1) then
   do i = 1, LNVart
      do k = 1, bas
         if( i.lt.LNOPP(k) ) exit
      end do
      k=k-1

      do j = 1, LNVart
         do l = 1, bas
            if( j.lt.LNOPP(l) ) exit
         end do
         l=l-1

         Jac(NOPP( globalNM(m,k) ) + i-LNOPP(k), NOPP( globalNM(m,l) ) + j-LNOPP(l) ) = &
            Jac(NOPP( globalNM(m,k) ) + i-LNOPP(k), NOPP( globalNM(m,l) ) + j-LNOPP(l) ) + locJac(i,j)

      end do

     Res( NOPP( globalNM(m,k) ) + i-LNOPP(k) ) = Res( NOPP( globalNM(m,k) ) + i-LNOPP(k) ) + locRes(i)

   end do
end if





  return
end subroutine assemble_local




