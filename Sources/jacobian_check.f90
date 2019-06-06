

subroutine jacobian_check(m, i, sj, sf, LNOPP, LNVar, id)
  use kind
  use data

  implicit none

  integer(kind=ik), intent(in):: m, i, LNVar, LNOPP(LNVar), id
  real(kind=rk), intent(in):: sj(LNVar, LNVar), sf(LNVar)

  integer(kind=ik):: k,p,q
  real(kind=rk):: sj_check, sf_check(LNVar), diffsj, rel, delta, TOL
  sf_check = 0.0_rk
  !write(*,*) 'secand method check for element m =',m, 'point i =', i

  delta = 10e-7_rk
  TOL = 10e-4_rk
 
  if(i.eq.1) then
     open(unit = 20, file = trim(folder)//'secand_method.dat', status = 'replace')
     write(20,'(A,i4)') 'element to check:', m
  else
     open(unit = 20, file = trim(folder)//'secand_method.dat', status = 'old', access = 'append')
  end if

  !change sol( NOPP( globalNM(m,k) ) + p )
  do k = 1, bas
     do p = 0, MDF( globalNM(m,k) )-1    !p = 0 --> si, p = MDF( globalNM(m,k) ) --> v or p

        !add delta to xj
        sol( NOPP( globalNM(m,k) ) + p ) = sol( NOPP( globalNM(m,k) ) + p ) + delta       
        !define soldot
        if (timestep.le.5) then
           soldot = ( sol - solp )/dt
        else
           soldot = 2.0_rk*( sol - solp )/dt - soldotp
        end if

        call values_in_an_element(m,id)

        call define_sf(m,i,sf_check, LNVar, LNOPP, id)
        
        !"do i = 1, bas" is in assemble subroutine
        do q = 2, MDF( globalNM(m,i) )-1    !q = 0 --> Rsi, q = MDF( globalNM(m,i) ) --> Rv or Rp
           
           ! if(q.eq.0) then
           !    if(i.ne.3 .and. i.ne.6 .and. i.ne.9) cycle
           ! end if
           ! if(q.eq.1) cycle

           sj_check = ( sf_check( LNOPP(i)+q ) - sf( LNOPP(i)+q ) )/delta
           diffsj = abs( sj_check - sj( LNOPP(i)+q, LNOPP(k)+p ) ) 
           rel = abs( diffsj / sj( LNOPP(i)+q, LNOPP(k)+p ) )

           if( rel.gt.TOL ) then   ! .and. sj( LNOPP(i)+q , LNOPP(k)+p ).gt.1.0e-8_rk
              if( diffsj.gt.TOL) then
                 ! if(rel.gt.0.1_rk) then  !write(*,*) '!!!'
                    write(*,*) 'ele =', m, ', i =', i,q+1, ', j =', k,p+1
                    ! write(*,*) 'diff =', diffsj, ', rel =', rel
                    write(*,*) 'sf =', sf( LNOPP(i)+q ), ', sf_check =', sf_check( LNOPP(i)+q )
                    write(*,*) 'sj =', sj( LNOPP(i)+q , LNOPP(k)+p ), ', sj_check =', sj_check
                 ! end if
                !  pause
              ! else if( sj( LNOPP(i)+q , LNOPP(k)+p ).gt.1.0e-7_rk ) then
              !    write(*,*) 'relative error large, but diffsj small'
              !    write(*,*) 'ele =', m, ', i =', i,q+1, ', j =', k,p+1
              end if
           ! else
           !    write(*,*) 'ele =', m, ', i =', i,q+1, ', j =', k,p+1
           end if

           write(20,'(A,i4,A,2i4,A,2i4)') 'ele =', m, ', i =', i,q+1, ', j =', k,p+1
           write(20,'(A,es13.6,A,es13.6)') 'sf =', sf( LNOPP(i)+q ), ', sf_check =', sf_check( LNOPP(i)+q )
           write(20,'(A,es13.6,A,es13.6)') 'sj =', sj( LNOPP(i)+q , LNOPP(k)+p ), ', sj_check =', sj_check
           write(20,'(A)') ' '

        end do
        

        !substract deltax to uj
        sol( NOPP( globalNM(m,k) ) + p ) = sol( NOPP( globalNM(m,k) ) + p ) - delta
        !define soldot
        if (timestep.le.5) then
           soldot = ( sol - solp )/dt
        else
           soldot = 2.0_rk*( sol - solp )/dt - soldotp
        end if
        call values_in_an_element(m,id)

     end do
  end do

  close(20)

  return

end subroutine jacobian_check
