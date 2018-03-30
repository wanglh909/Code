
subroutine initial_condition
  use kind
  use data!, only: NTN, NTE, NEL, NEM, NEV, NNX34, NNX3456, NNR12, globalNM, rNOP, rowNM, columnNM, RegN, Ngrid, rcoordinate, zcoordinate, usol, vsol, psol, sol, R, z, Nr, Nz, Hum, VN
  use NOP_mod, only: NOPP



  implicit none

  integer(kind=ik):: i, ii, j, m, n, NNL, NNM, rowN, rowN_alge, columnN, NEV1expand(vlayer+1), NEM_else
  real(kind=rk):: x, y, z, r1, r2, r3, z1, z2, z3, theta2, theta3, theta1, theta4, &
       Rtemp, k, R11expand(vlayer+1), d, dmin, dmax, Tper
  real(kind=rk):: x1, x2, y2, x3, y3, x4, x5, y5, x6, y6, theta

  integer(kind=ik):: NELp, NEMp, NEM_algep, NEVp
  real(kind=rk):: outerp, substratep

  x = 0.005_rk  !the length of box region
  !z = 0.1_rk   !the length inside algebraic mesh --> 1-x_alge
  !y = 0.5_rk   !the length of the base of region3
  !outer = 20.0_rk  !the length of the base of region 1&2&3, defined in AAdata
  Tper = 0.01_rk

  NNL = 2*NEL+1            !number of nodes of the less side
  NNM = 2*NEM+1            !number of nodes of the more side
  NEM_else = NEM - NEM_alge

  k_alge = 1.0_rk!1.0_rk!10.0_rk
  x_alge = 0.8_rk   !coordinate of algebraic mesh start point
  R11expand(1) = R
  R11expand(2:vlayer) = R11
  R11expand(vlayer+1) = outer*R
  NEV1expand(1) = 0
  NEV1expand(2:vlayer) = NEV1
  NEV1expand(vlayer+1) = NEV
!write(*,*) R11expand

  rowN_alge = 2*NEM_alge + 1    !node row number of nodes with algeN=1

  rcoordinate(:) = 0.0_rk
  zcoordinate(:) = 0.0_rk

  !judge if mesh condition same as 'Sources/elliptic_mesh.dat'
  open(10,file='Sources/elliptic_mesh.dat', Status ='old', action ='read') 
  read(10,'(4i4,2es13.6)') NELp, NEMp, NEM_algep, NEVp, outerp, substratep
  if(NEL.eq.NELp .and. NEM.eq.NEMp .and. NEM_alge.eq.NEM_algep .and. NEV.eq.NEVp .and. &
       outer.eq.outerp .and. substrate.eq.substratep ) then
     print *, 'mesh parameter same'
     read_coordinate_value = 1
     do i = 1, NTN
        read(10,'(2es15.7)') rcoordinate(i), zcoordinate(i)
        ! print *, rcoordinate(i), zcoordinate(i)
     end do
     

else !mesh parameter not same, redefine r&z initial values
  x1 = (1.0_rk-x)*R
  x2 = (1.0_rk-x)*R
  y2 = x*R
  x3 = sqrt(R**2 - (x*R)**2)
  y3 = x*R
  x4 = x_alge*R
  x6 = (k_alge**2*x_alge*R + sqrt(k_alge**4*x_alge**2*R**2 - &
       (1.0_rk+k_alge**2)*(k_alge**2*x_alge**2 - 1.0_rk)*R**2))/(1.0_rk+k_alge**2)
  y6 = sqrt(R**2 - x6**2)
  x5 = (x4+x6)/2.0_rk
  y5 = y6/2.0_rk
  theta1 = atan(x6/y6)
  theta4 = atan(x3/y3)



  do i = 1, NTN, 1            !globalN

     if( RegN( rNOP(i,1,1) ).eq.1 .or. RegN( rNOP(i,2,1) ).eq.1 .or. &
          RegN( rNOP(i,3,1) ).eq.1 .or. RegN( rNOP(i,4,1) ).eq.1 ) then       !nodes in region 1

        rowN = Ngrid(i,1)
        if(no_vapor.eq.0) rowN = rowN - 2*NEV
        columnN = Ngrid(i,2) - 2*NES

        zcoordinate(i) = x*R/real(2*NEL,rk)*real(columnN-1,rk)
        rcoordinate(i) = ( sqrt( R**2 - zcoordinate(i)**2 ) - (1.0_rk - x)*R ) /real(2*NEL,rk) *real(-rowN+NNL,rk) + &
             (1.0_rk - x)*R

     else if( RegN( rNOP(i,1,1) ).eq.2 .or. RegN( rNOP(i,2,1) ).eq.2 .or. &
          RegN( rNOP(i,3,1) ).eq.2 .or. RegN( rNOP(i,4,1) ).eq.2 ) then   !nodes in region 2
        
        rowN = Ngrid(i,1)
        columnN = Ngrid(i,2) - 2*NES
        theta = (pi/2.0_rk-theta4)/real(NNL-1,rk)*real(columnN-1, rk)

        m = 1
        do n = 1, vlayer-1
           if(layer(i).gt.2*NEV1(n)+1)   m = m + 1
        end do
        
        r1 = R11expand(m+1)*cos(theta)
        z1 = R11expand(m+1)*sin(theta)
        r2 = R11expand(m)*cos(theta)
        z2 = R11expand(m)*sin(theta)
        rcoordinate(i) = r2 + (r1-r2) /real(2*(NEV1expand(m+1)-NEV1expand(m)),rk) &
             *real(layer(i)-2*NEV1expand(m)-1,rk)
        zcoordinate(i) = z2 + (z1-z2) /real(2*(NEV1expand(m+1)-NEV1expand(m)),rk) &
             *real(layer(i)-2*NEV1expand(m)-1,rk)


        ! z2 = x*R/real(2*NEL,rk)*real(columnN-1,rk)
        ! r2 = sqrt( R11expand(m)**2 - z2**2 )

        ! ! if(m.eq.vlayer) then
        ! !    z3 = outer*R* cos( pi/2.0_rk* real(-columnN+2*(NEM+NEL)+1,rk) /real(2*(NEM+NEL),rk) )
        ! !    r3 = outer*R* sin( pi/2.0_rk* real(-columnN+2*(NEM+NEL)+1,rk) /real(2*(NEM+NEL),rk) )
        ! ! else
        !    r3 = sqrt( R11expand(m+1)**2 - z2**2 )
        !    z3 = z2
        ! ! end if

        ! zcoordinate(i) = z2 + (z3-z2) &
        !      /real(2*(NEV1expand(m+1)-NEV1expand(m)),rk) *real(layer(i)-2*NEV1expand(m)-1,rk) 
        ! rcoordinate(i) = r2 + (r3-r2) &
        !      /real(2*(NEV1expand(m+1)-NEV1expand(m)),rk) *real(layer(i)-2*NEV1expand(m)-1,rk) 



     else if( VN(i).ne.5 ) then! .and. BCflagN(i,2).ne.0 ) then         !nodes in region 3&4&5

        !calculate rowN&columnN of nodes in region 3&4&5
        ii = i-NNR1278
        if( mod( ii,NNX3456 ) .eq. 0 ) then
           rowN = ii/NNX3456         !no dble, rowN of nodes not elements (different from rowN in NOP)
           columnN = NNX3456
        else
           rowN = ii/NNX3456 + 1         !no dble
           columnN = mod( ii,NNX3456 )
        end if
        columnN = columnN - 2*NES
        
        ! if(columnN.le.NNX34) then   !nodes in region 3&4
        
        if(rowN.ge.rowN_alge) then !nodes outside algebraic corner in region 3&4&5
           rowN = rowN - 2*NEM_alge
           theta = theta1/real(2*NEM_else,rk)*real(2*NEM_else - rowN + 1, rk)
           if(columnN.le.NNX34) then   !nodes in region 3&4
              r1 = R*sin(theta)
              z1 = R*cos(theta)
              r2 = x4/real(2*NEM_else,rk)*real(2*NEM_else - rowN + 1, rk)
              rcoordinate(i) = (r1-r2)/real(NNX34-1,rk)*real(columnN-1,rk) + r2
              zcoordinate(i) = z1/real(NNX34-1,rk)*real(columnN-1,rk)
           ! else if(columnN.le.NNX34) then   !nodes in region 4


           else   !nodes in region 5
              
              m = 1
              do n = 1, vlayer-1
                 if(layer(i).gt.2*NEV1(n)+1)   m = m + 1

                 !write(*,*) layer(i), NEV1
              end do

              r1 = R11expand(m+1)*sin(theta)
              z1 = R11expand(m+1)*cos(theta)
              r2 = R11expand(m)*sin(theta)
              z2 = R11expand(m)*cos(theta)
              rcoordinate(i) = r2 + (r1-r2) /real(2*(NEV1expand(m+1)-NEV1expand(m)),rk) &
                   *real(layer(i)-2*NEV1expand(m)-1,rk)
              zcoordinate(i) = z2 + (z1-z2) /real(2*(NEV1expand(m+1)-NEV1expand(m)),rk) &
                   *real(layer(i)-2*NEV1expand(m)-1,rk)



           end if
        else     !nodes within algebraic corner
           if(columnN.le.NNL) then   !nodes in region 3 inside algebraic corner
              r1 = (x2-x5)/real(2*NEM_alge,rk)*real(2*NEM_alge - rowN + 1, rk) + x5
              z1 = (y2-y5)/real(2*NEM_alge,rk)*real(2*NEM_alge - rowN + 1, rk) + y5
              r2 = (x1-x4)/real(2*NEM_alge,rk)*real(2*NEM_alge - rowN + 1, rk) + x4
              rcoordinate(i) = (r1-r2)/real(2*NEL,rk)*real(columnN-1,rk) + r2
              zcoordinate(i) = z1/real(2*NEL,rk)*real(columnN-1,rk)
           else   !nodes in region 4&5 inside algebraic corner 
              columnN = columnN - NNL + 1
              theta = (theta4-theta1)/real(2*NEM_alge,rk)*real(2*NEM_alge - rowN + 1, rk) + theta1

              if(columnN.le.NNX34-NNL+1) then   !nodes in region 4
                 r1 = R*sin(theta)
                 z1 = R*cos(theta)
                 r2 = (x2-x5)/real(2*NEM_alge,rk)*real(2*NEM_alge - rowN + 1, rk) + x5
                 z2 = (y2-y5)/real(2*NEM_alge,rk)*real(2*NEM_alge - rowN + 1, rk) + y5
                 rcoordinate(i) = (r1-r2)/real(2*NEL,rk)*real(columnN-1,rk) + r2
                 zcoordinate(i) = (z1-z2)/real(2*NEL,rk)*real(columnN-1,rk) + z2
              else   !nodes in region 5
                 m = 1
                 do n = 1, vlayer-1
                    if(layer(i).gt.2*NEV1(n)+1)   m = m + 1

                    !write(*,*) layer(i), NEV1
                 end do

                 r1 = R11expand(m+1)*sin(theta)
                 z1 = R11expand(m+1)*cos(theta)
                 r2 = R11expand(m)*sin(theta)
                 z2 = R11expand(m)*cos(theta)
                 rcoordinate(i) = r2 + (r1-r2) /real(2*(NEV1expand(m+1)-NEV1expand(m)),rk) &
                      *real(layer(i)-2*NEV1expand(m)-1,rk)
                 zcoordinate(i) = z2 + (z1-z2) /real(2*(NEV1expand(m+1)-NEV1expand(m)),rk) &
                      *real(layer(i)-2*NEV1expand(m)-1,rk)

              end if
           end if


        ! else                       !nodes in region 5

        !    m = 1
        !    do n = 1, vlayer-1
        !       if(layer(i).gt.2*NEV1(n)+1)   m = m + 1

        !       !write(*,*) layer(i), NEV1
        !    end do

        !    theta2 = atan( sqrt( R11expand(m)**2 - x**2 ) / x )
        !    z2 = R11expand(m)* cos( theta2* real(-rowN+NNM,rk) /real(2*NEM,rk) )
        !    r2 = R11expand(m)* sin( theta2* real(-rowN+NNM,rk) /real(2*NEM,rk) )

           
        !    if(m.eq.vlayer) then
        !       theta3 = pi/2.0_rk / real(NEM+NEL,rk) * real(NEM,rk)
        !    else
        !       theta3 = atan( sqrt( R11expand(m+1)**2 - x**2 ) / x )
        !    end if

        !    z3 = R11expand(m+1)* cos( theta3* real(-rowN+NNM,rk) /real(2*NEM,rk) )
        !    r3 = R11expand(m+1)* sin( theta3* real(-rowN+NNM,rk) /real(2*NEM,rk) )

        !    zcoordinate(i) = z2+ (z3-z2) &
        !         /real(2*(NEV1expand(m+1)-NEV1expand(m)),rk) *real(layer(i)-2*NEV1expand(m)-1,rk)
        !    rcoordinate(i) = r2+ (r3-r2) &
        !         /real(2*(NEV1expand(m+1)-NEV1expand(m)),rk) *real(layer(i)-2*NEV1expand(m)-1,rk) 

        end if

     else  !nodes in region 6&7&8

        rowN = 2*rowNM( rNOP(i,1,1) )
        columnN = 2*columnNM( rNOP(i,1,1) )

        if( RegN( rNOP(i,1,1) ).eq.6 ) then
           rowN = rowN - ( rNOP(i,1,2) - 1 )/3 + 1
           if(mod(rNOP(i,1,2),3).eq.0) then
              columnN = columnN + 1                 
           else
              columnN = columnN + mod(rNOP(i,1,2),3) - 2
           end if

           if( rowN.le.2*(NEM-NEM_alge) ) then
              rcoordinate(i) = x_alge*R /real( 2*(NEM-NEM_alge),rk ) *real( rowN-1,rk )
           else 
              rcoordinate(i) = x_alge*R + &
                   (1.0_rk-x-x_alge)*R /real( 2*NEM_alge,rk ) *real( rowN-2*(NEM-NEM_alge)-1,rk )
           end if

        else   !nodes in region 7&8
           columnN = columnN + ( rNOP(i,1,2) - 1 )/3 - 1
           if(mod(rNOP(i,1,2),3).eq.0) then
              rowN = rowN + 1                 
           else
              rowN = rowN + mod(rNOP(i,1,2),3) - 2
           end if
           if( RegN( rNOP(i,1,1) ).eq.7 ) then
              rowN = rowN - 2*NEM
              rcoordinate(i) = (1.0_rk-x)*R + x*R /real( 2*NEL,rk ) * real( rowN-1,rk )
           else  !in region8
              rowN = rowN - 2*(NEM+NEL)

              if(no_vapor.eq.0) then
                 m = 1
                 do n = 1, vlayer-1
                    if(rowN.gt.2*NEV1(n)+1)   m = m + 1
                 end do
                 rcoordinate(i) = R11expand(m) + ( R11expand(m+1)-R11expand(m) )&
                      /real(2*(NEV1expand(m+1)-NEV1expand(m)),rk) *real(rowN-2*NEV1expand(m)-1,rk)

                 ! if( rowN.le.2*(NEM+NEL+NEV1(1)) ) then
                 !    rcoordinate(i) = R + (R11(1)-R) /real(2*NEV1(1),rk) *real(rowN-2*(NEL+NEM)-1,rk))
                 ! else if( rowN.le.2*(NEM+NEL+NEV1(2)) ) then
                 !    rcoordinate(i) = R11(1) + (R11(2)-R11(1)) &
                 !         /real(2*(NEV1(2)-NEV1(1)),rk) *real(rowN-2*(NEL+NEM+NEV1(1))-1,rk))
                 ! else
                 !    rcoordinate(i) = R11(2) + (outer*R-R11(2))
                 ! end if

              else   !no_vapor
                 rcoordinate(i) = R + (outer-R)/real(2*NEV,rk)*real(rowN-1,rk)
              end if

           end if
        end if
        zcoordinate(i) = -substrate/real(2*NES,rk) * real(2*NES-columnN+1,rk)



     end if  !regions

     if(VN(i).eq.1 .and. no_vapor.eq.1) then
        zcoordinate(i) = 0.0_rk
        rcoordinate(i) = R + (outer-R)/real(2*NEV,rk) *real(2*NEV-Ngrid(i,1)+1,rk)   
     end if

  end do  !nodes

end if  !judge if use data value or guess value for r&z
close(10)
              
  do i = 1, NTN
     if( VN(i).eq.0 .or. VN(i).eq.2 ) then
        usol(i) = 0.0_rk
        vsol(i) = 0.0_rk
        Tsol(i) = 0.0_rk!?
        psol(i) = 2.0_rk!/Ca    !include p at the node and p needed to be interpolated
     end if

     if( VN(i).eq.1 .or. VN(i).eq.2 ) then
        k = (1.0_rk - Hum) / real(NNXV-1,rk)
        csol(i) = -k*real(layer(i),rk) + 1.0_rk + k  !csat=1 is saturated concentration

        ! if( layer(i).eq.1 ) then
        !    dmin = sqrt(rcoordinate(i)**2+zcoordinate(i)**2)
        !    dmax = sqrt(rcoordinate(i+2*NEV)**2+zcoordinate(i+2*NEV)**2)
        !    csol(i) = 1.0_rk  !csat=1 is saturated concentration
        ! else
        !    d = sqrt(rcoordinate(i)**2+zcoordinate(i)**2)
        !    csol(i) = (Hum-1.0_rk)/real(2*NEV,rk)*real(layer(i)-1)+1.0_rk
        !    ! csol(i) = (Hum-1.0_rk)/(dmax-dmin)*(d-dmin)+1.0_rk
        ! end if
     end if
     
     if(VN(i).eq.5) then
        Tsol(i) = 0.0_rk
        if( BCflagN(i,5).ne.0 )  Tsol(i) = T_sub
     end if
     
        

  end do
  
  
  do i = 1, NTN, 1
     sol( NOPP(i) + Nr ) = rcoordinate(i)
     sol( NOPP(i) + Nz ) = zcoordinate(i)
     if( VN(i).eq.0 .or. VN(i).eq.2 ) then
        sol( NOPP(i) + Nu ) = usol(i)
        sol( NOPP(i) + Nv ) = vsol(i)
        sol( NOPP(i) + NT ) = Tsol(i)
        if ( PN(i).eq.1 )  then
           sol( NOPP(i) + Np ) = psol(i)
        end if
     end if
     if( VN(i).eq.1 .or. VN(i).eq.2 ) then
        sol( NOPP(i) + MDF(i)-1 ) = csol(i)
     end if
     if( VN(i).eq.1 .and. BCflagN(i,2).eq.1 ) sol( NOPP(i) + NT) = Tsol(i)
     if( VN(i).eq.5 )  sol( NOPP(i) + NTs) = Tsol(i)
  end do


! open(unit = 10, file = trim(folder)//'initial.dat', status = 'replace')
! do i = 1, NTN
!    if(VN(i).eq.0) then
!       write(10, '(i8, 5es15.7)', advance = 'no' ) i, sol( NOPP(i) + Nr ), sol( NOPP(i) + Nz ), &
!            sol( NOPP(i) + Nu ), sol( NOPP(i) + Nv ), sol( NOPP(i) + NT )
!       if(PN(i).eq.1) then 
!          write(10, '(es15.7)' ) sol( NOPP(i) + Np )
!       else
!          write(10, '(A)' ) ' '
!       end if
!    end if

!    if(VN(i).eq.1)   write(10, '(i8,2es15.7,A,A,A,es15.7)' ) i, sol( NOPP(i) + Nr ), sol( NOPP(i) + Nz ), &
!         '               ', '               ', '               ', '               ', sol( NOPP(i) + MDF(i)-1 )

!    if(VN(i).eq.2) then
!       write(10, '(i8, 5es15.7)', advance = 'no' ) i, sol( NOPP(i) + Nr ), sol( NOPP(i) + Nz ), &
!         sol( NOPP(i) + Nu ), sol( NOPP(i) + Nv ), sol( NOPP(i) + NT )
!       if(PN(i).eq.1) then
!          write(10,'(es15.7)', advance = 'no' ) sol( NOPP(i) + Np )
!       else
!          write(10,'(A)', advance = 'no' ) '               '
!       end if
!       write(10, '(es15.7)')  sol( NOPP(i) + MDF(i)-1 )
!    end if
! end do
! close(10)


  return
end subroutine initial_condition
