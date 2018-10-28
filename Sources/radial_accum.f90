subroutine radial_accumulation
  use kind
  use basis_f
  use data, only: NEL, NTE, folder, globalNM, rcoordinate, zcoordinate, cpsol, radial_cal_time, angle_c, angle_c_degree

  integer(kind=ik)::divide, k, eleN, j, jp
  real(kind=rk):: rlocation, dr1, dr2, dr3, dr4, rcut, z_surf, zlocation, dz(2*NEL), zaccum, &
       x(5), y(5), crossP(5), judge(4), jj(2,3), ff(2), delta(2), &
       sieta(2), rnode(9), znode(9), rcheck, zcheck, a, cp(2*NEL+1)
  real(kind=rk), parameter:: TOL = 1.0e-6_rk

  !print *, radial_cal_time
  !open data file
  if(radial_cal_time.eq.1) then
     open(unit=50, file = trim(folder)//'radial_cp.dat', status = 'replace')
  else
     open(unit=50, file = trim(folder)//'radial_cp.dat', status = 'old', access = 'append')
  end if
  write(50,'(A)') 'variables = "r", "cp" '
  write(50, '(A,f6.3,A)') 'Zone T = "theta=', angle_c_degree, '"'

  cp = 0.0_rk
  
  divide = 2*NEL
  rcut = 0.99_rk

  rlocation = 0.0_rk
  dr1 = 0.02_rk
  dr2 = 0.002_rk
  dr3 = 0.0002_rk
  dr4 = 0.00005_rk
  do while (rlocation.lt.1.0_rk-0.00005_rk)
     z_surf = sqrt( (1.0_rk/sin(angle_c))**2 - rlocation**2 ) - 1.0_rk/tan(angle_c)
     dzfix = z_surf/real(divide,rk)
     dz = dzfix
     dz(1) = 0.3_rk*dzfix
     dz(2) = 0.3_rk*dzfix
     dz(3) = 0.4_rk*dzfix
     dzrest = ( z_surf - 2.0_rk*( dz(1)+dz(2)+dz(3) ) ) / real(divide-6,rk)
     dz(4:2*NEL-3) = dzrest
     dz(2*NEL-2) = dz(3)
     dz(2*NEL-1) = dz(2)
     dz(2*NEL) = dz(1)*0.9_rk

     zaccum = 0.0_rk
     do k = 1, divide
        zaccum = zaccum + dz(k)
        ! print *, k, zaccum
        ! pause
     end do
     
     cp_radial = 0.0_rk
     zlocation = 0.0_rk
     !do or not do +1, zlocation+dz condition is different, change it at the end of this do loop
     !not+1: use cp_lower*dz to calculate accumulation
     !+1: use (cp_lower+cp_upper)/2*dz, but met error near the CL
     do k = 1, divide+1
        if( rlocation.ge.rcut .and. k.eq.divide+1 ) cycle
        
        ! print *, '(r,z)', rlocation,zlocation
        !(r,z) given

        if(zlocatin.gt.z_surf+TOL) then
           print *, 'z larger than z_surf, z', zlocation, 'z_surf', z_surf, 'divideNumk', k
           pause
        else if( abs(zlocation-z_surf).lt.TOL ) then
           zlocation = z_surf
        end if

        !judge which element is it in
        do eleN = 1, NTE
           do j = 1, 4
              if(j.eq.1) then
                 jp = 1
              else if(j.eq.2) then
                 jp = 3
              else if(j.eq.3) then
                 jp = 9
              else if(j.eq.4) then
                 jp = 7
              else
                 print *, 'error in j(which element, radial_accum)'
              end if
              x(j) = rcoordinate(globalNM(eleN,jp))
              y(j) = zcoordinate(globalNM(eleN,jp))
           end do
           x(5) = x(1)
           y(5) = y(1)
           
           do j = 1, 4
              crossP(j) = ( x(j+1)-x(j) )*( zlocation - y(j) ) - ( rlocation - x(j) )*( y(j+1)-y(j) )
           end do
           crossP(5) = crossP(1)
           
           flag = 1
           do j = 1, 4
              judge(j) = crossP(j)*crossP(j+1)
              if(.not.judge(j).ge.0.0_rk) flag = 0
           end do

           ! if(rlocation.eq.0.02_rk .and. k.eq.17 .and. (eleN.eq.NTE-divide) ) then
           !    print *, 'r', rlocation, x(1:4)
           !    print *, 'z', zlocation, y(1:4)
           !    print *, 'crossP', crossP
           !    print *, 'juedge', judge
           !    print *, 'element', eleN
           !    pause
           ! end if
           
           if(flag.eq.1) then
              ! print *, 'r', rlocation, x(1:4)
              ! print *, 'z', zlocation, y(1:4)
              ! print *, 'crossP', crossP
              ! print *, 'judge', judge
              ! print *, 'judge to be in this element', eleN
              !pause
              exit
           else if(eleN.eq.NTE) then
              print *, 'did not find an element for this r,z', rlocation, zlocation, 'zdivideNum k', k
              pause
           end if
           
        end do

        !interpolate with in element eleN

        !get r1-9, z1-9
        do j = 1, 9
           rnode(j) = rcoordinate( globalNM(eleN,j) )
           znode(j) = zcoordinate( globalNM(eleN,j) )
        end do

        !get sieta(1:2), viz.(si,eta), with Newton's method
        sieta(:) = 0.5_rk  !initial guess
        do
           jj(:,:) = 0.0_rk
           ff(:) = -(/rlocation,zlocation/)
           delta(:) = 0.0_rk
           do j = 1, 9
              jj(1,1) = jj(1,1) - rnode(j)*phiisi( sieta(1), sieta(2), j )
              jj(1,2) = jj(1,2) - rnode(j)*phiieta( sieta(1), sieta(2), j )
              jj(2,1) = jj(2,1) - znode(j)*phiisi( sieta(1), sieta(2), j )
              jj(2,2) = jj(2,2) - znode(j)*phiieta( sieta(1), sieta(2), j )
              ff(1) = ff(1) + rnode(j)*phii( sieta(1), sieta(2), j )
              ff(2) = ff(2) + znode(j)*phii( sieta(1), sieta(2), j )
           end do
           call FullGaussSolver(jj,ff,2,delta)
           sieta = sieta + delta
           ! print *, 'si, eta', sieta

           
           if( sqrt(delta(1)**2+delta(2)**2).le.TOL .and. sqrt(ff(1)**2+ff(2)**2).le.TOL ) then
              ! if( sieta(1).ge.-0.1_rk .and. sieta(1).le.1.0_rk+0.1_rk .and. &
              !      sieta(2).ge.-0.1_rk .and. sieta(2).le.1.0_rk+0.1_rk ) then  !si,eta are between 0 and 1
                 exit
              ! else
              !    print *, 'si, eta calculated to wrong value by Newtons method', sieta
              !    pause
              !    sieta(:) = 0.1_rk
              ! end if
           end if
        end do
        
        !check if rcoordinate(si,eta), zcoordinate(si,eta) are right
        rcheck = 0.0_rk
        zcheck = 0.0_rk
        do j = 1, 9
           rcheck = rcheck + rnode(j)*phii( sieta(1), sieta(2), j )
           zcheck = zcheck + znode(j)*phii( sieta(1), sieta(2), j )
        end do
        a = abs(rlocation-rcheck)
        if( abs(rlocation-rcheck).gt.abs(zlocation-zcheck) ) a = abs(zlocation-zcheck)
        if(a.gt.TOL) then
           print *, 'wrong si,eta calculated. (r,z)=', rlocation, zlocation
           stop
        end if

        !calculate cp(r,z) --> cp(k)
        cp(k) = 0.0_rk
        do j = 1, 9
           cp(k) = cp(k) + cpsol(globalNM(eleN,j))*phii( sieta(1), sieta(2), j )
        end do
        
        if(k.ne.divide+1) zlocation = zlocation + dz(k)
     end do  !for k, z value loop

     cp_radial = 0.0_rk
     do k = 1, divide
        if(rlocation.lt.rcut) then
           cp_radial = cp_radial + ( cp(k) + cp(k+1) )/2.0_rk *dz(k)
        else
           cp_radial = cp_radial + cp(k) *dz(k)
        end if
     end do
        

     ! print *, rlocation, cp_radial
     !write in file
     write(50,'(2es15.7)')  rlocation, cp_radial

     if(rlocation.le.0.7_rk) then
        rlocation = rlocation + dr1
     else if(rlocation.le.0.99_rk) then
        rlocation = rlocation + dr2
     else if(rlocation.le.0.995_rk) then
        rlocation = rlocation + dr3
     else
        rlocation = rlocation + dr4        
     end if
        
  end do  !loop of rlocation

  close(50)
  
  return

end subroutine radial_accumulation 




	Subroutine FullGaussSolver(JJ,FF,Nodes,SolutionMat)
!notice: J should be allocated as J(n,n+1)
!J gets changed after this program 

! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC														C
! C	FullGaussSolve is program used to solve full		C
! C	a matrix equation of the form:						C
! C														C
! C					JJ*u=FF								C
! C														C
! C	where JJ and FF are the Jacobian and RHS matrix		C
! C	sent into the program with dimensions (n,n) and		C
! C	n, respectively.  The program solves the matrices	C
! C	by full gaussian elimination and the resulting		C
! C	solution is stored and passed out as SolutionMat	C
! C														C
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	Implicit Integer (i-n)
	Implicit Real*8 (a-h,o-z)

	Real*8 rowhold(1,Nodes)
	Real*8 JJ(Nodes,Nodes+1), FF(Nodes)
	Real*8, INTENT(OUT):: SolutionMat(Nodes)

	!Append the response matrix to jacobian
	Do i = 1,Nodes
	  JJ(i,Nodes+1) = FF(i)
	End Do
	 
	Do j = 1,Nodes       

	  !Find the largest entry in the row, the pivot
	  entrymax = 0.0D0
	  Do i = j,Nodes
	    If (ABS(JJ(i,j)).gt.entrymax) then
	      entrymax = ABS(JJ(i,j))
		  imaxpos = i
	    End If
	  End Do

	  If (imaxpos.ne.j) then
	    Do m = 1,Nodes+1
		  c = JJ(j,m)
	
            JJ(j,m) = JJ(imaxpos,m)
            JJ(imaxpos,m) = c
	
          End Do
	 
	  End If



	  !Divide the pivot row by the current entry value
	  c = JJ(j,j) 
	  Do m = j,Nodes+1
	    JJ(j,m) = JJ(j,m)/c
	  End Do

	  !Do Elimination
	  Do i = j+1,Nodes
	    c = JJ(i,j)
		Do m = j,Nodes+1
	      JJ(i,m) = JJ(i,m) - c*JJ(j,m)
	    End Do
	  End Do

	End Do
	
	!Perform Back Substitution
	Do j = Nodes,1,-1
	  Do i = (j-1),1,-1
	    c = JJ(i,j)
	    Do m = j,Nodes+1
	      JJ(i,m) = JJ(i,m)-c*JJ(j,m)
	    End Do
	   End Do
	End Do

	Do i = 1,Nodes
	  SolutionMat(i) = JJ(i,Nodes+1)
	End Do

	Return
	End Subroutine FullGaussSolver
