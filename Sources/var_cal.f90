
subroutine variable_cal
  use kind
  use omp_lib
  use data
  use Ldata, only: dcdsi, dcdeta, rsi_right, zsi_right, reta_right, zeta_right, rintfac_right, rlocal, zlocal, cplocal, gammalocal
  use basis_f!, only: phii_1d, phiix_1d
  use NOP_mod, only: gaussian_quadrature, gaussian_quadrature_1d
  use special_points
  implicit none


  integer(kind=ik):: i,j,k,n1, n2, n3
  real(kind=rk):: rsi, reta, zsi, zeta, csi, ceta, rsi1, reta1, zsi1, zeta1, csi1, ceta1, rdot, zdot
  real(kind=rk), allocatable:: flux(:), flux1(:), Dflux(:), flux_gp(:)
  real(kind=rk):: J0, volume1, volume2
  ! real(kind=rk):: a, b, c, eta0, eta01, eta02, r_change
  real(kind=rk):: retap(3), zetap(3)
  real(kind=rk):: Rp, angle_c_sphe, err_sphe, z_sphe
  real(kind=rk):: MaranD, gradP, peta, Teta(3), gradT(3)
  real(kind=rk):: p0, Pe_change
  real(kind=rk):: t, cut = 1.0e-4_rk, delta_eta = 5e-2_rk
  real(kind=rk):: v_surf_p(3), h_surf, dPdr(3)
  
  real(kind=rk):: particle_m, intMass(3,3), intVol(3,3), cpintfac, rintfac, Jp
  real(kind=rk):: particle_sub, particle_free, intsurfMass(3), gammaintfac
  integer(kind=ik):: l, n, jpp
  real(kind=rk), allocatable:: cpsolp(:)

  t = REAL(omp_get_wtime(),rk)

  !----------------------------------contact angle------------------------------------
  if(final_size.eq.1) then
  angle_c = atan( zcoordinate(angle_c_node)/ ( rcoordinate(CL_node) - rcoordinate(angle_c_node) ) )
  !atan( solp( NOPP(angle_c_node)+Nz ) / ( solp( NOPP(1)+Nr ) - solp( NOPP(angle_c_node)+Nr ) ) )
  angle_c_degree = angle_c /pi*180.0_rk   !degree
  end if
  write(*,*) ' '
  write(*,*) 'contact angle', angle_c_degree

  if(timestep.eq.0) then
     open(unit = 11, file = trim(folder)//'angle_c.dat', status = 'replace')
     write(11, '(A)') 'variables = "time", "contact angle", "dt"'
  else
     open(unit = 11, file = trim(folder)//'angle_c.dat', status = 'old', access = 'append')
  end if

  write(11,'(es15.7,f9.3,es15.7)') time, angle_c_degree, dt
  close(11)

  
  !----------------substrate & surface , surface concentration---------------------------------
  !substrate
  if(timestep.gt.0 .and. solve_cp.eq.1 .and. sub_adsp.eq.1) then
     allocate( cpsolp(NTN) )
     if(timestep.eq.1) then
        open(unit = 32, file = trim(folder)//'sub_surf_concen.dat', status = 'replace')
     else
        open(unit = 32, file = trim(folder)//'sub_surf_concen.dat', status = 'old', access = 'append')
     end if
     write(32, '(A)') 'variables = "r", "<greek>Gamma</greek>"'
     write(32, '(A,f6.3,A)') 'Zone T = "theta=', angle_c_degree, '"'
     
     do i = 1, NTN
        if(BCflagN(i,2).ne.0 .and. VN(i).eq.0) then
           cpsolp(i) = solp( NOPP(i) + Ncp ) 
           gammasol(i) = gammasol(i) + Da_sub/Pep*( cpsol(i)+cpsolp(i) )/2.0_rk*dt
           write(32,'(2es15.7)')  rcoordinate(i), gammasol(i)
        end if
     end do
     
     close(32)
  end if

  !free surface
  if(timestep.gt.0 .and. solve_cp.eq.1 .and. surf_adsp.eq.1) then
     if(timestep.eq.1) then
        open(unit = 32, file = trim(folder)//'surf_surf_concen.dat', status = 'replace')
     else
        open(unit = 32, file = trim(folder)//'surf_surf_concen.dat', status = 'old', access = 'append')
     end if
     write(32, '(A)') 'variables = "r", "<greek>Gamma</greek>"'
     write(32, '(A,f6.3,A)') 'Zone T = "theta=', angle_c_degree, '"'
     
     do i = 1, NTN
        if(BCflagN(i,3).ne.0 .and. VN(i).eq.0) then
           write(32,'(2es15.7)')  rcoordinate(i), gammasol(i)
        end if
     end do
     
     close(32)
  end if
  !--------------------------------------------------------------------------------

  


  !--------------------------total particle mass & drop volume-------------------------
  particle_sub = 0.0_rk
  particle_free = 0.0_rk
  particle_m = 0.0_rk
  !substrtae particle
  if(final_size.eq.1 .and. solve_cp.eq.1 .and. sub_adsp.eq.1) then
     particle_sub = 0.0_rk
     do i = 1, NTE
        if(BCflagE(i,2).eq.0 .or. VE(i).ne.0) cycle
        do j = 1, 3
           if(BCflagE(i,2).eq.1) then
              jpp = j
           else! BCflagE=2
              jpp = j*3-2
           end if
           rlocal(j,1) = rcoordinate( globalNM(i,jpp) )
           gammalocal(j,1) = gammasol( globalNM(i,jpp) )
        end do
        do k = 1, Ng, 1
           rintfac = 0.0_rk
           rsi = 0.0_rk
           gammaintfac = 0.0_rk
           do n = 1, 3
              rintfac = rintfac + rlocal(n,1)*phi_1d(k,n)
              rsi = rsi + rlocal(n,1)*phix_1d(k,n)
              gammaintfac = gammaintfac + gammalocal(n,1)*phi_1d(k,n)
           end do
           intsurfMass(k) = gammaintfac * rintfac * abs(rsi)
        end do
        particle_sub = particle_sub + gaussian_quadrature_1d(intsurfMass)
     end do !NTE
     ! write(*,*) 'particle mass', particle_m
  end if

  !free surface particle
  if(final_size.eq.1 .and. solve_cp.eq.1 .and. surf_adsp.eq.1) then
     particle_free = 0.0_rk
     do i = 1, NTE
        if(BCflagE(i,3).ne.1 .or. VE(i).ne.0) cycle
        do j = 1, 3
           jpp = j*3
           rlocal(j,1) = rcoordinate( globalNM(i,jpp) )
           zlocal(j,1) = zcoordinate( globalNM(i,jpp) )
           gammalocal(j,1) = gammasol( globalNM(i,jpp) )
        end do
        do k = 1, Ng, 1
           rintfac = 0.0_rk
           reta = 0.0_rk
           zeta = 0.0_rk
           gammaintfac = 0.0_rk
           do n = 1, 3
              rintfac = rintfac + rlocal(n,1)*phi_1d(k,n)
              reta = reta + rlocal(n,1)*phix_1d(k,n)
              zeta = zeta + zlocal(n,1)*phix_1d(k,n)
              gammaintfac = gammaintfac + gammalocal(n,1)*phi_1d(k,n)
           end do
           intsurfMass(k) = gammaintfac * rintfac * sqrt(reta**2+zeta**2)
        end do
        particle_free = particle_free + gaussian_quadrature_1d(intsurfMass)
     end do !NTE
     ! write(*,*) 'particle mass', particle_m
  end if


  !bulk phase particle & drop volume
  if(final_size.eq.1) then
     particle_m = 0.0_rk
     volume1 = 0.0_rk
     do i = 1, NTE
        if(VE(i).ne.0) cycle
        do j = 1, 9
           rlocal(j,1) = rcoordinate( globalNM(i,j) )
           zlocal(j,1) = zcoordinate( globalNM(i,j) )
           if(solve_cp.eq.1) cplocal(j,1) = cpsol( globalNM(i,j) )
        end do
        do k = 1, Ng, 1
           do l = 1, Ng, 1
              rintfac = 0.0_rk
              rsi = 0.0_rk
              reta = 0.0_rk
              zsi = 0.0_rk
              zeta = 0.0_rk
              cpintfac = 0.0_rk
              do n = 1, 9, 1
                 rintfac = rintfac + rlocal(n,1)*phi(k,l,n)
                 if(solve_cp.eq.1) &
                      cpintfac = cpintfac + cplocal(n,1)*phi(k,l,n)
                 rsi = rsi + rlocal(n,1)*phisi(k,l,n)
                 reta = reta + rlocal(n,1)*phieta(k,l,n)
                 zsi = zsi + zlocal(n,1)*phisi(k,l,n)
                 zeta = zeta + zlocal(n,1)*phieta(k,l,n)
              end do
              !define Jp(3,3)
              Jp =  rsi*zeta - reta*zsi
              if(solve_cp.eq.1) &
                   intMass(k,l) = cpintfac * rintfac * abs( Jp )
              intVol(k,l) = rintfac * abs( Jp )
           end do
        end do
        if(solve_cp.eq.1) &
             particle_m = particle_m + gaussian_quadrature(intMass)
        volume1 = volume1 + gaussian_quadrature(intVol)
        ! if(timestep.eq.1) volume0 = volume1
     end do !NTE
     ! write(*,*) 'particle mass', particle_m

     !write particle mass
     if(solve_cp.eq.1) then
        if(timestep.eq.0) then
           open(unit = 31, file = trim(folder)//'particle_mass.dat', status = 'replace')
           write(31, '(A)') 'variables = "contact angle", "m_bulk", "m_sub", "m_free", "m_total", "time" '
        else
           open(unit = 31, file = trim(folder)//'particle_mass.dat', status = 'old', access = 'append')
        end if
        write(31, '(6es15.7)') angle_c_degree, particle_m, particle_sub, particle_free, &
             particle_m+particle_sub+particle_free, time
        close(31)
     end if

     ! !write drop volume
     ! write(*,*) 'drop volume', volume1, 'particle mass limit', volume1*(cpmax+1.0_rk)
     ! if(timestep.eq.0) then
     !    open(unit = 32, file = trim(folder)//'drop_volume.dat', status = 'replace')
     !    write(32, '(A)') 'variables = "time", "V", "contact angle" '
     ! else
     !    open(unit = 32, file = trim(folder)//'drop_volume.dat', status = 'old', access = 'append')
     ! end if
     ! write(32, '(3es15.7)') time, volume1, angle_c
     ! close(32)
  end if
  
  !---------------------------------------------------------------------------------------


  !in special_points module
  call turning_points(cut, delta_eta)
  call stagnation_and_extremum(cut, delta_eta)
  

  
  ! radial particle mass
  if(solve_cp.eq.1) then
     if(timestep.ge.10 .and. mod(timestep,graph_step).eq.0) then
        if(angle_c_degree.le.angle_int) then
           if(real(int(angle_c_degree),rk) .eq. angle_c_degree) then
              angle_int = angle_c_degree - 1.0_rk
           else
              angle_int = real(int(angle_c_degree),rk)
           end if
           ! if(timestep.ge.1) then
           radial_cal_time = radial_cal_time + 1
           call radial_accumulation
           !       pause
           ! end if
        end if
     end if
  end if
           
  ! ! packing
  ! if(solve_cp.eq.1) then
  !    !flag to judge if maximum packing everywhere
  !    pack_condition = 1.0_rk
  !    do i = 1, NTN
  !       if(pack_flag(i).eq.2) cycle
  !       pack_condition = pack_condition * real(pack_flag(i),rk)
  !       if(pack_condition.eq.0.0_rk) exit
  !    end do
  ! end if




  !----------------------------temperature&cp of free surface-------------------------
  if(timestep.ne.0) then
     
     if(solve_T.eq.1) then
        if(timestep.eq.1 ) then
           open(unit = 13, file = trim(folder)//'temp_surface.dat', status = 'replace')   
        else
           open(unit = 13, file = trim(folder)//'temp_surface.dat', status = 'old', access = 'append')
        end if
     end if

     if(solve_cp.eq.1) then
        if(timestep.eq.1) then 
           open(unit = 113, file = trim(folder)//'cp_surface.dat', status = 'replace')    
        else
           open(unit = 113, file = trim(folder)//'cp_surface.dat', status = 'old', access = 'append')
        end if

        if(timestep.eq.1) then 
           open(unit = 114, file = trim(folder)//'normalized_cp_surface.dat', status = 'replace')    
        else
           open(unit = 114, file = trim(folder)//'normalized_cp_surface.dat', status = 'old', access = 'append')
        end if
     end if


     if(timestep.le.5 .or. mod(timestep,graph_step).eq.0) then   !write data every several timestep
        if(solve_T.eq.1) then
           write(13, '(A)') 'variables = "r", "T"'
           write(13, '(A,f6.3,A)') 'Zone T = "theta=', angle_c_degree, '"'
           do i = 1, NTN
              if( ( ( BCflagN(i,3).eq.1 .or. BCflagN(i,3).eq.3 ) .and. PN(i).eq.1) .or. &
                   ( VN(i).eq.1 .and. BCflagN(i,2).eq.1 ) )  &
                   write(13,'(2es15.7)')  rcoordinate(i), Tsol(i)
           end do
        end if

        if(solve_cp.eq.1) then
           write(113, '(A)') 'variables = "r", "cp"'
           write(113, '(A,f6.3,A)') 'Zone T = "theta=', angle_c_degree, '"'
           do i = 1, NTN
              if( ( ( BCflagN(i,3).eq.1 .or. BCflagN(i,3).eq.3 ) .and. PN(i).eq.1) )  &
                   write(113,'(2es15.7)')  rcoordinate(i), cpsol(i)
           end do

           cp_average = particle_m/volume1
           print *, 'cp_average', cp_average
           write(114, '(A)') 'variables = "r", "cp"'
           write(114, '(A,f6.3,A)') 'Zone T = "theta=', angle_c_degree, '"'
           do i = 1, NTN
              if( ( ( BCflagN(i,3).eq.1 .or. BCflagN(i,3).eq.3 ) .and. PN(i).eq.1) )  &
                   write(114,'(2es15.7)')  rcoordinate(i), cpsol(i)/cp_average
           end do
        end if !solve_cp.eq.1
        
     end if!write data every several timestep

     
     if(solve_T.eq.1) close(13)
     if(solve_cp.eq.1) then
        close(113)
        close(114)
     end if

  end if  !timestep.ne.0


!   !-------------------------------------drop volume------------------------------------
!   if(timestep.eq.1) then
!      open(unit = 12, file = trim(folder)//'volume.dat', status = 'replace')
!      write(12, '(A)') 'variables = time, EvapSpeed, volume'!1, volume1+VolEvap1, volume1+VolEvap2, volume1+VolEvap3'
!   else
!      open(unit = 12, file = trim(folder)//'volume.dat', status = 'old', access = 'append')
!   end if

!   call drop_volume(volume1, volume2)
!   write(12,'(es13.6, f9.3, 2es13.6)') time, angle_c_degree, EvapSpeed, volume1!, volume1+VolEvap1, volume1+VolEvap2, volume1+VolEvap3

  !   close(12)
  

!   !--------------------------------------max value-----------------------------------
!   if(timestep.eq.1) then
!      open(unit = 19, file = trim(folder)//'max_v.dat', status = 'replace')
!      write(19,'(A)') 'variables = "contact angle", "umax", "vmax"'
!   else
!      open(unit = 19, file = trim(folder)//'max_v.dat', status = 'old', access = 'append')
!   end if

!   write(19,'(f9.3, 2es13.6)')  angle_c_degree, umax, vmax

!   close(19)


! !-------------------------------------peclet number---------------------------------
!   if(Maran_flow.eq.1) then
!      Pe_change = Pe*vmax*ztop

!      if(timestep.eq.1) then
!         open(unit = 22, file = trim(folder)//'Pe.dat', status = 'replace')
!         write(22,'(A)') 'variables = "contact angle", "Pe"'
!      else
!         open(unit = 22, file = trim(folder)//'Pe.dat', status = 'old', access = 'append')
!      end if
!      write(22,'(f6.3, es13.6)') angle_c_degree, Pe_change

!      close(22)
!   end if


!   !-----------------------------pressure of free surface-----------------------------
!   if(timestep.eq.1) then
!      open(unit = 20, file = trim(folder)//'pressure.dat', status = 'replace')
!   else
!      open(unit = 20, file = trim(folder)//'pressure.dat', status = 'old', access = 'append')
!   end if

!   if(timestep.le.5 .or. mod(timestep,graph_step).eq.0) then   !write data every several timestep
!      write(20, '(A)') 'variables = "r", "p"'
!      write(20, '(A,f6.3,A)') 'Zone T = "theta=', angle_c_degree, '"'
!      p0 = 2.0_rk*sin(angle_c_degree/180.0_rk*pi)
!      do i = 1, NTN
!         if( ( ( BCflagN(i,3).eq.1 .or. BCflagN(i,3).eq.3 ) .and. PN(i).eq.1) )  &
!              write(20,'(2es15.7)')  rcoordinate(i), psol(i)-p0
!      end do
!   end if

!   close(20)



  
!   !-------------------------------graph spherical cap--------------------------------
!   if(timestep.eq.1) then
!      open(unit = 15, file = trim(folder)//'sphe_cap.dat', status = 'replace')
!      open(unit = 16, file = trim(folder)//'err_sphe.dat', status = 'replace')
!      write(16, '(A)') 'variables = "spherical angle", "contact angle", "error"'
!   else
!      open(unit = 15, file = trim(folder)//'sphe_cap.dat', status = 'old', access = 'append')
!      open(unit = 16, file = trim(folder)//'err_sphe.dat', status = 'old', access = 'append')
!   end if

!   if(timestep.le.5 .or. mod(timestep,graph_step).eq.0) then   !write data every several timestep
!      Rp = (R**2 + ztop**2) / 2.0_rk/ztop
!      angle_c_sphe = atan(R/(Rp-ztop))/pi*180.0_rk
!      write(15, '(A)') 'variables = "r", "z"'
!      write(15,200) 'Zone T = "step:', timestep, '", STRANDID = 1, SOLUTIONTIME =', time, &
!           ', AUXDATA angle_sphe = "', angle_c_sphe, '"'
! 200  format(A,i8,A,es14.7,A,f7.3,A)

!      err_sphe = 0.0_rk
!      do i = 1, NTN
!         if( BCflagN(i,3).eq.1 .or. BCflagN(i,3).eq.3 ) then
!            z_sphe = -Rp + ztop + sqrt(Rp**2-rcoordinate(i)**2)
!            if( PN(i).eq.1 ) write(15,'(2es15.7)') rcoordinate(i), z_sphe
!            if( z_sphe - zcoordinate(i) .gt.err_sphe ) err_sphe = z_sphe - zcoordinate(i)
!            ! err_sphe = err_sphe + ( z_sphe - zcoordinate(i) )**2
!         end if
!      end do
!      ! write(15,'(A,f7.3,A,i7)') 'Text X=40, Y=90, F=Times, T= "contact angle =', angle_c_sphe , &
!      !      '", ZN= ', timestep
!      write(16,'(2f7.3,es15.7)') angle_c_sphe, angle_c_degree, err_sphe
     
!   end if

!   close(15)
!   close(16)


! !   !---------------------stress on the surface grad(p) & grad(sigma)--------------------
! !   if(timestep.eq.1) then
! !      open(unit = 17, file = trim(folder)//'surf_stress.dat', status = 'replace')
! !      open(unit = 21, file = trim(folder)//'Marangoni_stress.dat', status = 'replace')
! !      open(unit = 23, file = trim(folder)//'pressure_change.dat', status = 'replace')
! !      open(unit = 24, file = trim(folder)//'dTds.dat', status = 'replace')
! !      open(unit = 25, file = trim(folder)//'v_lubrication.dat', status = 'replace')
! !   else
! !      open(unit = 17, file = trim(folder)//'surf_stress.dat', status = 'old', access = 'append')
! !      open(unit = 21, file = trim(folder)//'Marangoni_stress.dat', status = 'old', access = 'append')
! !      open(unit = 23, file = trim(folder)//'pressure_change.dat', status = 'old', access = 'append')
! !      open(unit = 24, file = trim(folder)//'dTds.dat', status = 'old', access = 'append')
! !      open(unit = 25, file = trim(folder)//'v_lubrication.dat', status = 'old', access = 'append')
! !   end if

! !   if(timestep.le.5 .or. mod(timestep,graph_step).eq.0) then
! !      write(17, '(A)') 'variables = "r", "ratio"'
! !      write(21, '(A)') 'variables = "r", "Marangoni stress"'
! !      write(23, '(A)') 'variables = "r", "dP/ds"'
! !      write(24, '(A)') 'variables = "r", "dT/ds"'
! !      write(25, '(A)') 'variables = "r", "v_lubrication"'

! !      write(17, '(A,f6.3,A)') 'Zone T = "theta=', angle_c_degree, '"'
! !      write(21, '(A,f6.3,A)') 'Zone T = "theta=', angle_c_degree, '"'
! !      write(23, '(A,f6.3,A)') 'Zone T = "theta=', angle_c_degree, '"'
! !      write(24, '(A,f6.3,A)') 'Zone T = "theta=', angle_c_degree, '"'
! !      write(25, '(A,f6.3,A)') 'Zone T = "theta=', angle_c_degree, '"'

! !      do i = 1, NTE
! !         if(BCflagE(i,3).eq.1) then
! !            !at eta = 0
! !            retap(1) = -3.0_rk*rcoordinate(globalNM(i,3)) + 4.0_rk*rcoordinate(globalNM(i,6)) - rcoordinate(globalNM(i,9))
! !            zetap(1) = -3.0_rk*zcoordinate(globalNM(i,3)) + 4.0_rk*zcoordinate(globalNM(i,6)) - zcoordinate(globalNM(i,9))
! !            peta = -psol(globalNM(i,3)) + psol(globalNM(i,9))
! !            Teta(1) = -3.0_rk*Tsol(globalNM(i,3)) + 4.0_rk*Tsol(globalNM(i,6)) - Tsol(globalNM(i,9))
! !            gradP = - peta/sqrt( retap(1)**2+zetap(1)**2 )  != grad(p)
! !            gradT(1) = - Teta(1)/sqrt( retap(1)**2+zetap(1)**2 )
! !            MaranD = - beta*Teta(1)/sqrt( retap(1)**2+zetap(1)**2 )  != grad(sigma)
! ! ! write(*,*) presD, MaranD/beta, MaranD, Kdi
! !            h_surf = zcoordinate( globalNM(i,3) ) 
! !            dPdr(1) = peta/retap(1)
! !            v_surf_p(1) = ( -0.5_rk*dPdr(1)*h_surf + beta*gradT(1) ) *h_surf

! !            if( rcoordinate(globalNM(i,3)).ne.1.0_rk ) &
! !            write(17,'(2es15.7)') rcoordinate(globalNM(i,3)), -gradP*h_surf/2.0_rk /MaranD!abs(gradP*h_surf/2.0_rk /MaranD)
! !            write(23,'(2es15.7)') rcoordinate(globalNM(i,3)), gradP
! !            write(21,'(2es15.7)') rcoordinate(globalNM(i,3)), MaranD
! !            write(24,'(2es15.7)') rcoordinate(globalNM(i,3)), MaranD/beta
! !            write(25,'(2es15.7)') rcoordinate(globalNM(i,3)), v_surf_p(1)

! !         end if
! !      end do
! !   end if

! !   close(17)
! !   close(23)
! !   close(21)
! !   close(24)
! !   close(25)

  


!   !---------------------------flux on nodes of free surface---------------------------
!   if(timestep.eq.1) then
!      open(unit = 10, file = trim(folder)//'flux.dat', status = 'replace')
!   else
!      open(unit = 10, file = trim(folder)//'flux.dat', status = 'old', access = 'append')
!   end if

!   allocate(flux(NTN), flux1(NTN), Dflux(NTN))
!   J0 = 0.0_rk
!   flux(:) = 0.0_rk
!   Dflux(:) = 0.0_rk

!   if(timestep.le.5 .or. mod(timestep,graph_step).eq.0) then   !write data every several timestep

!      write(10, '(A)') 'variables = "r", "J"'!, "Deegan_flux"'
!      write(10, '(A,f6.3,A)') 'Zone T = "angle =', angle_c_degree, '"'

!      do i = 1, NTN
!         if( ( BCflagN(i,3).eq.1 .or. BCflagN(i,3).eq.3 ) .and. PN(i).eq.1) then

!            !(r,z,c)eta
!            if(i.eq.top_node) then      !use the below element
!               n1 = i
!               do j = i-1, 1, -1
!                  if(VN(j).eq.2) then
!                     n2 = j
!                     do k = j-1, 1, -1
!                        if(VN(k).eq.2) then
!                           n3 = k
!                           exit
!                        end if
!                     end do
!                     exit
!                  end if
!               end do
!               reta = 3.0_rk*rcoordinate(n1) - 4.0_rk*rcoordinate(n2) + rcoordinate(n3)
!               zeta = 3.0_rk*zcoordinate(n1) - 4.0_rk*zcoordinate(n2) + zcoordinate(n3)
!               ceta = 3.0_rk*csol(n1) - 4.0_rk*csol(n2) + csol(n3)

!            else                 !use the above element
!               n1 = i
!               do j = i+1, NTN
!                  if(BCflagN(j,3).eq.1) then
!                     n2 = j
!                     do k = j+1, NTN
!                        if(BCflagN(k,3).eq.1) then
!                           n3 = k
!                           exit
!                        end if
!                     end do
!                     exit
!                  end if
!               end do
!               ! write(*,*) n1,n2,n3
!               ! pause
!               reta = -3.0_rk*rcoordinate(n1) + 4.0_rk*rcoordinate(n2) - rcoordinate(n3)
!               zeta = -3.0_rk*zcoordinate(n1) + 4.0_rk*zcoordinate(n2) - zcoordinate(n3)
!               ! ceta = -3.0_rk*csol(n1) + 4.0_rk*csol(n2) - csol(n3)
!            end if
!            ! !(r,z,c)si
!            ! rsi = -3.0_rk*rcoordinate(i) + 4.0_rk*rcoordinate(i+1) - rcoordinate(i+2)
!            ! zsi = -3.0_rk*zcoordinate(i) + 4.0_rk*zcoordinate(i+1) - zcoordinate(i+2)
!            ! csi = -3.0_rk*csol(i) + 4.0_rk*csol(i+1) - csol(i+2)

!            !flux
!            ! flux(i) = 1.0_rk / ( rsi*zeta - reta*zsi ) * sqrt(reta**2 + zeta**2) *  (-csi)
!            rdot = soldot( NOPP(i) + Nr )
!            zdot = soldot( NOPP(i) + Nz )!??
!            flux(i) = ( zeta*(usol(i)-rdot) - reta*(vsol(i)-zdot) )/ sqrt(reta**2+zeta**2)
!            ! flux(i) = 1.0_rk / ( rsi*zeta - reta*zsi ) * sqrt(reta**2 + zeta**2) *  (-csi)
!            if(i.eq.top_node) J0 = flux(i)

!            Dflux(i) = (1.0_rk - rcoordinate(i)**2)**( -(0.5_rk-angle_c/pi) )

!         end if
!      end do


!      do i = 1, NTN
!         if( ( BCflagN(i,3).eq.1 .or. BCflagN(i,3).eq.3 ) .and. PN(i).eq.1) then 
!            ! flux(i) = flux(i)/J0
!            ! if(i.eq.1) cycle  !put infinity on hold

!            write(10,'(4es15.7)')  rcoordinate(i), flux(i)!, Dflux(i)!, abs(flux(i)-Dflux(i))

!         end if
!      end do

!   end if  !every 100 steps

!   close(10)

!   ! !flux on gausspoints of free surface
!   ! allocate( flux_gp(Ng) )
!   ! write(10, '(A)') 'variables = "r", "flux"'
!   ! write(10, '(A,f6.3,A)') 'Zone T = "t=', time, '"'
!   ! do i = 1, NTE
!   !    if(BCflagE(i,3).ne.2) cycle

!   !    call values_in_an_element(i,1)

!   !    do k = 1, Ng
!   !       flux_gp(k) = 1.0_rk / ( rsi_right(k,1)*zeta_right(k,1) - reta_right(k,1)*zsi_right(k,1) ) / &
!   !            sqrt(reta_right(k,1)**2 + zeta_right(k,1)**2) * &
!   !            ( -dcdsi(k,1) * ( reta_right(k,1)**2 + zeta_right(k,1)**2 ) + &
!   !            dcdeta(k,1) * ( rsi_right(k,1)*reta_right(k,1) + zsi_right(k,1)*zeta_right(k,1) ) )
!   !       write(10,'(2es15.7)')  rintfac_right(k,1), flux_gp(k)
!   !    end do

!   ! end do

!   deallocate(flux, flux1, Dflux)


  

  t = REAL(omp_get_wtime(),rk) - t
  open(unit = 200, file = trim(folder)//'cal_time.dat', status = 'old', access = 'append')
  write(200, '(A)') ' '
  write(200,'(A,es13.6)') 'post process:', t
  close(200)

  !pause

  return



contains






  






  subroutine drop_volume(vol1, vol2)
    use Ldata
    use NOP_mod, only: gaussian_quadrature, gaussian_quadrature_1d
    implicit none

    integer(kind=ik):: i,j, k,l, n, npp
    real(kind=rk), intent(out):: vol1, vol2
    real(kind=rk):: rlocal1(bas), zlocal1(bas), clocal1(bas), &
         rintfac1(Ng,Ng), rsi1(Ng,Ng), reta1(Ng,Ng), zsi1(Ng,Ng), zeta1(Ng,Ng), Jp1(Ng,Ng), intLocVol1(Ng,Ng), &
         rintfac_right1(Ng), zintfac_right1(Ng), reta_right1(Ng), intLocVol2(Ng), &
         zeta_right1(Ng), rsi_right1(Ng), zsi_right1(Ng), csi_right1(Ng), ceta_right1(Ng), &
         intLocEvap(Ng), Jp2(Ng)


    vol1 = 0.0_rk
    vol2 = 0.0_rk
    EvapSpeed = 0.0_rk

    do i = 1, NTE

       if(VE(i).eq.0) then

          do j = 1, bas   
             rlocal1(j) = sol( NOPP( globalNM(i,j) ) + Nr )
             zlocal1(j) = sol( NOPP( globalNM(i,j) ) + Nz )
          end do

          !for vol1
          rintfac1(:,:) = 0.0_rk
          rsi1(:,:) = 0.0_rk
          reta1(:,:) = 0.0_rk
          zsi1(:,:) = 0.0_rk
          zeta1(:,:) = 0.0_rk
          Jp1(:,:) = 0.0_rk
          do k = 1, Ng
             do l = 1, Ng
                !define rsi1, reta1, zsi1, zeta1
                do n = 1, bas
                   rintfac1(k,l) = rintfac1(k,l) + rlocal1(n)*phi(k,l,n)
                   rsi1(k,l) = rsi1(k,l) + rlocal1(n)*phisi(k,l,n)
                   reta1(k,l) = reta1(k,l) + rlocal1(n)*phieta(k,l,n)
                   zsi1(k,l) = zsi1(k,l) + zlocal1(n)*phisi(k,l,n)
                   zeta1(k,l) = zeta1(k,l) + zlocal1(n)*phieta(k,l,n)
                end do
                Jp1(k,l) =  rsi1(k,l)*zeta1(k,l) - reta1(k,l)*zsi1(k,l) 
                intLocVol1(k,l) = rintfac1(k,l) * abs( Jp1(k,l) )
             end do
          end do

          vol1 = vol1 + gaussian_quadrature(intLocVol1)

       end if


       !for vol2
       if(BCflagE(i,3).eq.1) then
          rintfac_right1(:) = 0.0_rk
          zintfac_right1(:) = 0.0_rk
          reta_right1(:) = 0.0_rk
          zeta_right1(:) = 0.0_rk
          do k = 1, Ng
             do n = 1, 3 !the summation of three terms, eg: rsi = sum( rlocal1(7,8,9)*phix_1d(1,2,3) )
                npp = 3*n  !only need r3, r6, r9 -->r(npp)
                rintfac_right1(k) = rintfac_right1(k) + rlocal1(npp) * phi_1d(k,n)
                zintfac_right1(k) = zintfac_right1(k) + zlocal1(npp) * phi_1d(k,n)
                reta_right1(k) = reta_right1(k) + rlocal1(npp) * phix_1d(k,n)
             end do
             intLocVol2(k) = -rintfac_right1(k)*zintfac_right1(k)*reta_right1(k)
          end do
          vol2 = vol2 + gaussian_quadrature_1d(intLocVol2)
       end if




       !for evaporating volume
       if(BCflagE(i,3).eq.2) then

          do j = 1, bas   
             rlocal1(j) = sol( NOPP( globalNM(i,j) ) + Nr )
             zlocal1(j) = sol( NOPP( globalNM(i,j) ) + Nz )
             clocal1(j) = sol( NOPP( globalNM(i,j) ) + MDF( globalNM(i,j) )-1 )
          end do

          rintfac_right1(:) = 0.0_rk
          reta_right1(:) = 0.0_rk 
          zeta_right1(:) = 0.0_rk
          ceta_right1(:) = 0.0_rk
          rsi_right1(:) = 0.0_rk
          zsi_right1(:) = 0.0_rk
          csi_right1(:) = 0.0_rk
          do k = 1, Ng
             do n = 1, 3 !the summation of three terms, eg: rsi = sum( rlocal1(7,8,9)*phix_1d(1,2,3) )
                npp = 3*n-2  !only need r1, r4, r7 -->r(npp)
                rintfac_right1(k) = rintfac_right1(k) + rlocal1(npp) * phi_1d(k,n)
                reta_right1(k) = reta_right1(k) + rlocal1(npp) * phix_1d(k,n)
                zeta_right1(k) = zeta_right1(k) + zlocal1(npp) * phix_1d(k,n)
                ceta_right1(k) = ceta_right1(k) + clocal1(npp) * phix_1d(k,n)
             end do
             do n = 1, 9              
                rsi_right1(k) = rsi_right1(k) + rlocal1(n) * phisi0_1d(k,n)
                zsi_right1(k) = zsi_right1(k) + zlocal1(n) * phisi0_1d(k,n)
                csi_right1(k) = csi_right1(k) + clocal1(n) * phisi0_1d(k,n)
             end do
             Jp2(k) =  rsi_right1(k)*zeta_right1(k) - reta_right1(k)*zsi_right1(k)
             intLocEvap(k) = ( -csi_right1(k)*( zeta_right1(k)**2 + reta_right1(k)**2 ) + &
                  ceta_right1(k)*( zsi_right1(k)*zeta_right1(k) + rsi_right1(k)*reta_right1(k) ) ) &
                  *rintfac_right1(k)/Jp2(k)
             !write(*,*)  csi_right1(k), zeta_right1(k), reta_right1(k),ceta_right1(k), zsi_right1(k), rsi_right1(k)
          end do
          EvapSpeed = EvapSpeed + gaussian_quadrature_1d(intLocEvap)/KBCgroup
       end if

    end do   !for i

    if(timestep.gt.0) then
       VolEvap1 = VolEvap1 + EvapSpeedp*dt
       VolEvap2 = VolEvap2 + EvapSpeed*dt
       VolEvap3 = VolEvap3 + (EvapSpeed+ EvapSpeedp)/2.0_rk*dt
    end if
    EvapSpeedp = EvapSpeed

    return
  end subroutine drop_volume





end subroutine variable_cal
