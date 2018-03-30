
subroutine variable_cal
  use kind
  use omp_lib
  use data
  use Ldata, only: dcdsi, dcdeta, rsi_right, zsi_right, reta_right, zeta_right, rintfac_right
  use basis_f!, only: phii_1d, phiix_1d
  implicit none


  integer(kind=ik):: i,j,k,n1, n2, n3, angle_c_node
  real(kind=rk):: rsi, reta, zsi, zeta, csi, ceta, rsi1, reta1, zsi1, zeta1, csi1, ceta1, rdot, zdot
  real(kind=rk), allocatable:: flux(:), flux1(:), Dflux(:), flux_gp(:)
  real(kind=rk):: J0, volume1, volume2
  ! real(kind=rk):: a, b, c, eta0, eta01, eta02, r_change
  real(kind=rk):: retap(3), zetap(3), v_surf(3), eta1, eta2, eta3, usolp, vsolp, r_change
  real(kind=rk):: Rp, angle_c_sphe, err_sphe, z_sphe
  real(kind=rk):: MaranD, gradP, peta, Teta(3), gradT(3)
  real(kind=rk):: p0, Pep
  real(kind=rk):: t
  real(kind=rk):: v_surf_p(3), h_surf, dPdr(3), zsolp


  t = REAL(omp_get_wtime(),rk)

  ! if(timestep.eq.1) then
  !    open(unit = 10, file = trim(folder)//'flux.dat', status = 'replace')
  !    open(unit = 13, file = trim(folder)//'temp_surface.dat', status = 'replace')
     
  !    open(unit = 14, file = trim(folder)//'surf_flow_dir.dat', status = 'replace')
  !    write(14, '(A)') 'variables = "contact angle", "r", "time" '

  !    open(unit = 18, file = trim(folder)//'surf_gradT_dir.dat', status = 'replace')
  !    write(18, '(A)') 'variables = "contact angle", "r", "time" '
     
  !    open(unit = 11, file = trim(folder)//'angle_c.dat', status = 'replace')
  !    write(11, '(A)') 'variables = "time", "contact angle", "dt"'

  !    open(unit = 12, file = trim(folder)//'volume.dat', status = 'replace')
  !    write(12, '(A)') 'variables = time, EvapSpeed, volume'!1, volume1+VolEvap1, volume1+VolEvap2, volume1+VolEvap3'

  !    open(unit = 15, file = trim(folder)//'sphe_cap.dat', status = 'replace')
  !    open(unit = 16, file = trim(folder)//'err_sphe.dat', status = 'replace')
  !    write(16, '(A)') 'variables = "spherical angle", "contact angle", "error"'

  !    open(unit = 17, file = trim(folder)//'surf_stress.dat', status = 'replace')
  !    open(unit = 21, file = trim(folder)//'Marangoni_stress.dat', status = 'replace')
  !    open(unit = 23, file = trim(folder)//'pressure_change.dat', status = 'replace')

  !    open(unit = 19, file = trim(folder)//'max_v.dat', status = 'replace')
  !    write(19,'(A)') 'variables = "contact angle", "umax", "vmax"'

  !    open(unit = 20, file = trim(folder)//'pressure.dat', status = 'replace')

  !    open(unit = 22, file = trim(folder)//'Pe.dat', status = 'replace')
  !    write(22,'(A)') 'variables = "contact angle", "Pe"'

  ! else
  !    open(unit = 10, file = trim(folder)//'flux.dat', status = 'old', access = 'append')
  !    open(unit = 13, file = trim(folder)//'temp_surface.dat', status = 'old', access = 'append')
  !    open(unit = 14, file = trim(folder)//'surf_flow_dir.dat', status = 'old', access = 'append')
  !    open(unit = 11, file = trim(folder)//'angle_c.dat', status = 'old', access = 'append')
  !    open(unit = 12, file = trim(folder)//'volume.dat', status = 'old', access = 'append')
  !    open(unit = 15, file = trim(folder)//'sphe_cap.dat', status = 'old', access = 'append')
  !    open(unit = 16, file = trim(folder)//'err_sphe.dat', status = 'old', access = 'append')
  !    open(unit = 17, file = trim(folder)//'surf_stress.dat', status = 'old', access = 'append')
  !    open(unit = 18, file = trim(folder)//'surf_gradT_dir.dat', status = 'old', access = 'append')
  !    open(unit = 19, file = trim(folder)//'max_v.dat', status = 'old', access = 'append')
  !    open(unit = 20, file = trim(folder)//'pressure.dat', status = 'old', access = 'append')
  !    open(unit = 21, file = trim(folder)//'Marangoni_stress.dat', status = 'old', access = 'append')
  !    open(unit = 22, file = trim(folder)//'Pe.dat', status = 'old', access = 'append')
  !    open(unit = 23, file = trim(folder)//'pressure_change.dat', status = 'old', access = 'append')
  ! end if


  !contact angle
  do i = 1, NTN
     if( BCflagN(i,3).eq.3 ) then
        exit
     end if
  end do   !i is the contact line node

  if(no_vapor.eq.0) then
     angle_c_node = i + 2*NEV + 2*(NES+1) + 1
  else !no_vapor=1
     angle_c_node = i + 2*(NES+1) + 1
  end if

  angle_c = atan( zcoordinate(angle_c_node)/ ( rcoordinate(i) - rcoordinate(angle_c_node) ) )
  !atan( solp( NOPP(angle_c_node)+Nz ) / ( solp( NOPP(1)+Nr ) - solp( NOPP(angle_c_node)+Nr ) ) )
  angle_c_degree = angle_c /pi*180.0_rk   !degree
  write(*,*) ' '
  write(*,*) 'contact angle', angle_c_degree

  if(timestep.eq.1) then
     open(unit = 11, file = trim(folder)//'angle_c.dat', status = 'replace')
     write(11, '(A)') 'variables = "time", "contact angle", "dt"'
  else
     open(unit = 11, file = trim(folder)//'angle_c.dat', status = 'old', access = 'append')
  end if

  write(11,'(es15.7,f9.3,es15.7)') time, angle_c_degree, dt
  close(11)

  !drop volume
  if(timestep.eq.1) then
     open(unit = 12, file = trim(folder)//'volume.dat', status = 'replace')
     write(12, '(A)') 'variables = time, EvapSpeed, volume'!1, volume1+VolEvap1, volume1+VolEvap2, volume1+VolEvap3'
  else
     open(unit = 12, file = trim(folder)//'volume.dat', status = 'old', access = 'append')
  end if

  call drop_volume(volume1, volume2)
  write(12,'(es13.6, f9.3, 2es13.6)') time, angle_c_degree, EvapSpeed, volume1!, volume1+VolEvap1, volume1+VolEvap2, volume1+VolEvap3

  close(12)



  !max value
  if(timestep.eq.1) then
     open(unit = 19, file = trim(folder)//'max_v.dat', status = 'replace')
     write(19,'(A)') 'variables = "contact angle", "umax", "vmax"'
  else
     open(unit = 19, file = trim(folder)//'max_v.dat', status = 'old', access = 'append')
  end if

  write(19,'(f9.3, 2es13.6)')  angle_c_degree, umax, vmax

  close(19)


  !peclet number
  Pep = Pe*vmax*ztop

  if(timestep.eq.1) then
     open(unit = 22, file = trim(folder)//'Pe.dat', status = 'replace')
     write(22,'(A)') 'variables = "contact angle", "Pe"'
  else
     open(unit = 22, file = trim(folder)//'Pe.dat', status = 'old', access = 'append')
  end if
  write(22,'(f6.3, es13.6)') angle_c_degree, Pep

  close(22)


  !pressure of free surface
  if(timestep.eq.1) then
     open(unit = 20, file = trim(folder)//'pressure.dat', status = 'replace')
  else
     open(unit = 20, file = trim(folder)//'pressure.dat', status = 'old', access = 'append')
  end if

  if(timestep.le.5 .or. mod(timestep,graph_step).eq.0) then   !write data every several timestep
     write(20, '(A)') 'variables = "r", "p"'
     write(20, '(A,f6.3,A)') 'Zone T = "theta=', angle_c_degree, '"'
     p0 = 2.0_rk*sin(angle_c_degree/180.0_rk*pi)
     do i = 1, NTN
        if( ( ( BCflagN(i,3).eq.1 .or. BCflagN(i,3).eq.3 ) .and. PN(i).eq.1) )  &
             write(20,'(2es15.7)')  rcoordinate(i), psol(i)-p0
     end do
  end if

  close(20)




  !temperature of free surface
  if(timestep.eq.1) then
     open(unit = 13, file = trim(folder)//'temp_surface.dat', status = 'replace')     
  else
     open(unit = 13, file = trim(folder)//'temp_surface.dat', status = 'old', access = 'append')
  end if

  if(timestep.le.5 .or. mod(timestep,graph_step).eq.0) then   !write data every several timestep
     write(13, '(A)') 'variables = "r", "T"'
     write(13, '(A,f6.3,A)') 'Zone T = "theta=', angle_c_degree, '"'
     do i = 1, NTN
        if( ( ( BCflagN(i,3).eq.1 .or. BCflagN(i,3).eq.3 ) .and. PN(i).eq.1) .or. &
             ( VN(i).eq.1 .and. BCflagN(i,2).eq.1 ) )  &
             write(13,'(2es15.7)')  rcoordinate(i), Tsol(i)
     end do
  end if

  close(13)


  !graph spherical cap
  if(timestep.eq.1) then
     open(unit = 15, file = trim(folder)//'sphe_cap.dat', status = 'replace')
     open(unit = 16, file = trim(folder)//'err_sphe.dat', status = 'replace')
     write(16, '(A)') 'variables = "spherical angle", "contact angle", "error"'
  else
     open(unit = 15, file = trim(folder)//'sphe_cap.dat', status = 'old', access = 'append')
     open(unit = 16, file = trim(folder)//'err_sphe.dat', status = 'old', access = 'append')
  end if

  if(timestep.le.5 .or. mod(timestep,graph_step).eq.0) then   !write data every several timestep
     Rp = (R**2 + ztop**2) / 2.0_rk/ztop
     angle_c_sphe = atan(R/(Rp-ztop))/pi*180.0_rk
     write(15, '(A)') 'variables = "r", "z"'
     write(15,200) 'Zone T = "step:', timestep, '", STRANDID = 1, SOLUTIONTIME =', time, &
          ', AUXDATA angle_sphe = "', angle_c_sphe, '"'
200  format(A,i8,A,es14.7,A,f7.3,A)

     err_sphe = 0.0_rk
     do i = 1, NTN
        if( BCflagN(i,3).eq.1 .or. BCflagN(i,3).eq.3 ) then
           z_sphe = -Rp + ztop + sqrt(Rp**2-rcoordinate(i)**2)
           if( PN(i).eq.1 ) write(15,'(2es15.7)') rcoordinate(i), z_sphe
           if( z_sphe - zcoordinate(i) .gt.err_sphe ) err_sphe = z_sphe - zcoordinate(i)
           ! err_sphe = err_sphe + ( z_sphe - zcoordinate(i) )**2
        end if
     end do
     ! write(15,'(A,f7.3,A,i7)') 'Text X=40, Y=90, F=Times, T= "contact angle =', angle_c_sphe , &
     !      '", ZN= ', timestep
     write(16,'(2f7.3,es15.7)') angle_c_sphe, angle_c_degree, err_sphe
     
  end if

  close(15)
  close(16)


  !stress on the surface grad(p) & grad(sigma)
  if(timestep.eq.1) then
     open(unit = 17, file = trim(folder)//'surf_stress.dat', status = 'replace')
     open(unit = 21, file = trim(folder)//'Marangoni_stress.dat', status = 'replace')
     open(unit = 23, file = trim(folder)//'pressure_change.dat', status = 'replace')
     open(unit = 24, file = trim(folder)//'dTds.dat', status = 'replace')
     open(unit = 25, file = trim(folder)//'v_lubrication.dat', status = 'replace')
  else
     open(unit = 17, file = trim(folder)//'surf_stress.dat', status = 'old', access = 'append')
     open(unit = 21, file = trim(folder)//'Marangoni_stress.dat', status = 'old', access = 'append')
     open(unit = 23, file = trim(folder)//'pressure_change.dat', status = 'old', access = 'append')
     open(unit = 24, file = trim(folder)//'dTds.dat', status = 'old', access = 'append')
     open(unit = 25, file = trim(folder)//'v_lubrication.dat', status = 'old', access = 'append')
  end if

  if(timestep.le.5 .or. mod(timestep,graph_step).eq.0) then
     write(17, '(A)') 'variables = "r", "ratio"'
     write(21, '(A)') 'variables = "r", "Marangoni stress"'
     write(23, '(A)') 'variables = "r", "dP/ds"'
     write(24, '(A)') 'variables = "r", "dT/ds"'
     write(25, '(A)') 'variables = "r", "v_lubrication"'

     write(17, '(A,f6.3,A)') 'Zone T = "theta=', angle_c_degree, '"'
     write(21, '(A,f6.3,A)') 'Zone T = "theta=', angle_c_degree, '"'
     write(23, '(A,f6.3,A)') 'Zone T = "theta=', angle_c_degree, '"'
     write(24, '(A,f6.3,A)') 'Zone T = "theta=', angle_c_degree, '"'
     write(25, '(A,f6.3,A)') 'Zone T = "theta=', angle_c_degree, '"'

     do i = 1, NTE
        if(BCflagE(i,3).eq.1) then
           !at eta = 0
           retap(1) = -3.0_rk*rcoordinate(globalNM(i,3)) + 4.0_rk*rcoordinate(globalNM(i,6)) - rcoordinate(globalNM(i,9))
           zetap(1) = -3.0_rk*zcoordinate(globalNM(i,3)) + 4.0_rk*zcoordinate(globalNM(i,6)) - zcoordinate(globalNM(i,9))
           peta = -psol(globalNM(i,3)) + psol(globalNM(i,9))
           Teta(1) = -3.0_rk*Tsol(globalNM(i,3)) + 4.0_rk*Tsol(globalNM(i,6)) - Tsol(globalNM(i,9))
           gradP = - peta/sqrt( retap(1)**2+zetap(1)**2 )  != grad(p)
           gradT(1) = - Teta(1)/sqrt( retap(1)**2+zetap(1)**2 )
           MaranD = - beta*Teta(1)/sqrt( retap(1)**2+zetap(1)**2 )  != grad(sigma)
! write(*,*) presD, MaranD/beta, MaranD, Kdi
           h_surf = zcoordinate( globalNM(i,3) ) 
           dPdr(1) = peta/retap(1)
           v_surf_p(1) = ( -0.5_rk*dPdr(1)*h_surf + beta*gradT(1) ) *h_surf

           if( rcoordinate(globalNM(i,3)).ne.1.0_rk ) &
           write(17,'(2es15.7)') rcoordinate(globalNM(i,3)), -gradP*h_surf/2.0_rk /MaranD!abs(gradP*h_surf/2.0_rk /MaranD)
           write(23,'(2es15.7)') rcoordinate(globalNM(i,3)), gradP
           write(21,'(2es15.7)') rcoordinate(globalNM(i,3)), MaranD
           write(24,'(2es15.7)') rcoordinate(globalNM(i,3)), MaranD/beta
           write(25,'(2es15.7)') rcoordinate(globalNM(i,3)), v_surf_p(1)

        end if
     end do
  end if

  close(17)
  close(23)
  close(21)
  close(24)
  close(25)


  !velocity direction change location on free surface & grad(T) direction

  if(timestep.eq.1) then
     open(unit = 18, file = trim(folder)//'surf_gradT_dir.dat', status = 'replace')
     write(18, '(A)') 'variables = "contact angle", "r", "time" '
  else
     open(unit = 18, file = trim(folder)//'surf_gradT_dir.dat', status = 'old', access = 'append')
  end if

  if(timestep.eq.1) then
     open(unit = 30, file = trim(folder)//'v0_lubrication.dat', status = 'replace')
     write(30, '(A)') 'variables = "contact angle", "r", "time" '
  else
     open(unit = 30, file = trim(folder)//'v0_lubrication.dat', status = 'old', access = 'append')
  end if


  if(timestep.eq.1) then
     open(unit = 14, file = trim(folder)//'surf_flow_dir.dat', status = 'replace')
     write(14, '(A)') 'variables = "contact angle", "r", "time" '
  else
     open(unit = 14, file = trim(folder)//'surf_flow_dir.dat', status = 'old', access = 'append')
  end if

  if(timestep.gt.0) then
     r_change = 0.0_rk 
     do i = 1, NTE
        if(BCflagE(i,3).eq.1) then  !surface element
           !the element where direction changes

           !at eta = 0
           retap(1) = -3.0_rk*rcoordinate(globalNM(i,3)) + 4.0_rk*rcoordinate(globalNM(i,6)) - rcoordinate(globalNM(i,9))
           zetap(1) = -3.0_rk*zcoordinate(globalNM(i,3)) + 4.0_rk*zcoordinate(globalNM(i,6)) - zcoordinate(globalNM(i,9))
           v_surf(1) = ( usol(globalNM(i,3))*retap(1) + vsol(globalNM(i,3))*zetap(1) ) /sqrt( retap(1)**2+zetap(1)**2 )
           Teta(1) = -3.0_rk*Tsol(globalNM(i,3)) + 4.0_rk*Tsol(globalNM(i,6)) - Tsol(globalNM(i,9))
           gradT(1) = -Teta(1)/sqrt( retap(1)**2+zetap(1)**2 )
           dPdr(1) = ( psol(globalNM(i,9)) - psol(globalNM(i,3)) ) /retap(1)
           v_surf_p(1) = -0.5_rk*dPdr(1)*zcoordinate(globalNM(i,3)) + beta*gradT(1)

           !at eta = 1
           retap(2) = rcoordinate(globalNM(i,3)) - 4.0_rk*rcoordinate(globalNM(i,6)) + 3.0_rk*rcoordinate(globalNM(i,9))
           zetap(2) = zcoordinate(globalNM(i,3)) - 4.0_rk*zcoordinate(globalNM(i,6)) + 3.0_rk*zcoordinate(globalNM(i,9))
           v_surf(2) = ( usol(globalNM(i,9))*retap(2) + vsol(globalNM(i,9))*zetap(2) ) /sqrt( retap(2)**2+zetap(2)**2 )
           Teta(2) = Tsol(globalNM(i,3)) - 4.0_rk*Tsol(globalNM(i,6)) + 3.0_rk*Tsol(globalNM(i,9))
           gradT(2) = -Teta(2)/sqrt( retap(2)**2+zetap(2)**2 )
           dPdr(2) = ( psol(globalNM(i,9)) - psol(globalNM(i,3)) ) /retap(2)
           v_surf_p(2) = -0.5_rk*dPdr(2)*zcoordinate(globalNM(i,9)) + beta*gradT(2)

           !lubrication velocity direction change
           if( angle_c_degree.lt.15.0_rk .and. &
                v_surf_p(1) * v_surf_p(2) .lt. 0.0_rk .and. &
                (i.ne.top_element .and. i.ne.CL_element) ) then 
 
              eta1 = 0.0_rk
              eta2 = 1.0_rk
              do while ( abs(eta1-eta2).gt.0.5e-1_rk )
                 eta3 = (eta1+eta2)/2.0_rk
                 retap(3) = 0.0_rk
                 zetap(3) = 0.0_rk
                 Teta(3) = 0.0_rk
                 zsolp = 0.0_rk
                 do j = 1, 3
                    retap(3) = retap(3) + rcoordinate(globalNM(i,3*j))*phiix_1d(eta3,j)
                    zetap(3) = zetap(3) + zcoordinate(globalNM(i,3*j))*phiix_1d(eta3,j)
                    Teta(3) = Teta(3) + Tsol(globalNM(i,3*j))*phiix_1d(eta3,j)
                    zsolp = zsolp + zcoordinate(globalNM(i,3*j))*phii_1d(eta3,j)
                 end do
                 gradT(3) = -Teta(3)/sqrt( retap(3)**2+zetap(3)**2 )
                 dPdr(3) = ( psol(globalNM(i,9)) - psol(globalNM(i,3)) ) /retap(3)
                 v_surf_p(3) = -0.5_rk*dPdr(3)*zsolp + beta*gradT(3)
                 if( v_surf_p(3) * v_surf_p(1) .lt. 0.0_rk ) then
                    eta2 = eta3
                 else  !v_surf(3) * v_surf(2) .le. 0.0_rk
                    eta1 = eta3
                 end if
              end do
              r_change = 0.0_rk
              do j = 1, 3
                 r_change = r_change + rcoordinate(globalNM(i,3*j))*phii_1d(eta1,j)
              end do

              write(30, '(f9.3,2es15.7)') angle_c_degree, r_change, time 
              write(*,*) 'r_change for lubrication velocity =', r_change, 'element:', i
              ! exit   !?not strict

           end if   !u change element



           !surface flow direction change
           if( v_surf(1) * v_surf(2) .lt. 0.0_rk .and. (i.ne.top_element .and. i.ne.CL_element) ) then
! write(*,*)  v_surf(1), v_surf(2)
              eta1 = 0.0_rk
              eta2 = 1.0_rk
              do while ( abs(eta1-eta2).gt.0.5e-1_rk )
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
                 end if
              end do
              r_change = 0.0_rk
              do j = 1, 3
                 r_change = r_change + rcoordinate(globalNM(i,3*j))*phii_1d(eta1,j)
              end do


              write(14, '(f9.3,2es15.7)') angle_c_degree, r_change, time 
              write(*,*) 'r_change for velocity =', r_change, 'element:', i
              ! exit   !?not strict

           end if   !u change element


           !grad(T) direction change
           if( gradT(1) * gradT(2) .lt. 0.0_rk .and. (i.ne.top_element .and. i.ne.CL_element) ) then  
              eta1 = 0.0_rk
              eta2 = 1.0_rk
              do while ( abs(eta1-eta2).gt.0.5e-1_rk )
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
                 end if
              end do
              r_change = 0.0_rk
              do j = 1, 3
                 r_change = r_change + rcoordinate(globalNM(i,3*j))*phii_1d(eta1,j)
              end do



              write(18, '(f9.3,2es15.7)') angle_c_degree, r_change, time 
              write(*,*) 'r_change for gradT =', r_change, 'element:', i
              ! exit   !?not strict


           end if   !gradT change element



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
  end if
  
  close(30)
  close(14)
  close(18)


  !flux on nodes of free surface
  if(timestep.eq.1) then
     open(unit = 10, file = trim(folder)//'flux.dat', status = 'replace')
  else
     open(unit = 10, file = trim(folder)//'flux.dat', status = 'old', access = 'append')
  end if

  allocate(flux(NTN), flux1(NTN), Dflux(NTN))
  J0 = 0.0_rk
  flux(:) = 0.0_rk
  Dflux(:) = 0.0_rk

  if(timestep.le.5 .or. mod(timestep,graph_step).eq.0) then   !write data every several timestep

     write(10, '(A)') 'variables = "r", "J"'!, "Deegan_flux"'
     write(10, '(A,f6.3,A)') 'Zone T = "angle =', angle_c_degree, '"'

     do i = 1, NTN
        if( ( BCflagN(i,3).eq.1 .or. BCflagN(i,3).eq.3 ) .and. PN(i).eq.1) then

           !(r,z,c)eta
           if(i.eq.top_node) then      !use the below element
              n1 = i
              do j = i-1, 1, -1
                 if(VN(j).eq.2) then
                    n2 = j
                    do k = j-1, 1, -1
                       if(VN(k).eq.2) then
                          n3 = k
                          exit
                       end if
                    end do
                    exit
                 end if
              end do
              reta = 3.0_rk*rcoordinate(n1) - 4.0_rk*rcoordinate(n2) + rcoordinate(n3)
              zeta = 3.0_rk*zcoordinate(n1) - 4.0_rk*zcoordinate(n2) + zcoordinate(n3)
              ceta = 3.0_rk*csol(n1) - 4.0_rk*csol(n2) + csol(n3)

           else                 !use the above element
              n1 = i
              do j = i+1, NTN
                 if(BCflagN(j,3).eq.1) then
                    n2 = j
                    do k = j+1, NTN
                       if(BCflagN(k,3).eq.1) then
                          n3 = k
                          exit
                       end if
                    end do
                    exit
                 end if
              end do
              ! write(*,*) n1,n2,n3
              ! pause
              reta = -3.0_rk*rcoordinate(n1) + 4.0_rk*rcoordinate(n2) - rcoordinate(n3)
              zeta = -3.0_rk*zcoordinate(n1) + 4.0_rk*zcoordinate(n2) - zcoordinate(n3)
              ! ceta = -3.0_rk*csol(n1) + 4.0_rk*csol(n2) - csol(n3)
           end if
           ! !(r,z,c)si
           ! rsi = -3.0_rk*rcoordinate(i) + 4.0_rk*rcoordinate(i+1) - rcoordinate(i+2)
           ! zsi = -3.0_rk*zcoordinate(i) + 4.0_rk*zcoordinate(i+1) - zcoordinate(i+2)
           ! csi = -3.0_rk*csol(i) + 4.0_rk*csol(i+1) - csol(i+2)

           !flux
           ! flux(i) = 1.0_rk / ( rsi*zeta - reta*zsi ) * sqrt(reta**2 + zeta**2) *  (-csi)
           rdot = soldot( NOPP(i) + Nr )
           zdot = soldot( NOPP(i) + Nz )!??
           flux(i) = ( zeta*(usol(i)-rdot) - reta*(vsol(i)-zdot) )/ sqrt(reta**2+zeta**2)
           ! flux(i) = 1.0_rk / ( rsi*zeta - reta*zsi ) * sqrt(reta**2 + zeta**2) *  (-csi)
           if(i.eq.top_node) J0 = flux(i)

           Dflux(i) = (1.0_rk - rcoordinate(i)**2)**( -(0.5_rk-angle_c/pi) )

        end if
     end do


     do i = 1, NTN
        if( ( BCflagN(i,3).eq.1 .or. BCflagN(i,3).eq.3 ) .and. PN(i).eq.1) then 
           ! flux(i) = flux(i)/J0
           ! if(i.eq.1) cycle  !put infinity on hold

           write(10,'(4es15.7)')  rcoordinate(i), flux(i)!, Dflux(i)!, abs(flux(i)-Dflux(i))

        end if
     end do

  end if  !every 100 steps

  close(10)

  ! !flux on gausspoints of free surface
  ! allocate( flux_gp(Ng) )
  ! write(10, '(A)') 'variables = "r", "flux"'
  ! write(10, '(A,f6.3,A)') 'Zone T = "t=', time, '"'
  ! do i = 1, NTE
  !    if(BCflagE(i,3).ne.2) cycle

  !    call values_in_an_element(i,1)

  !    do k = 1, Ng
  !       flux_gp(k) = 1.0_rk / ( rsi_right(k,1)*zeta_right(k,1) - reta_right(k,1)*zsi_right(k,1) ) / &
  !            sqrt(reta_right(k,1)**2 + zeta_right(k,1)**2) * &
  !            ( -dcdsi(k,1) * ( reta_right(k,1)**2 + zeta_right(k,1)**2 ) + &
  !            dcdeta(k,1) * ( rsi_right(k,1)*reta_right(k,1) + zsi_right(k,1)*zeta_right(k,1) ) )
  !       write(10,'(2es15.7)')  rintfac_right(k,1), flux_gp(k)
  !    end do

  ! end do

  deallocate(flux, flux1, Dflux)


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
