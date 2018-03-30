subroutine drop_volume(vol1, vol2)
  use kind
  use data
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



