
subroutine values_in_an_element(m,id)
  use kind
  use data
  use Ldata


  implicit none
  integer(kind=ik), intent(in):: m, id
  integer(kind=ik):: j, k, l, n, npp

  go to 120  !local initialization, make everything 0
  !give value to r1-r9 & z1-z9
  !define ulocal(9), vlocal(9), plocal(4)
121  do j = 1, 9, 1

     rlocal(j,id) = sol( NOPP( globalNM(m,j) ) + Nr )
     zlocal(j,id) = sol( NOPP( globalNM(m,j) ) + Nz )

     if(s_mode.eq.0) then

        if( VE(m).eq.0 ) then
           ulocal(j,id) = sol( NOPP( globalNM(m,j) ) + Nu )
           vlocal(j,id) = sol( NOPP( globalNM(m,j) ) + Nv )  
           Tlocal(j,id) = sol( NOPP( globalNM(m,j) ) + NT )
           if( PN( globalNM(m,j) ) .eq. 1 )  plocal(j,id) = sol( NOPP( globalNM(m,j) ) + Np ) 
           ! plocal(2, 4, 5, 6, 8) remain unchanged

           rdotlocal(j,id) = soldot( NOPP( globalNM(m,j) ) + Nr )
           zdotlocal(j,id) = soldot( NOPP( globalNM(m,j) ) + Nz )
           udotlocal(j,id) = soldot( NOPP( globalNM(m,j) ) + Nu )
           vdotlocal(j,id) = soldot( NOPP( globalNM(m,j) ) + Nv )
           Tdotlocal(j,id) = soldot( NOPP( globalNM(m,j) ) + NT )
        else if( VE(m).eq.1 ) then
           clocal(j,id) = sol( NOPP( globalNM(m,j) ) + MDF( globalNM(m,j) ) - 1 )
        else  !VE = 5
           rdotlocal(j,id) = soldot( NOPP( globalNM(m,j) ) + Nr )
           zdotlocal(j,id) = soldot( NOPP( globalNM(m,j) ) + Nz )
           if( VN(globalNM(m,j)).eq.5 ) then
              Tlocal(j,id) = sol( NOPP( globalNM(m,j) ) + NTs )
              Tdotlocal(j,id) = soldot( NOPP( globalNM(m,j) ) + NTs )
           else  !j is base node
              Tlocal(j,id) = sol( NOPP( globalNM(m,j) ) + NT )
              Tdotlocal(j,id) = soldot( NOPP( globalNM(m,j) ) + NT )
           end if
           !??
        end if
           
     end if !for s_mode=0

  end do


  rintfac(:,:,id) = 0.0_rk
  rsi(:,:,id) = 0.0_rk
  reta(:,:,id) = 0.0_rk
  zsi(:,:,id) = 0.0_rk
  zeta(:,:,id) = 0.0_rk
  do k = 1, 3, 1          !relate to a(3) (the value of si for gaussian_quadrature)
     do l = 1, 3, 1           !relate to a(3) (the value of eta for gaussian_quadrature)
        !define rsi, reta, zsi, zeta
        do n = 1, 9, 1
           rintfac(k,l,id) = rintfac(k,l,id) + rlocal(n,id)*phi(k,l,n)
           rsi(k,l,id) = rsi(k,l,id) + rlocal(n,id)*phisi(k,l,n)
           reta(k,l,id) = reta(k,l,id) + rlocal(n,id)*phieta(k,l,n)
           zsi(k,l,id) = zsi(k,l,id) + zlocal(n,id)*phisi(k,l,n)
           zeta(k,l,id) = zeta(k,l,id) + zlocal(n,id)*phieta(k,l,n)
        end do
        !define Jp(3,3)
        Jp(k,l,id) =  rsi(k,l,id)*zeta(k,l,id) - reta(k,l,id)*zsi(k,l,id) 
        Jpsign(k,l,id) = abs(Jp(k,l,id))/Jp(k,l,id)
        s_orth(k,l,id) = ( ( rsi(k,l,id)**2 + zsi(k,l,id)**2 )/( reta(k,l,id)**2 + zeta(k,l,id)**2 ) )**0.5_rk


        if(s_mode.eq.0) then
           !define phir, phiz
           do n = 1, 9, 1           !n is the notation of phi (phi1 - phi9)
              phir(k,l,n,id) = 1.0_rk/Jp(k,l,id) * &
                   ( zeta(k,l,id)*phisi(k,l,n) - zsi(k,l,id)*phieta(k,l,n) )
              phiz(k,l,n,id) = 1.0_rk/Jp(k,l,id) * &
                   ( -reta(k,l,id)*phisi(k,l,n) + rsi(k,l,id)*phieta(k,l,n) )
           end do
        end if !for s_mode=0

     end do
  end do


  if(s_mode.eq.0) then
     
     if( VE(m).eq.0 ) then

        !define uintfac, urintfac, uzintfac, vintfac, vrintfac, vzintfac
        uintfac(:,:,id) = 0.0_rk
        urintfac(:,:,id) = 0.0_rk
        uzintfac(:,:,id) = 0.0_rk
        vintfac(:,:,id) = 0.0_rk
        vrintfac(:,:,id) = 0.0_rk
        vzintfac(:,:,id) = 0.0_rk
        Trintfac(:,:,id) = 0.0_rk
        Tzintfac(:,:,id) = 0.0_rk
        pintfac(:,:,id) = 0.0_rk

        rdotintfac(:,:,id) = 0.0_rk
        zdotintfac(:,:,id) = 0.0_rk
        udotintfac(:,:,id) = 0.0_rk
        vdotintfac(:,:,id) = 0.0_rk
        Tdotintfac(:,:,id) = 0.0_rk

        do k = 1, 3
           do l = 1, 3

              do n = 1, 9, 1
                 uintfac(k,l,id) = uintfac(k,l,id) +  ulocal(n,id)*phi(k,l,n)
                 urintfac(k,l,id) = urintfac(k,l,id) +  ulocal(n,id)*phir(k,l,n,id)
                 uzintfac(k,l,id) = uzintfac(k,l,id) +  ulocal(n,id)*phiz(k,l,n,id)
                 vintfac(k,l,id) = vintfac(k,l,id) +  vlocal(n,id)*phi(k,l,n)
                 vrintfac(k,l,id) = vrintfac(k,l,id) +  vlocal(n,id)*phir(k,l,n,id)
                 vzintfac(k,l,id) = vzintfac(k,l,id) +  vlocal(n,id)*phiz(k,l,n,id)
                 Trintfac(k,l,id) = Trintfac(k,l,id) +  Tlocal(n,id)*phir(k,l,n,id)
                 Tzintfac(k,l,id) = Tzintfac(k,l,id) +  Tlocal(n,id)*phiz(k,l,n,id)
                 pintfac(k,l,id) =  pintfac(k,l,id) +  plocal(n,id)*psi(k,l,n)

                 rdotintfac(k,l,id) = rdotintfac(k,l,id) + rdotlocal(n,id)*phi(k,l,n)
                 zdotintfac(k,l,id) = zdotintfac(k,l,id) + zdotlocal(n,id)*phi(k,l,n)
                 udotintfac(k,l,id) = udotintfac(k,l,id) + udotlocal(n,id)*phi(k,l,n)
                 vdotintfac(k,l,id) = vdotintfac(k,l,id) + vdotlocal(n,id)*phi(k,l,n)
                 Tdotintfac(k,l,id) = Tdotintfac(k,l,id) + Tdotlocal(n,id)*phi(k,l,n)
              end do

           end do
        end do   !end loop for k&l, preparation for gaussian_quadrature

     else if( VE(m).eq.1 ) then

        crintfac(:,:,id) = 0.0_rk
        czintfac(:,:,id) = 0.0_rk
        do k = 1, 3
           do l = 1, 3

              do n = 1, 9, 1
                 crintfac(k,l,id) = crintfac(k,l,id) +  clocal(n,id)*phir(k,l,n,id)
                 czintfac(k,l,id) = czintfac(k,l,id) +  clocal(n,id)*phiz(k,l,n,id)
              end do
  
           end do
        end do   !end loop for k&l, preparation for gaussian_quadrature

     else  !VE = 5

        Trintfac(:,:,id) = 0.0_rk
        Tzintfac(:,:,id) = 0.0_rk

        rdotintfac(:,:,id) = 0.0_rk
        zdotintfac(:,:,id) = 0.0_rk
        Tdotintfac(:,:,id) = 0.0_rk

        do k = 1, 3
           do l = 1, 3

              do n = 1, 9, 1
                 Trintfac(k,l,id) = Trintfac(k,l,id) +  Tlocal(n,id)*phir(k,l,n,id)
                 Tzintfac(k,l,id) = Tzintfac(k,l,id) +  Tlocal(n,id)*phiz(k,l,n,id)

                 rdotintfac(k,l,id) = rdotintfac(k,l,id) + rdotlocal(n,id)*phi(k,l,n)
                 zdotintfac(k,l,id) = zdotintfac(k,l,id) + zdotlocal(n,id)*phi(k,l,n)
                 Tdotintfac(k,l,id) = Tdotintfac(k,l,id) + Tdotlocal(n,id)*phi(k,l,n)
              end do

           end do
        end do   !end loop for k&l, preparation for gaussian_quadrature

     end if  ! for VE

  end if !for s_mode=0


!---------------------------preparation for the SIs, define the temrs in SI---------------------------

  !other than free surface, SIs only exit in elliptic mesh equations
  rsi_down(:,id) = 0.0_rk
  zsi_down(:,id) = 0.0_rk
  rsi_left(:,id) = 0.0_rk
  zsi_left(:,id) = 0.0_rk
  reta_left(:,id) = 0.0_rk
  zeta_left(:,id) = 0.0_rk
  reta_right(:,id) = 0.0_rk
  zeta_right(:,id) = 0.0_rk
  Teta_right(:,id) = 0.0_rk

  rintfac_right(:,id) = 0.0_rk
  uintfac_right(:,id) = 0.0_rk
  vintfac_right(:,id) = 0.0_rk
  rdotintfac_right(:,id) = 0.0_rk
  zdotintfac_right(:,id) = 0.0_rk
  
  dTdsi(:,id) = 0.0_rk
  dcdsi(:,id) = 0.0_rk
  dTdeta(:,id) = 0.0_rk
  dcdeta(:,id) = 0.0_rk
  rsi_right(:,id) = 0.0_rk
  zsi_right(:,id) = 0.0_rk
  Jp_right(:,id) = 0.0_rk


!? change left to axis, change down to base, change right to free
  
  
  !base
  if( BCflagE(m,2).eq.1 ) then
     do k = 1, 3, 1  !the three gausspoints for SI
        do n = 1, 3, 1  !the summation of three terms, eg: rsi = sum( rlocal(7,8,9)*phix_1d(1,2,3) )
           npp = n  !only need r1, r2, r3 -->r(npp)
           rsi_left(k,id) = rsi_left(k,id) + rlocal(npp,id) * phix_1d(k,n)
           zsi_left(k,id) = zsi_left(k,id) + zlocal(npp,id) * phix_1d(k,n)
        end do  !end for n
     end do    !end for k
  end if
  if( BCflagE(m,2).eq.11 ) then
     do k = 1, 3, 1  !the three gausspoints for SI
        do n = 1, 3, 1  !the summation of three terms, eg: rsi = sum( rlocal(7,8,9)*phix_1d(1,2,3) )
           npp = n+6  !only need r7, r8, r9 -->r(npp)
           rsi_left(k,id) = rsi_left(k,id) + rlocal(npp,id) * phix_1d(k,n)
           zsi_left(k,id) = zsi_left(k,id) + zlocal(npp,id) * phix_1d(k,n)
        end do  !end for n
     end do    !end for k
  end if

  if( BCflagE(m,2).eq.2 ) then
     do k = 1, 3, 1  !the three gausspoints for SI
        do n = 1, 3, 1  !the summation of three terms, eg: rsi = sum( rlocal(7,8,9)*phix_1d(1,2,3) )
           npp = 3*n - 2  !only need r1, r4, r7 -->r(npp)
           reta_left(k,id) = reta_left(k,id) + rlocal(npp,id) * phix_1d(k,n)
           zeta_left(k,id) = zeta_left(k,id) + zlocal(npp,id) * phix_1d(k,n)
        end do  !end for n
     end do    !end for k
  end if
  if( BCflagE(m,2).eq.21 ) then
     do k = 1, 3, 1  !the three gausspoints for SI
        do n = 1, 3, 1  !the summation of three terms, eg: rsi = sum( rlocal(7,8,9)*phix_1d(1,2,3) )
           npp = 3*n  !only need r3, r6, r9 -->r(npp)
           reta_left(k,id) = reta_left(k,id) + rlocal(npp,id) * phix_1d(k,n)
           zeta_left(k,id) = zeta_left(k,id) + zlocal(npp,id) * phix_1d(k,n)
        end do  !end for n
     end do    !end for k
  end if

  !axis
  if( BCflagE(m,1).eq.1 ) then
     do k = 1, 3, 1  !the three gausspoints for SI
        do n = 1, 3, 1  !the summation of three terms, eg: rsi = sum( rlocal(7,8,9)*phix_1d(1,2,3) )
           npp = n + 6  !only need r7, r8, r9 -->r(npp)
           rsi_down(k,id) = rsi_down(k,id) + rlocal(npp,id)*phix_1d(k,n)
           zsi_down(k,id) = zsi_down(k,id) + zlocal(npp,id)*phix_1d(k,n)
        end do  !end for n
     end do    !end for k
  end if

  !outer surface
  if( BCflagE(m,4).eq.1 ) then
     do k = 1, 3, 1  !the three gausspoints for SI
        do n = 1, 3, 1  !the summation of three terms, eg: rsi = sum( rlocal(7,8,9)*phix_1d(1,2,3) )
           npp = 3*n  !only need r3, r6, r9 -->r(npp)
           reta_right(k,id) = reta_right(k,id) + rlocal(npp,id) * phix_1d(k,n)
           zeta_right(k,id) = zeta_right(k,id) + zlocal(npp,id) * phix_1d(k,n)
        end do  !end for n
     end do    !end for k
  end if
  
  
  !free surface  
  if( BCflagE(m,3).eq.1 )  then
     do k = 1, Ng, 1    !three gausspoints
        do n = 1, 3, 1 !the summation of three terms, eg: rsi = sum( rlocal(7,8,9)*phix_1d(1,2,3) )
           npp = 3*n  !only need r3, r6, r9 -->r(npp)
           reta_right(k,id) = reta_right(k,id) + rlocal(npp,id) * phix_1d(k,n)
           zeta_right(k,id) = zeta_right(k,id) + zlocal(npp,id) * phix_1d(k,n)

           if(s_mode.eq.0) then
              rintfac_right(k,id) = rintfac_right(k,id) + rlocal(npp,id) * phi_1d(k,n)

              uintfac_right(k,id) = uintfac_right(k,id) + ulocal(npp,id) * phi_1d(k,n)
              vintfac_right(k,id) = vintfac_right(k,id) + vlocal(npp,id) * phi_1d(k,n)
              rdotintfac_right(k,id) = rdotintfac_right(k,id) + rdotlocal(npp,id) * phi_1d(k,n)
              zdotintfac_right(k,id) = zdotintfac_right(k,id) + zdotlocal(npp,id) * phi_1d(k,n)

              Teta_right(k,id) = Teta_right(k,id) + Tlocal(npp,id) * phix_1d(k,n)

           end if  !for s_mode=0
        end do  !end for n
        
        if(s_mode.eq.0) then
           do n = 1, 9, 1 !the summation of nine terms, eg: csi = sum( clocal(1-9)*phisi_1d(1-9) )
              dTdsi(k,id) = dTdsi(k,id) + Tlocal(n,id) * phisi1_1d(k,n)
              rsi_right(k,id) = rsi_right(k,id) + rlocal(n,id) * phisi1_1d(k,n)
              zsi_right(k,id) = zsi_right(k,id) + zlocal(n,id) * phisi1_1d(k,n)
           end do  !end for n

           Jp_right(k,id) = rsi_right(k,id)*zeta_right(k,id) - reta_right(k,id)*zsi_right(k,id)

        end if  !for s_mode=0

     end do    !end for k
  end if

  if( BCflagE(m,3).eq.2 ) then
     do k = 1, Ng, 1    !three gausspoints

        do n = 1, 3, 1 !the summation of three terms, eg: rsi = sum( rlocal(7,8,9)*phix_1d(1,2,3) )
           npp = 3*n-2  !only need r1, r4, r7 -->r(npp)
           reta_right(k,id) = reta_right(k,id) + rlocal(npp,id) * phix_1d(k,n)
           zeta_right(k,id) = zeta_right(k,id) + zlocal(npp,id) * phix_1d(k,n)

           if(s_mode.eq.0) then
              rintfac_right(k,id) = rintfac_right(k,id) + rlocal(npp,id) * phi_1d(k,n)
              dcdeta(k,id) = dcdeta(k,id) + clocal(npp,id) * phix_1d(k,n)
           end if  !for s_mode=0

        end do  !end for n

        if(s_mode.eq.0) then
           do n = 1, 9, 1 !the summation of nine terms, eg: csi = sum( clocal(1-9)*phisi_1d(1-9) )
              dcdsi(k,id) = dcdsi(k,id) + clocal(n,id) * phisi0_1d(k,n)
              rsi_right(k,id) = rsi_right(k,id) + rlocal(n,id) * phisi0_1d(k,n)
              zsi_right(k,id) = zsi_right(k,id) + zlocal(n,id) * phisi0_1d(k,n)
           end do  !end for n

           Jp_right(k,id) = rsi_right(k,id)*zeta_right(k,id) - reta_right(k,id)*zsi_right(k,id)

        end if  !for s_mode=0

     end do    !end for k
  end if
go to 122


! if(m.eq.44 .and. s_mode.eq.0) then
!    write(*,*) Trintfac(:,:,id), Tzintfac(:,:,id), urintfac(:,:,id), uzintfac(:,:,id)
!    pause
! end if


120  rlocal(:,id) = 0.0_rk
  zlocal(:,id) = 0.0_rk
  ulocal(:,id) = 0.0_rk
  vlocal(:,id) = 0.0_rk
  Tlocal(:,id) = 0.0_rk
  plocal(:,id) = 0.0_rk
  rdotlocal(:,id) = 0.0_rk
  zdotlocal(:,id) = 0.0_rk
  udotlocal(:,id) = 0.0_rk
  vdotlocal(:,id) = 0.0_rk
  Tdotlocal(:,id) = 0.0_rk
  clocal(:,id) = 0.0_rk
  rintfac(:,:,id) = 0.0_rk
  rsi(:,:,id) = 0.0_rk
  reta(:,:,id) = 0.0_rk
  zsi(:,:,id) = 0.0_rk
  zeta(:,:,id) = 0.0_rk
  Jp(:,:,id) =  0.0_rk 
  Jpsign(:,:,id) = 0.0_rk
  s_orth(:,:,id) = 0.0_rk
  phir(:,:,:,id) = 0.0_rk
  phiz(:,:,:,id) = 0.0_rk
  uintfac(:,:,id) = 0.0_rk
  urintfac(:,:,id) = 0.0_rk
  uzintfac(:,:,id) = 0.0_rk
  vintfac(:,:,id) = 0.0_rk
  vrintfac(:,:,id) = 0.0_rk
  vzintfac(:,:,id) = 0.0_rk
  Trintfac(:,:,id) = 0.0_rk
  Tzintfac(:,:,id) = 0.0_rk
  pintfac(:,:,id) = 0.0_rk
  rdotintfac(:,:,id) = 0.0_rk
  zdotintfac(:,:,id) = 0.0_rk
  udotintfac(:,:,id) = 0.0_rk
  vdotintfac(:,:,id) = 0.0_rk
  Tdotintfac(:,:,id) = 0.0_rk
  crintfac(:,:,id) = 0.0_rk
  czintfac(:,:,id) = 0.0_rk
  rsi_down(:,id) = 0.0_rk
  zsi_down(:,id) = 0.0_rk
  rsi_left(:,id) = 0.0_rk
  zsi_left(:,id) = 0.0_rk
  reta_left(:,id) = 0.0_rk
  zeta_left(:,id) = 0.0_rk
  reta_right(:,id) = 0.0_rk
  zeta_right(:,id) = 0.0_rk
  Teta_right(:,id) = 0.0_rk
  rintfac_right(:,id) = 0.0_rk
  uintfac_right(:,id) = 0.0_rk
  vintfac_right(:,id) = 0.0_rk
  rdotintfac_right(:,id) = 0.0_rk
  zdotintfac_right(:,id) = 0.0_rk
  dTdsi(:,id) = 0.0_rk
  dcdsi(:,id) = 0.0_rk
  dTdeta(:,id) = 0.0_rk
  dcdeta(:,id) = 0.0_rk
  rsi_right(:,id) = 0.0_rk
  zsi_right(:,id) = 0.0_rk
  Jp_right(:,id) = 0.0_rk
go to 121

122  return
end subroutine values_in_an_element
