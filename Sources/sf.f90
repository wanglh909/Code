
subroutine define_sf(m,i, sf, LNVar, LNOPP,id)
  use kind
  use data
  use Ldata
  use NOP_mod, only: gaussian_quadrature, gaussian_quadrature_1d

  implicit none


  integer(kind=ik), intent(in):: m, i, LNVar, LNOPP(9), id
  real(kind=rk), intent(out):: sf(LNVar)

  integer(kind=ik):: k, l, ipp   !no i (i is the i in main program)
  real(kind=rk):: intRsi_V(Ng,Ng), intReta_V(Ng,Ng), &
       intRu_V(Ng,Ng), intRv_V(Ng,Ng), intRt_V(Ng,Ng), intRp(Ng,Ng), intRc(Ng,Ng)
  real(kind=rk):: intRsi_S(Ng), intReta_S(Ng), intRu_S(Ng), intRv_S(Ng), intRt_S(Ng)
  real(kind=rk):: conv, flux
  
  intRsi_V(:,:) = 0.0_rk 
  intReta_V(:,:) = 0.0_rk 
  intRu_V(:,:) = 0.0_rk 
  intRv_V(:,:) = 0.0_rk 
  intRt_V(:,:) = 0.0_rk 
  intRp(:,:) = 0.0_rk 
  intRc(:,:) = 0.0_rk 
  intRsi_S(:) = 0.0_rk
  intReta_S(:) = 0.0_rk
  intRu_S(:) = 0.0_rk
  intRv_S(:) = 0.0_rk
  intRt_S(:) = 0.0_rk

  do k = 1, MDF( globalNM(m,i) )
     sf(LNOPP(i) + k-1) = 0.0_rk
  end do


  do k = 1, Ng, 1
     do l = 1, Ng, 1

        !Aterm(k,l,i)
        Aterm(k,l,id) = ( zeta(k,l,id)**2 + reta(k,l,id)**2 )*phisi(k,l,i) - &
             ( zeta(k,l,id)*zsi(k,l,id) + reta(k,l,id)*rsi(k,l,id) )*phieta(k,l,i)
        Bterm(k,l,id) = -( zsi(k,l,id)*zeta(k,l,id) + rsi(k,l,id)*reta(k,l,id) )*phisi(k,l,i) &
             + ( zsi(k,l,id)**2 + rsi(k,l,id)**2 )*phieta(k,l,i)

        
        if( BCflagN( globalNM(m,i), 3 ).ne.1 .and. BCflagN( globalNM(m,i), 3 ).ne.3 )  &    !KBC on free surface has no volume integral
        intRsi_V(k,l) = ( s_orth(k,l,id) + epss )*Aterm(k,l,id)/Jp(k,l,id) * Jpsign(k,l,id) &
             - eps1*phisi(k,l,i)*fsi_size(m)*log( rsi(k,l,id)**2 + zsi(k,l,id)**2 )

        intReta_V(k,l) = ( 1.0_rk/s_orth(k,l,id) + epss )*Bterm(k,l,id)/Jp(k,l,id) * Jpsign(k,l,id) &
             - eps2*phieta(k,l,i)*geta_size(m)*log( reta(k,l,id)**2 + zeta(k,l,id)**2 )

        if(s_mode.eq.0) then

if(VE(m).eq.0) then
   
   if(NStrans.eq.1) then
      intRu_V(k,l) = intRu_V(k,l) + Re*phi(k,l,i)*( &
           udotintfac(k,l,id) - rdotintfac(k,l,id)*urintfac(k,l,id) - &
           zdotintfac(k,l,id)*uzintfac(k,l,id) ) *rintfac(k,l,id)*abs(Jp(k,l,id))

      intRv_V(k,l) = intRv_V(k,l) + Re*phi(k,l,i)*( &
           vdotintfac(k,l,id) - rdotintfac(k,l,id)*vrintfac(k,l,id) - &
           zdotintfac(k,l,id)*vzintfac(k,l,id) ) *rintfac(k,l,id)*abs(Jp(k,l,id))
   end if

   if(Inert.eq.1) then
      intRu_V(k,l) = intRu_V(k,l) + Re*phi(k,l,i)*( uintfac(k,l,id)*urintfac(k,l,id) +  &
           vintfac(k,l,id)*uzintfac(k,l,id) ) *rintfac(k,l,id)*abs(Jp(k,l,id))

      intRv_V(k,l) = intRv_V(k,l) + Re*phi(k,l,i)*( uintfac(k,l,id)*vrintfac(k,l,id) + &
           vintfac(k,l,id)*vzintfac(k,l,id) ) *rintfac(k,l,id)*abs(Jp(k,l,id))
   end if

   if(Capil.eq.1) then
      intRu_V(k,l) = intRu_V(k,l) - Kdi*pintfac(k,l,id)*Oh*&
           ( phir(k,l,i,id) + phi(k,l,i)/rintfac(k,l,id) ) *rintfac(k,l,id)*abs(Jp(k,l,id))

      intRv_V(k,l) = intRv_V(k,l) - Kdi*pintfac(k,l,id)*Oh*phiz(k,l,i,id) *rintfac(k,l,id)*abs(Jp(k,l,id))
   end if

   if(Viscous.eq.1) then
      intRu_V(k,l) = intRu_V(k,l) + ( &
           Oh*( 2.0_rk*urintfac(k,l,id)*phir(k,l,i,id) + &
           ( uzintfac(k,l,id) + vrintfac(k,l,id) )*phiz(k,l,i,id) + &
           phi(k,l,i)*2.0_rk/rintfac(k,l,id)**2 *uintfac(k,l,id) ) )  *rintfac(k,l,id)*abs(Jp(k,l,id))

      intRv_V(k,l) = intRv_V(k,l) + ( &
           Oh*( ( uzintfac(k,l,id) + vrintfac(k,l,id) ) *phir(k,l,i,id) + &
           2.0_rk*vzintfac(k,l,id)*phiz(k,l,i,id) ) ) *rintfac(k,l,id)*abs(Jp(k,l,id))
   end if

   if(GravI.eq.1) then
      intRv_V(k,l) = intRv_V(k,l) + Re*phi(k,l,i)* Grav  *rintfac(k,l,id)*abs(Jp(k,l,id))
   end if
     
   ! intRu_V(k,l) = ( Re*phi(k,l,i)*( udotintfac(k,l,id) + ( uintfac(k,l,id) - rdotintfac(k,l,id) )*urintfac(k,l,id) +  &
   !      ( vintfac(k,l,id) - zdotintfac(k,l,id) )*uzintfac(k,l,id) )  - &
   !      Kdi*pintfac(k,l,id)*Oh*( phir(k,l,i,id) + phi(k,l,i)/rintfac(k,l,id) ) + &
   !      Oh*( 2.0_rk*urintfac(k,l,id)*phir(k,l,i,id) + ( uzintfac(k,l,id) + vrintfac(k,l,id) )*phiz(k,l,i,id) + &
   !      phi(k,l,i)*2.0_rk/rintfac(k,l,id)**2 *uintfac(k,l,id) ) )  *rintfac(k,l,id)*abs(Jp(k,l,id))

   ! intRv_V(k,l) = ( Re*phi(k,l,i)*( vdotintfac(k,l,id) + ( uintfac(k,l,id) - rdotintfac(k,l,id) )*vrintfac(k,l,id) + &
   !      ( vintfac(k,l,id) - zdotintfac(k,l,id) )*vzintfac(k,l,id) + Grav ) - Kdi*pintfac(k,l,id)*Oh*phiz(k,l,i,id) + &
   !      Oh*( ( uzintfac(k,l,id) + vrintfac(k,l,id) ) *phir(k,l,i,id) + 2.0_rk*vzintfac(k,l,id)*phiz(k,l,i,id) ) ) &
   !      *rintfac(k,l,id)*abs(Jp(k,l,id))

   !evaporation cooling on free surface, no volume integral
   if( BCflagN( globalNM(m,i), 3 ).ne.1 .and. BCflagN( globalNM(m,i), 3 ).ne.3 )  then
      if(Ttime.eq.1) &
           intRt_V(k,l) = intRt_V(k,l) + Pe*phi(k,l,i)*( Tdotintfac(k,l,id) &
           - rdotintfac(k,l,id) *Trintfac(k,l,id)  - zdotintfac(k,l,id) *Tzintfac(k,l,id) ) &
           *rintfac(k,l,id)*abs(Jp(k,l,id))
      if(Tconv.eq.1) &
           intRt_V(k,l) = intRt_V(k,l) + Pe*phi(k,l,i)*(  &
           uintfac(k,l,id) *Trintfac(k,l,id) + vintfac(k,l,id) *Tzintfac(k,l,id) ) &
           *rintfac(k,l,id)*abs(Jp(k,l,id))
      if(Tdiff.eq.1) &
           intRt_V(k,l) = intRt_V(k,l) + ( &
           phir(k,l,i,id)*Trintfac(k,l,id) + phiz(k,l,i,id)*Tzintfac(k,l,id) ) &
           *rintfac(k,l,id)*abs(Jp(k,l,id))
   end if

   if( PN( globalNM(m,i) ).eq.1 ) then
      intRp(k,l) = psi(k,l,i)* ( urintfac(k,l,id) + uintfac(k,l,id)/rintfac(k,l,id) + vzintfac(k,l,id) ) &
           *rintfac(k,l,id)*abs(Jp(k,l,id))
   end if

else if( VE(m).eq.1 ) then

   intRc(k,l) = ( crintfac(k,l,id)*phir(k,l,i,id) + czintfac(k,l,id)*phiz(k,l,i,id) ) &
        *rintfac(k,l,id)*abs(Jp(k,l,id))

else  !VE = 5
   
   if(TtimeS.eq.1) &
        intRt_V(k,l) = intRt_V(k,l) + phi(k,l,i)*( Tdotintfac(k,l,id) &
        - rdotintfac(k,l,id) *Trintfac(k,l,id)  - zdotintfac(k,l,id) *Tzintfac(k,l,id) ) &
        *rintfac(k,l,id)*abs(Jp(k,l,id))
   if(TdiffS.eq.1) &
        intRt_V(k,l) = intRt_V(k,l) + F0*( &
        phir(k,l,i,id)*Trintfac(k,l,id) + phiz(k,l,i,id)*Tzintfac(k,l,id) ) &
        *rintfac(k,l,id)*abs(Jp(k,l,id))

end if    !for VE=0

        end if    !for s_mode=0

     end do
  end do


  if( BCflagN( globalNM(m,i), 3 ).ne.1 .and. BCflagN( globalNM(m,i), 3 ).ne.3 ) &  !KBC on free surface has no volume integral
       sf(LNOPP(i)+Nr) = gaussian_quadrature(intRsi_V)
  sf(LNOPP(i)+Nz) = gaussian_quadrature(intReta_V)

  if(s_mode.eq.0) then

     if(VE(m).eq.0) then
        sf(LNOPP(i)+Nu) = gaussian_quadrature(intRu_V)
        sf(LNOPP(i)+Nv) = gaussian_quadrature(intRv_V)

        if( BCflagN( globalNM(m,i), 3 ).eq.0 ) then  !evaporation cooling on free surface, no volume integral
           if( BCflagN( globalNM(m,i),2 ) .eq. 0 ) then
              sf(LNOPP(i)+NT) = gaussian_quadrature(intRt_V)
           else  !base nodes
              sf(LNOPP(i)+NT) = gaussian_quadrature(intRt_V)/kR
           end if
        end if

        if( PN( globalNM(m,i) ).eq.1 )  sf(LNOPP(i)+Np) = gaussian_quadrature(intRp)

     else if( VE(m).eq.1 ) then
        sf(LNOPP(i) + MDF(globalNM(m,i)) -1 ) = gaussian_quadrature(intRc)
     else !VE = 5
        if( VN(globalNM(m,i)) .eq. 5 ) then
           sf(LNOPP(i)+NTs) = gaussian_quadrature(intRt_V)
        else  !VN = base: if in drop
           sf(LNOPP(i)+NT) = gaussian_quadrature(intRt_V)/F0
        end if
     end if   !for VE

  end if  !for s_mode=0


  !*******************************************adding SI to sf********************************************************

!axis
if(BCflagN( globalNM(m,i), 1 ) .eq.1 ) then
   ipp = i - 6  !phix_1d(k,ipp)
   do k = 1, Ng, 1    !three gausspoints
      intRsi_S(k) = phix_1d(k,ipp) * fsi_size(m) * log( rsi_down(k,id)**2 + zsi_down(k,id)**2 )
   end do
   sf(LNOPP(i)+Nr) = sf(LNOPP(i)+Nr)  -  M1*gaussian_quadrature_1d(intRsi_S)
end if

!base
if( ( BCflagE(m,2).eq.1 .or. BCflagE(m,2).eq.11 ) &
     .and. (BCflagN(globalNM(m,i),2).eq.1 .or. BCflagN(globalNM(m,i),2).eq.3) ) then

     if(BCflagE(m,2).eq.1) then 
        ipp = i  !phix_1d(k,ipp)
        
        do k = 1, Ng, 1    !three gausspoints
           intRsi_S(k) = phix_1d(k,ipp) * fsi_size(m) * log( rsi_left(k,id)**2 + zsi_left(k,id)**2 )
        end do

        if( ( (VE(m).eq.1 .and. rowNM(m).gt.NEV1(1)).or.&
             (VE(m).eq.5 .and. rowNM(m).gt.NEV1(1)+NEM+NEL) )&
             .or. base_mesh.eq.0 ) then !elements far from drop
           sf(LNOPP(i)+Nr) = sf(LNOPP(i)+Nr)  -  M1*gaussian_quadrature_1d(intRsi_S)
        else
           sf(LNOPP(i)+Nr) = gaussian_quadrature_1d(intRsi_S)
        end if

     else !BCflagE(m,2)=11
        if( .not. ( no_vapor.eq.1 .and. VN(globalNM(m,i)).eq.1 ) ) sf(LNOPP(i)+Nr) = 0.0_rk
     end if

end if

if( ( BCflagE(m,2).eq.2 .or. BCflagE(m,2).eq.21 ) &
     .and. (BCflagN(globalNM(m,i),2).eq.2 .or. BCflagN(globalNM(m,i),2) .eq.3) ) then

     if(BCflagE(m,2).eq.2) then
        ipp = i/3 + 1  !phix_1d(k,ipp)
        
        do k = 1, Ng, 1    !three gausspoints
           intReta_S(k) = phix_1d(k,ipp)* geta_size(m)* log( reta_left(k,id)**2 + zeta_left(k,id)**2 )
        end do

        if( base_mesh.eq.0 ) then 
           sf(LNOPP(i)+Nz) = sf(LNOPP(i)+Nz)  -  M2*gaussian_quadrature_1d(intReta_S)
        else
           sf(LNOPP(i)+Nz) = gaussian_quadrature_1d(intReta_S)
        end if

     else !BCflagE(m,2)=21
        sf(LNOPP(i)+Nz) = 0.0_rk
     end if

end if


!outer surface
if(BCflagN( globalNM(m,i), 4 ) .ne.0 ) then
   ipp = i/3  !phix_1d(k,ipp)
   do k = 1, Ng, 1    !three gausspoints

intReta_S(k) = phix_1d(k,ipp) * geta_size(m) * log( reta_right(k,id)**2 + zeta_right(k,id)**2 )
   end do
   sf(LNOPP(i)+Nz) = sf(LNOPP(i)+Nz)  -  M2*gaussian_quadrature_1d(intReta_S)
end if


!free suface, decouple elliptic mesh inside and outside the drop
if(mesh_decouple.eq.1) then

   if(BCflagN( globalNM(m,i), 3 ).eq.1 .or. BCflagN( globalNM(m,i), 3 ).eq.3) then
      if(VE(m).eq.0) then
         ipp = i/3  !phix_1d(k,ipp)
      ! else
      !    ipp = i/3 + 1  !phix_1d(k,ipp)
      ! end if

         do k = 1, Ng, 1    !three gausspoints
intReta_S(k) = phix_1d(k,ipp) * geta_size(m) * log( reta_right(k,id)**2 + zeta_right(k,id)**2 )
         end do
         !directly replace Reta_V with Reta_S, viz SI of elliptic mesh
         sf(LNOPP(i)+Nz) = gaussian_quadrature_1d(intReta_S)

      else !VE = 1
         sf(LNOPP(i)+Nz) = 0.0_rk
      end if
   end if

end if  !for mesh_decouple


!free surface
if(s_mode.eq.0) then

!free suface from drop, KBC part1 & traction BC & evaporative cooling part1
if( (BCflagN( globalNM(m,i), 3 ).eq.1 .or. BCflagN( globalNM(m,i), 3 ).eq.3) .and. VE(m).eq.0 ) then
   ipp = i/3  !phix_1d(k,ipp)
   do k = 1, Ng, 1    !three gausspoints

!KBC
if(uniflux) then  !uniflux:  flux( angle_c,rintfac_right(k,id) ) 
   !KBC2
   intRsi_S(k) = phi_1d(k,ipp)* flux( angle_c,rintfac_right(k,id) ) * &
        ( reta_right(k,id)**2 + zeta_right(k,id)**2 )**0.5_rk *rintfac_right(k,id)
   !evaporation cooling 2
   intRt_S(k) = REH * intRsi_S(k)
   if(true_uniflux.eq.1) intRsi_S(k) = intRsi_S(k) / flux( angle_c,rintfac_right(k,id) )
end if
!KBC1
intRsi_S(k) = intRsi_S(k) + KBCgroup* ( phi_1d(k,ipp)* &
     ( -zeta_right(k,id)*( uintfac_right(k,id) - rdotintfac_right(k,id) ) + &
     reta_right(k,id)*( vintfac_right(k,id) - zdotintfac_right(k,id) ) )*rintfac_right(k,id) )


!traction BC
intRu_S(k) = ( reta_right(k,id)*phix_1d(k,ipp) / ( reta_right(k,id)**2 + zeta_right(k,id)**2 ) &
     + phi_1d(k,ipp)/rintfac_right(k,id) )*rintfac_right(k,id)* &
     ( reta_right(k,id)**2 + zeta_right(k,id)**2 )**0.5_rk
  intRu_S(k) = intRu_S(k) - phi_1d(k,ipp) *beta *Teta_right(k,id)* &
       ( reta_right(k,id)**2 + zeta_right(k,id)**2 )**(-0.5_rk) *reta_right(k,id)*rintfac_right(k,id)

intRv_S(k) = zeta_right(k,id)*phix_1d(k,ipp) / &
     ( reta_right(k,id)**2 + zeta_right(k,id)**2 )**0.5_rk *rintfac_right(k,id)
  intRv_S(k) = intRv_S(k) - phi_1d(k,ipp) *beta *Teta_right(k,id)* &
       ( reta_right(k,id)**2 + zeta_right(k,id)**2 )**(-0.5_rk) *zeta_right(k,id)*rintfac_right(k,id)

!evaporative cooling 1
intRt_S(k) = intRt_S(k) - phi_1d(k,ipp) * ( &
     -dTdsi(k,id)* ( reta_right(k,id)**2 + zeta_right(k,id)**2 ) + &
     Teta_right(k,id)* ( rsi_right(k,id)*reta_right(k,id) + zsi_right(k,id)*zeta_right(k,id) ) &
      ) *rintfac_right(k,id) /Jp_right(k,id)

   end do

   sf(LNOPP(i)+Nr) = sf(LNOPP(i)+Nr) + gaussian_quadrature_1d(intRsi_S)
   sf(LNOPP(i)+Nu) = sf(LNOPP(i)+Nu) + gaussian_quadrature_1d(intRu_S)!/Ca
   sf(LNOPP(i)+Nv) = sf(LNOPP(i)+Nv) + gaussian_quadrature_1d(intRv_S)!/Ca
   sf(LNOPP(i)+NT) = sf(LNOPP(i)+NT) + gaussian_quadrature_1d(intRt_S)

end if


!free surface from vapor, KBC part2 & evaporation cooling part2
if(.not.uniflux) then
if( (BCflagN( globalNM(m,i), 3 ).eq.1 .or. BCflagN( globalNM(m,i), 3 ).eq.3) .and. VE(m).eq.1 ) then
   ipp = i/3 + 1  !phi_1d(k,ipp)
   do k = 1, Ng, 1    !three gausspoints

!KBC2

! !dcdeta=0 so delete terms with dcdeta
! intRsi_S(k) = phi_1d(k,ipp)*( &
!      -dcdsi(k,id)*( reta_right(k,id)**2 + zeta_right(k,id)**2 )  &
!      ) *rintfac_right(k,id)/Jp_right(k,id)

intRsi_S(k) = phi_1d(k,ipp)*( &
     -dcdsi(k,id)*( reta_right(k,id)**2 + zeta_right(k,id)**2 ) + &
     dcdeta(k,id)*( zsi_right(k,id)*zeta_right(k,id) + rsi_right(k,id)*reta_right(k,id) ) &
     ) *rintfac_right(k,id)/Jp_right(k,id)

!evaporative cooling 2
intRt_S(k) = REH* intRsi_S(k)

   end do
      sf(LNOPP(i)+Nr) = sf(LNOPP(i)+Nr) + gaussian_quadrature_1d(intRsi_S)
      sf(LNOPP(i)+NT) = sf(LNOPP(i)+NT) + gaussian_quadrature_1d(intRt_S)
end if  !free surface nodes
end if  !not uniflux


end if  !for s_mode=0




  !*****************************************************************************************************************

if( ( WFLAG(m).eq.1 .and. ( (i.eq.1) .or. (i.eq.2) .or. (i.eq.3) ) ) .or. &
     ( WFLAG(m).eq.2 .and. ( (i.eq.1) .or. (i.eq.2) ) ) ) then
   conv = sf(LNOPP(i)+Nr)
   sf(LNOPP(i)+Nr) = -sf(LNOPP(i)+Nz)
   sf(LNOPP(i)+Nz) = conv
else if( WFLAG(m).eq.2 .and. i.eq.3 ) then
   sf(LNOPP(i)+Nr) = sf(LNOPP(i)+Nr) - sf(LNOPP(i)+Nz)
   sf(LNOPP(i)+Nz) = 0.0_rk
end if

return

end subroutine define_sf
