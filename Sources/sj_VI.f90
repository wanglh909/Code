subroutine VI_in_sj(m,i,j, sj, LNVar, LNOPP,id)
  use kind
  use data
  use Ldata
  use NOP_mod, only: gaussian_quadrature

  implicit none

  integer(kind=ik), intent(in):: m,i,j, LNVar, LNOPP(9), id
  real(kind=rk), intent(out):: sj(LNVar, LNVar)

  integer(kind=ik):: k, l
  real(kind=rk):: intRsi_r_V(Ng,Ng), intRsi_z_V(Ng,Ng), intReta_r_V(Ng,Ng), intReta_z_V(Ng,Ng)
  real(kind=rk):: intRu_r_V(Ng,Ng), intRu_z_V(Ng,Ng), intRu_u(Ng,Ng), intRuv(Ng,Ng), intRu_p(Ng,Ng)
  real(kind=rk):: intRv_r_V(Ng,Ng), intRv_z_V(Ng,Ng), intRv_u(Ng,Ng), intRvv(Ng,Ng),intRv_p(Ng,Ng) 
  real(kind=rk):: intRp_r(Ng,Ng), intRp_z(Ng,Ng), intRp_u(Ng,Ng), intRp_v(Ng,Ng)
  real(kind=rk):: intRc_r(Ng,Ng), intRc_z(Ng,Ng), intRc_c(Ng,Ng)
  real(kind=rk):: intRt_r_V(Ng,Ng), intRt_z_V(Ng,Ng), intRt_u(Ng,Ng), intRtv(Ng,Ng), intRt_T(Ng,Ng)
  intRsi_r_V(:,:) = 0.0_rk
  intRsi_z_V(:,:) = 0.0_rk
  intReta_r_V(:,:) = 0.0_rk
  intReta_z_V(:,:) = 0.0_rk
  intRu_r_V(:,:) = 0.0_rk
  intRu_z_V(:,:) = 0.0_rk
  intRu_u(:,:) = 0.0_rk
  intRuv(:,:) = 0.0_rk
  intRu_p(:,:) = 0.0_rk
  intRv_r_V(:,:) = 0.0_rk
  intRv_z_V(:,:) = 0.0_rk
  intRv_u(:,:) = 0.0_rk
  intRvv(:,:) = 0.0_rk
  intRv_p(:,:) = 0.0_rk
  intRp_r(:,:) = 0.0_rk
  intRp_z(:,:) = 0.0_rk
  intRp_u(:,:) = 0.0_rk
  intRp_v(:,:) = 0.0_rk
  intRc_r(:,:) = 0.0_rk
  intRc_z(:,:) = 0.0_rk
  intRc_c(:,:) = 0.0_rk
  intRt_r_V(:,:) = 0.0_rk
  intRt_z_V(:,:) = 0.0_rk
  intRt_u(:,:) = 0.0_rk
  intRtv(:,:) = 0.0_rk
  intRt_T(:,:) = 0.0_rk


  do k = 1, MDF( globalNM(m,i) ) 
     do l = 1, MDF( globalNM(m,j) )
        sj(LNOPP(i) + k-1, LNOPP(j) + l-1) = 0.0_rk
     end do
  end do

  !define integrand(k,l,p,q):
  do k = 1, Ng, 1          !relate to a(Ng) (the value of si for gaussian_quadrature)
     do l = 1, Ng, 1           !relate to a(Ng) (the value of eta for gaussian_quadrature)

        !intRsi_u,v,p = intReta_u,v,p = 0.0_rk
        
        !KBC on free surface, no volume integral
        if( BCflagN( globalNM(m,i), 3 ).ne.1 .and. BCflagN( globalNM(m,i), 3 ).ne.3 ) then

intRsi_r_V(k,l) = ( s_orth_r(k,l,id) * Aterm(k,l,id) / Jp(k,l,id)  +  &
     ( s_orth(k,l,id) + epss ) * Aterm_r(k,l,id) / Jp(k,l,id)  &
     -  ( s_orth(k,l,id) + epss ) * Aterm(k,l,id) / Jp(k,l,id)**2 * Jp_r(k,l,id) ) * Jpsign(k,l,id)  &
     -  eps1 * phisi(k,l,i) * fsi_size(m) / ( rsi(k,l,id)**2 + zsi(k,l,id)**2 ) &
     * 2.0_rk*rsi(k,l,id) * phisi(k,l,j)

intRsi_z_V(k,l) = ( s_orth_z(k,l,id) * Aterm(k,l,id) / Jp(k,l,id)  +  &
     ( s_orth(k,l,id) + epss ) * Aterm_z(k,l,id) / Jp(k,l,id)  &
     -  ( s_orth(k,l,id) + epss ) * Aterm(k,l,id) / (Jp(k,l,id)**2) * Jp_z(k,l,id) ) * Jpsign(k,l,id)  &
     -  eps1 * phisi(k,l,i) * fsi_size(m) / ( rsi(k,l,id)**2 + zsi(k,l,id)**2 ) &
     * 2.0_rk*zsi(k,l,id) * phisi(k,l,j) 

        end if

intReta_r_V(k,l) = (-1.0_rk/(s_orth(k,l,id)**2) * s_orth_r(k,l,id) * Bterm(k,l,id) / Jp(k,l,id)  &
     +  ( 1.0_rk/s_orth(k,l,id) + epss ) * Bterm_r(k,l,id) / Jp(k,l,id)  &
     -  ( 1.0_rk/s_orth(k,l,id) + epss ) * Bterm(k,l,id) / (Jp(k,l,id)**2) * Jp_r(k,l,id) ) * Jpsign(k,l,id)  &
     -  eps2 * phieta(k,l,i) * geta_size(m) / ( reta(k,l,id)**2 + zeta(k,l,id)**2 ) &
     * 2.0_rk*reta(k,l,id) * phieta(k,l,j)

intReta_z_V(k,l) =  ( -1.0_rk/s_orth(k,l,id)**2 * s_orth_z(k,l,id) * Bterm(k,l,id) / Jp(k,l,id)  &
     +  ( 1.0_rk/s_orth(k,l,id) + epss ) * Bterm_z(k,l,id) / Jp(k,l,id)  &
     -  ( 1.0_rk/s_orth(k,l,id) + epss ) * Bterm(k,l,id) / (Jp(k,l,id)**2) * Jp_z(k,l,id) ) * Jpsign(k,l,id)  &
     -  eps2 * phieta(k,l,i) * geta_size(m) / ( reta(k,l,id)**2 + zeta(k,l,id)**2 ) &
     * 2.0_rk*zeta(k,l,id) * phieta(k,l,j)


        if(s_mode.eq.0) then 
        if(VE(m).eq.0) then
!          Ruu, Ruv, Rup, Rvu, Rvv, Rvp
!NStran:     1    0    0    0    1    0 
!Inert:      1    1    0    1    1    0  
!Capil:      0    0    1    0    0    1 
!Viscous:    1    1    0    1    1    0 

           if(NStrans.eq.1) then
intRu_r_V(k,l) = intRu_r_V(k,l) + Re*phi(k,l,i)* ( -CTJ*phi(k,l,j)/dt*urintfac(k,l,id) - &
     rdotintfac(k,l,id)*urintfac_r(k,l,id) - zdotintfac(k,l,id)*uzintfac_r(k,l,id) ) &
     *rintfac(k,l,id)*abs(Jp(k,l,id)) + &
     
      Re*phi(k,l,i)*( udotintfac(k,l,id) - rdotintfac(k,l,id)*urintfac(k,l,id) -  &
      zdotintfac(k,l,id)*uzintfac(k,l,id) ) *rJp_r(k,l,id) 

intRu_z_V(k,l) = intRu_z_V(k,l) +  Re*phi(k,l,i)* ( -rdotintfac(k,l,id)*urintfac_z(k,l,id) - &
     CTJ*phi(k,l,j)/dt*uzintfac(k,l,id) - zdotintfac(k,l,id)*uzintfac_z(k,l,id) )  &
     *rintfac(k,l,id)*abs(Jp(k,l,id)) + &
     
     Re*phi(k,l,i)*( udotintfac(k,l,id) - rdotintfac(k,l,id)*urintfac(k,l,id) -  &
     zdotintfac(k,l,id)*uzintfac(k,l,id) ) *rintfac(k,l,id)* Jp_z(k,l,id)*Jpsign(k,l,id)


intRu_u(k,l) = intRu_u(k,l) + Re*phi(k,l,i)* ( CTJ*phi(k,l,j)/dt -  &
     rdotintfac(k,l,id)*phir(k,l,j,id) - zdotintfac(k,l,id)*phiz(k,l,j,id) )  &
     *rintfac(k,l,id)*abs(Jp(k,l,id))  

!intRuv(k,l) = intRuv(k,l) + 0 
!intRu_p(k,l) = intRu_p(k,l) + 0

intRv_r_V(k,l) = intRv_r_V(k,l) + Re*phi(k,l,i)* ( -CTJ*phi(k,l,j)/dt*vrintfac(k,l,id) - &
     rdotintfac(k,l,id)*vrintfac_r(k,l,id) - zdotintfac(k,l,id)*vzintfac_r(k,l,id) ) &
     *rintfac(k,l,id)*abs(Jp(k,l,id)) + &
     
      Re*phi(k,l,i)*( vdotintfac(k,l,id) - rdotintfac(k,l,id)*vrintfac(k,l,id) - &
      zdotintfac(k,l,id)*vzintfac(k,l,id) ) *rJp_r(k,l,id) 


intRv_z_V(k,l) = intRv_z_V(k,l) + Re*phi(k,l,i)* ( -rdotintfac(k,l,id)*vrintfac_z(k,l,id) - &
     CTJ*phi(k,l,j)/dt*vzintfac(k,l,id) - zdotintfac(k,l,id)*vzintfac_z(k,l,id) ) &
     *rintfac(k,l,id)*abs(Jp(k,l,id)) + &
     
      Re*phi(k,l,i)*( vdotintfac(k,l,id) - rdotintfac(k,l,id)*vrintfac(k,l,id) - &
     zdotintfac(k,l,id)*vzintfac(k,l,id) ) *rintfac(k,l,id)* Jp_z(k,l,id)*Jpsign(k,l,id)

!intRv_u(k,l) = intRv_u(k,l) + 0

intRvv(k,l) = intRvv(k,l) + Re*phi(k,l,i)*( CTJ*phi(k,l,j)/dt - rdotintfac(k,l,id)*phir(k,l,j,id) - &
     zdotintfac(k,l,id)*phiz(k,l,j,id) ) *rintfac(k,l,id)*abs(Jp(k,l,id)) 

!intRv_p(k,l) = intRv_p(k,l) + 0

           end if

           if(Inert.eq.1) then
intRu_r_V(k,l) = intRu_r_V(k,l) + Re*phi(k,l,i)* ( uintfac(k,l,id)*urintfac_r(k,l,id) + &
      vintfac(k,l,id)*uzintfac_r(k,l,id) ) *rintfac(k,l,id)*abs(Jp(k,l,id)) + &
     
      Re*phi(k,l,i)*(  uintfac(k,l,id)*urintfac(k,l,id) +  vintfac(k,l,id)*uzintfac(k,l,id) ) *rJp_r(k,l,id) 


intRu_z_V(k,l) = intRu_z_V(k,l) +  Re*phi(k,l,i)* ( uintfac(k,l,id)*urintfac_z(k,l,id) + &
     vintfac(k,l,id)*uzintfac_z(k,l,id) ) *rintfac(k,l,id)*abs(Jp(k,l,id)) + &
     
     Re*phi(k,l,i)*( uintfac(k,l,id) *urintfac(k,l,id) + vintfac(k,l,id) *uzintfac(k,l,id) ) *&
     rintfac(k,l,id)* Jp_z(k,l,id)*Jpsign(k,l,id)


intRu_u(k,l) = intRu_u(k,l) + Re*phi(k,l,i)* ( phi(k,l,j)*urintfac(k,l,id) +  &
     uintfac(k,l,id)*phir(k,l,j,id) + vintfac(k,l,id)*phiz(k,l,j,id) ) &
     *rintfac(k,l,id)*abs(Jp(k,l,id))  

intRuv(k,l) = intRuv(k,l) +  Re*phi(k,l,i)*phi(k,l,j)*uzintfac(k,l,id) *rintfac(k,l,id)*abs(Jp(k,l,id)) 

!intRu_p(k,l) = intRu_p(k,l) + 0


intRv_r_V(k,l) = intRv_r_V(k,l) + Re*phi(k,l,i)* ( &
     uintfac(k,l,id)*vrintfac_r(k,l,id) + vintfac(k,l,id)*vzintfac_r(k,l,id) ) &
     *rintfac(k,l,id)*abs(Jp(k,l,id)) + &
     
     Re*phi(k,l,i)*( uintfac(k,l,id)*vrintfac(k,l,id) + vintfac(k,l,id)*vzintfac(k,l,id) ) *rJp_r(k,l,id) 


intRv_z_V(k,l) = intRv_z_V(k,l) + Re*phi(k,l,i)* ( uintfac(k,l,id)*vrintfac_z(k,l,id)  &
      + vintfac(k,l,id)*vzintfac_z(k,l,id) )  *rintfac(k,l,id)*abs(Jp(k,l,id)) + &
     
     Re*phi(k,l,i)*( uintfac(k,l,id) *vrintfac(k,l,id) + &
     vintfac(k,l,id)*vzintfac(k,l,id) ) *rintfac(k,l,id)* Jp_z(k,l,id)*Jpsign(k,l,id)


intRv_u(k,l) = intRv_u(k,l) + Re*phi(k,l,i)*phi(k,l,j)*vrintfac(k,l,id) *rintfac(k,l,id)*abs(Jp(k,l,id)) 

intRvv(k,l) = intRvv(k,l) + Re*phi(k,l,i)*( uintfac(k,l,id)*phir(k,l,j,id) + &
     phi(k,l,j)*vzintfac(k,l,id) + vintfac(k,l,id)*phiz(k,l,j,id) ) *rintfac(k,l,id)*abs(Jp(k,l,id)) 

!intRv_p(k,l) = intRv_p(k,l) + 0

           end if

           if(Capil.eq.1) then
intRu_r_V(k,l) = intRu_r_V(k,l) + ( -Kdi*pintfac(k,l,id)*Oh* &
     ( phir_r(k,l,i,id) - phi(k,l,i)/rintfac(k,l,id)**2 *phi(k,l,j) ) ) *&
     rintfac(k,l,id)*abs(Jp(k,l,id))  &
     
     - Kdi*pintfac(k,l,id)*Oh*( phir(k,l,i,id) + phi(k,l,i)/rintfac(k,l,id) ) *rJp_r(k,l,id) 


intRu_z_V(k,l) = intRu_z_V(k,l) - Kdi*pintfac(k,l,id)*Oh*phir_z(k,l,i,id) *&
     rintfac(k,l,id)*abs(Jp(k,l,id))  &
     
      -Kdi*pintfac(k,l,id)*Oh*( phir(k,l,i,id) + phi(k,l,i)/rintfac(k,l,id) ) *&
      rintfac(k,l,id)* Jp_z(k,l,id)*Jpsign(k,l,id)

!intRu_u(k,l) = intRu_u(k,l) + 0
!intRuv(k,l) = intRuv(k,l) + 0

intRu_p(k,l) = intRu_p(k,l) - Kdi*psi(k,l,j)*Oh*&
     ( phir(k,l,i,id) + phi(k,l,i)/rintfac(k,l,id) ) *rintfac(k,l,id)*abs(Jp(k,l,id))


intRv_r_V(k,l) = intRv_r_V(k,l) - &
     Kdi*pintfac(k,l,id)*Oh*phiz_r(k,l,i,id) *rintfac(k,l,id)*abs(Jp(k,l,id)) - &
     
     Kdi*pintfac(k,l,id)*Oh*phiz(k,l,i,id) *rJp_r(k,l,id) 


intRv_z_V(k,l) = intRv_z_V(k,l) - &
     Kdi*pintfac(k,l,id)*Oh*phiz_z(k,l,i,id) *rintfac(k,l,id)*abs(Jp(k,l,id)) - &
     
      Kdi*pintfac(k,l,id)*Oh*phiz(k,l,i,id) *rintfac(k,l,id)* Jp_z(k,l,id)*Jpsign(k,l,id)

!intRv_u(k,l) = intRv_u(k,l) + 0
!intRvv(k,l) = intRvv(k,l) + 0

intRv_p(k,l) = intRv_p(k,l) - Kdi*psi(k,l,j)*Oh*phiz(k,l,i,id) *rintfac(k,l,id)*abs(Jp(k,l,id))  

           end if

           if(Viscous.eq.1) then
intRu_r_V(k,l) = intRu_r_V(k,l) + ( &
     Oh*( 2.0_rk*urintfac_r(k,l,id)*phir(k,l,i,id) + 2.0_rk*urintfac(k,l,id)*phir_r(k,l,i,id) + &
     ( uzintfac_r(k,l,id) + vrintfac_r(k,l,id) )*phiz(k,l,i,id) + &
     ( uzintfac(k,l,id) + vrintfac(k,l,id) )*phiz_r(k,l,i,id) - &
     4.0_rk*phi(k,l,i)*phi(k,l,j)/rintfac(k,l,id)**3 *uintfac(k,l,id) ) ) *rintfac(k,l,id)*abs(Jp(k,l,id)) + &
     
     ( Oh*( 2.0_rk*urintfac(k,l,id)*phir(k,l,i,id) + ( uzintfac(k,l,id) + vrintfac(k,l,id) )*phiz(k,l,i,id) + &
     phi(k,l,i)*2.0_rk/rintfac(k,l,id)**2 *uintfac(k,l,id) ) ) *rJp_r(k,l,id) 


intRu_z_V(k,l) = intRu_z_V(k,l) + ( &
     Oh*( 2.0_rk*urintfac_z(k,l,id)*phir(k,l,i,id) + 2.0_rk*urintfac(k,l,id)*phir_z(k,l,i,id) + &
     ( uzintfac_z(k,l,id) + vrintfac_z(k,l,id) )*phiz(k,l,i,id) + &
     ( uzintfac(k,l,id) + vrintfac(k,l,id) )*phiz_z(k,l,i,id) ) ) *rintfac(k,l,id)*abs(Jp(k,l,id)) + &
     
     ( Oh*( 2.0_rk*urintfac(k,l,id)*phir(k,l,i,id) + ( uzintfac(k,l,id) + vrintfac(k,l,id) )*phiz(k,l,i,id) + &
     phi(k,l,i)*2.0_rk/rintfac(k,l,id)**2 *uintfac(k,l,id) ) )  *rintfac(k,l,id)* Jp_z(k,l,id)*Jpsign(k,l,id)


intRu_u(k,l) = intRu_u(k,l) + Oh*( 2.0_rk*phir(k,l,i,id)*phir(k,l,j,id) + &
     phiz(k,l,i,id)*phiz(k,l,j,id) + 2.0_rk/rintfac(k,l,id)**2 *phi(k,l,i)*phi(k,l,j) )  &
     *rintfac(k,l,id)*abs(Jp(k,l,id))  

intRuv(k,l) = intRuv(k,l) + Oh*phir(k,l,j,id)*phiz(k,l,i,id) *rintfac(k,l,id)*abs(Jp(k,l,id)) 

!intRu_p(k,l) = intRu_p(k,l) + 0


intRv_r_V(k,l) = intRv_r_V(k,l) + Oh*( ( uzintfac_r(k,l,id) + vrintfac_r(k,l,id) )*phir(k,l,i,id) + &
     ( uzintfac(k,l,id) + vrintfac(k,l,id) )*phir_r(k,l,i,id) + &
     2.0_rk*vzintfac_r(k,l,id)*phiz(k,l,i,id) + 2.0_rk*vzintfac(k,l,id)*phiz_r(k,l,i,id) ) &
     *rintfac(k,l,id)*abs(Jp(k,l,id)) + &
     
     Oh*( ( uzintfac(k,l,id) + vrintfac(k,l,id) ) *phir(k,l,i,id) + 2.0_rk*vzintfac(k,l,id)*phiz(k,l,i,id) ) *rJp_r(k,l,id) 


intRv_z_V(k,l) = intRv_z_V(k,l) + Oh*( ( uzintfac_z(k,l,id) + vrintfac_z(k,l,id) )*phir(k,l,i,id) + &
     ( uzintfac(k,l,id) + vrintfac(k,l,id) )*phir_z(k,l,i,id) + &
     2.0_rk*vzintfac_z(k,l,id)*phiz(k,l,i,id) + 2.0_rk*vzintfac(k,l,id)*phiz_z(k,l,i,id) ) *rintfac(k,l,id)*abs(Jp(k,l,id)) + &
     
     Oh*( ( uzintfac(k,l,id) + vrintfac(k,l,id) ) *phir(k,l,i,id) + &
     2.0_rk*vzintfac(k,l,id)*phiz(k,l,i,id) ) *rintfac(k,l,id)* Jp_z(k,l,id)*Jpsign(k,l,id)


intRv_u(k,l) = intRv_u(k,l) + Oh*phiz(k,l,j,id)*phir(k,l,i,id) *rintfac(k,l,id)*abs(Jp(k,l,id)) 

intRvv(k,l) = intRvv(k,l) + Oh*( phir(k,l,i,id)*phir(k,l,j,id) + &
     2.0_rk*phiz(k,l,i,id)*phiz(k,l,j,id) ) *rintfac(k,l,id)*abs(Jp(k,l,id)) 

!intRv_p(k,l) = intRv_p(k,l) + 0

           end if
!j=2,4,5,6,8,  intRu_p(k,l) = 0.0_rk

           if(GravI.eq.1) then
intRv_r_V(k,l) = intRv_r_V(k,l) + Re*phi(k,l,i)* Grav *rJp_r(k,l,id) 

intRv_z_V(k,l) = intRv_z_V(k,l) + Re*phi(k,l,i)* Grav *rintfac(k,l,id)* Jp_z(k,l,id)*Jpsign(k,l,id)

!intRv_u(k,l) = intRv_u(k,l) + 0
!intRvv(k,l) = intRvv(k,l) + 0
!intRv_p(k,l) = intRv_p(k,l) + 0
           end if


        !evaporation cooling on free surface, no volume integral
        if( BCflagN( globalNM(m,i), 3 ).ne.1 .and. BCflagN( globalNM(m,i), 3 ).ne.3 ) then
!Rt_r_V
if(Ttime.eq.1) &
intRt_r_V(k,l) = intRt_r_V(k,l) + Pe*phi(k,l,i)* ( -CTJ*phi(k,l,j)/dt*Trintfac(k,l,id) - &
     rdotintfac(k,l,id)*Trintfac_r(k,l,id) - zdotintfac(k,l,id)*Tzintfac_r(k,l,id) ) &
     *rintfac(k,l,id)*abs(Jp(k,l,id)) + &
     
     Pe*phi(k,l,i)*( Tdotintfac(k,l,id) - rdotintfac(k,l,id) *Trintfac(k,l,id) - &
     zdotintfac(k,l,id) *Tzintfac(k,l,id) ) *rJp_r(k,l,id) 

if(Tconv.eq.1) &
intRt_r_V(k,l) = intRt_r_V(k,l) + Pe*phi(k,l,i)* ( &
     uintfac(k,l,id)*Trintfac_r(k,l,id) + vintfac(k,l,id)*Tzintfac_r(k,l,id) )&
     *rintfac(k,l,id)*abs(Jp(k,l,id)) + &
     
     Pe*phi(k,l,i)*( uintfac(k,l,id)*Trintfac(k,l,id) + vintfac(k,l,id)*Tzintfac(k,l,id) ) &
     *rJp_r(k,l,id) 

if(Tdiff.eq.1) &
intRt_r_V(k,l) = intRt_r_V(k,l) + ( &
     Trintfac_r(k,l,id)*phir(k,l,i,id) + Trintfac(k,l,id)*phir_r(k,l,i,id) + &
     Tzintfac_r(k,l,id)*phiz(k,l,i,id) + Tzintfac(k,l,id)*phiz_r(k,l,i,id) ) &
     *rintfac(k,l,id)*abs(Jp(k,l,id)) + &
     
     ( phir(k,l,i,id)*Trintfac(k,l,id) + phiz(k,l,i,id)*Tzintfac(k,l,id) ) *rJp_r(k,l,id) 

!intRt_z_V
if(Ttime.eq.1) &
intRt_z_V(k,l) = intRt_z_V(k,l) + Pe*phi(k,l,i)* ( -rdotintfac(k,l,id)*Trintfac_z(k,l,id) - &
     CTJ*phi(k,l,j)/dt*Tzintfac(k,l,id) - zdotintfac(k,l,id)*Tzintfac_z(k,l,id) ) &
     *rintfac(k,l,id)*abs(Jp(k,l,id)) + &
     
     Pe*phi(k,l,i)*( Tdotintfac(k,l,id) - &
     rdotintfac(k,l,id)*Trintfac(k,l,id) - zdotintfac(k,l,id)*Tzintfac(k,l,id) )  &
     *rintfac(k,l,id)* Jp_z(k,l,id)*Jpsign(k,l,id)

if(Tconv.eq.1) &
intRt_z_V(k,l) = intRt_z_V(k,l) + Pe*phi(k,l,i)* ( uintfac(k,l,id)*Trintfac_z(k,l,id) + &
     vintfac(k,l,id)*Tzintfac_z(k,l,id) ) *rintfac(k,l,id)*abs(Jp(k,l,id)) + &
     
     Pe*phi(k,l,i)*( uintfac(k,l,id)*Trintfac(k,l,id) + vintfac(k,l,id)*Tzintfac(k,l,id) ) &
     *rintfac(k,l,id)* Jp_z(k,l,id)*Jpsign(k,l,id)

if(Tdiff.eq.1) &
intRt_z_V(k,l) = intRt_z_V(k,l) + ( &
     Trintfac_z(k,l,id)*phir(k,l,i,id) + Trintfac(k,l,id)*phir_z(k,l,i,id) + &
     Tzintfac_z(k,l,id)*phiz(k,l,i,id) + Tzintfac(k,l,id)*phiz_z(k,l,i,id) ) &
     *rintfac(k,l,id)*abs(Jp(k,l,id)) + &
     
     ( phir(k,l,i,id)*Trintfac(k,l,id) + phiz(k,l,i,id)*Tzintfac(k,l,id) )  &
     *rintfac(k,l,id)* Jp_z(k,l,id)*Jpsign(k,l,id)


if(Tconv.eq.1) then
intRt_u(k,l) = Pe*phi(k,l,i)* phi(k,l,j)*Trintfac(k,l,id) *rintfac(k,l,id)*abs(Jp(k,l,id))  
intRtv(k,l) =  Pe*phi(k,l,i)* phi(k,l,j)*Tzintfac(k,l,id) *rintfac(k,l,id)*abs(Jp(k,l,id))  
end if

!Rt_T
if(Ttime.eq.1) &
intRt_T(k,l) = intRt_T(k,l) + Pe*phi(k,l,i)* ( CTJ*phi(k,l,j)/dt - &
     rdotintfac(k,l,id)*phir(k,l,j,id) - zdotintfac(k,l,id)*phiz(k,l,j,id) ) &
     *rintfac(k,l,id)*abs(Jp(k,l,id)) 
if(Tconv.eq.1) &
intRt_T(k,l) = intRt_T(k,l) + Pe*phi(k,l,i)* ( uintfac(k,l,id)*phir(k,l,j,id) + &
     vintfac(k,l,id)*phiz(k,l,j,id) ) *rintfac(k,l,id)*abs(Jp(k,l,id)) 
if(Tdiff.eq.1) &
intRt_T(k,l) = intRt_T(k,l) + ( phir(k,l,i,id)*phir(k,l,j,id) + phiz(k,l,i,id)*phiz(k,l,j,id) ) &
     *rintfac(k,l,id)*abs(Jp(k,l,id)) 

        end if

if( PN( globalNM(m,i) ).eq.1 ) then

intRp_r(k,l) = psi(k,l,i)* ( urintfac_r(k,l,id) - uintfac(k,l,id)/rintfac(k,l,id)**2 *phi(k,l,j) + vzintfac_r(k,l,id) ) &
     *rintfac(k,l,id)*abs(Jp(k,l,id)) + &
     
     psi(k,l,i)* ( urintfac(k,l,id) + uintfac(k,l,id)/rintfac(k,l,id) + vzintfac(k,l,id) ) *rJp_r(k,l,id)


intRp_z(k,l) = psi(k,l,i)* ( urintfac_z(k,l,id) + vzintfac_z(k,l,id) ) *rintfac(k,l,id)*abs(Jp(k,l,id)) + &
     
     psi(k,l,i)* ( urintfac(k,l,id) + uintfac(k,l,id)/rintfac(k,l,id) + vzintfac(k,l,id) ) *rintfac(k,l,id)* &
     Jp_z(k,l,id)*Jpsign(k,l,id)


intRp_u(k,l) = ( phi(k,l,j)/rintfac(k,l,id) + phir(k,l,j,id) )*psi(k,l,i) *rintfac(k,l,id)*abs(Jp(k,l,id))

intRp_v(k,l) = phiz(k,l,j,id)*psi(k,l,i) *rintfac(k,l,id)*abs(Jp(k,l,id)) 

end if

        else if ( VE(m).eq.1 ) then

intRc_r(k,l) = ( crintfac_r(k,l,id)*phir(k,l,i,id) + crintfac(k,l,id)*phir_r(k,l,i,id) + &
     czintfac_r(k,l,id)*phiz(k,l,i,id) + czintfac(k,l,id)*phiz_r(k,l,i,id) ) &
     *rintfac(k,l,id)*abs(Jp(k,l,id)) + &

     ( crintfac(k,l,id)*phir(k,l,i,id) + czintfac(k,l,id)*phiz(k,l,i,id) ) *rJp_r(k,l,id)

intRc_z(k,l) = ( crintfac_z(k,l,id)*phir(k,l,i,id) + crintfac(k,l,id)*phir_z(k,l,i,id) + &
     czintfac_z(k,l,id)*phiz(k,l,i,id) + czintfac(k,l,id)*phiz_z(k,l,i,id) ) &
     *rintfac(k,l,id)*abs(Jp(k,l,id)) + &

     ( crintfac(k,l,id)*phir(k,l,i,id) + czintfac(k,l,id)*phiz(k,l,i,id) ) *rintfac(k,l,id)*Jp_z(k,l,id)*Jpsign(k,l,id)

intRc_c(k,l) = ( phir(k,l,i,id)*phir(k,l,j,id) + phiz(k,l,i,id)*phiz(k,l,j,id) )*rintfac(k,l,id)*abs(Jp(k,l,id))


        else    !VE = 5

!Rt_r_V
if(TtimeS.eq.1) &
intRt_r_V(k,l) = intRt_r_V(k,l) + phi(k,l,i)* ( -CTJ*phi(k,l,j)/dt*Trintfac(k,l,id) - &
     rdotintfac(k,l,id)*Trintfac_r(k,l,id) - zdotintfac(k,l,id)*Tzintfac_r(k,l,id) ) &
     *rintfac(k,l,id)*abs(Jp(k,l,id)) + &
     
     phi(k,l,i)*( Tdotintfac(k,l,id) - rdotintfac(k,l,id) *Trintfac(k,l,id) - &
     zdotintfac(k,l,id) *Tzintfac(k,l,id) ) *rJp_r(k,l,id) 

if(TdiffS.eq.1) &
intRt_r_V(k,l) = intRt_r_V(k,l) + F0*( &
     Trintfac_r(k,l,id)*phir(k,l,i,id) + Trintfac(k,l,id)*phir_r(k,l,i,id) + &
     Tzintfac_r(k,l,id)*phiz(k,l,i,id) + Tzintfac(k,l,id)*phiz_r(k,l,i,id) ) &
     *rintfac(k,l,id)*abs(Jp(k,l,id)) + &
     
     F0*( phir(k,l,i,id)*Trintfac(k,l,id) + phiz(k,l,i,id)*Tzintfac(k,l,id) ) *rJp_r(k,l,id) 

!intRt_z_V
if(TtimeS.eq.1) &
intRt_z_V(k,l) = intRt_z_V(k,l) + phi(k,l,i)* ( -rdotintfac(k,l,id)*Trintfac_z(k,l,id) - &
     CTJ*phi(k,l,j)/dt*Tzintfac(k,l,id) - zdotintfac(k,l,id)*Tzintfac_z(k,l,id) ) &
     *rintfac(k,l,id)*abs(Jp(k,l,id)) + &
     
     phi(k,l,i)*( Tdotintfac(k,l,id) - &
     rdotintfac(k,l,id)*Trintfac(k,l,id) - zdotintfac(k,l,id)*Tzintfac(k,l,id) )  &
     *rintfac(k,l,id)* Jp_z(k,l,id)*Jpsign(k,l,id)!???

if(TdiffS.eq.1) &
intRt_z_V(k,l) = intRt_z_V(k,l) + F0*( &
     Trintfac_z(k,l,id)*phir(k,l,i,id) + Trintfac(k,l,id)*phir_z(k,l,i,id) + &
     Tzintfac_z(k,l,id)*phiz(k,l,i,id) + Tzintfac(k,l,id)*phiz_z(k,l,i,id) ) &
     *rintfac(k,l,id)*abs(Jp(k,l,id)) + &
     
     F0*( phir(k,l,i,id)*Trintfac(k,l,id) + phiz(k,l,i,id)*Tzintfac(k,l,id) )  &
     *rintfac(k,l,id)* Jp_z(k,l,id)*Jpsign(k,l,id)


!Rt_T
if(TtimeS.eq.1) &
intRt_T(k,l) = intRt_T(k,l) + phi(k,l,i)* ( CTJ*phi(k,l,j)/dt - &
     rdotintfac(k,l,id)*phir(k,l,j,id) - zdotintfac(k,l,id)*phiz(k,l,j,id) ) &
     *rintfac(k,l,id)*abs(Jp(k,l,id)) 
if(TdiffS.eq.1) &
intRt_T(k,l) = intRt_T(k,l) + F0*( phir(k,l,i,id)*phir(k,l,j,id) + phiz(k,l,i,id)*phiz(k,l,j,id) ) &
     *rintfac(k,l,id)*abs(Jp(k,l,id)) 

        end if   !for VE

        end if   !for s_mode=0


     end do
  end do   !end loop for k&l, preparation for gaussian_quadrature


  !KBC on free surface has no volume integral
  if( BCflagN( globalNM(m,i), 3 ).ne.1 .and. BCflagN( globalNM(m,i), 3 ).ne.3 ) then
     sj(LNOPP(i)+Nr,LNOPP(j)+Nr) = gaussian_quadrature(intRsi_r_V)            ! sjRsi_r(i,j)   
     sj(LNOPP(i)+Nr,LNOPP(j)+Nz) = gaussian_quadrature(intRsi_z_V)            ! sjRsi_z(i,j)  
     !can be skipped, cause locJac already have initial 0 in assemble
     ! if(s_mode.eq.0) then
     !    if(VE(m).eq.0) then
     !       sj(LNOPP(i)+Nr,LNOPP(j)+Nu) = 0.0_rk                                     ! sjRsi_u(i,j)   
     !       sj(LNOPP(i)+Nr,LNOPP(j)+Nv) = 0.0_rk                                     ! sjRsi_v(i,j)  
     !       if( PN( globalNM(m,j) ).eq.1 ) then
     !          sj(LNOPP(i)+Nr,LNOPP(j)+Np) = 0.0_rk                                  ! sjRsi_p(i,j)  
     !       end if
     !    else    !VE=1
     !       sj(LNOPP(i)+Nr,LNOPP(j)+MDF(globalNM(m,j))-1) = 0.0_rk                   ! sjRsi_c(i,j)   
     !    end if     !for VE
     ! end if   !for s_mode=0
  end if


  sj(LNOPP(i)+Nz,LNOPP(j)+Nr) = gaussian_quadrature(intReta_r_V)           ! sjReta_r(i,j)   
  sj(LNOPP(i)+Nz,LNOPP(j)+Nz) = gaussian_quadrature(intReta_z_V)           ! sjReta_z(i,j)  
  ! if(s_mode.eq.0) then
  !    if(VE(m).eq.0) then
  !       sj(LNOPP(i)+Nz,LNOPP(j)+Nu) = 0.0_rk                                     ! sjReta_u(i,j)   
  !       sj(LNOPP(i)+Nz,LNOPP(j)+Nv) = 0.0_rk                                     ! sjReta_v(i,j)  
  !       if( PN( globalNM(m,j) ).eq.1 ) then
  !          sj(LNOPP(i)+Nz,LNOPP(j)+Np) = 0.0_rk                                  ! sjReta_p(i,j)  
  !       end if
  !    else    !VE=1
  !       sj(LNOPP(i)+Nz,LNOPP(j)+MDF(globalNM(m,j))-1) = 0.0_rk                   ! sjReta_c(i,j)   
  !    end if     !for VE
  ! end if   !for s_mode=0


  if(s_mode.eq.0) then
     if(VE(m).eq.0) then
        
        sj(LNOPP(i)+Nu,LNOPP(j)+Nr) = gaussian_quadrature(intRu_r_V)           ! sjRur(i,j)   
        sj(LNOPP(i)+Nu,LNOPP(j)+Nz) = gaussian_quadrature(intRu_z_V)           ! sjRuz(i,j) 
        sj(LNOPP(i)+Nu,LNOPP(j)+Nu) = gaussian_quadrature(intRu_u)           ! sjRuu(i,j)   
        sj(LNOPP(i)+Nu,LNOPP(j)+Nv) = gaussian_quadrature(intRuv)           ! sjRuv(i,j)  
        if( PN( globalNM(m,j) ).eq.1 ) then
           sj(LNOPP(i)+Nu,LNOPP(j)+Np) = gaussian_quadrature(intRu_p)           ! sjRup(i,j) 
        end if

        sj(LNOPP(i)+Nv,LNOPP(j)+Nr) = gaussian_quadrature(intRv_r_V)           ! sjRvr(i,j)   
        sj(LNOPP(i)+Nv,LNOPP(j)+Nz) = gaussian_quadrature(intRv_z_V)           ! sjRvz(i,j) 
        sj(LNOPP(i)+Nv,LNOPP(j)+Nu) = gaussian_quadrature(intRv_u)           ! sjRvu(i,j)   
        sj(LNOPP(i)+Nv,LNOPP(j)+Nv) = gaussian_quadrature(intRvv)           ! sjRvv(i,j)  
        if( PN( globalNM(m,j) ).eq.1 ) then
           sj(LNOPP(i)+Nv,LNOPP(j)+Np) = gaussian_quadrature(intRv_p)           ! sjRvp(i,j) 
        end if

        if( BCflagN( globalNM(m,i), 3 ).eq.0 ) then  !evaporation cooling on free surface, no volume integral
           if( BCflagN( globalNM(m,i),2 ) .eq. 0 ) then
        sj(LNOPP(i)+NT,LNOPP(j)+Nr) = gaussian_quadrature(intRt_r_V)           ! sjRur(i,j)   
        sj(LNOPP(i)+NT,LNOPP(j)+Nz) = gaussian_quadrature(intRt_z_V)           ! sjRuz(i,j) 
        sj(LNOPP(i)+NT,LNOPP(j)+Nu) = gaussian_quadrature(intRt_u)           ! sjRuu(i,j)   
        sj(LNOPP(i)+NT,LNOPP(j)+Nv) = gaussian_quadrature(intRtv)           ! sjRuv(i,j)  
        sj(LNOPP(i)+NT,LNOPP(j)+NT) = gaussian_quadrature(intRt_T)           ! sjRup(i,j) 
           else  !base nodes 
        sj(LNOPP(i)+NT,LNOPP(j)+Nr) = gaussian_quadrature(intRt_r_V)/kR           ! sjRur(i,j)   
        sj(LNOPP(i)+NT,LNOPP(j)+Nz) = gaussian_quadrature(intRt_z_V)/kR           ! sjRuz(i,j) 
        sj(LNOPP(i)+NT,LNOPP(j)+Nu) = gaussian_quadrature(intRt_u)/kR           ! sjRuu(i,j)   
        sj(LNOPP(i)+NT,LNOPP(j)+Nv) = gaussian_quadrature(intRtv)/kR           ! sjRuv(i,j)  
        sj(LNOPP(i)+NT,LNOPP(j)+NT) = gaussian_quadrature(intRt_T)/kR           ! sjRup(i,j) 
           end if
        end if

        if( PN( globalNM(m,i) ).eq.1 ) then
           sj(LNOPP(i)+Np,LNOPP(j)+Nr) = gaussian_quadrature(intRp_r)           ! sjRpr(i,j)   
           sj(LNOPP(i)+Np,LNOPP(j)+Nz) = gaussian_quadrature(intRp_z)           ! sjRpz(i,j) 
           sj(LNOPP(i)+Np,LNOPP(j)+Nu) = gaussian_quadrature(intRp_u)           ! sjRpu(i,j)   
           sj(LNOPP(i)+Np,LNOPP(j)+Nv) = gaussian_quadrature(intRp_v)           ! sjRpv(i,j)  
           ! if( PN( globalNM(m,j) ).eq.1 ) then
           !    sj(LNOPP(i)+Np,LNOPP(j)+Np) = 0.0_rk                                    ! sjRpp(i,j) 
           ! end if
        end if

     else if( VE(m).eq.1 ) then
        
        sj(LNOPP(i)+MDF(globalNM(m,i))-1,LNOPP(j)+Nr) = gaussian_quadrature(intRc_r)!sjRcr(i,j)
        sj(LNOPP(i)+MDF(globalNM(m,i))-1,LNOPP(j)+Nz) = gaussian_quadrature(intRc_z)!sjRcr(i,j)
        sj(LNOPP(i)+MDF(globalNM(m,i))-1,LNOPP(j)+MDF(globalNM(m,j))-1) = gaussian_quadrature(intRc_c)!sjRcr(i,j)

     else !VE = 5
        if( VN(globalNM(m,i)) .eq. 5 ) then
           sj(LNOPP(i)+NTs,LNOPP(j)+Nr) = gaussian_quadrature(intRt_r_V)           ! sjRur(i,j)   
           sj(LNOPP(i)+NTs,LNOPP(j)+Nz) = gaussian_quadrature(intRt_z_V)           ! sjRuz(i,j) 
           if( VN(globalNM(m,j)) .eq. 5 ) then
              sj(LNOPP(i)+NTs,LNOPP(j)+NTs) = gaussian_quadrature(intRt_T)           ! sjRup(i,j) 
           else   !j is base node
              sj(LNOPP(i)+NTs,LNOPP(j)+NT) = gaussian_quadrature(intRt_T)           ! sjRup(i,j) 
           end if
        else  !i is base node
           sj(LNOPP(i)+NT,LNOPP(j)+Nr) = gaussian_quadrature(intRt_r_V)/F0           ! sjRur(i,j)   
           sj(LNOPP(i)+NT,LNOPP(j)+Nz) = gaussian_quadrature(intRt_z_V)/F0           ! sjRuz(i,j) 
           if( VN(globalNM(m,j)) .eq. 5 ) then
              sj(LNOPP(i)+NT,LNOPP(j)+NTs) = gaussian_quadrature(intRt_T)/F0           ! sjRup(i,j) 
           else   !j is base node
              sj(LNOPP(i)+NT,LNOPP(j)+NT) = gaussian_quadrature(intRt_T)/F0           ! sjRup(i,j) 
           end if
        end if
     end if   !for  VE

  end if   !for s_mode=0

! if(initial_vapor_solved.eq.1 .and. m.eq.28 .and. i.eq.4 .and. j.eq.2) then
!    write(*,*) sj(LNOPP(i)+NT,LNOPP(j)+Nr)
!    pause
! end if

  return
end subroutine VI_in_sj

