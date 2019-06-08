real(kind=rk) function flux_f(theta, r, m)
use kind
use data, only: Hum, pi


implicit none

real(kind=rk), intent(in):: theta, r
integer(kind=ik), intent(in):: m
real(kind=rk):: J0, lambda


! J0 = (1.0_rk-Hum) *( 0.27_rk*theta**2 + 1.3_rk ) *( 0.6381_rk - 0.2239_rk *( theta - pi/4.0_rk )**2 )
! lambda = 0.5_rk - theta/pi
      
! if(m.eq.1) then  !flux
!    if(r.eq.1.0_rk) then
!       flux_f = 1.0e5_rk
!    else if (r.gt.1.0_rk) then
!       write(*,*) '!error of r value for flux, larger than 1.0, r =', r
!       flux_f = 1.0e5_rk
!       ! pause
!    else
!       flux_f = J0 *( 1.0_rk - r**2 ) **(-lambda)
!    end if
   
! else if(m.eq.2) then  !d flux/ dr
!    if(r.eq.1.0_rk) then
!       flux_f = 1.0e5_rk
!    else if (r.gt.1.0_rk) then
!       write(*,*) '!error of r value for flux_r, larger than 1.0, r =', r
!       flux_f = 1.0e5_rk
!       ! pause
!    else
!       flux_f = 2.0_rk *J0 *lambda* r*( 1.0_rk - r**2 ) **(-lambda-1.0_rk)
!    end if
! end if

if(m.eq.1) then  !flux
   flux_f = 1.0_rk
else if(m.eq.2) then  !d flux/ dr
   flux_f = 0.0_rk
end if
   

return

end function flux_f



real(kind=rk) function Resisp_f(cp,m)
use kind
use data, only: cp_pack


implicit none

real(kind=rk), intent(in):: cp
integer(kind=ik), intent(in):: m
real(kind=rk):: steepness, drop_location, exponent_reciprocal1, exponent_reciprocal2

! steepness = 8.0e1_rk
! drop_location = 0.999_rk
! exponent_reciprocal1 = 1.0_rk / ( 1.0_rk + exp( -steepness * (cp/cp_pack - drop_location) ) )
! exponent_reciprocal2 = 1.0_rk / ( 1.0_rk/exp( -steepness * (cp/cp_pack - drop_location) ) + 1.0_rk )
! if(m.eq.1) then   !Resisp
!    Resisp_f = 1.0_rk - exponent_reciprocal1
! else if(m.eq.2) then   !d(Resisp)/d(cp)
!    Resisp_f = (-steepness) / cp_pack * exponent_reciprocal1 * exponent_reciprocal2
!    ! print *, 'Resisp function', Resisp_f
!    ! pause
! else if(m.eq.3) then   !d2(Resisp)/d(cp)2
!    Resisp_f = steepness**2 / (cp_pack**2) * exponent_reciprocal1* exponent_reciprocal2 * &
!         ( -2.0_rk*exponent_reciprocal2 + 1.0_rk )
! else
!    print *, 'Resisp_f(m) wrong number'
! end if


if(m.eq.1) then   !Resisp
   Resisp_f = 1.0_rk
else if(m.eq.2) then   !d(Resisp)/d(cp)
   Resisp_f = 0.0_rk
else if(m.eq.3) then   !d2(Resisp)/d(cp)2
   Resisp_f = 0.0_rk
else
   print *, 'Resisp_f(m) wrong number'
end if



return
end function Resisp_f
