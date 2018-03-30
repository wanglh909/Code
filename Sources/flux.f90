real(kind=rk) function flux(theta, r)
use kind
use data, only: Hum, pi, true_uniflux


implicit none

real(kind=rk), intent(in):: theta, r
real(kind=rk):: J0, lambda

if(r.eq.1.0_rk) then
   flux = 1.0e5_rk
else
   J0 = (1.0_rk-Hum) *( 0.27_rk*theta**2 + 1.3_rk ) *( 0.6381 - 0.2239 *( theta - pi/4.0_rk )**2 )
   lambda = 0.5_rk - theta/pi
   flux = J0 *( 1.0_rk - r**2 ) **(-lambda)
end if

! flux = (1.0_rk-Hum) *( 0.27_rk*theta**2 + 1.3_rk ) *( 0.6381 - 0.2239 *( theta - pi/4.0_rk )**2 )

! flux = 1.0_rk

return

end function flux

