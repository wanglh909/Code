module basis_f
use kind
contains

subroutine basis_function
  use data, only: phi, phisi, phieta, psi, psisi, psieta, phi_1d, phix_1d, phixx_1d, &
       phisi0_1d, phisi1_1d, phieta0_1d, phieta1_1d, phi0_1d, phi1_1d, phix_1d_flux, convert49, Ng


  implicit none

  integer(kind=ik):: j,k,l

  real(kind=rk), parameter, dimension(Ng):: gausspoint = &
       (/0.5_rk*(-0.7745966692414833770358531_rk + 1.0_rk), &
       0.5_rk*(0.000000000000000000000000000000_rk + 1.0_rk), &
       0.5_rk*(0.7745966692414833770358531_rk + 1.0_rk)/)


!*******************save phii into phi, phiisi into phisi, phiieta into phieta************************
!*******************save psii into psi, psiisi into psisi, psiieta into psieta************************
  do j = 1, 9, 1
     do k = 1, 3, 1
        do l = 1, 3, 1
           phi(k,l,j) = phii( gausspoint(k), gausspoint(l), j )
           phisi(k,l,j) = phiisi( gausspoint(k), gausspoint(l), j )
           phieta(k,l,j) = phiieta( gausspoint(k), gausspoint(l), j )
        end do
     end do
  end do
  ! write(*,*) 'phi:', phi

  psi(:,:,:) = 0.0_rk     !except 1, 3, 7, 9, other elements are 0.0_rk
  psisi(:,:,:) = 0.0_rk
  psieta(:,:,:) = 0.0_rk  
  do j = 1, 4, 1            !except psi(:,:,1) psi(:,:,3) psi(:,:,7) psi(:,:,9), psi(:,:,:) = 0.0_rk
     do k = 1, 3, 1
        do l = 1, 3, 1
           psi( k,l, convert49(j) ) = psii( gausspoint(k), gausspoint(l), j )
           psisi( k,l,convert49(j) ) = psiisi( gausspoint(l), j )
           psieta( k,l,convert49(j) ) = psiieta( gausspoint(k), j )
        end do
     end do
  end do
  !write(*,*) 'psi:', psi

  do j = 1, 3, 1
     do k = 1, 3, 1
        phi_1d(k,j) = phii_1d( gausspoint(k), j )
        phix_1d(k,j) = phiix_1d( gausspoint(k), j )
     end do
  end do
  phixx_1d(1) = 4.0_rk
  phixx_1d(2) = -8.0_rk
  phixx_1d(3) = 4.0_rk

  do j = 1, 9
     do k = 1, 3
        phisi0_1d(k,j) =  phiisi( 0.0_rk, gausspoint(k), j )
        phisi1_1d(k,j) =  phiisi( 1.0_rk, gausspoint(k), j )
        phieta0_1d(k,j) =  phiieta( 0.0_rk, gausspoint(k), j )
        phieta1_1d(k,j) =  phiieta( 1.0_rk, gausspoint(k), j )
        phi0_1d(k,j) = phii( 0.0_rk, gausspoint(k), j )
        phi1_1d(k,j) = phii( 1.0_rk, gausspoint(k), j )
     end do
  end do

  do j = 1, 3, 1
     phix_1d_flux(j) = phiix_1d( 0.0_rk, j )   !for calculating flux
  end do


!***************************************************************************************************

  return
end subroutine basis_function

! contains





real(kind=rk) function phii(si, eta, i)
  implicit none
  real(kind=rk):: si, eta
  integer(kind=ik):: i

  select case(i)
  case(1)
     phii = ( 1.0_rk - 3.0_rk*si + 2.0_rk*si**2 )*( 1.0_rk - 3.0_rk*eta + 2.0_rk*eta**2 )
  case(2)
     phii = ( 4.0_rk*si - 4.0_rk*si**2 )*( 1.0_rk - 3.0_rk*eta + 2.0_rk*eta**2 )
  case(3)
     phii = ( -si + 2.0_rk*si**2 )*( 1.0_rk - 3.0_rk*eta + 2.0_rk*eta**2 )
  case(4)
     phii = ( 1.0_rk - 3.0_rk*si + 2.0_rk*si**2 )*( 4.0_rk*eta - 4.0_rk*eta**2 )
  case(5)
     phii = ( 4.0_rk*si - 4.0_rk*si**2 )*( 4.0_rk*eta - 4.0_rk*eta**2 )
  case(6)
     phii = ( -si + 2.0_rk*si**2 )*( 4.0_rk*eta - 4.0_rk*eta**2 )
  case(7)
     phii = ( 1.0_rk - 3.0_rk*si + 2.0_rk*si**2 )*( -eta + 2.0_rk*eta**2 )
  case(8)
     phii = ( 4.0_rk*si - 4.0_rk*si**2 )*( -eta + 2.0_rk*eta**2 )
  case(9)
     phii = ( -si + 2.0_rk*si**2 )*( -eta + 2.0_rk*eta**2 )
  case default
     write(*,*) 'error in i, in phii'
  end select
  return
end function phii


real(kind=rk) function phiisi(si, eta, i)
  implicit none
  real(kind=rk):: si, eta
  integer(kind=ik):: i

  select case(i)
  case(1)
     phiisi = ( -3.0_rk + 4.0_rk*si )*( 1.0_rk - 3.0_rk*eta + 2.0_rk*eta**2 )
  case(2)
     phiisi = ( 4.0_rk - 8.0_rk*si )*( 1.0_rk - 3.0_rk*eta + 2.0_rk*eta**2 )
  case(3)
     phiisi = ( -1.0_rk + 4.0_rk*si )*( 1.0_rk - 3.0_rk*eta + 2.0_rk*eta**2 )
  case(4)
     phiisi = ( -3.0_rk + 4.0_rk*si )*( 4.0_rk*eta - 4.0_rk*eta**2 )
  case(5)
     phiisi = ( 4.0_rk - 8.0_rk*si )*( 4.0_rk*eta - 4.0_rk*eta**2 )
  case(6)
     phiisi = ( -1.0_rk + 4.0_rk*si )*( 4.0_rk*eta - 4.0_rk*eta**2 )
  case(7)
     phiisi = ( -3.0_rk + 4.0_rk*si )*( -eta + 2.0_rk*eta**2 )
  case(8)
     phiisi = ( 4.0_rk - 8.0_rk*si )*( -eta + 2.0_rk*eta**2 )
  case(9)
     phiisi = ( -1.0_rk + 4.0_rk*si )*( -eta + 2.0_rk*eta**2 )
  case default
     write(*,*) 'error in i, in phiisi'
  end select
  return
end function phiisi


real(kind=rk) function phiieta(si, eta, i) 
  implicit none
  real(kind=rk):: si, eta
  integer(kind=ik):: i

  select case(i)
  case(1)
     phiieta = ( 1.0_rk - 3.0_rk*si + 2.0_rk*si**2 )*( -3.0_rk + 4.0_rk*eta )
  case(2)
     phiieta = ( 4.0_rk*si - 4.0_rk*si**2 )*( -3.0_rk + 4.0_rk*eta )
  case(3)
     phiieta = ( -si + 2.0_rk*si**2 )*( -3.0_rk + 4.0_rk*eta )
  case(4)
     phiieta = ( 1.0_rk - 3.0_rk*si + 2.0_rk*si**2 )*( 4.0_rk - 8.0_rk*eta )
  case(5)
     phiieta = ( 4.0_rk*si - 4.0_rk*si**2 )*( 4.0_rk - 8.0_rk*eta )
  case(6)
     phiieta = ( -si + 2.0_rk*si**2 )*( 4.0_rk - 8.0_rk*eta )
  case(7)
     phiieta = ( 1.0_rk - 3.0_rk*si + 2.0_rk*si**2 )*( -1.0_rk + 4.0_rk*eta )
  case(8)
     phiieta = ( 4.0_rk*si - 4.0_rk*si**2 )*( -1.0_rk + 4.0_rk*eta )
  case(9)
     phiieta = ( -si + 2.0_rk*si**2 )*( -1.0_rk + 4.0_rk*eta )
  case default
     write(*,*) 'error in i, in phiieta'
  end select
  return
end function phiieta





real(kind=rk) function psii(si, eta, i)
  implicit none
  real(kind=rk):: si, eta
  integer(kind=ik):: i

  select case(i)
  case(1)
     psii = (1.0_rk-si)*(1.0_rk-eta)
  case(2)
     psii = si*(1.0_rk-eta)
  case(3)
     psii = (1.0_rk-si)*eta
  case(4)
     psii = si*eta
  case default
     write(*,*) 'error in i, in psii'
  end select
  return
end function psii


real(kind=rk) function psiisi(eta, i)
  implicit none
  real(kind=rk):: eta
  integer(kind=ik):: i

  select case(i)
  case(1)
     psiisi = (-1.0_rk)*(1.0_rk-eta)
  case(2)
     psiisi = 1.0_rk-eta
  case(3)
     psiisi = (-1.0_rk)*eta
  case(4)
     psiisi = eta
  case default
     write(*,*) 'error in i, in psiisi'
  end select
  return
end function psiisi


real(kind=rk) function psiieta(si, i)
  implicit none
  real(kind=rk):: si
  integer(kind=ik):: i

  select case(i)
  case(1)
     psiieta = (1.0_rk-si)*(-1.0_rk)
  case(2)
     psiieta = si*(-1.0_rk)
  case(3)
     psiieta = (1.0_rk-si)
  case(4)
     psiieta = si
  case default
     write(*,*) 'error in i, in psiieta'
  end select
  return
end function psiieta


real(kind=rk) function phii_1d(x,i) 
  implicit none
  real(kind=rk):: x
  integer(kind=ik) :: i

  select case (i)
  case(1)
     phii_1d = 1.0_rk - 3.0_rk*x + 2.0_rk*x**2
  case(2)
     phii_1d = 4.0_rk*( x - x**2 )
  case(3)
     phii_1d = -x + 2.0_rk*x**2
  case default
     write(*,*) 'error in i, in phii_1d'
  end select

  return
end function phii_1d



real(kind=rk) function phiix_1d(x,i)
  implicit none

  real(kind=rk):: x
  integer(kind=ik):: i
  select case(i)
  case(1)
     phiix_1d = -3.0_rk + 4.0_rk*x
  case(2)
     phiix_1d = 4.0_rk - 8.0_rk*x 
  case(3)
     phiix_1d = -1.0_rk + 4.0_rk*x

  case default
     write(*,*) 'error in i, in phiix_1d'
  end select
  return
end function phiix_1d








! end subroutine basis_function

end module basis_f
