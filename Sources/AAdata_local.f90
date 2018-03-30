module Ldata
  use kind
  use data, only:ths
  implicit none


  !for subroutine values_in_an_element
  real(kind=rk), allocatable:: rlocal(:,:), zlocal(:,:), &
       ulocal(:,:), vlocal(:,:), Tlocal(:,:), plocal(:,:), clocal(:,:)
  integer(kind=ik), parameter:: NNr = 1, NNz = 2, NNu = 3, NNv = 4, NNp = 6
  real(kind=rk), allocatable:: rintfac(:,:,:), rsi(:,:,:), reta(:,:,:), zsi(:,:,:), zeta(:,:,:)
  real(kind=rk), allocatable:: Jp(:,:,:), Jpsign(:,:,:), s_orth(:,:,:)
  real(kind=rk), allocatable:: phir(:,:,:,:), phiz(:,:,:,:)
  real(kind=rk), allocatable:: uintfac(:,:,:), urintfac(:,:,:), uzintfac(:,:,:), &
       vintfac(:,:,:), vrintfac(:,:,:), vzintfac(:,:,:), pintfac(:,:,:), &
       crintfac(:,:,:), czintfac(:,:,:), Trintfac(:,:,:), Tzintfac(:,:,:)
  real(kind=rk), allocatable:: rdotintfac(:,:,:), zdotintfac(:,:,:), &
       udotintfac(:,:,:), vdotintfac(:,:,:), Tdotintfac(:,:,:)
  real(kind=rk), allocatable:: udotlocal(:,:), vdotlocal(:,:), &
       rdotlocal(:,:), zdotlocal(:,:), Tdotlocal(:,:)
  real(kind=rk), allocatable:: rsi_down(:,:),  zsi_down(:,:), &
       rsi_left(:,:), zsi_left(:,:), reta_left(:,:), zeta_left(:,:), &
       reta_right(:,:), zeta_right(:,:), rsi_right(:,:), zsi_right(:,:), Teta_right(:,:)
  real(kind=rk), allocatable:: rintfac_right(:,:), uintfac_right(:,:), vintfac_right(:,:), &
       rdotintfac_right(:,:), zdotintfac_right(:,:), Jp_r_right(:,:), Jp_z_right(:,:), Jp_right(:,:)

  real(kind=rk), allocatable:: dcdsi(:,:), dcdeta(:,:) ,dTdsi(:,:) ,dTdeta(:,:)

  !for subroutine define_sf
  real(kind=rk), allocatable:: Aterm(:,:,:), Bterm(:,:,:)
  ! real(kind=rk), allocatable:: sf(:,:)
  real(kind=rk):: M1 = 0.0_rk, M2 = 0.0_rk
  real(kind=rk), parameter:: eps1 = 1.0_rk, eps2 = 1.0_rk, epss = 0.1_rk
  !real(kind=rk), parameter:: eps1 = 0.1_rk, eps2 = 0.1_rk, epss = 1.0_rk
  


  !for subroutine values_in_sj
  real(kind=rk), allocatable:: s_orth_r(:,:,:), s_orth_z(:,:,:)
  real(kind=rk), allocatable:: Aterm_r(:,:,:),Aterm_z(:,:,:), Bterm_r(:,:,:), Bterm_z(:,:,:)
  real(kind=rk), allocatable:: Jp_r(:,:,:), Jp_z(:,:,:), rJp_r(:,:,:)
  real(kind=rk), allocatable::  phir_r(:,:,:,:), phir_z(:,:,:,:), phiz_r(:,:,:,:), phiz_z(:,:,:,:)
  real(kind=rk), allocatable:: urintfac_r(:,:,:), urintfac_z(:,:,:), uzintfac_r(:,:,:), uzintfac_z(:,:,:), &
       vrintfac_r(:,:,:), vrintfac_z(:,:,:), vzintfac_r(:,:,:), vzintfac_z(:,:,:), &
       crintfac_r(:,:,:), crintfac_z(:,:,:), czintfac_r(:,:,:), czintfac_z(:,:,:), &
       Trintfac_r(:,:,:), Trintfac_z(:,:,:), Tzintfac_r(:,:,:), Tzintfac_z(:,:,:)


  ! !for subroutine values_in_an_element
  ! real(kind=rk):: rlocal(9,ths), zlocal(9,ths), ulocal(9,ths), vlocal(9,ths), Tlocal(9,ths), plocal(9,ths), clocal(9,ths)
  ! integer(kind=ik), parameter:: NNr = 1, NNz = 2, NNu = 3, NNv = 4, NNp = 6
  ! real(kind=rk):: rintfac(3,3,ths), rsi(3,3,ths), reta(3,3,ths), zsi(3,3,ths), zeta(3,3,ths)
  ! real(kind=rk):: Jp(3,3,ths), Jpsign(3,3,ths), s_orth(3,3,ths)
  ! real(kind=rk):: phir(3,3,9,ths), phiz(3,3,9,ths)
  ! real(kind=rk):: uintfac(3,3,ths), urintfac(3,3,ths), uzintfac(3,3,ths), &
  !      vintfac(3,3,ths), vrintfac(3,3,ths), vzintfac(3,3,ths), pintfac(3,3,ths), &
  !      crintfac(3,3,ths), czintfac(3,3,ths), Trintfac(3,3,ths), Tzintfac(3,3,ths)
  ! real(kind=rk):: rdotintfac(3,3,ths), zdotintfac(3,3,ths), udotintfac(3,3,ths), vdotintfac(3,3,ths), Tdotintfac(3,3,ths)
  ! real(kind=rk):: udotlocal(9,ths), vdotlocal(9,ths), rdotlocal(9,ths), zdotlocal(9,ths), Tdotlocal(9,ths)
  ! real(kind=rk):: rsi_down(3,ths),  zsi_down(3,ths), rsi_left(3,ths), zsi_left(3,ths), reta_left(3,ths), zeta_left(3,ths), &
  !      reta_right(3,ths), zeta_right(3,ths), rsi_right(3,ths), zsi_right(3,ths), Teta_right(3,ths)
  ! real(kind=rk):: rintfac_right(3,ths), uintfac_right(3,ths), vintfac_right(3,ths), &
  !      rdotintfac_right(3,ths), zdotintfac_right(3,ths), Jp_r_right(3,ths), Jp_z_right(3,ths), Jp_right(3,ths)

  ! real(kind=rk):: dcdsi(3,ths), dcdeta(3,ths) ,dTdsi(3,ths) ,dTdeta(3,ths)

  ! !for subroutine define_sf
  ! real(kind=rk):: Aterm(3,3,ths), Bterm(3,3,ths)
  ! ! real(kind=rk), allocatable:: sf(:,:)
  ! real(kind=rk):: M1 = 0.0_rk, M2 = 0.0_rk
  ! real(kind=rk), parameter:: eps1 = 1.0_rk, eps2 = 1.0_rk, epss = 0.1_rk
  ! !real(kind=rk), parameter:: eps1 = 0.1_rk, eps2 = 0.1_rk, epss = 1.0_rk
  


  ! !for subroutine values_in_sj
  ! real(kind=rk):: s_orth_r(3,3,ths), s_orth_z(3,3,ths)
  ! real(kind=rk):: Aterm_r(3,3,ths),Aterm_z(3,3,ths), Bterm_r(3,3,ths), Bterm_z(3,3,ths)
  ! real(kind=rk):: Jp_r(3,3,ths), Jp_z(3,3,ths), rJp_r(3,3,ths)
  ! real(kind=rk)::  phir_r(3,3,9,ths), phir_z(3,3,9,ths), phiz_r(3,3,9,ths), phiz_z(3,3,9,ths)
  ! real(kind=rk):: urintfac_r(3,3,ths), urintfac_z(3,3,ths), uzintfac_r(3,3,ths), uzintfac_z(3,3,ths), &
  !      vrintfac_r(3,3,ths), vrintfac_z(3,3,ths), vzintfac_r(3,3,ths), vzintfac_z(3,3,ths), &
  !      crintfac_r(3,3,ths), crintfac_z(3,3,ths), czintfac_r(3,3,ths), czintfac_z(3,3,ths), &
  !      Trintfac_r(3,3,ths), Trintfac_z(3,3,ths), Tzintfac_r(3,3,ths), Tzintfac_z(3,3,ths)

       


end module Ldata
