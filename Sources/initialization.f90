subroutine initialization
  use kind
  use data
  use Ldata
  use front_mod, only: init_front, piv


  implicit none


  
  allocate( sol(NVar), dsol(NVar) )
  allocate( solp(NVar), soldot(NVar), soldotp(NVar), soldotpp(NVar), solpred(NVar) )
  allocate( rcoordinate(NTN), zcoordinate(NTN), usol(NTN), vsol(NTN), Tsol(NTN), psol(NTN), csol(NTN) )
     allocate( fsi_size(NTE), geta_size(NTE) )
  sol = 0.0_rk
  dsol = 0.0_rk
  solp = 0.0_rk
  soldot = 0.0_rk
  soldotp = 0.0_rk
  soldotpp = 0.0_rk
  solpred = 0.0_rk
  rcoordinate = 0.0_rk
  zcoordinate = 0.0_rk
  usol = 0.0_rk
  vsol = 0.0_rk
  Tsol = 0.0_rk
  psol = 0.0_rk
  csol = 0.0_rk

  if(check_0_in_Jac.eq.1)  then
     allocate( Jac(NVar,NVar), Res(NVar) )
     Jac = 0.0_rk
     Res = 0.0_rk
  end if


  if(ths.eq.1) then
     call init_front(1, -1, 1)
  else if(ths.eq.2) then
     call init_front(1, -2, 1)
  else if(ths.ge.4) then
     ! call init_front(1, 9, 16)    !9 is fastest for ptn, always use 16 for ptn2
     call init_front(1, -3, 1)     !then use ptn = 8, ptn2 = 16, slower
  else
     write(*,*) 'error in ths'
     stop
  end if
  !NOP, NOPP, MDF, rNOP, NE etc must be defined before this step
  !subroutine init_front(init,ptn,ptn2)
  !init = 1 for first entry, 2 for changing setup after first, 0 for exit
  !ptn = size of sub rows in middle (-1 for overide single, -2 for overide double, -3 optimal nat/nest), -4 drops
  !ptn2 = size of sub columns in middle (if ptn<0 make ptn2 1)

  step = 0
  time = 0.0_rk
  timestep = 0
  size_function_change = 0
  final_size = 0
  initial_vapor_solved = 0   !start with 0, adjust the vapor concentration before solving for dynamics
  initial_vapor_solving = 0
  diverge = 0
  zone = 0

  angle_cd = 0.0_rk

  ! fluxuni = 0.0_rk
  ! if(uniflux) fluxuni = 1.0_rk


  !-------------------------------------allocate local data-------------------------------------
  allocate( rlocal(9,ths), zlocal(9,ths), ulocal(9,ths), vlocal(9,ths), &
       Tlocal(9,ths), plocal(9,ths), clocal(9,ths) )
  allocate( rintfac(3,3,ths), rsi(3,3,ths), reta(3,3,ths), zsi(3,3,ths), zeta(3,3,ths) )
  allocate( Jp(3,3,ths), Jpsign(3,3,ths), s_orth(3,3,ths) )
  allocate( phir(3,3,9,ths), phiz(3,3,9,ths) )
  allocate( uintfac(3,3,ths), urintfac(3,3,ths), uzintfac(3,3,ths), &
       vintfac(3,3,ths), vrintfac(3,3,ths), vzintfac(3,3,ths), pintfac(3,3,ths), &
       crintfac(3,3,ths), czintfac(3,3,ths), Trintfac(3,3,ths), Tzintfac(3,3,ths) )
  allocate( rdotintfac(3,3,ths), zdotintfac(3,3,ths), &
       udotintfac(3,3,ths), vdotintfac(3,3,ths), Tdotintfac(3,3,ths) )
  allocate( udotlocal(9,ths), vdotlocal(9,ths), rdotlocal(9,ths), zdotlocal(9,ths), Tdotlocal(9,ths) )
  allocate( rsi_down(3,ths),  zsi_down(3,ths), rsi_left(3,ths), zsi_left(3,ths), &
       reta_left(3,ths), zeta_left(3,ths), &
       reta_right(3,ths), zeta_right(3,ths), rsi_right(3,ths), zsi_right(3,ths), Teta_right(3,ths) )
  allocate( rintfac_right(3,ths), uintfac_right(3,ths), vintfac_right(3,ths), &
       rdotintfac_right(3,ths), zdotintfac_right(3,ths), &
       Jp_r_right(3,ths), Jp_z_right(3,ths), Jp_right(3,ths) )
  allocate( dcdsi(3,ths), dcdeta(3,ths) ,dTdsi(3,ths) ,dTdeta(3,ths) )
  allocate( Aterm(3,3,ths), Bterm(3,3,ths) )
  allocate( s_orth_r(3,3,ths), s_orth_z(3,3,ths) )
  allocate( Aterm_r(3,3,ths),Aterm_z(3,3,ths), Bterm_r(3,3,ths), Bterm_z(3,3,ths) )
  allocate( Jp_r(3,3,ths), Jp_z(3,3,ths), rJp_r(3,3,ths) )
  allocate( phir_r(3,3,9,ths), phir_z(3,3,9,ths), phiz_r(3,3,9,ths), phiz_z(3,3,9,ths) )
  allocate( urintfac_r(3,3,ths), urintfac_z(3,3,ths), uzintfac_r(3,3,ths), uzintfac_z(3,3,ths), &
       vrintfac_r(3,3,ths), vrintfac_z(3,3,ths), vzintfac_r(3,3,ths), vzintfac_z(3,3,ths), &
       crintfac_r(3,3,ths), crintfac_z(3,3,ths), czintfac_r(3,3,ths), czintfac_z(3,3,ths), &
       Trintfac_r(3,3,ths), Trintfac_z(3,3,ths), Tzintfac_r(3,3,ths), Tzintfac_z(3,3,ths) )
  !multifront module
  allocate( piv(ths) )


  return
end subroutine initialization
