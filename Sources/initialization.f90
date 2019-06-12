subroutine initialization
  use kind
  use data
  use Ldata
  use front_mod, only: init_front, &!piv, &
       assembler, associater, excluder, custom_order, var_finder, seed


  implicit none

!---------------------------------------FRONT SETUP------------------------------------
  !This setup must be done once
  
  !These must be declared external so front can point to them
  external assemble_local
  external associate_arrays
  external find_var_info
  !external get_custom_numbering
  !external exclude_check

  !Point internal front calls to external subroutine
  assembler => assemble_local               !REQUIRED
  associater => associate_arrays            !REQUIRED
  var_finder => find_var_info               !NOT REQUIRED, gives more info for NaN crash
  !excluder => exclude_check                !NOT REQUIRED, default include all
  !custom_order => get_custom_numbering     !NOT REQUIRED, default order 1 -> NE
!------------------------------------------------------------------------------------

  
  allocate( sol(NVar), dsol(NVar) )
  allocate( solp(NVar), soldot(NVar), soldotp(NVar), soldotpp(NVar), solpred(NVar) )
  allocate( rcoordinate(NTN), zcoordinate(NTN), usol(NTN), vsol(NTN), Tsol(NTN), psol(NTN), csol(NTN), cpsol(NTN), gammasol(NTN) )
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
  cpsol = 1.0_rk
  gammasol = 0.0_rk
  cp_average = 1.0_rk

  if(check_0_in_Jac.eq.1)  then
     allocate( Jac(NVar,NVar), Res(NVar) )
     Jac = 0.0_rk
     Res = 0.0_rk
  end if


  step = 0
  time = 0.0_rk
  timestep = 0
  read_coordinate_value = 0
  size_function_change = 0
  final_size = 0
  initial_vapor_solved = 0   !start with 0, adjust the vapor concentration before solving for dynamics
  initial_vapor_solving = 0
  diverge = 0
  pack_condition = 0.0_rk
  if(Maran_flow.eq.1) then  !change from 0 to 1 if solving for Marangoni cases
     init_stability = 1
  else
     init_stability = 1
  end if
  timestep_stable = 0
  radial_cal_time = 0
  
  angle_cd = 0.0_rk


  rmax = 1.0_rk
  zmax = 1.0_rk
  umax = 1.0_rk
  vmax = 1.0_rk
  Tmax = 1.0_rk
  cpmax = 1.0_rk
  pmax = 1.0_rk
  cmax = 1.0_rk

  !-------------------------------------allocate local data-------------------------------------
  allocate( rlocal(9,ths), zlocal(9,ths), ulocal(9,ths), vlocal(9,ths), &
       Tlocal(9,ths), plocal(9,ths), clocal(9,ths), cplocal(9,ths), gammalocal(9,ths) )
  allocate( rintfac(3,3,ths), rsi(3,3,ths), reta(3,3,ths), zsi(3,3,ths), zeta(3,3,ths) )
  allocate( Jp(3,3,ths), Jpsign(3,3,ths), s_orth(3,3,ths) )
  allocate( phir(3,3,9,ths), phiz(3,3,9,ths) )
  allocate( uintfac(3,3,ths), urintfac(3,3,ths), uzintfac(3,3,ths), &
       vintfac(3,3,ths), vrintfac(3,3,ths), vzintfac(3,3,ths), pintfac(3,3,ths), &
       crintfac(3,3,ths), czintfac(3,3,ths), Trintfac(3,3,ths), Tzintfac(3,3,ths), &
       cpintfac(3,3,ths), cprintfac(3,3,ths), cpzintfac(3,3,ths), mu(3,3,ths) )
  allocate( rdotintfac(3,3,ths), zdotintfac(3,3,ths), &
       udotintfac(3,3,ths), vdotintfac(3,3,ths), Tdotintfac(3,3,ths), cpdotintfac(3,3,ths) )
  allocate( udotlocal(9,ths), vdotlocal(9,ths), rdotlocal(9,ths), zdotlocal(9,ths), &
       Tdotlocal(9,ths), cpdotlocal(9,ths), gammadotlocal(9,ths) )
  allocate( rsi_down(3,ths),  zsi_down(3,ths), rsi_left(3,ths), zsi_left(3,ths), &
       reta_left(3,ths), zeta_left(3,ths), &
       reta_right(3,ths), zeta_right(3,ths), rsi_right(3,ths), zsi_right(3,ths), Teta_right(3,ths), cpeta_right(3,ths) )
  allocate( cpintfac_left(3,ths), rintfac_left(3,ths) )
  allocate( rintfac_right(3,ths), uintfac_right(3,ths), vintfac_right(3,ths), cpintfac_right(3,ths), &
       rdotintfac_right(3,ths), zdotintfac_right(3,ths), &
       Jp_r_right(3,ths), Jp_z_right(3,ths), Jp_right(3,ths) )
  allocate( gammaintfac(3,ths), gammadot(3,ths), gammaeta(3,ths), ueta(3,ths), veta(3,ths) )!, &
       !retaeta(3,ths), zetaeta(3,ths), gammaetaeta(3,ths) )
  
  allocate( dcdsi(3,ths), dcdeta(3,ths), dTdsi(3,ths), dcpdsi(3,ths) )!,dTdeta(3,ths) )
  allocate( Aterm(3,3,ths), Bterm(3,3,ths) )
  allocate( s_orth_r(3,3,ths), s_orth_z(3,3,ths) )
  allocate( Aterm_r(3,3,ths),Aterm_z(3,3,ths), Bterm_r(3,3,ths), Bterm_z(3,3,ths) )
  allocate( Jp_r(3,3,ths), Jp_z(3,3,ths), rJp_r(3,3,ths) )
  allocate( phir_r(3,3,9,ths), phir_z(3,3,9,ths), phiz_r(3,3,9,ths), phiz_z(3,3,9,ths) )
  allocate( urintfac_r(3,3,ths), urintfac_z(3,3,ths), uzintfac_r(3,3,ths), uzintfac_z(3,3,ths), &
       vrintfac_r(3,3,ths), vrintfac_z(3,3,ths), vzintfac_r(3,3,ths), vzintfac_z(3,3,ths), &
       crintfac_r(3,3,ths), crintfac_z(3,3,ths), czintfac_r(3,3,ths), czintfac_z(3,3,ths), &
       Trintfac_r(3,3,ths), Trintfac_z(3,3,ths), Tzintfac_r(3,3,ths), Tzintfac_z(3,3,ths), &
       cprintfac_r(3,3,ths), cprintfac_z(3,3,ths), cpzintfac_r(3,3,ths), cpzintfac_z(3,3,ths) )
  allocate( flux(3,ths), flux_r(3,ths) )
  allocate( SQr2z2(3,ths), dS(3,ths), dSQdr(3,ths), dSQdz(3,ths) )
  allocate( adsp_rate(3,ths), Rms1term(3,ths), ureandvze(3,ths), fourterms(3,ths), &
       Rms3_1(3,ths), rereeandzezee(3,ths), Rms3_2(3,ths) )
  allocate( rezeemzeree(3,ths), twoHterm(3,ths), rdzeandzdre(3,ths), &  
       dSdr(3,ths), dSdz(3,ths) )
  allocate( Rms5term(3,ths), reueandzeve(3,ths), Rms6term(3,ths) )  !rdoteta(3,ths), zdoteta(3,ths), 


!--------------------------------------multifront------------------------------------
  call omp_set_nested(.FALSE.) !Multifront is no longer nested

  !multifront module
  call init_front(SOLVER_MODE)
  seed = NTE


  
  
  return
end subroutine initialization
