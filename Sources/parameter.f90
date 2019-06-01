subroutine parameter_values
  use kind
  use data, only: Re, Ca, Kdi, Pe, KBCgroup, REH, beta, Oh, Grav, MaN, R, Hum, F0, kR, folder, substrate, outer, &
       NStrans, Inert, Capil, Viscous, GravI, Ttime, Tconv, Tdiff, TtimeS, TdiffS, NEM, NEL, NES, NEV, NEM_alge, T_sub, uniflux, &
       diameterp, Pep, kboltz, pi, solve_T, Maran_flow, fixed_Ma, cp_pack, Dp, no_vapor, solve_cp
  implicit none

  integer(kind=ik):: Ltype !, water, octane, hexanol
  real(kind=rk):: rho, mu, kT, cpl, alpha, beta0, Diff, csat, Hv, sigmac, lc, vc, Tc, ks, rhos, cpS, alphaS, Mmolar

  substrate = 1.0_rk!0.15_rk!   0.0_rk!    3.0_rk! !1.0_rk  !0.15_rk   
  outer = 1.3_rk    !10.0_rk!   20.0_rk, 1.3_rk

  Ltype = 1

  !solve_T
  if(Maran_flow.eq.1 .and. fixed_Ma.eq.0) then
     solve_T = 1
  else
     solve_T = 0
  end if

  !fluid
  select case(Ltype)
  case(1)   !water
     rho = 9.97e2_rk!(kg/m^3)
     mu =  8.90e-4_rk!(Pa.s)
     kT = 0.58_rk!(W/m/K)
     cpl = 4.18e3_rk!(J/kg/K)
     alpha = kT/rho/cpl     !1.3917e-7_rk  !
     beta0 = -1.657e-4_rk!(N/m/K)  -0.1657dyn/cm/K
     sigmac = 71.97e-3_rk!(N/m)(25C   71.97_rk dyn/cm  
     lc = 1.0e-3_rk!(m)
  case(2)  !octane
     Mmolar = 0.11423_rk  !(kg/mol)
     rho = 6.986e2_rk     !(kg/m^3)
     mu = 5.151e-4_rk  !(Pa.s)
     kT = 0.18_rk!(W/m/K)
     cpl = 255.68_rk*Mmolar  !(J/mol/K)*(kg/mol) = (J/kg/K)
     alpha = kT/rho/cpl     
     beta0 = -0.0935e-3_rk !(N/m/K)  
     sigmac = 21.14e-3_rk!(N/m)(25C   
     lc = 2.0e-3_rk
  case(3)  !hexanol
     Mmolar = 0.10217_rk  !(kg/mol)
     rho = 8.136e2_rk     !(kg/m^3)
     mu = 4.578e-3_rk  !(Pa.s)
     kT = 0.15_rk!(W/m/K)
     !cpl = 255.68*Mmolar  !(J/mol/K)*(kg/mol) = (J/kg/K)
     alpha = 7.84e-8_rk   !(m^2/s)
     beta0 = -8.0e-5_rk !(N/m/K) 
     sigmac = 2.581e-2_rk!(N/m)(25C 
     lc = 1.0e-3_rk
  end select

  !vapor
  select case(Ltype)
  case(1)   !water
     Diff = 2.82e-5_rk!(m^2/s)  !water
     csat = 2.32e-2_rk!(kg/m^3)   !water
     Hv = 2.265e6_rk!(J/kg)   2264.76 (kJ/kg)  
     Hum = 0.5_rk 
  case(2)  !octane
     Diff = 5.0e-6_rk!(m^2/s)  !octane
     csat = 2060.0_rk*Mmolar/8.314_rk/300.0_rk!(kg/m^3)   !octane: c = (Psat*M)/(R*T)   !0.0943kg/m^3
     !(Pa=kg/m/s^2)*(kg/mol)/( (J/mol/K=kg*m^2*s^-2/mol/K) * (K) ) = kg/m^3
     Hv = 41.49e3_rk/Mmolar   !(J/kg)    !3.632e5J/kg
     Hum = 0.0_rk
  case(3)  !hexanol
     Diff = 6.21e-6_rk!(m^2/s)  !octane
     csat = 6.55e-3_rk!(kg/m^3)   !octane: c = (Psat*M)/(R*T)   !0.0943kg/m^3
     !(Pa=kg/m/s^2)*(kg/mol)/( (J/mol/K=kg*m^2*s^-2/mol/K) * (K) ) = kg/m^3
     Hv = 6.03e5_rk  !(J/kg)  
     Hum = 0.0_rk
  end select

  !substrate
  ks = 0.80_rk!(W/m/K)
  rhos = 2.70e3_rk!(kg/m^3)
  cpS = 0.84e3_rk!(J/kg/K)
  alphaS = ks/rhos/cpS
  T_sub = 25.0_rk !     80 !(C)  !25 if not heated

  !particle
  diameterp = 1.0e-8_rk  !1.0e-7_rk  !  (m) particle diameter
  Dp = kboltz*(25.0_rk+273.15_rk)/(6.0_rk*pi*mu*diameterp)   !2.45e-12 (m^2/s)
  !cp0 = 2.5e-4  !(kg/m^3)  !??not used yet
  cp_pack = 5.0_rk
  

  !characteristic
  !sigmac & lc defined in liquid
  Tc = 1.0_rk!Hv*Diff*csat/kT  !2.54(K)
  vc = Diff*csat/rho/lc   !-beta0/mu*Tc   !6.56e-7  

  !dimensionless group
  Re = rho*vc*lc/mu   !6.8e-3_rk
  Ca = mu*vc/sigmac   !8.5e-8_rk
  Kdi = 1.0_rk/Ca       !Kdi: k, used before p as the dimensionless group
  Pe = vc*lc/alpha !4.36e-2_rk
  KBCgroup = rho*vc*lc/Diff/csat
  REH = Hv*diff*csat/kT/Tc   !1.0_rk    !8.55e-3_rk
  beta = Tc/sigmac*beta0 ! -5.87e-3_rk     0.0_rk!
  F0 = alphaS*(lc/vc)/(lc**2)!?
  kR = ks/kT  !relative thermal conductivity
  Pep = lc*vc/Dp    !267.35_rk
  if(solve_T.eq.0) then
     !beta = 0.0_rk
     substrate = 0.0_rk
  end if
  MaN = - beta0*Tc/(mu*vc)

  !change for equations
  Re = Re*Ca
  Oh = Ca  !1.0_rk    !Ca
  Grav = 0.0_rk  !input
  R = 1.0_rk

  T_sub = (T_sub-25.0_rk)/Tc

  open(unit = 10, file = trim(folder)//'parameter_values.dat', status = 'replace')
  write(10,'(A, es14.7)') 'lc =', lc
  write(10,'(A, es14.7)') 'vc =', vc
  write(10,'(A, es14.7)') 'Tc =', Tc
  write(10,'(A, es14.7)') 'Re =', Re/Ca
  write(10,'(A, es14.7)') 'Ca =', Ca
  write(10,'(A, es14.7)') 'Pe =', Pe
  write(10,'(A, es14.7)') 'Pep =', Pep
  write(10,'(A, es14.7)') 'F0 =', F0
  write(10,'(A, es14.7)') 'kR =', kR
  write(10,'(A, es14.7)') 'Grav =', Grav
  write(10,'(A, es14.7)') 'Ma =', MaN
  write(10,'(A, es14.7)') 'substrate =', substrate
  write(10,'(A, es14.7)') 'outer =', outer
  write(10,'(A, es14.7)') 'T_sub =', T_sub
  write(10,'(A, es14.7)') 'diameterp =', diameterp
  write(10,'(A)') ' '

  write(10,'(A, es14.7)') 'KBCgroup =', KBCgroup
  write(10,'(A, es14.7)') 'REH =', REH
  write(10,'(A, es14.7)') 'beta =', beta
  write(10,'(A)') ' '

  write(10,'(A, i8)') 'NStrans =', NStrans
  write(10,'(A, i8)') 'Inert =', Inert
  write(10,'(A, i8)') 'Capil =', Capil
  write(10,'(A, i8)') 'Viscous =', Viscous
  write(10,'(A, i8)') 'GravI =', GravI
  write(10,'(A, i8)') 'Ttime =', Ttime
  write(10,'(A, i8)') 'Tconv =', Tconv
  write(10,'(A, i8)') 'Tdiff =', Tdiff
  write(10,'(A, i8)') 'TtimeS =', TtimeS
  write(10,'(A, i8)') 'TdiffS =', TdiffS
  write(10,'(A, i8)') 'no_vapor =', no_vapor
  write(10,'(A, i8)') 'uniflux =', uniflux
  write(10,'(A, i8)') 'solve_T =', solve_T
  write(10,'(A, i8)') 'Maran_flow =', Maran_flow
  write(10,'(A, i8)') 'solve_cp =', solve_cp

  write(10,'(A)') ' '

  write(10,'(A, i8)') 'NEL = ', NEL
  write(10,'(A, i8)') 'NEM = ', NEM
  write(10,'(A, i8)') 'NEV = ', NEV
  write(10,'(A, i8)') 'NES = ', NES
  write(10,'(A, i8)') 'NEM_alge = ', NEM_alge
  if(no_vapor.eq.1) then
     write(10,'(A, i8)') 'rad_ele =', NEL*2+NES
     write(10,'(A, i8)') 'total elements: ', NEL**2+2*NEL*NEM+NES*(NEM+NEL+NEV)
  else ! no_vapor=0
     write(10,'(A, i8)') 'rad_ele =', NEL*2+NES+NEV
     ! write(10,'(A, i8)') 'total elements: '
  end if

  close(10)


  write(*,*) 'KBCgroup =', KBCgroup
  !write(*,*) 'REH =', REH
  !write(*,*) 'beta =', beta

  return
end subroutine parameter_values
