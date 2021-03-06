subroutine SI_in_sj(m,i,j, sj, LNVar, LNOPP, id)                     !adding SI to sj
  use kind
  use data
  use Ldata
  use NOP_mod, only: gaussian_quadrature_1d

  implicit none

  integer(kind=ik), intent(in):: m,i,j, LNVar, LNOPP(9), id
  real(kind=rk), intent(out):: sj(LNVar, LNVar)
  !real(kind=rk):: flux

  integer(kind=ik):: k, ipp, jpp
  real(kind=rk):: intRsi_r_S(Ng), intRsi_z_S(Ng), intReta_r_S(Ng), intReta_z_S(Ng)
  real(kind=rk):: intRsi_u_S(Ng), intRsi_v_S(Ng), intRsi_c_S(Ng), &
       intRu_r_S(Ng), intRu_z_S(Ng), intRv_r_S(Ng), intRv_z_S(Ng), &
       intRu_T_S(Ng), intRv_T_S(Ng)
  real(kind=rk):: intRt_r_S(Ng), intRt_z_S(Ng), intRt_c_S(Ng), intRt_T_S(Ng)
  real(kind=rk):: intRm_r_S(Ng), intRm_z_S(Ng), intRm_cp_S(Ng), intRm_gamma_S(Ng), intRm_c_S(Ng)
  real(kind=rk):: intRms_r_S(Ng), intRms_z_S(Ng), intRms_u_S(Ng), intRms_v_S(Ng), intRms_cp_S(Ng), intRms_gamma_S(Ng)
  real(kind=rk):: Rms1r, Rms2r, Rms31r, Rms32r, Rms4r, Rms5r, Rms6r, Rms1z, Rms2z, Rms31z, Rms32z, Rms4z, Rms5z, Rms6z, &
       Rms31u, Rms32u, Rms5u, Rms6u, Rms31v, Rms32v, Rms5v, Rms6v, &
       Rms2cp, Rms1gamma, Rms2gamma, Rms31gamma, Rms32gamma, Rms4gamma, Rms5gamma, Rms6gamma, &
       twoHterm_r, twoHterm_z, Rms5term_r, Rms5term_z, Rms32term
  intRsi_r_S(:) = 0.0_rk
  intRsi_z_S(:) = 0.0_rk
  intReta_r_S(:) = 0.0_rk
  intReta_z_S(:) = 0.0_rk
  intRsi_u_S(:) = 0.0_rk
  intRsi_v_S(:) = 0.0_rk
  intRsi_c_S(:) = 0.0_rk
  intRu_r_S(:) = 0.0_rk
  intRu_z_S(:) = 0.0_rk
  intRv_r_S(:) = 0.0_rk
  intRv_z_S(:) = 0.0_rk
  intRu_T_S(:) = 0.0_rk
  intRv_T_S(:) = 0.0_rk
  intRt_r_S(:) = 0.0_rk
  intRt_z_S(:) = 0.0_rk
  intRt_c_S(:) = 0.0_rk
  intRt_T_S(:) = 0.0_rk
  intRm_r_S(:) = 0.0_rk
  intRm_z_S(:) = 0.0_rk
  intRm_cp_S(:) = 0.0_rk
  intRm_gamma_S(:) = 0.0_rk
  intRm_c_S(:) = 0.0_rk
  intRms_r_S(:) = 0.0_rk
  intRms_z_S(:) = 0.0_rk
  intRms_u_S(:) = 0.0_rk
  intRms_v_S(:) = 0.0_rk
  intRms_cp_S(:) = 0.0_rk
  intRms_gamma_S(:) = 0.0_rk

  ! Rms1r = 0.0_rk
  ! Rms2r = 0.0_rk
  ! Rms3r = 0.0_rk
  ! Rms4r = 0.0_rk
  ! Rms1z = 0.0_rk
  ! Rms2z = 0.0_rk
  ! Rms3z = 0.0_rk
  ! Rms4z = 0.0_rk
  ! Rms1gamma = 0.0_rk
  ! Rms2gamma = 0.0_rk
  ! Rms3gamma = 0.0_rk
  ! Rms4gamma = 0.0_rk
  


  !axis
  if( BCflagN( globalNM(m,i),1 ).eq.1 .and. BCflagN( globalNM(m,j),1 ).eq.1 ) then
     ipp = i - 6  !phix_1d(k,ipp)
     jpp = j - 6  !phix_1d(k,jpp)
     do k = 1, Ng, 1    !three gausspoints
        intRsi_r_S(k) = phix_1d(k,ipp) * fsi_size(m) / ( rsi_down(k,id)**2 + zsi_down(k,id)**2 ) &
             * 2.0_rk*rsi_down(k,id) * phix_1d(k,jpp)
        intRsi_z_S(k) = phix_1d(k,ipp) * fsi_size(m) / ( rsi_down(k,id)**2 + zsi_down(k,id)**2 ) &
             * 2.0_rk*zsi_down(k,id) * phix_1d(k,jpp)
     end do
     sj(LNOPP(i)+Nr,LNOPP(j)+Nr) = sj(LNOPP(i)+Nr,LNOPP(j)+Nr) - M1*gaussian_quadrature_1d(intRsi_r_S)
     sj(LNOPP(i)+Nr,LNOPP(j)+Nz) = sj(LNOPP(i)+Nr,LNOPP(j)+Nz) - M1*gaussian_quadrature_1d(intRsi_z_S)
  end if

  !base
  if( ( BCflagE(m,2).eq.1 .or. BCflagE(m,2).eq.11 ) &
     .and. (BCflagN(globalNM(m,i),2).eq.1 .or. BCflagN(globalNM(m,i),2).eq.3) ) then

     if(BCflagE(m,2).eq.1) then
        ipp = i  !phix_1d(k,ipp)
        
        if( BCflagN( globalNM(m,j),2 ).eq.1 .or. BCflagN( globalNM(m,j),2 ).eq.3 ) then

           if(BCflagE(m,2).eq.1) then
              jpp = j  !phix_1d(k,jpp)
           end if

           do k = 1, Ng, 1    !three gausspoints
              intRsi_r_S(k) = phix_1d(k,ipp) * fsi_size(m) / ( rsi_left(k,id)**2 + zsi_left(k,id)**2 ) &
                   * 2.0_rk*rsi_left(k,id) * phix_1d(k,jpp)
              intRsi_z_S(k) = phix_1d(k,ipp) * fsi_size(m) / ( rsi_left(k,id)**2 + zsi_left(k,id)**2 ) &
                   * 2.0_rk*zsi_left(k,id) * phix_1d(k,jpp)

           end do

           if(((VE(m).eq.1 .and. rowNM(m).gt.NEV1(1)).or.(VE(m).eq.5 .and. rowNM(m).gt.NEV1(1)+NEM+NEL))&
                .or. base_mesh.eq.0 ) then!elements farfrom drop
              sj(LNOPP(i)+Nr,LNOPP(j)+Nr) = sj(LNOPP(i)+Nr,LNOPP(j)+Nr) - M1*gaussian_quadrature_1d(intRsi_r_S)
              sj(LNOPP(i)+Nr,LNOPP(j)+Nz) = sj(LNOPP(i)+Nr,LNOPP(j)+Nz) - M1*gaussian_quadrature_1d(intRsi_z_S)
           else
              sj(LNOPP(i)+Nr,LNOPP(j)+Nr) = gaussian_quadrature_1d(intRsi_r_S)
              sj(LNOPP(i)+Nr,LNOPP(j)+Nz) = gaussian_quadrature_1d(intRsi_z_S)
           end if

        else
           if( (.not.( (VE(m).eq.1 .and. rowNM(m).gt.NEV1(1)).or.&
                (VE(m).eq.5 .and. rowNM(m).gt.NEV1(1)+NEM+NEL) ) ) .and. &!.not.(elements far from drop)
                base_mesh.eq.1 ) then
              do k = Nr, Nz
                 sj(LNOPP(i)+Nr,LNOPP(j)+k) = 0.0_rk    !in SI, d Rsi/ d xj = 0.0_rk other than j = 1,2,3. so replace original volume integral with 0.0
              end do
           end if  !for j = 1,2,3 on the free surface or not
        end if

     else !BCflagE(m,2) = 11
        if( .not. ( no_vapor.eq.1 .and. VN(globalNM(m,i)).eq.1 ) ) then
           do k = Nr, Nz
              sj(LNOPP(i)+Nr,LNOPP(j)+k) = 0.0_rk    
           end do
        end if
     end if

  end if   !base of region 1&2 & 7&8

  if( ( BCflagE(m,2).eq.2 .or. BCflagE(m,2).eq.21 ) &
      .and. (BCflagN(globalNM(m,i),2).eq.2 .or. BCflagN(globalNM(m,i),2).eq.3) ) then

     if(BCflagE(m,2).eq.2) then
        ipp = i/3 + 1  !phix_1d(k,ipp)
        
        if( BCflagN( globalNM(m,j),2 ).eq.2 .or. BCflagN( globalNM(m,j),2 ).eq.3 )  then

           if(BCflagE(m,2).eq.2) then
              jpp = j/3 + 1  !phix_1d(k,jpp)
           end if

           do k = 1, Ng, 1    !three gausspoints
              intReta_r_S(k) = phix_1d(k,ipp) * geta_size(m) / &
                   ( reta_left(k,id)**2 + zeta_left(k,id)**2 ) * 2.0_rk*reta_left(k,id) * phix_1d(k,jpp)
              intReta_z_S(k) = phix_1d(k,ipp) * geta_size(m) / &
                   ( reta_left(k,id)**2 + zeta_left(k,id)**2 ) * 2.0_rk*zeta_left(k,id) * phix_1d(k,jpp)
           end do
           if(base_mesh.eq.0) then
              ! sj(LNOPP(i)+Nz,LNOPP(j)+Nr) = sj(LNOPP(i)+Nz,LNOPP(j)+Nr) - M1*gaussian_quadrature_1d(intReta_r_S)
              ! sj(LNOPP(i)+Nz,LNOPP(j)+Nz) = sj(LNOPP(i)+Nz,LNOPP(j)+Nz) - M1*gaussian_quadrature_1d(intReta_z_S)
           else
              sj(LNOPP(i)+Nz,LNOPP(j)+Nr) = gaussian_quadrature_1d(intReta_r_S)
              sj(LNOPP(i)+Nz,LNOPP(j)+Nz) = gaussian_quadrature_1d(intReta_z_S)
           end if

        else
           if(base_mesh.eq.1) then
              do k = Nr, Nz
                 sj(LNOPP(i)+Nz,LNOPP(j)+k) = 0.0_rk  !in SI, d Reta/ d xj = 0.0_rk other than j = 1,4,7. so replace original volume integral with 0.0
              end do
           end if
        end if  !for j = 1,4,7 on the free surface or not
        
     else  !BCflagE(m,2) = 21
        do k = Nr, Nz
           sj(LNOPP(i)+Nz,LNOPP(j)+k) = 0.0_rk    
        end do
     end if

  end if   !base of region 3


  
  !base, substrate adsorption
  if(s_mode.eq.0 .and. solve_cp.eq.1 .and. sub_adsp.eq.1) then
     
     if( ( BCflagE(m,2).eq.1 .and. VE(m).eq.0 ) &
          .and. (BCflagN(globalNM(m,i),2).eq.1 .or. BCflagN(globalNM(m,i),2).eq.3) ) then
        ipp = i  !phix_1d(k,ipp)
        if( BCflagN( globalNM(m,j),2 ).eq.1 .or. BCflagN( globalNM(m,j),2 ).eq.3 ) then
           jpp = j  !phix_1d(k,jpp)
           do k = 1, Ng, 1    !three gausspoints
              intRm_r_S(k) = phi_1d(k,ipp) *cpintfac_left(k,id) * &
                   ( phi_1d(k,jpp)*rsi_left(k,id) + rintfac_left(k,id)*phix_1d(k,jpp) ) &
                   *abs(rsi_left(k,id))/rsi_left(k,id)
              intRm_cp_S(k) = phi_1d(k,ipp)*phi_1d(k,jpp)*rintfac_left(k,id)*rsi_left(k,id)&
                   *abs(rsi_left(k,id))/rsi_left(k,id)
           end do
           sj(LNOPP(i)+Ncp,LNOPP(j)+Nr) = sj(LNOPP(i)+Ncp,LNOPP(j)+Nr) + &
                Da_sub *gaussian_quadrature_1d(intRm_r_S)
           sj(LNOPP(i)+Ncp,LNOPP(j)+Ncp) = sj(LNOPP(i)+Ncp,LNOPP(j)+Ncp) + &
                Da_sub *gaussian_quadrature_1d(intRm_cp_S)
        end if !j
     end if ! m element on base, i node on base

     if( BCflagE(m,2).eq.2  &
          .and. (BCflagN(globalNM(m,i),2).eq.2 .or. BCflagN(globalNM(m,i),2) .eq.3) ) then
        ipp = i/3 + 1  !phix_1d(k,ipp)
        if( BCflagN( globalNM(m,j),2 ).eq.2 .or. BCflagN( globalNM(m,j),2 ).eq.3 ) then
           jpp = j/3 + 1  !phix_1d(k,ipp)
           do k = 1, Ng, 1    !three gausspoints
              intRm_r_S(k) = phi_1d(k,ipp) *cpintfac_left(k,id) * &
                   ( phi_1d(k,jpp)*reta_left(k,id) + rintfac_left(k,id)*phix_1d(k,jpp) )&
                   *abs(reta_left(k,id))/reta_left(k,id)
              intRm_cp_S(k) = phi_1d(k,ipp)*phi_1d(k,jpp)*rintfac_left(k,id)*reta_left(k,id)&
                   *abs(reta_left(k,id))/reta_left(k,id)
           end do
           sj(LNOPP(i)+Ncp,LNOPP(j)+Nr) = sj(LNOPP(i)+Ncp,LNOPP(j)+Nr) + &
                Da_sub *gaussian_quadrature_1d(intRm_r_S)
           sj(LNOPP(i)+Ncp,LNOPP(j)+Ncp) = sj(LNOPP(i)+Ncp,LNOPP(j)+Ncp) + &
                Da_sub *gaussian_quadrature_1d(intRm_cp_S)
        end if!j
     end if !m & i

  end if  !s_mode=0 and solve_cp=1 and sub_adsp=1


  !outer surface
  if( BCflagN( globalNM(m,i),4 ).ne.0 .and. BCflagN( globalNM(m,j),4 ).ne.0 ) then
     ipp = i/3  !phix_1d(k,l)
     jpp = j/3  !phix_1d(k,l)
     do k = 1, Ng, 1    !three gausspoints

intReta_r_S(k) = phix_1d(k,ipp) * geta_size(m) / ( reta_right(k,id)**2 + zeta_right(k,id)**2 ) &
     * 2.0_rk*reta_right(k,id) * phix_1d(k,jpp)
intReta_z_S(k) = phix_1d(k,ipp) * geta_size(m) / ( reta_right(k,id)**2 + zeta_right(k,id)**2 ) &
     * 2.0_rk*zeta_right(k,id) * phix_1d(k,jpp)

     end do
     sj(LNOPP(i)+Nz,LNOPP(j)+Nr) = sj(LNOPP(i)+Nz,LNOPP(j)+Nr) - M1*gaussian_quadrature_1d(intReta_r_S)
     sj(LNOPP(i)+Nz,LNOPP(j)+Nz) = sj(LNOPP(i)+Nz,LNOPP(j)+Nz) - M1*gaussian_quadrature_1d(intReta_z_S)!?should be M2, doesn't matter because M1=M2=0

  end if


!free suface, decouple elliptic mesh inside and outside the drop
if(mesh_decouple.eq.1) then

  if( BCflagN( globalNM(m,i),3 ).eq.1 .or. BCflagN( globalNM(m,i),3 ).eq.3 ) then
     if(VE(m).eq.0) then
        ipp = i/3  !phix_1d(k,ipp)
     ! else
     !    ipp = i/3 + 1  !phix_1d(k,ipp)
     ! end if
        
        if( BCflagN( globalNM(m,j),3 ).eq.1 .or.BCflagN( globalNM(m,j),3 ).eq.3 ) then
           ! if(VE(m).eq.0) then
              jpp = j/3  !phix_1d(k,jpp)
           ! else
           !    jpp = j/3 + 1  !phix_1d(k,jpp)
           ! end if

           do k = 1, Ng, 1    !three gausspoints

intReta_r_S(k) = phix_1d(k,ipp) * geta_size(m) / SQr2z2(k,id) * 2.0_rk*reta_right(k,id) * phix_1d(k,jpp)
intReta_z_S(k) = phix_1d(k,ipp) * geta_size(m) / SQr2z2(k,id) * 2.0_rk*zeta_right(k,id) * phix_1d(k,jpp)

           end do

           !directly replace Reta_V with Reta_S, viz SI of elliptic mesh
           sj(LNOPP(i)+Nz,LNOPP(j)+Nr) = gaussian_quadrature_1d(intReta_r_S)
           sj(LNOPP(i)+Nz,LNOPP(j)+Nz) = gaussian_quadrature_1d(intReta_z_S)

        else   ! j not on boundary
           do k = Nr, Nz
              sj(LNOPP(i)+Nz,LNOPP(j)+k) = 0.0_rk  
              !in SI, d Reta/ d xj = 0.0_rk other than j = 1,4,7. so replace original volume integral with 0.0
           end do
        end if  !for j = 1,4,7 on the free surface or not

     else  !VE = 1
        do k = Nr, Nz
           sj(LNOPP(i)+Nz,LNOPP(j)+k) = 0.0_rk  !no SI from vapor element
        end do
     end if
        
  end if  !for i on the free surface

end if  !for mesh_decouple



if(s_mode.eq.0) then

   !free suface from drop, KBC part1 & traction BC & evaporation cooling part1
   if( ( BCflagN( globalNM(m,i),3 ).eq.1 .or. BCflagN( globalNM(m,i),3 ).eq.3 ).and. VE(m).eq.0) then
      ipp = i/3  !phix_1d(k,l)

      if( BCflagN( globalNM(m,j),3 ).eq.1 .or. BCflagN( globalNM(m,j),3 ).eq.3 ) then
         jpp = j/3  !phix_1d(k,l)
     
         do k = 1, Ng, 1    !three gausspoints

dSdr(k,id) = phi_1d(k,jpp)*sqrt(SQr2z2(k,id)) + rintfac_right(k,id)/sqrt(SQr2z2(k,id))*reta_right(k,id)*phix_1d(k,jpp)
            
dSdz(k,id) = rintfac_right(k,id)/sqrt(SQr2z2(k,id))*zeta_right(k,id)*phix_1d(k,jpp)

if(no_vapor.eq.1) then  !flux: flux( angle_c,rintfac_right(k,id) ) 
   
   !KBC2
   intRsi_r_S(k) =  phi_1d(k,ipp)* flux(k,id) *&
        dSQdr(k,id) *phix_1d(k,jpp)* rintfac_right(k,id) + &
        
        phi_1d(k,ipp)* ( flux(k,id)+ flux_r(k,id)*rintfac_right(k,id) ) *&
        SQr2z2(k,id)**0.5_rk *phi_1d(k,jpp)


   intRsi_z_S(k) = phi_1d(k,ipp)* flux(k,id) *&
        dSQdz(k,id) *phix_1d(k,jpp)* rintfac_right(k,id)


   !evaporation cooling 2
   if(solve_T.eq.1) then
      intRt_r_S(k) = REH* intRsi_r_S(k)
      intRt_z_S(k) = REH* intRsi_z_S(k)
   end if

   
   !KBC2 with uniflux: flux = 1.0_rk, only apply in KBC & accumulation: Rsi&Rm
   if(uniflux.eq.1) then
      intRsi_r_S(k) = intRsi_r_S(k) / flux(k,id)
      intRsi_z_S(k) = intRsi_z_S(k) / flux(k,id)
   end if

   
   !particle accumulation 2
   if(solve_cp.eq.1) then
      intRm_r_S(k) = intRsi_r_S(k) * cpintfac_right(k,id)
      intRm_z_S(k) = intRsi_z_S(k) * cpintfac_right(k,id)
      
      intRm_cp_S(k) = phi_1d(k,ipp)* flux(k,id) *phi_1d(k,jpp)* dS(k,id)
      if(uniflux.eq.1) intRm_cp_S(k) = intRm_cp_S(k) / flux(k,id)

      if(surf_adsp.eq.1) then
intRm_r_S(k) = intRm_r_S(k) - KBCgroup* phi_1d(k,ipp)* &
     adsp_rate(k,id) * ( phi_1d(k,jpp)*SQr2z2(k,id)**0.5 + &
     rintfac_right(k,id)* dSQdr(k,id) *phix_1d(k,jpp) )

intRm_z_S(k) = intRm_z_S(k) - KBCgroup* phi_1d(k,ipp)* &
     adsp_rate(k,id) *rintfac_right(k,id)* dSQdz(k,id) *phix_1d(k,jpp)

intRm_cp_S(k) = intRm_cp_S(k) - KBCgroup* phi_1d(k,ipp)*( Da_surf1*gammaintfac(k,id) + Da_surf2) *phi_1d(k,jpp)*dS(k,id)

intRm_gamma_S(k) = -KBCgroup* phi_1d(k,ipp)* Da_surf1* phi_1d(k,jpp) *cpintfac_right(k,id) *dS(k,id)
      end if  !surf_adsp

   end if  !solve_cp
   
end if   !no_vapor=1


!KBC1
intRsi_r_S(k) = intRsi_r_S(k) + KBCgroup* ( phi_1d(k,ipp)* &
     phix_1d(k,jpp)* vintfac_right(k,id) *rintfac_right(k,id) + &
     
     phi_1d(k,ipp)*( -zeta_right(k,id)* uintfac_right(k,id) + &
     reta_right(k,id)* vintfac_right(k,id) )*phi_1d(k,jpp) )

if(NStrans.eq.1) &
intRsi_r_S(k) = intRsi_r_S(k) + KBCgroup* ( phi_1d(k,ipp)*( zeta_right(k,id)*CTJ/dt*phi_1d(k,jpp) - &
     phix_1d(k,jpp)* zdotintfac_right(k,id) ) *rintfac_right(k,id) + &
     
     phi_1d(k,ipp)*( zeta_right(k,id)* rdotintfac_right(k,id)  - &
     reta_right(k,id)* zdotintfac_right(k,id) )*phi_1d(k,jpp) )

! intRsi_r_S(k) = intRsi_r_S(k) + KBCgroup* ( phi_1d(k,ipp)*( zeta_right(k,id)*CTJ/dt*phi_1d(k,jpp) + &
!      phix_1d(k,jpp)*( vintfac_right(k,id) - zdotintfac_right(k,id) ) )*rintfac_right(k,id) + &
     
!      phi_1d(k,ipp)*( -zeta_right(k,id)*( uintfac_right(k,id) - rdotintfac_right(k,id) ) + &
!      reta_right(k,id)*( vintfac_right(k,id) - zdotintfac_right(k,id) ) )*phi_1d(k,jpp) )


intRsi_z_S(k) = intRsi_z_S(k) + KBCgroup* ( phi_1d(k,ipp)*( -phix_1d(k,jpp)*&
      uintfac_right(k,id) )*rintfac_right(k,id) )

if(NStrans.eq.1) &
intRsi_z_S(k) = intRsi_z_S(k) + KBCgroup* ( phi_1d(k,ipp)* &
     ( phix_1d(k,jpp)*rdotintfac_right(k,id) - reta_right(k,id)*CTJ/dt*phi_1d(k,jpp) )*rintfac_right(k,id) )

! intRsi_z_S(k) = intRsi_z_S(k) + KBCgroup* ( phi_1d(k,ipp)*( -phix_1d(k,jpp)*&
!      ( uintfac_right(k,id) - rdotintfac_right(k,id) ) - &
!      reta_right(k,id)*CTJ/dt*phi_1d(k,jpp) )*rintfac_right(k,id) )


intRsi_u_S(k) = intRsi_u_S(k) + KBCgroup* ( -phi_1d(k,jpp)*phi_1d(k,ipp)*zeta_right(k,id)*rintfac_right(k,id) )

intRsi_v_S(k) = intRsi_v_S(k) + KBCgroup* ( phi_1d(k,jpp)*phi_1d(k,ipp)*reta_right(k,id)*rintfac_right(k,id) )


!traction BC
intRu_r_S(k) = ( phix_1d(k,ipp)*phix_1d(k,jpp) / SQr2z2(k,id) - &
     reta_right(k,id)*phix_1d(k,ipp) / SQr2z2(k,id)**2 * 2.0_rk*reta_right(k,id)*phix_1d(k,jpp) - &
     1.0_rk/rintfac_right(k,id)**2 *phi_1d(k,jpp)*phi_1d(k,ipp) ) * dS(k,id) + &
     
     ( reta_right(k,id)*phix_1d(k,ipp) / SQr2z2(k,id) + phi_1d(k,ipp)/rintfac_right(k,id) )&
     *( phi_1d(k,jpp)*SQr2z2(k,id)**0.5 + &
     rintfac_right(k,id) / SQr2z2(k,id)**0.5 *reta_right(k,id)*phix_1d(k,jpp) )


if(Maran_flow.eq.1) then
   if(fixed_Ma.eq.0) then
   intRu_r_S(k) = intRu_r_S(k) - phi_1d(k,ipp) *beta *Teta_right(k,id)* ( &
        -SQr2z2(k,id)**(-1.5) *reta_right(k,id)**2 *phix_1d(k,jpp) *rintfac_right(k,id) + &
        SQr2z2(k,id)**(-0.5) *&
        ( phix_1d(k,jpp) *rintfac_right(k,id) + reta_right(k,id) *phi_1d(k,jpp) ) )
   else !fixed_Ma.eq.1
   intRu_r_S(k) = intRu_r_S(k) + Ca*MaN*phi_1d(k,ipp)* &
        ( phix_1d(k,jpp) *rintfac_right(k,id) + reta_right(k,id) *phi_1d(k,jpp) )
   end if !fixed_Ma
end if  !Maran_flow

intRu_z_S(k) = -reta_right(k,id)*phix_1d(k,ipp)*SQr2z2(k,id)**(-1.5) *&
     2.0_rk*zeta_right(k,id)*phix_1d(k,jpp)*rintfac_right(k,id) + &
     
     ( reta_right(k,id)*phix_1d(k,ipp) / SQr2z2(k,id) + phi_1d(k,ipp)/rintfac_right(k,id) ) *&
     rintfac_right(k,id) / SQr2z2(k,id)**0.5 *zeta_right(k,id)*phix_1d(k,jpp)

if(Maran_flow.eq.1 .and. fixed_Ma.eq.0) then
   intRu_z_S(k) = intRu_z_S(k) + phi_1d(k,ipp) *beta *Teta_right(k,id)* &
        SQr2z2(k,id)**(-1.5)* &
        zeta_right(k,id) *phix_1d(k,jpp) *reta_right(k,id) *rintfac_right(k,id) 

   intRu_T_S(k) = -phi_1d(k,ipp) *beta *phix_1d(k,jpp) * &
        SQr2z2(k,id)**(-0.5)* reta_right(k,id) *rintfac_right(k,id)
end if   !Maran_flow.eq.1 .and. fixed_Ma.eq.0


intRv_r_S(k) = zeta_right(k,id)*phix_1d(k,ipp)*&
     ( -SQr2z2(k,id)**(-1.5) *reta_right(k,id)*phix_1d(k,jpp)*rintfac_right(k,id) + &
     SQr2z2(k,id)**(-0.5) *phi_1d(k,jpp) )

 if(Maran_flow.eq.1) then
   if(fixed_Ma.eq.0) then
   intRv_r_S(k) = intRv_r_S(k) - phi_1d(k,ipp)*beta*Teta_right(k,id)*zeta_right(k,id)* ( &
      -SQr2z2(k,id)**(-1.5) *reta_right(k,id) *phix_1d(k,jpp) *rintfac_right(k,id) + &
      SQr2z2(k,id)**(-0.5) * phi_1d(k,jpp) )
   else 
   intRv_r_S(k) = intRv_r_S(k) + Ca*MaN*phi_1d(k,ipp)*zeta_right(k,id)*phi_1d(k,jpp)
   end if !fixed_Ma
end if  !Maran_flow
      

intRv_z_S(k) = phix_1d(k,ipp)*phix_1d(k,jpp) / sqrt( SQr2z2(k,id) ) *rintfac_right(k,id) - &
     zeta_right(k,id)**2 *phix_1d(k,ipp)*phix_1d(k,jpp) / SQr2z2(k,id)**1.5 *rintfac_right(k,id)

 if(Maran_flow.eq.1) then
    if(fixed_Ma.eq.0) then
       intRv_z_S(k) = intRv_z_S(k) - phi_1d(k,ipp) *beta *Teta_right(k,id)* ( &
      phix_1d(k,jpp) *SQr2z2(k,id)**(-0.5) - &
      SQr2z2(k,id)**(-1.5)*zeta_right(k,id)**2 *phix_1d(k,jpp) ) &
      *rintfac_right(k,id) 

       intRv_T_S(k) = -phi_1d(k,ipp) *beta *phix_1d(k,jpp) * &
            SQr2z2(k,id)**(-0.5)* zeta_right(k,id) *rintfac_right(k,id)
    else  !fixed_Ma.eq.1
       intRv_z_S(k) = intRv_z_S(k) + Ca*MaN*phi_1d(k,ipp)*phix_1d(k,jpp)*rintfac_right(k,id)
    end if !fixed_Ma
 end if !Maran_flow
       

!adsorbed particle
 if(solve_cp.eq.1 .and. surf_adsp.eq.1) then
    !Rms1, temporal change
    Rms1r = phi_1d(k,ipp)*gammadot(k,id)*dSdr(k,id) - &
         phi_1d(k,ipp)*gammaeta(k,id)*( &
         (CTJ/dt*phi_1d(k,jpp)*reta_right(k,id) )*rintfac_right(k,id)*SQr2z2(k,id)**(-0.5) + &
         Rms1term(k,id)*phi_1d(k,jpp)*SQr2z2(k,id)**(-0.5) - &
         Rms1term(k,id)*rintfac_right(k,id)*SQr2z2(k,id)**(-1.5)*reta_right(k,id)*phix_1d(k,jpp) )

    Rms1z = phi_1d(k,ipp)*gammadot(k,id)*dSdz(k,id) &
         + phi_1d(k,ipp)*( &
         - gammaeta(k,id)*( CTJ/dt*phi_1d(k,jpp) *zeta_right(k,id) + &
         zdotintfac_right(k,id)*phix_1d(k,jpp) )*rintfac_right(k,id)/sqrt(SQr2z2(k,id)) + &
         gammaeta(k,id)*Rms1term(k,id)*rintfac_right(k,id)&
         *SQr2z2(k,id)**(-1.5_rk)*zeta_right(k,id)*phix_1d(k,jpp)  & 
         )

    Rms1gamma = phi_1d(k,ipp)*rintfac_right(k,id)*CTJ/dt*phi_1d(k,jpp) *sqrt( SQr2z2(k,id) ) &
         - phi_1d(k,ipp)*rintfac_right(k,id)*phix_1d(k,jpp)*Rms1term(k,id)/sqrt( SQr2z2(k,id) ) 

    !Rms2, adsorption
    Rms2r = -phi_1d(k,ipp)*adsp_rate(k,id)*&
         ( phi_1d(k,jpp)*sqrt(SQr2z2(k,id)) + rintfac_right(k,id)*dSQdr(k,id)*phix_1d(k,jpp) )

    Rms2z = -phi_1d(k,ipp)*adsp_rate(k,id)*rintfac_right(k,id)*dSQdz(k,id)*phix_1d(k,jpp)
    
    Rms2cp = -phi_1d(k,ipp)*( Da_surf1*gammaintfac(k,id) + Da_surf2 )*phi_1d(k,jpp)*dS(k,id)
    
    Rms2gamma = -phi_1d(k,ipp)*Da_surf1*phi_1d(k,jpp)*cpintfac_right(k,id)*dS(k,id)
    
    !Rms31, convection
    Rms31r = phi_1d(k,ipp)* ( &
         phi_1d(k,jpp)* gammaeta(k,id)*ureandvze(k,id) /sqrt(SQr2z2(k,id)) + &
         rintfac_right(k,id)* gammaeta(k,id)*uintfac_right(k,id)*phix_1d(k,jpp) /sqrt(SQr2z2(k,id)) - &
         rintfac_right(k,id)* gammaeta(k,id)*ureandvze(k,id) *SQr2z2(k,id)**(-1.5)*reta_right(k,id)*phix_1d(k,jpp)  )

    Rms31z = phi_1d(k,ipp)*rintfac_right(k,id)*( &
         gammaeta(k,id)*vintfac_right(k,id)*phix_1d(k,jpp)/sqrt(SQr2z2(k,id)) -&
         
         gammaeta(k,id)*ureandvze(k,id) / SQr2z2(k,id)**(1.5_rk) *zeta_right(k,id)*phix_1d(k,jpp)  )

    Rms31u = phi_1d(k,ipp)*rintfac_right(k,id)*( &
         gammaeta(k,id)*phi_1d(k,jpp)*reta_right(k,id)/sqrt(SQr2z2(k,id)) )

    Rms31v = phi_1d(k,ipp)*rintfac_right(k,id)*(gammaeta(k,id)*phi_1d(k,jpp)*zeta_right(k,id) / sqrt(SQr2z2(k,id)) )

    Rms31gamma = phi_1d(k,ipp)*rintfac_right(k,id)/sqrt(SQr2z2(k,id)) *phix_1d(k,jpp)*ureandvze(k,id) 
    
    !Rms32, stretching
    Rms32term = phix_1d(k,ipp)*gammaintfac(k,id) + phi_1d(k,ipp)*gammaeta(k,id)
    Rms32r = -Rms32term*( uintfac_right(k,id)*phix_1d(k,jpp)/sqrt(SQr2z2(k,id))*rintfac_right(k,id) &
         + ureandvze(k,id)*( -SQr2z2(k,id)**(-1.5)*reta_right(k,id)*phix_1d(k,jpp)*rintfac_right(k,id) &
         + phi_1d(k,jpp)/sqrt(SQr2z2(k,id)) ) )

    Rms32z = -Rms32term*rintfac_right(k,id)*( vintfac_right(k,id)*phix_1d(k,jpp)/sqrt(SQR2z2(k,id)) &
         - ureandvze(k,id)/SQr2z2(k,id)**1.5*zeta_right(k,id)*phix_1d(k,jpp) )
    
    Rms32u = -Rms32term*phi_1d(k,jpp)*reta_right(k,id)/sqrt(SQr2z2(k,id))*rintfac_right(k,id)
    
    Rms32v = -Rms32term*phi_1d(k,jpp)*zeta_right(k,id)/sqrt(SQr2z2(k,id))*rintfac_right(k,id)
    
    Rms32gamma = -( phix_1d(k,ipp)*phi_1d(k,jpp) + phi_1d(k,ipp)*phix_1d(k,jpp) ) * &
         ureandvze(k,id)/sqrt(SQr2z2(k,id))*rintfac_right(k,id)
    
    !Rms4, diffusion
    Rms4r = 1.0_rk/Pep_surf*phix_1d(k,ipp)*gammaeta(k,id)*( &
         -2.0_rk/SQr2z2(k,id)**2*reta_right(k,id)*phix_1d(k,jpp)*dS(k,id) + dSdr(k,id)/SQr2z2(k,id) )

    Rms4z = 1.0_rk/Pep_surf*phix_1d(k,ipp)*gammaeta(k,id)*( &
         -2.0_rk/SQr2z2(k,id)**2*zeta_right(k,id)*phix_1d(k,jpp)*dS(k,id) + dSdz(k,id)/SQr2z2(k,id))
    
    Rms4gamma = 1.0_rk/Pep_surf*phix_1d(k,ipp)*phix_1d(k,jpp)/SQr2z2(k,id)*dS(k,id)

    !Rms5, dilatation
    twoHterm_r = -2.0_rk*sin(angle_c)/R/SQr2z2(k,id)**1.5*reta_right(k,id)*phix_1d(k,jpp)
    Rms5term_r = (CTJ/dt*phi_1d(k,jpp)*zeta_right(k,id) - zdotintfac_right(k,id)*phix_1d(k,jpp))*twoHterm(k,id) + &
         rdzeandzdre(k,id)*twoHterm_r
    Rms5r = phi_1d(k,ipp)*gammaintfac(k,id)*(Rms5term_r*dS(k,id) + Rms5term(k,id)*dSdr(k,id))

    twoHterm_z = -2.0_rk*sin(angle_c)/R/SQr2z2(k,id)**1.5*zeta_right(k,id)*phix_1d(k,jpp)
    Rms5term_z = ( rdotintfac_right(k,id)*phix_1d(k,jpp) - CTJ/dt*phi_1d(k,jpp)*reta_right(k,id) )*twoHterm(k,id) + &
         rdzeandzdre(k,id)*twoHterm_z
    Rms5z = phi_1d(k,ipp)*gammaintfac(k,id)*( Rms5term_z*dS(k,id) + Rms5term(k,id)*dSdz(k,id) )

    Rms5gamma = phi_1d(k,ipp)*phi_1d(k,jpp)*Rms5term(k,id)*dS(k,id)

    !Rms6, stretching&dilatation   
    Rms6r = phi_1d(k,ipp)*gammaintfac(k,id)* ( &
         dSdr(k,id)*Rms6term(k,id) + dS(k,id)* ( &
         phix_1d(k,jpp)*ueta(k,id)/SQr2z2(k,id) &
         - reueandzeve(k,id)/SQr2z2(k,id)**2*2.0_rk*reta_right(k,id)*phix_1d(k,jpp) &
         - phi_1d(k,jpp)/rintfac_right(k,id)**2*uintfac_right(k,id) - flux_r(k,id)*phi_1d(k,jpp)* mtwoH )       )

    Rms6z = phi_1d(k,ipp)*gammaintfac(k,id)* ( &
         dSdz(k,id)*Rms6term(k,id) + dS(k,id)* ( phix_1d(k,jpp)*veta(k,id)/SQr2z2(k,id) &
         - reueandzeve(k,id)/SQr2z2(k,id)**2*2.0_rk*zeta_right(k,id)*phix_1d(k,jpp)  ) )

    Rms6u = phi_1d(k,ipp)*gammaintfac(k,id)*dS(k,id)*&
          ( reta_right(k,id)*phix_1d(k,jpp)/SQr2z2(k,id) + &
         phi_1d(k,jpp)/rintfac_right(k,id) )

    Rms6v = phi_1d(k,ipp)*gammaintfac(k,id)*dS(k,id)* zeta_right(k,id)*phix_1d(k,jpp)/SQr2z2(k,id) 

    Rms6gamma = phi_1d(k,ipp)*phi_1d(k,jpp)*dS(k,id)* Rms6term(k,id)

    
    ! intRms_r_S(k) = Rms1r + Rms31r + Rms6r
    ! intRms_z_S(k) = Rms1z + Rms31z + Rms6z
    ! intRms_u_S(k) = Rms31u + Rms6u !
    ! intRms_v_S(k) = Rms31v + Rms6v!
    ! intRms_cp_S(k) = 0.0_rk!Rms2cp
    ! intRms_gamma_S(k) = Rms1gamma + Rms31gamma + Rms6gamma
    
    intRms_r_S(k) = Rms1r + Rms2r + Rms31r + Rms4r + Rms6r !+ Rms5r + Rms32r !
    intRms_z_S(k) = Rms1z + Rms2z + Rms31z + Rms4z + Rms6z !+ Rms5z + Rms32z !
    !only Rms3u&Rms6u
    intRms_u_S(k) = Rms31u + Rms6u !+ Rms32u !
    !only Rms3v&Rms6v    
    intRms_v_S(k) = Rms31v + Rms6v !+ Rms32v !
    !only Rms2cp
    intRms_cp_S(k) = Rms2cp
    intRms_gamma_S(k) = Rms1gamma + Rms2gamma + Rms31gamma + Rms4gamma + Rms6gamma !+ Rms5gamma + Rms32gamma !
    
 end if !solve_cp.eq.1 .and. surf_adsp.eq.1
 
         end do ! k: three gausspoints
         
     sj(LNOPP(i)+Nr,LNOPP(j)+Nr) = sj(LNOPP(i)+Nr,LNOPP(j)+Nr) + gaussian_quadrature_1d(intRsi_r_S)
     sj(LNOPP(i)+Nr,LNOPP(j)+Nz) = sj(LNOPP(i)+Nr,LNOPP(j)+Nz) + gaussian_quadrature_1d(intRsi_z_S)
     sj(LNOPP(i)+Nr,LNOPP(j)+Nu) = sj(LNOPP(i)+Nr,LNOPP(j)+Nu) + gaussian_quadrature_1d(intRsi_u_S)
     sj(LNOPP(i)+Nr,LNOPP(j)+Nv) = sj(LNOPP(i)+Nr,LNOPP(j)+Nv) + gaussian_quadrature_1d(intRsi_v_S)

     sj(LNOPP(i)+Nu,LNOPP(j)+Nr) = sj(LNOPP(i)+Nu,LNOPP(j)+Nr) + gaussian_quadrature_1d(intRu_r_S)!/Ca
     sj(LNOPP(i)+Nu,LNOPP(j)+Nz) = sj(LNOPP(i)+Nu,LNOPP(j)+Nz) + gaussian_quadrature_1d(intRu_z_S)!/Ca
     if(solve_T.eq.1) &
       sj(LNOPP(i)+Nu,LNOPP(j)+NT) = sj(LNOPP(i)+Nu,LNOPP(j)+NT) + gaussian_quadrature_1d(intRu_T_S)!/Ca

     sj(LNOPP(i)+Nv,LNOPP(j)+Nr) = sj(LNOPP(i)+Nv,LNOPP(j)+Nr) + gaussian_quadrature_1d(intRv_r_S)!/Ca
     sj(LNOPP(i)+Nv,LNOPP(j)+Nz) = sj(LNOPP(i)+Nv,LNOPP(j)+Nz) + gaussian_quadrature_1d(intRv_z_S)!/Ca
     if(solve_T.eq.1) &
       sj(LNOPP(i)+Nv,LNOPP(j)+NT) = sj(LNOPP(i)+Nv,LNOPP(j)+NT) + gaussian_quadrature_1d(intRv_T_S)!/Ca

     if( solve_cp.eq.1 .and. surf_adsp.eq.1 ) then
     sj(LNOPP(i)+MDF(globalNM(m,i))-1,LNOPP(j)+Nr) = gaussian_quadrature_1d(intRms_r_S)
     sj(LNOPP(i)+MDF(globalNM(m,i))-1,LNOPP(j)+Nz) = gaussian_quadrature_1d(intRms_z_S)
     sj(LNOPP(i)+MDF(globalNM(m,i))-1,LNOPP(j)+Nu) = gaussian_quadrature_1d(intRms_u_S)
     sj(LNOPP(i)+MDF(globalNM(m,i))-1,LNOPP(j)+Nv) = gaussian_quadrature_1d(intRms_v_S)
     sj(LNOPP(i)+MDF(globalNM(m,i))-1,LNOPP(j)+Ncp) = gaussian_quadrature_1d(intRms_cp_S)
     sj(LNOPP(i)+MDF(globalNM(m,i))-1,LNOPP(j)+MDF(globalNM(m,j))-1) = gaussian_quadrature_1d(intRms_gamma_S)

     
     end if  !solve_cp.eq.1 .and. surf_adsp.eq.1
     
      end if  !for j = 1,4,7 on the free surface
      
      !evaporation cooling 1
      if(solve_T.eq.1) then
         
      do k = 1, Ng, 1    !three gausspoints
intRt_r_S(k) =  intRt_r_S(k) - phi_1d(k,ipp)*( (-dTdsi(k,id)*2.0_rk*reta_right(k,id)* phieta1_1d(k,j) +   &
     Teta_right(k,id)*( rsi_right(k,id)*phieta1_1d(k,j) + reta_right(k,id)*phisi1_1d(k,j) ) &
     ) *rintfac_right(k,id)/Jp_right(k,id) + &

     ( -dTdsi(k,id)*SQr2z2(k,id) + &
     Teta_right(k,id)*( zsi_right(k,id)*zeta_right(k,id) + rsi_right(k,id)*reta_right(k,id) ) ) *&
     (phi1_1d(k,j)/Jp_right(k,id) - rintfac_right(k,id)/Jp_right(k,id)**2 * Jp_r_right(k,id) ) )

intRt_z_S(k) = intRt_z_S(k) - phi_1d(k,ipp)*( (-dTdsi(k,id)*2.0_rk*zeta_right(k,id)* phieta1_1d(k,j) + &
     Teta_right(k,id)*( zeta_right(k,id)*phisi1_1d(k,j) + zsi_right(k,id)*phieta1_1d(k,j) ) &
     ) *rintfac_right(k,id)/Jp_right(k,id) + &

     ( -dTdsi(k,id)*SQr2z2(k,id) + &
     Teta_right(k,id)*( zsi_right(k,id)*zeta_right(k,id) + rsi_right(k,id)*reta_right(k,id) ) ) *&
     (-rintfac_right(k,id))/Jp_right(k,id)**2 * Jp_z_right(k,id) )

intRt_T_S(k) = intRt_T_S(k) - phi_1d(k,ipp)*( -phisi1_1d(k,j)*SQr2z2(k,id) + &
     phieta1_1d(k,j)*( zsi_right(k,id)*zeta_right(k,id) + rsi_right(k,id)*reta_right(k,id) ) ) * &
     rintfac_right(k,id)/Jp_right(k,id)
      end do
      sj(LNOPP(i)+NT,LNOPP(j)+Nr) = sj(LNOPP(i)+NT,LNOPP(j)+Nr) + gaussian_quadrature_1d(intRt_r_S)
      sj(LNOPP(i)+NT,LNOPP(j)+Nz) = sj(LNOPP(i)+NT,LNOPP(j)+Nz) + gaussian_quadrature_1d(intRt_z_S)
      sj(LNOPP(i)+NT,LNOPP(j)+NT ) = sj(LNOPP(i)+NT,LNOPP(j)+NT ) + gaussian_quadrature_1d(intRt_T_S)
      
      end if   !solve_T.eq.1

      !particle accumulation 1
      if(solve_cp.eq.1) then
         
      do k = 1, Ng, 1    !three gausspoints
intRm_r_S(k) =  intRm_r_S(k) + KBCgroup/Pep* phi_1d(k,ipp)*( (-dcpdsi(k,id)*2.0_rk*reta_right(k,id)* phieta1_1d(k,j) +   &
     cpeta_right(k,id)*( rsi_right(k,id)*phieta1_1d(k,j) + reta_right(k,id)*phisi1_1d(k,j) ) &
     ) *rintfac_right(k,id)/Jp_right(k,id) + &

     ( -dcpdsi(k,id)*SQr2z2(k,id) + &
     cpeta_right(k,id)*( zsi_right(k,id)*zeta_right(k,id) + rsi_right(k,id)*reta_right(k,id) ) ) *&
     (phi1_1d(k,j)/Jp_right(k,id) - rintfac_right(k,id)/Jp_right(k,id)**2 * Jp_r_right(k,id) ) )

intRm_z_S(k) = intRm_z_S(k) + KBCgroup/Pep* phi_1d(k,ipp)*( (-dcpdsi(k,id)*2.0_rk*zeta_right(k,id)* phieta1_1d(k,j) + &
     cpeta_right(k,id)*( zeta_right(k,id)*phisi1_1d(k,j) + zsi_right(k,id)*phieta1_1d(k,j) ) &
     ) *rintfac_right(k,id)/Jp_right(k,id) + &

     ( -dcpdsi(k,id)*SQr2z2(k,id) + &
     cpeta_right(k,id)*( zsi_right(k,id)*zeta_right(k,id) + rsi_right(k,id)*reta_right(k,id) ) ) *&
     (-rintfac_right(k,id))/Jp_right(k,id)**2 * Jp_z_right(k,id) )

intRm_cp_S(k) = intRm_cp_S(k) + KBCgroup/Pep* phi_1d(k,ipp)*( -phisi1_1d(k,j)*SQr2z2(k,id) + &
     phieta1_1d(k,j)*( zsi_right(k,id)*zeta_right(k,id) + rsi_right(k,id)*reta_right(k,id) ) ) * &
     rintfac_right(k,id)/Jp_right(k,id)
      end do
      sj(LNOPP(i)+Ncp,LNOPP(j)+Nr) = sj(LNOPP(i)+Ncp,LNOPP(j)+Nr) + gaussian_quadrature_1d(intRm_r_S)
      sj(LNOPP(i)+Ncp,LNOPP(j)+Nz) = sj(LNOPP(i)+Ncp,LNOPP(j)+Nz) + gaussian_quadrature_1d(intRm_z_S)
      sj(LNOPP(i)+Ncp,LNOPP(j)+Ncp ) = sj(LNOPP(i)+Ncp,LNOPP(j)+Ncp ) + gaussian_quadrature_1d(intRm_cp_S)
      if(surf_adsp.eq.1) sj( LNOPP(i) + Ncp, LNOPP(j)+MDF(globalNM(m,j))-1 ) = &
           sj(LNOPP(i)+Ncp,LNOPP(j)+MDF(globalNM(m,j))-1 ) + gaussian_quadrature_1d(intRm_gamma_S)

      end if   !solve_cp.eq.1

   end if   !for i = 1,4,7 on the free surface


  !free surface from vapor, KBC part2 & evaporation cooling part2
  if(no_vapor.eq.0) then
  if( (BCflagN( globalNM(m,i), 3 ).eq.1 .or. BCflagN( globalNM(m,i), 3 ).eq.3) .and. VE(m).eq.1 ) then
     ipp = i/3 + 1  !phi_1d(k,l)

     ! if( BCflagN( globalNM(m,j),3 ).eq.1 .or. BCflagN( globalNM(m,j),3 ).eq.3 ) then
     !    jpp = j/3 + 1  !phi_1d(k,l)
     
     do k = 1, Ng, 1    !three gausspoints

!KBC2
intRsi_r_S(k) =  phi_1d(k,ipp)*( (-dcdsi(k,id)*2.0_rk*reta_right(k,id)* phieta0_1d(k,j) +   &
     dcdeta(k,id)*( rsi_right(k,id)*phieta0_1d(k,j) + reta_right(k,id)*phisi0_1d(k,j) ) &
     ) *rintfac_right(k,id)/Jp_right(k,id) + &

     ( -dcdsi(k,id)*( reta_right(k,id)**2 + zeta_right(k,id)**2 ) + &
     dcdeta(k,id)*( zsi_right(k,id)*zeta_right(k,id) + rsi_right(k,id)*reta_right(k,id) ) ) *&
     (phi0_1d(k,j)/Jp_right(k,id) - rintfac_right(k,id)/Jp_right(k,id)**2 * Jp_r_right(k,id) ) )

intRsi_z_S(k) = phi_1d(k,ipp)*( (-dcdsi(k,id)*2.0_rk*zeta_right(k,id)* phieta0_1d(k,j) + &
     dcdeta(k,id)*( zeta_right(k,id)*phisi0_1d(k,j) + zsi_right(k,id)*phieta0_1d(k,j) ) &
     ) *rintfac_right(k,id)/Jp_right(k,id) + &

     ( -dcdsi(k,id)*( reta_right(k,id)**2 + zeta_right(k,id)**2 ) + &
     dcdeta(k,id)*( zsi_right(k,id)*zeta_right(k,id) + rsi_right(k,id)*reta_right(k,id) ) ) *&
     (-rintfac_right(k,id))/Jp_right(k,id)**2 * Jp_z_right(k,id) )

intRsi_c_S(k) = phi_1d(k,ipp)*( -phisi0_1d(k,j)*( reta_right(k,id)**2 + zeta_right(k,id)**2 ) + &
     phieta0_1d(k,j)*( zsi_right(k,id)*zeta_right(k,id) + rsi_right(k,id)*reta_right(k,id) ) ) * &
     rintfac_right(k,id)/Jp_right(k,id)


!evaporation cooling 2
if(solve_T.eq.1) then
   intRt_r_S(k) = REH* intRsi_r_S(k)
   intRt_z_S(k) = REH* intRsi_z_S(k)
   intRt_c_S(k) = REH* intRsi_c_S(k)
end if

!particle accumulation 2  !??
if(solve_cp.eq.1) then
   intRm_r_S(k) = cpintfac_right(k,id) * intRsi_r_S(k)
   intRm_z_S(k) = cpintfac_right(k,id) * intRsi_z_S(k)
   intRm_c_S(k) = cpintfac_right(k,id) * intRsi_c_S(k)
   
   if( BCflagN( globalNM(m,j),3 ).eq.1 .or. BCflagN( globalNM(m,j),3 ).eq.3 ) then
      jpp = j/3 + 1  !phi_1d(k,l)
      intRm_cp_S(k) = phi_1d(k,ipp)*( &
           -dcdsi(k,id)*( reta_right(k,id)**2 + zeta_right(k,id)**2 ) + &
           dcdeta(k,id)*( zsi_right(k,id)*zeta_right(k,id) + rsi_right(k,id)*reta_right(k,id) ) &
           ) *phi_1d(k,jpp) *rintfac_right(k,id)/Jp_right(k,id)
   end if
   
end if  !solve_cp.eq.1


     end do
     
     sj(LNOPP(i)+Nr,LNOPP(j)+Nr) = sj(LNOPP(i)+Nr,LNOPP(j)+Nr) + gaussian_quadrature_1d(intRsi_r_S)
     sj(LNOPP(i)+Nr,LNOPP(j)+Nz) = sj(LNOPP(i)+Nr,LNOPP(j)+Nz) + gaussian_quadrature_1d(intRsi_z_S)
     sj(LNOPP(i)+Nr,LNOPP(j) + MDF( globalNM(m,j) ) -1 ) = &
          sj(LNOPP(i)+Nr,LNOPP(j) + MDF( globalNM(m,j) ) -1 ) + gaussian_quadrature_1d(intRsi_c_S)
     
     if(solve_T.eq.1) then
     sj(LNOPP(i)+NT,LNOPP(j)+Nr) = sj(LNOPP(i)+NT,LNOPP(j)+Nr) + gaussian_quadrature_1d(intRt_r_S)
     sj(LNOPP(i)+NT,LNOPP(j)+Nz) = sj(LNOPP(i)+NT,LNOPP(j)+Nz) + gaussian_quadrature_1d(intRt_z_S)
     sj(LNOPP(i)+NT,LNOPP(j)+ MDF( globalNM(m,j) ) -1 ) = &
          sj(LNOPP(i)+NT,LNOPP(j)+ MDF( globalNM(m,j) ) -1 ) + gaussian_quadrature_1d(intRt_c_S)
     end if

     if(solve_cp.eq.1) then
     sj(LNOPP(i)+Ncp,LNOPP(j)+Nr) = sj(LNOPP(i)+Ncp,LNOPP(j)+Nr) + gaussian_quadrature_1d(intRm_r_S)
     sj(LNOPP(i)+Ncp,LNOPP(j)+Nz) = sj(LNOPP(i)+Ncp,LNOPP(j)+Nz) + gaussian_quadrature_1d(intRm_z_S)
     sj(LNOPP(i)+Ncp,LNOPP(j)+ MDF( globalNM(m,j) ) -1 ) = &
          sj(LNOPP(i)+Ncp,LNOPP(j)+ MDF( globalNM(m,j) ) -1 ) + gaussian_quadrature_1d(intRm_c_S)
     
     if( BCflagN( globalNM(m,j),3 ).eq.1 .or. BCflagN( globalNM(m,j),3 ).eq.3 ) &
          sj(LNOPP(i)+Ncp,LNOPP(j)+Ncp ) = &
          sj(LNOPP(i)+Ncp,LNOPP(j)+Ncp ) + gaussian_quadrature_1d(intRm_cp_S)
     
     end if !solve_cp.eq.1

  end if   !for i = 1,4,7 on the free surface
  end if   !solve for vapor phase
  
end if   !for s_mode=0


! if(initial_vapor_solved.eq.1 .and. m.eq.28 .and. i.eq.4 .and. j.eq.2) then
!    write(*,*) sj(LNOPP(i)+NT,LNOPP(j)+Nr)
!    pause
! end if

  return
end subroutine SI_in_sj
