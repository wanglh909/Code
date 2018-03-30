subroutine SI_in_sj(m,i,j, sj, LNVar, LNOPP, id)                     !adding SI to sj
  use kind
  use data
  use Ldata
  use NOP_mod, only: gaussian_quadrature_1d

  implicit none

  integer(kind=ik), intent(in):: m,i,j, LNVar, LNOPP(9), id
  real(kind=rk), intent(out):: sj(LNVar, LNVar)
  real(kind=rk):: flux

  integer(kind=ik):: k, ipp, jpp
  real(kind=rk):: intRsi_r_S(Ng), intRsi_z_S(Ng), intReta_r_S(Ng), intReta_z_S(Ng)
  real(kind=rk):: intRsi_u_S(Ng), intRsi_v_S(Ng), intRsi_c_S(Ng), &
       intRu_r_S(Ng), intRu_z_S(Ng), intRv_r_S(Ng), intRv_z_S(Ng), &
       intRu_T_S(Ng), intRv_T_S(Ng)

  real(kind=rk):: intRt_r_S(Ng), intRt_z_S(Ng), intRt_c_S(Ng), intRt_T_S(Ng)
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

intReta_r_S(k) = phix_1d(k,ipp) * geta_size(m) / ( reta_right(k,id)**2 + zeta_right(k,id)**2 ) &
                   * 2.0_rk*reta_right(k,id) * phix_1d(k,jpp)
intReta_z_S(k) = phix_1d(k,ipp) * geta_size(m) / ( reta_right(k,id)**2 + zeta_right(k,id)**2 ) &
                   * 2.0_rk*zeta_right(k,id) * phix_1d(k,jpp)

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

if(uniflux) then  !uniflux( flux( angle_c,rintfac_right(k,id) ) )
   
   !KBC2
   intRsi_r_S(k) =  phi_1d(k,ipp)* flux( angle_c,rintfac_right(k,id) ) *&
        ( reta_right(k,id)**2 + zeta_right(k,id)**2 )**(-0.5_rk) *&
        reta_right(k,id) *phix_1d(k,jpp)* rintfac_right(k,id) + &
        
        phi_1d(k,ipp)* flux( angle_c,rintfac_right(k,id) ) *&
        ( reta_right(k,id)**2 + zeta_right(k,id)**2 )**0.5_rk *phi_1d(k,jpp)

   intRsi_z_S(k) = phi_1d(k,ipp)* flux( angle_c,rintfac_right(k,id) ) *&
        ( reta_right(k,id)**2 + zeta_right(k,id)**2 )**(-0.5_rk) *&
        zeta_right(k,id) *phix_1d(k,jpp)* rintfac_right(k,id)


   !evaporation cooling 2
   intRt_r_S(k) = REH* intRsi_r_S(k)
   intRt_z_S(k) = REH* intRsi_z_S(k)
   if(true_uniflux.eq.1) then
      intRsi_r_S(k) = intRsi_r_S(k) / flux( angle_c,rintfac_right(k,id) )
      intRsi_z_S(k) = intRsi_z_S(k) / flux( angle_c,rintfac_right(k,id) )
   end if
end if


!KBC1
intRsi_r_S(k) = intRsi_r_S(k) + KBCgroup* ( phi_1d(k,ipp)*( zeta_right(k,id)*CTJ/dt*phi_1d(k,jpp) + &
     phix_1d(k,jpp)*( vintfac_right(k,id) - zdotintfac_right(k,id) ) )*rintfac_right(k,id) + &
     
     phi_1d(k,ipp)*( -zeta_right(k,id)*( uintfac_right(k,id) - rdotintfac_right(k,id) ) + &
     reta_right(k,id)*( vintfac_right(k,id) - zdotintfac_right(k,id) ) )*phi_1d(k,jpp) )

intRsi_z_S(k) = intRsi_z_S(k) + KBCgroup* ( phi_1d(k,ipp)*( -phix_1d(k,jpp)*&
     ( uintfac_right(k,id) - rdotintfac_right(k,id) ) - &
     reta_right(k,id)*CTJ/dt*phi_1d(k,jpp) )*rintfac_right(k,id) )

intRsi_u_S(k) = intRsi_u_S(k) + KBCgroup* ( -phi_1d(k,jpp)*phi_1d(k,ipp)*zeta_right(k,id)*rintfac_right(k,id) )

intRsi_v_S(k) = intRsi_v_S(k) + KBCgroup* ( phi_1d(k,jpp)*phi_1d(k,ipp)*reta_right(k,id)*rintfac_right(k,id) )


!traction BC
intRu_r_S(k) = ( phix_1d(k,ipp)*phix_1d(k,jpp) / ( reta_right(k,id)**2 + zeta_right(k,id)**2 ) - &
     reta_right(k,id)*phix_1d(k,ipp) / ( reta_right(k,id)**2 + zeta_right(k,id)**2 )**2 * 2.0_rk*reta_right(k,id)*phix_1d(k,jpp) - &
     1.0_rk/rintfac_right(k,id)**2 *phi_1d(k,jpp)*phi_1d(k,ipp) ) * &
     rintfac_right(k,id)*( reta_right(k,id)**2 + zeta_right(k,id)**2 )**0.5_rk + &
     
     ( reta_right(k,id)*phix_1d(k,ipp) / ( reta_right(k,id)**2 + zeta_right(k,id)**2 ) + phi_1d(k,ipp)/rintfac_right(k,id) )&
     *( phi_1d(k,jpp)*( reta_right(k,id)**2 + zeta_right(k,id)**2 )**0.5_rk + &
     rintfac_right(k,id) / ( reta_right(k,id)**2 + zeta_right(k,id)**2 )**0.5_rk *reta_right(k,id)*phix_1d(k,jpp) )

 intRu_r_S(k) = intRu_r_S(k) - phi_1d(k,ipp) *beta *Teta_right(k,id)* ( &
      -( reta_right(k,id)**2 + zeta_right(k,id)**2 )**(-1.5_rk) *reta_right(k,id)**2 *phix_1d(k,jpp) *rintfac_right(k,id) + &
      ( reta_right(k,id)**2 + zeta_right(k,id)**2 )**(-0.5_rk) *&
      ( phix_1d(k,jpp) *rintfac_right(k,id) + reta_right(k,id) *phi_1d(k,jpp) ) )

intRu_z_S(k) = -reta_right(k,id)*phix_1d(k,ipp)*( reta_right(k,id)**2 + zeta_right(k,id)**2 )**(-1.5_rk) *&
     2.0_rk*zeta_right(k,id)*phix_1d(k,jpp)*rintfac_right(k,id) + &
     
     ( reta_right(k,id)*phix_1d(k,ipp) / ( reta_right(k,id)**2 + zeta_right(k,id)**2 ) + phi_1d(k,ipp)/rintfac_right(k,id) ) *&
     rintfac_right(k,id) / ( reta_right(k,id)**2 + zeta_right(k,id)**2 )**0.5_rk *zeta_right(k,id)*phix_1d(k,jpp)

 intRu_z_S(k) = intRu_z_S(k) + phi_1d(k,ipp) *beta *Teta_right(k,id)* &
      ( reta_right(k,id)**2 + zeta_right(k,id)**2 )**(-1.5_rk)* &
      zeta_right(k,id) *phix_1d(k,jpp) *reta_right(k,id) *rintfac_right(k,id) 

intRu_T_S(k) = -phi_1d(k,ipp) *beta *phix_1d(k,jpp) * &
      ( reta_right(k,id)**2 + zeta_right(k,id)**2 )**(-0.5_rk)* reta_right(k,id) *rintfac_right(k,id) 


intRv_r_S(k) = zeta_right(k,id)*phix_1d(k,ipp)*&
     ( -( reta_right(k,id)**2 + zeta_right(k,id)**2 )**(-1.5_rk) *reta_right(k,id)*phix_1d(k,jpp)*rintfac_right(k,id) + &
     ( reta_right(k,id)**2 + zeta_right(k,id)**2 )**(-0.5_rk) *phi_1d(k,jpp) )

 intRv_r_S(k) = intRv_r_S(k) - phi_1d(k,ipp)*beta*Teta_right(k,id)*zeta_right(k,id)* ( &
      -( reta_right(k,id)**2 + zeta_right(k,id)**2 )**(-1.5_rk) *reta_right(k,id) *phix_1d(k,jpp) *rintfac_right(k,id) + &
      ( reta_right(k,id)**2 + zeta_right(k,id)**2 )**(-0.5_rk) * phi_1d(k,jpp) )


intRv_z_S(k) = phix_1d(k,ipp)*phix_1d(k,jpp) / sqrt( reta_right(k,id)**2 + zeta_right(k,id)**2 ) *rintfac_right(k,id) - &
     zeta_right(k,id)**2 *phix_1d(k,ipp)*phix_1d(k,jpp) / ( reta_right(k,id)**2 + zeta_right(k,id)**2 )**1.5_rk *rintfac_right(k,id)

 intRv_z_S(k) = intRv_z_S(k) - phi_1d(k,ipp) *beta *Teta_right(k,id)* ( &
      phix_1d(k,jpp) *( reta_right(k,id)**2 + zeta_right(k,id)**2 )**(-0.5_rk) - &
      ( reta_right(k,id)**2 + zeta_right(k,id)**2 )**(-1.5_rk)*zeta_right(k,id)**2 *phix_1d(k,jpp) ) &
      *rintfac_right(k,id) 

intRv_T_S(k) = -phi_1d(k,ipp) *beta *phix_1d(k,jpp) * &
      ( reta_right(k,id)**2 + zeta_right(k,id)**2 )**(-0.5_rk)* zeta_right(k,id) *rintfac_right(k,id) 

         end do
     
     sj(LNOPP(i)+Nr,LNOPP(j)+Nr) = sj(LNOPP(i)+Nr,LNOPP(j)+Nr) + gaussian_quadrature_1d(intRsi_r_S)
     sj(LNOPP(i)+Nr,LNOPP(j)+Nz) = sj(LNOPP(i)+Nr,LNOPP(j)+Nz) + gaussian_quadrature_1d(intRsi_z_S)
     sj(LNOPP(i)+Nr,LNOPP(j)+Nu) = sj(LNOPP(i)+Nr,LNOPP(j)+Nu) + gaussian_quadrature_1d(intRsi_u_S)
     sj(LNOPP(i)+Nr,LNOPP(j)+Nv) = sj(LNOPP(i)+Nr,LNOPP(j)+Nv) + gaussian_quadrature_1d(intRsi_v_S)

     sj(LNOPP(i)+Nu,LNOPP(j)+Nr) = sj(LNOPP(i)+Nu,LNOPP(j)+Nr) + gaussian_quadrature_1d(intRu_r_S)!/Ca
     sj(LNOPP(i)+Nu,LNOPP(j)+Nz) = sj(LNOPP(i)+Nu,LNOPP(j)+Nz) + gaussian_quadrature_1d(intRu_z_S)!/Ca
     sj(LNOPP(i)+Nu,LNOPP(j)+NT) = sj(LNOPP(i)+Nu,LNOPP(j)+NT) + gaussian_quadrature_1d(intRu_T_S)!/Ca

     sj(LNOPP(i)+Nv,LNOPP(j)+Nr) = sj(LNOPP(i)+Nv,LNOPP(j)+Nr) + gaussian_quadrature_1d(intRv_r_S)!/Ca
     sj(LNOPP(i)+Nv,LNOPP(j)+Nz) = sj(LNOPP(i)+Nv,LNOPP(j)+Nz) + gaussian_quadrature_1d(intRv_z_S)!/Ca
     sj(LNOPP(i)+Nv,LNOPP(j)+NT) = sj(LNOPP(i)+Nv,LNOPP(j)+NT) + gaussian_quadrature_1d(intRv_T_S)!/Ca

      end if  !for j = 1,4,7 on the free surface
      
      !evaporation cooling 1
      do k = 1, Ng, 1    !three gausspoints
intRt_r_S(k) =  intRt_r_S(k) - phi_1d(k,ipp)*( (-dTdsi(k,id)*2.0_rk*reta_right(k,id)* phieta1_1d(k,j) +   &
     Teta_right(k,id)*( rsi_right(k,id)*phieta1_1d(k,j) + reta_right(k,id)*phisi1_1d(k,j) ) &
     ) *rintfac_right(k,id)/Jp_right(k,id) + &

     ( -dTdsi(k,id)*( reta_right(k,id)**2 + zeta_right(k,id)**2 ) + &
     Teta_right(k,id)*( zsi_right(k,id)*zeta_right(k,id) + rsi_right(k,id)*reta_right(k,id) ) ) *&
     (phi1_1d(k,j)/Jp_right(k,id) - rintfac_right(k,id)/Jp_right(k,id)**2 * Jp_r_right(k,id) ) )

intRt_z_S(k) = intRt_z_S(k) - phi_1d(k,ipp)*( (-dTdsi(k,id)*2.0_rk*zeta_right(k,id)* phieta1_1d(k,j) + &
     Teta_right(k,id)*( zeta_right(k,id)*phisi1_1d(k,j) + zsi_right(k,id)*phieta1_1d(k,j) ) &
     ) *rintfac_right(k,id)/Jp_right(k,id) + &

     ( -dTdsi(k,id)*( reta_right(k,id)**2 + zeta_right(k,id)**2 ) + &
     Teta_right(k,id)*( zsi_right(k,id)*zeta_right(k,id) + rsi_right(k,id)*reta_right(k,id) ) ) *&
     (-rintfac_right(k,id))/Jp_right(k,id)**2 * Jp_z_right(k,id) )

intRt_T_S(k) = intRt_T_S(k) - phi_1d(k,ipp)*( -phisi1_1d(k,j)*( reta_right(k,id)**2 + zeta_right(k,id)**2 ) + &
     phieta1_1d(k,j)*( zsi_right(k,id)*zeta_right(k,id) + rsi_right(k,id)*reta_right(k,id) ) ) * &
     rintfac_right(k,id)/Jp_right(k,id)
      end do
      sj(LNOPP(i)+NT,LNOPP(j)+Nr) = sj(LNOPP(i)+NT,LNOPP(j)+Nr) + gaussian_quadrature_1d(intRt_r_S)
      sj(LNOPP(i)+NT,LNOPP(j)+Nz) = sj(LNOPP(i)+NT,LNOPP(j)+Nz) + gaussian_quadrature_1d(intRt_z_S)
      sj(LNOPP(i)+NT,LNOPP(j)+NT ) = sj(LNOPP(i)+NT,LNOPP(j)+NT ) + gaussian_quadrature_1d(intRt_T_S)
      

   end if   !for i = 1,4,7 on the free surface


  !free surface from vapor, KBC part2 & evaporation cooling part2
  if(.not.uniflux) then
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
intRt_r_S(k) = REH* intRsi_r_S(k)
intRt_z_S(k) = REH* intRsi_z_S(k)
intRt_c_S(k) = REH* intRsi_c_S(k)


     end do
     
     sj(LNOPP(i)+Nr,LNOPP(j)+Nr) = sj(LNOPP(i)+Nr,LNOPP(j)+Nr) + gaussian_quadrature_1d(intRsi_r_S)
     sj(LNOPP(i)+Nr,LNOPP(j)+Nz) = sj(LNOPP(i)+Nr,LNOPP(j)+Nz) + gaussian_quadrature_1d(intRsi_z_S)
     sj(LNOPP(i)+Nr,LNOPP(j) + MDF( globalNM(m,j) ) -1 ) = &
          sj(LNOPP(i)+Nr,LNOPP(j) + MDF( globalNM(m,j) ) -1 ) + gaussian_quadrature_1d(intRsi_c_S)

     sj(LNOPP(i)+NT,LNOPP(j)+Nr) = sj(LNOPP(i)+NT,LNOPP(j)+Nr) + gaussian_quadrature_1d(intRt_r_S)
     sj(LNOPP(i)+NT,LNOPP(j)+Nz) = sj(LNOPP(i)+NT,LNOPP(j)+Nz) + gaussian_quadrature_1d(intRt_z_S)
     sj(LNOPP(i)+NT,LNOPP(j)+ MDF( globalNM(m,j) ) -1 ) = &
          sj(LNOPP(i)+NT,LNOPP(j)+ MDF( globalNM(m,j) ) -1 ) + gaussian_quadrature_1d(intRt_c_S)
     

  end if   !for i = 1,4,7 on the free surface
  end if   !for uniflux = 0
  
end if   !for s_mode=0


! if(initial_vapor_solved.eq.1 .and. m.eq.28 .and. i.eq.4 .and. j.eq.2) then
!    write(*,*) sj(LNOPP(i)+NT,LNOPP(j)+Nr)
!    pause
! end if

  return
end subroutine SI_in_sj
