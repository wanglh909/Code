subroutine values_in_sj(m,i,j,id)
  use kind
  use data
  use Ldata

  implicit none
  integer(kind=ik), intent(in):: m,i,j, id
  integer(kind=ik):: k, l, n

  do k = 1, Ng, 1          !relate to a(Ng) (the value of si for gaussian_quadrature)
     do l = 1, Ng, 1           !relate to a(Ng) (the value of eta for gaussian_quadrature)

        !Actually, these terms are (k,l,i,j), but can be refurbished in each cycle of i&j, so just need (k,l)
        !s_orth_r(k,l,i,j)
        s_orth_r(k,l,id) = 1.0_rk / ( reta(k,l,id)**2 + zeta(k,l,id)**2 ) * ( rsi(k,l,id)*phisi(k,l,j)/s_orth(k,l,id) &
             - reta(k,l,id)*phieta(k,l,j)*s_orth(k,l,id) )

        !s_orth_z(k,l,i,j)
        s_orth_z(k,l,id) = 1.0_rk / ( reta(k,l,id)**2 + zeta(k,l,id)**2 ) * ( zsi(k,l,id)*phisi(k,l,j)/s_orth(k,l,id) &
             - zeta(k,l,id)*phieta(k,l,j)*s_orth(k,l,id) )

        !Aterm_r(k,l,i,j)
        Aterm_r(k,l,id) = 2.0_rk*reta(k,l,id)*phieta(k,l,j)*phisi(k,l,i) &
             - ( rsi(k,l,id)*phieta(k,l,j) + reta(k,l,id)*phisi(k,l,j) )*phieta(k,l,i)

        !Aterm_z(k,l,i,j)
        Aterm_z(k,l,id) = 2.0_rk*zeta(k,l,id)*phieta(k,l,j)*phisi(k,l,i) &
             - ( zsi(k,l,id)*phieta(k,l,j) + zeta(k,l,id)*phisi(k,l,j) )*phieta(k,l,i)

        !Bterm_r(k,l,i,j)
        Bterm_r(k,l,id) = - ( reta(k,l,id)*phisi(k,l,j) + rsi(k,l,id)*phieta(k,l,j) )*phisi(k,l,i) &
             + 2.0_rk*rsi(k,l,id)*phisi(k,l,j)*phieta(k,l,i) 

        !Bterm_z(k,l,i,j)
        Bterm_z(k,l,id) = - ( zeta(k,l,id)*phisi(k,l,j) + zsi(k,l,id)*phieta(k,l,j) )*phisi(k,l,i) &
             + 2.0_rk*zsi(k,l,id)*phisi(k,l,j)*phieta(k,l,i) 

        !Jp_r((k,l,i,j)
        Jp_r(k,l,id) = ( phisi(k,l,j)*zeta(k,l,id) - phieta(k,l,j)*zsi(k,l,id) ) 

        !Jp_z((k,l,i,j)
        Jp_z(k,l,id) = ( rsi(k,l,id)*phieta(k,l,j) - reta(k,l,id)*phisi(k,l,j) )

        !r|Jp|_r(k,l,i,j)
        rJp_r(k,l,id) = ( phi(k,l,j)*(Jp(k,l,id)) + rintfac(k,l,id)*Jp_r(k,l,id) )*Jpsign(k,l,id)


        if(s_mode.eq.0) then
           
           do n = 1, 9            !phir_r(k,l,n,j)
phir_r(k,l,n,id) = -1.0_rk/Jp(k,l,id)**2 *Jp_r(k,l,id)*( phisi(k,l,n)*zeta(k,l,id) - phieta(k,l,n)*zsi(k,l,id) )     !d phir(n) / d rj

phir_z(k,l,n,id) = -1.0_rk/Jp(k,l,id)**2 *Jp_z(k,l,id)*( phisi(k,l,n)*zeta(k,l,id) - phieta(k,l,n)*zsi(k,l,id) ) + &
                   1.0_rk/Jp(k,l,id)*( phisi(k,l,n)*phieta(k,l,j) - phieta(k,l,n)*phisi(k,l,j) )

phiz_r(k,l,n,id) = -1.0_rk/Jp(k,l,id)**2 *Jp_r(k,l,id)*( -phisi(k,l,n)*reta(k,l,id) + phieta(k,l,n)*rsi(k,l,id) ) + &
                   1.0_rk/Jp(k,l,id)*( -phisi(k,l,n)*phieta(k,l,j) + phieta(k,l,n)*phisi(k,l,j) )

phiz_z(k,l,n,id) = -1.0_rk/Jp(k,l,id)**2 *Jp_z(k,l,id)*( -phisi(k,l,n)*reta(k,l,id) + phieta(k,l,n)*rsi(k,l,id) )
           end do

           if(VE(m).eq.0) then

              !urintfac_r(k,l,j)
              urintfac_r(k,l,id) = 0.0_rk
              urintfac_z(k,l,id) = 0.0_rk
              uzintfac_r(k,l,id) = 0.0_rk
              uzintfac_z(k,l,id) = 0.0_rk
              vrintfac_r(k,l,id) = 0.0_rk
              vrintfac_z(k,l,id) = 0.0_rk
              vzintfac_r(k,l,id) = 0.0_rk
              vzintfac_z(k,l,id) = 0.0_rk
              Trintfac_r(k,l,id) = 0.0_rk
              Trintfac_z(k,l,id) = 0.0_rk
              Tzintfac_r(k,l,id) = 0.0_rk
              Tzintfac_z(k,l,id) = 0.0_rk
              do n = 1, 9
                 urintfac_r(k,l,id) = urintfac_r(k,l,id) + ulocal(n,id)*phir_r(k,l,n,id)
                 urintfac_z(k,l,id) = urintfac_z(k,l,id) + ulocal(n,id)*phir_z(k,l,n,id)
                 uzintfac_r(k,l,id) = uzintfac_r(k,l,id) + ulocal(n,id)*phiz_r(k,l,n,id)
                 uzintfac_z(k,l,id) = uzintfac_z(k,l,id) + ulocal(n,id)*phiz_z(k,l,n,id)
                 vrintfac_r(k,l,id) = vrintfac_r(k,l,id) + vlocal(n,id)*phir_r(k,l,n,id)
                 vrintfac_z(k,l,id) = vrintfac_z(k,l,id) + vlocal(n,id)*phir_z(k,l,n,id)
                 vzintfac_r(k,l,id) = vzintfac_r(k,l,id) + vlocal(n,id)*phiz_r(k,l,n,id)
                 vzintfac_z(k,l,id) = vzintfac_z(k,l,id) + vlocal(n,id)*phiz_z(k,l,n,id)
                 Trintfac_r(k,l,id) = Trintfac_r(k,l,id) + Tlocal(n,id)*phir_r(k,l,n,id)
                 Trintfac_z(k,l,id) = Trintfac_z(k,l,id) + Tlocal(n,id)*phir_z(k,l,n,id)
                 Tzintfac_r(k,l,id) = Tzintfac_r(k,l,id) + Tlocal(n,id)*phiz_r(k,l,n,id)
                 Tzintfac_z(k,l,id) = Tzintfac_z(k,l,id) + Tlocal(n,id)*phiz_z(k,l,n,id)
              end do

           else if(VE(m).eq.1) then
              crintfac_r(k,l,id) = 0.0_rk
              crintfac_z(k,l,id) = 0.0_rk
              czintfac_r(k,l,id) = 0.0_rk
              czintfac_z(k,l,id) = 0.0_rk
              do n = 1, 9
                 crintfac_r(k,l,id) = crintfac_r(k,l,id) + clocal(n,id)*phir_r(k,l,n,id)
                 crintfac_z(k,l,id) = crintfac_z(k,l,id) + clocal(n,id)*phir_z(k,l,n,id)
                 czintfac_r(k,l,id) = czintfac_r(k,l,id) + clocal(n,id)*phiz_r(k,l,n,id)
                 czintfac_z(k,l,id) = czintfac_z(k,l,id) + clocal(n,id)*phiz_z(k,l,n,id)
              end do

           else  !VE = 5
              Trintfac_r(k,l,id) = 0.0_rk
              Trintfac_z(k,l,id) = 0.0_rk
              Tzintfac_r(k,l,id) = 0.0_rk
              Tzintfac_z(k,l,id) = 0.0_rk
              do n = 1, 9
                 Trintfac_r(k,l,id) = Trintfac_r(k,l,id) + Tlocal(n,id)*phir_r(k,l,n,id)
                 Trintfac_z(k,l,id) = Trintfac_z(k,l,id) + Tlocal(n,id)*phir_z(k,l,n,id)
                 Tzintfac_r(k,l,id) = Tzintfac_r(k,l,id) + Tlocal(n,id)*phiz_r(k,l,n,id)
                 Tzintfac_z(k,l,id) = Tzintfac_z(k,l,id) + Tlocal(n,id)*phiz_z(k,l,n,id)
              end do

           end if  !for VE

        end if  !for s_mode=0


     end do
  end do            !end loop for k&l

!----------------------------------------for sj_SI-------------------------------------

  if( BCflagN( globalNM(m,i),3 ).eq.1 .and. VE(m).eq.0) then
     do k = 1, Ng
        !Jp_r_right(k,j)
        Jp_r_right(k,id) = ( phisi1_1d(k,j)*zeta_right(k,id) - phieta1_1d(k,j)*zsi_right(k,id) ) 

        !Jp_z_right(k,j)
        Jp_z_right(k,id) = ( rsi_right(k,id)*phieta1_1d(k,j) - reta_right(k,id)*phisi1_1d(k,j) )
     end do
  end if
  if( BCflagN( globalNM(m,i),3 ).eq.1 .and. VE(m).eq.1 ) then
     do k = 1, Ng
        !Jp_r_right(k,j)
        Jp_r_right(k,id) = ( phisi0_1d(k,j)*zeta_right(k,id) - phieta0_1d(k,j)*zsi_right(k,id) ) 

        !Jp_z_right(k,j)
        Jp_z_right(k,id) = ( rsi_right(k,id)*phieta0_1d(k,j) - reta_right(k,id)*phisi0_1d(k,j) )
     end do
  end if


  return
end subroutine values_in_sj
