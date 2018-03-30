subroutine reverse_sj(m,i,j, sj, LNVar, LNOPP, id)
  use kind
  use data, only: WFLAG, globalNM, Nr, Nz, Np

  implicit none

  integer(kind=ik), intent(in):: m,i,j, LNVar, LNOPP(9), id
  real(kind=rk), intent(out):: sj(LNVar, LNVar)

  integer(kind=ik):: k
  real(kind=rk):: conv(5)

  if( ( WFLAG(m).eq.1 .and. ( (i.eq.1) .or. (i.eq.2) .or. (i.eq.3) ) ) .or. &
       ( WFLAG(m).eq.2 .and. ( (i.eq.1) .or. (i.eq.2) ) ) ) then

     do k = Nr, Nz
        conv(k+1) = sj(LNOPP(i)+Nr,LNOPP(j)+k)  !k = Nr, k+1 = NNr
        sj(LNOPP(i)+Nr,LNOPP(j)+k) = -sj(LNOPP(i)+Nz,LNOPP(j)+k)
        sj(LNOPP(i)+Nz,LNOPP(j)+k) = conv(k+1)
     end do
     ! conv = sj(i,j,NNr,:)
     ! sj(i,j,NNr,:) = -sj(i,j,NNz,:)
     ! sj(i,j,NNz,:) = conv

  else if( WFLAG(m).eq.2 .and. i.eq.3 ) then

     do k = Nr, Nz
        sj(LNOPP(i)+Nr,LNOPP(j)+k) = sj(LNOPP(i)+Nr,LNOPP(j)+k) - sj(LNOPP(i)+Nz,LNOPP(j)+k)
        sj(LNOPP(i)+Nz,LNOPP(j)+k) = 0.0_rk
     end do
     ! sj(i,j,NNr,:) = sj(i,j,NNr,:) - sj(i,j,NNz,:)
     ! sj(i,j,NNz,:) = 0.0_rk

  end if
  return
end subroutine reverse_sj
