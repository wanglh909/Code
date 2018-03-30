subroutine jac_check_0
  use kind
  use data, only: Jac, NOPP, NTN, Res, NVar, MDF, folder, globalNM
  use NOP_mod, only: NOPP_define
  implicit none
  
  integer(kind=ik):: i, j, k
  real(kind=rk):: sum

  ! call NOPP_define
  open(unit = 10, file = trim(folder)//'JacSum.dat', status = 'replace')

  do i = 1, NTN
     do k = 1, MDF(i)

  ! do i = 1, NVar
        sum = 0.0_rk
        do j = 1, NVar
           sum = sum + Jac(NOPP(i)+k-1,j)**2
        end do
        ! if(sum.eq.0.0_rk) then
        ! do j = 1, NTN
        !    if( i.lt.NOPP(j) ) exit
        ! end do
        ! j = j-1

        write(10,'(A,i4,A,es14.7,A,i4,A,i4)') 'Jac, i =', NOPP(i)+k-1, ', Jac row =', sum, ', node=', i, ', variable:', k
        write(10,'(A,es14.7)') 'Res', Res(NOPP(i)+k-1)

        if(sum.eq.0.0_rk .or. sum.ne.sum) then
           write(*,*) 'Jac, i =', NOPP(i)+k-1, ', Jac row =', sum, ', node=', i, ', variable:', k
           write(*,*) 'Res', Res(NOPP(i)+k-1)
        end if


        ! end if
  ! end do

     end do
  end do
  close(10)
  
  open(unit = 10, file = trim(folder)//'Jac.dat', status = 'replace')
  do i = 1, NVar  !NOPP(globalNM(61,1))+5, NOPP(globalNM(61,1))+5 !
     write(10, '(i8)',ADVANCE='NO') i
     do j = 1, NTN
        do k = 1, MDF(j)
           write(10, '(es14.7)',ADVANCE='NO') Jac(NOPP(j)+k-1,i)
           ! if(Jac(NOPP(j)+k-1,i).ne.0.0_rk) write(*,*) j,k, Jac(NOPP(j)+k-1,i)
        end do
     end do
     write(10, '(es14.7)',ADVANCE='yes')
  end do
  close(10)
  pause


  return
end subroutine jac_check_0
