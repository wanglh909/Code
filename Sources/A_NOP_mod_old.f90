
module NOP_mod
  use kind
  use data

  implicit none

contains



  subroutine NOP
    implicit none
    integer(kind=ik):: eleN, locN, eleNp
    integer(kind=ik):: i,j,k,l,n
    integer(kind=ik):: LNrowN, LNcolumnN,  NrowNp,NcolumnNp, NrowN, NcolumnN
    !LNrowN: local rowN of the node, LNcolumnN: local columnN of the node
    !NrowN: rowN of the node, NcolumnN: columnN of the node
    !( NrowNp: NrowN just for Region 1 )
    integer(kind=ik):: ErowN  !ErowN: rowN of the element ?can be deleted

    integer(kind=ik):: globalNb, globalNbp  !globalN basis for each element

    ! NEM_alge = NEM/3*2

    NNR1278 = ( 2*(NEL+NES) )*( 2*(NEV+NEL) )    !number of nodes in region 1&2, excluding the sharing line with region 3&4&5
    NNR1278p = (2*(NEL+NES)+1)*( 2*(NEV+NEL)+1 )    !number of nodes in region 1&2, including the sharing line with region 3&4&5
    NNX34 = 4*NEL + 1
    NNX3456 = NNX34 + 2*(NEV+NES)
    NNXV = 2*NEV + 1
    NTN = NNX3456*( 2*NEM + 1 ) + NNR1278
    NTE = (2*NEL+NEV+NES)*NEM + (NEL+NES)*(NEV+NEL)
    top_node = NTN-2*NEV

    NND = (2*NEL+1)**2 + 2*NEM*NNX34
    NNV = (2*NEV+1) * (2*(NEL+NEM)+1)
    NNS = (2*NES+1) * (2*(NEL+NEM+NEV)+1)
    ED = NEL**2 + NEM*2*NEL
    EV = NEV*(NEL+NEM)
    ES = NES*(NEM+NEL+NEV)
    
    open(unit=10, file = trim(folder)//'element_distribution.dat', status = 'replace')
200 format(A,i4)
    write(10,200) 'NEL =', NEL
    write(10,200) 'NEM =', NEM
    write(10,200) 'NEM_alge =', NEM_alge
    write(10,200) 'NES =', NES
    write(10,200) 'NEV =', NEV
    write(10,200) 'rad_ele =', NEL*2+NEV+NES
    close(10)


    allocate( globalNM(NTE,9), rowNM(NTE), columnNM(NTE), WFLAG(NTE), RegN(NTE), VE(NTE) )
    !VE: vapor element
    allocate( NOPDV(NTN,3), layer(NTN), BCflagE(NTE,5), BCflagN(NTN,5) )
    allocate( Ngrid(NNR1278p,2) )!?
    !Ngrid: rowN & columnN of the node (only for Reg1&2, count together)
    allocate( VN(NTN) ) !VN = 0: droplet, 1: vapor phase, 2: free surface, 5: substrate only
    allocate( algeN(NTN), algeS(NTN) )
    globalNM(:,:) = 0
    NOPDV(:,:) = 0
    rowNM(:) = 0
    columnNM(:) = 0
    WFLAG(:) = 0
    RegN(:) = 0
    Ngrid(:,:) = 0
    VE(:) = 0
    VN(:) = 0
    layer(:) = 0
    algeN(:) = 0
    algeS(:) = 0

! if(NEV.lt.NES) then
!    write(*,*) 'NEV < NES, need a new way to number elements'
!    stop
! end if

do eleN = 1, NTE

   if ( eleN.le.(NEL+NES)*(NEL+NEV) ) then    !the element is in region 1,2,7 or 8

      !calculate rowN & columnN
      do i = 1, NEL+NES         !solve for the elements in region1&2 round by round
         if ( eleN.le.(i**2 + i*abs(NEV-NES) ) ) exit
      end do
      if(NEV.ge.NES) then
         rowNM(eleN) = i+NEV-NES
         columnNM(eleN) = i
         if( eleN .lt. ( i**2 + (i-1)*(NEV-NES-1) ) ) then
            columnNM(eleN) = eleN - ( (i-1)**2 + (i-1)*(NEV-NES) )
         else if( eleN .gt. ( i**2 + (i-1)*(NEV-NES-1) ) ) then
            rowNM(eleN) = i**2 + i*(NEV-NES) - eleN + 1
         end if
      else
         rowNM(eleN) = i
         columnNM(eleN) = i+NES-NEV
         if( eleN .lt. ( i**2 + i*(NES-NEV-1) + 1 ) ) then
            columnNM(eleN) = eleN - ( (i-1)**2 + (i-1)*(NES-NEV) )
         else if( eleN .gt. ( i**2 + i*(NES-NEV-1) + 1 ) ) then
            rowNM(eleN) = i**2 + i*(NES-NEV) - eleN + 1
         end if
      end if
           !  exit
      !    end if
      ! end do

      do locN = 1, 9

         !calculate globalN
         select case (locN)
         case(1)
            LNrowN = 3
            LNcolumnN = 1
         case(2)
            LNrowN = 2
            LNcolumnN = 1
         case(3)
            LNrowN = 1
            LNcolumnN = 1
         case(4)
            LNrowN = 3
            LNcolumnN = 2
         case(5)
            LNrowN = 2
            LNcolumnN = 2
         case(6)
            LNrowN = 1
            LNcolumnN = 2
         case(7)
            LNrowN = 3
            LNcolumnN = 3
         case(8)
            LNrowN = 2
            LNcolumnN = 3
         case(9)
            LNrowN = 1
            LNcolumnN = 3
         end select

         NrowN = 2*rowNM(eleN) + LNrowN - 2
         NcolumnN = 2*columnNM(eleN) + LNcolumnN - 2

         if(NEV.ge.NES) then
            if( (NrowN-2*(NEV-NES)) .gt. NcolumnN ) then
               NrowNp = NrowN - 2*(NEV-NES)
               globalNM(eleN, locN) = NrowNp**2 - NrowNp + 1 + 2*(NEV-NES)*(NrowNp-1) - (NrowNp-NcolumnN)
            else
               globalNM(eleN, locN) = NcolumnN**2 - NcolumnN + 1 + 2*(NEV-NES)*(NcolumnN-1) - NrowN + &
                    (2*(NEV-NES)+NcolumnN)
            end if
         else!???
            if( (NcolumnN-2*(NES-NEV)) .lt. NrowN ) then
               NrowNp = NrowN - 1
               globalNM(eleN, locN) = NrowNp**2 + 2*(NES-NEV)*(NrowNp) + (NcolumnN)
              
            else
               NcolumnNp = NcolumnN - 2*(NES-NEV)
               globalNM(eleN, locN) = NcolumnNp**2 + 2*(NES-NEV)*(NcolumnNp) - NrowN+1
            end if


         end if

         Ngrid( globalNM(eleN, locN), 1 ) = NrowN
         Ngrid( globalNM(eleN, locN), 2 ) = NcolumnN

      end do

      
      if( columnNM(eleN).le.NES ) then   !this element is in region 7 or 8 
         if( rowNM(eleN).le.NEV ) then
            RegN(eleN) = 8
         else
            RegN(eleN) = 7
         end if
      else   !this element is in region 1 or 2 
         if( rowNM(eleN).le.NEV ) then
            RegN(eleN) = 2
            ! rowNM(eleN) = rowNM(eleN) - NEV
         else
            RegN(eleN) = 1
         end if
         ! columnNM(eleN) = columnNM(eleN) - NES
      end if


   else   !in region 3, 4, 5 or 6
      
      eleNp = eleN - (NEL+NES)*(NEL+NEV)

      !calculate rowN & columnN
      if( mod( eleNp, 2*NEL+NEV+NES ) .eq. 0 ) then
         rowNM(eleN) = eleNp/(2*NEL+NEV+NES)              !no dble, rowN of elements
         columnNM(eleN) = (2*NEL+NEV+NES)
      else
         rowNM(eleN) = eleNp/(2*NEL+NEV+NES) + 1          !no dble
         columnNM(eleN) = mod( eleNp,(2*NEL+NEV+NES) )
      end if
      

      if( columnNM(eleN).le.NES ) then
         RegN(eleN) = 6
      else if( columnNM(eleN).le.NEL+NES ) then
         RegN(eleN) = 3
      else if( columnNM(eleN).le.NEL*2+NES ) then
         RegN(eleN) = 4
      else
         RegN(eleN) = 5
      end if

      if( rowNM(eleN).eq.1 .and. columnNM(eleN).le.NEL+NES ) then
         WFLAG(eleN) = WFLAG(eleN) + 1
         if( columnNM(eleN).eq.NEL+NES )  WFLAG(eleN) = WFLAG(eleN) + 1
      end if

      !calculate globalN
      globalNb = 2*NNX3456*(rowNM(eleN)-1) + 2*(columnNM(eleN)-1) + NNR1278

      do locN = 1, 9

         select case (locN)
         case(1)
            globalNM(eleN, locN) = globalNb + 1
         case(2)
            globalNM(eleN, locN) = globalNb + 2
         case(3)
            globalNM(eleN, locN) = globalNb + 3
         case(4)
            globalNM(eleN, locN) = globalNb + NNX3456 + 1
         case(5)
            globalNM(eleN, locN) = globalNb + NNX3456 + 2
         case(6)
            globalNM(eleN, locN) = globalNb + NNX3456 + 3
         case(7)
            globalNM(eleN, locN) = globalNb + 2*NNX3456 + 1
         case(8)
            globalNM(eleN, locN) = globalNb + 2*NNX3456 + 2
         case(9)
            globalNM(eleN, locN) = globalNb + 2*NNX3456 + 3
         end select
      end do

   end if

   !VN&VE
   !vapor
   if(RegN(eleN).eq.2 .or. RegN(eleN).eq.5) then
      VE(eleN) = 1
      do locN = 1, 9
         VN( globalNM(eleN, locN) ) = 1
         if( ( RegN(eleN).eq.2 .and. rowNM(eleN).eq.NEV ) .or. &
              ( RegN(eleN).eq.5 .and. columnNM(eleN).eq.NES+2*NEL+1 ) ) then
            if(locN.eq.1 .or. locN.eq.4 .or. locN.eq.7) then
               VN( globalNM(eleN, locN) ) = 2
            end if
         end if
      end do
   end if
   !substrate
   if(RegN(eleN).eq.6 .or. RegN(eleN).eq.7 .or. RegN(eleN).eq.8 ) then
      VE(eleN) = 5
      do locN = 1, 9
         if(columnNM(eleN).eq.NES) then
            if( ( RegN(eleN).eq.6 .and. (locN.eq.3.or.locN.eq.6.or.locN.eq.9 ) ) .or. &
                 ( (RegN(eleN).eq.7.or.RegN(eleN).eq.8).and.(locN.eq.7.or.locN.eq.8.or.locN.eq.9) ) )&
                 cycle
         end if
         VN( globalNM(eleN, locN) ) = 5
      end do
   end if

end do  !eleN



!BCflagE: axis, base, free, outer
BCflagE(:,:) = 0
do i = 1, NTE
   !axis
   if( (RegN(i).eq.3.or. RegN(i).eq.4.or. RegN(i).eq.5.or. RegN(i).eq.6) .and. rowNM(i).eq.NEM ) then
      BCflagE(i,1) = 1
   end if
   !base
   if( (RegN(i).eq.1 .or. RegN(i).eq.2) .and. columnNM(i).eq.NES+1 ) then
      BCflagE(i,2) = 1
   end if
   if( RegN(i).eq.3 .and. columnNM(i).eq.NES+1 ) then
      BCflagE(i,2) = 2
   end if
   if( (RegN(i).eq.7 .or. RegN(i).eq.8) .and. columnNM(i).eq.NES ) then
      BCflagE(i,2) = 11
   end if
   if( RegN(i).eq.6 .and. columnNM(i).eq.NES ) then
      BCflagE(i,2) = 21
   end if
   !free
   if( (RegN(i).eq.1 .and. rowNM(i).eq.NEV+1).or. (RegN(i).eq.4 .and. columnNM(i).eq.2*NEL+NES) ) then
      BCflagE(i,3) = 1
   end if
   if( (RegN(i).eq.2 .and. rowNM(i).eq.NEV) .or.(RegN(i).eq.5 .and. columnNM(i).eq.2*NEL+NES+1) ) then
      BCflagE(i,3) = 2
   end if
   !outer circle
   if( ( (RegN(i).eq.2.or.RegN(i).eq.8) .and. rowNM(i).eq.1).or. (RegN(i).eq.5 .and. columnNM(i).eq.2*NEL+NEV+NES) ) then
      BCflagE(i,4) = 1
   end if
   !substrate base
   if( (RegN(i).eq.6 .or. RegN(i).eq.7 .or. RegN(i).eq.8) .and. columnNM(i).eq.1 ) then
      if(RegN(i).eq.6) then 
         BCflagE(i,5) = 2
      else 
         BCflagE(i,5) = 1
      end if
   end if
end do

!BCflagN: axis, base, free, outer
BCflagN(:,:) = 0
do i = 1, NTE
   do j = 1, 9
      
   !axis
   if( (RegN(i).eq.3.or. RegN(i).eq.4.or. RegN(i).eq.5.or. RegN(i).eq.6) .and. rowNM(i).eq.NEM .and. &
        ( j.eq.7 .or. j.eq.8 .or. j.eq.9 ) ) then
      BCflagN( globalNM(i,j), 1) = 1
   end if

   !base
   if( (RegN(i).eq.1 .or. RegN(i).eq.2) .and. columnNM(i).eq.NES+1 .and. ( j.eq.1 .or. j.eq.2 .or. j.eq.3 ) ) then
      if(rowNM(i).eq.NEL+NEV .and. j.eq.1) cycle
      BCflagN( globalNM(i,j), 2 ) = 1
   else if( RegN(i).eq.3 .and. columnNM(i).eq.NES+1 .and. ( j.eq.1 .or. j.eq.4 .or. j.eq.7 ) ) then
      BCflagN( globalNM(i,j), 2 ) = 2
      if( WFLAG(i).eq.1 .and. j.eq.1 )   BCflagN( globalNM(i,j), 2 ) = 3
   end if
   
   !free
   if( VN(globalNM(i,j)).eq.2 ) then
      BCflagN( globalNM(i,j), 3) = 1
      if( RegN(i).eq.2 .and. columnNM(i).eq.NES+1 .and. j.eq.1 )  BCflagN( globalNM(i,j), 3) = 3!???
   else if( RegN(i).eq.8 .and. rowNM(i).eq.NEV .and. ( j.eq.1 .or. j.eq.4 .or. j.eq.7 ) ) then
      BCflagN( globalNM(i,j), 3) = 2 
   end if

   !outer circle
   if(( (RegN(i).eq.2.and.rowNM(i).eq.1) .or. ( RegN(i).eq.5 .and. columnNM(i).eq.2*NEL+NEV+NES ) ) &
        .and. ( j.eq.3 .or. j.eq.6 .or. j.eq.9 ) )then
      if(columnNM(i).eq.NES+1 .and. j.eq.3) cycle
      BCflagN( globalNM(i,j), 4) = 1
   else if(RegN(i).eq.8 .and. rowNM(i).eq.1 .and. ( j.eq.3 .or. j.eq.6 .or. j.eq.9) ) then
      BCflagN( globalNM(i,j), 4) = 2
      if( columnNM(i).eq.NES .and. j.eq.9 ) BCflagN( globalNM(i,j), 4) = 3
   end if

   !substrate base
   if( (RegN(i).eq.7 .or. RegN(i).eq.8) .and. columnNM(i).eq.1 .and. ( j.eq.1 .or. j.eq.2 .or. j.eq.3 ) ) then
      if(rowNM(i).eq.NEL+NEV .and. j.eq.1) cycle
      BCflagN( globalNM(i,j), 5 ) = 1
   else if( RegN(i).eq.6 .and. columnNM(i).eq.1 .and. ( j.eq.1 .or. j.eq.4 .or. j.eq.7 ) ) then
      BCflagN( globalNM(i,j), 5 ) = 2
      if( WFLAG(i).eq.1 .and. j.eq.1 )   BCflagN( globalNM(i,j), 5 ) = 3
   end if

   end do
end do


!NOPD( NOPDV(:,1) ) & NOPV( NOPDV(:,2) ) & NOPS( NOPDV(:,3) )
do i = 1, NTN
   !NOPD
   if(VN(i).ne.1 .and. VN(i).ne.5) then
      NOPDV(i,1) = i
      do j = 1, i
         if(VN(j).eq.1.or.VN(j).eq.5) NOPDV(i,1) = NOPDV(i,1) - 1
      end do
   end if
   !NOPV
   if(VN(i).ne.0 .and. VN(i).ne.5) then
      NOPDV(i,2) = i
      do j = 1, i
         if(VN(j).eq.0 .or. VN(j).eq.5) NOPDV(i,2) = NOPDV(i,2) - 1
      end do
      layer(i) = mod(NOPDV(i,2),NNXV)
      if(layer(i).eq.0) layer(i) = NNXV
   end if
   !NOPS
   if(VN(i).eq.5 .or. BCflagN(i,2).ne.0) then
      NOPDV(i,3) = i
      do j = 1, i
         if(VN(j).ne.5.and.BCflagN(j,2).eq.0) NOPDV(i,3) = NOPDV(i,3) - 1
      end do
   end if

end do


!algeN
algeN(:) = 0
do i = 1, NTE
   if(RegN(i).eq.3 .or. RegN(i).eq.4 .or. RegN(i).eq.5) then
      if(rowNM(i).eq.NEM_alge) then
         do j = 7, 9
            algeN(globalNM(i,j)) = 1
         end do
      end if
   end if
end do

!algeS
algeS(:) = 0
do i = 1, NTE
   if(RegN(i).eq.6) then
      if( rowNM(i).eq.NEM_alge ) then
         do j = 7,9
            algeS(globalNM(i,j)) = 2
         end do
      ! else if( rowNM(i).eq.1 ) then
      !    do j = 1, 3
      !       algeS(globalNM(i,j)) = 1
      !    end do
      end if
   else if(RegN(i).eq.8) then
      if( rowNM(i).eq.NEV ) then
         do j = 1,7,3
            algeS(globalNM(i,j)) = 1
         end do
      end if
      do j = 1, 9
         do n = 1, vlayer-1
            if( Ngrid( globalNM(i,j),1 ) .eq. 2*(NEV-NEV1(n))+1 ) then
               algeS( globalNM(i,j) ) = 1
            end if
         end do
      end do
   end if
end do

! open( unit = 10, file = trim(folder)//'rowN.dat', status = 'replace' )
! do eleN = 1, NTE
!    write(10, '(3i8)') eleN, rowNM(eleN), columnNM(eleN)
! end do
! close(10)
! pause

!adjust rowNM&columnNM
do eleN = 1, NTE
   if(RegN(eleN).eq.6 )  rowNM(eleN) = -rowNM(eleN) + NEM+1
   if(RegN(eleN).eq.3 .or. RegN(eleN).eq.4 .or. RegN(eleN).eq.1 .or. RegN(eleN).eq.2 ) &
        columnNM(eleN) = columnNM(eleN) - NES
   if(RegN(eleN).eq.1) rowNM(eleN) = rowNM(eleN) - NEV
   if(RegN(eleN).eq.5) columnNM(eleN) = columnNM(eleN) - NES - 2*NEL
   if(RegN(eleN).eq.7 .or. RegN(eleN).eq.8 )  rowNM(eleN) = -rowNM(eleN) + NEM+NEL+NEV+1
   if(RegN(eleN).eq.2 ) rowNM(eleN) = -rowNM(eleN) + NEV+1
end do

! open(unit = 10, file = trim(folder)//'algeN.dat', status = 'replace' )
! do i = 1, NTN
!    write(10, '(2i8)' ) i, algeN(i)
! end do
! close(10)
! ! pause

! open(unit = 10, file = trim(folder)//'algeS.dat', status = 'replace' )
! do i = 1, NTN
!    write(10, '(2i8)' ) i, algeS(i)
! end do
! close(10)
! pause

! open(unit = 10, file = trim(folder)//'BCflagN.dat', status = 'replace' )
! do i = 1, NTN
!    write(10, '(6i8)' ) i, BCflagN(i,1), BCflagN(i,2), BCflagN(i,3), BCflagN(i,4), BCflagN(i,5)
! end do
! close(10)
! pause

! open(unit = 10, file = trim(folder)//'BCflagE.dat', status = 'replace' )
! do i = 1, NTE
!    write(10, '(6i8)' ) i, BCflagE(i,1), BCflagE(i,2), BCflagE(i,3), BCflagE(i,4), BCflagE(i,5)
! end do
! close(10)
! !pause

! open(unit = 10, file = trim(folder)//'NOPD.dat', status = 'replace' )
! do i = 1, NTN
!    write(10, '(5i8)' ) i, NOPDV(i,1), NOPDV(i,2), NOPDV(i,3), layer(i)
! end do
! close(10)


! open(unit = 10, file = trim(folder)//'VN.dat', status = 'replace' )
! do i = 1, NTN
!    write(10, '(2i8)' ) i, VN(i)
! end do
! close(10)

! open(unit = 10, file = trim(folder)//'VE.dat', status = 'replace' )
! do i = 1, NTE
!    write(10, '(2i8)' ) i, VE(i)
! end do
! close(10)


! open(unit = 10, file = trim(folder)//'globalNM.dat', status = 'replace' )
! do eleN = 1, NTE
!    do locN = 1, 9
!       write(10, '(i8)', advance = 'no' ) globalNM(eleN, locN) 
!    end do
!    write(10,'(A)') ' '
! end do
! close(10)
! pause

! open( unit = 10, file = trim(folder)//'RegN.dat', status = 'replace' )
! do eleN = 1, NTE
!    write(10, '(2i8)') eleN, RegN(eleN)
! end do
! close(10)
! pause

! open( unit = 10, file = trim(folder)//'rowN.dat', status = 'replace' )
! do eleN = 1, NTE
!    write(10, '(3i8)') eleN, rowNM(eleN), columnNM(eleN)
! end do
! close(10)
! pause

! open( unit = 10, file = trim(folder)//'WFLAG.dat', status = 'replace' )
! do eleN = 1, NTE
!    write(10, '(2i8)') eleN, WFLAG(eleN)
! end do
! close(10)

! open( unit = 10, file = trim(folder)//'NrowN.dat', status = 'replace' )
! do i = 1, NNR1278p
!    write(10, '(3i8)') i, Ngrid(i,1), Ngrid(i,2)
! end do
! close(10)

! pause


    allocate( rNOP(NTN,4,2) )
    rNOP(:,:,:) = 0
    Do i = 1, NTE, 1
       Do j = 1, 9, 1
          l = 0
          k = 1
152       if (l.eq.0.and.k.le.4) then
             if (rNOP(globalNM(i,j),k,1).eq.0) then
                rNOP(globalNM(i,j),k,1) = i
                rNOP(globalNM(i,j),k,2) = j
                l = 1
             end if
             k = k + 1
             goto 152

          end if
       end do
    end do
    ! do i = 1, NTN
    !    do j = 1, 4
    !       if(rNOP(i,j,1).eq.0) then
    !          rNOP(i,j,1) = rNOP(i,j-1,1)
    !          rNOP(i,j,2) = rNOP(i,j-1,2)
    !       end if
    !    end do
    ! end do

    ! open(unit = 10, file = trim(folder)//'rNOP.dat', status = 'replace')
    ! do i = 1, NTN, 1
    !    write(10,'(A, i8, A)') '----------------------globalN:', i, '------------------------'
    !    do j = 1, 4
    !       write(10,'(2i8)') rNOP(i,j,1), rNOP(i,j,2)
    !    end do
    ! end do
    ! close(10)
    ! pause

    return
  end subroutine NOP



subroutine variableN
  implicit none

integer(kind=ik):: i, j,k,l


  allocate( MDF(NTN), NOPP(NTN), MDFd(NTN), PN(NTN) )
  MDF = 0
  MDFd = 0
  PN = 0
  
! MDF = 2
    do i = 1, NTE
       do j = 1, 9

          if(VE(i).eq.0) then

             if ( j.eq.1 .or. j.eq.3 .or. j.eq.7 .or. j.eq.9 ) then
                MDF( globalNM(i,j) ) = NNVar1
                PN( globalNM(i,j) ) = 1
             else
                MDF( globalNM(i,j) ) = NNVar1-1
             end if
             
             ! if( VN( globalNM(i,j) ).eq.2 ) then
             !    MDF( globalNM(i,j) ) = MDF( globalNM(i,j) ) + 1
             ! end if

          end if


          ! if( VN(globalNM(i,j)).eq.1 ) then
          !    MDF( globalNM(i,j) ) = NNvar2
          ! end if

       end do
    end do

    do i = 1, NTN
       !vapor
       if(VN(i).eq.2) MDF(i) = MDF(i) + 1 !c
       if(VN(i).eq.1) MDF(i) = NNVar2 !r,z,c
       !substrate
       if(BCflagN(i,2).eq.1 .and. VN(i).eq.1) MDF(i) = MDF(i) + 1 !T
       if(VN(i).eq.5) MDF(i) = NNVar3 !r,T
    end do


    MDFd(:) = MDF(:)

    call NOPP_define

    !define iBW
    iBW1 = NOPP( (NNR1278p + 5) + NNX3456*5 ) + MDF( (NNR1278p + 5) + NNX3456*5 ) - NOPP(NNR1278p + NNX3456 + 1) + 1
    !iBW2 = NOPP(NNR1278p) + MDF(NNR1278p) - NOPP( globalNM( (NEL-1)**2+(NEL-2)*(NEV-1) ,3) + 2*(NEV+NEL-4) ) + 1
    iBW = iBW1 !max(iBW1, iBW2)


! open(unit = 10, file = trim(folder)//'MDF.dat', status = 'replace')
! do i = 1, NTN
!    write(10, '(3i8)') i, MDF(i), PN(i)
! end do
! close(10)
! pause


! open(unit = 10, file = trim(folder)//'NOPP.dat', status = 'replace')
! do i = 1, NTN
!    write(10, '(2i8)') i, NOPP(i)
! end do
! close(10)
! pause
    
    return

  end subroutine variableN
  


  subroutine NOPP_define!(x)
    implicit none
    integer(kind=ik):: i,j
    !integer(kind=ik),allocatable:: x(:)

    ! allocate( x(NTN) )
    NOPP(:) = 0
    do i = 1, NTN, 1
       do j = 1, i-1, 1
          NOPP(i) = NOPP(i) + MDF(j)
       end do
       NOPP(i) = NOPP(i) + 1
    end do


    NVar = NOPP(NTN) + MDF(NTN) - 1
    
    return
  end subroutine NOPP_DEFINE





  real(kind=rk) function damfac(x)
    implicit none

    real(kind=rk):: x

    ! if(x.gt.50.0_rk) then
    !    damfac = 0.05_rk

    ! else 
    if(x.gt.5.0_rk) then
       damfac = 0.1_rk

    else if(x.gt.3.0_rk) then
       damfac = 0.3_rk

    else if(x.gt.0.5_rk) then
       damfac = 0.6_rk

    else
       damfac = 1.0_rk
    end if
    return
  end function damfac




real(kind=rk) function gaussian_quadrature(f) 
  implicit none

  !gaussian quadrature integration in 0-1 with 3 mesh points.
  ! f(3,3)(integrand) needs to be defined:
  !use "do" loop to give value to f(3): eg, f(i) = a(i)**2 (This stands for f(x) = x**2 )

  integer(kind=ik):: i, j
  real(kind=rk):: f(3,3)
  real(kind=rk), parameter, dimension(3):: w = (/0.2777777777777777777777777777778_rk,&
       0.44444444444444444444444444444444_rk, 0.2777777777777777777777777777778_rk/)

  gaussian_quadrature = 0.0_rk
  do i = 1,3,1
     do j = 1,3,1
     gaussian_quadrature = gaussian_quadrature + w(i)*w(j)*f(i,j)
     end do
  end do
  return
end function gaussian_quadrature

real(kind=rk) function gaussian_quadrature_1d(f)
  implicit none

  !gaussian quadrature integration in 0-1 with 3 mesh points.
  ! f(3)(integrand) needs to be defined
  !use "do" loop to give value to f(3): eg, f(i) = a(i)**2 (This stands for f(x) = x**2 )

  integer(kind=ik):: i
  real(kind=rk):: f(3)
  real(kind=rk), parameter, dimension(3):: w = (/0.2777777777777777777777777777778_rk,&
       0.44444444444444444444444444444444_rk, 0.2777777777777777777777777777778_rk/)

  gaussian_quadrature_1d = 0.0_rk
  do i = 1,3,1
     gaussian_quadrature_1d = gaussian_quadrature_1d + w(i)*f(i)
  end do

  return
end function gaussian_quadrature_1d






end module NOP_mod

