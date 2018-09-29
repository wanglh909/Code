!-*- mode: f90;-*-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module front_mod
  use kind, only: rk, ik, ths
  use data, only: rmax, zmax, umax, vmax, Tmax, cpmax, pmax, cmax, &
       VN, initial_vapor_solved, BCflagN, PN, no_vapor, &
       Nr, Nz, NT, Nu, Nv, Ncp, Np, NTs
  use omp_lib

  integer(kind=ik):: check_pivot=0
  
  !Search for "USER_SPECIFIED" to find parts in need of modification
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!USER_SPECIFIED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  !Must be suuplied from a user module
  !rk: real precision, parameter
  !ik: integer size, parameter
  !ths: Number of threads, parameter
  !NVar: Number of unknowns
  !NN: Number of Nodes
  !NE: Number of elements
  !iBW: Bandwidth (warning will be changed by init_front
  !threads: Number of threads used
  !MDF(NN): Unkowns at given node
  !NOP(NE,9): Local to global nodes
  !NOPP(NN): first var at node
  !load(NVar): Solution goes here
  !s_mode: Set to 0 if solving full dynamics and 1 if only mesh 
  !(must be zero if solve mesh by putting BC's on all velocity and pressure vars)
  !bas: number of nodes in element (9)

  
  !rNOP(NN,4,2): Reverse NOP, create if you don't have it using below code
  
  !Section to create rNOP
  !rNOP(:,:,:) = 0
  ! Do i = 1, NE, 1
  !      Do j = 1, 9, 1
  !         l = 0
  !         k = 1
  ! 152     if (l.eq.0.and.k.le.4) then
  !            if (rNOP(NOP(i,j),k,1).eq.0) then
  !               rNOP(NOP(i,j),k,1) = i
  !               rNOP(NOP(i,j),k,2) = j
  !               l = 1
  !            end if
  !            k = k + 1
  !            goto 152
  
  !         end if
  !      end do
  !   end do

  !Must have a callable assembly subroutine:
  !call assemble_local(ele,local,loc,NOPPl,NB,id)
  !ele: element to be assembled
  !local(NB,NB): local matrix (Warning assemble as (j,i), we normally assemble i,j i = row, j = col)
  !NB: Integer, defines sizes
  !loc(NB): local rhs
  !NOPPl(bas): Local NOPPl, must be filled and returned
!!!!!!!!!!!!!!!!!!!!USER_SPECIFIED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  abstract interface
     subroutine assembler_type(ele,amat_local,loc,noppl_o,NB,id,dum)
       use kind
       integer(kind=ik) :: ele, NB, id, dum, noppl_o(9)
       real(kind=rk) :: amat_local(NB,NB), loc(NB), alcl_o(NB,NB)

     end subroutine assembler_type
  end interface

  abstract interface
     subroutine void_type()
       
     end subroutine void_type
  end interface
  
  abstract interface
     subroutine exclude_type(ele,inc)
       use kind
       integer(kind=ik) :: ele
       logical :: inc

     end subroutine exclude_type
  end interface

  abstract interface
     subroutine find_var_type(global_var,element,proc_id)
       use kind
       integer(kind=ik) :: global_var, element, proc_id
     end subroutine find_var_type
  end interface

  procedure (assembler_type), pointer :: assembler => null () !Gets an assembled element
  procedure (void_type), pointer :: associater => null ()     !Associates pointers to external data 
  procedure (void_type), pointer :: custom_order => null ()   !If left null order is 1 to NE
  procedure (exclude_type), pointer :: excluder => null () !Excluder, default include all
  procedure (find_var_type), pointer :: var_finder => null () !If NaN execute this
  logical :: debug_NAN = .false.		
  logical :: load_balance_flag = .TRUE., SWAP_LOCAL_IJ = .TRUE., REVERSE_FLAG = .TRUE.
  integer(kind=ik) :: seed = 1, THREADS_FRONT = -1
  
  real(kind=rk), allocatable :: loadc(:), LT(:,:), load_dum(:)
  integer(kind=ik), allocatable ::  IT(:,:), lt_i(:), NODf(:,:), NOA(:), NCN(:),&
       ele_list(:,:), lt_i2(:), dm_list(:,:), dm_type(:), Etest(:), eprs(:,:,:), dprs(:,:,:), piv(:),&
       BNOP(:), sBNOP(:)
  integer(kind=ik) :: NB, LT_l, m2, m1, mDM, rad_ele, z_ele, iDM, ths2, ths3, ptn_g, ptn2_g, NSi,&
       prm, prm2, max_prs, s_frac, f_cut, fs_max, nDM(2), sDM(2)
  real(kind=rk) :: tr

  !Pointers
  integer(kind=ik), pointer :: NVar, NN, NE, s_mode, bas, MDF(:), NOPP(:), DNOP(:), NOP(:,:), rNOP(:,:,:)
  real(kind=rk), pointer :: load(:)

  integer(kind=ik), parameter :: DUAL=1, SINGLE=2, DOMAINS=3

  !LT: Stores the pivotal rows
  !IT: Stores the information on how to use LT to back substitute
  !loadc: b vector used during solution, but not where the answer goes
  !iPg: global pointers i
  !jPg: global pointers j
  !rhead: Where the local matrix goes rowwise
  !chead: Where the local goes column wise
  !NB: Number of variables in a local block
  !front: the front
  !Lt_l: Max size of any given fully summed pivotal row to be held
  
  !Layout of sub domains
  !Arrows indicate the direction that elimination  occurs in
  !Middle domains may also be subdivided into smaller domains and 
  !eliminated via nested dissection
  !----------------------------------------------------------------------------------------------|
  !               X|--->      |           |            |            |--->        |       * <-----|
  !                |    7...  |           |            |            |   iDM-1    |               |
  !                |          |           |            |            |            |               |
  !                |----------|-----------|------------|------------|------------|               |
  !                |          |           |            |            |            |               |
  !                |    6     |           |            |            |            |               |
  !                |          |           |            |            |            |               |
  !                |----------|---------- |------------|------------|------------| -             |
  !                |          |           |            |            |            | |             |
  !                |    5     |           |            |            |            | |=ptn2        |
  !                |          |           |            |            |            | |             |
  !       d1       |----------|-----------|------------|------------|------------| -             |
  !                |          |           |            |            |            |       iDM     |
  !                |    4     |           |            |            |            |               |
  !                |          |           |            |            |            |               |
  !                |----------|-----------|------------|------------|------------|               |
  !                |          |           |            |            |            |               |
  !                |    3     |           |            |            |            |               |
  !                |          |           |            |            |            |               |
  !                |----------|-----------|------------|------------|------------|               |
  !                |          |           |            |            |            |               |
  !                |    2     |           |            |            |            |               |
  !                |          |           |            |            |            |               |
  !------>*________|Z--->_____|--->______ |____________|____________|____________|Y______________|
  !                                       |------------|
  !                                             =ptn
  !X = end_ele_d1 
  !Y = end_ele_dn
  !Z = first_ele_di
  !*Default setup is lazy domain 1 and iDM numbering (Number from bottom up)
  !Must manully create these ordering for max efficiency

  
  
  
  
  
contains

  subroutine init_front(init)
    !init = 1 for first entry, 2 for changing setup after first, 0 for exit
    !ptn = size of sub rows in middle (-1 for overide single, -2 for overide double, -3 optimal nat/nest)
    !ptn2 = size of sub columns in middle

    !Recommended ptn2 = rad_ele/3ish (want three regions close to equal, 
    !so for 28 pick 10 for 2 regions of 10 and 1 of 8)
    !ptn = 4-7, 5 in general but test it out
    !If you intend to use nested ptn = ptn2 is better
    implicit none
    integer(kind=ik) :: init, h, i, j, k, l, m, n, o, p, q, r, g, last_e, last_n, flag1, l_i,&
         ele, ptn, rp, zp, ptn2, end_ele_d1, end_ele_dn, first_ele_di, e_move, iRG, elcheck
    logical :: inc

    !Deallocate arrays
    if (allocated(NODF)) deallocate(loadc,LT,IT,lt_i,load_dum,NODf,NOA,ele_list,dm_list,lt_i2,dm_type,&
         Etest,NCN,piv,BNOP,sBNOP)

    if (init.eq.0) return

    allocate(piv(ths))

    if(.not.ASSOCIATED(associater)) then
       write(*,*) 'Error in multifront: associater not set'
       stop
    end if
    call nullifier()
    call associater()
    call check_associated()
    
    call omp_set_num_threads(THREADS_FRONT)
    tr = 0.0_rk

    !Number of vars in an element
    allocate(NCN(NE))
    NCN(:) = 0
    do j = 1, NE, 1
       do i = 1, bas, 1
          NCN(j) = NCN(j) + MDF(NOP(j,i))

       end do
    end do




    !Find max NCN => NB
    NB = 0
    do i = 1, NE, 1
       if (NCN(i).gt.NB) NB = INT(NCN(i),ik)
    end do



    !Override (2 for dual domain, 1 for single domain)
    if(ths.le.0) then
       write(*,*) 'Error in multifront: ths must be greater than 0'
       stop
    end if
    
    select case(init)
    case(SINGLE)
       iDM = 1
    case(DUAL)
       iDM = 2
       if(2.gt.ths) then
          write(*,*) 'Error in multifront: For DUAL mode ths must be >= 2'
          stop
       end if
    case(DOMAINS)
       iDM = 3
       if(THREADS_FRONT.gt.ths) then
          write(*,*) 'Error in multifront: Front has more threads specified than ths'
          stop
       else if(THREADS_FRONT.le.0) then
          write(*,*) 'Error in multifront: Number of threads has not been set'
          stop
       end if
    case default
       write(*,*) 'Error in multifront: Invalid solver type'
    end select

    !Single domain standard order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ORDERING FOR SINGLE FRONT!!!!!!!!!!!!!!!!!!!!!!!!

    !Default takes the order in which your elements are already ordered
    !But this can be changed by reordering ele_list
    !ele_list(ele. #, domain #)
    !ele. # referes to the order in which you number your elements already
    !domain # is always 1 here since single front solves one solution domain
    !with all the elements

    !Example: Assembling backwards
    ! k = 0
    ! do i = NE, 1, -1
    !    call exclude_check(i,inc)
    !    if (inc) then
    !       k = k + 1
    !       ele_list(k,1) = i
    !    end if
    ! end do

    !Notice the use of exclude check, always include this check
    !It will return .TRUE. if element should be assembled
    !What elements get turned off is determined by that subroutine
    !Default is everything is true

    !Default numbering
    if (iDM.eq.1.or.iDM.eq.2) then

       allocate(loadc(NVar),load_dum(NVar),lt_i(iDM+1),dm_type(iDM),Etest(NE),&
            NODf(NVar,iDM+1),NOA(NVar),ele_list(NE+1,iDM),dm_list(iDM,2),lt_i2(2+1),&
            BNOP(0:NE),sBNOP(0:NE))
       ele_list = 0
       k = 0
       BNOP(0) = 0
       sBNOP(0) = 0
       BNOP(1:NE) = 1
       sBNOP(1:NE) = 1

       if(ASSOCIATED(custom_order)) then
          call custom_order()
       else
          k = 0
          do i = 1, NE, 1
             call excluder(i,inc)
             if (inc) then
                k = k + 1
                ele_list(k,1) = i
             end if
          end do         
       end if

    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ORDERING FOR DUAL FRONT!!!!!!!!!!!!!!!!!!!!!!!!
    !This will take the order which is specified for single front and evenly split the elements
    !between the two solution domains
    !Don't modify if you single front is correct, multifront will automatically load
    !from here if the regions have disparate solve times
    if (iDM.eq.2) then !Dual domain standard order, and reversed standard



       !Loop that counts how many elements are in th single front domain (p)
       p = 0
       do i = 1, NE, 1
          if (ele_list(i,1).eq.0) exit
          call excluder(ele_list(i,1),inc)
          if (inc) p = p + 1
       end do

       e_move = p/2 !number of elements to move
       j = 0        !number of elements currently in domain o (empty)
       k = p        !number of elements currently in domain p
       o = 2        !Domain receiving elements
       p = 1        !Domain giving elements

       do i = 1, e_move, 1
          if (k.lt.10) then
             !If giving region has less than 10 stop
             write(*,*) 'Less than 10 elements'
             exit
          end if
          ele_list(j+1,o) = ele_list(k,p)
          BNOP(ele_list(j+1,o)) = o
          sBNOP(ele_list(j+1,o)) = o
          ele_list(k,p) = 0
          k = k - 1
          j = j + 1


       end do

    end if
    if (iDM.gt.2) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!USER_SPECIFIED!!!!!!!!!!!!!
       call multifront2_domains()

    end if

    call ETESTer()

    !NOA stores number of appearance per variable
    NOA(:) = 0

    do i = 1, NN, 1
       do j = 1, 4, 1
          if (rNOP(i,j,1).eq.0) cycle
          call excluder(rNOP(i,j,1),inc)
          if (.not.inc) cycle
          do k = 0, MDF(i)-1, 1
             NOA(NOPP(i)+k) = NOA(NOPP(i)+k) + 1
          end do
       end do
    end do

    call determine_offsets()
    !OVERRIDE FS_MAX HERE
    !fs_max =
    write(*,*) 'FS_MAX:', fs_max
    allocate(LT(fs_max,NVar),IT(3+fs_max,NVar))


  end subroutine init_front

   subroutine multifront2_domains()
    implicit none

    abstract interface
       function x_next (ele)
         use kind
         integer(kind=ik) :: x_next
         integer(kind=ik), intent (in) :: ele
       end function x_next
    end interface

    procedure (x_next), pointer :: r_next => null ()
    procedure (x_next), pointer :: r_prev => null ()
    procedure (x_next), pointer :: z_next => null ()
    procedure (x_next), pointer :: z_prev => null ()
    procedure (x_next), pointer :: block_next1 => null ()
    procedure (x_next), pointer :: block_next2 => null ()
    procedure (x_next), pointer :: temp => null ()
    integer(kind=ik), allocatable :: block_list(:,:), sub_blocks(:,:), block_dms(:,:), temp_elist(:)
    integer(kind=ik) :: blocks
    integer(kind=ik) :: blk, ele, xi_ele, eta_ele, xi_size, eta_size, rp, zp, hold_block, hold_eta,&
         xi_count, eta_count,r_size, z_size, block_size1, block_size2, step_ii, step_jj
    integer(kind=ik) :: i, j, k, DM, next, ii, jj, eta_num, xi_num, hold_block_eta, num_r
    logical, allocatable  :: REVERSED(:)
    
    call subdivide_domain(blocks,block_list)
    block_list(blocks+1,1) = NE+1
    i = 0
    do blk = 1, blocks, 1
       i = i + block_list(blk,2)*block_list(blk,3)
    end do
    write(*,*) 'NE in blocks:', i, 'NE total:', NE
 
    allocate(sub_blocks(blocks,4),block_dms(blocks,2),REVERSED(blocks))
    iDM = 0
    do blk = 1, blocks, 1

       xi_ele = block_list(blk,2)
       eta_ele = block_list(blk,3)
       !if(eta_ele < xi_ele) then
       !   sub_blocks(blk,1) = -1
       !   iDM = iDM + 1
       !   cycle
       !end if

       num_r = 3
       if(xi_ele.le.eta_ele.or.eta_ele.lt.6) then
          REVERSED(blk) = .FALSE.
          
          xi_size = xi_ele/(num_r-1)
          !if (mod(xi_ele,2).ne.0) then
          do
             xi_size = xi_size - 1
             !write(*,*) xi_size, xi_ele/xi_size, mod(xi_ele,xi_size)
             if (xi_ele/xi_size.ge.num_r) exit
          end do
          if (mod(xi_ele,num_r).ne.0) xi_size = xi_size + 1
          !end if
          eta_size = (xi_ele - 2* xi_size)/2 + mod((xi_ele - 2* xi_size),2)
          if (mod((xi_ele - 2* xi_size),2).eq.0) eta_size = eta_size + 1
          !write(*,*) 'SELECTED PTNS FOR BLOCK', blk, ':', eta_size, xi_size
          !write(*,*) 'From (radial, axial):', '(', xi_ele, ',', eta_ele, ')'

          !Size of radial and axial partitions
          rp = xi_ele/xi_size
          zp = eta_ele/eta_size
       else
          REVERSED(blk) = .TRUE.
          eta_size = eta_ele/(num_r-1)
          !if (mod(xi_ele,2).ne.0) then
          do
             eta_size = eta_size - 1
             !write(*,*) xi_size, xi_ele/xi_size, mod(xi_ele,xi_size)
             if (eta_ele/eta_size.ge.num_r) exit
          end do
          if (mod(eta_ele,num_r).ne.0) eta_size = eta_size + 1
          !end if
          xi_size = (eta_ele - 2* eta_size)/2 + mod((eta_ele - 2* eta_size),2)
          if (mod((eta_ele - 2* eta_size),2).eq.0) xi_size = xi_size + 1
          !write(*,*) 'SELECTED PTNS FOR BLOCK', blk, ':', eta_size, xi_size
          !write(*,*) 'From (radial, axial):', '(', xi_ele, ',', eta_ele, ')'

          !Size of radial and axial partitions
          rp = xi_ele/xi_size
          zp = eta_ele/eta_size
       end if
       
       
       !Determine number of domains
       ! if (mod(xi_ele,xi_size).ne.0) rp = rp + 1
       ! iDM = rp*zp
       ! if (mod(eta_ele,eta_size).ne.0) then
       !    iDM = iDM + rp
       ! end if
       ! if (mod(xi_ele,xi_size).ne.0) rp = rp - 1
        if (mod(xi_ele,xi_size).ne.0) rp = rp + 1
        if (mod(eta_ele,eta_size).ne.0) zp = zp + 1

        block_dms(blk,1) = iDM
        block_dms(blk,2) = rp*zp
        iDM = iDM + rp*zp
        
       
       hold_block = ele
       hold_eta = ele
       
       sub_blocks(blk,1) = xi_size
       sub_blocks(blk,2) = eta_size
       sub_blocks(blk,3) = rp
       sub_blocks(blk,4) = zp
       !write(*,*) 'Num. Blocks:', rp, zp
    end do
    !pause
    allocate(loadc(NVar),load_dum(NVar),lt_i(iDM+1),dm_type(iDM),Etest(NE),&
         NODf(NVar,iDM+1),NOA(NVar),ele_list(NE+1,iDM),dm_list(iDM,2),lt_i2(2+1),&
         BNOP(NE),sBNOP(0:NE))
    BNOP(:) = 0
    sBNOP(:) = 0
    
    ele_list = 0
    DM = 1
    do blk = 1, blocks, 1

       if(sub_blocks(blk,1).eq.-1) then
          k = 0
          write(*,*) iDM
          do i = block_list(blk,1), block_list(blk+1,1)-1, 1
             write(*,*) i
             call add_ele(k,DM,i,blk)
          end do
          DM = DM + 1
          cycle
       end if

       hold_block = block_list(blk,1)
       hold_block_eta = hold_block
       !xi_count = 1
       !eta_count = 0
       xi_ele = block_list(blk,2)
       eta_ele = block_list(blk,3)
       !xi_size = sub_blocks(blk,1)
       !eta_size = sub_blocks(blk,2)
       xi_num = sub_blocks(blk,3)
       eta_num = sub_blocks(blk,4)

       if(.not.REVERSED(blk)) then
          block_next1 => xi_next
          block_next2 => eta_next
          step_ii = eta_num
          step_jj = xi_num
       else
          block_next1 => eta_next
          block_next2 => xi_next
          step_ii = xi_num
          step_jj = eta_num
       end if

       !write(*,*) 'blk:', blk
       do ii = 1, step_ii, 1
          !write(*,*) 'ii:', ii
          do jj = 1, step_jj, 1
             !write(*,*) 'jj:', jj
!!!!!!!!!!!!!!!!!!!!DO BLOCK!!!!!!!!!!!!!!!!!
!!!!TBD tweak for optimal order
            

             xi_size = sub_blocks(blk,1)
             eta_size = sub_blocks(blk,2)

             
             
             if(.not.REVERSED(blk)) then

                if(jj.eq.xi_num) then
                   xi_size = xi_ele - (xi_num-1)*xi_size
                end if

                if(ii.eq.eta_num) then
                   eta_size = eta_ele - (eta_num-1)*eta_size
                end if

                block_size1 = xi_size
                block_size2 = eta_size

             else
                
                if(ii.eq.xi_num) then
                   xi_size = xi_ele - (xi_num-1)*xi_size
                end if

                if(jj.eq.eta_num) then
                   eta_size = eta_ele - (eta_num-1)*eta_size
                end if

                block_size1 = eta_size
                block_size2 = xi_size
                
             end if


            
             
             if(jj.eq.xi_num) then
                ele = hold_block
                do i = 1, xi_size - 1, 1
                   ele = xi_next(ele)
                end do
                               
                r_next => xi_prev
                !r_prev => xi_next
                r_size = xi_size
             else
                ele = hold_block
                hold_eta = ele
                r_next => xi_next
                !r_prev => xi_prev
                r_size = xi_size
                
             end if

             if(ii.eq.eta_num) then
                do i = 1, eta_size - 1, 1
                   ele = eta_next(ele)
                end do
                z_next => eta_prev
                !z_prev => eta_next    
                z_size = eta_size
             else
                z_next => eta_next
                !z_prev => eta_prev    
                z_size = eta_size
             end if

             if(xi_size.gt.eta_size) then
                temp => z_next
                z_next => r_next
                r_next => temp
                i = z_size
                z_size = r_size
                r_size = i
             end if
             
             k = 0
             hold_eta = ele
             do i = 1, z_size, 1
                call add_ele(k,DM,ele,blk)
                do j = 1, r_size-1, 1

                   next = r_next(ele)
                   if(next.eq.0) exit

                   ele = next
                   call add_ele(k,DM,ele,blk)

                end do

                
                next = z_next(hold_eta)
                !write(*,*) 'eta_next', next
                if(next.eq.0) exit
                ele = next
                hold_eta = ele
             end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


             do i = 1, block_size1, 1
                hold_block = block_next1(hold_block)
             end do
             DM = DM + 1
             k = 0
          end do !xi_block
          do i = 1, block_size2, 1
             hold_block_eta = block_next2(hold_block_eta)
          end do
         
          hold_block = hold_block_eta         
       end do !eta_block
          
    end do !block

    !Make sure domains connected
    if(REVERSE_FLAG) then
       do blk = 1, blocks-1, 1
          ii = block_dms(blk,1)+block_dms(blk,2)+1
          jj = 0
          !write(*,*) 'Block ', blk, 'searching for ', ii
          do i = block_dms(blk,1)+1, block_dms(blk,1)+block_dms(blk,2), 1
             j = 0
             do
                j = j + 1
                if(ele_list(j,i).eq.0) exit
                ele = ele_list(j,i)
                !write(*,*) 'check:', sBNOP(xi_next(ele)), sBNOP(xi_prev(ele)),sBNOP(eta_next(ele)), sBNOP(eta_prev(ele))
                if(sBNOP(xi_next(ele)).eq.ii.or.sBNOP(xi_prev(ele)).eq.ii&
                     .or.sBNOP(eta_next(ele)).eq.ii.or.sBNOP(eta_prev(ele)).eq.ii) then
                   jj = i
                   exit
                end if
             end do

             if(jj.ne.0) then
                exit
             end if
          end do
          if(jj.eq.0) cycle
          if(jj-block_dms(blk,1).lt.block_dms(blk,2)/2) then
             write(*,*) 'Reversing domain order in block #', blk, NE
             allocate(temp_elist(NE+1))
             j = block_dms(blk,1)+block_dms(blk,2)
             i = block_dms(blk,1)+1
             do


                temp_elist(:) = ele_list(:,j)
                ele_list(:,j) = ele_list(:,i)
                ele_list(:,i) = temp_elist(:)
                ii = 0
                do
                   ii = ii + 1
                   if(ele_list(ii,i).ne.0) sBNOP(ele_list(ii,i)) = i
                   if(ele_list(ii,j).ne.0) sBNOP(ele_list(ii,j)) = j
                   if(ele_list(ii,i).eq.0.and.ele_list(ii,j).eq.0) exit
                end do

                i = i + 1
                j = j - 1
                if(i.ge.j) exit
             end do
             deallocate(temp_elist)
          end if
       end do
    end if
       
       dm_list = 0
       k = 0
       do i = 1, iDM/2, 1
          k = k + 1
          dm_list(i,1) = i
       end do
       sDM(1) = 1
       nDM(1) = k  
       nDM(2) = k+1
       sDM(2) = iDM

       k = 0
       do i = iDM, iDM/2+1, -1
          k = k + 1
          dm_list(k,2) = i
       end do

    deallocate(sub_blocks,block_dms,REVERSED)
    
    !pause
  end subroutine multifront2_domains
  
  subroutine add_ele(k,DM,ele,blk)
    !use data, only: sol, n2, NOPP
    !use locals
    implicit none
    integer(kind=ik) :: k
    integer(kind=ik), intent(in) :: DM, ele, blk
    k = k + 1
    ele_list(k,DM) = ele
    BNOP(ele) = blk
    sBNOP(ele) = DM
    !write(*,*) 'ele, DM:', ele, DM
    if (ele.eq.0) then
       write(*,*) 'error in add_ele', ele, DM!, sol(ra+NOPP(NOP(ele,1)),n2),sol(za+NOPP(NOP(ele,1)),n2)
       pause
    end if
  end subroutine add_ele
  
  ! subroutine get_default_numbering(ele_list)
  !   implicit none
  !   integer(kind=ik) :: i, k
  !   integer(kind=ik) :: ele_list(:)
  !   logical :: inc

  !   k = 0
  !   do i = 1, NE, 1
  !      call excluder(i,inc)
  !      if (inc) then
  !         k = k + 1
  !         ele_list(k) = i
  !      end if
  !   end do

  ! end subroutine get_default_numbering

  integer(kind=ik) function find_start()
    implicit none
    integer(kind=ik) :: ele, next
    ele = seed
    
    do
       next = eta_prev(ele)
       if (next.eq.0) exit
       ele = next
    end do
    
    do
       next = xi_prev(ele)
       if (next.eq.0) exit
       if (NOP(ele,4).ne.NOP(next,6)) exit !weird bound jump
       ele = next
    end do

    find_start = ele

  end function find_start
  
  subroutine subdivide_domain(blocks, block_list)
    !use data, only: n2, sol, NOPP
    !use locals
    implicit none
    integer(kind=ik) :: xi_count_prev, xi_count, eta_count, ele, hold_eta, ele_count, block_start,&
         next
    integer(kind=ik), allocatable :: block_list(:,:)
    integer(kind=ik) :: blocks

    xi_count_prev = -1
    xi_count = 1
    eta_count = 0
    
    ele_count = 1
    block_start = find_start()
    ele = block_start
    hold_eta = ele
    blocks = 0
    !write(*,*) ele, BNOP(1,1,1)
    do

       !Check if there is a xi neighbor
       next = xi_next(ele)
       !write(*,*) ele, next, sol(ra+NOPP(NOP(ele,1)),n2),sol(za+NOPP(NOP(ele,1)),n2)
       if (next.gt.0) then
          ele = next
          xi_count = xi_count + 1
       else
          !Check if this was a block
          if (xi_count.ne.xi_count_prev.and.xi_count_prev.gt.0) then
             call add_block(blocks,block_start,xi_count_prev,eta_count,block_list)            
             block_start = hold_eta
          end if

          !increment eta
          xi_count_prev = xi_count
          eta_count = eta_count+1
          
          !find next column
          next = eta_next(hold_eta)
          ele = hold_eta
          !write(*,*) 'eta_next', next, 'from', hold_eta
          if(next.gt.0) then
             !Weird region
             if(NOP(ele,8).ne.NOP(next,2)) then
                !pause
                do
                   ele = xi_next(ele)
                   next = eta_next(ele)
                   if(NOP(ele,8).eq.NOP(next,2)) exit
                end do
                call add_block(blocks,block_start,xi_count_prev,eta_count,block_list)
             end if
             
             ele = next
             hold_eta = ele
             do
                next = xi_prev(ele)
                if (next.eq.0) exit
                if (NOP(ele,4).ne.NOP(next,6)) exit !weird bound jump
                ele = next
                !write(*,*) ele
                !pause
             end do

             if(ele.ne.hold_eta) then
                call add_block(blocks,block_start,xi_count_prev,eta_count,block_list)
                hold_eta = ele
             end if
          !No next column end of search   
          else
             call add_block(blocks,block_start,xi_count_prev,eta_count,block_list)
             exit
          end if
          if(block_start.eq.0) block_start = ele
          xi_count = 1
          !write(*,*) 'New column at:', ele, hold_eta
          !pause
       end if

    end do
    
  end subroutine subdivide_domain

  subroutine add_block(blocks,block_start,xi_size,eta_size,block_list)
    implicit none
    integer(kind=ik) :: blocks, block_start, xi_size, eta_size
    integer(kind=ik), allocatable :: block_list(:,:), temp(:,:)
    
    write(*,*) 'Found Block #', blocks+1, xi_size, 'x', eta_size, 'ele_start:', block_start
    if (allocated(block_list)) then
       allocate(temp(blocks+1,3))
       temp = block_list
       deallocate(block_list)
       allocate(block_list(blocks+2,3))
       block_list(1:blocks,:) = temp(:,:)
       deallocate(temp)
    else
       allocate(block_list(2,3))
    end if
    
    blocks = blocks + 1
    block_list(blocks,1) = block_start
    block_list(blocks,2) = xi_size
    block_list(blocks,3) = eta_size
    
    
    block_start = 0
    eta_size = 0
    xi_size = 0
    
    !pause
    
  end subroutine add_block

  integer(kind=ik) function xi_next(ele)
    implicit none
    integer(kind=ik), intent(in) :: ele

    if(rNOP(NOP(ele,6),1,1).ne.ele) then
       xi_next = rNOP(NOP(ele,6),1,1)
    else
       xi_next = rNOP(NOP(ele,6),2,1)
    end if

    return

  end function xi_next

  integer(kind=ik) function xi_prev(ele)
    implicit none
    integer(kind=ik), intent(in) :: ele
    
    if(rNOP(NOP(ele,4),1,1).ne.ele) then
       xi_prev = rNOP(NOP(ele,4),1,1)
    else
       xi_prev = rNOP(NOP(ele,4),2,1)
    end if

    return

  end function xi_prev

  integer(kind=ik) function eta_next(ele)
    implicit none
    integer(kind=ik), intent(in) :: ele
    
    if(rNOP(NOP(ele,8),1,1).ne.ele) then
       eta_next = rNOP(NOP(ele,8),1,1)
    else
       eta_next = rNOP(NOP(ele,8),2,1)
    end if

    return

  end function eta_next

  integer(kind=ik) function eta_prev(ele)
    implicit none
    integer(kind=ik), intent(in) :: ele

    if(rNOP(NOP(ele,2),1,1).ne.ele) then
       eta_prev = rNOP(NOP(ele,2),1,1)
    else
       eta_prev = rNOP(NOP(ele,2),2,1)
    end if

    return

  end function eta_prev
  
  subroutine determine_offsets()
    implicit none
    integer(kind=ik) :: i, j, k, ele, l, l_i, o, p, fs
    logical :: inc

    !write(*,*) 'Determining offsets'

    NODf = 0
    !Use NODf to find the amount of variables that will be
    !deleted in each level and sublevels
    !if (iDM.ge.1) then
    !Determine l_i offsets LEVEL 1
    l_i = 0
    lt_i(1) = 0
    fs_max = 0
   
    do l = 1, iDM, 1
       j = 0
       fs = 0
       do 
          j = j + 1
          ele = ele_list(j,l)
          if (ele.eq.0) exit
          do i = 1, 9, 1

             o = NOP(ele,i)
             p = NOPP(o)
             !Add new vars
             do k = 0, MDF(o)-1, 1
                if(NODf(p+k,l).eq.0) then
                   fs = fs + 1
                end if
                NODf(p+k,l) = NODf(p+k,l) + 1
             end do
          end do
          
          if (fs.gt.fs_max) fs_max = fs
          
          do i = 1, 9, 1
             o = NOP(ele,i)
             p = NOPP(o)
             do k = 0, MDF(o)-1, 1
                if (NODf(p+k,l).eq.NOA(p+k)) then
                   l_i = l_i + 1
                   NODf(p+k,l) = -NODf(p+k,l)
                   fs = fs - 1
                end if
             end do
          end do
       end do
       
       lt_i(l+1) = l_i
       !write(*,*) 'l_i1:', l_i
    end do
    !write(*,*) 'l_i 1:', l_i

    if (iDM.gt.2) then
       lt_i2(1) = lt_i(iDM+1)
       l = 1
       j = 1

       do
          j = j + 1
          ele = dm_list(j,l)
          if (ele.eq.0) exit
          NODf(:,ele) = NODf(:,ele) + NODf(:,ele-1)
          call check_NODF(NODf(:,ele),fs)
         
          if(fs.gt.fs_max) fs_max = fs
          do i = 1, NVar, 1
             if (NODf(i,ele).eq.NOA(i)) then
                l_i = l_i + 1
                NODf(i,ele) = -NODf(i,ele)
             end if
          end do
          !write(*,*) 'l_i:', l_i
       end do
       lt_i2(2) = l_i
       !write(*,*) 'l_i2:', l_i
       j = 1
       l = 2
       do
          j = j + 1
          ele = dm_list(j,l)
          if (ele.eq.0) exit
          NODf(:,ele) = NODf(:,ele) + NODf(:,ele+1)
          call check_NODF(NODf(:,ele),fs)
         
          if(fs.gt.fs_max) fs_max = fs
          do i = 1, NVar, 1
             if (NODf(i,ele).eq.NOA(i)) then
                l_i = l_i + 1
                NODf(i,ele) = -NODf(i,ele)
             end if
          end do
       end do
       lt_i2(3) = l_i
       !write(*,*) 'l_i2:', l_i
    end if
    lt_i2(2) = NVar+1
    
    !write(*,*) 'l_i:', l_i
    fs_max = fs_max
    !write(*,*) 'Calculated FS_MAX:', fs_max
  end subroutine determine_offsets
  
  subroutine check_NODF(NODf,fs)
    implicit none
    integer(kind=ik) :: NODf(NVar), i, fs
    fs = 0
    do i = 1, NVar, 1
       if(NODf(i).gt.0) then      
          fs = fs + 1
       end if
    end do

  end subroutine check_NODF

  subroutine exclude_default(ele,inc)
    implicit none
    integer(kind=ik) :: ele
    logical inc

    inc = .TRUE.
  end subroutine exclude_default
  
  subroutine ETESTer()
    implicit none
    integer(kind=ik) :: E(NE), i,j, ele
    logical :: frr
    E = 0
    do ele = 1, NE, 1
       do i = 1, iDM, 1
          do j = 1, NE, 1
             if (ele_list(j,i).eq.0) exit
             if (ele.eq.ele_list(j,i)) E(ele) = E(ele)+1
          end do
       end do
    end do

    do ele = 1, NE, 1
       if (E(ele).ne.1) then
          call excluder(ele,frr)
          if (frr) then
             write(*,'(A,i8,A,i8)',ADVANCE='NO') 'failed', ele,'with',E(ele)
             if(ASSOCIATED(DNOP)) then
                write(*,*) 'in Domain', DNOP(ele)
             else
                write(*,*) ' '
             end if
          end if
       end if
    end do
    write(*,*) 'ETEST DONE'

  end subroutine ETESTER

  integer(kind=ik) function grid(i,j,s_ele,r_ele)
    implicit none
    integer(kind=ik) :: i, j, s_ele, r_ele

    grid = s_ele + (i-1)*(r_ele) + j
    return

  end function grid

  subroutine multifront(L2res,t)
    implicit none
    real(kind=rk) :: L2res, t

    select case(iDM)
    case(1)
       call multifront_single(L2res,t)
    case(2)
       call multifront_dual2(L2res,t)
    case default
       call multifront3(L2res,t)
    end select
  end subroutine multifront

  subroutine multifront_single(L2res,t)
    !L2res: returns this value
    !t: returns solve time, useful for comparing different params
    !mode: determines natural(1) or nested(2) sub ordering
    !If you overrode init_front with -1 or -2 ptn then no purpose

    !DO NOT CALCULATE L2RES outside of this subroutine
    !Just feed in L2res and it will be returned
    !If you must use the variable load_dum to calculate yourself
    implicit none
    real(kind=rk) :: front_m(fs_max,fs_max), loadf_m(fs_max), loadf_dumm(fs_max), L2res, t1, t
    integer(kind=ik) :: e_n(4), iPg_m(fs_max), jPg_m(fs_max), fs_m, i, j, l_i, id, last_sides(2) = 0

    if (iDM.ne.1) then
       write(*,*) 'Error iDM != 1'
       return
    end if
    !mode = 1 is nested natural ordering
    !mode = 2 is nested dissection (ptn should equal ptn2)
    t1 = REAL(omp_get_wtime(),rk)
    piv = 0


    !Blank global arrays
    IT = 0              !Back sub array
    loadc = 0.0_rk      !working array
    load_dum = 0.0_rk   !dummy right hand side
    load = 0.0_rk       !Solution vector

    front_m(:,:) = 0.0_rk  !Single front matrix
    loadf_m(:) = 0.0_rk    !single right hand
    loadf_dumm(:) = 0.0_rk !single dummy
    iPg_m = 0              !local to global rows
    jPg_m = 0              !local to glbal col
    fs_m = 0               !Size of front
    l_i = 0                !Variables removed from front
    NODf = 0               !Records variable appearances



    !Single front
    call Front_solve(fs_m,front_m,loadf_m,loadf_dumm,iPg_m,jPg_m,l_i,1,1,1)

    !Perform back substitution
    call back_sub(1,l_i, last_sides)

    !Calculate L2res
    L2res = 0.0_rk
    do i = 1, NN   
       L2res = L2res + ( load_dum(NOPP(i)+Nr)/rmax )**2
       if(NTs.eq.2 .or. VN(i).ne.5)  L2res = L2res + ( load_dum(NOPP(i)+Nz)/zmax )**2
       if(s_mode.eq.0) then
          if(initial_vapor_solved.eq.1) then
             if( VN(i).eq.0 .or. VN(i).eq.2 .or. BCflagN(i,2).eq.1) then
                L2res = L2res + ( load_dum(NOPP(i)+NT)/Tmax )**2
                if( VN(i).eq.0 .or. VN(i).eq.2 ) then
                   L2res = L2res + ( load_dum(NOPP(i)+Nu)/umax )**2
                   L2res = L2res + ( load_dum(NOPP(i)+Nv)/vmax )**2
                   L2res = L2res + ( load_dum(NOPP(i)+Ncp)/cpmax )**2
                   if(PN(i).eq.1) L2res = L2res + ( load_dum(NOPP(i)+Np)/pmax )**2
                end if
             end if
             if(VN(i).eq.5) L2res = L2res + ( load_dum(NOPP(i)+NTs)/Tmax )**2
          end if
          if(no_vapor.eq.0 .and. ( VN(i).eq.1 .or. VN(i).eq.2 ) ) &
               L2res = L2res + ( load_dum(NOPP(i)+MDF(i)-1)/cmax )**2
       end if
    end do!??
    ! do i = 1, NVar, 1
    !    ! if(sol(i).ne.0.0_rk) then
    !    !    L2res = L2res + ( load_dum(i)/sol(i) )**2
    !    ! else
    !    !    L2res = L2res + load_dum(i)**2
    !    ! end if
    !    L2res = L2res + (load_dum(i))**2
    ! end do
    L2res = sqrt(L2res)

    !Save total runtime
    t = REAL(omp_get_wtime(),rk) - t1
    j = 0
    do i = 1, ths, 1
       j = piv(i) + j
    end do
    !write(*,*) 'Small Pivots:', j



  end subroutine multifront_single

  subroutine multifront_dual2(L2res,t)
    !L2res: returns this value
    !t: returns solve time, useful for comparing different params
    !mode: determines natural(1) or nested(2) sub ordering
    !If you overrode init_front with -1 or -2 ptn then no purpose

    !DO NOT CALCULATE L2RES outside of this subroutine
    !Just feed in L2res and it will be returned
    !If you must use the variable load_dum to calculate yourself
    use omp_lib
    implicit none
    real(kind=rk) :: front(fs_max,fs_max,iDM),  loadf(fs_max,iDM), loadf_dum(fs_max,iDM), L2res,&
         t0, t1, t2, t, t_diff, t_ele, t_i, t_d(iDM)
    integer(kind=ik) :: iPg(fs_max,iDM), jPg(fs_max,iDM), fs(iDM), e_n(4), i, j, k, l_i, id, s, &
         e_move, o, p, last_sides(2) = 0

    if (iDM.ne.2) then
       write(*,*) 'Error iDM != 2'
       return
    end if
    !mode = 1 is nested natural ordering
    !mode = 2 is nested dissection (ptn should equal ptn2)

    piv = 0


    !Blank global arrays
    IT = 0
    loadc = 0.0_rk      !working array
    load_dum = 0.0_rk   !dummy right hand side
    load = 0.0_rk       !Solution vector
    l_i = 0                !Variables removed from front
    NODf = 0               !Records variable appearances

    front = 0.0_rk
    loadf = 0.0_rk
    loadf_dum = 0.0_rk
    iPg = 0
    jPg = 0
    fs = 0

    t1 = REAL(omp_get_wtime(),rk)


    !PHASE 1: Solve each half
    !$omp parallel private(id,j,l_i,t_i) num_threads(ths2)
    !$omp do
    do j = 1, iDM, 1
       id = omp_get_thread_num()+1  !Thread id used to determine assembly arrays

       l_i = lt_i(j)                !initilaize var counter to predetermined starting point


       !Solve domain j and save front
       call Front_solve(fs(j),front(:,:,j),loadf(:,j),loadf_dum(:,j),iPg(:,j),jPg(:,j),l_i,id,j,1)
       t_d(j) = REAL(omp_get_wtime(),rk)
    end do
    !$omp end do
    !$omp end parallel


    !write(*,*) "fs's:", fs(1), fs(2)
    

    !PHASE 3: Stitch halves
    l_i = lt_i(3)
    id = 1
    call add_front(front(:,:,iDM),front(:,:,1),loadf(:,iDM),loadf(:,1),loadf_dum(:,iDM),&
         loadf_dum(:,1),iPg(:,iDM),iPg(:,1),jPg(:,iDM),jPg(:,1),fs(iDM),fs(1),iDM,1,l_i,1,1)
    !write(*,*) "fs's:", fs(1), fs(2)
    !call elim_stitch(front(:,:,1),loadf(:,1),loadf_dum(:,1),&
    !     iPg(:,1),jPg(:,1),l_i,fs(1),id)

    !write(*,*) 'fs:', fs(1), l_i

    !Perform back substitution
    call back_sub(2,l_i,last_sides)

    !Calculate L2res
    L2res = 0.0_rk
    !$omp parallel private(i) reduction(+:L2res)
    !$omp do
    do i = 1, NN   
       L2res = L2res + ( load_dum(NOPP(i)+Nr)/rmax )**2
       if(NTs.eq.2 .or. VN(i).ne.5)  L2res = L2res + ( load_dum(NOPP(i)+Nz)/zmax )**2
       if(s_mode.eq.0) then
          if(initial_vapor_solved.eq.1) then
             if( VN(i).eq.0 .or. VN(i).eq.2 .or. BCflagN(i,2).eq.1) then
                L2res = L2res + ( load_dum(NOPP(i)+NT)/Tmax )**2
                if( VN(i).eq.0 .or. VN(i).eq.2 ) then
                   L2res = L2res + ( load_dum(NOPP(i)+Nu)/umax )**2
                   L2res = L2res + ( load_dum(NOPP(i)+Nv)/vmax )**2
                   L2res = L2res + ( load_dum(NOPP(i)+Ncp)/cpmax )**2
                   if(PN(i).eq.1) L2res = L2res + ( load_dum(NOPP(i)+Np)/pmax )**2
                end if
             end if
             if(VN(i).eq.5) L2res = L2res + ( load_dum(NOPP(i)+NTs)/Tmax )**2
          end if
          if(no_vapor.eq.0 .and. ( VN(i).eq.1 .or. VN(i).eq.2 ) ) &
               L2res = L2res + ( load_dum(NOPP(i)+MDF(i)-1)/cmax )**2
       end if
    end do
    ! do i = 1, NVar, 1
    !    ! if(sol(i).ne.0.0_rk) then
    !    !    L2res = L2res + ( load_dum(i)/sol(i) )**2
    !    ! else
    !    !    L2res = L2res + load_dum(i)**2
    !    ! end if
    !    L2res = L2res + (load_dum(i))**2
    ! end do
    !$omp end do
    !$omp end parallel
    L2res = sqrt(L2res)    

    !Save total runtime
    t = REAL(omp_get_wtime(),rk) - t1
    j = 0
    do i = 1, ths, 1
       j = piv(i) + j
    end do
    !write(*,*) 'Small Pivots:', j



    t_diff = t_d(2) - t_d(1)
    !t_diff = 1.0_rk
    if (ABS(t_d(1) - t_d(2)).gt.0.01_rk*(t_d(1)+t_d(2)-2.0_rk*t1).and.load_balance_flag) then
       !write(*,*) 'Load balancing'
       !write(*,*) 'Times:', t_d(1)-t1, t_d(2)-t1
       !write(*,*) (t_d(1) - t_d(2)), 0.01_rk*(t_d(1)+t_d(2)-2.0_rk*t1)

       if (t_diff.gt.0) then
          o = 1 !receives
          p = 2 !gives
       else
          o = 2
          p = 1
       end if

       j = 0
       k = 0
       do i = 1, NE, 1
          if (j.eq.0) then
             if (ele_list(i,o).eq.0) j = i - 1
          end if
          if (k.eq.0) then
             if (ele_list(i,p).eq.0) k = i - 1
          end if
       end do
       !write(*,*) j, k, j+k, NE



       t_ele = (t_d(2) - t1)/REAL(k,rk)
       e_move = INT(ABS(t_diff/(2.0_rk*t_ele)),ik)

       do i = 1, e_move, 1
          if (k.lt.10) exit
          ele_list(j+1,o) = ele_list(k,p)
          ele_list(k,p) = 0
          k = k - 1
          j = j + 1


       end do


       !write(*,*) 'First half:', (t_d(1)-t1)
       !write(*,*) 'Second half:', (t_d(2)-t1)
       !write(*,*) 'Diff and per ele:', t_diff, t_ele
       !write(*,*) 'Move:', e_move, k+j, NE
       k = fs_max
       call determine_offsets()
       if(fs_max.gt.k) then
          deallocate(LT,IT)
          allocate(LT(fs_max,NVar),IT(3+fs_max,NVar))
       end if

    end if



  end subroutine multifront_dual2
  
  subroutine multifront3(L2res,t)
    !L2res: returns this value
    !t: returns solve time, useful for comparing different params
    !mode: determines natural(1) or nested(2) sub ordering
    !If you overrode init_front with -1 or -2 ptn then no purpose

    !DO NOT CALCULATE L2RES outside of this subroutine
    !Just feed in L2res and it will be returned
    !If you must use the variable load_dum to calculate yourself
    implicit none
    real(kind=rk) :: front(fs_max,fs_max,iDM), loadf(fs_max,iDM), loadf_dum(fs_max,iDM),&
         L2res, t1, t2, t3, t, ti, t0,  t_d(2), t_diff, t_ele
    integer(kind=ik) :: iPg(fs_max,iDM), jPg(fs_max,iDM), fs(iDM), e_n(4,ths), &
         i, j, k, l_i, id, o, p, f1, f2, f3, l_i2, l_i3, s_stop, thsb,&
         omp_sec, e_move, inc
    logical :: is_master, ASSEMBLED(0:iDM+1), MERGED(0:iDM+1), ASSEMBLING(0:iDM+1), MERGING(0:iDM+1), e
    integer(kind=ik) :: Domain_state(iDM), next_merge(2),&
         next_assemble(2), num_assembled(2), job_id, side, num_merged(2), merge_id, master_id,&
         assemble_id, l_i_master, thsl, l_i_inc, l_i_end(2), MC(iDM), AC(iDM), LC(NVar), NEW_THREADS
    integer(kind=ik), parameter :: WAIT = 0, ASSEMBLE = 1, MERGE = 2, J_CYCLE = 3, J_EXIT = 4

    if (iDM.lt.3) then
       write(*,*) 'Error iDM lt 3'
       return
    end if

    ti = REAL(omp_get_wtime(),rk)

    INQUIRE(FILE='THREADS_FRONT.dat',EXIST=e)
    if (e) then
       !write(*,*) 'EXISTS'
       Open(unit = 77,file = 'THREADS_FRONT.dat', status = 'old', action='read')
       read(77,'(i2)',IOSTAT=i) NEW_THREADS
       if(i.eq.0) then
          if(NEW_THREADS.le.12.and.NEW_THREADS.gt.0) then
             if(NEW_THREADS.ne.THREADS_FRONT) then
                if(NEW_THREADS.gt.ths) then
                   write(*,*) 'Cannot change number of threads to ', NEW_THREADS, ' as it is higher than ths'
                else              
                   write(*,*) 'CHANGING THREADS TO:', NEW_THREADS 
                   THREADS_FRONT = NEW_THREADS
                end if
             end if
          end if
       end if
       close(77)
    end if

    i = 0
    
    !$omp parallel num_threads(THREADS_FRONT)
    !$omp critical (test)
    i = i + 1
    !$omp end critical (test)
    !$omp end parallel
    if (i.ne.THREADS_FRONT) then
       write(*,*) 'CODE APPEARS TO BE COMPILED WITHOUT OpenMP, RUNNING IN SERIAL', i, THREADS_FRONT
       THREADS_FRONT = 1
       
    end if

    


    !Blank global arrays
    IT = 0
    loadc = 0.0_rk      !working array
    load_dum = 0.0_rk   !dummy right hand side
    load = 0.0_rk       !Solution vector
    l_i = 0                !Variables removed from front
    NODf = 0               !Records variable appearances
    piv = 0
    front = 0.0_rk
    loadf = 0.0_rk
    loadf_dum = 0.0_rk
    iPg = 0
    jPg = 0
    fs = 0

    !Set up assembly/merging flags
    next_merge(1) = 1
    next_assemble(1) = 1
    next_merge(2) = iDM
    next_assemble(2) = iDM
    ASSEMBLED = .FALSE.
    ASSEMBLED(0) = .TRUE.
    ASSEMBLED(iDM+1) = .TRUE.
    MERGED = .FALSE.
    MERGED(0) = .TRUE.
    MERGED(iDM+1) = .TRUE. 
    ASSEMBLING = ASSEMBLED
    MERGING = MERGED

    !If not load balancing this will stop both groups at the middle
    if (.not.load_balance_flag.and.THREADS_FRONT.gt.1) MERGING(iDM/2) = .TRUE.
    
    !$omp parallel private(l_i_inc,id,is_master,side,inc,l_i,l_i_master,job_id,merge_id,assemble_id,&
    !$omp master_id) num_threads(THREADS_FRONT)

    !Get thread id
    id = omp_get_thread_num()+1

    !Set masters
    if (id.eq.1.or.id.eq.THREADS_FRONT) then
       is_master = .TRUE.     
    else
       is_master = .FALSE.
    end if

    !Set group info
    if(id.le.THREADS_FRONT/2.or.THREADS_FRONT.eq.1) then
       side = 1
       inc = 1
       l_i_inc = 1
       if(is_master) master_id = 1
    else
       side = 2
       inc = -1
       l_i_inc = -1
       if(is_master) master_id = iDM
    end if
    if(is_master) l_i_master = lt_i2(side)
    
    !Do nothing on entry
    job_id = J_CYCLE 
    do
       !$omp critical (get_job)
       !Record what I just did
       select case(job_id)
       case(ASSEMBLE) 
          ASSEMBLED(assemble_id) = .TRUE. 
       case(MERGE)
          MERGED(merge_id) = .TRUE.
       end select

       !Get next job
       if(is_master.and.ASSEMBLED(next_merge(side))) then
          job_id = MERGE
          merge_id = next_merge(side)
          !If this is already being Merged I found the other group and should exit
          if(.not.MERGING(merge_id)) then 
             MERGING(merge_id) = .TRUE.
             next_merge(side) = next_merge(side) + inc
          else
             job_id = J_EXIT
          end if
       else            
          job_id = ASSEMBLE
          assemble_id = next_assemble(side)        
          if(.not.ASSEMBLING(assemble_id)) then
             ASSEMBLING(assemble_id) = .TRUE.
             next_assemble(side) = next_assemble(side) + inc
          else
             if(is_master) then !Masters may need to merge still
                job_id = J_CYCLE
             else               !Workers can quit here
                job_id = J_EXIT
             end if
          end if
       end if
       !$omp end critical (get_job)

       !Go and do my job
       select case(job_id)
       case(ASSEMBLE)       
          l_i = lt_i(assemble_id)
          call Front_solve(fs(assemble_id),&
                           front(:,:,assemble_id),&
                           loadf(:,assemble_id),&
                           loadf_dum(:,assemble_id),&
                           iPg(:,assemble_id),&
                           jPg(:,assemble_id),&
                           l_i,&
                           id,&
                           assemble_id,&
                           1)
       case(MERGE)
          call add_front(front(:,:,merge_id),&
                         front(:,:,master_id),&
                         loadf(:,merge_id),&
                         loadf(:,master_id),&
                         loadf_dum(:,merge_id),&
                         loadf_dum(:,master_id),&
                         iPg(:,merge_id),&
                         iPg(:,master_id),&
                         jPg(:,merge_id),&
                         jPg(:,master_id),&
                         fs(merge_id),&
                         fs(master_id),&
                         merge_id,&
                         master_id,&
                         l_i_master,&
                         l_i_inc,&
                         1)
       case(J_CYCLE)
          cycle
       case(J_EXIT)
          exit
       end select

    end do

    !$omp critical (end)
    if(is_master) then
       l_i_end(side) = l_i_master
    end if 
    !$omp end critical (end)

    if (.not.load_balance_flag.and.id.eq.1.and.THREADS_FRONT.gt.1) then
       call add_front(front(:,:,iDM/2),front(:,:,1),loadf(:,iDM/2),&
         loadf(:,1),loadf_dum(:,iDM/2),loadf_dum(:,1),iPg(:,iDM/2),&
         iPg(:,1),jPg(:,iDM/2),jPg(:,1),fs(iDM/2),fs(1),iDM/2,1,l_i_end(1),1,1)
       MERGED(iDM/2) = .TRUE.
    end if
!$omp end parallel

    
    
    !Quick error check
    do i = 1, iDM, 1
       if(.not.MERGED(i)) write(*,*) 'Merge    error:', i
       if(.not.ASSEMBLED(i)) write(*,*) 'Assemble error:', I
    end do

   

    if(THREADS_FRONT.ne.1) then !Single thread doesn't need to merge groups
       id = 1
       l_i = l_i_end(1)
       !Combine and eliminate final fronts
       call add_front(front(:,:,iDM),front(:,:,1),loadf(:,iDM),loadf(:,1),loadf_dum(:,iDM),&
            loadf_dum(:,1),iPg(:,iDM),iPg(:,1),jPg(:,iDM),jPg(:,1),fs(iDM),fs(1),iDM,1,l_i,1,1)
       do j = l_i_end(1)+1, l_i, 1
          LC(j) = LC(j) + 1
       end do
    else
       l_i = l_i_end(1)
    end if

    !Perform back substitution
    call back_sub(3,l_i, l_i_end)

    
    !Calculate L2res
    L2res = 0.0_rk
    !$omp parallel private(i) reduction(+:L2res)
    !$omp do
    do i = 1, NN   
       L2res = L2res + ( load_dum(NOPP(i)+Nr)/rmax )**2
       if(NTs.eq.2 .or. VN(i).ne.5)  L2res = L2res + ( load_dum(NOPP(i)+Nz)/zmax )**2
       if(s_mode.eq.0) then
          if(initial_vapor_solved.eq.1) then
             if( VN(i).eq.0 .or. VN(i).eq.2 .or. BCflagN(i,2).eq.1) then
                L2res = L2res + ( load_dum(NOPP(i)+NT)/Tmax )**2
                if( VN(i).eq.0 .or. VN(i).eq.2 ) then
                   L2res = L2res + ( load_dum(NOPP(i)+Nu)/umax )**2
                   L2res = L2res + ( load_dum(NOPP(i)+Nv)/vmax )**2
                   L2res = L2res + ( load_dum(NOPP(i)+Ncp)/cpmax )**2
                   if(PN(i).eq.1) L2res = L2res + ( load_dum(NOPP(i)+Np)/pmax )**2
                end if
             end if
             if(VN(i).eq.5) L2res = L2res + ( load_dum(NOPP(i)+NTs)/Tmax )**2
          end if
          if(no_vapor.eq.0 .and. ( VN(i).eq.1 .or. VN(i).eq.2 ) ) &
               L2res = L2res + ( load_dum(NOPP(i)+MDF(i)-1)/cmax )**2
       end if
    end do
    ! do i = 1, NVar, 1
    !    ! if(sol(i).ne.0.0_rk) then
    !    !    L2res = L2res + ( load_dum(i)/sol(i) )**2
    !    ! else
    !    !    L2res = L2res + load_dum(i)**2
    !    ! end if
    !    L2res = L2res + (load_dum(i))**2
    ! end do
    !$omp end do
    !$omp end parallel
    L2res = sqrt(L2res)

    !Save total runtime
    t = REAL(omp_get_wtime(),rk) - ti

    
  end subroutine multifront3

  subroutine back_sub(mode,last, last_sides)
    implicit none
    integer(kind=ik) :: i, j, k, l, mode, last, last_sides(2), id, l_i_inc(2),lvl2
    !integer(kind=ik) :: ITT(NVar)
    l_i_inc(1) = 1
    l_i_inc(2) = -1


    !For single, dual or nested single thread front
    !ITT(:) = 0
    if (mode.eq.1) then
       
       load(IT(3,last)) = loadc(IT(2,last))
       !ITT( IT(3,last) ) = ITT( IT(3,last) ) + 1
       do i = last-1, 1, -1
          do j = 1, IT(1,i), 1
             loadc(IT(2,i)) = loadc(IT(2,i)) - LT(j,i)*load(IT(j+3,i))
          end do
          load(IT(3,i)) = loadc(IT(2,i))
          !ITT( IT(3,i) ) = ITT( IT(3,i) ) + 1
       end do

    else if (mode.eq.2) then
       load(IT(3,last)) = loadc(IT(2,last))
       do i = last-1, lt_i(3)+1, -1
          do j = 1, IT(1,i), 1
             loadc(IT(2,i)) = loadc(IT(2,i)) - LT(j,i)*load(IT(j+3,i))
          end do
          load(IT(3,i)) = loadc(IT(2,i))
       end do

       !$omp parallel private(k,i,j) num_threads(2)
       !$omp do
       do k = iDM, 1, -1
          do i = lt_i(k+1), lt_i(k)+1, -1
             do j = 1, IT(1,i), 1
                loadc(IT(2,i)) = loadc(IT(2,i)) - LT(j,i)*load(IT(j+3,i))
             end do
             load(IT(3,i)) = loadc(IT(2,i))
          end do
       end do
       !$omp end do
       !$omp end parallel




    else if (mode.eq.3) then!For mode 1 parallel implementation
       ! write(*,*) last_sides(2)-1, 'to', last_sides(1)+1, 'by', -1
       ! write(*,*) last_sides(1), 'to', lt_i2(1)+l_i_inc(1), 'by',SIGN(1,lt_i2(1)+l_i_inc(1) - last_sides(1))
       ! write(*,*) last_sides(2), 'to', lt_i2(2)+l_i_inc(2), 'by',SIGN(1,lt_i2(2)+l_i_inc(2) - last_sides(2))   
      

      
       
       if(THREADS_FRONT.eq.1) then
          lvl2 = 1
          load(IT(3,last_sides(1))) = loadc(IT(2,last_sides(1)))
          !ITT( last_sides(1) ) = ITT( last_sides(1) ) + 1
          last_sides(1) = last_sides(1) - 1
          
       else
          lvl2 = 2
          !write(*,*) last_sides(2)-1, last_sides(1)+1, IT(1,i)
          do i = last_sides(2)-1, last_sides(1)+1, -1
             do j = 1, IT(1,i), 1
                loadc(IT(2,i)) = loadc(IT(2,i)) - LT(j,i)*load(IT(j+3,i))
             end do
             load(IT(3,i)) = loadc(IT(2,i))
             !ITT( i ) = ITT( i ) + 1
          end do
       end if
       !$omp parallel private(k,i,j,id) num_threads(lvl2)

       id = omp_get_thread_num() + 1

       !write(*,*) last_sides(id), lt_i2(id)+l_i_inc(id)
       do i = last_sides(id), lt_i2(id)+l_i_inc(id), SIGN(1,lt_i2(id)+l_i_inc(id) - last_sides(id))
          do j = 1, IT(1,i), 1
             loadc(IT(2,i)) = loadc(IT(2,i)) - LT(j,i)*load(IT(j+3,i))
          end do
          load(IT(3,i)) = loadc(IT(2,i))
          !$omp critical (ITT)
          !ITT( i ) = ITT( i ) + 1
          !$omp end critical (ITT)
       end do


       !$omp end parallel


       !$omp parallel private(k,i,j) num_threads(THREADS_FRONT)
       !$omp do
       do k = iDM, 1, -1
          !if (127899.le.lt_i(k+1).and.127899.ge.lt_i(k)+1) then
          !   write(*,*) lt_i(k+1), lt_i(k)+1, k
          !end if
          do i = lt_i(k+1), lt_i(k)+1, -1
             do j = 1, IT(1,i), 1
                loadc(IT(2,i)) = loadc(IT(2,i)) - LT(j,i)*load(IT(j+3,i))
             end do
             load(IT(3,i)) = loadc(IT(2,i))
             !$omp critical (ITT)
             !ITT( i ) = ITT( i ) + 1
             !$omp end critical (ITT)
          end do
       end do
       !$omp end do
       !$omp end parallel

      
       
    end if
    ! k = 0
    ! do i = 1, NVar, 1
    !    if(ITT(i).ne.1) then
    !       write(*,*) 'Error for var:', i, 'with', ITT(i)
    !       k = k + 1
    !    end if
    ! end do
    ! write(*,*) k, 'of', NVar, 'errors'

  end subroutine back_sub

  !Adds two partial fronts together
  subroutine add_front(front,front_m,loadf,loadf_m,loadf_dum,loadf_dumm,iPg,iPg_m,jPg,jPg_m,fs,fs_m,id2,id1,l_i,l_i_inc,th)
    implicit none
    real(kind=rk) :: front(fs_max,fs_max), front_m(fs_max,fs_max), loadf(fs_max), loadf_m(fs_max),&
         loadf_dum(fs_max), loadf_dumm(fs_max)
    integer(kind=ik) :: iPg(fs_max), jPg(fs_max), fs, l_i_inc,&
         iPg_m(fs_max), jPg_m(fs_max), fs_m, i, j, k, flag1, id2, id1, mode, th, rhead(fs), chead(fs),&
         l_i, rvs, cvs, RSUM(fs_max), CSUM(fs_max)
    !Should never be called for single domain

    if(id1.eq.id2) return

    !Modes determine which NDF to use
    !Combine domain id2 and id1
    
    NODf(:,id1) = NODf(:,id2) + NODf(:,id1)
    
    do i = 1, fs, 1
       !Find matching row
       flag1 = 0
       do k = 1, fs_m, 1
          if (iPg(i).eq.iPg_m(k)) then
             rhead(i) = k
             flag1 = 1

          end if
       end do
       if (flag1.eq.0) then
          fs_m = fs_m + 1      
          iPg_m(fs_m) = iPg(i)
          jPg_m(fs_m) = iPg(i)
          rhead(i) = fs_m
       end if
       
       flag1 = 0
       do k = 1, fs_m, 1
          if (jPg(i).eq.jPg_m(k)) then
             chead(i) = k
             flag1 = 1
          end if
       end do
       if (flag1.eq.0) then
          fs_m = fs_m + 1
          !rhead = fs_m
          iPg_m(fs_m) = jPg(i)
          jPg_m(fs_m) = jPg(i)
          chead(i) = fs_m
       end if    
    end do

    do i = 1, fs, 1
       do j = 1, fs, 1
          front_m(chead(j),rhead(i)) = front(j,i) + front_m(chead(j),rhead(i))
       end do
       loadf_m(rhead(i)) = loadf_m(rhead(i)) + loadf(i)
       loadf_dumm(rhead(i)) = loadf_dumm(rhead(i)) + loadf_dum(i)     
    end do
    
    
    rvs = 0
    cvs = 0
    do i =1, fs_m, 1
       if( NODf(iPg_m(i),id1).eq.NOA(iPg_m(i))) then
          !iPg_m(i) = -iPg_m(i)
          rvs = rvs + 1
          RSUM(rvs) = i
       end if
       if( NODf(jPg_m(i),id1).eq.NOA(jPg_m(i))) then
          !jPg_m(i) = -jPg_m(i)
          cvs = cvs + 1
          CSUM(cvs) = i
       end if
    end do
        !elim_summed3(front  ,loadf  ,loadf_dum ,iPg  ,jPg  ,rvs,cvs,RSUM,CSUM,l_i,l_i_inc,fs  ,id ,ele
    call elim_summed3(front_m,loadf_m,loadf_dumm,iPg_m,jPg_m,rvs,cvs,RSUM,CSUM,l_i,l_i_inc,fs_m,id1,id2)

    ! do i =1, fs_m, 1
    !    if( NODf(iPg_m(i),id1).eq.NOA(iPg_m(i))) iPg_m(i) = -iPg_m(i)
    !    if( NODf(jPg_m(i),id1).eq.NOA(jPg_m(i))) jPg_m(i) = -jPg_m(i)
    ! end do
    !call elim_stitch(front_m,loadf_m,loadf_dumm,iPg_m,jPg_m,l_i,fs_m,id1)
  end subroutine add_front


  !Solves front from element list
  Subroutine Front_solve(fs,front,loadf,loadf_dum,iPg,jPg,l_i,id,dm,l_i_inc)
    implicit none
    integer(kind=ik) :: ele1, ele2, fs, NOPf(NE,bas), iPg(fs_max),jPg(fs_max),e1, e2, es
    real(kind=rk) :: front(fs_max,fs_max), loadf(fs_max), loadf_dum(fs_max)
    integer(kind=ik) :: ele, i, j, k, rhead(NB), chead(NB), NOPPl(bas), flag1,&
         RSUM(fs_max), CSUM(fs_max), rvs, cvs, CPIV, RPIV, l_i, elim, j_max, k_max, &
         last_e, last_n, o, p, l, NK(NB), id, dm, l_i_inc, e_n(4)
    real(kind=rk) ::  local(NB,NB), loc(NB), PIVOT, FAC, t2, t1

    !Initialize local to global
    CPIV = 0 !Pivotal column
    RPIV = 0 !Pivotal row
    cvs = 0 !#summed col
    rvs = 0 !#summed rows (same as cvs)
    !piv = 0 !#small pivots (global now)

    !Find last appearance of node using NODf
    NOPf = NOP
    flag1 = 0


    !if (iDM.ne.dm) then



    e1 = 0
    do
       e1 = e1 + 1
       ele = ele_list(e1,dm)
       if (ele.eq.0) exit

       do i = 1, 9, 1
          o = NOP(ele,i)
          p = NOPP(o)
          NODf(p,dm) = NODf(p,dm) + 1
          if (NODf(p,dm).eq.NOA(p)) then
             NOPf(ele,i) = -NOP(ele,i)

          end if
          do k = 1, MDF(o)-1, 1
             NODf(p+k,dm) = NODf(p+k,dm) + 1
          end do
       end do
    end do


    e1 = 0
    do
       !Extract next element from proper list
       e1 = e1 + 1
       ele = ele_list(e1,dm)
       if (ele.eq.0) exit


       !USER_SPECIFIED
       !Assemble local using proper thread (id)
       if(.not.ASSOCIATED(assembler)) then
          write(*,*) 'Error in multifront: assembler not set'
          stop
       end if
       
       call assembler(ele,local,loc,NOPPl,NB,id,0)
       !USER_SPECIFIED

       !Create NK vector
       do i = 1, bas, 1
          j = NOPP(NOP(ele,i))
          o = NOPPl(i)
          do k = 0, MDF(NOP(ele,i))-1, 1
             NK(k+o) = k+j
             !write(*,*) NK(k+o)
          end do
       end do


       !Create heading vectors
       do o = 1, bas, 1
          do k = 0, MDF(NOP(ele,o))-1, 1
             i = NOPPl(o)+k
             flag1 = 0
             do j = 1, fs, 1
                if (NK(i).eq.iPg(j)) then
                   rhead(i) = j
                   flag1 = 1
                   exit
                end if
             end do
             !Create new row and column
             if (flag1.eq.0) then
                fs = fs + 1
                iPg(fs) = NK(i)
                jPg(fs) = NK(i)
                rhead(i) = fs
                chead(i) = fs
             else !find column
                do j = 1, fs, 1
                   if (NK(i).eq.jPg(j)) then                
                      chead(i) = j
                      exit
                   end if

                end do
             end if

          end do
       end do

       !   !Create heading vectors
       ! do i = 1, NCN(ele), 1
       !    flag1 = 0
       !    do j = 1, fs, 1
       !       if (NK(i).eq.iPg(j)) then
       !          rhead(i) = j
       !          flag1 = 1
       !          exit
       !       end if
       !    end do
       !    !Create new row and column
       !    if (flag1.eq.0) then
       !       fs = fs + 1
       !       iPg(fs) = NK(i)
       !       jPg(fs) = NK(i)
       !       rhead(i) = fs
       !       chead(i) = fs
       !    else !find column
       !       do j = 1, fs, 1
       !          if (NK(i).eq.jPg(j)) then                
       !             chead(i) = j
       !             exit
       !          end if

       !       end do
       !    end if

       ! end do
       !IF (S_MODE.EQ.0) write(*,*) ELE, FS

       !Warning front width exceeded, will probably seg fault
       if (fs.gt.fs_max) then
          write(*,*) 'Error front exceeds max size, fs:', fs, fs_max
          t1 = 0.0_rk
          loadf(1) = 1.0_rk/t1
          pause
          return
       end if

       ! !Add local to front
       ! do  o = 1, bas, 1
       !    do k = 0, MDF(NOP(ele,o))-1, 1
       !       i = NOPPl(o)+k

       !       do  p = 1, bas, 1
       !          do l = 0, MDF(NOP(ele,p))-1, 1
       !             j = NOPPl(p) + l
       !             front(chead(j),rhead(i)) = local(j,i) + front(chead(j),rhead(i))
       !          end do

       !       end do
       !       loadf(rhead(i)) = loc(i) + loadf(rhead(i))
       !       loadf_dum(rhead(i)) = loc(i) + loadf_dum(rhead(i))
       !    end do
       ! end do
       ! !   do i = 1, NCN(ele), 1

       ! !    do j = 1, NCN(ele), 1
       ! !       front(chead(j),rhead(i)) = local(j,i) + front(chead(j),rhead(i))
       ! !    end do
       ! !    loadf(rhead(i)) = loc(i) + loadf(rhead(i))
       ! !    loadf_dum(rhead(i)) = loc(i) + loadf_dum(rhead(i))
       ! ! end do


       !Add local to front
       if(SWAP_LOCAL_IJ) then
          do  o = 1, bas, 1
             do k = 0, MDF(NOP(ele,o))-1, 1
                i = NOPPl(o)+k

                do  p = 1, bas, 1
                   do l = 0, MDF(NOP(ele,p))-1, 1
                      j = NOPPl(p) + l
                      front(chead(j),rhead(i)) = local(i,j) + front(chead(j),rhead(i))
                   end do

                end do
                loadf(rhead(i)) = loc(i) + loadf(rhead(i))
                loadf_dum(rhead(i)) = loc(i) + loadf_dum(rhead(i))
             end do
          end do
       else
          do  o = 1, bas, 1
             do k = 0, MDF(NOP(ele,o))-1, 1
                i = NOPPl(o)+k

                do  p = 1, bas, 1
                   do l = 0, MDF(NOP(ele,p))-1, 1
                      j = NOPPl(p) + l
                      front(chead(j),rhead(i)) = local(j,i) + front(chead(j),rhead(i))
                   end do

                end do
                loadf(rhead(i)) = loc(i) + loadf(rhead(i))
                loadf_dum(rhead(i)) = loc(i) + loadf_dum(rhead(i))
             end do
          end do
       end if


       !Find fully summed using NOPf
       do i = 1, bas, 1
          if (NOPf(ele,i).lt.0) then
             do k = 0, MDF(NOP(ele,i))-1, 1

                rvs = rvs + 1
                RSUM(rvs) = rhead(k+NOPPl(i))
                cvs = cvs + 1
                CSUM(cvs) = chead(k+NOPPl(i))

             end do
          end if
       end do

       !Call elimination phase
       call elim_summed3(front,loadf,loadf_dum,iPg,jPg,rvs,cvs,RSUM,CSUM,l_i,l_i_inc,fs,id,ele)



    end do

    !if (piv.ne.0) write(*,*) piv, 'small pivots'
    !write(*,*) 'rem. fs:', fs

  End Subroutine Front_solve

 subroutine elim_summed3(front,loadf,loadf_dum,iPg,jPg,rvs,cvs,RSUM,CSUM,l_i,l_i_inc,fs,id,ele)
    implicit none
    integer(kind=ik) :: fs
    integer(kind=ik) :: iPg(fs_max),jPg(fs_max)
    real(kind=rk) :: front(fs_max,fs_max)
    real(kind=rk) :: loadf(fs_max), loadf_dum(fs_max)
    integer(kind=ik) ::  i, k, j, id, l_i_inc,&
         RSUM(fs_max), CSUM(fs_max), rvs, cvs, CPIV, RPIV, l_i, elim, j_max, k_max, o, ele
    real(kind=rk) ::  PIVOT, c_piv, GG(fs_max)
    !write(*,*) 'Called'
    !eliminate fully summed full pivotal choice
    elim = rvs
    if (l_i+l_i_inc.lt.1.or.l_i+l_i_inc.gt.NVar) write(*,*) 'Error l_i', l_i, NVar, l_i_inc
    do i = 1, elim, 1

       !Find sufficiently large pivot
       j_max = 1
       k_max = 1
       PIVOT = ABS(front(CSUM(k_max),RSUM(j_max)))
       if (PIVOT.lt.0.0001_rk) then
          do k = 1, rvs, 1

             do j = 1, cvs, 1
                if (PIVOT.lt.ABS(front(CSUM(k),RSUM(j)))) then
                   j_max = j
                   k_max = k
                   PIVOT = ABS(front(CSUM(k),RSUM(j)))
                end if
             end do
          end do
       end if
       
       RPIV = RSUM(j_max)
       RSUM(j_max) = RSUM(rvs)
       do j = 1, rvs
          if(RSUM(j).eq.fs) then         
             RSUM(j) = RPIV
          end if
       end do
       
       CPIV = CSUM(k_max)
       CSUM(k_max) = CSUM(cvs)
       do j = 1, cvs
          if(CSUM(j).eq.fs) then         
             CSUM(j) = CPIV
          end if
       end do

       !Decrement rows and columns to delete
       rvs = rvs - 1
       cvs = cvs - 1

       !Set PIVOT
       PIVOT = front(CPIV,RPIV)
       !Check if pivot too small
       if (ABS(PIVOT).le.1.0e-18_rk) then
          if(check_pivot.eq.1) write(*,*) 'Small Pivot:', PIVOT, id, ele, rvs
           if(IT(1,l_i+l_i_inc).ne.0) write(*,*) 'o is not 0', IT(1, l_i+l_i_inc), id, ele
          !call find_var_info(iPG(RPIV),ele,id)
          piv(id) = piv(id)+1
       end if

       !Normalize Pivotal Row
       do j = 1, fs, 1
          front(j,RPIV) = front(j,RPIV)/PIVOT
       end do
       loadf(RPIV) = loadf(RPIV)/PIVOT
       
       o = fs - 1
       
        !Cogy out pivotal row and column
       l_i = l_i + l_i_inc
       do j = 1, o, 1
          LT(j,l_i) = front(j,RPIV)
          GG(j) = front(CPIV,j)                                
       end do
       GG(RPIV) = front(CPIV,fs)               
       LT(CPIV,l_i) = front(fs,RPIV)     
       c_piv = loadf(RPIV) !Save rhs
       
       !Swap row and col
       do j = 1, fs, 1
          front(j,RPIV) = front(j,fs)
          front(j,fs) = 0.0_rk
       end do

       do j = 1, fs, 1
          front(CPIV,j) = front(fs,j)
          front(fs,j) = 0.0_rk
       end do

       PIVOT = loadf(RPIV)
       loadf(RPIV) = loadf(fs)
       loadf(FS) = PIVOT

       PIVOT = loadf_dum(RPIV)
       loadf_dum(RPIV) = loadf_dum(fs)
       loadf_dum(FS) = PIVOT

       j = iPg(RPIV)
       iPg(RPIV) = iPg(fs)
       iPg(fs) = j

       j = jPg(CPIV)
       jPg(CPIV) = jPg(fs)
       jPg(fs) = j
       
       !Save pivotal row
                             !increment eliminated var count
           !(fs-1)!Save number in row
       if(IT(1,l_i).ne.0) then
          write(*,*) 'o is not 0', IT(1, l_i), id, ele
          !pause
       end if
       IT(1,l_i) = o  
       IT(2,l_i) = iPg(fs)!(RPIV) !Save global headers
       IT(3,l_i) = jPg(fs)!(CPIV) !Save column header for pivot
       loadc(iPg(fs)) = loadf(fs)!RPIV,RPIV 
       load_dum(iPg(fs)) = loadf_dum(fs)

       do j = 1, o, 1        
          IT(3+j,l_i) = jPg(j)
       end do
     
       !Perform actual elimination
       ![                            |           ]
       ![                            |           ]
       ![              1             |     2     ]
       ![                            |           ]
       ![                            |           ]
       ![                            |           ]
       ![                            |           ]
       ![----------------------------------------]
       ![                            |           ]
       ![              3             |     4     ]
       ![                            |           ]

       ! if (ths.eq.1) then
       do k = 1, o, 1
          
          if (GG(k).ne.0.0_rk) then
             do j = 1, o, 1 !Doing Sector 1 update
                front(j,k) = front(j,k) - LT(j,l_i)*GG(k)
             end do
             loadf(k) = loadf(k) - c_piv*GG(k)
          end if

       end do

       loadf(fs) = 0.0_rk
       loadf_dum(fs) = 0.0_rk

       fs  = fs - 1

       if(debug_NAN) then
          do k = 1, fs, 1
             do j = 1, fs, 1
                if (front(k,j).ne.front(k,j)) then
                   write(*,*) 'NaN in front:', id, l_i, iPG(RPIV)
                   if(ASSOCIATED(var_finder)) call var_finder(iPG(RPIV),ele,id)                  
                end if
             end do
          end do
       end if
          

    end do

  end subroutine elim_summed3

  subroutine nullifier()
    implicit none
    NULLIFY(NVar)
    NULLIFY(NN)
    NULLIFY(NE)
    NULLIFY(s_mode)
    NULLIFY(bas)
    NULLIFY(MDF)
    NULLIFY(NOPP)
    NULLIFY(DNOP)
    NULLIFY(NOP)
    NULLIFY(rNOP)
    NULLIFY(load)
  end subroutine nullifier

  subroutine check_associated()
    implicit none

    if(.not.ASSOCIATED(NVar)) then
       write(*,*) 'Error in multifront: NVar not set'
       stop
    end if

    if(.not.ASSOCIATED(NN)) then
       write(*,*) 'Error in multifront: NN not set'
       stop
    end if

    if(.not.ASSOCIATED(NE)) then
       write(*,*) 'Error in multifront: NE not set'
       stop
    end if

    if(.not.ASSOCIATED(s_mode)) then
       write(*,*) 'Error in multifront: s_mode not set'
       stop
    end if

    if(.not.ASSOCIATED(bas)) then
       write(*,*) 'Error in multifront: bas not set'
       stop
    end if

    if(.not.ASSOCIATED(MDF)) then
       write(*,*) 'Error in multifront: MDF not set'
       stop
    end if

    if(.not.ASSOCIATED(NOPP)) then
       write(*,*) 'Error in multifront: NOPP not set'
       stop
    end if

    if(.not.ASSOCIATED(DNOP)) then
       write(*,*) 'Error in multifront: DNOP not set'
       stop
    end if

    if(.not.ASSOCIATED(NOP)) then
       write(*,*) 'Error in multifront: NOP not set'
       stop
    end if

    if(.not.ASSOCIATED(rNOP)) then
       write(*,*) 'Error in multifront: rNOP not set'
       stop
    end if

    if(.not.ASSOCIATED(load)) then
       write(*,*) 'Error in multifront: load not set'
       stop
    end if

    if(.not.ASSOCIATED(excluder)) then
       write(*,*) 'By default including all elements in multifront'
       excluder => exclude_default
    end if


  end subroutine check_associated

end module front_mod
