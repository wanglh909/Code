module data
  use kind
  implicit none
  !?split into data & Ldatas
  
  integer(kind=ik), target :: NVar, NTN, NTE, s_mode, bas = 9
  !s_mode = 0: solve for dynamics, 1: mesh only
  !s_mode: Set to 0 if solving full dynamics and 1 if only mesh (must be 0 if solve mesh by putting BC's on all velocity and pressure vars)
  !bas: number of nodes in element (9)
  integer(kind=ik), allocatable, target :: rNOP(:,:,:), globalNM(:,:), MDF(:), NOPP(:), RegN(:)
  real(kind=rk), allocatable, target :: dsol(:)
  integer(kind=ik) :: SOLVER_MODE

  
  real(kind=rk):: R, outer, substrate
  real(kind=rk):: Re, Ca, Oh, Grav, Kdi, KBCgroup, Pe, REH, beta, F0, kR, MaN, Da_sub, Da_surf1, Da_surf2 !Marangoni number (MaN)

  integer(kind=ik):: NEL, NEM, NEV, NES, NNX34, NNX3456, NNXV, NNR1278, NNR1278p, iBW,iBW1,iBW2
  !NVar: Number of Variables   !NNX: number of nodes in r (abscissa)   
  integer(kind=ik):: NND, NNV, NNS, ED, EV, ES
  !NND: number of nodes of droplet,   NNV: number of nodes of vapor
  !ED: number of elements of droplet, EV: number of elements of vapor
  integer(kind=ik):: NNVar1 , NNVar2 = 3, NNVar3 = 3   !NNVar: node number of variable (r,z,u,v,p --> 5)
  !NNVar1: r,z,T,u,v,cp,p; define in program main(=7 if solve for T, =6 if not)
  !NNVar2: r,z,c; Nnvar3: r,z,T

  !integer(kind=ik), parameter:: Nr = 0, Nz = 1, Nu = 3, Nv = 4, NT = 2, Ncp = 5, Np = 6, NTs = 2 !Nc = MDF(i)-1
  integer(kind=ik):: Nr = 0, Nz = 1, NT = 2, Nu, Nv, Ncp, Np, Ngamma, NTs = 2 !Nc = MDF(i)-1

  integer(kind=ik):: timestep, step

  real(kind=rk):: error1, error2

  integer(kind=ik), allocatable:: rowNM(:), columnNM(:), WFLAG(:)!, RegNN(:)
  !RegN: tell which region the element is in
  !RegNN: tell which region the node is in
  integer(kind=ik),allocatable:: Ngrid(:,:)
  integer(kind=ik), allocatable:: MDFd(:), PN(:)
  !PN = 1, node with pressure; 0, node without pressure


  real(kind=rk):: time
  real(kind=rk), allocatable:: rcoordinate(:), zcoordinate(:), usol(:), vsol(:), Tsol(:), psol(:), csol(:), cpsol(:), gammasol(:)
  !solp & dtp means previous solution & dt, solpred means the predicted solution
  real(kind=rk), allocatable:: sol(:), solp(:), soldot(:), soldotp(:), soldotpp(:), solpred(:)

  real(kind=rk), allocatable:: fsi_size(:), geta_size(:)

  ! real(kind=rk), parameter, dimension(3):: gausspoint = (/0.5_rk*(-0.7745966692414833770358531_rk + 1.0_rk),&
  !      0.5_rk*(0.000000000000000000000000000000_rk + 1.0_rk),0.5_rk*(0.7745966692414833770358531_rk + 1.0_rk)/)



  !for subroutine basis_function
  real(kind=rk):: phi(3,3,9), phisi(3,3,9), phieta(3,3,9), psi(3,3,9), psisi(3,3,9), psieta(3,3,9), &
       phi_1d(3,3), phix_1d(3,3), phixx_1d(3), phisi0_1d(3,9), phisi1_1d(3,9), phieta0_1d(3,9), phieta1_1d(3,9), &
       phi0_1d(3,9), phi1_1d(3,9), phix_1d_flux(3)
  !phisi_1d: 3 is for three gausspoints of eta, si = 1.0, 9 is for the number of basis functions

  integer(kind=ik),parameter:: convert49(4) = (/1, 3, 7, 9/)


  integer, parameter:: Ng = 3 !number of gausspoints


  !for subroutine prediction
  real(kind=rk):: dt, dtp, CTJ, change, trunerr
  real(kind=rk), parameter:: eps = 1.0e-3_rk

  
  !for size_function change
  integer(kind=ik):: size_function_change, final_size

  !for spherical cap 
  real(kind=rk):: angle_c, angle_cd, angle_c_degree, mtwoH
  real(kind=rk),parameter:: pi=3.141592653589793238462643383279_rk 

  !fixed timestep change
  integer(kind=ik):: FTS  !fixed time steps

  !vapor phase
  integer(kind=ik), allocatable:: VE(:), VN(:)  !vapor element, vapor node
  !VE = 0 for drop element, 1 for vapor element, 5 for substrate element
  !VN = 0 for drop, 1 for vapor, 2 for free surface, 5 for substrate( exclude base nodes(judge by BCflagN) )
  real(kind=rk):: Hum         !Hum: humidity
  integer(kind=ik),allocatable:: NOPDV(:,:), layer(:)   
  !NOPDV(:,1): globalNM in droplet, NOPDV(:,2): globalNM in vapor
  integer(kind=ik),allocatable:: BCflagE(:,:), BCflagN(:,:) 
  !BCflag(:,1) is for axis, (:,2) for base, (:,3) for free surface, (:,4) for outer circle, (:,4) for substrate base; 
  !BCflag(1,2) = 0 means node1 is not on base

  integer(kind=ik):: mesh_decouple = 1, base_mesh = 1

  !vapor_mesh_size 
  integer(kind=ik):: vlayer
  integer(kind=ik), allocatable:: NEV1(:)
  real(kind=rk), allocatable:: R11(:)


  !frontal solver debug
  real(kind=rk), allocatable:: Jac(:,:), Res(:)
  integer(kind=ik):: check_0_in_Jac 

  integer(kind=ik):: graph_mode !1: graph each step; 0: graph each timestep

  !temperature convection control
  integer(kind=ik):: Ttime, Tconv, Tdiff, TtimeS, TdiffS, cptime, cpconv, cpdiff
  integer(kind=ik):: NStrans, Inert, Capil, Viscous, GravI

  real(kind=rk):: EvapSpeed, EvapSpeedp=0.0_rk, VolEvap1=0.0_rk, VolEvap2=0.0_rk, VolEvap3=0.0_rk
  integer(kind=ik):: initial_vapor_solved, initial_vapor_solving
  !indicator for calculating and graphing initial vapor adjustment

  integer(kind=ik):: diverge

  !algebraic mesh for contact corner
  real(kind=rk):: k_alge, x_alge  !x_alge: coordinate of algebraic mesh start point
  integer(kind=ik):: NEM_alge
  integer(kind=ik), allocatable:: algeN(:), algeS(:)

  !substrate phase
  ! integer(kind=ik):: SN

  real(kind=rk):: rmax, zmax, umax, vmax, Tmax, cpmax, pmax, cmax
  integer(kind=ik):: graph_step   !graph every how many timestep

  character(LEN=40):: folder
  ! logical:: exist

  integer:: uniflux    !determine if the imposed flux is uniform. Notice: If flux is uniform, still use divergent heat flux.
  
  integer(kind=ik):: simple_mesh

  integer(kind=ik):: alge_corner, no_vapor

  real(kind=rk):: ztop

  integer(kind=ik):: top_node, top_element, CL_element, CL_node, angle_c_node
  integer(kind=ik):: read_coordinate_value

  real(kind=rk):: T_sub


  real(kind=rk):: diameterp, Dp, Pep, Pep_surf, kboltz = 1.3806485279e-23 !(J/K)    


  integer(kind=ik):: Maran_flow, solve_T, fixed_Ma, solve_cp, surf_adsp, sub_adsp, muchange, fixed_dt

  integer(kind=ik):: init_stability, timestep_stable, radial_cal_time
  
  !maximum packing
  integer(kind=ik), allocatable:: pack_flag(:)
  real(kind=rk):: cp_pack, pack_condition,cp_average  !, volume0 


  real(kind=rk):: angle_int

end module data
