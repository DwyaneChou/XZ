module spatial_operators_mod
  use constants_mod
  use mesh_mod
  use parameters_mod
  use stat_mod
  use tend_mod
  use reconstruction_mod
  use test_case_mod
  implicit none
  
  private
  
  public init_spatial_operator,spatial_operator
  
  logical, dimension(:,:), allocatable :: inDomain
  
  integer(i_kind), dimension(:,:), allocatable :: iCenCell ! center cell index on reconstruction stencil
  
  integer(i_kind), dimension(:,:,:), allocatable :: iRecCell ! x index of reconstruction cells
  integer(i_kind), dimension(:,:,:), allocatable :: kRecCell ! k index of reconstruction cells
  
  real   (r_kind), dimension(:,:,:,:), allocatable :: xRel   ! relative x coordinate of reconstruction cells
  real   (r_kind), dimension(:,:,:,:), allocatable :: etaRel ! relative eta coordinate of reconstruction cells
  
  real   (r_kind), dimension(:,:,:,:), allocatable :: recMatrixL
  real   (r_kind), dimension(:,:,:,:), allocatable :: recMatrixR
  real   (r_kind), dimension(:,:,:,:), allocatable :: recMatrixB
  real   (r_kind), dimension(:,:,:,:), allocatable :: recMatrixT
  real   (r_kind), dimension(:,:,:,:), allocatable :: recMatrixG ! Reconstruction matrix on Gaussian points
  
  real   (r_kind) :: dV
  
  real   (r_kind), dimension(:,:,:  ), allocatable :: qC
  real   (r_kind), dimension(:,:,:,:), allocatable :: qL
  real   (r_kind), dimension(:,:,:,:), allocatable :: qR
  real   (r_kind), dimension(:,:,:,:), allocatable :: qB
  real   (r_kind), dimension(:,:,:,:), allocatable :: qT
  real   (r_kind), dimension(:,:,:,:), allocatable :: qG

  real(r_kind), dimension(:,:,:), allocatable :: F
  real(r_kind), dimension(:,:,:), allocatable :: H
  real(r_kind), dimension(  :,:), allocatable :: P
  
  real(r_kind), dimension(:,:,:,:), allocatable :: FL    ! Reconstructed F_(i-1/2,k)
  real(r_kind), dimension(:,:,:,:), allocatable :: FR    ! Reconstructed F_(i+1/2,k)
  real(r_kind), dimension(:,:,:,:), allocatable :: FB    ! Reconstructed F_(i,k-1/2)
  real(r_kind), dimension(:,:,:,:), allocatable :: FT    ! Reconstructed F_(i,k+1/2)
      
  real(r_kind), dimension(:,:,:,:), allocatable :: HL    ! Reconstructed H_(i-1/2,k)
  real(r_kind), dimension(:,:,:,:), allocatable :: HR    ! Reconstructed H_(i+1/2,k)
  real(r_kind), dimension(:,:,:,:), allocatable :: HB    ! Reconstructed H_(i,k-1/2)
  real(r_kind), dimension(:,:,:,:), allocatable :: HT    ! Reconstructed H_(i,k+1/2)
  
  real(r_kind), dimension(  :,:), allocatable :: PC
  real(r_kind), dimension(:,:,:), allocatable :: PL    ! Reconstructed P_(i-1/2,k)
  real(r_kind), dimension(:,:,:), allocatable :: PR    ! Reconstructed P_(i+1/2,k)
  real(r_kind), dimension(:,:,:), allocatable :: PB    ! Reconstructed P_(i,k-1/2)
  real(r_kind), dimension(:,:,:), allocatable :: PT    ! Reconstructed P_(i,k+1/2)
  
  real(r_kind), dimension(:,:,:  ), allocatable :: Fe    ! F on edges of each cell
  real(r_kind), dimension(:,:,:  ), allocatable :: He    ! H on edges of each cell
  
  real(r_kind), dimension(:,:,:,:), allocatable :: FeP   ! F on points on edges of each cell
  real(r_kind), dimension(:,:,:,:), allocatable :: HeP   ! H on points on edges of each cell
  
  real(r_kind), dimension(:,:,:), allocatable :: src   ! source term
  
  real(r_kind), dimension(  :,:), allocatable :: rho_p ! density perturbation
  
  real(r_kind), dimension(:,:,:,:), allocatable :: qL_ref
  real(r_kind), dimension(:,:,:,:), allocatable :: qR_ref
  real(r_kind), dimension(:,:,:,:), allocatable :: qB_ref
  real(r_kind), dimension(:,:,:,:), allocatable :: qT_ref
  real(r_kind), dimension(:,:,:,:), allocatable :: qG_ref
  
  real(r_kind), dimension(  :,:), allocatable :: PC_ref
  real(r_kind), dimension(:,:,:), allocatable :: PL_ref! Reconstructed P_ref_(i-1/2,k)
  real(r_kind), dimension(:,:,:), allocatable :: PR_ref! Reconstructed P_ref_(i+1/2,k)
  real(r_kind), dimension(:,:,:), allocatable :: PB_ref! Reconstructed P_ref_(i,k-1/2)
  real(r_kind), dimension(:,:,:), allocatable :: PT_ref! Reconstructed P_ref_(i,k+1/2)
  
  real(r_kind), dimension(:,:), allocatable :: relax_coef ! Relax coefficient of Rayleigh damping
  
  real(r_kind), dimension(:,:,:), allocatable :: q_diff ! u wind, for viscosity terms only
  
  real(r_kind), dimension(:,:,:,:), allocatable :: CG_coef
  
  real(r_kind), dimension(:,:), allocatable :: eigen_mtx_x
  real(r_kind), dimension(:,:), allocatable :: eigen_mtx_z
  
contains
  subroutine init_spatial_operator
    integer(i_kind) :: i,j,k,iR,kR,iVar,iPOE,iEOC
    integer(i_kind) :: iRec,kRec
    
    real(r_kind) :: uB_ref,uT_ref
    real(r_kind) :: nx,nz,nv(2),pm(2,2)
    real(r_kind) :: xG  (nQuadPointsOnCell)
    real(r_kind) :: etaG(nQuadPointsOnCell)
    
    allocate(inDomain  (ics:ice,kcs:kce))
    
    allocate(iCenCell  (ids:ide,kds:kde))
    
    allocate(iRecCell  (maxRecCells,ids:ide,kds:kde))
    allocate(kRecCell  (maxRecCells,ids:ide,kds:kde))
    
    allocate(xRel  (4,maxRecCells,ids:ide,kds:kde))
    allocate(etaRel(4,maxRecCells,ids:ide,kds:kde))
    
    allocate(recMatrixL(nPointsOnEdge,maxRecTerms,ids:ide,kds:kde))
    allocate(recMatrixR(nPointsOnEdge,maxRecTerms,ids:ide,kds:kde))
    allocate(recMatrixB(nPointsOnEdge,maxRecTerms,ids:ide,kds:kde))
    allocate(recMatrixT(nPointsOnEdge,maxRecTerms,ids:ide,kds:kde))
    
    allocate(recMatrixG(nQuadPointsOnCell,maxRecTerms,ids:ide,kds:kde))
    
    allocate(qC(nVar,              ics:ice,kcs:kce))
    allocate(qL(nVar,nPointsOnEdge,ics:ice,kcs:kce))
    allocate(qR(nVar,nPointsOnEdge,ics:ice,kcs:kce))
    allocate(qB(nVar,nPointsOnEdge,ics:ice,kcs:kce))
    allocate(qT(nVar,nPointsOnEdge,ics:ice,kcs:kce))
    
    allocate(qG(nVar,nQuadPointsOnCell,ics:ice,kcs:kce))
  
    allocate(F (nVar,ids:ide,kds:kde))
    allocate(H (nVar,ids:ide,kds:kde))
    allocate(P (     ids:ide,kds:kde))
    
    allocate(FL(nVar,nPointsOnEdge,ics:ice,kcs:kce))
    allocate(FR(nVar,nPointsOnEdge,ics:ice,kcs:kce))
    allocate(HB(nVar,nPointsOnEdge,ics:ice,kcs:kce))
    allocate(HT(nVar,nPointsOnEdge,ics:ice,kcs:kce))
    
    allocate(PC(                   ics:ice,kcs:kce))
    allocate(PL(     nPointsOnEdge,ics:ice,kcs:kce))
    allocate(PR(     nPointsOnEdge,ics:ice,kcs:kce))
    allocate(PB(     nPointsOnEdge,ics:ice,kcs:kce))
    allocate(PT(     nPointsOnEdge,ics:ice,kcs:kce))
    
    allocate(Fe(nVar,ids:ide+1,kds:kde  ))
    allocate(He(nVar,ids:ide  ,kds:kde+1))
    
    allocate(FeP(nVar,nPointsOnEdge,ids:ide+1,kds:kde  ))
    allocate(HeP(nVar,nPointsOnEdge,ids:ide  ,kds:kde+1))
    
    allocate(src  (nVar,ids:ide,kds:kde))
    
    allocate(rho_p (    ids:ide,kds:kde))
  
    allocate(qL_ref  (nVar,nPointsOnEdge,ics:ice,kcs:kce))
    allocate(qR_ref  (nVar,nPointsOnEdge,ics:ice,kcs:kce))
    allocate(qB_ref  (nVar,nPointsOnEdge,ics:ice,kcs:kce))
    allocate(qT_ref  (nVar,nPointsOnEdge,ics:ice,kcs:kce))
    
    allocate(qG_ref  (nVar,nQuadPointsOnCell,ics:ice,kcs:kce))
    
    allocate(PC_ref(              ics:ice,kcs:kce))
    allocate(PL_ref(nPointsOnEdge,ics:ice,kcs:kce))
    allocate(PR_ref(nPointsOnEdge,ics:ice,kcs:kce))
    allocate(PB_ref(nPointsOnEdge,ics:ice,kcs:kce))
    allocate(PT_ref(nPointsOnEdge,ics:ice,kcs:kce))
    
    allocate(relax_coef(ics:ice,kcs:kce))
    
    allocate(q_diff(nVar,ics:ice,kcs:kce))
    
    allocate(CG_coef(nVar,maxRecTerms,ids:ide,kds:kde))
    
    !allocate(eigen_mtx_x(nVar,nVar))
    !allocate(eigen_mtx_z(nVar,nVar))
    allocate(eigen_mtx_x(4,4))
    allocate(eigen_mtx_z(4,4))
  
    src   = 0
    
    dV = dx * deta
    
    CG_coef = 1.
    
    ! set reconstruction cells on each stencil
    print*,'Set reconstruction cells on each stencil'
    print*,''
    inDomain(ics:ice,kcs:kce) = .false. 
    inDomain(ids:ide,kds:kde) = .true.
    ! Scheme 1
    do k = kds,kde
      do i = ids,ide
        j = 0
        do kRec = -recBdy,recBdy
          do iRec = -recBdy,recBdy
        !do kRec = -1,1
        !  do iRec = -1,1
            if(inDomain(i+iRec,k+kRec))then
              j = j + 1
              iRecCell(j,i,k) = i+iRec
              kRecCell(j,i,k) = k+kRec
            endif
          enddo
        enddo
        nRecCells(i,k) = j
        
        locPolyDegree(i,k) = min( maxval(iRecCell(1:j,i,k)) - minval(iRecCell(1:j,i,k)), maxval(kRecCell(1:j,i,k)) - minval(kRecCell(1:j,i,k)) )
        
        !if( nRecCells(i,k) >= 3           ) locPolyDegree(i,k) = 1
        !if( nRecCells(i,k) >= 9           ) locPolyDegree(i,k) = 2
        !if( nRecCells(i,k) >= 16          ) locPolyDegree(i,k) = 3
        !if( nRecCells(i,k) >= 25          ) locPolyDegree(i,k) = 4
        !if( nRecCells(i,k) == maxRecCells ) locPolyDegree(i,k) = recPolyDegree
      enddo
    enddo
    
    ! special treatment on low boundary
    do i = ibs,ibe
      k = kds
      j = 0
      do kRec = 0,2
        do iRec = -1,1
          j = j + 1
          iRecCell(j,i,k) = i + iRec
          kRecCell(j,i,k) = k + kRec
        enddo
      enddo
      nRecCells    (i,k) = j
      locPolyDegree(i,k) = 2
      
      k = kds + 1
      j = 0
      do kRec = -1,3
        do iRec = -recBdy,recBdy
          j = j + 1
          iRecCell(j,i,k) = i + iRec
          kRecCell(j,i,k) = k + kRec
        enddo
      enddo
      nRecCells    (i,k) = j
      locPolyDegree(i,k) = 3
      
      !do k = kds, kds + recBdy - 1
      !  j = 0
      !  do kRec = - k + 1, - k + stencil_width
      !    do iRec = -recBdy,recBdy
      !      j = j + 1
      !      iRecCell(j,i,k) = i + iRec
      !      kRecCell(j,i,k) = k + kRec
      !    enddo
      !  enddo
      !  nRecCells    (i,k) = j
      !  locPolyDegree(i,k) = 4!k + 1
      !enddo
      
      !locPolyDegree(i,k) = min( maxval(iRecCell(1:j,i,k)) - minval(iRecCell(1:j,i,k)), maxval(kRecCell(1:j,i,k)) - minval(kRecCell(1:j,i,k)) )
    enddo
    
    !!! Scheme 2
    !!do i = ibs,ibe
    !!  do k = kds,kds+1
    !do k = kds,kde
    !  do i = ids,ide
    !    j = 0
    !    do kRec = -recBdy,recBdy
    !      do iRec = -recBdy,recBdy
    !        if(abs(iRec)+abs(kRec)<=recBdy)then
    !          if(inDomain(i+iRec,k+kRec))then
    !            j = j + 1
    !            iRecCell(j,i,k) = i + iRec
    !            kRecCell(j,i,k) = k + kRec
    !          endif
    !        endif
    !      enddo
    !    enddo
    !    nRecCells(i,k) = j
    !    if(nRecCells(i,k)>6 )locPolyDegree(i,k) = 2
    !    if(nRecCells(i,k)>10)locPolyDegree(i,k) = 3
    !    if(nRecCells(i,k)>16)locPolyDegree(i,k) = 4
    !    if(nRecCells(i,k)>21)locPolyDegree(i,k) = 5
    !    !print*,k,i,j,locPolyDegree(i,k)
    !  enddo
    !enddo
    
    !! special treatment on low boundary
    !do k = kds,kds+1
    !  do i = ibs,ibe
    !    j = 0
    !    do kRec = 1,5 ! stencil_width=5
    !      do iRec = 1,5
    !        if(inDomain(i-2+iRec-1,k-2+kRec-1))then! recBdy=2
    !          j = j + 1
    !          iRecCell(j,i,k) = i-2+iRec-1
    !          kRecCell(j,i,k) = k-2+kRec-1
    !        endif
    !      enddo
    !    enddo
    !    nRecCells(i,k) = j
    !      
    !    if( nRecCells(i,k) >= 3           ) locPolyDegree(i,k) = 1
    !    if( nRecCells(i,k) >= 9           ) locPolyDegree(i,k) = 2
    !    if( nRecCells(i,k) >= 16          ) locPolyDegree(i,k) = 3
    !    if( nRecCells(i,k) >= 25          ) locPolyDegree(i,k) = 4
    !    if( nRecCells(i,k) == maxRecCells ) locPolyDegree(i,k) = recPolyDegree
    !  enddo
    !enddo

    do k = kds,kde
      do i = ids,ide
        nRecTerms(i,k) = ( locPolyDegree(i,k) + 1 ) * ( locPolyDegree(i,k) + 2 ) / 2
      enddo
    enddo
    
    ! initilize polyCoordCoef
    print*,'Initilize polyCoordCoef'
    print*,''
    
    recCoef = 1.
    recdx   = 1. / ( dx   * recCoef )
    recdeta = 1. / ( deta * recCoef )
    recdV   = 1. / ( recCoef**2 )
    
    !$OMP PARALLEL DO PRIVATE(i,j,iRec,kRec,xG,etaG)
    do k = kds,kde
      do i = ids,ide
        do j = 1,nRecCells(i,k)
          iRec = iRecCell(j,i,k)
          kRec = kRecCell(j,i,k)
          
          if(iRec==i.and.kRec==k)iCenCell(i,k) = j
          
          xRel  (:,j,i,k) = ( xCorner  (:,iRec,kRec) - xCenter  (i,k) ) * recdx
          etaRel(:,j,i,k) = ( etaCorner(:,iRec,kRec) - etaCenter(i,k) ) * recdeta
          
          call calc_polynomial_square_integration(locPolyDegree(i,k),xRel  (1,j,i,k),xRel  (2,j,i,k),&
                                                                     etaRel(1,j,i,k),etaRel(4,j,i,k),polyCoordCoef(j,:,i,k))
        enddo
    
        ! Calculate reconstruction matrix on edge
        call calc_polynomial_matrix(locPolyDegree(i,k),nPointsOnEdge,nRecTerms(i,k),xL*recdx,etaL*recdeta,recMatrixL(:,1:nRecTerms(i,k),i,k))
        call calc_polynomial_matrix(locPolyDegree(i,k),nPointsOnEdge,nRecTerms(i,k),xR*recdx,etaR*recdeta,recMatrixR(:,1:nRecTerms(i,k),i,k))
        call calc_polynomial_matrix(locPolyDegree(i,k),nPointsOnEdge,nRecTerms(i,k),xB*recdx,etaB*recdeta,recMatrixB(:,1:nRecTerms(i,k),i,k))
        call calc_polynomial_matrix(locPolyDegree(i,k),nPointsOnEdge,nRecTerms(i,k),xT*recdx,etaT*recdeta,recMatrixT(:,1:nRecTerms(i,k),i,k))
        
        xG   = ( ( x  (:,i,k) - xCenter  (i,k) ) / dx   ) / recCoef
        etaG = ( ( eta(:,i,k) - etaCenter(i,k) ) / deta ) / recCoef
        call calc_polynomial_matrix(locPolyDegree(i,k),nQuadPointsOnCell,nRecTerms(i,k),xG,etaG,recMatrixG(:,1:nRecTerms(i,k),i,k))
      enddo
    enddo
    !$OMP END PARALLEL DO
    
    !! Reconstruct metric function
    !call reconstruct_field(sqrtGC*recdV,&
    !                       sqrtGL      ,&
    !                       sqrtGR      ,&
    !                       sqrtGB      ,&
    !                       sqrtGT)
    !
    !call reconstruct_field(G13C*recdV,&
    !                       G13L      ,&
    !                       G13R      ,&
    !                       G13B      ,&
    !                       G13T)
    
        
    ! Set mertric functions by analytical value
    do k = kds,kde
      do i = ids,ide
        sqrtGL(:,i,k) = sqrtG_ext(:,4,i,k)
        sqrtGR(:,i,k) = sqrtG_ext(:,2,i,k)
        sqrtGB(:,i,k) = sqrtG_ext(:,1,i,k)
        sqrtGT(:,i,k) = sqrtG_ext(:,3,i,k)
        
        G13L  (:,i,k) = G13_ext  (:,4,i,k)
        G13R  (:,i,k) = G13_ext  (:,2,i,k)
        G13B  (:,i,k) = G13_ext  (:,1,i,k)
        G13T  (:,i,k) = G13_ext  (:,3,i,k)
      enddo
    enddo
    
    ! Calculate reference pressure
    qC = FillValue
    qC(:,ids:ide,kds:kde) = ref%q(:,ids:ide,kds:kde)
    do iVar = 1,nVar
      call reconstruct_field(qC(iVar  ,:,:)*recdV,&
                             qL(iVar,:,:,:)      ,&
                             qR(iVar,:,:,:)      ,&
                             qB(iVar,:,:,:)      ,&
                             qT(iVar,:,:,:)      ,&
                             qG(iVar,:,:,:))
    enddo
    
    qL_ref = qL
    qR_ref = qR
    qB_ref = qB
    qT_ref = qT
    qG_ref = qG
  
    ! Fill outside boundary
    ! left boundary
    qR_ref(:,:,ids-1,kds:kde) = qL_ref(:,:,ids,kds:kde)
    sqrtGR(  :,ids-1,kds:kde) = sqrtGL(  :,ids,kds:kde)
    G13R  (  :,ids-1,kds:kde) = G13L  (  :,ids,kds:kde)
    
    ! right boundary
    qL_ref(:,:,ide+1,kds:kde) = qR_ref(:,:,ide,kds:kde)
    sqrtGL(  :,ide+1,kds:kde) = sqrtGR(  :,ide,kds:kde)
    G13L  (  :,ide+1,kds:kde) = G13R  (  :,ide,kds:kde)
    
    ! bottom boundary
    qT_ref(:,:,ids:ide,kds-1) = qB_ref(:,:,ids:ide,kds)
    sqrtGT(  :,ids:ide,kds-1) = sqrtGB(  :,ids:ide,kds)
    G13T  (  :,ids:ide,kds-1) = G13B  (  :,ids:ide,kds)
    
    ! top boundary
    qB_ref(:,:,ids:ide,kde+1) = qT_ref(:,:,ids:ide,kde)
    sqrtGB(  :,ids:ide,kde+1) = sqrtGT(  :,ids:ide,kde)
    G13B  (  :,ids:ide,kde+1) = G13T  (  :,ids:ide,kde)
    
    ! Calculate projection matrix for no-flux boundary
    k    = kds
    iEOC = 1
    do i = ids,ide
      do iPOE = 1,nPointsOnEdge
        if(G13B(iPOE,i,k)/=0)then
          nx = 1
          nz = 1. / G13B(iPOE,i,k)
        else
          nx = 0
          nz = 1
        endif
        nv(1) = nx
        nv(2) = nz
        nv    = nv / sqrt( dot_product( nv, nv ) ) ! calc unit norm vector
        nx    = nv(1)
        nz    = nv(2)
        
        pm(1,1) = 1. - nx**2
        pm(1,2) = -nx * nz
        pm(2,1) = pm(1,2)
        pm(2,2) = 1. - nz**2
        
        nvec(:  ,iPOE,iEOC,i,k) = nv
        pmtx(:,:,iPOE,iEOC,i,k) = pm
      enddo
    enddo
    
    k    = kde
    iEOC = 3
    do i = ids,ide
      do iPOE = 1,nPointsOnEdge
        if(G13T(iPOE,i,k)/=0)then
          nx = 1
          nz = 1. / G13T(iPOE,i,k)
        else
          nx = 0
          nz = 1
        endif
        nv(1) = nx
        nv(2) = nz
        nv    = nv / sqrt( dot_product( nv, nv ) ) ! calc unit norm vector
        nx    = nv(1)
        nz    = nv(2)
        
        pm(1,1) = 1. - nx**2
        pm(1,2) = -nx * nz
        pm(2,1) = pm(1,2)
        pm(2,2) = 1. - nz**2
        
        nvec(:  ,iPOE,iEOC,i,k) = nv
        pmtx(:,:,iPOE,iEOC,i,k) = pm
      enddo
    enddo
    
    !$OMP PARALLEL DO PRIVATE(i,iPOE)
    do k = kds-1,kde+1
      do i = ids-1,ide+1
        do iPOE = 1,nPointsOnEdge
          PL_ref(iPOE,i,k) = calc_pressure(sqrtGL(iPOE,i,k),qL_ref(:,iPOE,i,k))
          PR_ref(iPOE,i,k) = calc_pressure(sqrtGR(iPOE,i,k),qR_ref(:,iPOE,i,k))
          PB_ref(iPOE,i,k) = calc_pressure(sqrtGB(iPOE,i,k),qB_ref(:,iPOE,i,k))
          PT_ref(iPOE,i,k) = calc_pressure(sqrtGT(iPOE,i,k),qT_ref(:,iPOE,i,k))
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    
    ! Calculate Rayleigh damping coef
    call Rayleigh_coef(relax_coef)
    
    ! Fill out qC
    qC = FillValue
    
  end subroutine init_spatial_operator
  
  subroutine spatial_operator(stat,tend)
    type(stat_field), target, intent(inout) :: stat
    type(tend_field), target, intent(inout) :: tend
    
    real(r_kind)                  :: sqrtGe
    real(r_kind)                  :: G13e
    real(r_kind), dimension(nVar) :: qe
    real(r_kind)                  :: rho
    real(r_kind)                  :: u
    real(r_kind)                  :: w
    real(r_kind)                  :: w_eta
    real(r_kind)                  :: theta
    
    integer(i_kind) :: i,j,k,iR,kR,iVar,iPOE,iEOC
    real   (r_kind) :: wind_vector(2)
    
    integer(i_kind) ip1,im1
    integer(i_kind) kp1,km1
    integer(i_kind) ip2,im2
    integer(i_kind) kp2,km2
    
    ! Attension stat is changed here!
    if(case_num==2)call Rayleigh_damping(stat%q,ref%q)
    
    ! copy stat
    !qC(:,ids:ide,kds:kde) = stat%q(:,ids:ide,kds:kde)
    !$OMP PARALLEL DO PRIVATE(i,iVar) COLLAPSE(3)
    do k = kds,kde
      do i = ids,ide
        do iVar = 1,nVar
          qC(iVar,i,k) = stat%q(iVar,i,k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    
    do iVar = 1,nVar
      call reconstruct_field(qC(iVar  ,:,:)*recdV,&
                             qL(iVar,:,:,:)      ,&
                             qR(iVar,:,:,:)      ,&
                             qB(iVar,:,:,:)      ,&
                             qT(iVar,:,:,:)      ,&
                             qG(iVar,:,:,:))
      !call reconstruct_field(qC(iVar  ,:,:)*recdV,&
      !                       qL(iVar,:,:,:)      ,&
      !                       qR(iVar,:,:,:)      ,&
      !                       qB(iVar,:,:,:)      ,&
      !                       qT(iVar,:,:,:))
      !call reconstruct_field(qC(iVar  ,:,:)*recdV,&
      !                       qL(iVar,:,:,:)      ,&
      !                       qR(iVar,:,:,:)      ,&
      !                       qB(iVar,:,:,:)      ,&
      !                       qT(iVar,:,:,:)      ,&
      !                       iVar=iVar)
    enddo
  
    ! Boundary Condition
    ! Correct wind on bottom and top boundaries
    k    = kds
    iEOC = 1
    !$OMP PARALLEL DO PRIVATE(iPOE,wind_vector) COLLAPSE(2)
    do i = ids,ide
      do iPOE = 1,nPointsOnEdge
        wind_vector(1) = qB(2,iPOE,i,k)
        wind_vector(2) = qB(3,iPOE,i,k)
        wind_vector    = matmul(pmtx(:,:,iPOE,iEOC,i,k),wind_vector)
        qB(2,iPOE,i,k) = wind_vector(1)
        qB(3,iPOE,i,k) = wind_vector(2)
      enddo
    enddo
    !$OMP END PARALLEL DO
    
    k    = kde
    iEOC = 3
    !$OMP PARALLEL DO PRIVATE(iPOE,wind_vector) COLLAPSE(2)
    do i = ids,ide
      do iPOE = 1,nPointsOnEdge
        wind_vector(1) = qT(2,iPOE,i,k)
        wind_vector(2) = qT(3,iPOE,i,k)
        wind_vector    = matmul(pmtx(:,:,iPOE,iEOC,i,k),wind_vector)
        qT(2,iPOE,i,k) = wind_vector(1)
        qT(3,iPOE,i,k) = wind_vector(2)
      enddo
    enddo
    !$OMP END PARALLEL DO
      
    ! Fill lateral boundary
    if(case_num==1.or.case_num==3)then
      qL(2,:,ids,kds:kde) = 0
      qR(2,:,ide,kds:kde) = 0
    elseif(case_num==2)then
      !$OMP PARALLEL DO PRIVATE(iPOE) COLLAPSE(2)
      do k = kds,kde
        do iPOE = 1,nPointsOnEdge
          qL(2,iPOE,ids,k) = ref%q(2,ids,k)
          qR(2,iPOE,ide,k) = ref%q(2,ide,k)
        enddo
      enddo
      !$OMP END PARALLEL DO
    endif
    
    ! Start OpenMP zone
    !$OMP PARALLEL
    
    ! Fill outside boundary
    !$OMP DO
    do k = kds,kde
      ! left boundary
      qR(:,:,ids-1,k) = qL(:,:,ids,k)
      
      ! right boundary
      qL(:,:,ide+1,k) = qR(:,:,ide,k)
    enddo
    !$OMP ENDDO
    
    !$OMP DO
    do i = ids,ide
      ! bottom boundary
      qT(:,:,i,kds-1) = qB(:,:,i,kds)
      
      ! top boundary
      qB(:,:,i,kde+1) = qT(:,:,i,kde)
    enddo
    !$OMP ENDDO
  
    ! Calculate pressure X
    !$OMP DO PRIVATE(i,iPOE) COLLAPSE(3)
    do k = kds,kde
      do i = ids,ide
        do iPOE = 1,nPointsOnEdge
          PL(iPOE,i,k) = calc_pressure(sqrtGL(iPOE,i,k),qL(:,iPOE,i,k))
          PR(iPOE,i,k) = calc_pressure(sqrtGR(iPOE,i,k),qR(:,iPOE,i,k))
        enddo
      enddo
    enddo
    !$OMP END DO NOWAIT
    !$OMP DO PRIVATE(iPOE,i) COLLAPSE(2)
    do k = kds,kde
      do iPOE = 1,nPointsOnEdge
        i = ide + 1
        PL(iPOE,i,k) = calc_pressure(sqrtGL(iPOE,i,k),qL(:,iPOE,i,k))
        i = ids - 1
        PR(iPOE,i,k) = calc_pressure(sqrtGR(iPOE,i,k),qR(:,iPOE,i,k))
      enddo
    enddo
    !$OMP ENDDO NOWAIT
    
    ! Calculate pressure Z
    !$OMP DO PRIVATE(i,iPOE) COLLAPSE(3)
    do k = kds,kde
      do i = ids,ide
        do iPOE = 1,nPointsOnEdge
          PB(iPOE,i,k) = calc_pressure(sqrtGB(iPOE,i,k),qB(:,iPOE,i,k))
          PT(iPOE,i,k) = calc_pressure(sqrtGT(iPOE,i,k),qT(:,iPOE,i,k))
        enddo
      enddo
    enddo
    !$OMP END DO NOWAIT
    !$OMP DO PRIVATE(k,iPOE) COLLAPSE(2)
    do i = ids,ide
      do iPOE = 1,nPointsOnEdge
        k = kde + 1
        PB(iPOE,i,k) = calc_pressure(sqrtGB(iPOE,i,k),qB(:,iPOE,i,k))
        k = kds - 1
        PT(iPOE,i,k) = calc_pressure(sqrtGT(iPOE,i,k),qT(:,iPOE,i,k))
      enddo
    enddo
    !$OMP END DO NOWAIT
    
    !$OMP BARRIER
    
    ! Calculate F
    !$OMP DO PRIVATE(i,iPOE) COLLAPSE(3)
    do k = kds,kde
      do i = ids,ide
        do iPOE = 1,nPointsOnEdge
          FL(:,iPOE,i,k) = calc_F(sqrtGL(iPOE,i,k),qL(:,iPOE,i,k),PL(iPOE,i,k))
          FR(:,iPOE,i,k) = calc_F(sqrtGR(iPOE,i,k),qR(:,iPOE,i,k),PR(iPOE,i,k))
        enddo
      enddo
    enddo
    !$OMP END DO NOWAIT
    !$OMP DO PRIVATE(i,iPOE) COLLAPSE(2)
    do k = kds,kde
      do iPOE = 1,nPointsOnEdge
        i = ide + 1
        FL(:,iPOE,i,k) = calc_F(sqrtGL(iPOE,i,k),qL(:,iPOE,i,k),PL(iPOE,i,k))
        i = ids - 1
        FR(:,iPOE,i,k) = calc_F(sqrtGR(iPOE,i,k),qR(:,iPOE,i,k),PR(iPOE,i,k))
      enddo
    enddo
    !$OMP END DO NOWAIT
    
    ! Calculate pressure H
    !$OMP DO PRIVATE(i,iPOE) COLLAPSE(3)
    do k = kds,kde
      do i = ids,ide
        do iPOE = 1,nPointsOnEdge
          HB(:,iPOE,i,k) = calc_H(sqrtGB(iPOE,i,k),G13B(iPOE,i,k),qB(:,iPOE,i,k),pB(iPOE,i,k),pB_ref(iPOE,i,k))
          HT(:,iPOE,i,k) = calc_H(sqrtGT(iPOE,i,k),G13T(iPOE,i,k),qT(:,iPOE,i,k),pT(iPOE,i,k),pT_ref(iPOE,i,k))
        enddo
      enddo
    enddo
    !$OMP END DO NOWAIT
    !$OMP DO PRIVATE(k,iPOE) COLLAPSE(2)
    do i = ids,ide
      do iPOE = 1,nPointsOnEdge
        k = kde + 1
        HB(:,iPOE,i,k) = calc_H(sqrtGB(iPOE,i,k),G13B(iPOE,i,k),qB(:,iPOE,i,k),pB(iPOE,i,k),pB_ref(iPOE,i,k))
        k = kds - 1
        HT(:,iPOE,i,k) = calc_H(sqrtGT(iPOE,i,k),G13T(iPOE,i,k),qT(:,iPOE,i,k),pT(iPOE,i,k),pT_ref(iPOE,i,k))
      enddo
    enddo
    !$OMP END DO NOWAIT
    
    !$OMP BARRIER
    
    ! calc x flux
    !$OMP DO PRIVATE(i,im1,iPOE,eigen_mtx_x,qe,u,w,theta,sqrtGe,iVar) COLLAPSE(2)
    do k = kds,kde
      do i = ids,ide+1
        im1 = i - 1
        do iPOE = 1,nPointsOnEdge
          ! Scheme 1
          qe = 0.5 * ( qL(:,iPOE,i,k) + qR(:,iPOE,im1,k) )
          !u  = qe(2) / qe(1)
          !qe = qe - 0.5 * sign(1._r_kind,u) * ( qL(:,iPOE,i,k) - qR(:,iPOE,im1,k) )
          
          !! Scheme 2
          !rho   = ( 0.5 * ( sqrt(qL(1,iPOE,i,k)) + sqrt(qR(1,iPOE,im1,k)) ) )**2
          !u     = ( qL(2,iPOE,i,k) / sqrt(qL(1,iPOE,i,k)) + qR(2,iPOE,im1,k) / sqrt(qR(1,iPOE,im1,k)) ) / ( sqrt(qL(1,iPOE,i,k)) + sqrt(qR(1,iPOE,im1,k)) )
          !w     = ( qL(3,iPOE,i,k) / sqrt(qL(1,iPOE,i,k)) + qR(3,iPOE,im1,k) / sqrt(qR(1,iPOE,im1,k)) ) / ( sqrt(qL(1,iPOE,i,k)) + sqrt(qR(1,iPOE,im1,k)) )
          !theta = ( qL(4,iPOE,i,k) / sqrt(qL(1,iPOE,i,k)) + qR(4,iPOE,im1,k) / sqrt(qR(1,iPOE,im1,k)) ) / ( sqrt(qL(1,iPOE,i,k)) + sqrt(qR(1,iPOE,im1,k)) )
          !
          !qe(1) = rho
          !qe(2) = rho * u
          !qe(3) = rho * w
          !qe(4) = rho * theta
          
          sqrtGe = 0.5 * ( sqrtGL(iPOE,i,k) + sqrtGR(iPOE,im1,k) )
          
          eigen_mtx_x = calc_eigen_matrix_x( qe(1:4), sqrtGe )
          
          FeP(:,iPOE,i,k) = 0.5 * ( FL(1:4,iPOE,i,k) + FR(1:4,iPOE,im1,k) - matmul( eigen_mtx_x, ( qL(1:4,iPOE,i,k) - qR(1:4,iPOE,im1,k) ) ) )
        enddo
        do iVar = 1,4
          Fe(iVar,i,k) = Gaussian_quadrature_1d(FeP(iVar,:,i,k))
        enddo
      enddo
    enddo
    !$OMP END DO NOWAIT
    
    ! calc z flux
    !$OMP DO PRIVATE(k,km1,iPOE,u,w,theta,w_eta,eigen_mtx_z,qe,sqrtGe,G13e,iVar) COLLAPSE(2)
    do i = ids,ide
      do k = kds,kde+1
        km1 = k - 1
        do iPOE = 1,nPointsOnEdge
          ! Scheme 1
          qe          = 0.5 * ( qB    (:,iPOE,i,k) + qT    (:,iPOE,i,km1) )
          sqrtGe      = 0.5 * ( sqrtGB(  iPOE,i,k) + sqrtGT(  iPOE,i,km1) )
          G13e        = 0.5 * ( G13B  (  iPOE,i,k) + G13T  (  iPOE,i,km1) )
          !w_eta       = calc_w_eta(sqrtGe,G13e,qe)
          !qe          = qe - 0.5 * sign(1._r_kind,w_eta) * ( qB(:,iPOE,i,k) - qT(:,iPOE,i,km1) )
          
          !! Scheme 2
          !rho   = ( 0.5 * ( sqrt(qB(1,iPOE,i,k)) + sqrt(qT(1,iPOE,i,km1)) ) )**2
          !u     = ( qB(2,iPOE,i,k) / sqrt(qB(1,iPOE,i,k)) + qT(2,iPOE,i,km1) / sqrt(qT(1,iPOE,i,km1)) ) / ( sqrt(qB(1,iPOE,i,k)) + sqrt(qT(1,iPOE,i,km1)) )
          !w     = ( qB(3,iPOE,i,k) / sqrt(qB(1,iPOE,i,k)) + qT(3,iPOE,i,km1) / sqrt(qT(1,iPOE,i,km1)) ) / ( sqrt(qB(1,iPOE,i,k)) + sqrt(qT(1,iPOE,i,km1)) )
          !theta = ( qB(4,iPOE,i,k) / sqrt(qB(1,iPOE,i,k)) + qT(4,iPOE,i,km1) / sqrt(qT(1,iPOE,i,km1)) ) / ( sqrt(qB(1,iPOE,i,k)) + sqrt(qT(1,iPOE,i,km1)) )
          !
          !qe(1) = rho
          !qe(2) = rho * u
          !qe(3) = rho * w
          !qe(4) = rho * theta
          !
          !sqrtGe = 0.5 * ( sqrtGB(iPOE,i,k) + sqrtGT(iPOE,i,km1) )
          !G13e   = 0.5 * ( G13B  (iPOE,i,k) + G13T  (iPOE,i,km1) )
          
          eigen_mtx_z = calc_eigen_matrix_z( qe(1:4), sqrtGe, G13e )
          
          HeP(1:4,iPOE,i,k) = 0.5 * ( HB(1:4,iPOE,i,k) + HT(1:4,iPOE,i,km1) - matmul( eigen_mtx_z, ( qB(1:4,iPOE,i,k) - qT(1:4,iPOE,i,km1) ) ) )
        enddo
        do iVar = 1,4
          He(iVar,i,k) = Gaussian_quadrature_1d(HeP(iVar,:,i,k))
        enddo
      enddo
    enddo
    !$OMP END DO NOWAIT
    
    !$OMP DO PRIVATE(i) COLLAPSE(2)
    do k = kds,kde
      do i = ids,ide
        !rho_p(i,k) = ( qC(1,i,k) + qC(5,i,k) - ref%q(1,i,k) - ref%q(5,i,k) ) / sqrtGC(i,k)  ! hydrostatic
        !rho_p(i,k) = ( qC(1,i,k) + qC(5,i,k) ) / sqrtGC(i,k)  ! nonhydrostatic
        rho_p(i,k) = Gaussian_quadrature_2d( ( qG(1,:,i,k) + qG(5,:,i,k) - qG_ref(1,:,i,k) - qG_ref(5,:,i,k) ) ) ! hydrostatic
        if( abs( rho_p(i,k) ) <= 1.e-13 ) rho_p(i,k) = 0 ! hydrostatic
      enddo
    enddo
    !$OMP END DO
    
    !$OMP DO PRIVATE(i,iVar) COLLAPSE(2)
    do k = kds,kde
      do i = ids,ide
        do iVar = 1,nVar
          src(iVar,i,k) = 0
        enddo
        !src(3,i,k) = src(3,i,k) - sqrtGC(i,k) * rho_p(i,k) * gravity
        src(3,i,k) = src(3,i,k) - rho_p(i,k) * gravity
      enddo
    enddo
    !$OMP END DO
    
    !$OMP END PARALLEL
    ! End OpenMP zone
  
    ! Viscosity terms for Density Current case only
    if(case_num==3)then
      do iVar = 2,4
        q_diff(iVar,ids:ide,kds:kde) = qC(iVar,ids:ide,kds:kde) / qC(1,ids:ide,kds:kde)
        
        q_diff(iVar,ids:ide,kds-1) = q_diff(iVar,ids:ide,kds)
        q_diff(iVar,ids:ide,kde+1) = q_diff(iVar,ids:ide,kde)
        q_diff(iVar,ids-1,kds:kde) = q_diff(iVar,ids,kds:kde)
        q_diff(iVar,ide+1,kds:kde) = q_diff(iVar,ide,kds:kde)
      enddo
      
      !$OMP PARALLEL DO PRIVATE(i,iVar)
      do k = kds,kde
        do i = ids,ide
          do iVar = 2,4
            src(iVar,i,k) = src(iVar,i,k) + viscosity_coef * qC(1,i,k) * ( ( q_diff(iVar,i+1,k) - 2. * q_diff(iVar,i,k) + q_diff(iVar,i-1,k) ) / dx  **2 &
                                                                         + ( q_diff(iVar,i,k+1) - 2. * q_diff(iVar,i,k) + q_diff(iVar,i,k-1) ) / deta**2 / sqrtGC(i,k)**2 )
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
    endif
    
    ! Calculate tend
    !$OMP PARALLEL DO PRIVATE(i,ip1,kp1,iVar) COLLAPSE(2)
    do k = kds,kde
      do i = ids,ide
        ip1 = i + 1
        kp1 = k + 1
        do iVar = 1,nVar
          tend%q(iVar,i,k) = - ( ( Fe(iVar,ip1,k) - Fe(iVar,i,k) ) / dx + ( He(iVar,i,kp1) - He(iVar,i,k) ) / deta ) + src(iVar,i,k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    
  end subroutine spatial_operator

  subroutine reconstruct_field(qC,qL,qR,qB,qT,qG,iVar)
    real   (r_kind), dimension(                  ics:ice,kcs:kce),intent(in )           :: qC
    real   (r_kind), dimension(nPointsOnEdge    ,ics:ice,kcs:kce),intent(out)           :: qL
    real   (r_kind), dimension(nPointsOnEdge    ,ics:ice,kcs:kce),intent(out)           :: qR
    real   (r_kind), dimension(nPointsOnEdge    ,ics:ice,kcs:kce),intent(out)           :: qB
    real   (r_kind), dimension(nPointsOnEdge    ,ics:ice,kcs:kce),intent(out)           :: qT
    real   (r_kind), dimension(nQuadPointsOnCell,ics:ice,kcs:kce),intent(out), optional :: qG
    integer(i_kind)                                                          , optional :: iVar ! Specify variable for speed up CG method
  
    integer(i_kind) :: i,j,k
    integer(i_kind) :: iRec,kRec
    integer(i_kind) :: ic
    integer(i_kind) :: m,n
    
    real(r_kind)                                     :: h
    real(r_kind), dimension(maxRecCells            ) :: u
    real(r_kind), dimension(maxRecCells,maxRecTerms) :: A
    real(r_kind), dimension(            maxRecTerms) :: polyCoef
    
    h = dx
    
    ! Full WLS-ENO
    !$OMP PARALLEL DO PRIVATE(i,j,m,n,iRec,kRec,iC,u,A,polyCoef) collapse(2)
    do k = kds,kde
      do i = ids,ide
        m = nRecCells(i,k)
        n = nRecTerms(i,k)
        ! Set variable for reconstruction
        do j = 1,m
          iRec = iRecCell(j,i,k)
          kRec = kRecCell(j,i,k)
          
          u(j    ) = qC(iRec,kRec)
          A(j,1:n) = polyCoordCoef(j,1:n,i,k)
        enddo
        ic = iCenCell(i,k)
        
        !if(present(iVar))then
        !  polyCoef(1:n) = WLS_ENO(A(1:m,1:n),u(1:m),h,m,n,ic,CG_coef(iVar,1:n,i,k))
        !  CG_coef(iVar,1:n,i,k) = polyCoef(1:n)
        !else
        !  polyCoef(1:n) = WLS_ENO(A(1:m,1:n),u(1:m),h,m,n,ic)
        !endif
        polyCoef(1:n) = WLS_ENO(A(1:m,1:n),u(1:m),h,m,n,ic)
        
        qL(:,i,k) = matmul(recMatrixL(:,1:n,i,k),polyCoef(1:n))
        qR(:,i,k) = matmul(recMatrixR(:,1:n,i,k),polyCoef(1:n))
        qB(:,i,k) = matmul(recMatrixB(:,1:n,i,k),polyCoef(1:n))
        qT(:,i,k) = matmul(recMatrixT(:,1:n,i,k),polyCoef(1:n))
        
        if(present(qG))qG(:,i,k) = matmul(recMatrixG(:,1:n,i,k),polyCoef(1:n))
      enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine reconstruct_field
  
  function calc_pressure(sqrtG,q)
    real(r_kind) calc_pressure
    real(r_kind) sqrtG
    real(r_kind) q(5)
    
    real(r_kind) w1
    real(r_kind) w4
    real(r_kind) w5
    
    real(r_kind) rho
    real(r_kind) R
    real(r_kind) gamma
    real(r_kind) theta
    real(r_kind) cp
    real(r_kind) cv
    real(r_kind) sh
    real(r_kind) kappa
    
    w1 = q(1)
    w4 = q(4)
    w5 = q(5)
    
    rho      = ( w1 + w5 ) / sqrtG
    gamma    = w5 / w1
    sh       = gamma / ( 1. + gamma )
    R        = ( 1. + eq * sh ) * Rd
    theta    = w4 / w1
    cp       = cpd + ( cpv - cpd ) * sh
    cv       = cvd + ( cvv - cvd ) * sh
    kappa    = cp / cv
    
    calc_pressure = p0 * ( rho * theta * R / p0 )**kappa
    
  end function calc_pressure
    
  function calc_F(sqrtG,q,P)
    real(r_kind),dimension(nVar) :: calc_F
    real(r_kind)                 :: sqrtG
    real(r_kind),dimension(nVar) :: q
    real(r_kind)                 :: p      ! pressure
    
    real(r_kind) w1
    real(r_kind) w2
    real(r_kind) w3
    real(r_kind) w4
    
    real(r_kind) sqrtGrho
    real(r_kind) u
    
    w1 = q(1)
    w2 = q(2)
    w3 = q(3)
    w4 = q(4)
    
    sqrtGrho = w1
    u        = w2 / sqrtGrho
    
    calc_F(1) = w1 * u
    calc_F(2) = w2 * u + sqrtG * p
    calc_F(3) = w3 * u
    calc_F(4) = w4 * u
    
  end function calc_F
  
  function calc_H(sqrtG,G13,q,p,p_ref)
    real(r_kind),dimension(nVar) :: calc_H
    real(r_kind)                 :: sqrtG
    real(r_kind)                 :: G13
    real(r_kind),dimension(nVar) :: q
    real(r_kind)                 :: p      ! pressure
    real(r_kind)                 :: p_ref  ! reference pressure
    
    real(r_kind)                 :: p_pert ! pressure  perturbation
    
    real(r_kind) w1
    real(r_kind) w2
    real(r_kind) w3
    real(r_kind) w4
    real(r_kind) ww
    
    real(r_kind) sqrtGrho
    real(r_kind) u
    real(r_kind) w
    
    w1 = q(1)
    w2 = q(2)
    w3 = q(3)
    w4 = q(4)
    
    sqrtGrho = w1
    u        = w2 / sqrtGrho
    w        = w3 / sqrtGrho
    p_pert   = p - p_ref
    if(abs(p_pert)/p_ref<1.e-13)p_pert=0
    
    ww = w / sqrtG + G13 * u
    
    calc_H(1) = w1 * ww
    calc_H(2) = w2 * ww + sqrtG * G13 * p
    calc_H(3) = w3 * ww + p_pert
    calc_H(4) = w4 * ww
    
  end function calc_H
    
  function calc_w_eta(sqrtG,G13,q)
    real(r_kind)                 :: calc_w_eta
    real(r_kind)                 :: sqrtG
    real(r_kind)                 :: G13
    real(r_kind),dimension(nVar) :: q
    
    real(r_kind) :: u
    real(r_kind) :: w

    u = q(2) / q(1)
    w = q(3) / q(1)
    
    calc_w_eta = w / sqrtG + G13 * u
  
  end function calc_w_eta
  
  subroutine Rayleigh_damping(q,q_ref)
    real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(inout) :: q
    real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(in   ) :: q_ref
    
    integer(i_kind), parameter :: vs = 1
    integer(i_kind), parameter :: ve = 5
    
    integer i,k,iVar
    
    !$OMP PARALLEL DO PRIVATE(i,iVar) COLLAPSE(3)
    do k = kds,kde
      do i = ids,ide
        do iVar = vs,ve
          q(iVar,i,k) = q(iVar,i,k) - relax_coef(i,k) * ( q(iVar,i,k) - q_ref(iVar,i,k) )
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    
  end subroutine Rayleigh_damping
  
  ! Rayleigh damping ( Wong and Stull, MWR, 2015 )
  subroutine Rayleigh_coef(mu)
    real(r_kind), dimension(ics:ice,kcs:kce), intent(out) :: mu
    
    real(r_kind), dimension(ics:ice,kcs:kce) :: muT
    real(r_kind), dimension(ics:ice,kcs:kce) :: muL
    real(r_kind), dimension(ics:ice,kcs:kce) :: muR
    
    real(r_kind), parameter :: topSpongeThickness   = 9000   ! 12000 for "best" result
    real(r_kind), parameter :: leftSpongeThickness  = 10000  ! 10000 for "best" result
    real(r_kind), parameter :: rightSpongeThickness = 10000  ! 10000 for "best" result
    
    real(r_kind), parameter :: mu_max_top = 0.2!0.28
    real(r_kind), parameter :: mu_max_lat = 0.2!0.15
    
    real(r_kind) zd, zt
    
    integer i,k
    
    muT = 0
    muL = 0
    MuR = 0
    
    ! Top
    zt = z_max
    zd = zt - topSpongeThickness
    where( zC > zd )
      !muT = mu_max_top * sin( pi / 2. * ( zC - zd ) / ( zt - zd ) )**2  !( Wong and Stull, MWR, 2015 )
      muT = mu_max_top * ( ( zC - zd ) / ( zt - zd ) )**4 ! ( Li Xingliang, MWR, 2013 )
    elsewhere
      muT = 0.
    endwhere
    
    ! Left
    zt = -x_min
    zd = zt - leftSpongeThickness
    where( abs(xC) > zd )
      !muT = mu_max_lat * sin( pi / 2. * ( abs(xC) - zd ) / ( zt - zd ) )**2  !( Wong and Stull, MWR, 2015 )
      muT = mu_max_lat * ( ( abs(xC) - zd ) / ( zt - zd ) )**4 ! ( Li Xingliang, MWR, 2013 )
    elsewhere
      muL = 0.
    endwhere
    
    ! Right
    zt = x_max
    zd = zt - rightSpongeThickness
    where( xC > zd )
      !muT = mu_max_lat * sin( pi / 2. * ( abs(xC) - zd ) / ( zt - zd ) )**2  !( Wong and Stull, MWR, 2015 )
      muT = mu_max_lat * ( ( abs(xC) - zd ) / ( zt - zd ) )**4 ! ( Li Xingliang, MWR, 2013 )
    elsewhere
      muR = 0.
    endwhere
    
    do k = kds,kde
      do i = ids,ide
        mu(i,k) = max( muT(i,k), muL(i,k), muR(i,k) )
      enddo
    enddo
  
  end subroutine Rayleigh_coef
  
  function calc_eigenvalue_x(sqrtG,q)
    !real(r_kind),dimension(nVar) :: calc_eigenvalue_x
    real(r_kind),dimension(4) :: calc_eigenvalue_x
    real(r_kind)                 :: sqrtG
    !real(r_kind),dimension(nVar) :: q
    real(r_kind),dimension(4) :: q
    
    real(r_kind) w1
    real(r_kind) w2
    real(r_kind) w3
    real(r_kind) w4
    
    !real(r_kind) eig(nVar)
    real(r_kind) eig(4)
    
    real(r_kind) coef1,coef2,coef3
    
    if(any(q==FillValue))then
      eig = 0
    else
      w1 = q(1)
      w2 = q(2)
      w3 = q(3)
      w4 = q(4)
      
      coef1 = w2
      
      coef2 = sqrt( cpd * p0 * sqrtG * w1 / cvd ) * ((Rd*w4)/(p0*sqrtG))**(cpd/(2*cvd))
      
      coef3 = w1
      
      eig(1) = w2 / w1
      eig(2) = eig(1)
      eig(3) = ( coef1 - coef2 ) / coef3
      eig(4) = ( coef1 + coef2 ) / coef3
    endif
    
    calc_eigenvalue_x = eig
    
  end function calc_eigenvalue_x
  
  function calc_eigenvalue_z(sqrtG,G13,q)
    !real(r_kind),dimension(nVar) :: calc_eigenvalue_z
    real(r_kind),dimension(4) :: calc_eigenvalue_z
    real(r_kind)                 :: sqrtG
    real(r_kind)                 :: G13
    !real(r_kind),dimension(nVar) :: q
    real(r_kind),dimension(4) :: q
    
    real(r_kind) w1
    real(r_kind) w2
    real(r_kind) w3
    real(r_kind) w4
    
    !real(r_kind) eig(nVar)
    real(r_kind) eig(4)
    
    real(r_kind) sqrtGrho
    real(r_kind) u
    real(r_kind) w
    
    real(r_kind) drhoetadt
    
    real(r_kind) coef1,coef2,coef3
    
    if(any(q==FillValue))then
      eig = 0
    else
      w1 = q(1)
      w2 = q(2)
      w3 = q(3)
      w4 = q(4)
      
      coef1 = cvd*sqrtG**2*w1**3*(G13*sqrtG*w2 + w3)*w4
      
      coef2 = sqrt( cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7*w4**2&
                  *((Rd*w4)/(p0*sqrtG))**(cpd/cvd) )
      
      coef3 = cvd*sqrtG**3*w1**4*w4
      
      drhoetadt = (G13*sqrtG*w2 + w3)/ (sqrtG*w1)
      
      eig(1) = drhoetadt
      eig(2) = drhoetadt
      eig(3) = ( coef1 - coef2 ) / coef3
      eig(4) = ( coef1 + coef2 ) / coef3
      
      !eig(3) = (sqrtG**2*w1**3*(G13*sqrtG*w2 + w3) -                           &
      !            Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7*w4**2*  &
      !                  ((Rd*w4)/(p0*sqrtG))**(cpd/cvd))/(cvd*w4))/            &
      !         (sqrtG**3*w1**4)
      !eig(4) = (cvd*sqrtG**2*w1**3*(G13*sqrtG*w2 + w3)*w4 +                    &
      !             Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7*w4**2* &
      !                 ((Rd*w4)/(p0*sqrtG))**(cpd/cvd)))/                      &
      !          (cvd*sqrtG**3*w1**4*w4)
    endif
    
    calc_eigenvalue_z = eig
    
  end function calc_eigenvalue_z
  
  ! A_x = R_x \lambda_x L_x, L_x = R_x^-1
  function calc_eigen_matrix_x(q,sqrtG)
    !real(r_kind), dimension(nVar,nVar)              :: calc_eigen_matrix_x
    real(r_kind), dimension(4,4)              :: calc_eigen_matrix_x
    !real(r_kind), dimension(nVar     ), intent(in ) :: q
    real(r_kind), dimension(4     ), intent(in ) :: q
    real(r_kind),                       intent(in ) :: sqrtG
    
    !real(r_kind), dimension(nVar,nVar)              :: mtx
    real(r_kind), dimension(4,4)              :: mtx
    
    real(r_kind) w1
    real(r_kind) w2
    real(r_kind) w3
    real(r_kind) w4
    
    real(r_kind) a,b,c
    
    !real(r_kind), dimension(nVar) :: eigen_value
    real(r_kind), dimension(4) :: eigen_value
    
    w1 = q(1)
    w2 = q(2)
    w3 = q(3)
    w4 = q(4)
    
    eigen_value = calc_eigenvalue_x(sqrtG,q)
    a = abs( eigen_value(1) )
    b = abs( eigen_value(3) )
    c = abs( eigen_value(4) )
    
    mtx(1,1) = (b*Sqrt(cvd)*w2 - c*Sqrt(cvd)*w2 + 2.*a*Sqrt(cpd)*Sqrt(p0)* &
                 Sqrt(sqrtG)*Sqrt(w1)*((Rd*w4)/(p0*sqrtG))**               &
                   (cpd/(2.*cvd)))/(((Rd*w4)/(p0*sqrtG))**(cpd/(2.*cvd))*  &
               (2.*Sqrt(cpd)*Sqrt(p0)*Sqrt(sqrtG)*Sqrt(w1)))
    
    mtx(1,2) = ((-b + c)*Sqrt(cvd)*Sqrt(w1))/         &
               (((Rd*w4)/(p0*sqrtG))**(cpd/(2.*cvd))* &
                  (2.*Sqrt(cpd)*Sqrt(p0)*Sqrt(sqrtG)))
    
    mtx(1,3) = 0
    
    mtx(1,4) = ((-2.*a + b + c)*w1)/(2.*w4)
    
    mtx(2,1) = (w2*(2.*a*Sqrt(w1) - b*Sqrt(w1) - c*Sqrt(w1) +           &
               (b*Sqrt(cvd)*w2)/(((Rd*w4)/(p0*sqrtG))**(cpd/(2.*cvd))*  &
                    (Sqrt(cpd)*Sqrt(p0)*Sqrt(sqrtG))) -                 &
               (c*Sqrt(cvd)*w2)/(((Rd*w4)/(p0*sqrtG))**(cpd/(2.*cvd))*  &
                    (Sqrt(cpd)*Sqrt(p0)*Sqrt(sqrtG)))))/(2.*w1**(3./2.))
    
    mtx(2,2) = (1./2.)*(b + c - (b*Sqrt(cvd)*w2)/                      &
                 (((Rd*w4)/(p0*sqrtG))**(cpd/(2.*cvd))*                &
                    (Sqrt(cpd)*Sqrt(p0)*Sqrt(sqrtG)*Sqrt(w1))) +       &
               (c*Sqrt(cvd)*w2)/(((Rd*w4)/(p0*sqrtG))**(cpd/(2.*cvd))* &
                    (Sqrt(cpd)*Sqrt(p0)*Sqrt(sqrtG)*Sqrt(w1))))
    
    mtx(2,3) = 0
    
    mtx(2,4) = -((2.*a*w2 - b*w2 - c*w2 + (b*Sqrt(cpd)*Sqrt(p0)*         &
                       Sqrt(sqrtG)*Sqrt(w1)*((Rd*w4)/(p0*sqrtG))**       &
                         (cpd/(2.*cvd)))/Sqrt(cvd) -                     &
                  (c*Sqrt(cpd)*Sqrt(p0)*Sqrt(sqrtG)*Sqrt(w1)*            &
                       ((Rd*w4)/(p0*sqrtG))**(cpd/(2.*cvd)))/Sqrt(cvd))/ &
               (2.*w4))
    
    mtx(3,1) = ((b - c)*Sqrt(cvd)*w2*w3)/                         &
               (((Rd*w4)/(p0*sqrtG))**(cpd/(2.*cvd))*             &
                  (2.*Sqrt(cpd)*Sqrt(p0)*Sqrt(sqrtG)*w1**(3./2.)))
    
    mtx(3,2) = -(((b - c)*Sqrt(cvd)*w3)/(((Rd*w4)/(p0*sqrtG))**   &
               (cpd/(2.*cvd))*(2.*Sqrt(cpd)*Sqrt(p0)*Sqrt(sqrtG)* &
                Sqrt(w1))))
    
    mtx(3,3) = a
    
    mtx(3,4) = ((-2.*a + b + c)*w3)/(2.*w4)
    
    mtx(4,1) = ((b - c)*Sqrt(cvd)*w2*w4)/                         &
               (((Rd*w4)/(p0*sqrtG))**(cpd/(2.*cvd))*             &
                  (2.*Sqrt(cpd)*Sqrt(p0)*Sqrt(sqrtG)*w1**(3./2.)))
    
    mtx(4,2) = -(((b - c)*Sqrt(cvd)*w4)/(((Rd*w4)/(p0*sqrtG))**   &
               (cpd/(2.*cvd))*(2.*Sqrt(cpd)*Sqrt(p0)*Sqrt(sqrtG)* &
                Sqrt(w1))))
          
    mtx(4,3) = 0
    
    mtx(4,4) = ( b + c ) / 2.
    
    calc_eigen_matrix_x = mtx
  end function calc_eigen_matrix_x
  
  ! A_z = R_z \lambda_z L_z, L_z = R_z^-1
  function calc_eigen_matrix_z(q,sqrtG,G13)
    !real(r_kind), dimension(nVar,nVar)              :: calc_eigen_matrix_z
    real(r_kind), dimension(4,4)              :: calc_eigen_matrix_z
    !real(r_kind), dimension(nVar     ), intent(in ) :: q
    real(r_kind), dimension(4     ), intent(in ) :: q
    real(r_kind),                       intent(in ) :: sqrtG
    real(r_kind),                       intent(in ) :: G13
    
    !real(r_kind), dimension(nVar,nVar)              :: mtx
    real(r_kind), dimension(4,4)              :: mtx
    
    real(r_kind) w1
    real(r_kind) w2
    real(r_kind) w3
    real(r_kind) w4
    
    real(r_kind) a,b,c
    
    !real(r_kind), dimension(nVar) :: eigen_value
    real(r_kind), dimension(4) :: eigen_value
    
    w1 = q(1)
    w2 = q(2)
    w3 = q(3)
    w4 = q(4)
    
    eigen_value = calc_eigenvalue_z(sqrtG,G13,q)
    a = abs( eigen_value(1) )
    b = abs( eigen_value(3) )
    c = abs( eigen_value(4) )
    
    mtx(1,1) = (b*cvd*sqrtG**2*w1**3*(G13*sqrtG*w2 + w3)*w4 -                &
                  c*cvd*sqrtG**2*w1**3*(G13*sqrtG*w2 + w3)*w4 +              &
                  2*a*Sqrt(cpd*cvd*p0*sqrtG**5*(1. + G13**2*sqrtG**2)*w1**7* &
                        w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)))/             &
               (2.*Sqrt(cpd*cvd*p0*sqrtG**5*(1. + G13**2*sqrtG**2)*w1**7*    &
                      w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)))
        
    mtx(1,2) = -(((b - c)*cvd*G13*sqrtG**3*w1**4*w4)/                      &
               (2.*Sqrt(cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)*w1**7*   &
                      w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd))))
    
    mtx(1,3) = -(((b - c)*cvd*sqrtG**2*w1**4*w4)/                       &
              (2.*Sqrt(cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)*w1**7* &
                     w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd))))
    
    mtx(1,4) = ((-2.*a + b + c)*w1)/(2.*w4)
    
    mtx(2,1) = -((sqrtG*(G13*sqrtG*w2 + w3)*                                 &
                 ((-b)*cvd*sqrtG*(1.+ G13**2*sqrtG**2)*w1**3*w2*w4 +         &
                    c*cvd*sqrtG*(1.+ G13**2*sqrtG**2)*w1**3*w2*w4 -          &
                    2.*a*G13*Sqrt(cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)* &
                          w1**7*w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)) +     &
                    b*G13*Sqrt(cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)*    &
                          w1**7*w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)) +     &
                    c*G13*Sqrt(cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)*    &
                          w1**7*w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd))))/    &
              (2.*(w1 + G13**2*sqrtG**2*w1)*                                 &
                 Sqrt(cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)*w1**7*w4**2* &
                     ((Rd*w4)/(p0*sqrtG))**(cpd/cvd))))
    
    mtx(2,2) = (2.*a*Sqrt(                                                     &
                cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)*w1**7*               &
                 w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)) -                      &
                 b*G13*sqrtG**2*(cvd*sqrtG*(1.+ G13**2*sqrtG**2)*w1**3*w2*w4 - &
                 G13*Sqrt(                                                     &
                   cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)*w1**7*            &
                    w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd))) +                  &
                 c*G13*sqrtG**2*(cvd*sqrtG*(1.+ G13**2*sqrtG**2)*w1**3*w2*w4 + &
                 G13*Sqrt(                                                     &
                   cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)*w1**7*            &
                    w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd))))/                  &
              (2.*(1.+ G13**2*sqrtG**2)*                                       &
              Sqrt(cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)*w1**7*            &
                w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)))
    
    mtx(2,3) = 0.5*sqrtG*(-((2.*a*G13)/(1.+ G13**2*sqrtG**2)) + (b*G13)/(1.+&
                G13**2*sqrtG**2) + (c*G13)/(1.+ G13**2*sqrtG**2) -          &
                (b*cvd*sqrtG*w1**3*w2*w4)/                                  &
              Sqrt(cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)*w1**7*         &
                w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)) +                    &
                (c*cvd*sqrtG*w1**3*w2*w4)/                                  &
              Sqrt(cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)*w1**7*         &
                w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)))
    
    mtx(2,4) = (-2.*a*cvd*sqrtG*(1.+ G13**2*sqrtG**2)*w1**3*w2*w4 +    &
                b*cvd*sqrtG*(1.+ G13**2*sqrtG**2)*w1**3*w2*w4 +        &
                   c*cvd*sqrtG*(1.+ G13**2*sqrtG**2)*w1**3*w2*w4 -     &
                   b*G13*                                              &
                 Sqrt(cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)*w1**7* &
                   w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)) +            &
                   c*G13*                                              &
                 Sqrt(cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)*w1**7* &
                   w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)))/            &
                (2.*cvd*sqrtG*(1.+ G13**2*sqrtG**2)*w1**3*w4**2)
    
    mtx(3,1) = -(((G13*sqrtG*w2 + w3)*                                         &
                 (-2.*a*Sqrt(cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)*        &
                          w1**7*w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)) +       &
                    b*((-cvd)*sqrtG**2*(1.+ G13**2*sqrtG**2)*w1**3*w3*w4 +     &
                         Sqrt(cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)*w1**7* &
                             w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd))) +         &
                    c*(cvd*sqrtG**2*(1.+ G13**2*sqrtG**2)*w1**3*w3*w4 +        &
                         Sqrt(cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)*w1**7* &
                             w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)))))/        &
              (2.*(w1 + G13**2*sqrtG**2*w1)*                                   &
                 Sqrt(cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)*w1**7*w4**2*   &
                     ((Rd*w4)/(p0*sqrtG))**(cpd/cvd))))
    
    mtx(3,2) = (G13*sqrtG*(-2.*a*Sqrt(cpd*cvd*p0*sqrtG**5*                        &
                         (1.+ G13**2*sqrtG**2)*w1**7*w4**2*((Rd*w4)/(p0*sqrtG))** &
                           (cpd/cvd)) +                                           &
                   b*((-cvd)*sqrtG**2*(1.+ G13**2*sqrtG**2)*w1**3*w3*w4 +         &
                        Sqrt(cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)*w1**7*     &
                            w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd))) +             &
                   c*(cvd*sqrtG**2*(1.+ G13**2*sqrtG**2)*w1**3*w3*w4 +            &
                        Sqrt(cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)*w1**7*     &
                            w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)))))/            &
             (2.*(1.+ G13**2*sqrtG**2)*Sqrt(cpd*cvd*p0*sqrtG**5*                  &
                    (1.+ G13**2*sqrtG**2)*w1**7*w4**2*((Rd*w4)/(p0*sqrtG))**      &
                      (cpd/cvd)))
    
    mtx(3,3) = (2.*a*G13**2*sqrtG**2*Sqrt(cpd*cvd*p0*sqrtG**5*                 &
                   (1.+ G13**2*sqrtG**2)*w1**7*w4**2*((Rd*w4)/(p0*sqrtG))**    &
                     (cpd/cvd)) + b*((-cvd)*sqrtG**2*(1.+ G13**2*sqrtG**2)*    &
                    w1**3*w3*w4 + Sqrt(cpd*cvd*p0*sqrtG**5*                    &
                      (1.+ G13**2*sqrtG**2)*w1**7*w4**2*((Rd*w4)/(p0*sqrtG))** &
                        (cpd/cvd))) + c*(cvd*sqrtG**2*(1.+ G13**2*sqrtG**2)*   &
                    w1**3*w3*w4 + Sqrt(cpd*cvd*p0*sqrtG**5*                    &
                      (1.+ G13**2*sqrtG**2)*w1**7*w4**2*((Rd*w4)/(p0*sqrtG))** &
                        (cpd/cvd))))/(2.*(1.+ G13**2*sqrtG**2)*                &
             Sqrt(cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)*w1**7*w4**2*       &
                 ((Rd*w4)/(p0*sqrtG))**(cpd/cvd)))
    
    mtx(3,4) = (-2.*a*w3*w4 + b*w3*w4 + c*w3*w4 -                        &
               (b*Sqrt(cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)*w1**7*  &
                        w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)))/         &
                 (sqrtG**2*(cvd + cvd*G13**2*sqrtG**2)*w1**3) +          &
               (c*Sqrt(cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)*w1**7*  &
                        w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)))/         &
                 (sqrtG**2*(cvd + cvd*G13**2*sqrtG**2)*w1**3))/(2.*w4**2)
    
    mtx(4,1) = ((b - c)*cvd*sqrtG**2*w1**2*(G13*sqrtG*w2 + w3)*w4**2)/   &
               (2.*Sqrt(cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)*w1**7* &
                      w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)))
    
    mtx(4,2) = -(((b - c)*cvd*G13*sqrtG**3*w1**3*w4**2)/                 &
               (2.*Sqrt(cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)*w1**7* &
                      w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd))))
          
    mtx(4,3) = -(((b - c)*cvd*sqrtG**2*w1**3*w4**2)/                    &
              (2.*Sqrt(cpd*cvd*p0*sqrtG**5*(1.+ G13**2*sqrtG**2)*w1**7* &
                     w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd))))
    
    mtx(4,4) = ( b + c ) / 2.
    
    calc_eigen_matrix_z = mtx
  end function calc_eigen_matrix_z
    
end module spatial_operators_mod

