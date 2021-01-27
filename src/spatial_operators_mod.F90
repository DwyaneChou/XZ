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
  
  real   (r_kind) :: dV
  
  real   (r_kind), dimension(:,:,:  ), allocatable :: qC
  real   (r_kind), dimension(:,:,:,:), allocatable :: qL
  real   (r_kind), dimension(:,:,:,:), allocatable :: qR
  real   (r_kind), dimension(:,:,:,:), allocatable :: qB
  real   (r_kind), dimension(:,:,:,:), allocatable :: qT

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

  real(r_kind), dimension(:,:,:), allocatable :: rhoL_ref ! reference density
  real(r_kind), dimension(:,:,:), allocatable :: rhoR_ref ! reference density
  real(r_kind), dimension(:,:,:), allocatable :: rhoB_ref ! reference density
  real(r_kind), dimension(:,:,:), allocatable :: rhoT_ref ! reference density
      
  real(r_kind), dimension(:,:,:), allocatable :: cL_ref ! sound speed
  real(r_kind), dimension(:,:,:), allocatable :: cR_ref ! sound speed
  real(r_kind), dimension(:,:,:), allocatable :: cB_ref ! sound speed
  real(r_kind), dimension(:,:,:), allocatable :: cT_ref ! sound speed
  
  real(r_kind), dimension(  :,:), allocatable :: PC_ref
  real(r_kind), dimension(:,:,:), allocatable :: PL_ref! Reconstructed P_ref_(i-1/2,k)
  real(r_kind), dimension(:,:,:), allocatable :: PR_ref! Reconstructed P_ref_(i+1/2,k)
  real(r_kind), dimension(:,:,:), allocatable :: PB_ref! Reconstructed P_ref_(i,k-1/2)
  real(r_kind), dimension(:,:,:), allocatable :: PT_ref! Reconstructed P_ref_(i,k+1/2)
  real(r_kind), dimension(:,:,:), allocatable :: Pz_ref! Reference pressure on bottom and top boundary
  
  real(r_kind), dimension(:,:), allocatable :: relax_coef ! Relax coefficient of Rayleigh damping
  
  real(r_kind), dimension(:,:,:), allocatable :: q_diff ! u wind, for viscosity terms only
  
  real(r_kind), dimension(:,:,:,:), allocatable :: CG_coef
  
contains
  subroutine init_spatial_operator
    integer(i_kind) :: i,j,k,iR,kR,iVar,iPOE
    integer(i_kind) :: iRec,kRec
    
    real(r_kind) :: uB_ref,uT_ref
    
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
    
    allocate(qC(nVar,              ics:ice,kcs:kce))
    allocate(qL(nVar,nPointsOnEdge,ics:ice,kcs:kce))
    allocate(qR(nVar,nPointsOnEdge,ics:ice,kcs:kce))
    allocate(qB(nVar,nPointsOnEdge,ics:ice,kcs:kce))
    allocate(qT(nVar,nPointsOnEdge,ics:ice,kcs:kce))
  
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
    
    allocate(rhoL_ref(     nPointsOnEdge,ics:ice,kcs:kce)) ! reference density
    allocate(rhoR_ref(     nPointsOnEdge,ics:ice,kcs:kce)) ! reference density
    allocate(rhoB_ref(     nPointsOnEdge,ics:ice,kcs:kce)) ! reference density
    allocate(rhoT_ref(     nPointsOnEdge,ics:ice,kcs:kce)) ! reference density
    
    allocate(cL_ref  (     nPointsOnEdge,ics:ice,kcs:kce)) ! sound speed
    allocate(cR_ref  (     nPointsOnEdge,ics:ice,kcs:kce)) ! sound speed
    allocate(cB_ref  (     nPointsOnEdge,ics:ice,kcs:kce)) ! sound speed
    allocate(cT_ref  (     nPointsOnEdge,ics:ice,kcs:kce)) ! sound speed
    
    allocate(PC_ref(              ics:ice,kcs:kce))
    allocate(PL_ref(nPointsOnEdge,ics:ice,kcs:kce))
    allocate(PR_ref(nPointsOnEdge,ics:ice,kcs:kce))
    allocate(PB_ref(nPointsOnEdge,ics:ice,kcs:kce))
    allocate(PT_ref(nPointsOnEdge,ics:ice,kcs:kce))
    
    allocate(Pz_ref(nPointsOnEdge,ids:ide,kds:kde+1))
    
    allocate(relax_coef(ics:ice,kcs:kce))
    
    allocate(q_diff(nVar,ics:ice,kcs:kce))
    
    allocate(CG_coef(nVar,maxRecTerms,ids:ide,kds:kde))
  
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
        do kRec = 1,stencil_width
          do iRec = 1,stencil_width
            if(inDomain(i-recBdy+iRec-1,k-recBdy+kRec-1))then
              j = j + 1
              iRecCell(j,i,k) = i-recBdy+iRec-1
              kRecCell(j,i,k) = k-recBdy+kRec-1
            endif
          enddo
        enddo
        nRecCells(i,k) = j
          
        if( nRecCells(i,k) >= 3           ) locPolyDegree(i,k) = 1
        if( nRecCells(i,k) >= 9           ) locPolyDegree(i,k) = 2
        if( nRecCells(i,k) >= 16          ) locPolyDegree(i,k) = 3
        if( nRecCells(i,k) >= 25          ) locPolyDegree(i,k) = 4
        if( nRecCells(i,k) == maxRecCells ) locPolyDegree(i,k) = recPolyDegree
      enddo
    enddo
    
    !do i = ibs,ibe
    !  ! special treatment on low boundary
    !  k = kds
    !  iRecCell(11,i,k) = i - 1
    !  kRecCell(11,i,k) = k + 2
    !  iRecCell(12,i,k) = i
    !  kRecCell(12,i,k) = k + 2
    !  iRecCell(13,i,k) = i + 1
    !  kRecCell(13,i,k) = k + 2
    !  iRecCell(14,i,k) = i
    !  kRecCell(14,i,k) = k + 3
    !  
    !  nRecCells    (i,k) = 14
    !  locPolyDegree(i,k) = 3
    !  
    !  k = kds + 1
    !  iRecCell(16,i,k) = i - 1
    !  kRecCell(16,i,k) = k + 2
    !  iRecCell(17,i,k) = i
    !  kRecCell(17,i,k) = k + 2
    !  iRecCell(18,i,k) = i + 1
    !  kRecCell(18,i,k) = k + 2
    !  iRecCell(19,i,k) = i - 1
    !  kRecCell(19,i,k) = k + 3
    !  iRecCell(20,i,k) = i
    !  kRecCell(20,i,k) = k + 3
    !  iRecCell(21,i,k) = i + 1
    !  kRecCell(21,i,k) = k + 3
    !  
    !  nRecCells    (i,k) = 21
    !  locPolyDegree(i,k) = 3
    !enddo
    
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
    !
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
    !
    !do i = ibs,ibe
    !  ! special treatment on low boundary
    !  k = kds
    !  iRecCell(11,i,k) = i - 1
    !  kRecCell(11,i,k) = k + 2
    !  iRecCell(12,i,k) = i
    !  kRecCell(12,i,k) = k + 2
    !  iRecCell(13,i,k) = i + 1
    !  kRecCell(13,i,k) = k + 2
    !  iRecCell(14,i,k) = i
    !  kRecCell(14,i,k) = k + 3
    !  
    !  nRecCells    (i,k) = 14
    !  locPolyDegree(i,k) = 3
    !  
    !  k = kds + 1
    !  iRecCell(16,i,k) = i - 1
    !  kRecCell(16,i,k) = k + 2
    !  iRecCell(17,i,k) = i
    !  kRecCell(17,i,k) = k + 2
    !  iRecCell(18,i,k) = i + 1
    !  kRecCell(18,i,k) = k + 2
    !  iRecCell(19,i,k) = i - 1
    !  kRecCell(19,i,k) = k + 3
    !  iRecCell(20,i,k) = i
    !  kRecCell(20,i,k) = k + 3
    !  iRecCell(21,i,k) = i + 1
    !  kRecCell(21,i,k) = k + 3
    !  
    !  nRecCells    (i,k) = 21
    !  locPolyDegree(i,k) = 3
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
    
    !$OMP PARALLEL DO PRIVATE(i,j,iRec,kRec)
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
      enddo
    enddo
    !$OMP END PARALLEL DO
    
    !! Reconstruct metric function
    !call reconstruct_field(sqrtGL,&
    !                       sqrtGR,&
    !                       sqrtGB,&
    !                       sqrtGT,&
    !                       sqrtGC*recdV)
    !
    !call reconstruct_field(G13L,&
    !                       G13R,&
    !                       G13B,&
    !                       G13T,&
    !                       G13C*recdV)
    
        
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
      call reconstruct_field(qL(iVar,:,:,:),&
                             qR(iVar,:,:,:),&
                             qB(iVar,:,:,:),&
                             qT(iVar,:,:,:),&
                             qC(iVar  ,:,:)*recdV)
    enddo
    
    qL_ref = qL
    qR_ref = qR
    qB_ref = qB
    qT_ref = qT
  
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
    
    !$OMP PARALLEL DO PRIVATE(i,iPOE)
    do k = kds-1,kde+1
      do i = ids-1,ide+1
        do iPOE = 1,nPointsOnEdge
          rhoL_ref(iPOE,i,k) = ( qL_ref(1,iPOE,i,k) + qL_ref(5,iPOE,i,k) ) / sqrtGL(iPOE,i,k)
          rhoR_ref(iPOE,i,k) = ( qR_ref(1,iPOE,i,k) + qR_ref(5,iPOE,i,k) ) / sqrtGR(iPOE,i,k)
          rhoB_ref(iPOE,i,k) = ( qB_ref(1,iPOE,i,k) + qB_ref(5,iPOE,i,k) ) / sqrtGB(iPOE,i,k)
          rhoT_ref(iPOE,i,k) = ( qT_ref(1,iPOE,i,k) + qT_ref(5,iPOE,i,k) ) / sqrtGT(iPOE,i,k)
          
          cL_ref(iPOE,i,k) = calc_sound_speed_x(sqrtGL(iPOE,i,k)               ,qL_ref(:,iPOE,i,k))
          cR_ref(iPOE,i,k) = calc_sound_speed_x(sqrtGR(iPOE,i,k)               ,qR_ref(:,iPOE,i,k))
          cB_ref(iPOE,i,k) = calc_sound_speed_z(sqrtGB(iPOE,i,k),G13B(iPOE,i,k),qB_ref(:,iPOE,i,k)) * sqrtGB(iPOE,i,k) ! unit: m/s
          cT_ref(iPOE,i,k) = calc_sound_speed_z(sqrtGT(iPOE,i,k),G13T(iPOE,i,k),qT_ref(:,iPOE,i,k)) * sqrtGT(iPOE,i,k) ! unit: m/s
      
          PL_ref(iPOE,i,k) = calc_pressure(sqrtGL(iPOE,i,k),qL_ref(:,iPOE,i,k))
          PR_ref(iPOE,i,k) = calc_pressure(sqrtGR(iPOE,i,k),qR_ref(:,iPOE,i,k))
          PB_ref(iPOE,i,k) = calc_pressure(sqrtGB(iPOE,i,k),qB_ref(:,iPOE,i,k))
          PT_ref(iPOE,i,k) = calc_pressure(sqrtGT(iPOE,i,k),qT_ref(:,iPOE,i,k))
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    
    do k = kds,kde+1
      do i = ids,ide
        do iPOE = 1,nPointsOnEdge
          uB_ref = ( qB_ref(3,iPOE,i,k  ) + sqrtGB(iPOE,i,k  ) * G13B(iPOE,i,k  ) * qB_ref(2,iPOE,i,k  ) ) / ( qB_ref(1,iPOE,i,k  ) + qB_ref(5,iPOE,i,k  ) )
          uT_ref = ( qT_ref(3,iPOE,i,k-1) + sqrtGT(iPOE,i,k-1) * G13T(iPOE,i,k-1) * qT_ref(2,iPOE,i,k-1) ) / ( qT_ref(1,iPOE,i,k-1) + qT_ref(5,iPOE,i,k-1) )
          call AUSM_p(Pz_ref(iPOE,i,k), rhoT_ref(iPOE,i,k-1), rhoB_ref(iPOE,i,k),&
                                        uT_ref              , uB_ref            ,&
                                        cT_ref  (iPOE,i,k-1), cB_ref  (iPOE,i,k),&
                                        pT_ref  (iPOE,i,k-1), pB_ref  (iPOE,i,k))
        enddo
      enddo
    enddo
    
    ! Calculate Rayleigh damping coef
    call Rayleigh_coef(relax_coef)
    
    ! Fill out qC
    qC = FillValue
    
  end subroutine init_spatial_operator
  
  subroutine spatial_operator(stat,tend)
    type(stat_field), target, intent(inout) :: stat
    type(tend_field), target, intent(inout) :: tend
  
    real(r_kind) maxeigen_x
    real(r_kind) maxeigen_z
    
    integer(i_kind) :: i,j,k,iR,kR,iVar,iPOE
  
    integer(i_kind) ip1,im1
    integer(i_kind) kp1,km1
    integer(i_kind) ip2,im2
    integer(i_kind) kp2,km2
  
    ! Attension stat is changed here!
    if(case_num==2)call Rayleigh_damping(stat%q,ref%q)
    
    ! copy stat
    !qC(:,ids:ide,kds:kde) = stat%q(:,ids:ide,kds:kde)
    !$OMP PARALLEL DO PRIVATE(i,iVar)
    do k = kds,kde
      do i = ids,ide
        do iVar = 1,nVar
          qC(iVar,i,k) = stat%q(iVar,i,k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    
    do iVar = 1,nVar
      call reconstruct_field(qL(iVar,:,:,:),&
                             qR(iVar,:,:,:),&
                             qB(iVar,:,:,:),&
                             qT(iVar,:,:,:),&
                             qC(iVar  ,:,:)*recdV,&
                             iVar)
    enddo
  
    ! Boundary Condition
    ! Fill boundary
    if(case_num==1.or.case_num==3)then
      qL(2,:,ids,kds:kde) = 0
      qR(2,:,ide,kds:kde) = 0
      qB(3,:,ids:ide,kds) = 0
      qT(3,:,ids:ide,kde) = 0
    elseif(case_num==2)then
      !$OMP PARALLEL DO PRIVATE(iPOE)
      do k = kds,kde
        do iPOE = 1,nPointsOnEdge
          qL(2,iPOE,ids,k) = ref%q(2,ids,k)
          qR(2,iPOE,ide,k) = ref%q(2,ide,k)
        enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO PRIVATE(iPOE)
      do i = ids,ide
        do iPOE = 1,nPointsOnEdge
          qB(3,iPOE,i,kds) = -sqrtGB(iPOE,i,kds) * G13B(iPOE,i,kds) * qB(2,iPOE,i,kds)
          qT(3,iPOE,i,kde) = -sqrtGT(iPOE,i,kde) * G13T(iPOE,i,kde) * qT(2,iPOE,i,kde)
        enddo
      enddo
      !$OMP END PARALLEL DO
    endif
    
    ! Fill outside boundary
    !$OMP PARALLEL DO
    do k = kds,kde
      ! left boundary
      qR(:,:,ids-1,k) = qL(:,:,ids,k)
      
      ! right boundary
      qL(:,:,ide+1,k) = qR(:,:,ide,k)
    enddo
    !$OMP END PARALLEL DO
    
    !$OMP PARALLEL DO
    do i = ids,ide
      ! bottom boundary
      qT(:,:,i,kds-1) = qB(:,:,i,kds)
      
      ! top boundary
      qB(:,:,i,kde+1) = qT(:,:,i,kde)
    enddo
    !$OMP END PARALLEL DO
  
    ! Calculate pressure X
    !$OMP PARALLEL DO PRIVATE(i,iPOE)
    do k = kds,kde
      do i = ids,ide
        do iPOE = 1,nPointsOnEdge
          PL(iPOE,i,k) = calc_pressure(sqrtGL(iPOE,i,k),qL(:,iPOE,i,k))
          PR(iPOE,i,k) = calc_pressure(sqrtGR(iPOE,i,k),qR(:,iPOE,i,k))
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO PRIVATE(iPOE,i)
    do k = kds,kde
      do iPOE = 1,nPointsOnEdge
        i = ide + 1
        PL(iPOE,i,k) = calc_pressure(sqrtGL(iPOE,i,k),qL(:,iPOE,i,k))
        i = ids - 1
        PR(iPOE,i,k) = calc_pressure(sqrtGR(iPOE,i,k),qR(:,iPOE,i,k))
      enddo
    enddo
    !$OMP END PARALLEL DO
    
    ! Calculate pressure Z
    !$OMP PARALLEL DO PRIVATE(i,iPOE)
    do k = kds,kde
      do i = ids,ide
        do iPOE = 1,nPointsOnEdge
          PB(iPOE,i,k) = calc_pressure(sqrtGB(iPOE,i,k),qB(:,iPOE,i,k))
          PT(iPOE,i,k) = calc_pressure(sqrtGT(iPOE,i,k),qT(:,iPOE,i,k))
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO PRIVATE(k,iPOE)
    do i = ids,ide
      do iPOE = 1,nPointsOnEdge
        k = kde + 1
        PB(iPOE,i,k) = calc_pressure(sqrtGB(iPOE,i,k),qB(:,iPOE,i,k))
        k = kds - 1
        PT(iPOE,i,k) = calc_pressure(sqrtGT(iPOE,i,k),qT(:,iPOE,i,k))
      enddo
    enddo
    !$OMP END PARALLEL DO
    
    ! Calculate x directional flux on each points on edge
    !$OMP PARALLEL DO PRIVATE(i,im1,iPOE,iVar)
    do k = kds,kde
      do i = ids,ide+1
        do iPOE = 1,nPointsOnEdge
          im1 = i - 1
          FeP(:,iPOE,i,k) = calc_F(sqrtGR(iPOE,im1,k),sqrtGL(iPOE,i,k),qR(:,iPOE,im1,k),qL(:,iPOE,i,k),pR(iPOE,im1,k),pL(iPOE,i,k))
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    ! Calculate x directional flux on edges
    !$OMP PARALLEL DO PRIVATE(i,iVar)
    do k = kds,kde
      do i = ids,ide+1
        do iVar = 1,nVar
          Fe(iVar,i,k) = Gaussian_quadrature_1d(FeP(iVar,:,i,k))
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    
    ! Calculate z directional flux on each points on edge
    !$OMP PARALLEL DO PRIVATE(k,km1,iPOE,iVar)
    do i = ids,ide
      do k = kds,kde+1
        do iPOE = 1,nPointsOnEdge
          km1 = k - 1
          HeP(:,iPOE,i,k) = calc_H(sqrtGT(  iPOE,i,km1),sqrtGB(  iPOE,i,k),G13T(iPOE,i,km1),G13B(iPOE,i,k),&
                                   qT    (:,iPOE,i,km1),qB    (:,iPOE,i,k),pT  (iPOE,i,km1),pB  (iPOE,i,k),&
                                   Pz_ref(  iPOE,i,k))
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    ! Calculate z directional flux on edges
    !$OMP PARALLEL DO PRIVATE(k,iVar)
    do i = ids,ide
      do k = kds,kde+1
        do iVar = 1,nVar
          He(iVar,i,k) = Gaussian_quadrature_1d(HeP(iVar,:,i,k))
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    
    !$OMP PARALLEL DO PRIVATE(i)
    do k = kds,kde
      do i = ids,ide
        rho_p(i,k) = ( qC(1,i,k) + qC(5,i,k) - ref%q(1,i,k) - ref%q(5,i,k) ) / sqrtGC(i,k)  ! hydrostatic
        !rho_p(i,k) = ( qC(1,i,k) + qC(5,i,k) ) / sqrtGC(i,k)  ! nonhydrostatic
      enddo
    enddo
    !$OMP END PARALLEL DO
    where(abs(rho_p)<=1.E-13)rho_p=0.! hydrostatic
    
    !$OMP PARALLEL DO PRIVATE(i,iVar)
    do k = kds,kde
      do i = ids,ide
        do iVar = 1,nVar
          src(iVar,i,k) = 0
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO PRIVATE(i)
    do k = kds,kde
      do i = ids,ide
        src(3,i,k) = src(3,i,k) - sqrtGC(i,k) * rho_p(i,k) * gravity
      enddo
    enddo
    !$OMP END PARALLEL DO
  
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
    !$OMP PARALLEL DO PRIVATE(i,ip1,kp1,iVar)
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

  subroutine reconstruct_field(qL,qR,qB,qT,qC,iVar)
    real(r_kind), dimension(nPointsOnEdge,ics:ice,kcs:kce),intent(out) :: qL
    real(r_kind), dimension(nPointsOnEdge,ics:ice,kcs:kce),intent(out) :: qR
    real(r_kind), dimension(nPointsOnEdge,ics:ice,kcs:kce),intent(out) :: qB
    real(r_kind), dimension(nPointsOnEdge,ics:ice,kcs:kce),intent(out) :: qT
    real(r_kind), dimension(              ics:ice,kcs:kce),intent(in ) :: qC
    
    integer(i_kind), optional :: iVar ! Specify variable for speed up CG method
  
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
    !$OMP PARALLEL DO PRIVATE(i,j,m,n,iRec,kRec,iC,u,A,polyCoef)
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
        
        
        if(present(iVar))then
          polyCoef(1:n) = WLS_ENO(A(1:m,1:n),u(1:m),h,m,n,ic,CG_coef(iVar,1:n,i,k))
          CG_coef(iVar,1:n,i,k) = polyCoef(1:n)
        else
          polyCoef(1:n) = WLS_ENO(A(1:m,1:n),u(1:m),h,m,n,ic)
        endif
        
        qL(:,i,k) = matmul(recMatrixL(:,1:n,i,k),polyCoef(1:n))
        qR(:,i,k) = matmul(recMatrixR(:,1:n,i,k),polyCoef(1:n))
        qB(:,i,k) = matmul(recMatrixB(:,1:n,i,k),polyCoef(1:n))
        qT(:,i,k) = matmul(recMatrixT(:,1:n,i,k),polyCoef(1:n))
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
  
  function calc_F(sqrtGL,sqrtGR,qL,qR,pL,pR)
    real(r_kind) :: calc_F(nVar)
    real(r_kind) :: sqrtGL
    real(r_kind) :: sqrtGR
    real(r_kind) :: qL(nVar)
    real(r_kind) :: qR(nVar)
    real(r_kind) :: pL      ! pressure
    real(r_kind) :: pR      ! pressure
    
    real(r_kind) sqrtGrhoL
    real(r_kind) sqrtGrhoR
    
    real(r_kind) rhoL
    real(r_kind) rhoR
    
    real(r_kind) uL
    real(r_kind) uR
    
    real(r_kind) cL
    real(r_kind) cR
    
    real(r_kind) m
    real(r_kind) p ! sqrtG * p
    
    sqrtGrhoL = qL(1) + qL(5)
    sqrtGrhoR = qR(1) + qR(5)
    
    rhoL = sqrtGrhoL / sqrtGL
    rhoR = sqrtGrhoR / sqrtGR
    uL   = qL(2) / sqrtGrhoL
    uR   = qR(2) / sqrtGrhoR
    cL   = calc_sound_speed_x(sqrtGL,qL)
    cR   = calc_sound_speed_x(sqrtGR,qR)
    
    call AUSM_up_x(m,p,sqrtGL,sqrtGR,rhoL,rhoR,uL,uR,pL,pR,cL,cR)
    
    calc_F = 0.5 * m * ( qL + qR - sign(1._r_kind,m) * ( qR - qL ) )
    calc_F(2) = calc_F(2) + p
    
  end function calc_F
  
  function calc_H(sqrtGL,sqrtGR,G13L,G13R,qL,qR,pL,pR,p_ref)
    real(r_kind) :: calc_H(nVar)
    real(r_kind) :: sqrtGL
    real(r_kind) :: sqrtGR
    real(r_kind) :: G13L
    real(r_kind) :: G13R
    real(r_kind) :: qL(nVar)
    real(r_kind) :: qR(nVar)
    real(r_kind) :: pL        ! pressure
    real(r_kind) :: pR        ! pressure
    real(r_kind) :: p_ref     ! reference pressure
    
    real(r_kind) sqrtGrhoL
    real(r_kind) sqrtGrhoR
    
    real(r_kind) rhoL
    real(r_kind) rhoR
    
    real(r_kind) uL
    real(r_kind) uR
    
    real(r_kind) cL
    real(r_kind) cR
    
    real(r_kind) m
    real(r_kind) p      ! sqrtG * G13 * p
    real(r_kind) p_pert ! p - p_ref
    
    sqrtGrhoL = qL(1) + qL(5)
    sqrtGrhoR = qR(1) + qR(5)
    
    rhoL = sqrtGrhoL / sqrtGL
    rhoR = sqrtGrhoR / sqrtGR
    uL   = ( qL(3) + sqrtGL * G13L * qL(2) ) / sqrtGrhoL ! a sqrtG has been multipled here
    uR   = ( qR(3) + sqrtGR * G13R * qR(2) ) / sqrtGrhoR ! a sqrtG has been multipled here
    cL   = calc_sound_speed_z(sqrtGL,G13L,qL) * sqrtGL
    cR   = calc_sound_speed_z(sqrtGR,G13R,qR) * sqrtGR
    
    call AUSM_up_z(m, p, p_pert,                                                  &
                   sqrtGL, sqrtGR, G13L, G13R, rhoL, rhoR, uL, uR, cL, cR, pL, pR,&
                   p_ref)
    
    calc_H = 0.5 * m * ( qL + qR - sign(1._r_kind,m) * ( qR - qL ) )
    calc_H(2) = calc_H(2) + p
    calc_H(3) = calc_H(3) + p_pert
    
  end function calc_H
      
  function calc_sound_speed_x(sqrtG,q)
    real(r_kind) :: calc_sound_speed_x
    real(r_kind) :: sqrtG
    real(r_kind) :: q(nVar)
    
    real(r_kind) w1
    real(r_kind) w2
    real(r_kind) w3
    real(r_kind) w4
    real(r_kind) w5
    
    real(r_kind) coef1,coef2

    real(r_kind) c1,c2,c3,c4,c5
    
    if(any(q==FillValue).or.sqrtG==FillValue)then
      calc_sound_speed_x = 0
    else
      w1 = q(1)
      w2 = q(2)
      w3 = q(3)
      w4 = q(4)
      w5 = q(5)

      c1 = w1 + w5
      c2 = cpd*w1 + cpv*w5
      c3 = cvd*w1 + cvv*w5
      c4 = c1 + eq*w5
      c5 = p0*sqrtG
      
      coef1 = sqrt( c5*w1**2*w4**2*c1**3*c2*&
      c3**3*c4**2*((Rd*w4*c4)/(c5*w1))**&
      (c2/c3) )
      
      coef2 = w1*w4*c1**2*c3**2*c4
      
      calc_sound_speed_x = coef1 / coef2
    endif
    
  end function calc_sound_speed_x
  
  function calc_sound_speed_z(sqrtG,G13,q)
    real(r_kind) :: calc_sound_speed_z
    real(r_kind) :: sqrtG
    real(r_kind) :: G13
    real(r_kind) :: q(nVar)
    
    real(r_kind) w1
    real(r_kind) w2
    real(r_kind) w3
    real(r_kind) w4
    real(r_kind) w5
    
    real(r_kind) coef1,coef2

    real(r_kind) c1,c2,c3,c4,c5
    
    if(any(q==FillValue).or.sqrtG==FillValue)then
      calc_sound_speed_z = 0
    else
      w1 = q(1)
      w2 = q(2)
      w3 = q(3)
      w4 = q(4)
      w5 = q(5)

      c1 = w1 + w5
      c2 = cpd*w1 + cpv*w5
      c3 = cvd*w1 + cvv*w5
      c4 = c1 + eq*w5
      c5 = p0*sqrtG
      
      coef1 = sqrt( c5*(1. + G13**2*sqrtG**2)*w1**2* &
            w4**2*c1**3*c2*c3**3*                    &
            c4**2*((Rd*w4*c4)/(c5*w1))**             &
            (c2/c3) )
      
      coef2 = sqrtG*w1*w4*c1**2*c3**2*c4
      
      calc_sound_speed_z = coef1 / coef2
    endif
    
  end function calc_sound_speed_z
  
  subroutine AUSM_up_x(m,p,sqrtGL,sqrtGR,rhoL,rhoR,uL,uR,pL,pR,cL,cR)
    real(r_kind),intent(out) :: m
    real(r_kind),intent(out) :: p
    real(r_kind),intent(in ) :: sqrtGL
    real(r_kind),intent(in ) :: sqrtGR
    real(r_kind),intent(in ) :: rhoL
    real(r_kind),intent(in ) :: rhoR
    real(r_kind),intent(in ) :: uL
    real(r_kind),intent(in ) :: uR
    real(r_kind),intent(in ) :: pL
    real(r_kind),intent(in ) :: pR
    real(r_kind),intent(in ) :: cL ! Sound speed
    real(r_kind),intent(in ) :: cR ! Sound speed
    
    real(r_kind),parameter :: Ku    = 0.75
    real(r_kind),parameter :: Kp    = 0.25
    real(r_kind),parameter :: sigma = 1.
    real(r_kind),parameter :: sp    = 1.
    real(r_kind),parameter :: sn    = -1.
    
    real(r_kind) :: rho
    real(r_kind) :: a
    real(r_kind) :: ML
    real(r_kind) :: MR
    real(r_kind) :: Mbar2
    real(r_kind) :: Mh
    
    real(r_kind) :: P5MLsp
    real(r_kind) :: P5MRsn
    
    rho = 0.5 * ( rhoL + rhoR )
    a   = 0.5 * ( cL + cR )
    
    ML = uL / a
    MR = uR / a
    
    Mbar2 = ( uL**2 + uR**2 ) / ( 2. * a**2 )
    
    Mh = M4( ML, sp ) + M4( MR, sn ) - Kp * max( 1. - sigma * Mbar2, 0. ) * ( PR - PL ) / ( rho * a**2 )
    m  = a * Mh
    
    P5MLsp = P5(ML,sp)
    P5MRsn = P5(MR,sn)
    
    p = P5MLsp * sqrtGL * PL + P5MRsn * sqrtGR * PR - Ku * P5MLsp * P5MRsn * ( sqrtGL * rhoL + sqrtGR * rhoR ) * a * ( uR - uL )
    
  end subroutine AUSM_up_x
  
  subroutine AUSM_up_z(m, p, p_pert,                                                  & ! output
                       sqrtGL, sqrtGR, G13L, G13R, rhoL, rhoR, uL, uR, cL, cR, pL, pR,& ! input state
                                                                                 P_ref) ! input reference state
    real(r_kind),intent(out) :: m
    real(r_kind),intent(out) :: p        ! Actully sqrtG * G13 * p
    real(r_kind),intent(out) :: p_pert   ! Actully p - p_ref
    real(r_kind),intent(in ) :: sqrtGL
    real(r_kind),intent(in ) :: sqrtGR
    real(r_kind),intent(in ) :: G13L
    real(r_kind),intent(in ) :: G13R
    real(r_kind),intent(in ) :: rhoL
    real(r_kind),intent(in ) :: rhoR
    real(r_kind),intent(in ) :: uL
    real(r_kind),intent(in ) :: uR
    real(r_kind),intent(in ) :: cL       ! Sound speed
    real(r_kind),intent(in ) :: cR       ! Sound speed
    real(r_kind),intent(in ) :: pL
    real(r_kind),intent(in ) :: pR
    real(r_kind),intent(in ) :: p_ref   ! reference pressure (hydrostatic pressure)
    
    real(r_kind),parameter :: Ku    = 0.75
    real(r_kind),parameter :: Kp    = 0.25
    real(r_kind),parameter :: sigma = 1.
    real(r_kind),parameter :: sp    = 1.
    real(r_kind),parameter :: sn    = -1.
    
    real(r_kind) :: rho
    real(r_kind) :: a
    real(r_kind) :: ML
    real(r_kind) :: MR
    real(r_kind) :: Mbar2
    real(r_kind) :: Mh
    
    real(r_kind) :: p_diff
    
    real(r_kind) :: coefL
    real(r_kind) :: coefR
    
    real(r_kind) :: P5MLsp
    real(r_kind) :: P5MRsn
    
    rho = 0.5 * ( rhoL + rhoR )
    a   = 0.5 * ( cL + cR )
    
    ML = uL / a
    MR = uR / a
    
    Mbar2 = ( uL**2 + uR**2 ) / ( 2. * a**2 )
    
    P5MLsp = P5(ML,sp)
    P5MRsn = P5(MR,sn)
    
    p_diff = PR - PL
    
    Mh     = M4( ML, sp ) + M4( MR, sn ) - Kp * max( 1. - sigma * Mbar2, 0. ) * p_diff / ( rho * a**2 )
    m      = 0.5 * ( cL / sqrtGL + cR / sqrtGR ) * Mh
    
    coefL = sqrtGL * G13L
    coefR = sqrtGR * G13R
    
    p = P5MLsp * coefL * PL + P5MRsn * coefR * PR &
      - Ku * P5MLsp * P5MRsn * ( coefL * rhoL + coefR * rhoR ) * a * ( uR - uL )
    
    p_pert = P5MLsp * pL + P5MRsn * pR - 2. * Ku * P5MLsp * P5MRsn * ( rho * a ) * ( uR - uL ) ! pressure
    p_pert = p_pert - p_ref ! pressure perturbation, hydrostatic
    
    if( abs(p_pert/p_ref) < 1.e-14 )p_pert = 0! hydrostatic
    
  end subroutine AUSM_up_z

  ! AUSM calculate pressure only
  subroutine AUSM_p(p, rhoL, rhoR, uL, uR, cL, cR, pL, pR)
    real(r_kind),intent(out) :: p
    real(r_kind),intent(in ) :: rhoL
    real(r_kind),intent(in ) :: rhoR
    real(r_kind),intent(in ) :: uL
    real(r_kind),intent(in ) :: uR
    real(r_kind),intent(in ) :: cL       ! Sound speed
    real(r_kind),intent(in ) :: cR       ! Sound speed
    real(r_kind),intent(in ) :: pL
    real(r_kind),intent(in ) :: pR
    
    real(r_kind),parameter :: Ku    = 0.75
    real(r_kind),parameter :: Kp    = 0.25
    real(r_kind),parameter :: sp    = 1.
    real(r_kind),parameter :: sn    = -1.
    
    real(r_kind) :: rho
    real(r_kind) :: a
    real(r_kind) :: ML
    real(r_kind) :: MR
    
    real(r_kind) :: P5MLsp
    real(r_kind) :: P5MRsn
    
    rho = 0.5 * ( rhoL + rhoR )
    a   = 0.5 * ( cL + cR )
    
    ML = uL / a
    MR = uR / a
    
    P5MLsp = P5(ML,sp)
    P5MRsn = P5(MR,sn)
    
    p = P5MLsp * PL + P5MRsn * PR &
      - Ku * P5MLsp * P5MRsn * ( rhoL + rhoR ) * a * ( uR - uL )
    
  end subroutine AUSM_p
  
  function M2(M,signal)
    real(r_kind) :: M2
    real(r_kind) :: M
    real(r_kind) :: signal ! must be 1 or -1
    
    M2 = signal * 0.25 * ( M + signal )**2
    
  end function M2
  
  function M4(M,signal)
    real(r_kind) :: M4
    real(r_kind) :: M
    real(r_kind) :: signal ! must be 1 or -1
    
    real(r_kind),parameter :: beta  = 0.125
    
    if(abs(M)>=1)then
      M4 = 0.5 * ( M + signal * abs(M) )
    else
      M4 = M2( M, signal ) * ( 1. - signal * 16. * beta * M2( M, -signal ) )
    endif
    
  end function M4
  
  function P5(M,signal)
    real(r_kind) :: P5
    real(r_kind) :: M
    real(r_kind) :: signal ! must be 1 or -1
    
    real(r_kind),parameter :: alpha = 0.1875
    
    if(abs(M)>=1)then
      P5 = 0.5 * ( 1. + signal * sign(1._r_kind,M) )
    else
      P5 = M2( M, signal ) * ( ( 2. * signal - M ) - signal * 16.*alpha * M * M2( M, -signal ) )
    endif
    
  end function P5
  
  subroutine Rayleigh_damping(q,q_ref)
    real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(inout) :: q
    real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(in   ) :: q_ref
    
    integer(i_kind), parameter :: vs = 1
    integer(i_kind), parameter :: ve = 5
    
    integer i,k,iVar
    
    !$OMP PARALLEL DO PRIVATE(i,iVar)
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
    
    real(r_kind), parameter :: mu_max_top = 0.2
    real(r_kind), parameter :: mu_max_lat = 0.2
    
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
end module spatial_operators_mod

