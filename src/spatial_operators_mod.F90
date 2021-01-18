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
      
      integer(i_kind), dimension(:,:), allocatable :: iCenCell ! center cell index on reconstruction stencil
      
      integer(i_kind), dimension(:,:,:,:), allocatable :: iRecCell ! x index of reconstruction cells
      integer(i_kind), dimension(:,:,:,:), allocatable :: kRecCell ! k index of reconstruction cells
      
      real   (r_kind), dimension(:,:,:,:,:), allocatable :: xRel   ! relative x coordinate of reconstruction cells
      real   (r_kind), dimension(:,:,:,:,:), allocatable :: etaRel ! relative eta coordinate of reconstruction cells
      
      real   (r_kind), dimension(:,:,:), allocatable :: disCenter ! distance between adjacent cells and center cell on stencil
      
      real   (r_kind), dimension(:,:), allocatable :: recMatrixL
      real   (r_kind), dimension(:,:), allocatable :: recMatrixR
      real   (r_kind), dimension(:,:), allocatable :: recMatrixB
      real   (r_kind), dimension(:,:), allocatable :: recMatrixT
      
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
      
      real(r_kind), dimension(:,:), allocatable :: relax_coef ! Relax coefficient of Rayleigh damping
      
      real(r_kind), dimension(:,:,:), allocatable :: q_diff ! u wind, for viscosity terms only
      
    contains
      subroutine init_spatial_operator
        integer(i_kind) :: i,j,k,iR,kR,iVar,iPOE
        integer(i_kind) :: iRec,kRec
        
        allocate(iCenCell  (ids:ide,kds:kde))
        
        allocate(iRecCell  (stencil_width,stencil_width,ids:ide,kds:kde))
        allocate(kRecCell  (stencil_width,stencil_width,ids:ide,kds:kde))
        
        allocate(xRel  (4,stencil_width,stencil_width,ids:ide,kds:kde))
        allocate(etaRel(4,stencil_width,stencil_width,ids:ide,kds:kde))
        
        allocate(disCenter(nRecCells,ids:ide,kds:kde))
        
        allocate(recMatrixL(nPointsOnEdge,nRecTerms))
        allocate(recMatrixR(nPointsOnEdge,nRecTerms))
        allocate(recMatrixB(nPointsOnEdge,nRecTerms))
        allocate(recMatrixT(nPointsOnEdge,nRecTerms))
        
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
        
        allocate(relax_coef(ics:ice,kcs:kce))
        
        allocate(q_diff(nVar,ics:ice,kcs:kce))
      
        src   = 0
        
        dV = dx * deta
        
        ! set reconstruction cells on each stencil
        print*,'Set reconstruction cells on each stencil'
        print*,''
        ! x-dir
        do k = kds,kde
          ! middle
          do i = ibs,ibe
            do kR = 1,stencil_width
              do iR = 1,stencil_width
                iRecCell(iR,kR,i,k) = iR + i - 1 - recBdy
              enddo
            enddo            
          enddo
          ! left
          do i = ids,ibs-1
            do kR = 1,stencil_width
              do iR = 1,stencil_width
                iRecCell(iR,kR,i,k) = ids - 1 + iR
              enddo
            enddo
          enddo
          ! right
          do i = ibe+1,ide
            do kR = 1,stencil_width
              do iR = 1,stencil_width
                iRecCell(iR,kR,i,k) = ide - stencil_width + iR
              enddo
            enddo
          enddo
        enddo
        
        ! z-dir
        do i = ids,ide
          ! middle
          do k = kbs,kbe
            do kR = 1,stencil_width
              do iR = 1,stencil_width
                kRecCell(iR,kR,i,k) = kR + k - 1 - recBdy
              enddo
            enddo            
          enddo
          ! bottom
          do k = kds,kbs-1
            do kR = 1,stencil_width
              do iR = 1,stencil_width
                kRecCell(iR,kR,i,k) = kds - 1 + kR
              enddo
            enddo
          enddo
          ! top
          do k = kbe+1,kde
            do kR = 1,stencil_width
              do iR = 1,stencil_width
                kRecCell(iR,kR,i,k) = kde - stencil_width  + kR
              enddo
            enddo
          enddo
        enddo
        
        ! initilize polyCoordCoef
        print*,'Initilize polyCoordCoef'
        print*,''
        
        recCoef = 1.
        recdx   = 1. / ( dx   * recCoef )
        recdeta = 1. / ( deta * recCoef )
        recdV   = 1. / ( recCoef**2 )
        
        do k = kds,kde
          do i = ids,ide
            j = 0 ! Reset cell number on stencil
            do kR = 1,stencil_width
              do iR = 1,stencil_width
                iRec = iRecCell(iR,kR,i,k)
                kRec = kRecCell(iR,kR,i,k)
                
                j = j + 1
                
                if(iRec==i.and.kRec==k)iCenCell(i,k) = j
                
                xRel  (:,iR,kR,i,k) = ( xCorner  (:,iRec,kRec) - xCenter  (i,k) ) * recdx
                etaRel(:,iR,kR,i,k) = ( etaCorner(:,iRec,kRec) - etaCenter(i,k) ) * recdeta
                
                disCenter(j,i,k) = sqrt( ( xCenter(iRec,kRec) - xCenter(i,k) )**2 + ( etaCenter(iRec,kRec) - etaCenter(i,k) )**2 )
                
                call calc_polynomial_square_integration(recPolyDegree,xRel  (1,iR,kR,i,k),xRel  (2,iR,kR,i,k),&
                                                                      etaRel(1,iR,kR,i,k),etaRel(4,iR,kR,i,k),polyCoordCoef(j,:,i,k))
              enddo
            enddo
          enddo
        enddo
        
        ! Calculate reconstruction matrix on edge
        call calc_polynomial_matrix(recPolyDegree,nPointsOnEdge,nRecTerms,xL*recdx,etaL*recdeta,recMatrixL)
        call calc_polynomial_matrix(recPolyDegree,nPointsOnEdge,nRecTerms,xR*recdx,etaR*recdeta,recMatrixR)
        call calc_polynomial_matrix(recPolyDegree,nPointsOnEdge,nRecTerms,xB*recdx,etaB*recdeta,recMatrixB)
        call calc_polynomial_matrix(recPolyDegree,nPointsOnEdge,nRecTerms,xT*recdx,etaT*recdeta,recMatrixT)
        
        ! Reconstruct metric function
        call reconstruct_field(sqrtGL,&
                               sqrtGR,&
                               sqrtGB,&
                               sqrtGT,&
                               sqrtGC*recdV)
        
        call reconstruct_field(G13L,&
                               G13R,&
                               G13B,&
                               G13T,&
                               G13C*recdV)
        
        ! Calculate reference pressure
        qC = ref%q
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
        qC(:,ids:ide,kds:kde) = stat%q(:,ids:ide,kds:kde)
        
        do iVar = 1,nVar
          call reconstruct_field(qL(iVar,:,:,:),&
                                 qR(iVar,:,:,:),&
                                 qB(iVar,:,:,:),&
                                 qT(iVar,:,:,:),&
                                 qC(iVar  ,:,:)*recdV)
        enddo
      
        ! Boundary Condition
        ! Fill boundary
        if(case_num==1.or.case_num==3)then
          qL(2,:,ids,:) = 0
          qR(2,:,ide,:) = 0
          qB(3,:,:,kds) = 0
          qT(3,:,:,kde) = 0
        elseif(case_num==2)then
          do iPOE = 1,nPointsOnEdge
            qL(2,iPOE,ids,kds:kde) = ref%q(2,ids,kds:kde)
            qR(2,iPOE,ide,kds:kde) = ref%q(2,ide,kds:kde)
            qB(3,iPOE,ids:ide,kds) = -sqrtGB(iPOE,ids:ide,kds) * G13B(iPOE,ids:ide,kds) * qB(2,iPOE,ids:ide,kds)
            qT(3,iPOE,ids:ide,kde) = -sqrtGT(iPOE,ids:ide,kde) * G13T(iPOE,ids:ide,kde) * qT(2,iPOE,ids:ide,kde)
          enddo
        endif
        
        ! Fill outside boundary
        ! left boundary
        qR(:,:,ids-1,kds:kde) = qL(:,:,ids,kds:kde)
        FR(:,:,ids-1,kds:kde) = FL(:,:,ids,kds:kde)
        
        ! right boundary
        qL(:,:,ide+1,kds:kde) = qR(:,:,ide,kds:kde)
        FL(:,:,ide+1,kds:kde) = FR(:,:,ide,kds:kde)
        
        ! bottom boundary
        qT(:,:,ids:ide,kds-1) = qB(:,:,ids:ide,kds)
        HT(:,:,ids:ide,kds-1) = HB(:,:,ids:ide,kds)
        
        ! top boundary
        qB(:,:,ids:ide,kde+1) = qT(:,:,ids:ide,kde)
        HB(:,:,ids:ide,kde+1) = HT(:,:,ids:ide,kde)
      
        ! Calculate functions X
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
        
        ! Calculate functions Z
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
        
        !$OMP PARALLEL DO PRIVATE(i,im1,iPOE,iVar)
        do k = kds,kde
          do i = ids,ide+1
            do iPOE = 1,nPointsOnEdge
              im1 = i - 1
              FeP(:,iPOE,i,k) = calc_F(sqrtGR(iPOE,im1,k),sqrtGL(iPOE,i,k),qR(:,iPOE,im1,k),qL(:,iPOE,i,k),pR(iPOE,im1,k),pL(iPOE,i,k))
            enddo
            do iVar = 1,nVar
              Fe(iVar,i,k) = Gaussian_quadrature_1d(FeP(iVar,:,i,k))
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
        
        !$OMP PARALLEL DO PRIVATE(k,km1,iPOE,iVar)
        do i = ids,ide
          do k = kds,kde+1
            do iPOE = 1,nPointsOnEdge
              km1 = k - 1
              HeP(:,iPOE,i,k) = calc_H(sqrtGT  (  iPOE,i,km1),sqrtGB  (  iPOE,i,k),G13T  (iPOE,i,km1),G13B  (iPOE,i,k),&
                                       qT      (:,iPOE,i,km1),qB      (:,iPOE,i,k),pT    (iPOE,i,km1),pB    (iPOE,i,k),&
                                       rhoT_ref(  iPOE,i,km1),rhoB_ref(  iPOE,i,k),cT_ref(iPOE,i,km1),cB_ref(iPOE,i,k),&
                                       pT_ref  (  iPOE,i,km1),pB_ref  (  iPOE,i,k))
            enddo
            do iVar = 1,nVar
              He(iVar,i,k) = Gaussian_quadrature_1d(HeP(iVar,:,i,k))
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
        
        !rho_p = ( qC   (1,ids:ide,kds:kde) + qC   (5,ids:ide,kds:kde) &
        !        - ref%q(1,ids:ide,kds:kde) - ref%q(5,ids:ide,kds:kde) ) / sqrtGC(ids:ide,kds:kde)
        rho_p = ( qC   (1,ids:ide,kds:kde) + qC   (5,ids:ide,kds:kde) ) / sqrtGC(ids:ide,kds:kde)
        
        where(abs(rho_p)<=1.E-13)rho_p=0.
        
        src = 0
        src(3,ids:ide,kds:kde) = src(3,ids:ide,kds:kde) - sqrtGC(ids:ide,kds:kde) * rho_p(ids:ide,kds:kde) * gravity
      
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
            tend%q(:,i,k) = - ( ( Fe(:,ip1,k) - Fe(:,i,k) ) / dx + ( He(:,i,kp1) - He(:,i,k) ) / deta ) + src(:,i,k)
          enddo
        enddo
        !$OMP END PARALLEL DO
      end subroutine spatial_operator
    
      subroutine reconstruct_field(qL,qR,qB,qT,qC)
        real(r_kind), dimension(nPointsOnEdge,ics:ice,kcs:kce),intent(out) :: qL
        real(r_kind), dimension(nPointsOnEdge,ics:ice,kcs:kce),intent(out) :: qR
        real(r_kind), dimension(nPointsOnEdge,ics:ice,kcs:kce),intent(out) :: qB
        real(r_kind), dimension(nPointsOnEdge,ics:ice,kcs:kce),intent(out) :: qT
        real(r_kind), dimension(              ics:ice,kcs:kce),intent(in ) :: qC
      
        integer(i_kind) :: i,j,k,iR,kR
        integer(i_kind) :: iRec,kRec
        integer(i_kind) :: ic
        integer(i_kind) :: m,n
        
        real(r_kind), dimension(nRecCells          ) :: u
        real(r_kind)                                 :: h
        real(r_kind), dimension(nRecCells,nRecTerms) :: A
        real(r_kind), dimension(          nRecTerms) :: polyCoef
        
        m = nRecCells
        n = nRecTerms
        
        h = dx
        
        !$OMP PARALLEL DO PRIVATE(i,j,kR,iR,iRec,kRec,iC,u,A,polyCoef)
        !do k = kds,kde
        !  do i = ids,ide
        do k = kbs,kbe
          do i = ibs,ibe
            ! Set variable for reconstruction
            do kR = 1,stencil_width
              do iR = 1,stencil_width
                iRec = iRecCell(iR,kR,i,k)
                kRec = kRecCell(iR,kR,i,k)
                
                j = stencil_width * ( kR - 1 ) + iR
                u(j  ) = qC(iRec,kRec)
                A(j,:) = polyCoordCoef(j,:,i,k)
              enddo
            enddo
            ic = iCenCell(i,k)
            
            polyCoef = WLS_ENO(A,u,h,m,n,ic)
            
            qL(:,i,k) = matmul(recMatrixL,polyCoef)
            qR(:,i,k) = matmul(recMatrixR,polyCoef)
            qB(:,i,k) = matmul(recMatrixB,polyCoef)
            qT(:,i,k) = matmul(recMatrixT,polyCoef)
          enddo
        enddo
        !$OMP END PARALLEL DO
        
        ! bottom
        call reconstruct_bdy(ids,ide,kds,kbs-1,qL,qR,qB,qT,qC)
        ! top
        call reconstruct_bdy(ids,ide,kbe+1,kde,qL,qR,qB,qT,qC)
        ! left
        call reconstruct_bdy(ids,ibs-1,kbs,kbe,qL,qR,qB,qT,qC)
        ! right
        call reconstruct_bdy(ibe+1,ide,kbs,kbe,qL,qR,qB,qT,qC)
        
      end subroutine reconstruct_field
      
      ! reconstruct boundary by weno
      subroutine reconstruct_bdy(is,ie,ks,ke,qL,qR,qB,qT,qC)
        integer(i_kind)                                          , intent(in ) :: is
        integer(i_kind)                                          , intent(in ) :: ie
        integer(i_kind)                                          , intent(in ) :: ks
        integer(i_kind)                                          , intent(in ) :: ke
        real   (r_kind), dimension(nPointsOnEdge,ics:ice,kcs:kce), intent(out) :: qL
        real   (r_kind), dimension(nPointsOnEdge,ics:ice,kcs:kce), intent(out) :: qR
        real   (r_kind), dimension(nPointsOnEdge,ics:ice,kcs:kce), intent(out) :: qB
        real   (r_kind), dimension(nPointsOnEdge,ics:ice,kcs:kce), intent(out) :: qT
        real   (r_kind), dimension(              ics:ice,kcs:kce), intent(in ) :: qC
        
        real   (r_kind) :: q_weno(5)
        integer(i_kind) :: dir
        
        integer(i_kind) :: i,k
      
        integer(i_kind) :: ip1,im1
        integer(i_kind) :: kp1,km1
        integer(i_kind) :: ip2,im2
        integer(i_kind) :: kp2,km2
        
        !$OMP PARALLEL DO PRIVATE(kp1,km1,kp2,km2,i,ip1,im1,ip2,im2,q_weno,dir)
        do k = ks,ke
          kp1 = k + 1
          km1 = k - 1
          kp2 = k + 2
          km2 = k - 2
          do i = is,ie
            ip1 = i + 1
            im1 = i - 1
            ip2 = i + 2
            im2 = i - 2
            
            ! x-dir
            q_weno = qC(im2:ip2,k)
            
            dir = -1
            call WENO_limiter(qL(1,i,k),q_weno,dir)
            dir = 1
            call WENO_limiter(qR(1,i,k),q_weno,dir)
            
            ! z-dir
            q_weno = qC(i,km2:kp2)
            
            dir = -1
            call WENO_limiter(qB(1,i,k),q_weno,dir)
            dir = 1
            call WENO_limiter(qT(1,i,k),q_weno,dir)
          
            qL(2:nPointsOnEdge,i,k) = qL(1,i,k)
            qR(2:nPointsOnEdge,i,k) = qR(1,i,k)
            qB(2:nPointsOnEdge,i,k) = qB(1,i,k)
            qT(2:nPointsOnEdge,i,k) = qT(1,i,k)
          enddo
        enddo
        !$OMP END PARALLEL DO
      
      end subroutine reconstruct_bdy
      
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
        
        real(r_kind) rhoL
        real(r_kind) rhoR
        
        real(r_kind) uL
        real(r_kind) uR
        
        real(r_kind) cL
        real(r_kind) cR
        
        real(r_kind) m
        real(r_kind) p ! sqrtG * p
        
        rhoL = ( qL(1) + qL(5) ) / sqrtGL
        rhoR = ( qR(1) + qR(5) ) / sqrtGR
        uL   = qL(2) / ( qL(1) + qL(5) )
        uR   = qR(2) / ( qR(1) + qR(5) )
        cL   = calc_sound_speed_x(sqrtGL,qL)
        cR   = calc_sound_speed_x(sqrtGR,qR)
        
        call AUSM_up_x(m,p,sqrtGL,sqrtGR,rhoL,rhoR,uL,uR,pL,pR,cL,cR)
        
        calc_F = 0.5 * m * ( qL + qR - sign(1.,m) * ( qR - qL ) )
        calc_F(2) = calc_F(2) + p
        
      end function calc_F
      
      function calc_H(sqrtGL,sqrtGR,G13L,G13R,qL,qR,pL,pR,rhoL_ref,rhoR_ref,cL_ref,cR_ref,pL_ref,pR_ref)
        real(r_kind) :: calc_H(nVar)
        real(r_kind) :: sqrtGL
        real(r_kind) :: sqrtGR
        real(r_kind) :: G13L
        real(r_kind) :: G13R
        real(r_kind) :: qL(nVar)
        real(r_kind) :: qR(nVar)
        real(r_kind) :: pL        ! pressure
        real(r_kind) :: pR        ! pressure
        real(r_kind) :: rhoL_ref  ! reference density
        real(r_kind) :: rhoR_ref  ! reference density
        real(r_kind) :: cL_ref    ! reference sound speed
        real(r_kind) :: cR_ref    ! reference sound speed
        real(r_kind) :: pL_ref    ! reference pressure
        real(r_kind) :: pR_ref    ! reference pressure
      
        real(r_kind) :: pL_pert ! pressure  perturbation
        real(r_kind) :: pR_pert ! pressure  perturbation
        
        
        real(r_kind) rhoL
        real(r_kind) rhoR
        
        real(r_kind) uL
        real(r_kind) uR
        
        real(r_kind) cL
        real(r_kind) cR
        
        real(r_kind) m
        real(r_kind) p      ! sqrtG * G13 * p
        real(r_kind) p_pert ! p - p_ref
        
        rhoL = ( qL(1) + qL(5) ) / sqrtGL
        rhoR = ( qR(1) + qR(5) ) / sqrtGR
        uL   = ( qL(3) + sqrtGL * G13L * qL(2) ) / ( qL(1) + qL(5) ) ! a sqrtG has been multipled here
        uR   = ( qR(3) + sqrtGR * G13R * qR(2) ) / ( qR(1) + qR(5) ) ! a sqrtG has been multipled here
        cL   = calc_sound_speed_z(sqrtGL,G13L,qL) * sqrtGL
        cR   = calc_sound_speed_z(sqrtGR,G13R,qR) * sqrtGR
        
        call AUSM_up_z(m, p, p_pert,                                                  &
                       sqrtGL, sqrtGR, G13L, G13R, rhoL, rhoR, uL, uR, cL, cR, pL, pR,&
                       rhoL_ref, rhoR_ref, cL_ref, cR_ref, pL_ref, pR_ref)
        
        calc_H = 0.5 * m * ( qL + qR - sign(1.,m) * ( qR - qL ) )
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
        
        if(any(q==FillValue).or.sqrtG==FillValue)then
          calc_sound_speed_x = 0
        else
          w1 = q(1)
          w2 = q(2)
          w3 = q(3)
          w4 = q(4)
          w5 = q(5)
          
          coef1 = sqrt( p0*sqrtG*w1**2*w4**2*(w1 + w5)**3*(cpd*w1 + cpv*w5)*&
                (cvd*w1 + cvv*w5)**3*(w1 + w5 + &
                eq*w5)**2*((Rd*w4*(w1 + w5 + eq*w5))/(p0*sqrtG*w1))**&
                ((cpd*w1 + cpv*w5)/(cvd*w1 + cvv*w5)) )
          
          coef2 = w1*w4*(w1 + w5)**2*(cvd*w1 + cvv*w5)**2*(w1 + w5 + eq*w5)
          
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
        
        if(any(q==FillValue))then
          calc_sound_speed_z = 0
        else
          w1 = q(1)
          w2 = q(2)
          w3 = q(3)
          w4 = q(4)
          w5 = q(5)
          
          coef1 = sqrt( p0*sqrtG*(1 + G13**2*sqrtG**2)*w1**2*                    &
                w4**2*(w1 + w5)**3*(cpd*w1 + cpv*w5)*(cvd*w1 + cvv*w5)**3*       &
                (w1 + w5 + eq*w5)**2*((Rd*w4*(w1 + w5 + eq*w5))/(p0*sqrtG*w1))** &
                ((cpd*w1 + cpv*w5)/(cvd*w1 + cvv*w5)) )
          
          coef2 = sqrtG*w1*w4*(w1 + w5)**2*(cvd*w1 + cvv*w5)**2*(w1 + w5 + eq*w5)
          
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
        
        rho = 0.5 * ( rhoL + rhoR )
        a   = 0.5 * ( cL + cR )
        
        ML = uL / a
        MR = uR / a
        
        Mbar2 = ( uL**2 + uR**2 ) / ( 2. * a**2 )
        
        Mh = M4( ML, sp ) + M4( MR, sn ) - Kp * max( 1. - sigma * Mbar2, 0. ) * ( PR - PL ) / ( rho * a**2 )
        m  = a * Mh
        
        p = P5(ML,sp) * sqrtGL * PL + P5(MR,sn) * sqrtGR * PR - Ku * P5(ML,sp) * P5(MR,sn) * ( sqrtGL * rhoL + sqrtGR * rhoR ) * a * ( uR - uL )
        
      end subroutine AUSM_up_x
      
      subroutine AUSM_up_z(m, p, p_pert,                                                  & ! output
                           sqrtGL, sqrtGR, G13L, G13R, rhoL, rhoR, uL, uR, cL, cR, pL, pR,& ! input state
                           rhoL_ref, rhoR_ref, cL_ref, cR_ref, pL_ref, pR_ref)              ! input reference state
        real(r_kind),intent(out) :: m
        real(r_kind),intent(out) :: p
        real(r_kind),intent(out) :: p_pert
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
        real(r_kind),intent(in ) :: rhoL_ref ! reference state
        real(r_kind),intent(in ) :: rhoR_ref ! reference state
        real(r_kind),intent(in ) :: cL_ref   ! reference state
        real(r_kind),intent(in ) :: cR_ref   ! reference state
        real(r_kind),intent(in ) :: pL_ref   ! reference state
        real(r_kind),intent(in ) :: pR_ref   ! reference state
        
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
        
        real(r_kind) :: rho_ref
        real(r_kind) :: a_ref
        
        real(r_kind) :: pL_pert
        real(r_kind) :: pR_pert
        
        real(r_kind) :: p_diff
        
        real(r_kind) :: coefL
        real(r_kind) :: coefR
        
        rho = 0.5 * ( rhoL + rhoR )
        a   = 0.5 * ( cL + cR )
        
        rho_ref = 0.5 * ( rhoL_ref + rhoR_ref )
        a_ref   = 0.5 * ( cL_ref + cR_ref )
        
        pL_pert = pL - pL_ref
        pR_pert = pR - pR_ref
        
        ML = uL / a
        MR = uR / a
        
        Mbar2 = ( uL**2 + uR**2 ) / ( 2. * a**2 )
        
        !if( pL_pert/pL_ref<1.e-13 .and. pR_pert/pR_ref<1.e-13 )then
        !  p_diff = 0
        !else
        !  p_diff = PR - PL
        !endif
        
        p_diff = PR - PL
        
        Mh     = M4( ML, sp ) + M4( MR, sn ) - Kp * max( 1. - sigma * Mbar2, 0. ) * p_diff / ( rho * a**2 )
        m      = 0.5 * ( cL / sqrtGL + cR / sqrtGR ) * Mh
        
        coefL = sqrtGL * G13L
        coefR = sqrtGR * G13R
        
        p = P5(ML,sp) * coefL * PL + P5(MR,sn) * coefR * PR &
          - Ku * P5(ML,sp) * P5(MR,sn) * ( coefL * rhoL + coefR * rhoR ) * a * ( uR - uL )
        
        !p_pert = P5(ML,sp) * pL_pert + P5(MR,sn) * pR_pert &
        !       - 2. * Ku * P5(ML,sp) * P5(MR,sn) * ( rho * a - rho_ref * a_ref  ) * ( uR - uL )
        p_pert = P5(ML,sp) * pL + P5(MR,sn) * pR &
               - 2. * Ku * P5(ML,sp) * P5(MR,sn) * rho * a * ( uR - uL )
        
        !if( abs(2.*p_pert/(pL+pR)) < 1.e-13 )p_pert = 0
        
      end subroutine AUSM_up_z
      
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
          P5 = 0.5 * ( 1. + signal * sign(1.,M) )
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
        
        do iVar = vs,ve
          q(iVar,ids:ide,kds:kde) = q(iVar,ids:ide,kds:kde) - relax_coef(ids:ide,kds:kde) * ( q(iVar,ids:ide,kds:kde) - q_ref(iVar,ids:ide,kds:kde) )
        enddo
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

