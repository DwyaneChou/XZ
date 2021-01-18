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
      
      real(r_kind), dimension(  :,:), allocatable :: PC_ref
      real(r_kind), dimension(:,:,:), allocatable :: PL_ref! Reconstructed P_ref_(i-1/2,k)
      real(r_kind), dimension(:,:,:), allocatable :: PR_ref! Reconstructed P_ref_(i+1/2,k)
      real(r_kind), dimension(:,:,:), allocatable :: PB_ref! Reconstructed P_ref_(i,k-1/2)
      real(r_kind), dimension(:,:,:), allocatable :: PT_ref! Reconstructed P_ref_(i,k+1/2)
  
      real(r_kind), dimension(:,:), allocatable :: eig_x
      real(r_kind), dimension(:,:), allocatable :: eig_z
      
    contains
      subroutine init_spatial_operator
        integer(i_kind) :: i,j,k,iR,kR,iVar,iPOE
        integer(i_kind) :: ibs,ibe,kbs,kbe
        integer(i_kind) :: recBdy
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
        
        allocate(PC_ref(              ics:ice,kcs:kce))
        allocate(PL_ref(nPointsOnEdge,ics:ice,kcs:kce))
        allocate(PR_ref(nPointsOnEdge,ics:ice,kcs:kce))
        allocate(PB_ref(nPointsOnEdge,ics:ice,kcs:kce))
        allocate(PT_ref(nPointsOnEdge,ics:ice,kcs:kce))
      
        allocate(eig_x(ics:ice,kcs:kce))
        allocate(eig_z(ics:ice,kcs:kce))
      
        src   = 0
        eig_x = 0
        eig_z = 0
        
        dV = dx * deta
        
        ! set reconstruction cells on each stencil
        print*,'Set reconstruction cells on each stencil'
        print*,''
        recBdy = ( stencil_width - 1 ) / 2
        ibs    = ids + recBdy
        ibe    = ide - recBdy
        kbs    = kds + recBdy
        kbe    = kde - recBdy
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
        
        !$OMP PARALLEL DO PRIVATE(i,iPOE)
        do k = kds,kde
          do i = ids,ide
            do iPOE = 1,nPointsOnEdge
              PL_ref(iPOE,i,k) = calc_pressure(sqrtGL(iPOE,i,k),qL(:,iPOE,i,k))
              PR_ref(iPOE,i,k) = calc_pressure(sqrtGR(iPOE,i,k),qR(:,iPOE,i,k))
              PB_ref(iPOE,i,k) = calc_pressure(sqrtGB(iPOE,i,k),qB(:,iPOE,i,k))
              PT_ref(iPOE,i,k) = calc_pressure(sqrtGT(iPOE,i,k),qT(:,iPOE,i,k))
              
              ! Reset pressure at domain boundary for limiting oscillation
              if(k==kds.or.k==kde)then
                PB_ref(iPOE,i,k) = calc_pressure(sqrtGC(i,k),ref%q(:,i,k))
                PT_ref(iPOE,i,k) = PB_ref(iPOE,i,k)!calc_pressure(sqrtGC(i,k),ref%q(:,i,k))
              endif
              if(i==ids.or.i==ide)then
                PL_ref(iPOE,i,k) = calc_pressure(sqrtGC(i,k),ref%q(:,i,k))
                PR_ref(iPOE,i,k) = PL_ref(iPOE,i,k)!calc_pressure(sqrtGC(i,k),ref%q(:,i,k))
              endif
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
        
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
        
        ! Calculate functions
        !$OMP PARALLEL DO PRIVATE(i,iPOE)
        do k = kds,kde
          do i = ids,ide
            do iPOE = 1,nPointsOnEdge
              PL(iPOE,i,k) = calc_pressure(sqrtGL(iPOE,i,k),qL(:,iPOE,i,k))
              PR(iPOE,i,k) = calc_pressure(sqrtGR(iPOE,i,k),qR(:,iPOE,i,k))
              PB(iPOE,i,k) = calc_pressure(sqrtGB(iPOE,i,k),qB(:,iPOE,i,k))
              PT(iPOE,i,k) = calc_pressure(sqrtGT(iPOE,i,k),qT(:,iPOE,i,k))
              
              ! Reset pressure at domain boundary for limiting oscillation
              if(k==kds.or.k==kde)then
                PB(iPOE,i,k) = calc_pressure(sqrtGC(i,k),qC(:,i,k))
                PT(iPOE,i,k) = PB(iPOE,i,k)!calc_pressure(sqrtGC(i,k),qC(:,i,k))
              endif
              if(i==ids.or.i==ide)then
                PL(iPOE,i,k) = calc_pressure(sqrtGC(i,k),qC(:,i,k))
                PR(iPOE,i,k) = PL(iPOE,i,k)!calc_pressure(sqrtGC(i,k),qC(:,i,k))
              endif
              
              FL(:,iPOE,i,k) = calc_F(sqrtGL(iPOE,i,k),qL(:,iPOE,i,k),PL(iPOE,i,k))
              FR(:,iPOE,i,k) = calc_F(sqrtGR(iPOE,i,k),qR(:,iPOE,i,k),PR(iPOE,i,k))
              HB(:,iPOE,i,k) = calc_H(sqrtGB(iPOE,i,k),G13B(iPOE,i,k),qB(:,iPOE,i,k),PB(iPOE,i,k),PB_ref(iPOE,i,k))
              HT(:,iPOE,i,k) = calc_H(sqrtGT(iPOE,i,k),G13T(iPOE,i,k),qT(:,iPOE,i,k),PT(iPOE,i,k),PT_ref(iPOE,i,k))
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
        
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
        
        rho_p = ( qC   (1,ids:ide,kds:kde) + qC   (5,ids:ide,kds:kde) &
                - ref%q(1,ids:ide,kds:kde) - ref%q(5,ids:ide,kds:kde) ) / sqrtGC(ids:ide,kds:kde)
        
        !where(abs(rho_p)<=1.E-13)rho_p=0.
        
        src = 0
        src(3,ids:ide,kds:kde) = src(3,ids:ide,kds:kde) - sqrtGC(ids:ide,kds:kde) * rho_p(ids:ide,kds:kde) * gravity
        
        ! Calculate eigenvalues x-dir
        !$OMP PARALLEL DO PRIVATE(i)
        do k = kds,kde
          do i = ids,ide
            eig_x(i,k) = calc_eigenvalue_x(sqrtGC(i,k),qC(:,i,k))
          enddo
        enddo
        !$OMP END PARALLEL DO
        
        ! Calculate eigenvalues z-dir
        !$OMP PARALLEL DO PRIVATE(k)
        do i = ids,ide
          do k = kds,kde
            eig_z(i,k) = calc_eigenvalue_z(sqrtGC(i,k),G13C(i,k),qC(:,i,k))
          enddo
        enddo
        !$OMP END PARALLEL DO

        ! calc x flux
        !$OMP PARALLEL DO PRIVATE(i,im1,maxeigen_x,iVar,iPOE)
        do k = kds,kde
          do i = ids,ide+1
            im1 = i - 1
            
            maxeigen_x = max(eig_x(i,k),eig_x(im1,k))
            !maxeigen_x = max(eig_x(i,k),eig_x(im1,k),eig_x(i,k+1),eig_x(i,k-1),eig_x(im1,k+1),eig_x(im1,k-1))
            
            do iVar = 1,nVar
              do iPOE = 1,nPointsOnEdge
                if(abs(FL(iVar,iPOE,i,k) + FR(iVar,iPOE,im1,k))<=1.E-15)then
                  FeP(iVar,iPOE,i,k) = 0
                else
                  FeP(iVar,iPOE,i,k) = 0.5 * ( FL(iVar,iPOE,i,k) + FR(iVar,iPOE,im1,k) - maxeigen_x * ( qL(iVar,iPOE,i,k) - qR(iVar,iPOE,im1,k) ) )
                endif
              enddo
              Fe(iVar,i,k) = Gaussian_quadrature_1d(FeP(iVar,:,i,k))
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
        
        ! calc z flux
        !$OMP PARALLEL DO PRIVATE(k,km1,maxeigen_z,iVar,iPOE)
        do i = ids,ide
          do k = kds,kde+1
            km1 = k - 1
            
            maxeigen_z = max(eig_z(i,k),eig_z(i,km1))
            !maxeigen_z = max(eig_z(i,k),eig_z(i,km1),eig_z(i+1,k),eig_z(i-1,k),eig_z(i+1,km1),eig_z(i-1,km1))
            
            do iVar = 1,nVar
              do iPOE = 1,nPointsOnEdge
                if(abs(HB(iVar,iPOE,i,k) + HT(iVar,iPOE,i,km1))<=1.E-15)then
                  HeP(iVar,iPOE,i,k) = 0
                else
                  HeP(iVar,iPOE,i,k) = 0.5 * ( HB(iVar,iPOE,i,k) + HT(iVar,iPOE,i,km1) - maxeigen_z * ( qB(iVar,iPOE,i,k) - qT(iVar,iPOE,i,km1) ) )
                endif
              enddo
              He(iVar,i,k) = Gaussian_quadrature_1d(HeP(iVar,:,i,k))
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
        
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
        do k = kds,kde
          do i = ids,ide
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
      
      function calc_F(sqrtG,q,p)
        real(r_kind),dimension(5) :: calc_F
        real(r_kind)              :: sqrtG
        real(r_kind),dimension(5) :: q(5)
        real(r_kind)              :: p      ! pressure
        
        real(r_kind) w1
        real(r_kind) w2
        real(r_kind) w3
        real(r_kind) w4
        real(r_kind) w5
        
        real(r_kind) sqrtGrho
        real(r_kind) u
        
        w1 = q(1)
        w2 = q(2)
        w3 = q(3)
        w4 = q(4)
        w5 = q(5)
        
        sqrtGrho = w1 + w5
        u        = w2 / sqrtGrho
        
        calc_F(1) = w1 * u
        calc_F(2) = w2 * u + sqrtG * p
        calc_F(3) = w3 * u
        calc_F(4) = w4 * u
        calc_F(5) = w5 * u
        
      end function calc_F
      
      function calc_H(sqrtG,G13,q,p,p_ref)
        real(r_kind),dimension(5) :: calc_H
        real(r_kind)              :: sqrtG
        real(r_kind)              :: G13
        real(r_kind),dimension(5) :: q(5)
        real(r_kind)              :: p      ! pressure
        real(r_kind)              :: p_ref  ! reference pressure
        
        real(r_kind)              :: p_pert ! pressure  perturbation
        
        real(r_kind) w1
        real(r_kind) w2
        real(r_kind) w3
        real(r_kind) w4
        real(r_kind) w5
        real(r_kind) ww
        
        real(r_kind) sqrtGrho
        real(r_kind) u
        real(r_kind) w
        
        w1 = q(1)
        w2 = q(2)
        w3 = q(3)
        w4 = q(4)
        w5 = q(5)
        
        sqrtGrho = w1 + w5
        u        = w2 / sqrtGrho
        w        = w3 / sqrtGrho
        p_pert   = p - p_ref
        
        ww = w / sqrtG + G13 * u
        
        calc_H(1) = w1 * ww
        calc_H(2) = w2 * ww + sqrtG * G13 * p
        calc_H(3) = w3 * ww + p_pert
        calc_H(4) = w4 * ww
        calc_H(5) = w5 * ww
        
      end function calc_H
      
      function calc_eigenvalue_x(sqrtG,q)
        real(r_kind) :: calc_eigenvalue_x
        real(r_kind) :: sqrtG
        real(r_kind) :: q(5)
        real(r_kind) :: eig(5)
        
        real(r_kind) w1
        real(r_kind) w2
        real(r_kind) w3
        real(r_kind) w4
        real(r_kind) w5
        
        real(r_kind) coef1,coef2,coef3
        
        if(any(q==FillValue).or.sqrtG==FillValue)then
          eig = 0
        else
          w1 = q(1)
          w2 = q(2)
          w3 = q(3)
          w4 = q(4)
          w5 = q(5)
          
          coef1 = cvd**2*w1**3*w2*w4*(w1 + w5)*(w1 + w5 + eq*w5) + &
                2.*cvd*cvv*w1**2*w2*w4*w5*(w1 + w5)*(w1 + w5 + eq*w5) + &
                cvv**2*w1*w2*w4*w5**2*(w1 + w5)*(w1 + w5 + eq*w5)
          
          coef2 = sqrt( p0*sqrtG*w1**2*w4**2*(w1 + w5)**3*(cpd*w1 + cpv*w5)*&
                (cvd*w1 + cvv*w5)**3*(w1 + w5 + &
                eq*w5)**2*((Rd*w4*(w1 + w5 + eq*w5))/(p0*sqrtG*w1))**&
                ((cpd*w1 + cpv*w5)/(cvd*w1 + cvv*w5)) )
          
          coef3 = w1*w4*(w1 + w5)**2*(cvd*w1 + cvv*w5)**2*(w1 + w5 + eq*w5)
          
          eig(1) = w2 / ( w1 + w5 )
          eig(2) = eig(1)
          eig(3) = eig(1)
          eig(4) = ( coef1 - coef2 ) / coef3
          eig(5) = ( coef1 + coef2 ) / coef3
        endif
        
        calc_eigenvalue_x = maxval(abs(eig))
        
      end function calc_eigenvalue_x
      
      function calc_eigenvalue_z(sqrtG,G13,q)
        real(r_kind) :: calc_eigenvalue_z
        real(r_kind) :: sqrtG
        real(r_kind) :: G13
        real(r_kind) :: q(5)
        real(r_kind) :: eig(5)
        
        real(r_kind) w1
        real(r_kind) w2
        real(r_kind) w3
        real(r_kind) w4
        real(r_kind) w5
        
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
          w5 = q(5)
          
          !sqrtGrho = w1 + w5
          !u        = w2 / sqrtGrho
          !w        = w3 / sqrtGrho
          
          coef1 = cvd**2*w1**3*(G13*sqrtG*w2 + w3)*w4*(w1 + w5)*(w1 + w5 + eq*w5) + &
                2.*cvd*cvv*w1**2*(G13*sqrtG*w2 + w3)*w4*                            &
                w5*(w1 + w5)*(w1 + w5 + eq*w5) +                                    &
                cvv**2*w1*(G13*sqrtG*w2 + w3)*w4*w5**2*(w1 + w5)*(w1 + w5 + eq*w5)
          
          coef2 = sqrt( p0*sqrtG*(1 + G13**2*sqrtG**2)*w1**2*                    &
                w4**2*(w1 + w5)**3*(cpd*w1 + cpv*w5)*(cvd*w1 + cvv*w5)**3*       &
                (w1 + w5 + eq*w5)**2*((Rd*w4*(w1 + w5 + eq*w5))/(p0*sqrtG*w1))** &
                ((cpd*w1 + cpv*w5)/(cvd*w1 + cvv*w5)) )
          
          coef3 = sqrtG*w1*w4*(w1 + w5)**2*(cvd*w1 + cvv*w5)**2*(w1 + w5 + eq*w5)
          
          drhoetadt = (G13*sqrtG*w2 + w3)/(sqrtG*w1 + sqrtG*w5)
          
          eig(1) = drhoetadt
          eig(2) = drhoetadt
          eig(3) = drhoetadt
          eig(4) = ( coef1 - coef2 ) / coef3
          eig(5) = ( coef1 + coef2 ) / coef3
        endif
        
        calc_eigenvalue_z = maxval(abs(eig))
        
      end function calc_eigenvalue_z
    end module spatial_operators_mod

