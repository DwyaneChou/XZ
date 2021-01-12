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
      
      public init_spatial_operator!,spatial_operator
      
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
      
    contains
      subroutine init_spatial_operator
        integer(i_kind) :: i,j,k,iR,kR
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
        do k = kds,kde
          do i = ids,ide
            j = 0 ! Reset cell number on stencil
            do kR = 1,stencil_width
              do iR = 1,stencil_width
                iRec = iRecCell(iR,kR,i,k)
                kRec = kRecCell(iR,kR,i,k)
                
                j    = j + 1
                
                if(iRec==i.and.kRec==k)iCenCell(i,k) = j
                
                xRel  (:,iR,kR,i,k) = xCorner  (:,iRec,kRec) - xCenter  (i,k)
                etaRel(:,iR,kR,i,k) = etaCorner(:,iRec,kRec) - etaCenter(i,k)
                
                disCenter(j,i,k) = sqrt( ( xCenter(iRec,kRec) - xCenter(i,k) )**2 + ( etaCenter(iRec,kRec) - etaCenter(i,k) )**2 )
                
                call calc_polynomial_square_integration(recPolyDegree,xRel(1,iR,kR,i,k),xRel(2,iR,kR,i,k),etaRel(1,iR,kR,i,k),etaRel(4,iR,kR,i,k),polyCoordCoef(j,:,i,k))
              enddo
            enddo
          enddo
        enddo
        
        ! Calculate reconstruction matrix on edge
        call calc_polynomial_matrix(recPolyDegree,nPointsOnEdge,nRecTerms,xL,etaL,recMatrixL)
        call calc_polynomial_matrix(recPolyDegree,nPointsOnEdge,nRecTerms,xR,etaR,recMatrixR)
        call calc_polynomial_matrix(recPolyDegree,nPointsOnEdge,nRecTerms,xB,etaB,recMatrixB)
        call calc_polynomial_matrix(recPolyDegree,nPointsOnEdge,nRecTerms,xT,etaT,recMatrixT)
        
        ! Reconstruct metric function
        call reconstruct_field(sqrtGL,sqrtGR,sqrtGB,sqrtGT,sqrtGC)
        
      end subroutine init_spatial_operator
    
      subroutine reconstruct_field(qL,qR,qB,qT,qC)
        real(r_kind), dimension(:,:,:),intent(out) :: qL
        real(r_kind), dimension(:,:,:),intent(out) :: qR
        real(r_kind), dimension(:,:,:),intent(out) :: qB
        real(r_kind), dimension(:,:,:),intent(out) :: qT
        real(r_kind), dimension(  :,:),intent(in ) :: qC
      
        integer(i_kind) :: i,j,k,iR,kR
        integer(i_kind) :: iRec,kRec
        integer(i_kind) :: ic
        integer(i_kind) :: m,n
        
        real(r_kind), dimension(nRecCells          ) :: u
        real(r_kind), dimension(nRecCells          ) :: h
        real(r_kind), dimension(nRecCells,nRecTerms) :: A
        real(r_kind), dimension(          nRecTerms) :: polyCoef
        
        real(r_kind) dV
        
        m = nRecCells
        n = nRecTerms
        
        dV = dx * deta
        
        do k = kds,kde
          do i = ids,ide
            ! Set variable for reconstruction
            j = 0
            do kR = 1,stencil_width
              do iR = 1,stencil_width
                iRec = iRecCell(iR,kR,i,k)
                kRec = kRecCell(iR,kR,i,k)
                j = j + 1
                u(j  ) = qC(iRec,kRec)
                h(j  ) = disCenter(j,i,k)
                A(j,:) = polyCoordCoef(j,:,i,k)
              enddo
            enddo
            ic = iCenCell(i,k)
            
            polyCoef = WLS_ENO(A,u,h,m,n,ic)
            
            qL(:,i,k) = matmul(recMatrixL,polyCoef) * dV
            qR(:,i,k) = matmul(recMatrixR,polyCoef) * dV
            qB(:,i,k) = matmul(recMatrixB,polyCoef) * dV
            qT(:,i,k) = matmul(recMatrixT,polyCoef) * dV
            
            !print*,polyCoef
            !print*,''
            
            print*,qL(:,i,k)
            print*,''
            print*,qR(:,i,k)
            print*,''
            print*,qB(:,i,k)
            print*,''
            print*,qT(:,i,k)
            print*,''
          enddo
        enddo
        
      end subroutine reconstruct_field
      
    end module spatial_operators_mod

