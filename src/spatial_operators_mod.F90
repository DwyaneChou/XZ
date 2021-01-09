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
      
      integer(i_kind), dimension(:,:,:,:), allocatable :: iRecCell ! x index of reconstruction cells
      integer(i_kind), dimension(:,:,:,:), allocatable :: kRecCell ! k index of reconstruction cells
      
      real   (r_kind), dimension(:,:,:,:,:), allocatable :: xRel   ! relative x coordinate of reconstruction cells
      real   (r_kind), dimension(:,:,:,:,:), allocatable :: etaRel ! relative eta coordinate of reconstruction cells
      
      real   (r_kind), dimension(:,:,:), allocatable :: disCenter ! distance between adjacent cells and center cell on stencil
      
    contains
      subroutine init_spatial_operator
        integer(i_kind) :: i,j,k,iR,kR
        integer(i_kind) :: ibs,ibe,kbs,kbe
        integer(i_kind) :: recBdy
        integer(i_kind) :: iRec,kRec
        
        allocate(iRecCell  (stencil_width,stencil_width,ids:ide,kds:kde))
        allocate(kRecCell  (stencil_width,stencil_width,ids:ide,kds:kde))
        
        allocate(xRel  (4,stencil_width,stencil_width,ids:ide,kds:kde))
        allocate(etaRel(4,stencil_width,stencil_width,ids:ide,kds:kde))
        
        allocate(disCenter(nRecCells,ids:ide,kds:kde))
        
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
                
                xRel  (:,iR,kR,i,k) = xCorner  (:,iRec,kRec) - xCenter  (i,k)
                etaRel(:,iR,kR,i,k) = etaCorner(:,iRec,kRec) - etaCenter(i,k)
                
                disCenter(j,i,k) = sqrt( ( xCenter(iRec,kRec) - xCenter(i,k) )**2 + ( etaCenter(iRec,kRec) - etaCenter(i,k) )**2 )
                
                call calc_polynomial_square_integration(recPolyDegree,xRel(1,iR,kR,i,k),xRel(2,iR,kR,i,k),etaRel(1,iR,kR,i,k),etaRel(4,iR,kR,i,k),polyCoordCoef(j,:,i,k))
              enddo
            enddo
          enddo
        enddo
        
      end subroutine init_spatial_operator
        
    end module spatial_operators_mod

