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
      
      integer(i_kind), dimension(:,:,:,:), allocatable :: iRec ! x index of reconstruction cells
      integer(i_kind), dimension(:,:,:,:), allocatable :: kRec ! k index of reconstruction cells
      
    contains
      subroutine init_spatial_operator
        integer(i_kind) :: i,j,k,iR,kR
        integer(i_kind) :: ibs,ibe,kbs,kbe
        integer(i_kind) :: recBdy
        
        allocate(iRec(stencil_width,stencil_width,ids:ide,kds:kde))
        allocate(kRec(stencil_width,stencil_width,ids:ide,kds:kde))
        
        ! set reconstruction cells on each cell
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
                iRec(iR,kR,i,k) = iR + i - 1 - recBdy
              enddo
            enddo            
          enddo
          ! left
          do i = ids,ibs-1
            do kR = 1,stencil_width
              do iR = 1,stencil_width
                iRec(iR,kR,i,k) = ids - 1 + iR
              enddo
            enddo
          enddo
          ! right
          do i = ibe+1,ide
            do kR = 1,stencil_width
              do iR = 1,stencil_width
                iRec(iR,kR,i,k) = ide - stencil_width + iR
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
                kRec(iR,kR,i,k) = kR + k - 1 - recBdy
              enddo
            enddo            
          enddo
          ! bottom
          do k = kds,kbs-1
            do kR = 1,stencil_width
              do iR = 1,stencil_width
                kRec(iR,kR,i,k) = kds - 1 + kR
              enddo
            enddo
          enddo
          ! top
          do k = kbe+1,kde
            do kR = 1,stencil_width
              do iR = 1,stencil_width
                kRec(iR,kR,i,k) = kde - stencil_width  + kR
              enddo
            enddo
          enddo
        enddo
        
        !i = ide
        !k = kde-1
        !do kR = 1,stencil_width
        !  print*,kR,(kRec(iR,kR,i,k),iR=1,stencil_width)
        !enddo
        !do iR = 1,stencil_width
        !  print*,iR,(iRec(iR,kR,i,k),kR=1,stencil_width)
        !enddo
        
        ! initilize polyCoordCoef
        do k = kcs,kce
          do i = ics,ice
            call calc_polynomial_square_integration(recPolyDegree,xCorner(1,i,k),xCorner(2,i,k),etaCorner(1,i,k),etaCorner(4,i,k),polyCoordCoef(:,i,k))
          enddo
        enddo
        
      end subroutine init_spatial_operator
        
    end module spatial_operators_mod

