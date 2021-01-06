    module mesh_mod
      use constants_mod
      use parameters_mod
      use math_mod
      implicit none
      
      integer(i_kind), parameter :: nEdgesOnCell      = 4
      integer(i_kind), parameter :: nPointsOnEdge     = 5
      integer(i_kind), parameter :: nQuadPointsOnCell = nPointsOnEdge**2 !n Quadrature Points On Cell
      
      real(r_kind) deta
      real(r_kind) dxi
      
      real(r_kind), dimension(:,:,:), allocatable :: xCorner   ! 1: low left, 2: low right, 3: up right, 4:  up left
      real(r_kind), dimension(:,:,:), allocatable :: etaCorner ! 1: low left, 2: low right, 3: up right, 4:  up left
      real(r_kind), dimension(:,:,:), allocatable :: zCorner   ! 1: low left, 2: low right, 3: up right, 4:  up left
      
      real(r_kind), dimension(:,:,:), allocatable :: xQP   ! x coordinate of Quadrature Points
      real(r_kind), dimension(:,:,:), allocatable :: etaQP ! eta coordinate of Quadrature Points
      real(r_kind), dimension(:,:,:), allocatable :: xiQP  ! xi coordinate of Quadrature Points
      real(r_kind), dimension(:,:,:), allocatable :: zQP   ! z coordinate of Quadrature Points
      
      real   (r_kind), dimension(:,:,:), allocatable :: xL
      real   (r_kind), dimension(:,:,:), allocatable :: xR
      real   (r_kind), dimension(:,:,:), allocatable :: xB
      real   (r_kind), dimension(:,:,:), allocatable :: xT
      
      real   (r_kind), dimension(:,:,:), allocatable :: zL
      real   (r_kind), dimension(:,:,:), allocatable :: zR
      real   (r_kind), dimension(:,:,:), allocatable :: zB
      real   (r_kind), dimension(:,:,:), allocatable :: zT
      
      real(r_kind), dimension(:,:,:), allocatable :: zsQP ! topography on Quadrature Points
      
      real(r_kind), dimension(:,:,:), allocatable :: sqrtGQP ! sqrtG on Quadrature Points
      real(r_kind), dimension(:,:,:), allocatable :: G13QP   ! G13 on Quadrature Points
      
      real(r_kind), dimension(:,:), allocatable :: sqrtG ! sqrtG on Cell
      real(r_kind), dimension(:,:), allocatable :: G13   ! G13 on Cell
      
      real   (r_kind), dimension(:,:,:), allocatable :: sqrtGL ! sqrtGL on RiemannPoints
      real   (r_kind), dimension(:,:,:), allocatable :: sqrtGR ! sqrtGR on RiemannPoints
      real   (r_kind), dimension(:,:,:), allocatable :: sqrtGB ! sqrtGB on RiemannPoints
      real   (r_kind), dimension(:,:,:), allocatable :: sqrtGT ! sqrtGT on RiemannPoints
      
      real   (r_kind), dimension(:,:,:), allocatable :: G13L   ! G13L on RiemannPoints
      real   (r_kind), dimension(:,:,:), allocatable :: G13R   ! G13R on RiemannPoints
      real   (r_kind), dimension(:,:,:), allocatable :: G13B   ! G13B on RiemannPoints
      real   (r_kind), dimension(:,:,:), allocatable :: G13T   ! G13T on RiemannPoints
      
      real   (r_kind), dimension(nPointsOnEdge              ) :: quad_pos_1d
      real   (r_kind), dimension(nPointsOnEdge              ) :: quad_wts_1d
      real   (r_kind), dimension(nPointsOnEdge,nPointsOnEdge) :: quad_wts_2d
      
      ! For dzdeta ( sqrt(G) )
      real(r_kind), dimension(:,:,:), allocatable :: dzdxi
      real(r_kind), dimension(:,:,:), allocatable :: dxideta
      real(r_kind), dimension(:,:,:), allocatable :: dzdeta
      
      ! For detadx ( G13 )
      real(r_kind), dimension(:,:,:), allocatable :: detadxi
      real(r_kind), dimension(:,:,:), allocatable :: dxidx
      real(r_kind), dimension(:,:,:), allocatable :: dzsdx
      real(r_kind), dimension(:,:,:), allocatable :: detadx
    contains
    
    subroutine init_mesh
      integer(i_kind) :: i,j,k
      integer(i_kind) :: iQP,jQP
      
      real(r_kind), dimension(nPointsOnEdge,1) :: quad_wts_tmp
      
      quad_pos_1d = (/-0.90617984593866399279762687829939,&
                      -0.53846931010568309103631442070021,&
                       0.                                ,&
                       0.53846931010568309103631442070021,&
                       0.90617984593866399279762687829939/)
      quad_pos_1d = ( quad_pos_1d + 1. ) / 2.
      
      quad_wts_1d = (/0.23692688505618908751426404071992,&
	                  0.47862867049936646804129151483564,&
	                  0.56888888888888888888888888888889,&
	                  0.47862867049936646804129151483564,&
	                  0.23692688505618908751426404071992/)
      quad_wts_1d  = quad_wts_1d / 2.
      
      quad_wts_tmp(:,1) = quad_wts_1d
      
      quad_wts_2d = matmul(quad_wts_tmp,transpose(quad_wts_tmp))
      
      allocate(xCorner  (nEdgesOnCell,ics:ice,kcs:kce))
      allocate(etaCorner(nEdgesOnCell,ics:ice,kcs:kce))
      allocate(zCorner  (nEdgesOnCell,ics:ice,kcs:kce))
      
      allocate(xQP    (nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(etaQP  (nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(xiQP   (nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(zQP    (nQuadPointsOnCell,ics:ice,kcs:kce))
      
      allocate(zsQP   (nQuadPointsOnCell,ics:ice,kcs:kce)) ! topography on Quadrature Points
      
      allocate(sqrtGQP(nQuadPointsOnCell,ics:ice,kcs:kce)) ! sqrtG on Quadrature Points
      allocate(G13QP  (nQuadPointsOnCell,ics:ice,kcs:kce)) ! G13 on Quadrature Points
      
      allocate(xL(nPointsOnEdge,ics:ice,kcs:kce))
      allocate(xR(nPointsOnEdge,ics:ice,kcs:kce))
      allocate(xB(nPointsOnEdge,ics:ice,kcs:kce))
      allocate(xT(nPointsOnEdge,ics:ice,kcs:kce))
      
      allocate(zL(nPointsOnEdge,ics:ice,kcs:kce))
      allocate(zR(nPointsOnEdge,ics:ice,kcs:kce))
      allocate(zB(nPointsOnEdge,ics:ice,kcs:kce))
      allocate(zT(nPointsOnEdge,ics:ice,kcs:kce))
      
      allocate(sqrtGL(nPointsOnEdge,ics:ice,kcs:kce)) ! sqrtGL on RiemannPoints
      allocate(sqrtGR(nPointsOnEdge,ics:ice,kcs:kce)) ! sqrtGR on RiemannPoints
      allocate(sqrtGB(nPointsOnEdge,ics:ice,kcs:kce)) ! sqrtGB on RiemannPoints
      allocate(sqrtGT(nPointsOnEdge,ics:ice,kcs:kce)) ! sqrtGT on RiemannPoints
      
      allocate(G13L(nPointsOnEdge,ics:ice,kcs:kce))   ! G13L on RiemannPoints
      allocate(G13R(nPointsOnEdge,ics:ice,kcs:kce))   ! G13R on RiemannPoints
      allocate(G13B(nPointsOnEdge,ics:ice,kcs:kce))   ! G13B on RiemannPoints
      allocate(G13T(nPointsOnEdge,ics:ice,kcs:kce))   ! G13T on RiemannPoints
      
      ! For dzdeta ( sqrt(G) )
      allocate(dzdxi  (nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(dxideta(nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(dzdeta (nQuadPointsOnCell,ics:ice,kcs:kce))
      
      ! For detadx ( G13 )
      allocate(detadxi(nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(dxidx  (nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(dzsdx  (nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(detadx (nQuadPointsOnCell,ics:ice,kcs:kce))
      
      deta = 1. / real(nz,r_kind)
      dxi  = ( z_max - z_min ) / nz
      
      dxideta = z_max - z_min
      detadxi = 1. / dxideta
      
      do k = kcs,kce
        do i = ics,ice
          xCorner(1,i,k) = x_min + ( real(i) - 1. ) * dx
          xCorner(2,i,k) = xCorner(1,i,k) + dx
          xCorner(3,i,k) = xCorner(2,i,k)
          xCorner(4,i,k) = xCorner(1,i,k)
          
          etaCorner(1,i,k) = ( real(k) - 1. ) * deta
          etaCorner(2,i,k) = etaCorner(1,i,k)
          etaCorner(3,i,k) = etaCorner(1,i,k) + deta
          etaCorner(4,i,k) = etaCorner(3,i,k)
          
          j = 0
          do jQP = 1,nPointsOnEdge
            do iQP = 1,nPointsOnEdge
              j = j + 1
              xQP  (j,i,k) = xCorner  (1,i,k) + dx   * quad_pos_1d(iQP)
              etaQP(j,i,k) = etaCorner(1,i,k) + deta * quad_pos_1d(jQP)
              xiQP (j,i,k) = etaQP(j,i,k) * detadxi(j,i,k)
            enddo
          enddo
        enddo
      enddo
      
    end subroutine init_mesh
    
    subroutine init_vertical_distribiution
      integer(r_kind) i,j,k
    
      ! For Schar, 2002
      real(r_kind) :: H
      real(r_kind) :: s
      
      if( vertical_coordinate == 1 )then
        ! Gal-Chen, 1975
        do k = kcs,kce
          do i = ics,ice
            do j = 1,nQuadPointsOnCell
              dzdxi(j,i,k) = ( z_max - zsQP(j,i,k) ) / z_max
              dxidx(j,i,k) = ( xiQP(j,i,k) - z_max ) / ( z_max - zsQP(j,i,k) ) * dzsdx(j,i,k)
              
              dzdeta(j,i,k) = dzdxi  (j,i,k) * dxideta(j,i,k)
              detadx(j,i,k) = detadxi(j,i,k) * dxidx  (j,i,k)
              
              zQP(j,i,k) = ( z_max - zsQP(j,i,k) ) / z_max * xiQP(j,i,k)  + zsQP(j,i,k)
            enddo
          enddo
        enddo
      elseif( vertical_coordinate == 2 )then
        ! Schar, 2002
        H = z_max - z_min
        if(case_num==2)then
          s = 3000.
        else
          s = H / exp(1.)
        endif
        
        do k = kcs,kce
          do i = ics,ice
            do j = 1,nQuadPointsOnCell
              dzdxi(j,i,k) = 1. - zsQP(j,i,k) * cosh( ( H - xiQP(j,i,k) ) / s ) / ( s * sinh( H / s ) )
              dxidx(j,i,k) = -sinh( ( z_max - xiQP(j,i,k) ) / s ) / sinh( z_max / s ) &
                           / ( 1. - zsQP(j,i,k) / s * cosh( ( z_max - xiQP(j,i,k) ) / s ) / sinh( z_max / s ) ) * dzsdx(j,i,k)
              
              dzdeta(j,i,k) = dzdxi  (j,i,k) * dxideta(j,i,k)
              detadx(j,i,k) = detadxi(j,i,k) * dxidx  (j,i,k)
              
              zQP(j,i,k) = xiQP(j,i,k) + zsQP(j,i,k) * sinh( ( H - xiQP(j,i,k) ) / s ) / sinh( H / s )
            enddo
          enddo
        enddo
      endif
      
    end subroutine init_vertical_distribiution
      
    end module mesh_mod
    
