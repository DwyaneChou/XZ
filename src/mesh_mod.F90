    module mesh_mod
      use constants_mod
      use parameters_mod
      use reconstruction_mod
      use math_mod
      implicit none
      
      real(r_kind) deta
      real(r_kind) dxi
      
      real(r_kind), dimension(:,:,:), allocatable :: xCorner   ! 1: low left, 2: low right, 3: up right, 4:  up left
      real(r_kind), dimension(:,:,:), allocatable :: etaCorner ! 1: low left, 2: low right, 3: up right, 4:  up left
      real(r_kind), dimension(:,:,:), allocatable :: zCorner   ! 1: low left, 2: low right, 3: up right, 4:  up left
      
      real(r_kind), dimension(:,:,:), allocatable :: x   ! x coordinate of Quadrature Points
      real(r_kind), dimension(:,:,:), allocatable :: eta ! eta coordinate of Quadrature Points
      real(r_kind), dimension(:,:,:), allocatable :: xi  ! xi coordinate of Quadrature Points
      real(r_kind), dimension(:,:,:), allocatable :: z   ! z coordinate of Quadrature Points
      
      real   (r_kind), dimension(:), allocatable :: xL
      real   (r_kind), dimension(:), allocatable :: xR
      real   (r_kind), dimension(:), allocatable :: xB
      real   (r_kind), dimension(:), allocatable :: xT
      
      real   (r_kind), dimension(:), allocatable :: etaL
      real   (r_kind), dimension(:), allocatable :: etaR
      real   (r_kind), dimension(:), allocatable :: etaB
      real   (r_kind), dimension(:), allocatable :: etaT
      
      real(r_kind), dimension(:,:,:), allocatable :: zs ! topography on Quadrature Points
      
      real(r_kind), dimension(:,:,:), allocatable :: sqrtG ! sqrtG on Quadrature Points
      real(r_kind), dimension(:,:,:), allocatable :: G13   ! G13 on Quadrature Points
      
      real(r_kind), dimension(  :,:), allocatable :: sqrtGC ! sqrtG on Cell
      real(r_kind), dimension(:,:,:), allocatable :: sqrtGL ! sqrtG on RiemannPoints on Left edge of Cell
      real(r_kind), dimension(:,:,:), allocatable :: sqrtGR ! sqrtG on RiemannPoints on Right edge of Cell
      real(r_kind), dimension(:,:,:), allocatable :: sqrtGB ! sqrtG on RiemannPoints on Bottom edge of Cell
      real(r_kind), dimension(:,:,:), allocatable :: sqrtGT ! sqrtG on RiemannPoints on Top edge of Cell
      
      real(r_kind), dimension(  :,:), allocatable :: G13C   ! G13 on Cell
      real(r_kind), dimension(:,:,:), allocatable :: G13L   ! G13 on RiemannPoints on Left edge of Cell
      real(r_kind), dimension(:,:,:), allocatable :: G13R   ! G13 on RiemannPoints on Right edge of Cell
      real(r_kind), dimension(:,:,:), allocatable :: G13B   ! G13 on RiemannPoints on Bottom edge of Cell
      real(r_kind), dimension(:,:,:), allocatable :: G13T   ! G13 on RiemannPoints on Top edge of Cell
      
      real(r_kind), dimension(:,:), allocatable :: xC   ! x coordinate on Cell by Gaussian quadrature
      real(r_kind), dimension(:,:), allocatable :: zC   ! z coordinate on Cell by Gaussian quadrature
      real(r_kind), dimension(:,:), allocatable :: zsC ! topography on Cell by Gaussian quadrature
        
      real(r_kind), dimension(:,:), allocatable :: xCenter   ! center coordinate
      real(r_kind), dimension(:,:), allocatable :: etaCenter ! center coordinate
      
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
        
        print*,'Initialize mesh'
        print*,''
        
        allocate(xCorner  (nEdgesOnCell,ics:ice,kcs:kce))
        allocate(etaCorner(nEdgesOnCell,ics:ice,kcs:kce))
        allocate(zCorner  (nEdgesOnCell,ics:ice,kcs:kce))
        
        allocate(x  (nQuadPointsOnCell,ics:ice,kcs:kce))
        allocate(eta(nQuadPointsOnCell,ics:ice,kcs:kce))
        allocate(xi (nQuadPointsOnCell,ics:ice,kcs:kce))
        allocate(z  (nQuadPointsOnCell,ics:ice,kcs:kce))
        
        allocate(xL(nPointsOnEdge))
        allocate(xR(nPointsOnEdge))
        allocate(xB(nPointsOnEdge))
        allocate(xT(nPointsOnEdge))
        
        allocate(etaL(nPointsOnEdge))
        allocate(etaR(nPointsOnEdge))
        allocate(etaB(nPointsOnEdge))
        allocate(etaT(nPointsOnEdge))
        
        allocate(zs   (nQuadPointsOnCell,ics:ice,kcs:kce)) ! topography on Quadrature Points
        
        allocate(sqrtG(nQuadPointsOnCell,ics:ice,kcs:kce)) ! sqrtG on Quadrature Points
        allocate(G13  (nQuadPointsOnCell,ics:ice,kcs:kce)) ! G13 on Quadrature Points
        
        allocate(sqrtGC(              ics:ice,kcs:kce)) ! sqrtG on Cell
        allocate(sqrtGL(nPointsOnEdge,ics:ice,kcs:kce)) ! sqrtG on RiemannPoints on Left edge of Cell
        allocate(sqrtGR(nPointsOnEdge,ics:ice,kcs:kce)) ! sqrtG on RiemannPoints on Right edge of Cell
        allocate(sqrtGB(nPointsOnEdge,ics:ice,kcs:kce)) ! sqrtG on RiemannPoints on Bottom edge of Cell
        allocate(sqrtGT(nPointsOnEdge,ics:ice,kcs:kce)) ! sqrtG on RiemannPoints on Top edge of Cell
        
        allocate(G13C(              ics:ice,kcs:kce))   ! G13 on Cell
        allocate(G13L(nPointsOnEdge,ics:ice,kcs:kce))   ! G13 on RiemannPoints on Left edge of Cell
        allocate(G13R(nPointsOnEdge,ics:ice,kcs:kce))   ! G13 on RiemannPoints on Right edge of Cell
        allocate(G13B(nPointsOnEdge,ics:ice,kcs:kce))   ! G13 on RiemannPoints on Bottom edge of Cell
        allocate(G13T(nPointsOnEdge,ics:ice,kcs:kce))   ! G13 on RiemannPoints on Top edge of Cell
        
        allocate(xC (ics:ice,kcs:kce))
        allocate(zC (ics:ice,kcs:kce))
        allocate(zsC(ics:ice,kcs:kce))
        
        allocate(xCenter  (ics:ice,kcs:kce))
        allocate(etaCenter(ics:ice,kcs:kce))
        
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
                x  (j,i,k) = xCorner  (1,i,k) + dx   * quad_pos_1d(iQP)
                eta(j,i,k) = etaCorner(1,i,k) + deta * quad_pos_1d(jQP)
                xi (j,i,k) = eta(j,i,k) * dxideta(j,i,k)
              enddo
            enddo
            xC(i,k) = Gaussian_quadrature_2d(x(:,i,k))
          enddo
        enddo
        
        do k = kds,kde
          do i = ids,ide
            xCenter   = ( xCorner  (1,i,k) + xCorner  (2,i,k) ) / 2.
            etaCenter = ( etaCorner(1,i,k) + etaCorner(4,i,k) ) / 2.
          enddo
        enddo
        
        xL = - dx / 2.
        xR =   dx / 2.
        do iQP = 1,nPointsOnEdge
          xB(iQP) = xL(1) + dx * quad_pos_1d(iQP)
          xT(iQP) = xB(iQP)
        enddo
        
        etaB = - deta / 2.
        etaT =   deta / 2.
        do jQP = 1,nPointsOnEdge
          etaL(jQP) = etaB(1) + deta * quad_pos_1d(jQP)
          etaR(jQP) = etaL(jQP)
        enddo
        
      end subroutine init_mesh
      
      subroutine init_vertical_coordinate
        integer(r_kind) i,j,k
      
        ! For Schar, 2002
        real(r_kind) :: H
        real(r_kind) :: s
        
        print*,'Initialize vertical coordinate'
        print*,''
        
        if( vertical_coordinate == 1 )then
          ! Gal-Chen, 1975
          do k = kcs,kce
            do i = ics,ice
              do j = 1,nQuadPointsOnCell
                dzdxi(j,i,k) = ( z_max - zs(j,i,k) ) / z_max
                dxidx(j,i,k) = ( xi(j,i,k) - z_max ) / ( z_max - zs(j,i,k) ) * dzsdx(j,i,k)
                
                dzdeta(j,i,k) = dzdxi  (j,i,k) * dxideta(j,i,k)
                detadx(j,i,k) = detadxi(j,i,k) * dxidx  (j,i,k)
                
                z(j,i,k) = ( z_max - zs(j,i,k) ) / z_max * xi(j,i,k)  + zs(j,i,k)
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
                dzdxi(j,i,k) = 1. - zs(j,i,k) * cosh( ( H - xi(j,i,k) ) / s ) / ( s * sinh( H / s ) )
                dxidx(j,i,k) = -sinh( ( z_max - xi(j,i,k) ) / s ) / sinh( z_max / s ) &
                             / ( 1. - zs(j,i,k) / s * cosh( ( z_max - xi(j,i,k) ) / s ) / sinh( z_max / s ) ) * dzsdx(j,i,k)
                
                dzdeta(j,i,k) = dzdxi  (j,i,k) * dxideta(j,i,k)
                detadx(j,i,k) = detadxi(j,i,k) * dxidx  (j,i,k)
                
                z(j,i,k) = xi(j,i,k) + zs(j,i,k) * sinh( ( H - xi(j,i,k) ) / s ) / sinh( H / s )
              enddo
            enddo
          enddo
        endif
        
        ! Reconstruct mertric tensor by analytical value
        do k = kcs,kce
          do i = ics,ice
            do j = 1,nQuadPointsOnCell
              sqrtG(j,i,k) = dzdeta(j,i,k)
              G13  (j,i,k) = detadx(j,i,k)
            enddo
          enddo
        enddo
        
        ! Calculate cell value
        do k = kcs,kce
          do i = ics,ice
            zC    (i,k) = Gaussian_quadrature_2d(z    (:,i,k))
            zsC   (i,k) = Gaussian_quadrature_2d(zs   (:,i,k))
            sqrtGC(i,k) = Gaussian_quadrature_2d(sqrtG(:,i,k))
            G13C  (i,k) = Gaussian_quadrature_2d(G13  (:,i,k))
          enddo
        enddo
        
      end subroutine init_vertical_coordinate
      
    end module mesh_mod
    
