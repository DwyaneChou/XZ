    module mesh_mod
      use constants_mod
      use parameters_mod
      use reconstruction_mod
      use math_mod
      implicit none
      
      real(r_kind) deta
      
      real(r_kind), dimension(:,:,:), allocatable :: xCorner   ! 1: low left, 2: low right, 3: up right, 4:  up left
      real(r_kind), dimension(:,:,:), allocatable :: etaCorner ! 1: low left, 2: low right, 3: up right, 4:  up left
      real(r_kind), dimension(:,:,:), allocatable :: zCorner   ! 1: low left, 2: low right, 3: up right, 4:  up left
      
      real(r_kind), dimension(:,:,:), allocatable :: x   ! x coordinate of Points on cell
      real(r_kind), dimension(:,:,:), allocatable :: eta ! eta coordinate of Points on cell
      real(r_kind), dimension(:,:,:), allocatable :: z   ! z coordinate of Points on cell
      
      real   (r_kind), dimension(:), allocatable :: xL
      real   (r_kind), dimension(:), allocatable :: xR
      real   (r_kind), dimension(:), allocatable :: xB
      real   (r_kind), dimension(:), allocatable :: xT
      
      real   (r_kind), dimension(:), allocatable :: etaL
      real   (r_kind), dimension(:), allocatable :: etaR
      real   (r_kind), dimension(:), allocatable :: etaB
      real   (r_kind), dimension(:), allocatable :: etaT
      
      real(r_kind), dimension(:,:,:), allocatable :: zs ! topography on Points on cell
      
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
      real(r_kind), dimension(:,:,:), allocatable :: dzdeta
      
      ! For detadx ( G13 )
      real(r_kind), dimension(:,:,:), allocatable :: dzsdx
      real(r_kind), dimension(:,:,:), allocatable :: detadx
      
      ! For calculate analytical metric term
      real(r_kind), dimension(:,:,:,:), allocatable :: x_ext     ! (nPointsOnEdge,nEdgesOnCell,ids:ide,kds:kde) &
      real(r_kind), dimension(:,:,:,:), allocatable :: eta_ext   ! edge index: 1 for bottom
      real(r_kind), dimension(:,:,:,:), allocatable :: z_ext     !             2 for right
      real(r_kind), dimension(:,:,:,:), allocatable :: zs_ext    !             3 for top
      real(r_kind), dimension(:,:,:,:), allocatable :: sqrtG_ext !             4 for left
      real(r_kind), dimension(:,:,:,:), allocatable :: G13_ext
      real(r_kind), dimension(:,:,:,:), allocatable :: dzdeta_ext
      real(r_kind), dimension(:,:,:,:), allocatable :: dzsdx_ext
      real(r_kind), dimension(:,:,:,:), allocatable :: detadx_ext
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
        allocate(dzdeta (nQuadPointsOnCell,ics:ice,kcs:kce))
        
        ! For detadx ( G13 )
        allocate(dzsdx  (nQuadPointsOnCell,ics:ice,kcs:kce))
        allocate(detadx (nQuadPointsOnCell,ics:ice,kcs:kce))
        
        ! For calculate analytical metric term
        allocate(x_ext     (nPointsOnEdge,nEdgesOnCell,ids:ide,kds:kde))
        allocate(eta_ext   (nPointsOnEdge,nEdgesOnCell,ids:ide,kds:kde))
        allocate(z_ext     (nPointsOnEdge,nEdgesOnCell,ids:ide,kds:kde))
        allocate(zs_ext    (nPointsOnEdge,nEdgesOnCell,ids:ide,kds:kde))
        allocate(sqrtG_ext (nPointsOnEdge,nEdgesOnCell,ids:ide,kds:kde))
        allocate(G13_ext   (nPointsOnEdge,nEdgesOnCell,ids:ide,kds:kde))
        allocate(dzdeta_ext(nPointsOnEdge,nEdgesOnCell,ids:ide,kds:kde))
        allocate(dzsdx_ext (nPointsOnEdge,nEdgesOnCell,ids:ide,kds:kde))
        allocate(detadx_ext(nPointsOnEdge,nEdgesOnCell,ids:ide,kds:kde))
        
        deta = ( z_max - z_min ) / nz
        
        do k = kcs,kce
          do i = ics,ice
            xCorner(1,i,k) = x_min + ( real(i) - 1. ) * dx
            xCorner(2,i,k) = xCorner(1,i,k) + dx
            xCorner(3,i,k) = xCorner(2,i,k)
            xCorner(4,i,k) = xCorner(1,i,k)
            
            etaCorner(1,i,k) = z_min + ( real(k) - 1. ) * deta
            etaCorner(2,i,k) = etaCorner(1,i,k)
            etaCorner(3,i,k) = etaCorner(1,i,k) + deta
            etaCorner(4,i,k) = etaCorner(3,i,k)
            
            j = 0
            do jQP = 1,nPointsOnEdge
              do iQP = 1,nPointsOnEdge
                j = j + 1
                x  (j,i,k) = xCorner  (1,i,k) + dx   * quad_pos_1d(iQP)
                eta(j,i,k) = etaCorner(1,i,k) + deta * quad_pos_1d(jQP)
              enddo
            enddo
            xC(i,k) = Gaussian_quadrature_2d(x(:,i,k))
          enddo
        enddo
        
        do k = kds,kde
          do i = ids,ide
            xCenter  (i,k) = ( xCorner  (1,i,k) + xCorner  (2,i,k) ) / 2.
            etaCenter(i,k) = ( etaCorner(1,i,k) + etaCorner(4,i,k) ) / 2.
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
        
        do k = kds,kde
          do i = ids,ide
            x_ext  (:,1,i,k) = xCenter  (i,k) + xB
            x_ext  (:,2,i,k) = xCenter  (i,k) + xR
            x_ext  (:,3,i,k) = xCenter  (i,k) + xT
            x_ext  (:,4,i,k) = xCenter  (i,k) + xL
            eta_ext(:,1,i,k) = etaCenter(i,k) + etaB
            eta_ext(:,2,i,k) = etaCenter(i,k) + etaR
            eta_ext(:,3,i,k) = etaCenter(i,k) + etaT
            eta_ext(:,4,i,k) = etaCenter(i,k) + etaL
          enddo
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
          dzdeta = ( z_max - zs ) / z_max
          detadx = ( eta - z_max ) / ( z_max - zs ) * dzsdx
          z      = ( z_max - zs ) / z_max * eta + zs
          
          dzdeta_ext = ( z_max - zs_ext ) / z_max
          detadx_ext = ( eta_ext - z_max ) / ( z_max - zs_ext ) * dzsdx_ext
          z_ext      = ( z_max - zs_ext ) / z_max * eta_ext + zs_ext
        elseif( vertical_coordinate == 2 )then
          ! Schar, 2002
          H = z_max - z_min
          if(case_num==2)then
            s = 3000.
          else
            s = H / exp(1.)
          endif
          
          dzdeta = 1. - zs * cosh( ( H - eta ) / s ) / ( s * sinh( H / s ) )
          detadx = -dzsdx*sinh((H-eta)/s)/sinh(H/s) / dzdeta
          z      = eta + zs * sinh( ( H - eta ) / s ) / sinh( H / s )
          
          dzdeta_ext = 1. - zs_ext * cosh( ( H - eta_ext ) / s ) / ( s * sinh( H / s ) )
          detadx_ext = -dzsdx_ext*sinh((H-eta_ext)/s)/sinh(H/s) / dzdeta_ext
          z_ext      = eta_ext + zs_ext * sinh( ( H - eta_ext ) / s ) / sinh( H / s )
        endif
        
        ! Reconstruct mertric tensor by analytical value
        sqrtG = dzdeta
        G13   = detadx
        
        sqrtG_ext = dzdeta_ext
        G13_ext   = detadx_ext
        
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
    
