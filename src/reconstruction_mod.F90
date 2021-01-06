    module reconstruction_mod
      use constants_mod
      use parameters_mod
      implicit none
      
      real   (r_kind), dimension(nPointsOnEdge    ) :: quad_pos_1d
      real   (r_kind), dimension(nPointsOnEdge    ) :: quad_wts_1d
      real   (r_kind), dimension(nQuadPointsOnCell) :: quad_wts_2d
      
      integer(i_kind) :: nRecCells
      
      integer(i_kind) :: nRecTerms
      
      real   (r_kind), dimension(:,:,:), allocatable :: polyCoordCoef
      
    contains
      subroutine init_reconstruction
        integer(i_kind) :: i,j,k
        
        real(r_kind), dimension(nPointsOnEdge,1            ) :: quad_wts_tmp_1d
        real(r_kind), dimension(nPointsOnEdge,nPointsOnEdge) :: quad_wts_tmp_2d
        
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
        
        quad_wts_tmp_1d(:,1) = quad_wts_1d
        
        quad_wts_tmp_2d = matmul(quad_wts_tmp_1d,transpose(quad_wts_tmp_1d))
        
        k = 0
        do j = 1,nPointsOnEdge
          do i = 1,nPointsOnEdge
            k = k + 1
            quad_wts_2d(k) = quad_wts_tmp_2d(i,j)
          enddo
        enddo
        
        nRecCells = stencil_width**2
        nRecTerms = (recPolyDegree+1) * (recPolyDegree+2) / 2
        
        allocate( polyCoordCoef(nRecTerms,ics:ice,kcs:kce) )
        
      end subroutine init_reconstruction
      
      function Gaussian_quadrature_1d(q)
        real(r_kind) :: Gaussian_quadrature_1d
        real(r_kind) :: q(nPointsOnEdge)
        
        Gaussian_quadrature_1d = dot_product(quad_wts_1d,q)
        
      end function Gaussian_quadrature_1d
      
      function Gaussian_quadrature_2d(q)
        real(r_kind) :: Gaussian_quadrature_2d
        real(r_kind) :: q(nQuadPointsOnCell)
        
        Gaussian_quadrature_2d = dot_product(quad_wts_2d,q)
        
      end function Gaussian_quadrature_2d
    
    end module reconstruction_mod
    