    module reconstruction_mod
      use constants_mod
      use parameters_mod
      use qr_solver_mod
      implicit none
      
      real   (r_kind), dimension(nPointsOnEdge    ) :: quad_pos_1d
      real   (r_kind), dimension(nPointsOnEdge    ) :: quad_wts_1d
      real   (r_kind), dimension(nQuadPointsOnCell) :: quad_wts_2d
      
      integer(i_kind) :: nRecCells
      
      integer(i_kind) :: nRecTerms
      
      real   (r_kind), dimension(:,:,:,:), allocatable :: polyCoordCoef
      
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
        
        ! Check known and unknown values
        if(nRecCells<nRecTerms)then
          print*,'Error! nRecCells<nRecTerms, nRecCells =',nRecCells,' nRecTerm =',nRecTerms
          print*,'Choose larger stencil_width or smaller recPolyDegree'
          stop
        endif
        
        allocate( polyCoordCoef(nRecCells,nRecTerms,ids:ide,kds:kde) )
        
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
    
      function WLS_ENO(A,u,h,m,n,ic)
        ! WLS_ENO
        ! Ax = b -> WAx=Wb -> min|| WAx - Wb || -> x
        ! where A->A, x->WLS_ENO, b->u
        ! WLS_ENO : unknown vector x (the vector of reconstruction polynomial coefficients)
        ! A       : matrix of coordinate coefficient
        ! u       : vector of known values
        ! h       : vector of distence between adjacent cells and center cell
        ! m       : number of known values (volumn integration value on cell, or in other words, number of cells in a stencil)
        ! n       : number of coeficients of reconstruction polynomial
        ! ic      : index of center cell on stencil
        real   (r_kind), dimension(n  )             :: WLS_ENO
        integer(i_kind)                , intent(in) :: m
        integer(i_kind)                , intent(in) :: n
        real   (r_kind), dimension(m,n), intent(in) :: A
        real   (r_kind), dimension(m  ), intent(in) :: u
        real   (r_kind), dimension(m  ), intent(in) :: h
        integer(i_kind)                , intent(in) :: ic
        
        real(r_kind),parameter :: alpha   = 1.5
        real(r_kind),parameter :: epsilon = 1.e-10
        
        real(r_kind), dimension(m,n) :: WA
        real(r_kind), dimension(m  ) :: Wu
        real(r_kind), dimension(m  ) :: W ! weights on each cells
        real(r_kind), dimension(m  ) :: beta
        
        integer(i_kind) :: i,j,k
        
        do j = 1,m
          beta(j) = ( u(j) - u(ic) )**2 + epsilon * h(j)**2
        enddo
        beta(ic) = minval(beta,abs(beta)>1.e-15)
        
        W = 1./beta
        W(ic) = alpha * W(ic)
        
        do j = 1,m
          WA(j,:) = W(j) * A(j,:)
          Wu(j  ) = W(j) * u(j  )
        enddo
        
        call qr_solver(WA,Wu,WLS_ENO,M,N)
      end function WLS_ENO
      
    end module reconstruction_mod
    