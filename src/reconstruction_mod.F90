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
      
      real   (r_kind) :: recCoef
      real   (r_kind) :: recdx
      real   (r_kind) :: recdeta
      real   (r_kind) :: recdV
      
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
        ! h       : average edge length of cells
        ! m       : number of known values (volumn integration value on cell, or in other words, number of cells in a stencil)
        ! n       : number of coeficients of reconstruction polynomial
        ! ic      : index of center cell on stencil
        real   (r_kind), dimension(n  )             :: WLS_ENO
        integer(i_kind)                , intent(in) :: m
        integer(i_kind)                , intent(in) :: n
        real   (r_kind), dimension(m,n), intent(in) :: A
        real   (r_kind), dimension(m  ), intent(in) :: u
        real   (r_kind)                , intent(in) :: h
        integer(i_kind)                , intent(in) :: ic
        
        real(r_kind),parameter :: alpha   = 10
        real(r_kind),parameter :: epsilon = 1.e-2
        
        real(r_kind), dimension(m,n) :: WA
        real(r_kind), dimension(m  ) :: Wu
        real(r_kind), dimension(m  ) :: W ! weights on each cells
        real(r_kind), dimension(m  ) :: beta
        
        integer(i_kind) :: i,j,k
        
        do j = 1,m
          beta(j) = ( u(j) - u(ic) )**2 + epsilon * h*h
        enddo
        beta(ic) = minval(beta,abs(beta)>1.e-15)
        
        W = 1./beta
        W(ic) = alpha * W(ic)
        
        do j = 1,m
          WA(j,:) = W(j) * A(j,:)
          Wu(j  ) = W(j) * u(j  )
        enddo
        
        call qr_solver(WA,Wu,WLS_ENO,M,N)
        
        !if(abs(WLS_ENO(1)-1.)>1.e-15)then
        !  print*,W
        !  print*,''
        !  !print*,u
        !  !print*,''
        !  print*,A
        !  print*,''
        !  print*,WLS_ENO
        !endif
        
      end function WLS_ENO
      
      subroutine WENO_limiter(Qrec,Q,dir)
        real   (r_kind)              , intent(out) :: Qrec
        real   (r_kind), dimension(5), intent(in ) :: Q
        integer(i_kind)              , intent(in ) :: dir
        
        integer(i_kind), parameter :: nStencil = 3
        real   (r_kind), parameter :: weno_coef(3)  = [0.1, 0.6, 0.3]
        real   (r_kind), parameter :: eps           = 1.E-16
        
        real(r_kind) Qim(nStencil-1)
        real(r_kind) Qip(nStencil-1)
        real(r_kind) Qi
        
        integer(i_kind) iStencil
        
        Qim(2) = Q(1)
        Qim(1) = Q(2)
        Qi     = Q(3)
        Qip(1) = Q(4)
        Qip(2) = Q(5)
        
        if( .not. any(Q==FillValue) )then
          call WENO5(Qrec,Q,dir)
        elseif( any(Q==FillValue) .and. Qi/=FillValue )then
          call WENO3(Qrec,Q(2:4),dir)
        else
          print*,'Center point is FillValue, WENO cannot be implemented.'
        endif
        
      end subroutine WENO_limiter
      
      ! 1D WENO 5th order slope limiter, according to Sun,2015
      ! "A Slope Constrained 4th OrderMulti-Moment Finite Volume Method with WENO Limiter"
      ! and Jiang and Shu, 1996
      subroutine WENO5(Qrec,Q,dir)
        real   (r_kind)              , intent(out) :: Qrec
        real   (r_kind), dimension(5), intent(in ) :: Q
        integer(i_kind)              , intent(in ) :: dir
        
        integer(i_kind), parameter :: nStencil = 3
        real   (r_kind), parameter :: weno_coef(3)  = [0.1, 0.6, 0.3]
        real   (r_kind), parameter :: eps           = 1.E-16
        
        real(r_kind) Qim(nStencil-1)
        real(r_kind) Qip(nStencil-1)
        real(r_kind) Qi
        
        real(r_kind), dimension(nStencil) :: stencil
        real(r_kind), dimension(nStencil) :: coefA
        real(r_kind), dimension(nStencil) :: coefB
        real(r_kind), dimension(nStencil) :: alpha
        real(r_kind), dimension(nStencil) :: beta
        real(r_kind), dimension(nStencil) :: omega
        
        real(r_kind) tau40
        real(r_kind) tau41
        real(r_kind) tau5
        
        integer(i_kind) iStencil
        
        Qim(2) = Q(1)
        Qim(1) = Q(2)
        Qi     = Q(3)
        Qip(1) = Q(4)
        Qip(2) = Q(5)
        
        if(dir>0)then
          stencil (1) =  Qim(2)/3. - 7./6. * Qim(1) + 11./6. * Qi     
          stencil (2) = -Qim(1)/6. + 5./6. * Qi     +  1./3. * Qip(1) 
          stencil (3) =  Qi    /3. + 5./6. * Qip(1) -  1./6. * Qip(2)
          
          coefA(1) = Qim(2) - 2. * Qim(1) + Qi
          coefA(2) = Qim(1) - 2. * Qi     + Qip(1)
          coefA(3) = Qi     - 2. * Qip(1) + Qip(2)
          
          coefB(1) =      Qim(2) - 4. * Qim(1) + 3. * Qi
          coefB(2) =      Qim(1) -      Qip(1)
          coefB(3) = 3. * Qi     - 4. * Qip(1) +      Qip(2)
        elseif(dir<0)then
          stencil (1) =  Qip(2)/3. - 7./6. * Qip(1) + 11./6. * Qi     
          stencil (2) = -Qip(1)/6. + 5./6. * Qi     +  1./3. * Qim(1) 
          stencil (3) =  Qi    /3. + 5./6. * Qim(1) -  1./6. * Qim(2)
          
          coefA(1) = Qip(2) - 2. * Qip(1) + Qi
          coefA(2) = Qip(1) - 2. * Qi     + Qim(1)
          coefA(3) = Qi     - 2. * Qim(1) + Qim(2)
          
          coefB(1) =      Qip(2) - 4. * Qip(1) + 3. * Qi
          coefB(2) =      Qip(1) -      Qim(1)
          coefB(3) = 3. * Qi     - 4. * Qim(1) +      Qim(2)
        endif
        
        beta = coefA**2 * 13. / 12. + coefB**2 * 0.25
        
        ! WENO-Z
        tau40 = abs( beta(1) - beta(2) )
        tau41 = abs( beta(2) - beta(3) )
        tau5  = abs( beta(3) - beta(1) )
        
        if( tau40<=minval(beta) .and. tau41>minval(beta) )then
          omega = [1./3.,2./3.,0.]
        elseif( tau40>minval(beta) .and. tau41<minval(beta) )then
          omega = [0.,2./3.,1./3.]
        else
          alpha = weno_coef * ( 1. + tau5 / ( eps + beta ) )
          omega = alpha / sum(alpha)
        endif
        ! WENO-Z
        
        !! origin WENO
        !alpha = weno_coef / ( eps + beta )**2
        !omega = alpha / sum(alpha)
        !! origin WENO
        
        Qrec = dot_product( stencil, omega )
        
      end subroutine WENO5
      
      ! 1D 3th order WENO limiter, accroding to 
      ! T. Lunnet et al., 2017, Monthly Weather Review
      subroutine WENO3(Qrec,Q,dir)
        real   (r_kind)              , intent(out) :: Qrec
        real   (r_kind), dimension(3), intent(in ) :: Q
        integer(i_kind)              , intent(in ) :: dir
        
        integer(i_kind), parameter :: nStencil = 2
        real   (r_kind), parameter :: weno_coef(2)  = [1./3., 2./3.]
        real   (r_kind), parameter :: eps           = 1.E-16
        
        real(r_kind) Qim(nStencil-1)
        real(r_kind) Qip(nStencil-1)
        real(r_kind) Qi
        
        real(r_kind), dimension(nStencil) :: stencil
        real(r_kind), dimension(nStencil) :: coefA
        real(r_kind), dimension(nStencil) :: coefB
        real(r_kind), dimension(nStencil) :: alpha
        real(r_kind), dimension(nStencil) :: beta
        real(r_kind), dimension(nStencil) :: omega
        
        integer(i_kind) iStencil
        
        Qim(1) = Q(1)
        Qi     = Q(2)
        Qip(1) = Q(3)
        
        if(dir>0)then
          stencil (1) = -Qim(1)/2. + 3./2. * Qi
          stencil (2) =  Qi / 2. + Qip(1) / 2.
          
          coefA(1) =-Qim(1) + Qi
          coefA(2) = Qip(1) - Qi
        elseif(dir<0)then
          stencil (1) = -Qip(1)/2. + 3./2. * Qi
          stencil (2) = Qi / 2. + Qim(1) / 2.
          
          coefA(1) =-Qip(1) + Qi
          coefA(2) = Qim(1) - Qi
        endif
        
        beta = coefA**2
        
        ! origin WENO
        alpha = weno_coef / ( eps + beta )**2
        omega = alpha / sum(alpha)
        ! origin WENO
        
        Qrec = dot_product( stencil, omega )
        
      end subroutine WENO3
      
      ! Reconstruct on the left most boundary
      function left_side_recon3(q)
        real(r_kind) :: left_side_recon3
        real(r_kind),dimension(3),intent(in) :: q
        
        real(r_kind) Qip(2)
        real(r_kind) Qi
        
        Qi     = Q(1)
        Qip(1) = Q(2)
        Qip(2) = Q(3)
        
        left_side_recon3 = Qip(2)/3. - 7./6. * Qip(1) + 11./6. * Qi
        !left_side_recon3 = 1.5 * Qi - 0.5 * Qip(1)
        !left_side_recon3 = 3./8. * Qip(2) - 5./4. * Qip(1) + 15./8. * Qi
      
      end function left_side_recon3
      
      ! Reconstruct on the right most boundary
      function right_side_recon3(q)
        real(r_kind) :: right_side_recon3
        real(r_kind),dimension(3),intent(in) :: q
        
        real(r_kind) Qim(2)
        real(r_kind) Qi
        
        Qim(2) = Q(1)
        Qim(1) = Q(2)
        Qi     = Q(3)
        
        right_side_recon3 = Qim(2)/3. - 7./6. * Qim(1) + 11./6. * Qi
        !right_side_recon3 = 1.5 * Qi - 0.5 * Qim(1)
        !right_side_recon3 = 3./8. * Qim(2) - 5./4. * Qim(1) + 15./8. * Qi
      
      end function right_side_recon3
    end module reconstruction_mod
    