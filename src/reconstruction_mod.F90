  module reconstruction_mod
    use constants_mod
    contains
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
    
    real(r_kind) function dqdx(q,dx)
      real(r_kind), dimension(4), intent(in ) :: q
      real(r_kind)              , intent(in ) :: dx
      real(r_kind) :: dq
      
      dq = ( q(1) / 12. - 5./4. * q(2) + 5./4. * q(3) - q(4) / 12.  ) / dx
      
      dqdx = dq
    
    end function dqdx
    
    real(r_kind) function dqdxC(q,dx)
      real(r_kind), dimension(2), intent(in ) :: q
      real(r_kind)              , intent(in ) :: dx
      real(r_kind) :: dq
      
      dqdxC = ( q(2) - q(1)  ) / dx
      
      dqdxC = dq
    
    end function dqdxC
    
    real(r_kind) function dqdxL(q,dx)
      real(r_kind), dimension(3), intent(in ) :: q
      real(r_kind)              , intent(in ) :: dx
      real(r_kind) :: dq
      
      dq = - ( 2. * q(1) - 3. * q(2) + q(3) ) / dx
      
      if(maxval(abs(q))/=0)then
        if( abs( dq/maxval(abs(q)) * dx ) < 1.e-14 ) dq = 0
      endif
      
      dqdxL = dq
    
    end function dqdxL
    
    real(r_kind) function dqdxR(q,dx)
      real(r_kind), dimension(3), intent(in ) :: q
      real(r_kind)              , intent(in ) :: dx
      real(r_kind) :: dq
      
      dq = ( 2. * q(3) - 3. * q(2) + q(1) ) / dx
    
      dqdxR = dq
    end function dqdxR
  end module reconstruction_mod
    