MODULE spatial_operators_mod
  use constants_mod
  use mesh_mod
  use parameters_mod
  use stat_mod
  use tend_mod
  
  implicit none
  
  private
  
  public init_spatial_operator, spatial_operator
  
  integer(i_kind), parameter :: extPts = 3
  
  integer(i_kind) :: pos
  integer(i_kind) :: neg
      
  integer(i_kind) :: ics
  integer(i_kind) :: ice
  integer(i_kind) :: kcs
  integer(i_kind) :: kce
  
  real(r_kind), dimension(:,:,:), allocatable :: q_ext ! Extended forecast variables
  real(r_kind), dimension(:,:,:), allocatable :: qL    ! Reconstructed q_(i-1/2,k)
  real(r_kind), dimension(:,:,:), allocatable :: qR    ! Reconstructed q_(i+1/2,k)
  real(r_kind), dimension(:,:,:), allocatable :: qB    ! Reconstructed q_(i,k-1/2)
  real(r_kind), dimension(:,:,:), allocatable :: qT    ! Reconstructed q_(i,k+1/2)
  
  real(r_kind), dimension(:,:,:), allocatable :: FL    ! Reconstructed F_(i-1/2,k)
  real(r_kind), dimension(:,:,:), allocatable :: FR    ! Reconstructed F_(i+1/2,k)
  real(r_kind), dimension(:,:,:), allocatable :: FB    ! Reconstructed F_(i,k-1/2)
  real(r_kind), dimension(:,:,:), allocatable :: FT    ! Reconstructed F_(i,k+1/2)
      
  real(r_kind), dimension(:,:,:), allocatable :: HL    ! Reconstructed H_(i-1/2,k)
  real(r_kind), dimension(:,:,:), allocatable :: HR    ! Reconstructed H_(i+1/2,k)
  real(r_kind), dimension(:,:,:), allocatable :: HB    ! Reconstructed H_(i,k-1/2)
  real(r_kind), dimension(:,:,:), allocatable :: HT    ! Reconstructed H_(i,k+1/2)
  
  real(r_kind), dimension(:,:  ), allocatable :: PL    ! Reconstructed P_(i-1/2,k)
  real(r_kind), dimension(:,:  ), allocatable :: PR    ! Reconstructed P_(i+1/2,k)
  real(r_kind), dimension(:,:  ), allocatable :: PB    ! Reconstructed P_(i,k-1/2)
  real(r_kind), dimension(:,:  ), allocatable :: PT    ! Reconstructed P_(i,k+1/2)
  
    contains
    subroutine init_spatial_operator
      pos = 1
      neg = -1
      
      ics = ids - extPts
      ice = ide + extPts
      kcs = kds - extPts
      kce = kde + extPts
      
      allocate(q_ext(nVar,ics:ice,kcs:kce))
      allocate(qL   (nVar,ids:ide,kds:kde))
      allocate(qR   (nVar,ids:ide,kds:kde))
      allocate(qB   (nVar,ids:ide,kds:kde))
      allocate(qT   (nVar,ids:ide,kds:kde))
      
      allocate(FL   (nVar,ids:ide,kds:kde))
      allocate(FR   (nVar,ids:ide,kds:kde))
      allocate(FB   (nVar,ids:ide,kds:kde))
      allocate(FT   (nVar,ids:ide,kds:kde))
      
      allocate(HL   (nVar,ids:ide,kds:kde))
      allocate(HR   (nVar,ids:ide,kds:kde))
      allocate(HB   (nVar,ids:ide,kds:kde))
      allocate(HT   (nVar,ids:ide,kds:kde))
      
      allocate(PL   (     ids:ide,kds:kde))
      allocate(PR   (     ids:ide,kds:kde))
      allocate(PB   (     ids:ide,kds:kde))
      allocate(PT   (     ids:ide,kds:kde))
    end subroutine init_spatial_operator
    
    subroutine spatial_operator(stat,tend)
      type(stat_field), target, intent(in   ) :: stat
      type(tend_field), target, intent(inout) :: tend
      
      real(r_kind), dimension(5) :: q_weno
      
      character(30) :: bdy_type
      
      real(r_kind) eigenvalue_x(2)
      real(r_kind) eigenvalue_z(2)
      real(r_kind) maxeigen_x
      real(r_kind) maxeigen_z
      
      integer(i_kind) dir
      
      integer(i_kind) i,k,iVar
      
      integer(i_kind) ip2,im2
      integer(i_kind) kp2,km2
      
      ! Fill ghost
      bdy_type = 'slip_boundary'
      call bdy_condition(q_ext,stat%q,bdy_type)
      
      !$OMP PARALLEL DO PRIVATE(i,iVar,ip2,im2,q_weno,dir,kp2,km2)
      do k = kds,kde
        do i = ids,ide
          do iVar = 1,nVar
            ! x-dir
            ip2 = i + 2
            im2 = i - 2
            q_weno = q_ext(iVar,im2:ip2,k)
            
            dir = -1
            call WENO_limiter(qL(iVar,i,k),q_weno,dir)
            dir = 1
            call WENO_limiter(qR(iVar,i,k),q_weno,dir)
            
            ! z-dir
            kp2 = k + 2
            km2 = k - 2
            q_weno = q_ext(iVar,i,km2:kp2)
            
            dir = -1
            call WENO_limiter(qB(iVar,i,k),q_weno,dir)
            dir = 1
            call WENO_limiter(qT(iVar,i,k),q_weno,dir)
          enddo
          
          PL(i,k) = calc_pressure(sqrtG(i,k),qL(:,i,k))
          PR(i,k) = calc_pressure(sqrtG(i,k),qR(:,i,k))
          PB(i,k) = calc_pressure(sqrtG(i,k),qB(:,i,k))
          PT(i,k) = calc_pressure(sqrtG(i,k),qT(:,i,k))
          
          FL(:,i,k) = calc_F(sqrtG(i,k),qL(:,i,k),PL(i,k))
          FR(:,i,k) = calc_F(sqrtG(i,k),qR(:,i,k),PR(i,k))
          FB(:,i,k) = calc_F(sqrtG(i,k),qB(:,i,k),PB(i,k))
          FT(:,i,k) = calc_F(sqrtG(i,k),qT(:,i,k),PT(i,k))
          
          HL(:,i,k) = calc_H(sqrtG(i,k),G13(i,k),qL(:,i,k),PL(i,k))
          HR(:,i,k) = calc_H(sqrtG(i,k),G13(i,k),qR(:,i,k),PR(i,k))
          HB(:,i,k) = calc_H(sqrtG(i,k),G13(i,k),qB(:,i,k),PB(i,k))
          HT(:,i,k) = calc_H(sqrtG(i,k),G13(i,k),qT(:,i,k),PT(i,k))
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      !do k = kds,kde
      !  do i = ids,ide
      !    tend%q(:,i,k) = 
      !  enddo
      !enddo
      
    end subroutine spatial_operator
    
    ! 1D WENO slope limiter, according to Sun,2015
    ! "A Slope Constrained 4th OrderMulti-Moment Finite Volume Method with WENO Limiter"
    ! and Jiang and Shu, 1996
    subroutine WENO_limiter(Qrec,Q,dir)
      real   (r_kind)              , intent(out) :: Qrec
      real   (r_kind), dimension(5), intent(in ) :: Q
      integer(i_kind)              , intent(in ) :: dir
      
      integer(i_kind), parameter :: nStencil = 3
      real   (r_kind), parameter :: weno_coef(3)  = [0.1, 0.6, 0.3]
      real   (r_kind), parameter :: eps           = 1.E-2
      
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
      
      !! WENO-Z
      !tau40 = abs( beta(1) - beta(2) )
      !tau41 = abs( beta(2) - beta(3) )
      !tau5  = abs( beta(3) - beta(1) )
      !
      !! xi
      !if( tau40<=minval(beta) .and. tau41>minval(beta) )then
      !  omega = [1./3.,2./3.,0.]
      !elseif( tau40>minval(beta) .and. tau41<minval(beta) )then
      !  omega = [0.,2./3.,1./3.]
      !else
      !  alpha = weno_coef * ( 1. + tau5 / ( eps + beta ) )
      !  omega = alpha / sum(alpha)
      !endif
      !! WENO-Z
      
      ! origin WENO
      alpha = weno_coef / ( eps + beta )**2
      omega = alpha / sum(alpha)
      ! origin WENO
      
      Qrec = dot_product( stencil, omega )
      
    end subroutine WENO_limiter
    
    subroutine bdy_condition(q_ext,q,bdy_type)
      real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(out) :: q_ext
      real(r_kind), dimension(nVar,ids:ide,kds:kde), intent(in ) :: q
      character(30)                                , intent(in ) :: bdy_type
      
      integer(i_kind) dir
      
      integer(i_kind) i,k
      
      q_ext(:,ids:ide,kds:kde) = q
      
      if(trim(adjustl(bdy_type))=='slip_boundary')then
        dir = 1
        do i = 1,extPts
          call fill_ghost(q_ext(1,:,:),q(1,:,:),dir,pos)
          call fill_ghost(q_ext(2,:,:),q(2,:,:),dir,neg)
          call fill_ghost(q_ext(3,:,:),q(3,:,:),dir,pos)
          call fill_ghost(q_ext(4,:,:),q(4,:,:),dir,pos)
          call fill_ghost(q_ext(5,:,:),q(5,:,:),dir,pos)
        enddo
        
        dir = 2
        do k = 1,extPts
          call fill_ghost(q_ext(1,:,:),q(1,:,:),dir,pos)
          call fill_ghost(q_ext(2,:,:),q(2,:,:),dir,pos)
          call fill_ghost(q_ext(3,:,:),q(3,:,:),dir,neg)
          call fill_ghost(q_ext(4,:,:),q(4,:,:),dir,pos)
          call fill_ghost(q_ext(5,:,:),q(5,:,:),dir,pos)
        enddo
      elseif(trim(adjustl(bdy_type))=='nonreflecting')then
        
      endif
    end subroutine bdy_condition
    
    subroutine fill_ghost(q_ext,q,dir,sign)
      real   (r_kind), dimension(ics:ice,kcs:kce), intent(out) :: q_ext
      real   (r_kind), dimension(ids:ide,kds:kde), intent(in ) :: q
      integer(i_kind)                            , intent(in ) :: dir
      integer(i_kind)                            , intent(in ) :: sign
      
      integer(i_kind) i,k
      
      if(dir == 1)then
        ! x-dir
        do i = 1,extPts
          q_ext(ids-i,kds:kde) = sign * q(ids+i-1,kds:kde)
          q_ext(ide+i,kds:kde) = sign * q(ide-i+1,kds:kde)
        enddo
      elseif(dir == 2)then
        ! z-dir
        do k = 1,extPts
          q_ext(ids:ide,kds-k) = sign * q(ids:ide,kds+k-1)
          q_ext(ids:ide,kde+k) = sign * q(ids:ide,kde-k+1)
        enddo
      endif
      
    end subroutine fill_ghost
    
    function calc_pressure(sqrtG,q)
      real(r_kind) calc_pressure
      real(r_kind) sqrtG
      real(r_kind) q(5)
      
      real(r_kind) w1
      real(r_kind) w4
      real(r_kind) w5
      
      w1 = q(1)
      w4 = q(4)
      w5 = q(5)
      
      calc_pressure = p0 * ( ( 1. + eq*w5/w1 ) * Rd * w4 / sqrtG / p0 )**( ( Cpd + ( Cpv - Cpd ) * w5 / w1) / ( Cvd + ( Cvv - Cvd ) * w5 / w1 ) )
    end function calc_pressure
    
    function calc_F(sqrtG,q,P)
      real(r_kind),dimension(5) :: calc_F
      real(r_kind)              :: sqrtG
      real(r_kind),dimension(5) :: q(5)
      real(r_kind)              :: P   ! pressure
      
      real(r_kind) w1
      real(r_kind) w2
      real(r_kind) w3
      real(r_kind) w4
      real(r_kind) w5
      
      w1 = q(1)
      w2 = q(2)
      w3 = q(3)
      w4 = q(4)
      w5 = q(5)
      
      calc_F(1) = w2
      calc_F(2) = w2**2 / w1 + sqrtG * P
      calc_F(3) = w2 * w3 / w1
      calc_F(4) = w4 * w2 / w1
      calc_F(5) = w5 * w2 / w1
      
    end function calc_F
    
    function calc_H(sqrtG,G13,q,P)
      real(r_kind),dimension(5) :: calc_H
      real(r_kind)              :: sqrtG
      real(r_kind)              :: G13
      real(r_kind),dimension(5) :: q(5)
      real(r_kind)              :: P   ! pressure
      
      real(r_kind) w1
      real(r_kind) w2
      real(r_kind) w3
      real(r_kind) w4
      real(r_kind) w5
      real(r_kind) ww
      
      w1 = q(1)
      w2 = q(2)
      w3 = q(3)
      w4 = q(4)
      w5 = q(5)
      
      ww = w3 / sqrtG + G13 * w2
      
      calc_H(1) = ww
      calc_H(2) = ww * w2 / w1 + sqrtG * G13 * P
      calc_H(3) = ww * w3 / w1 + P
      calc_H(4) = ww * w4 / w1
      calc_H(5) = ww * w5 / w1
      
    end function calc_H
    
    function calc_eigenvalue_x(sqrtG,q)
      real(r_kind),dimension(5) :: calc_eigenvalue_x
      real(r_kind)              :: sqrtG
      real(r_kind),dimension(5) :: q(5)
      
      real(r_kind) w1
      real(r_kind) w2
      real(r_kind) w3
      real(r_kind) w4
      real(r_kind) w5
      
      w1 = q(1)
      w2 = q(2)
      w3 = q(3)
      w4 = q(4)
      w5 = q(5)
      
      real coef1,coef2,coef3
      
      coef1 = cvv**2 * w1**2 * w2 * w4 * w5**2                             &
            + cvv**2 * eq * w1 * w2 * w4 * w5**3                           &
            + cvd**2 * w1 * w2 * w4 * (w1 - w5)**2 * (w1 + eq*w5)          &
            + 2. * cvd * cvv * w1 * w2 * w4 * (w1 - w5) * w5 * (w1 + eq*w5)
      
      coef2 = sqrt( Rd * w1**2 * w4**3 * ( cpd * (w1 - w5) + cpv * w5 ) * ( cvd * (w1 - w5) + cvv * w5 )**3 ( w1 + eq * w5 )**3 * ( ( Rd * w4 * ( w1 + eq * w5 ) )&
                  / ( sqrtG * p0 * w1 ) )**( -1. + ( cpd * w1 - cpd * w5 + cpv * w5 ) / ( cvd * w1 - cvd * w5 + cvv * w5 ) ) )
      
      coef3 = w1**2 * w4 * ( cvd*( w1 - w5 ) + cvv * w5 )**2 ( w1 + eq * w5 )
      
      calc_eigenvalue_x(1) = w2 / w1
      calc_eigenvalue_x(2) = w2 / w1
      calc_eigenvalue_x(3) = w2 / w1
      calc_eigenvalue_x(4) = ( coef1 + coef2 ) / coef3
      calc_eigenvalue_x(5) = ( coef1 - coef2 ) / coef3
      
    end function calc_eigenvalue_x
    
END MODULE spatial_operators_mod

