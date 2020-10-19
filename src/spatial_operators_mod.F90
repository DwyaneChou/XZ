MODULE spatial_operators_mod
  use constants_mod
  use mesh_mod
  use parameters_mod
  use stat_mod
  use tend_mod
  use test_case_mod
  implicit none
  
  private
  
  public init_spatial_operator, spatial_operator
  
  integer(i_kind) :: pos
  integer(i_kind) :: neg
  
  integer(i_kind) :: ils
  integer(i_kind) :: ile
  integer(i_kind) :: kls
  integer(i_kind) :: kle
  integer(i_kind) :: irs
  integer(i_kind) :: ire
  integer(i_kind) :: krs
  integer(i_kind) :: kre
  integer(i_kind) :: ibs
  integer(i_kind) :: ibe
  integer(i_kind) :: kbs
  integer(i_kind) :: kbe
  integer(i_kind) :: its
  integer(i_kind) :: ite
  integer(i_kind) :: kts
  integer(i_kind) :: kte
  
  real(r_kind), dimension(:,:,:), allocatable :: q_ext ! Extended forecast variables
  
  real(r_kind), dimension(:,:,:), allocatable :: qL    ! Reconstructed q_(i-1/2,k)
  real(r_kind), dimension(:,:,:), allocatable :: qR    ! Reconstructed q_(i+1/2,k)
  real(r_kind), dimension(:,:,:), allocatable :: qB    ! Reconstructed q_(i,k-1/2)
  real(r_kind), dimension(:,:,:), allocatable :: qT    ! Reconstructed q_(i,k+1/2)
  
  real(r_kind), dimension(:,:,:), allocatable :: F
  real(r_kind), dimension(:,:,:), allocatable :: H
  real(r_kind), dimension(  :,:), allocatable :: P
  
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
  
  real(r_kind), dimension(:,:,:), allocatable :: Fe    ! F on edges of each cell
  real(r_kind), dimension(:,:,:), allocatable :: He    ! H on edges of each cell
  
  real(r_kind), dimension(:,:,:), allocatable :: src   ! source term
  
  real(r_kind), dimension(:,:  ), allocatable :: sqrtG_ext
  real(r_kind), dimension(:,:  ), allocatable :: G13_ext
  
    contains
    subroutine init_spatial_operator
      integer(i_kind) dir
      
      pos = 1
      neg = -1
      
      ils = ids
      ile = ide + 1
      kls = kds
      kle = kde
      
      irs = ids - 1
      ire = ide
      krs = kds
      kre = kde
      
      ibs = ids
      ibe = ide
      kbs = kds
      kbe = kde + 1
      
      its = ids
      ite = ide
      kts = kds - 1
      kte = kde
      
      allocate(q_ext(nVar,ics:ice,kcs:kce))
      
      allocate(qL   (nVar,ils:ile,kls:kle))
      allocate(qR   (nVar,irs:ire,krs:kre))
      allocate(qB   (nVar,ibs:ibe,kbs:kbe))
      allocate(qT   (nVar,its:ite,kts:kte))
      
      allocate(F    (nVar,ids:ide,kds:kde))
      allocate(H    (nVar,ids:ide,kds:kde))
      allocate(P    (     ids:ide,kds:kde))
      
      allocate(FL   (nVar,ils:ile,kls:kle))
      allocate(FR   (nVar,irs:ire,krs:kre))
      
      allocate(HB   (nVar,ibs:ibe,kbs:kbe))
      allocate(HT   (nVar,its:ite,kts:kte))
      
      allocate(PL   (     ils:ile,kls:kle))
      allocate(PR   (     irs:ire,krs:kre))
      allocate(PB   (     ibs:ibe,kbs:kbe))
      allocate(PT   (     its:ite,kts:kte))
      
      allocate(Fe   (nVar,ids:ide+1,kds:kde  ))
      allocate(He   (nVar,ids:ide  ,kds:kde+1))
      
      allocate(src  (nVar,ids:ide,kds:kde))
    end subroutine init_spatial_operator
    
    subroutine spatial_operator(stat,tend)
      type(stat_field), target, intent(in   ) :: stat
      type(tend_field), target, intent(inout) :: tend
      
      real(r_kind), dimension(5) :: q_weno
      
      real(r_kind) eigenvalue_x(5,2)
      real(r_kind) eigenvalue_z(5,2)
      real(r_kind) maxeigen_x
      real(r_kind) maxeigen_z
      
      integer(i_kind) dir
      
      integer(i_kind) i,k,iVar
      
      integer(i_kind) ip1,im1
      integer(i_kind) kp1,km1
      integer(i_kind) ip2,im2
      integer(i_kind) kp2,km2
      
      ! initialize source terms
      src = 0.
      
      ! Set no flux and nonreflecting boundary condition
      call bdy_condition(q_ext,stat%q(:,ids:ide,kds:kde),q_ref,src)
      
      !$OMP PARALLEL DO PRIVATE(i,iVar,ip1,im1,ip2,im2,q_weno,dir,kp1,km1,kp2,km2)
      do k = kds,kde
        kp1 = k + 1
        km1 = k - 1
        kp2 = k + 2
        km2 = k - 2
        do i = ids,ide
          ip1 = i + 1
          im1 = i - 1
          ip2 = i + 2
          im2 = i - 2
          do iVar = 1,nVar
            ! x-dir
            q_weno = q_ext(iVar,im2:ip2,k)
            if(im1<ids             )q_weno(1:2) = FillValue
            if(im2<ids.and.im1>=ids)q_weno(1  ) = FillValue
            if(ip1>ide             )q_weno(4:5) = FillValue
            if(ip2>ide.and.ip1<=ide)q_weno(5  ) = FillValue
            
            dir = -1
            call WENO_limiter(qL(iVar,i,k),q_weno,dir)
            dir = 1
            call WENO_limiter(qR(iVar,i,k),q_weno,dir)
            
            ! z-dir
            q_weno = q_ext(iVar,i,km2:kp2)
            if(km1<kds             )q_weno(1:2) = FillValue
            if(km2<kds.and.km1>=kds)q_weno(1  ) = FillValue
            if(kp1>kde             )q_weno(4:5) = FillValue
            if(kp2>kde.and.kp1<=kde)q_weno(5  ) = FillValue
            
            dir = -1
            call WENO_limiter(qB(iVar,i,k),q_weno,dir)
            dir = 1
            call WENO_limiter(qT(iVar,i,k),q_weno,dir)
          enddo
          if(i==ids)qL(2,i,k) = 0
          if(i==ide)qR(2,i,k) = 0
          if(k==kds)qB(3,i,k) = 0
          if(k==kde)qT(3,i,k) = 0
          
          PL(i,k) = calc_pressure(sqrtG(i,k),qL(:,i,k))
          PR(i,k) = calc_pressure(sqrtG(i,k),qR(:,i,k))
          PB(i,k) = calc_pressure(sqrtG(i,k),qB(:,i,k))
          PT(i,k) = calc_pressure(sqrtG(i,k),qT(:,i,k))
          
          FL(:,i,k) = calc_F(sqrtG(i,k),qL(:,i,k),PL(i,k))
          FR(:,i,k) = calc_F(sqrtG(i,k),qR(:,i,k),PR(i,k))
          
          HB(:,i,k) = calc_H(sqrtG(i,k),G13(i,k),qB(:,i,k),PB(i,k))
          HT(:,i,k) = calc_H(sqrtG(i,k),G13(i,k),qT(:,i,k),PT(i,k))
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! Fill boundary
      ! left
      qR(:,irs,:) = qL(:,ids,:)
      FR(:,irs,:) = FL(:,ids,:)
      
      ! right
      qL(:,ile,:) = qR(:,ide,:)
      FL(:,ile,:) = FR(:,ide,:)
      
      ! bottom
      qT(:,:,kts) = qB(:,:,kds)
      HT(:,:,kts) = HB(:,:,kds)
      
      ! top
      qB(:,:,kbe) = qT(:,:,kde)
      HB(:,:,kbe) = HT(:,:,kde)
      
      src(3,ids:ide,kds:kde) = src(3,ids:ide,kds:kde) - stat%q(1,ids:ide,kds:kde) * gravity
      
      ! calc x flux
      !$OMP PARALLEL DO PRIVATE(i,im1,eigenvalue_x,maxeigen_x,iVar)
      do k = kds,kde
        do i = ids+1,ide
          im1 = i - 1
          eigenvalue_x(:,1) = calc_eigenvalue_x(sqrtG(im1,k),q_ext(:,im1,k))
          eigenvalue_x(:,2) = calc_eigenvalue_x(sqrtG(i  ,k),q_ext(:,i  ,k))
          
          maxeigen_x = maxval(abs(eigenvalue_x))
          
          Fe(:,i,k) = 0.5 * ( FL(:,i,k) + FR(:,im1,k) - maxeigen_x * ( qL(:,i,k) - qR(:,im1,k) ) )
          !do iVar = 1,nVar
          !  if(abs(FL(iVar,i,k) + FR(iVar,im1,k))<=1.E-15)then
          !    Fe(iVar,i,k) = 0
          !  else
          !    Fe(iVar,i,k) = 0.5 * ( FL(iVar,i,k) + FR(iVar,im1,k) - maxeigen_x * ( qL(iVar,i,k) - qR(iVar,im1,k) ) )
          !  endif
          !enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      Fe(1,ids  ,kds:kde) = 0
      Fe(1,ide+1,kds:kde) = 0
      Fe(2,ids  ,kds:kde) = sqrtG(ids,kds:kde) * PL(ids,kds:kde)
      Fe(2,ide+1,kds:kde) = sqrtG(ide,kds:kde) * PR(ide,kds:kde)
      Fe(3,ids  ,kds:kde) = 0
      Fe(3,ide+1,kds:kde) = 0
      Fe(4,ids  ,kds:kde) = 0
      Fe(4,ide+1,kds:kde) = 0
      Fe(5,ids  ,kds:kde) = 0
      Fe(5,ide+1,kds:kde) = 0
      
      ! calc z flux
      !$OMP PARALLEL DO PRIVATE(k,km1,eigenvalue_z,maxeigen_z,iVar)
      do i = ids,ide
        do k = kds,kde+1
          km1 = k - 1
          eigenvalue_z(:,1) = calc_eigenvalue_z(sqrtG(i,km1),G13(i,km1),q_ext(:,i,km1))
          eigenvalue_z(:,2) = calc_eigenvalue_z(sqrtG(i  ,k),G13(i  ,k),q_ext(:,i  ,k))
          
          maxeigen_z = maxval(abs(eigenvalue_z))
          
          He(:,i,k) = 0.5 * ( HB(:,i,k) + HT(:,i,km1) - maxeigen_z * ( qB(:,i,k) - qT(:,i,km1) ) )
          !do iVar = 1,nVar
          !  if(abs(HB(iVar,i,k) + HT(iVar,i,km1))<=1.E-15)then
          !    He(iVar,i,k) = 0
          !  else
          !    He(iVar,i,k) = 0.5 * ( HB(iVar,i,k) + HT(iVar,i,km1) - maxeigen_z * ( qB(iVar,i,k) - qT(iVar,i,km1) ) )
          !  endif
          !enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      He(1,ids:ide,kds  ) = 0
      He(1,ids:ide,kde+1) = 0
      He(2,ids:ide,kds  ) = sqrtG(ids:ide,kds) * G13(ids:ide,kds) * PB(ids:ide,kds)
      He(2,ids:ide,kde+1) = sqrtG(ids:ide,kde) * G13(ids:ide,kde) * PT(ids:ide,kde)
      He(3,ids:ide,kds  ) = PB(ids:ide,kds)
      He(3,ids:ide,kde+1) = PT(ids:ide,kde)
      He(4,ids:ide,kds  ) = 0
      He(4,ids:ide,kde+1) = 0
      He(5,ids:ide,kds  ) = 0
      He(5,ids:ide,kde+1) = 0
      
      !$OMP PARALLEL DO PRIVATE(i,ip1,kp1)
      do k = kds,kde
        do i = ids,ide
          ip1 = i + 1
          kp1 = k + 1
          tend%q(:,i,k) = - ( ( Fe(:,ip1,k) - Fe(:,i,k) ) / dx + ( He(:,i,kp1) - He(:,i,k) ) / deta ) + src(:,i,k)
        enddo
      enddo
      !$OMP END PARALLEL DO
      
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
      
      !! WENO-Z
      !tau40 = abs( beta(1) - beta(2) )
      !tau41 = abs( beta(2) - beta(3) )
      !tau5  = abs( beta(3) - beta(1) )
      !
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
      
      !print*,Q
      !print*,alpha,beta
      
      if(any(Q==FillValue))then
        if(dir>0)then
          if(Q(1)==FillValue)alpha(1  ) = 0
          if(Q(2)==FillValue)alpha(1:2) = 0
          if(Q(4)==FillValue)alpha(2:3) = 0
          if(Q(5)==FillValue)alpha(3  ) = 0
        elseif(dir<0)then
          if(Q(5)==FillValue)alpha(1  ) = 0
          if(Q(4)==FillValue)alpha(1:2) = 0
          if(Q(2)==FillValue)alpha(2:3) = 0
          if(Q(1)==FillValue)alpha(3  ) = 0
        endif
        omega = alpha / sum(alpha)
      endif
      
      Qrec = dot_product( stencil, omega )
      
    end subroutine WENO_limiter
    
    subroutine bdy_condition(q_ext,q,q_ref,src)
      real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(out  ) :: q_ext
      real(r_kind), dimension(nVar,ids:ide,kds:kde), intent(in   ) :: q
      real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(in   ) :: q_ref
      real(r_kind), dimension(nVar,ids:ide,kds:kde), intent(inout) :: src
      
      real(r_kind), dimension(nVar,kds:kde) :: dqx
      real(r_kind), dimension(nVar,ids:ide) :: dqz
      
      integer(i_kind), parameter :: vs = 2
      integer(i_kind), parameter :: ve = 4
      integer(i_kind), parameter :: bdy_width = 5
      real   (r_kind), parameter :: exp_ceof  = 2
      
      integer(i_kind) dir
      integer(i_kind) i,k,iVar
      
      integer(i_kind) il,ir
      integer(i_kind) kt
      
      integer(i_kind) kls,kle ! pure lateral boundary layer indices
      integer(i_kind) its,ite ! pure top boundary layer indices
      
      real(r_kind) :: relax_coef(bdy_width)
      real(r_kind) :: max_exp
      
      q_ext(:,ids:ide,kds:kde) = q
      
      ! Nonreflecting condition
      kls = kds
      kle = kde-bdy_width
      its = ids+bdy_width
      ite = ide-bdy_width
      ! calculate relax coefficients
      !max_exp = exp( ( real( bdy_width - 1 ) / real(bdy_width) )**exp_ceof ) - 1.
      do i = 1,bdy_width
        relax_coef(i) = ( real( bdy_width - i + 1 ) / real( bdy_width ) )**4 / dt
        !relax_coef(i) = ( exp( ( real( bdy_width - i ) / real(bdy_width) )**exp_ceof ) - 1. ) / ( max_exp * dt )
      enddo
      
      !! pure zone
      !do i = 1,bdy_width
      !  il = i
      !  ir = ide-i+1
      !  kt = kde-i+1
      !  do iVar = 1,4
      !    src(iVar,ids:ide,kt     ) = - relax_coef(i) * ( q(iVar,ids:ide,kt) - q_ref(iVar,ids:ide,kt) )
      !  enddo
      !  !iVar = 4
      !  !src(iVar,ids:ide,kt     ) = - relax_coef(i) * ( q(iVar,ids:ide,kt) - q_ref(iVar,ids:ide,kt) / q_ref(1,ids:ide,kt) * q(1,ids:ide,kt) )
      !enddo
      
      ! pure zone
      do i = 1,bdy_width
        il = i
        ir = ide-i+1
        kt = kde-i+1
        do iVar = vs,ve
          src(iVar,il     ,kls:kle) = - relax_coef(i) * ( q(iVar,il,kls:kle) - q_ref(iVar,il,kls:kle) / q(1,il,kls:kle) * q(1,il,kls:kle) )
          src(iVar,ir     ,kls:kle) = - relax_coef(i) * ( q(iVar,ir,kls:kle) - q_ref(iVar,ir,kls:kle) / q(1,ir,kls:kle) * q(1,ir,kls:kle) )
          src(iVar,its:ite,kt     ) = - relax_coef(i) * ( q(iVar,its:ite,kt) - q_ref(iVar,its:ite,kt) / q(1,its:ite,kt) * q(1,its:ite,kt) )
        enddo
      enddo
      
      !overlap zone
      do k = 1,bdy_width
        do i = 1,bdy_width
          il = i
          ir = ide-i+1
          kt = kde-i+1
          do iVar = vs,ve
            src(iVar,il,kt) = - max( relax_coef(i), relax_coef(k) ) * ( q(iVar,il,kt) - q_ref(iVar,il,kt) / q(1,il,kt) * q(1,il,kt) )
            src(iVar,ir,kt) = - max( relax_coef(i), relax_coef(k) ) * ( q(iVar,ir,kt) - q_ref(iVar,ir,kt) / q(1,ir,kt) * q(1,ir,kt) )
          enddo
        enddo
      enddo
    end subroutine bdy_condition
    
    subroutine fill_ghost(q_ext,q,dir,sign)
      real   (r_kind), dimension(ics:ice,kcs:kce), intent(out) :: q_ext
      real   (r_kind), dimension(ids:ide,kds:kde), intent(in ) :: q
      integer(i_kind)                            , intent(in ) :: dir
      integer(i_kind)                            , intent(in ) :: sign
      
      integer(i_kind) i,k
      
      q_ext(ids:ide,kds:kde) = q
      
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
      
      real(r_kind) coef1,coef2,coef3
      
      w1 = q(1)
      w2 = q(2)
      w3 = q(3)
      w4 = q(4)
      w5 = q(5)
      
      coef1 = cvv**2 * w1**2 * w2 * w4 * w5**2                             &
            + cvv**2 * eq * w1 * w2 * w4 * w5**3                           &
            + cvd**2 * w1 * w2 * w4 * (w1 - w5)**2 * (w1 + eq*w5)          &
            + 2. * cvd * cvv * w1 * w2 * w4 * (w1 - w5) * w5 * (w1 + eq*w5)
      
      coef2 = sqrt( Rd * w1**2 * w4**3 * ( cpd * ( w1 - w5 ) + cpv * w5 )   &
                  * ( cvd * ( w1 - w5 ) + cvv * w5 )**3                     &
                  * ( w1 + eq * w5 )**3                                     &
                  * ( ( Rd * w4 * ( w1 + eq * w5 ) ) / ( sqrtG * p0 * w1 ) )&
                  **( -1. + ( cpd * w1 - cpd * w5 + cpv * w5 ) / ( cvd * w1 - cvd * w5 + cvv * w5 ) ) )
      
      coef3 = w1**2 * w4 * ( cvd * ( w1 - w5 ) + cvv * w5 )**2 * ( w1 + eq * w5 )
      
      calc_eigenvalue_x(1) = w2 / w1
      calc_eigenvalue_x(2) = calc_eigenvalue_x(1)
      calc_eigenvalue_x(3) = calc_eigenvalue_x(1)
      calc_eigenvalue_x(4) = ( coef1 - coef2 ) / coef3
      calc_eigenvalue_x(5) = ( coef1 + coef2 ) / coef3
      
    end function calc_eigenvalue_x
    
    function calc_eigenvalue_z(sqrtG,G13,q)
      real(r_kind),dimension(5) :: calc_eigenvalue_z
      real(r_kind)              :: sqrtG
      real(r_kind)              :: G13
      real(r_kind),dimension(5) :: q(5)
      
      real(r_kind) w1
      real(r_kind) w2
      real(r_kind) w3
      real(r_kind) w4
      real(r_kind) w5
      
      real(r_kind) drhoetadt
      
      real(r_kind) coef1,coef2,coef3
      
      w1 = q(1)
      w2 = q(2)
      w3 = q(3)
      w4 = q(4)
      w5 = q(5)
      
      coef1 = cvd**2 * w1 * ( sqrtG * G13 * w2 + w3 ) * w4 * ( w1 - w5 )**2 * ( w1 + eq * w5 )&
            + 2. * cvd * cvv * w1 * ( sqrtG * G13 * w2 + w3 ) * w4 * ( w1 - w5 ) * w5 * ( w1 + eq * w5 )&
            + cvv**2 * w1 * ( sqrtG * G13 * w2 + w3 ) * w4 * w5**2 * ( w1 + eq * w5 )
      
      coef2 = sqrt( ( 1. + ( sqrtG * G13 )**2 ) * Rd * w1**2 * w4**3        &
                  * ( cpd * ( w1 - w5 ) + cpv * w5 )                        &
                  * ( cvd * ( w1 - w5 ) + cvv * w5 )**3                     &
                  * ( w1 + eq * w5 )**3                                     &
                  * ( ( Rd * w4 * ( w1 + eq * w5 ) ) / ( sqrtG * p0 * w1 ) )&
                  **( -1. + ( cpd * w1 - cpd * w5 + cpv * w5 ) / ( cvd * w1 - cvd * w5 + cvv * w5 ) ) )
      
      coef3 = sqrtG * w1**2 * w4 * ( cvd * ( w1 - w5 ) + cvv * w5 )**2 * ( w1 + eq * w5 )
      
      drhoetadt = ( G13 * w2 + w3 / sqrtG ) / w1
      
      calc_eigenvalue_z(1) = drhoetadt
      calc_eigenvalue_z(2) = drhoetadt
      calc_eigenvalue_z(3) = drhoetadt
      calc_eigenvalue_z(4) = ( coef1 - coef2 ) / coef3
      calc_eigenvalue_z(5) = ( coef1 + coef2 ) / coef3
      
    end function calc_eigenvalue_z
    
END MODULE spatial_operators_mod

