MODULE spatial_operators_mod
  use constants_mod
  use mesh_mod
  use parameters_mod
  use stat_mod
  use tend_mod
  use reconstruction_mod, only: weno_limiter
  implicit none
  
  private
  
  public init_spatial_operator, &
         spatial_operator
  
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
  
  real(r_kind), dimension(  :,:), allocatable :: rho_p ! density perturbation
  
  real(r_kind), dimension(:,:  ), allocatable :: PL_ref! Reconstructed P_ref_(i-1/2,k)
  real(r_kind), dimension(:,:  ), allocatable :: PR_ref! Reconstructed P_ref_(i+1/2,k)
  real(r_kind), dimension(:,:  ), allocatable :: PB_ref! Reconstructed P_ref_(i,k-1/2)
  real(r_kind), dimension(:,:  ), allocatable :: PT_ref! Reconstructed P_ref_(i,k+1/2)
  
  real(r_kind), dimension(:,:), allocatable :: eig_x
  real(r_kind), dimension(:,:), allocatable :: eig_z
      
  real(r_kind), dimension(:,:), allocatable :: relax_coef ! Relax coefficient of Rayleigh damping
  
    contains
    subroutine init_spatial_operator
      integer(i_kind) dir
      integer(i_kind) i,k,iVar
      
      real(r_kind), dimension(5) :: q_weno
      
      integer(i_kind) ip1,im1
      integer(i_kind) kp1,km1
      integer(i_kind) ip2,im2
      integer(i_kind) kp2,km2
      
      allocate(q_ext(nVar,ics:ice,kcs:kce))
      
      allocate(qL   (nVar,ics:ice,kcs:kce))
      allocate(qR   (nVar,ics:ice,kcs:kce))
      allocate(qB   (nVar,ics:ice,kcs:kce))
      allocate(qT   (nVar,ics:ice,kcs:kce))
      
      allocate(F    (nVar,ids:ide,kds:kce))
      allocate(H    (nVar,ids:ide,kds:kce))
      allocate(P    (     ids:ide,kds:kce))
      
      allocate(FL   (nVar,ics:ice,kcs:kce))
      allocate(FR   (nVar,ics:ice,kcs:kce))
      
      allocate(HB   (nVar,ics:ice,kcs:kce))
      allocate(HT   (nVar,ics:ice,kcs:kce))
      
      allocate(PL   (     ics:ice,kcs:kce))
      allocate(PR   (     ics:ice,kcs:kce))
      allocate(PB   (     ics:ice,kcs:kce))
      allocate(PT   (     ics:ice,kcs:kce))
      
      allocate(Fe   (nVar,ids:ide+1,kds:kde  ))
      allocate(He   (nVar,ids:ide  ,kds:kde+1))
      
      allocate(src  (nVar,ids:ide,kds:kde))
      
      allocate(rho_p (    ids:ide,kds:kde))
      
      allocate(PL_ref(    ids:ide,kds:kde))
      allocate(PR_ref(    ids:ide,kds:kde))
      allocate(PB_ref(    ids:ide,kds:kde))
      allocate(PT_ref(    ids:ide,kds:kde))
      
      allocate(eig_x(     ics:ice,kcs:kce))
      allocate(eig_z(     ics:ice,kcs:kce))
      
      allocate(relax_coef(ics:ice,kcs:kce))
      
      ! Set reference pressure
      q_ext = ref%q
      call fill_ghost(q_ext,ref%q)
      !$OMP PARALLEL DO PRIVATE(kp1,km1,kp2,km2,i,ip1,im1,ip2,im2,iVar,q_weno,dir)
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
            
            dir = -1
            call WENO_limiter(qL(iVar,i,k),q_weno,dir)
            dir = 1
            call WENO_limiter(qR(iVar,i,k),q_weno,dir)
            
            ! z-dir
            q_weno = q_ext(iVar,i,km2:kp2)
            
            dir = -1
            call WENO_limiter(qB(iVar,i,k),q_weno,dir)
            dir = 1
            call WENO_limiter(qT(iVar,i,k),q_weno,dir)
          enddo
          
          PL_ref(i,k) = calc_pressure(sqrtGL(i,k),qL(:,i,k))
          PR_ref(i,k) = calc_pressure(sqrtGR(i,k),qR(:,i,k))
          PB_ref(i,k) = calc_pressure(sqrtGB(i,k),qB(:,i,k))
          PT_ref(i,k) = calc_pressure(sqrtGT(i,k),qT(:,i,k))
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! Calculate Rayleigh damping coef
      call Rayleigh_coef(relax_coef)
      
    end subroutine init_spatial_operator
    
    subroutine spatial_operator(stat,tend)
      type(stat_field), target, intent(inout) :: stat
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
      
      ! Attension stat is changed here!
      if(case_num==2)call Rayleigh_damping(stat%q,ref%q)
      
      ! copy stat
      q_ext = stat%q
      
      ! Fill ghost cells
      call fill_ghost(q_ext,ref%q)
      
      ! Reconstruct X
      !$OMP PARALLEL DO PRIVATE(i,iVar,ip1,im1,ip2,im2,q_weno,dir)
      do k = kds,kde
        !do i = ids-1,ide+1
        do i = ids,ide
          ip1 = i + 1
          im1 = i - 1
          ip2 = i + 2
          im2 = i - 2
          do iVar = 1,nVar
            ! x-dir
            q_weno = q_ext(iVar,im2:ip2,k)
            
            dir = -1
            call WENO_limiter(qL(iVar,i,k),q_weno,dir)
            dir = 1
            call WENO_limiter(qR(iVar,i,k),q_weno,dir)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! Reconstruct Z
      !$OMP PARALLEL DO PRIVATE(i,iVar,q_weno,dir,kp1,km1,kp2,km2)
      do k = kds,kde
        kp1 = k + 1
        km1 = k - 1
        kp2 = k + 2
        km2 = k - 2
        do i = ids,ide
          do iVar = 1,nVar
            ! z-dir
            q_weno = q_ext(iVar,i,km2:kp2)
            
            dir = -1
            call WENO_limiter(qB(iVar,i,k),q_weno,dir)
            dir = 1
            call WENO_limiter(qT(iVar,i,k),q_weno,dir)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! Boundary Condition
      ! Fill boundary
      if(case_num==1)then
        qL(2,ids,:) = 0
        qR(2,ide,:) = 0
        qB(3,:,kds) = 0
        qT(3,:,kde) = 0
      elseif(case_num==2)then
        qL(2,ids,kds:kde) = ref%q(2,ids,kds:kde)
        qR(2,ide,kds:kde) = ref%q(2,ide,kds:kde)
        !qL(2,ids,kds:kde) = ref%q(2,ids,kds:kde) / ref%q(1,ids,kds:kde) * q_ext(1,ids,kds:kde)
        !qR(2,ide,kds:kde) = ref%q(2,ide,kds:kde) / ref%q(1,ide,kds:kde) * q_ext(1,ids,kds:kde)
        qB(3,:,kds) = -sqrtGB(:,kds) * G13B(:,kds) * qB(2,:,kds)
        qT(3,:,kde) = -sqrtGT(:,kde) * G13T(:,kde) * qT(2,:,kde)
      endif
      
      ! left boundary
      qR(:,ids-1,:) = qL(:,ids,:)
      
      ! right boundary
      qL(:,ide+1,:) = qR(:,ide,:)
      
      ! bottom boundary
      qT(:,:,kds-1) = qB(:,:,kds)
      
      ! top boundary
      qB(:,:,kde+1) = qT(:,:,kde)
      
      ! Calculate functions X
      !$OMP PARALLEL DO PRIVATE(i)
      do k = kds,kde
        !do i = ids-1,ide+1
        do i = ids,ide
          PL(i,k) = calc_pressure(sqrtGL(i,k),qL(:,i,k))
          PR(i,k) = calc_pressure(sqrtGR(i,k),qR(:,i,k))
          
          FL(:,i,k) = calc_F(sqrtGL(i,k),qL(:,i,k),PL(i,k))
          FR(:,i,k) = calc_F(sqrtGR(i,k),qR(:,i,k),PR(i,k))
        enddo
      enddo
      !$OMP END PARALLEL DO
      FR(:,ids-1,:) = FL(:,ids,:)
      FL(:,ide+1,:) = FR(:,ide,:)
      
      ! Calculate functions Z
      !$OMP PARALLEL DO PRIVATE(i)
      do k = kds,kde
        do i = ids,ide
          PB(i,k) = calc_pressure(sqrtGB(i,k),qB(:,i,k))
          PT(i,k) = calc_pressure(sqrtGT(i,k),qT(:,i,k))
          
          HB(:,i,k) = calc_H(sqrtGB(i,k),G13B(i,k),qB(:,i,k),PB(i,k),PB_ref(i,k))
          HT(:,i,k) = calc_H(sqrtGT(i,k),G13T(i,k),qT(:,i,k),PT(i,k),PT_ref(i,k))
        enddo
      enddo
      !$OMP END PARALLEL DO
      HT(:,:,kds-1) = HB(:,:,kds)
      HB(:,:,kde+1) = HT(:,:,kde)
      
      ! initialize source terms
      src = 0.
      
      rho_p = ( stat%q(1,ids:ide,kds:kde) + stat%q(5,ids:ide,kds:kde) &
              - ref%q (1,ids:ide,kds:kde) - ref%q (5,ids:ide,kds:kde) ) / sqrtG(ids:ide,kds:kde)
      
      where(abs(rho_p)<=1.E-13)rho_p=0.
      
      src(3,ids:ide,kds:kde) = src(3,ids:ide,kds:kde) - sqrtG(ids:ide,kds:kde) * rho_p(ids:ide,kds:kde) * gravity
      
      ! Calculate eigenvalues x-dir
      !$OMP PARALLEL DO PRIVATE(i)
      do k = kds,kde
        do i = ids,ide
          eig_x(i,k) = calc_eigenvalue_x(sqrtG(i,k),q_ext(:,i,k))
        enddo
        i = ids - 1
        eig_x(i,k) = 0
        i = ide + 1
        eig_x(i,k) = 0
      enddo
      !$OMP END PARALLEL DO
      
      ! Calculate eigenvalues z-dir
      !$OMP PARALLEL DO PRIVATE(k)
      do i = ids,ide
        do k = kds,kde
          eig_z(i,k) = calc_eigenvalue_z(sqrtG(i,k),G13(i,k),q_ext(:,i,k))
        enddo
        k = kds - 1
        eig_z(i,k) = 0
        k = kde + 1
        eig_z(i,k) = 0
      enddo
      !$OMP END PARALLEL DO
      
      ! calc x flux
      !$OMP PARALLEL DO PRIVATE(i,im1,maxeigen_x,iVar)
      do k = kds,kde
        do i = ids,ide+1
          im1 = i - 1
          
          maxeigen_x = max(abs(eig_x(i,k)),abs(eig_x(im1,k)))
          
          !Fe(:,i,k) = 0.5 * ( FL(:,i,k) + FR(:,im1,k) - maxeigen_x * ( qL(:,i,k) - qR(:,im1,k) ) )
          
          do iVar = 1,nVar
            if(abs(FL(iVar,i,k) + FR(iVar,im1,k))<=1.E-15)then
              Fe(iVar,i,k) = 0
            else
              Fe(iVar,i,k) = 0.5 * ( FL(iVar,i,k) + FR(iVar,im1,k) - maxeigen_x * ( qL(iVar,i,k) - qR(iVar,im1,k) ) )
            endif
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! calc z flux
      !$OMP PARALLEL DO PRIVATE(k,km1,maxeigen_z,iVar)
      do i = ids,ide
        do k = kds,kde+1
          km1 = k - 1
          
          maxeigen_z = max(abs(eig_z(i,k)),abs(eig_z(i,km1)))
          
          !He(:,i,k) = 0.5 * ( HB(:,i,k) + HT(:,i,km1) - maxeigen_z * ( qB(:,i,k) - qT(:,i,km1) ) )
          
          do iVar = 1,nVar
            if(abs(HB(iVar,i,k) + HT(iVar,i,km1))<=1.E-15)then
              He(iVar,i,k) = 0
            else
              He(iVar,i,k) = 0.5 * ( HB(iVar,i,k) + HT(iVar,i,km1) - maxeigen_z * ( qB(iVar,i,k) - qT(iVar,i,km1) ) )
            endif
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      
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
    
    subroutine fill_ghost(q,q_ref)
      real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(inout) :: q
      real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(in   ) :: q_ref
      
      integer(i_kind) dir
      integer(i_kind) i,k,iVar
      
      if(case_num==1)then
        ! No-flux
        ! left
        q(:,ics:ids-1,:) = FillValue
        
        ! right
        q(:,ide+1:ice,:) = FillValue
        
        ! bottom
        q(:,:,kcs:kds-1) = FillValue
        
        ! top
        q(:,:,kde+1:kce) = FillValue
        
      elseif(case_num==2)then
        ! left
        q(:,ics:ids-1,:) = FillValue
        !q(:,ids:ids+2,:) = q_ref(:,ids:ids+2,:)
        
        ! right
        q(:,ide+1:ice,:) = FillValue
        !q(:,ide-2:ide,:) = q_ref(:,ide-2:ide,:)
        
        !! left
        !q(:,ics:ids-1,:) = q_ref(:,ics:ids-1,:)
        !
        !! right
        !q(:,ide+1:ice,:) = q_ref(:,ide+1:ice,:)
        
        ! bottom
        q(:,:,kcs:kds-1) = FillValue
        
        ! top
        q(:,:,kde+1:kce) = FillValue
        
      endif
      
    end subroutine fill_ghost
    
    subroutine Rayleigh_damping(q,q_ref)
      real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(inout) :: q
      real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(in   ) :: q_ref
      
      integer(i_kind), parameter :: vs = 1
      integer(i_kind), parameter :: ve = 5
      
      integer i,k,iVar
      
      do iVar = vs,ve
        q(iVar,ids:ide,kds:kde) = q(iVar,ids:ide,kds:kde) - relax_coef(ids:ide,kds:kde) * ( q(iVar,ids:ide,kds:kde) - q_ref(iVar,ids:ide,kds:kde) )
      enddo
    end subroutine Rayleigh_damping
    
    ! Rayleigh damping ( Wong and Stull, MWR, 2015 )
    subroutine Rayleigh_coef(mu)
      real(r_kind), dimension(ics:ice,kcs:kce), intent(out) :: mu
      
      real(r_kind), dimension(ics:ice,kcs:kce) :: muT
      real(r_kind), dimension(ics:ice,kcs:kce) :: muL
      real(r_kind), dimension(ics:ice,kcs:kce) :: muR
      
      real(r_kind), parameter :: topSpongeThickness   = 12000
      real(r_kind), parameter :: leftSpongeThickness  = 10000
      real(r_kind), parameter :: rightSpongeThickness = 10000
      
      real(r_kind), parameter :: mu_max_top = 0.5
      real(r_kind), parameter :: mu_max_lat = 0.5
      
      real(r_kind) zd, zt
      
      integer i,k
      
      muT = 0
      muL = 0
      MuR = 0
      
      ! Top
      zt = z_max
      zd = zt - topSpongeThickness
      where( z > zd )
        !muT = mu_max_top * sin( pi / 2. * ( z - zd ) / ( zt - zd ) )**2  !( Wong and Stull, MWR, 2015 )
        muT = mu_max_top * ( ( z - zd ) / ( zt - zd ) )**4 ! ( Li Xingliang, MWR, 2013 )
      elsewhere
        muT = 0.
      endwhere
      
      ! Left
      zt = -x_min
      zd = zt - leftSpongeThickness
      where( abs(x) > zd )
        !muT = mu_max_lat * sin( pi / 2. * ( abs(x) - zd ) / ( zt - zd ) )**2  !( Wong and Stull, MWR, 2015 )
        muT = mu_max_lat * ( ( abs(x) - zd ) / ( zt - zd ) )**4 ! ( Li Xingliang, MWR, 2013 )
      elsewhere
        muL = 0.
      endwhere
      
      ! Right
      zt = x_max
      zd = zt - rightSpongeThickness
      where( x > zd )
        !muT = mu_max_lat * sin( pi / 2. * ( abs(x) - zd ) / ( zt - zd ) )**2  !( Wong and Stull, MWR, 2015 )
        muT = mu_max_lat * ( ( abs(x) - zd ) / ( zt - zd ) )**4 ! ( Li Xingliang, MWR, 2013 )
      elsewhere
        muR = 0.
      endwhere
      
      do k = kds,kde
        do i = ids,ide
          mu(i,k) = max( muT(i,k), muL(i,k), muR(i,k) )
        enddo
      enddo
    
    end subroutine Rayleigh_coef
    
    function calc_pressure(sqrtG,q)
      real(r_kind) calc_pressure
      real(r_kind) sqrtG
      real(r_kind) q(5)
      
      real(r_kind) w1
      real(r_kind) w4
      real(r_kind) w5
      
      real(r_kind) rho
      real(r_kind) R
      real(r_kind) gamma
      real(r_kind) theta
      real(r_kind) cp
      real(r_kind) cv
      real(r_kind) sh
      real(r_kind) kappa
      
      w1 = q(1)
      w4 = q(4)
      w5 = q(5)
      
      rho      = ( w1 + w5 ) / sqrtG
      gamma    = w5 / w1
      sh       = gamma / ( 1. + gamma )
      R        = ( 1. + eq * sh ) * Rd
      theta    = w4 / w1
      cp       = cpd + ( cpv - cpd ) * sh
      cv       = cvd + ( cvv - cvd ) * sh
      kappa    = cp / cv
      
      calc_pressure = p0 * ( rho * theta * R / p0 )**kappa
      
    end function calc_pressure
    
    function calc_F(sqrtG,q,p)
      real(r_kind),dimension(5) :: calc_F
      real(r_kind)              :: sqrtG
      real(r_kind),dimension(5) :: q(5)
      real(r_kind)              :: p      ! pressure
      
      real(r_kind) w1
      real(r_kind) w2
      real(r_kind) w3
      real(r_kind) w4
      real(r_kind) w5
      
      real(r_kind) sqrtGrho
      real(r_kind) u
      
      w1 = q(1)
      w2 = q(2)
      w3 = q(3)
      w4 = q(4)
      w5 = q(5)
      
      sqrtGrho = w1 + w5
      u        = w2 / sqrtGrho
      
      calc_F(1) = w1 * u
      calc_F(2) = w2 * u + sqrtG * p
      calc_F(3) = w3 * u
      calc_F(4) = w4 * u
      calc_F(5) = w5 * u
      
    end function calc_F
    
    function calc_H(sqrtG,G13,q,p,p_ref)
      real(r_kind),dimension(5) :: calc_H
      real(r_kind)              :: sqrtG
      real(r_kind)              :: G13
      real(r_kind),dimension(5) :: q(5)
      real(r_kind)              :: p      ! pressure
      real(r_kind)              :: p_ref  ! reference pressure
      
      real(r_kind)              :: p_pert ! pressure  perturbation
      
      real(r_kind) w1
      real(r_kind) w2
      real(r_kind) w3
      real(r_kind) w4
      real(r_kind) w5
      real(r_kind) ww
      
      real(r_kind) sqrtGrho
      real(r_kind) u
      real(r_kind) w
      
      w1 = q(1)
      w2 = q(2)
      w3 = q(3)
      w4 = q(4)
      w5 = q(5)
      
      sqrtGrho = w1 + w5
      u        = w2 / sqrtGrho
      w        = w3 / sqrtGrho
      p_pert   = p - p_ref
      
      ww = w / sqrtG + G13 * u
      
      calc_H(1) = w1 * ww
      calc_H(2) = w2 * ww + sqrtG * G13 * p
      calc_H(3) = w3 * ww + p_pert
      calc_H(4) = w4 * ww
      calc_H(5) = w5 * ww
      
    end function calc_H
    
    function calc_eigenvalue_x(sqrtG,q)
      real(r_kind) :: calc_eigenvalue_x
      real(r_kind) :: sqrtG
      real(r_kind) :: q(5)
      real(r_kind) :: eig(5)
      
      real(r_kind) w1
      real(r_kind) w2
      real(r_kind) w3
      real(r_kind) w4
      real(r_kind) w5
      
      real(r_kind) coef1,coef2,coef3
      
      if(any(q==FillValue).or.sqrtG==FillValue)then
        eig = 0
      else
        w1 = q(1)
        w2 = q(2)
        w3 = q(3)
        w4 = q(4)
        w5 = q(5)
        
        coef1 = cvd**2*w1**3*w2*w4*(w1 + w5)*(w1 + w5 + eq*w5) + &
              2.*cvd*cvv*w1**2*w2*w4*w5*(w1 + w5)*(w1 + w5 + eq*w5) + &
              cvv**2*w1*w2*w4*w5**2*(w1 + w5)*(w1 + w5 + eq*w5)
        
        coef2 = sqrt( p0*sqrtG*w1**2*w4**2*(w1 + w5)**3*(cpd*w1 + cpv*w5)*&
              (cvd*w1 + cvv*w5)**3*(w1 + w5 + &
              eq*w5)**2*((Rd*w4*(w1 + w5 + eq*w5))/(p0*sqrtG*w1))**&
              ((cpd*w1 + cpv*w5)/(cvd*w1 + cvv*w5)) )
        
        coef3 = w1*w4*(w1 + w5)**2*(cvd*w1 + cvv*w5)**2*(w1 + w5 + eq*w5)
        
        eig(1) = w2 / ( w1 + w5 )
        eig(2) = eig(1)
        eig(3) = eig(1)
        eig(4) = ( coef1 - coef2 ) / coef3
        eig(5) = ( coef1 + coef2 ) / coef3
      endif
      
      calc_eigenvalue_x = maxval(eig)
      
    end function calc_eigenvalue_x
    
    function calc_eigenvalue_z(sqrtG,G13,q)
      real(r_kind) :: calc_eigenvalue_z
      real(r_kind) :: sqrtG
      real(r_kind) :: G13
      real(r_kind) :: q(5)
      real(r_kind) :: eig(5)
      
      real(r_kind) w1
      real(r_kind) w2
      real(r_kind) w3
      real(r_kind) w4
      real(r_kind) w5
      
      real(r_kind) sqrtGrho
      real(r_kind) u
      real(r_kind) w
      
      real(r_kind) drhoetadt
      
      real(r_kind) coef1,coef2,coef3
      
      if(any(q==FillValue))then
        eig = 0
      else
        w1 = q(1)
        w2 = q(2)
        w3 = q(3)
        w4 = q(4)
        w5 = q(5)
        
        !sqrtGrho = w1 + w5
        !u        = w2 / sqrtGrho
        !w        = w3 / sqrtGrho
        
        coef1 = cvd**2*w1**3*(G13*sqrtG*w2 + w3)*w4*(w1 + w5)*(w1 + w5 + eq*w5) + &
              2.*cvd*cvv*w1**2*(G13*sqrtG*w2 + w3)*w4*                            &
              w5*(w1 + w5)*(w1 + w5 + eq*w5) +                                    &
              cvv**2*w1*(G13*sqrtG*w2 + w3)*w4*w5**2*(w1 + w5)*(w1 + w5 + eq*w5)
        
        coef2 = sqrt( p0*sqrtG*(1 + G13**2*sqrtG**2)*w1**2*                    &
              w4**2*(w1 + w5)**3*(cpd*w1 + cpv*w5)*(cvd*w1 + cvv*w5)**3*       &
              (w1 + w5 + eq*w5)**2*((Rd*w4*(w1 + w5 + eq*w5))/(p0*sqrtG*w1))** &
              ((cpd*w1 + cpv*w5)/(cvd*w1 + cvv*w5)) )
        
        coef3 = sqrtG*w1*w4*(w1 + w5)**2*(cvd*w1 + cvv*w5)**2*(w1 + w5 + eq*w5)
        
        drhoetadt = (G13*sqrtG*w2 + w3)/(sqrtG*w1 + sqrtG*w5)
        
        eig(1) = drhoetadt
        eig(2) = drhoetadt
        eig(3) = drhoetadt
        eig(4) = ( coef1 - coef2 ) / coef3
        eig(5) = ( coef1 + coef2 ) / coef3
      endif
      
      calc_eigenvalue_z = maxval(eig)
      
    end function calc_eigenvalue_z
    
END MODULE spatial_operators_mod

