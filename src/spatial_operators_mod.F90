MODULE spatial_operators_mod
  use constants_mod
  use mesh_mod
  use parameters_mod
  use stat_mod
  use tend_mod
  use reconstruction_mod
  implicit none
  
  private
  
  public init_spatial_operator, &
         spatial_operator
  
  real(r_kind), dimension(:,:,:), allocatable :: qC ! Extended forecast variables
  
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
      
  real(r_kind), dimension(:,:,:), allocatable :: qL_ref
  real(r_kind), dimension(:,:,:), allocatable :: qR_ref
  real(r_kind), dimension(:,:,:), allocatable :: qB_ref
  real(r_kind), dimension(:,:,:), allocatable :: qT_ref
  
  real(r_kind), dimension(  :,:), allocatable :: rhoL_ref ! reference density
  real(r_kind), dimension(  :,:), allocatable :: rhoR_ref ! reference density
  real(r_kind), dimension(  :,:), allocatable :: rhoB_ref ! reference density
  real(r_kind), dimension(  :,:), allocatable :: rhoT_ref ! reference density
      
  real(r_kind), dimension(  :,:), allocatable :: cL_ref ! sound speed
  real(r_kind), dimension(  :,:), allocatable :: cR_ref ! sound speed
  real(r_kind), dimension(  :,:), allocatable :: cB_ref ! sound speed
  real(r_kind), dimension(  :,:), allocatable :: cT_ref ! sound speed
  
  real(r_kind), dimension(:,:  ), allocatable :: PL_ref! Reconstructed P_ref_(i-1/2,k)
  real(r_kind), dimension(:,:  ), allocatable :: PR_ref! Reconstructed P_ref_(i+1/2,k)
  real(r_kind), dimension(:,:  ), allocatable :: PB_ref! Reconstructed P_ref_(i,k-1/2)
  real(r_kind), dimension(:,:  ), allocatable :: PT_ref! Reconstructed P_ref_(i,k+1/2)
      
  real(r_kind), dimension(:,:), allocatable :: relax_coef ! Relax coefficient of Rayleigh damping
  
  real(r_kind), dimension(:,:,:), allocatable :: q_diff ! u wind, for viscosity terms only
  
    contains
    subroutine init_spatial_operator
      integer(i_kind) dir
      integer(i_kind) i,k,iVar,iEOC
      
      real(r_kind) :: nx,nz,nv(2),pm(2,2)
      
      real(r_kind), dimension(5) :: q_weno
      
      integer(i_kind) ip1,im1
      integer(i_kind) kp1,km1
      integer(i_kind) ip2,im2
      integer(i_kind) kp2,km2
      
      allocate(qC(nVar,ics:ice,kcs:kce))
      
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
      
      allocate(qL_ref(nVar,ics:ice,kcs:kce))
      allocate(qR_ref(nVar,ics:ice,kcs:kce))
      allocate(qB_ref(nVar,ics:ice,kcs:kce))
      allocate(qT_ref(nVar,ics:ice,kcs:kce))
      
      allocate(rhoL_ref(  ics:ice,kcs:kce)) ! reference density
      allocate(rhoR_ref(  ics:ice,kcs:kce)) ! reference density
      allocate(rhoB_ref(  ics:ice,kcs:kce)) ! reference density
      allocate(rhoT_ref(  ics:ice,kcs:kce)) ! reference density
      
      allocate(cL_ref(    ics:ice,kcs:kce)) ! sound speed
      allocate(cR_ref(    ics:ice,kcs:kce)) ! sound speed
      allocate(cB_ref(    ics:ice,kcs:kce)) ! sound speed
      allocate(cT_ref(    ics:ice,kcs:kce)) ! sound speed
      
      allocate(PL_ref(    ics:ice,kcs:kce))
      allocate(PR_ref(    ics:ice,kcs:kce))
      allocate(PB_ref(    ics:ice,kcs:kce))
      allocate(PT_ref(    ics:ice,kcs:kce))
      
      allocate(relax_coef(ics:ice,kcs:kce))
      
      allocate(q_diff(nVar,ics:ice,kcs:kce))
      
      src   = 0
      
      ! Set reference pressure
      qC = ref%q
      call fill_ghost(qC,ref%q)
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
            q_weno = qC(iVar,im2:ip2,k)
            
            dir = -1
            call WENO_limiter(qL_ref(iVar,i,k),q_weno,dir)
            dir = 1
            call WENO_limiter(qR_ref(iVar,i,k),q_weno,dir)
            
            ! z-dir
            q_weno = qC(iVar,i,km2:kp2)
            
            dir = -1
            call WENO_limiter(qB_ref(iVar,i,k),q_weno,dir)
            dir = 1
            call WENO_limiter(qT_ref(iVar,i,k),q_weno,dir)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! Fill outside boundary
      ! left boundary
      qR_ref(:,ids-1,kds:kde) = qL_ref(:,ids,kds:kde)
      sqrtGR(  ids-1,kds:kde) = sqrtGL(  ids,kds:kde)
      G13R  (  ids-1,kds:kde) = G13L  (  ids,kds:kde)
      
      ! right boundary
      qL_ref(:,ide+1,kds:kde) = qR_ref(:,ide,kds:kde)
      sqrtGL(  ide+1,kds:kde) = sqrtGR(  ide,kds:kde)
      G13L  (  ide+1,kds:kde) = G13R  (  ide,kds:kde)
      
      ! bottom boundary
      qT_ref(:,ids:ide,kds-1) = qB_ref(:,ids:ide,kds)
      sqrtGT(  ids:ide,kds-1) = sqrtGB(  ids:ide,kds)
      G13T  (  ids:ide,kds-1) = G13B  (  ids:ide,kds)
      
      ! top boundary
      qB_ref(:,ids:ide,kde+1) = qT_ref(:,ids:ide,kde)
      sqrtGB(  ids:ide,kde+1) = sqrtGT(  ids:ide,kde)
      G13B  (  ids:ide,kde+1) = G13T  (  ids:ide,kde)
    
      ! Calculate projection matrix for no-flux boundary
      k    = kds
      iEOC = 1
      do i = ids,ide
        nx = G13B(i,k) * sqrtGB(i,k)
        nz = 1
        
        nv(1) = nx
        nv(2) = nz
        nv    = nv / sqrt( dot_product( nv, nv ) ) ! calc unit norm vector
        nx    = nv(1)
        nz    = nv(2)
        
        pm(1,1) = 1. - nx**2
        pm(1,2) = -nx * nz
        pm(2,1) = pm(1,2)
        pm(2,2) = 1. - nz**2
        
        nvec(:  ,iEOC,i,k) = nv
        pmtx(:,:,iEOC,i,k) = pm
      enddo
      
      k    = kde
      iEOC = 3
      do i = ids,ide
        nx = G13T(i,k) * sqrtGT(i,k)
        nz = 1
        
        nv(1) = nx
        nv(2) = nz
        nv    = nv / sqrt( dot_product( nv, nv ) ) ! calc unit norm vector
        nx    = nv(1)
        nz    = nv(2)
        
        pm(1,1) = 1. - nx**2
        pm(1,2) = -nx * nz
        pm(2,1) = pm(1,2)
        pm(2,2) = 1. - nz**2
        
        nvec(:  ,iEOC,i,k) = nv
        pmtx(:,:,iEOC,i,k) = pm
      enddo
          
      do k = kds-1,kde+1
        do i = ids-1,ide+1
          rhoL_ref(i,k) = ( qL_ref(1,i,k) + qL_ref(5,i,k) ) / sqrtGL(i,k)
          rhoR_ref(i,k) = ( qR_ref(1,i,k) + qR_ref(5,i,k) ) / sqrtGR(i,k)
          rhoB_ref(i,k) = ( qB_ref(1,i,k) + qB_ref(5,i,k) ) / sqrtGB(i,k)
          rhoT_ref(i,k) = ( qT_ref(1,i,k) + qT_ref(5,i,k) ) / sqrtGT(i,k)
          
          cL_ref(i,k) = calc_sound_speed_x(sqrtGL(i,k)          ,qL_ref(:,i,k))
          cR_ref(i,k) = calc_sound_speed_x(sqrtGR(i,k)          ,qR_ref(:,i,k))
          cB_ref(i,k) = calc_sound_speed_z(sqrtGB(i,k),G13B(i,k),qB_ref(:,i,k)) * sqrtGB(i,k) ! unit: m/s
          cT_ref(i,k) = calc_sound_speed_z(sqrtGT(i,k),G13T(i,k),qT_ref(:,i,k)) * sqrtGT(i,k) ! unit: m/s
          
          PL_ref(i,k) = calc_pressure(sqrtGL(i,k),qL_ref(:,i,k))
          PR_ref(i,k) = calc_pressure(sqrtGR(i,k),qR_ref(:,i,k))
          PB_ref(i,k) = calc_pressure(sqrtGB(i,k),qB_ref(:,i,k))
          PT_ref(i,k) = calc_pressure(sqrtGT(i,k),qT_ref(:,i,k))
        enddo
      enddo
      
      ! Calculate Rayleigh damping coef
      call Rayleigh_coef(relax_coef)
      
      ! Fill out qC
      qC = FillValue
      
    end subroutine init_spatial_operator
    
    subroutine spatial_operator(stat,tend)
      type(stat_field), target, intent(inout) :: stat
      type(tend_field), target, intent(inout) :: tend
      
      real(r_kind), dimension(5) :: q_weno
      
      real(r_kind) maxeigen_x
      real(r_kind) maxeigen_z
      
      integer(i_kind) dir
      
      integer(i_kind) i,k,iVar,iEOC
      real   (r_kind) :: wind_vector(2)
      
      integer(i_kind) ip1,im1
      integer(i_kind) kp1,km1
      integer(i_kind) ip2,im2
      integer(i_kind) kp2,km2
      
      ! Attension stat is changed here!
      if(case_num==2)call Rayleigh_damping(stat%q,ref%q)
      
      ! copy stat
      qC(:,ids:ide,kds:kde) = stat%q(:,ids:ide,kds:kde)
      
      !! Fill ghost cells
      !call fill_ghost(qC,ref%q)
      
      ! Reconstruct X
      !$OMP PARALLEL
      !$OMP DO PRIVATE(i,iVar,q_weno,dir) COLLAPSE(3)
      do k = kds,kde
        do i = ids,ide
          do iVar = 1,nVar
            ! x-dir
            q_weno = qC(iVar,i-2:i+2,k)
            
            dir = -1
            !call WENO_limiter(qL(iVar,i,k),q_weno,dir)
            call WENO5(qL(iVar,i,k),q_weno,dir)
            dir = 1
            !call WENO_limiter(qR(iVar,i,k),q_weno,dir)
            call WENO5(qR(iVar,i,k),q_weno,dir)
          enddo
        enddo
      enddo
      !$OMP END DO NOWAIT
      
      ! Reconstruct Z
      !$OMP DO PRIVATE(i,iVar,q_weno,dir) COLLAPSE(3)
      do k = kds,kde
        do i = ids,ide
          do iVar = 1,nVar
            ! z-dir
            q_weno = qC(iVar,i,k-2:k+2)
            
            dir = -1
            !call WENO_limiter(qB(iVar,i,k),q_weno,dir)
            call WENO5(qB(iVar,i,k),q_weno,dir)
            dir = 1
            !call WENO_limiter(qT(iVar,i,k),q_weno,dir)
            call WENO5(qT(iVar,i,k),q_weno,dir)
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      
      ! Boundary Condition
      ! Correct wind on bottom and top boundaries
      k    = kds
      iEOC = 1
      !$OMP PARALLEL DO PRIVATE(wind_vector)
      do i = ids,ide
        wind_vector(1) = qB(2,i,k)
        wind_vector(2) = qB(3,i,k)
        wind_vector    = matmul(pmtx(:,:,iEOC,i,k),wind_vector)
        qB(2,i,k) = wind_vector(1)
        qB(3,i,k) = wind_vector(2)
      enddo
      !$OMP END PARALLEL DO
      
      k    = kde
      iEOC = 3
      !$OMP PARALLEL DO PRIVATE(wind_vector)
      do i = ids,ide
        wind_vector(1) = qT(2,i,k)
        wind_vector(2) = qT(3,i,k)
        wind_vector    = matmul(pmtx(:,:,iEOC,i,k),wind_vector)
        qT(2,i,k) = wind_vector(1)
        qT(3,i,k) = wind_vector(2)
      enddo
      !$OMP END PARALLEL DO
        
      ! Fill lateral boundary
      if(case_num==1.or.case_num==3)then
        qL(2,ids,kds:kde) = 0
        qR(2,ide,kds:kde) = 0
      elseif(case_num==2)then
        !$OMP PARALLEL DO
        do k = kds,kde
          qL(2,ids,k) = ref%q(2,ids,k)
          qR(2,ide,k) = ref%q(2,ide,k)
        enddo
        !$OMP END PARALLEL DO
      endif
      
      ! Fill outside boundary
      ! left boundary
      qR(:,ids-1,kds:kde) = qL(:,ids,kds:kde)
      
      ! right boundary
      qL(:,ide+1,kds:kde) = qR(:,ide,kds:kde)
      
      ! bottom boundary
      qT(:,ids:ide,kds-1) = qB(:,ids:ide,kds)
      
      ! top boundary
      qB(:,ids:ide,kde+1) = qT(:,ids:ide,kde)
      
      ! Calculate functions X
      !$OMP PARALLEL DO PRIVATE(i)
      do k = kds,kde
        do i = ids,ide
          PL(i,k) = calc_pressure(sqrtGL(i,k),qL(:,i,k))
          PR(i,k) = calc_pressure(sqrtGR(i,k),qR(:,i,k))
        enddo
      enddo
      !$OMP END PARALLEL DO
      do k = kds,kde
        i = ide + 1
        PL(i,k) = calc_pressure(sqrtGL(i,k),qL(:,i,k))
        i = ids - 1
        PR(i,k) = calc_pressure(sqrtGR(i,k),qR(:,i,k))
      enddo
      
      ! Calculate functions Z
      !$OMP PARALLEL DO PRIVATE(i)
      do k = kds,kde
        do i = ids,ide
          PB(i,k) = calc_pressure(sqrtGB(i,k),qB(:,i,k))
          PT(i,k) = calc_pressure(sqrtGT(i,k),qT(:,i,k))
        enddo
      enddo
      !$OMP END PARALLEL DO
      do i = ids,ide
        k = kde + 1
        PB(i,k) = calc_pressure(sqrtGB(i,k),qB(:,i,k))
        k = kds - 1
        PT(i,k) = calc_pressure(sqrtGT(i,k),qT(:,i,k))
      enddo
      
      !$OMP PARALLEL
      !$OMP DO PRIVATE(i,im1) COLLAPSE(2)
      do k = kds,kde
        do i = ids,ide+1
          im1 = i - 1
          Fe(:,i,k) = calc_F(sqrtGR(im1,k),sqrtGL(i,k),qR(:,im1,k),qL(:,i,k),pR(im1,k),pL(i,k))
        enddo
      enddo
      !$OMP END DO NOWAIT
      
      !$OMP DO PRIVATE(k,km1) COLLAPSE(2)
      do i = ids,ide
        do k = kds,kde+1
          km1 = k - 1
          He(:,i,k) = calc_H(sqrtGT(i,km1),sqrtGB(i,k),G13T(i,km1),G13B(i,k),qT(:,i,km1),qB(:,i,k),pT(i,km1),pB(i,k),&
                             rhoT_ref(i,km1),rhoB_ref(i,k),cT_ref(i,km1),cB_ref(i,k),pT_ref(i,km1),pB_ref(i,k))
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      
      ! Calculate source term
      !rho_p = ( qC   (1,ids:ide,kds:kde) + qC   (5,ids:ide,kds:kde) &
      !        - ref%q(1,ids:ide,kds:kde) - ref%q(5,ids:ide,kds:kde) ) / sqrtG(ids:ide,kds:kde)
      rho_p = ( qC   (1,ids:ide,kds:kde) + qC   (5,ids:ide,kds:kde) ) / sqrtG(ids:ide,kds:kde)
      
      where(abs(rho_p)<=1.E-13)rho_p=0.
      
      src = 0
      src(3,ids:ide,kds:kde) = src(3,ids:ide,kds:kde) - sqrtG(ids:ide,kds:kde) * rho_p(ids:ide,kds:kde) * gravity
      
      ! Viscosity terms for Density Current case only
      if(case_num==3)then
        do iVar = 2,4
          q_diff(iVar,ids:ide,kds:kde) = qC(iVar,ids:ide,kds:kde) / qC(1,ids:ide,kds:kde)
          
          q_diff(iVar,ids:ide,kds-1) = q_diff(iVar,ids:ide,kds)
          q_diff(iVar,ids:ide,kde+1) = q_diff(iVar,ids:ide,kde)
          q_diff(iVar,ids-1,kds:kde) = q_diff(iVar,ids,kds:kde)
          q_diff(iVar,ide+1,kds:kde) = q_diff(iVar,ide,kds:kde)
        enddo
        
        !$OMP PARALLEL DO PRIVATE(i,iVar)
        do k = kds,kde
          do i = ids,ide
            do iVar = 2,4
              src(iVar,i,k) = src(iVar,i,k) + viscosity_coef * qC(1,i,k) * ( ( q_diff(iVar,i+1,k) - 2. * q_diff(iVar,i,k) + q_diff(iVar,i-1,k) ) / dx  **2 &
                                                                           + ( q_diff(iVar,i,k+1) - 2. * q_diff(iVar,i,k) + q_diff(iVar,i,k-1) ) / deta**2 / sqrtG(i,k)**2 )
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
      endif
      
      ! Calculate tend
      !$OMP PARALLEL DO PRIVATE(i) COLLAPSE(2)
      do k = kds,kde
        do i = ids,ide
          tend%q(:,i,k) = - ( ( Fe(:,i+1,k) - Fe(:,i,k) ) / dx + ( He(:,i,k+1) - He(:,i,k) ) / deta ) + src(:,i,k)
        enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine spatial_operator
    
    subroutine fill_ghost(q,q_ref)
      real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(inout) :: q
      real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(in   ) :: q_ref
      
      integer(i_kind) dir
      integer(i_kind) i,k,iVar
      
      ! No-flux
      ! left
      q(:,ics:ids-1,:) = FillValue
      
      ! right
      q(:,ide+1:ice,:) = FillValue
      
      ! bottom
      q(:,:,kcs:kds-1) = FillValue
      
      ! top
      q(:,:,kde+1:kce) = FillValue
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
      
      real(r_kind), parameter :: topSpongeThickness   = 12000  ! 12000 for "best" result
      real(r_kind), parameter :: leftSpongeThickness  = 10000  ! 10000 for "best" result
      real(r_kind), parameter :: rightSpongeThickness = 10000  ! 10000 for "best" result
      
      real(r_kind), parameter :: mu_max_top = 0.1
      real(r_kind), parameter :: mu_max_lat = 0.1
      
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
      real(r_kind) q(nVar)
      
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
    
    function calc_F(sqrtGL,sqrtGR,qL,qR,pL,pR)
      real(r_kind) :: calc_F(nVar)
      real(r_kind) :: sqrtGL
      real(r_kind) :: sqrtGR
      real(r_kind) :: qL(nVar)
      real(r_kind) :: qR(nVar)
      real(r_kind) :: pL      ! pressure
      real(r_kind) :: pR      ! pressure
      
      real(r_kind) rhoL
      real(r_kind) rhoR
      
      real(r_kind) uL
      real(r_kind) uR
      
      real(r_kind) cL
      real(r_kind) cR
      
      real(r_kind) m
      real(r_kind) p ! sqrtG * p
      
      rhoL = ( qL(1) + qL(5) ) / sqrtGL
      rhoR = ( qR(1) + qR(5) ) / sqrtGR
      uL   = qL(2) / ( qL(1) + qL(5) )
      uR   = qR(2) / ( qR(1) + qR(5) )
      cL   = calc_sound_speed_x(sqrtGL,qL)
      cR   = calc_sound_speed_x(sqrtGR,qR)
      
      call AUSM_up_x(m,p,sqrtGL,sqrtGR,rhoL,rhoR,uL,uR,pL,pR,cL,cR)
      
      calc_F = 0.5 * m * ( qL + qR - sign(1.,m) * ( qR - qL ) )
      calc_F(2) = calc_F(2) + p
      
    end function calc_F
    
    function calc_H(sqrtGL,sqrtGR,G13L,G13R,qL,qR,pL,pR,rhoL_ref,rhoR_ref,cL_ref,cR_ref,pL_ref,pR_ref)
      real(r_kind) :: calc_H(nVar)
      real(r_kind) :: sqrtGL
      real(r_kind) :: sqrtGR
      real(r_kind) :: G13L
      real(r_kind) :: G13R
      real(r_kind) :: qL(nVar)
      real(r_kind) :: qR(nVar)
      real(r_kind) :: pL        ! pressure
      real(r_kind) :: pR        ! pressure
      real(r_kind) :: rhoL_ref  ! reference density
      real(r_kind) :: rhoR_ref  ! reference density
      real(r_kind) :: cL_ref    ! reference sound speed
      real(r_kind) :: cR_ref    ! reference sound speed
      real(r_kind) :: pL_ref    ! reference pressure
      real(r_kind) :: pR_ref    ! reference pressure

      real(r_kind) :: pL_pert ! pressure  perturbation
      real(r_kind) :: pR_pert ! pressure  perturbation
      
      
      real(r_kind) rhoL
      real(r_kind) rhoR
      
      real(r_kind) uL
      real(r_kind) uR
      
      real(r_kind) cL
      real(r_kind) cR
      
      real(r_kind) m
      real(r_kind) p      ! sqrtG * G13 * p
      real(r_kind) p_pert ! p - p_ref
      
      rhoL = ( qL(1) + qL(5) ) / sqrtGL
      rhoR = ( qR(1) + qR(5) ) / sqrtGR
      uL   = ( qL(3) + sqrtGL * G13L * qL(2) ) / ( qL(1) + qL(5) ) ! a sqrtG has been multipled here
      uR   = ( qR(3) + sqrtGR * G13R * qR(2) ) / ( qR(1) + qR(5) ) ! a sqrtG has been multipled here
      cL   = calc_sound_speed_z(sqrtGL,G13L,qL) * sqrtGL
      cR   = calc_sound_speed_z(sqrtGR,G13R,qR) * sqrtGR
      
      call AUSM_up_z(m, p, p_pert,                                                  &
                     sqrtGL, sqrtGR, G13L, G13R, rhoL, rhoR, uL, uR, cL, cR, pL, pR,&
                     rhoL_ref, rhoR_ref, cL_ref, cR_ref, pL_ref, pR_ref)
      
      calc_H = 0.5 * m * ( qL + qR - sign(1.,m) * ( qR - qL ) )
      calc_H(2) = calc_H(2) + p
      calc_H(3) = calc_H(3) + p_pert
      
    end function calc_H
    
    function calc_sound_speed_x(sqrtG,q)
      real(r_kind) :: calc_sound_speed_x
      real(r_kind) :: sqrtG
      real(r_kind) :: q(nVar)
      
      real(r_kind) w1
      real(r_kind) w2
      real(r_kind) w3
      real(r_kind) w4
      real(r_kind) w5
      
      real(r_kind) coef1,coef2
      
      if(any(q==FillValue).or.sqrtG==FillValue)then
        calc_sound_speed_x = 0
      else
        w1 = q(1)
        w2 = q(2)
        w3 = q(3)
        w4 = q(4)
        w5 = q(5)
        
        coef1 = sqrt( p0*sqrtG*w1**2*w4**2*(w1 + w5)**3*(cpd*w1 + cpv*w5)*&
              (cvd*w1 + cvv*w5)**3*(w1 + w5 + &
              eq*w5)**2*((Rd*w4*(w1 + w5 + eq*w5))/(p0*sqrtG*w1))**&
              ((cpd*w1 + cpv*w5)/(cvd*w1 + cvv*w5)) )
        
        coef2 = w1*w4*(w1 + w5)**2*(cvd*w1 + cvv*w5)**2*(w1 + w5 + eq*w5)
        
        calc_sound_speed_x = coef1 / coef2
      endif
      
    end function calc_sound_speed_x
    
    function calc_sound_speed_z(sqrtG,G13,q)
      real(r_kind) :: calc_sound_speed_z
      real(r_kind) :: sqrtG
      real(r_kind) :: G13
      real(r_kind) :: q(nVar)
      
      real(r_kind) w1
      real(r_kind) w2
      real(r_kind) w3
      real(r_kind) w4
      real(r_kind) w5
      
      real(r_kind) coef1,coef2
      
      if(any(q==FillValue))then
        calc_sound_speed_z = 0
      else
        w1 = q(1)
        w2 = q(2)
        w3 = q(3)
        w4 = q(4)
        w5 = q(5)
        
        coef1 = sqrt( p0*sqrtG*(1 + G13**2*sqrtG**2)*w1**2*                    &
              w4**2*(w1 + w5)**3*(cpd*w1 + cpv*w5)*(cvd*w1 + cvv*w5)**3*       &
              (w1 + w5 + eq*w5)**2*((Rd*w4*(w1 + w5 + eq*w5))/(p0*sqrtG*w1))** &
              ((cpd*w1 + cpv*w5)/(cvd*w1 + cvv*w5)) )
        
        coef2 = sqrtG*w1*w4*(w1 + w5)**2*(cvd*w1 + cvv*w5)**2*(w1 + w5 + eq*w5)
        
        calc_sound_speed_z = coef1 / coef2
      endif
      
    end function calc_sound_speed_z
    
    subroutine AUSM_up_x(m,p,sqrtGL,sqrtGR,rhoL,rhoR,uL,uR,pL,pR,cL,cR)
      real(r_kind),intent(out) :: m
      real(r_kind),intent(out) :: p
      real(r_kind),intent(in ) :: sqrtGL
      real(r_kind),intent(in ) :: sqrtGR
      real(r_kind),intent(in ) :: rhoL
      real(r_kind),intent(in ) :: rhoR
      real(r_kind),intent(in ) :: uL
      real(r_kind),intent(in ) :: uR
      real(r_kind),intent(in ) :: pL
      real(r_kind),intent(in ) :: pR
      real(r_kind),intent(in ) :: cL ! Sound speed
      real(r_kind),intent(in ) :: cR ! Sound speed
      
      real(r_kind),parameter :: Ku    = 0.75
      real(r_kind),parameter :: Kp    = 0.25
      real(r_kind),parameter :: sigma = 1.
      real(r_kind),parameter :: sp    = 1.
      real(r_kind),parameter :: sn    = -1.
      
      real(r_kind) :: rho
      real(r_kind) :: a
      real(r_kind) :: ML
      real(r_kind) :: MR
      real(r_kind) :: Mbar2
      real(r_kind) :: Mh
      
      rho = 0.5 * ( rhoL + rhoR )
      a   = 0.5 * ( cL + cR )
      
      ML = uL / a
      MR = uR / a
      
      Mbar2 = ( uL**2 + uR**2 ) / ( 2. * a**2 )
      
      Mh = M4( ML, sp ) + M4( MR, sn ) - Kp * max( 1. - sigma * Mbar2, 0. ) * ( PR - PL ) / ( rho * a**2 )
      m  = a * Mh
      
      p = P5(ML,sp) * sqrtGL * PL + P5(MR,sn) * sqrtGR * PR - Ku * P5(ML,sp) * P5(MR,sn) * ( sqrtGL * rhoL + sqrtGR * rhoR ) * a * ( uR - uL )
      
    end subroutine AUSM_up_x
    
    subroutine AUSM_up_z(m, p, p_pert,                                                  & ! output
                         sqrtGL, sqrtGR, G13L, G13R, rhoL, rhoR, uL, uR, cL, cR, pL, pR,& ! input state
                         rhoL_ref, rhoR_ref, cL_ref, cR_ref, pL_ref, pR_ref)              ! input reference state
      real(r_kind),intent(out) :: m
      real(r_kind),intent(out) :: p
      real(r_kind),intent(out) :: p_pert
      real(r_kind),intent(in ) :: sqrtGL
      real(r_kind),intent(in ) :: sqrtGR
      real(r_kind),intent(in ) :: G13L
      real(r_kind),intent(in ) :: G13R
      real(r_kind),intent(in ) :: rhoL
      real(r_kind),intent(in ) :: rhoR
      real(r_kind),intent(in ) :: uL
      real(r_kind),intent(in ) :: uR
      real(r_kind),intent(in ) :: cL       ! Sound speed
      real(r_kind),intent(in ) :: cR       ! Sound speed
      real(r_kind),intent(in ) :: pL
      real(r_kind),intent(in ) :: pR
      real(r_kind),intent(in ) :: rhoL_ref ! reference state
      real(r_kind),intent(in ) :: rhoR_ref ! reference state
      real(r_kind),intent(in ) :: cL_ref   ! reference state
      real(r_kind),intent(in ) :: cR_ref   ! reference state
      real(r_kind),intent(in ) :: pL_ref   ! reference state
      real(r_kind),intent(in ) :: pR_ref   ! reference state
      
      real(r_kind),parameter :: Ku    = 0.75
      real(r_kind),parameter :: Kp    = 0.25
      real(r_kind),parameter :: sigma = 1.
      real(r_kind),parameter :: sp    = 1.
      real(r_kind),parameter :: sn    = -1.
      
      real(r_kind) :: rho
      real(r_kind) :: a
      real(r_kind) :: ML
      real(r_kind) :: MR
      real(r_kind) :: Mbar2
      real(r_kind) :: Mh
      
      real(r_kind) :: rho_ref
      real(r_kind) :: a_ref
      
      real(r_kind) :: pL_pert
      real(r_kind) :: pR_pert
      
      real(r_kind) :: p_diff
      
      real(r_kind) :: coefL
      real(r_kind) :: coefR
      
      rho = 0.5 * ( rhoL + rhoR )
      a   = 0.5 * ( cL + cR )
      
      rho_ref = 0.5 * ( rhoL_ref + rhoR_ref )
      a_ref   = 0.5 * ( cL_ref + cR_ref )
      
      pL_pert = pL - pL_ref
      pR_pert = pR - pR_ref
      
      ML = uL / a
      MR = uR / a
      
      Mbar2 = ( uL**2 + uR**2 ) / ( 2. * a**2 )
      
      !if( pL_pert/pL_ref<1.e-13 .and. pR_pert/pR_ref<1.e-13 )then
      !  p_diff = 0
      !else
      !  p_diff = PR - PL
      !endif
      
      p_diff = PR - PL
      
      Mh     = M4( ML, sp ) + M4( MR, sn ) - Kp * max( 1. - sigma * Mbar2, 0. ) * p_diff / ( rho * a**2 )
      m      = 0.5 * ( cL / sqrtGL + cR / sqrtGR ) * Mh
      
      coefL = sqrtGL * G13L
      coefR = sqrtGR * G13R
      
      p = P5(ML,sp) * coefL * PL + P5(MR,sn) * coefR * PR &
        - Ku * P5(ML,sp) * P5(MR,sn) * ( coefL * rhoL + coefR * rhoR ) * a * ( uR - uL )
      
      !p_pert = P5(ML,sp) * pL_pert + P5(MR,sn) * pR_pert &
      !       - 2. * Ku * P5(ML,sp) * P5(MR,sn) * ( rho * a - rho_ref * a_ref  ) * ( uR - uL )
      p_pert = P5(ML,sp) * pL + P5(MR,sn) * pR &
             - 2. * Ku * P5(ML,sp) * P5(MR,sn) * rho * a * ( uR - uL )
      
      !if( abs(2.*p_pert/(pL+pR)) < 1.e-13 )p_pert = 0
      
    end subroutine AUSM_up_z
    
    function M2(M,signal)
      real(r_kind) :: M2
      real(r_kind) :: M
      real(r_kind) :: signal ! must be 1 or -1
      
      M2 = signal * 0.25 * ( M + signal )**2
      
    end function M2
    
    function M4(M,signal)
      real(r_kind) :: M4
      real(r_kind) :: M
      real(r_kind) :: signal ! must be 1 or -1
      
      real(r_kind),parameter :: beta  = 0.125
      
      if(abs(M)>=1)then
        M4 = 0.5 * ( M + signal * abs(M) )
      else
        M4 = M2( M, signal ) * ( 1. - signal * 16. * beta * M2( M, -signal ) )
      endif
      
    end function M4
    
    function P5(M,signal)
      real(r_kind) :: P5
      real(r_kind) :: M
      real(r_kind) :: signal ! must be 1 or -1
      
      real(r_kind),parameter :: alpha = 0.1875
      
      if(abs(M)>=1)then
        P5 = 0.5 * ( 1. + signal * sign(1.,M) )
      else
        P5 = M2( M, signal ) * ( ( 2. * signal - M ) - signal * 16.*alpha * M * M2( M, -signal ) )
      endif
      
    end function P5
END MODULE spatial_operators_mod

