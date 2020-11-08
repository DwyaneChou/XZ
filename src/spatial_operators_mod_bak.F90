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
  
  real(r_kind), dimension(  :,:), allocatable :: rho_p ! density perturbation
  
  real(r_kind), dimension(:,:  ), allocatable :: PL_ref! Reconstructed P_ref_(i-1/2,k)
  real(r_kind), dimension(:,:  ), allocatable :: PR_ref! Reconstructed P_ref_(i+1/2,k)
  real(r_kind), dimension(:,:  ), allocatable :: PB_ref! Reconstructed P_ref_(i,k-1/2)
  real(r_kind), dimension(:,:  ), allocatable :: PT_ref! Reconstructed P_ref_(i,k+1/2)
  
  real(r_kind), dimension(:,:), allocatable :: eig_x
  real(r_kind), dimension(:,:), allocatable :: eig_z
  
    contains
    subroutine init_spatial_operator
      integer(i_kind) dir
      integer(i_kind) i,k,iVar
      
      real(r_kind), dimension(5) :: q_weno
      
      integer(i_kind) ip1,im1
      integer(i_kind) kp1,km1
      integer(i_kind) ip2,im2
      integer(i_kind) kp2,km2
      
      real(r_kind) dpdx
      real(r_kind) dpdeta
      
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
      
      allocate(rho_p (    ids:ide,kds:kde))
      
      allocate(PL_ref(    ids:ide,kds:kde))
      allocate(PR_ref(    ids:ide,kds:kde))
      allocate(PB_ref(    ids:ide,kds:kde))
      allocate(PT_ref(    ids:ide,kds:kde))
      
      allocate(eig_x(     ics:ice,kcs:kce))
      allocate(eig_z(     ics:ice,kcs:kce))
      
      ! Set reference pressure
      q_ext = ref%q
      call bdy_condition(q_ext,q_ext,ref%q,src)
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
      
      ! copy stat
      q_ext = stat%q
      
      ! initialize source terms
      src = 0.
      
      ! Set no flux and nonreflecting boundary condition
      call bdy_condition(q_ext,stat%q,ref%q,src)
      
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
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      if(case_num==1)then
        qL(2,ids,:) = 0
        qR(2,ide,:) = 0
        qB(3,:,kds) = 0
        qT(3,:,kde) = 0
      elseif(case_num==2)then
        qL(2,ids,:) = ref%q(2,ids,kds:kde)
        qR(2,ide,:) = ref%q(2,ide,kds:kde)
        !qL(2,ids,:) = ref%q(2,ids,kds:kde) / ref%q(1,ids,kds:kde) * q_ext(1,ids,kds:kde)
        !qR(2,ide,:) = ref%q(2,ide,kds:kde) / ref%q(1,ide,kds:kde) * q_ext(1,ids,kds:kde)
        !qB(3,:,kds) = -sqrtGB(:,kds) * G13B(:,kds) * qB(2,:,kds)
        !qT(3,:,kde) = -sqrtGT(:,kde) * G13T(:,kde) * qT(2,:,kde)
      endif
      
      !$OMP PARALLEL DO PRIVATE(i)
      do k = kds,kde
        do i = ids,ide
          PL(i,k) = calc_pressure(sqrtGL(i,k),qL(:,i,k))
          PR(i,k) = calc_pressure(sqrtGR(i,k),qR(:,i,k))
          PB(i,k) = calc_pressure(sqrtGB(i,k),qB(:,i,k))
          PT(i,k) = calc_pressure(sqrtGT(i,k),qT(:,i,k))
          
          FL(:,i,k) = calc_F(sqrtGL(i,k),qL(:,i,k),PL(i,k))
          FR(:,i,k) = calc_F(sqrtGR(i,k),qR(:,i,k),PR(i,k))
          
          HB(:,i,k) = calc_H(sqrtGB(i,k),G13B(i,k),qB(:,i,k),PB(i,k),PB_ref(i,k))
          HT(:,i,k) = calc_H(sqrtGT(i,k),G13T(i,k),qT(:,i,k),PT(i,k),PT_ref(i,k))
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! Fill boundary
      ! left boundary
      qR(:,irs,:) = qL(:,ids,:)
      FR(:,irs,:) = FL(:,ids,:)
      
      ! right boundary
      qL(:,ile,:) = qR(:,ide,:)
      FL(:,ile,:) = FR(:,ide,:)
      
      ! bottom boundary
      qT(:,:,kts) = qB(:,:,kds)
      HT(:,:,kts) = HB(:,:,kds)
      
      ! top boundary
      qB(:,:,kbe) = qT(:,:,kde)
      HB(:,:,kbe) = HT(:,:,kde)
      
      rho_p = ( stat%q(1,ids:ide,kds:kde) + stat%q(5,ids:ide,kds:kde) &
              - ref%q (1,ids:ide,kds:kde) - ref%q (5,ids:ide,kds:kde) ) / sqrtG(ids:ide,kds:kde)
      
      where(abs(rho_p)<=1.E-13)rho_p=0.
      
      src(3,ids:ide,kds:kde) = src(3,ids:ide,kds:kde) - sqrtG(ids:ide,kds:kde) * rho_p(ids:ide,kds:kde) * gravity
      
      ! Calculate eigenvalues x-dir
      !$OMP PARALLEL DO PRIVATE(i)
      do k = kds,kde
        do i = ids-1,ide+1
          eig_x(i,k) = calc_eigenvalue_x(sqrtG(i,k),q_ext(:,i,k))
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! Calculate eigenvalues z-dir
      !$OMP PARALLEL DO PRIVATE(i)
      do k = kds-1,kde+1
        do i = ids,ide
          eig_z(i,k) = calc_eigenvalue_z(sqrtG(i,k),G13(i,k),q_ext(:,i,k))
        enddo
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
      ! No flux boundary
      He(1,ids:ide,kds) = 0
      He(2,ids:ide,kds) = sqrtGB(ids:ide,kds) * G13B(ids:ide,kds) *  PB(ids:ide,kds)
      He(3,ids:ide,kds) = PB(ids:ide,kds) - PB_ref(ids:ide,kds)
      He(4,ids:ide,kds) = 0
      He(5,ids:ide,kds) = 0
      
      He(1,ids:ide,kde) = 0
      He(2,ids:ide,kde) = sqrtGT(ids:ide,kde) * G13T(ids:ide,kde) *  PT(ids:ide,kde)
      He(3,ids:ide,kde) = PT(ids:ide,kde) - PT_ref(ids:ide,kde)
      He(4,ids:ide,kde) = 0
      He(5,ids:ide,kde) = 0
      
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
    
    subroutine bdy_condition(q_ext,q,q_ref,src)
      real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(out  ) :: q_ext
      real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(in   ) :: q
      real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(in   ) :: q_ref
      real(r_kind), dimension(nVar,ids:ide,kds:kde), intent(inout) :: src
      
      integer(i_kind), parameter :: vs = 1
      integer(i_kind), parameter :: ve = 5
      integer(i_kind), parameter :: bdy_widthT = 40
      integer(i_kind), parameter :: bdy_widthL = 60
      integer(i_kind), parameter :: bdy_widthR = 60
      real   (r_kind), parameter :: exp_ceof  = 2
      
      integer(i_kind) dir
      integer(i_kind) i,k,iVar
      
      integer(i_kind) il,ir
      integer(i_kind) kt
      
      integer(i_kind) kls,kle ! pure lateral boundary layer indices
      integer(i_kind) its,ite ! pure top boundary layer indices
      
      real(r_kind) :: relax_coefT(bdy_widthT)
      real(r_kind) :: relax_coefL(bdy_widthL)
      real(r_kind) :: relax_coefR(bdy_widthR)
      real(r_kind) :: max_exp
      
      if(case_num==1)then
        ! No-flux
        ! left
        q_ext(:,ics:ids-1,:) = FillValue
        
        ! right
        q_ext(:,ide+1:ice,:) = FillValue
        
        ! bottom
        q_ext(:,:,kcs:kds-1) = FillValue
        
        ! top
        q_ext(:,:,kde+1:kce) = FillValue
        
      elseif(case_num==2)then
        ! left
        q_ext(:,ics:ids-1,:) = FillValue
        
        ! right
        q_ext(:,ide+1:ice,:) = FillValue
        
        !! left
        !q_ext(:,ics:ids-1,:) = q_ref(:,ics:ids-1,:)
        !
        !! right
        !q_ext(:,ide+1:ice,:) = q_ref(:,ide+1:ice,:)
        
        ! bottom
        q_ext(:,:,kcs:kds-1) = FillValue
        
        ! top
        q_ext(:,:,kde+1:kce) = FillValue
      
        ! Nonreflecting condition
        kls = kds
        kle = kde-bdy_widthL
        its = ids+bdy_widthT
        ite = ide-bdy_widthT
        ! calculate relax coefficients
        max_exp = exp( ( real( bdy_widthT - 1 ) / real(bdy_widthT) )**exp_ceof ) - 1.
        do i = 1,bdy_widthT
          !relax_coefT(i) = ( real( bdy_widthT - i + 1 ) / real( bdy_widthT ) )**4 / dt
          relax_coefT(i) = ( exp( ( real( bdy_widthT - i ) / real(bdy_widthT) )**exp_ceof ) - 1. ) / ( max_exp * dt )
        enddo
        
        max_exp = exp( ( real( bdy_widthL - 1 ) / real(bdy_widthL) )**exp_ceof ) - 1.
        do i = 1,bdy_widthL
          !relax_coefL(i) = ( real( bdy_widthL - i + 1 ) / real( bdy_widthL ) )**4 / dt
          relax_coefL(i) = ( exp( ( real( bdy_widthL - i ) / real(bdy_widthL) )**exp_ceof ) - 1. ) / ( max_exp * dt )
        enddo
        
        max_exp = exp( ( real( bdy_widthR - 1 ) / real(bdy_widthR) )**exp_ceof ) - 1.
        do i = 1,bdy_widthR
          !relax_coefR(i) = ( real( bdy_widthR - i + 1 ) / real( bdy_widthR ) )**4 / dt
          relax_coefR(i) = ( exp( ( real( bdy_widthR - i ) / real(bdy_widthR) )**exp_ceof ) - 1. ) / ( max_exp * dt )
        enddo
        
        do i = 1,bdy_widthL
          il = i
          ir = ide-i+1
          kt = kde-i+1
          do iVar = vs,ve
            src(iVar,il     ,kds:kde) = - relax_coefL(i) * ( q(iVar,il,kds:kde) - q_ref(iVar,il,kds:kde) )
            src(iVar,ir     ,kds:kde) = - relax_coefR(i) * ( q(iVar,ir,kds:kde) - q_ref(iVar,ir,kds:kde) )
          enddo
        enddo
        
        do i = 1,bdy_widthT
          il = i
          ir = ide-i+1
          kt = kde-i+1
          do iVar = vs,ve
            src(iVar,its:ite,kt     ) = - relax_coefT(i) * ( q(iVar,its:ite,kt) - q_ref(iVar,its:ite,kt) )
          enddo
        enddo
        
        !overlap zone
        do k = 1,bdy_widthT
          do i = 1,bdy_widthL
            il = i
            ir = ide-i+1
            kt = kde-i+1
            do iVar = vs,ve
              src(iVar,il,kt) = - max( relax_coefL(i), relax_coefT(k) ) * ( q(iVar,il,kt) - q_ref(iVar,il,kt) )
              src(iVar,ir,kt) = - max( relax_coefR(i), relax_coefT(k) ) * ( q(iVar,ir,kt) - q_ref(iVar,ir,kt) )
            enddo
          enddo
        enddo
      endif
      
    end subroutine bdy_condition
    
    !subroutine bdy_condition(q_ext,q,q_ref,src)
    !  real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(out  ) :: q_ext
    !  real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(in   ) :: q
    !  real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(in   ) :: q_ref
    !  real(r_kind), dimension(nVar,ids:ide,kds:kde), intent(inout) :: src
    !  
    !  integer(i_kind), parameter :: vs = 1
    !  integer(i_kind), parameter :: ve = 5
    !  integer(i_kind), parameter :: bdy_width = 40
    !  real   (r_kind), parameter :: exp_ceof  = 2
    !  
    !  integer(i_kind) dir
    !  integer(i_kind) i,k,iVar
    !  
    !  integer(i_kind) il,ir
    !  integer(i_kind) kt
    !  
    !  integer(i_kind) kls,kle ! pure lateral boundary layer indices
    !  integer(i_kind) its,ite ! pure top boundary layer indices
    !  
    !  real(r_kind) :: relax_coef(bdy_width)
    !  real(r_kind) :: max_exp
    !  
    !  if(case_num==1)then
    !    ! No-flux
    !    ! left
    !    q_ext(:,ics:ids-1,:) = FillValue
    !    
    !    ! right
    !    q_ext(:,ide+1:ice,:) = FillValue
    !    
    !    ! bottom
    !    q_ext(:,:,kcs:kds-1) = FillValue
    !    
    !    ! top
    !    q_ext(:,:,kde+1:kce) = FillValue
    !    
    !  elseif(case_num==2)then
    !    ! left
    !    q_ext(:,ics:ids-1,:) = FillValue
    !    
    !    ! right
    !    q_ext(:,ide+1:ice,:) = FillValue
    !    
    !    !! left
    !    !q_ext(:,ics:ids-1,:) = q_ref(:,ics:ids-1,:)
    !    !
    !    !! right
    !    !q_ext(:,ide+1:ice,:) = q_ref(:,ide+1:ice,:)
    !    
    !    ! bottom
    !    q_ext(:,:,kcs:kds-1) = FillValue
    !    
    !    ! top
    !    q_ext(:,:,kde+1:kce) = FillValue
    !  
    !    ! Nonreflecting condition
    !    kls = kds
    !    kle = kde-bdy_width
    !    its = ids+bdy_width
    !    ite = ide-bdy_width
    !    ! calculate relax coefficients
    !    max_exp = exp( ( real( bdy_width - 1 ) / real(bdy_width) )**exp_ceof ) - 1.
    !    do i = 1,bdy_width
    !      !relax_coef(i) = ( real( bdy_width - i + 1 ) / real( bdy_width ) )**4 / dt
    !      relax_coef(i) = ( exp( ( real( bdy_width - i ) / real(bdy_width) )**exp_ceof ) - 1. ) / ( max_exp * dt )
    !    enddo
    !    
    !    do i = 1,bdy_width
    !      il = i
    !      ir = ide-i+1
    !      kt = kde-i+1
    !      do iVar = vs,ve
    !        !src(iVar,its:ite,kt     ) = - relax_coef(i) * ( q(iVar,its:ite,kt) - q_ref(iVar,its:ite,kt) / q_ref(1,its:ite,kt) * q(1,its:ite,kt) )
    !        src(iVar,its:ite,kt     ) = - relax_coef(i) * ( q(iVar,its:ite,kt) - q_ref(iVar,its:ite,kt) )
    !        src(iVar,il     ,kds:kde) = - relax_coef(i) * ( q(iVar,il,kds:kde) - q_ref(iVar,il,kds:kde) )
    !        src(iVar,ir     ,kds:kde) = - relax_coef(i) * ( q(iVar,ir,kds:kde) - q_ref(iVar,ir,kds:kde) )
    !      enddo
    !    enddo
    !    
    !    !! lateral only
    !    !do i = 1,bdy_width
    !    !  il = i
    !    !  ir = ide-i+1
    !    !  kt = kde-i+1
    !    !  do iVar = vs,ve
    !    !    src(iVar,il,kds:kde) = - relax_coef(i) * ( q(iVar,il,kds:kde) - q_ref(iVar,il,kds:kde) )
    !    !    src(iVar,ir,kds:kde) = - relax_coef(i) * ( q(iVar,ir,kds:kde) - q_ref(iVar,ir,kds:kde) )
    !    !  enddo
    !    !enddo
    !    
    !    !! pure zone
    !    !do i = 1,bdy_width
    !    !  il = i
    !    !  ir = ide-i+1
    !    !  kt = kde-i+1
    !    !  do iVar = vs,ve
    !    !    !src(iVar,il     ,kds:kde) = - relax_coef(i) * ( q(iVar,il,kds:kde) - q_ref(iVar,il,kds:kde) )
    !    !    src(iVar,ir     ,kds:kde) = - relax_coef(i) * ( q(iVar,ir,kds:kde) - q_ref(iVar,ir,kds:kde) )
    !    !    src(iVar,its:ite,kt     ) = - relax_coef(i) * ( q(iVar,its:ite,kt) - q_ref(iVar,its:ite,kt) )
    !    !  enddo
    !    !  !src(2,il,kds:kde) = 0
    !    !enddo
    !    
    !    !overlap zone
    !    do k = 1,bdy_width
    !      do i = 1,bdy_width
    !        il = i
    !        ir = ide-i+1
    !        kt = kde-i+1
    !        do iVar = vs,ve
    !          src(iVar,il,kt) = - max( relax_coef(i), relax_coef(k) ) * ( q(iVar,il,kt) - q_ref(iVar,il,kt) )
    !          src(iVar,ir,kt) = - max( relax_coef(i), relax_coef(k) ) * ( q(iVar,ir,kt) - q_ref(iVar,ir,kt) )
    !        enddo
    !      enddo
    !    enddo
    !  endif
    !  
    !end subroutine bdy_condition
    
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

