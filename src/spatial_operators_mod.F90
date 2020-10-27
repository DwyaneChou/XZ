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
  
  real(r_kind), dimension(:,:,:), allocatable :: src_ref ! reference source term
  
  real(r_kind), dimension(:,:  ), allocatable :: P_ref ! Reference pressure
  real(r_kind), dimension(:,:  ), allocatable :: PL_ref! Reconstructed P_ref_(i-1/2,k)
  real(r_kind), dimension(:,:  ), allocatable :: PR_ref! Reconstructed P_ref_(i+1/2,k)
  real(r_kind), dimension(:,:  ), allocatable :: PB_ref! Reconstructed P_ref_(i,k-1/2)
  real(r_kind), dimension(:,:  ), allocatable :: PT_ref! Reconstructed P_ref_(i,k+1/2)
  
    real(r_kind), dimension(:,:  ), allocatable :: sqrtGP ! sqrtG * P
    real(r_kind), dimension(:,:  ), allocatable :: sqrtGG13P ! sqrtG * G13 * P
  
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
      
      allocate(src    (nVar,ids:ide,kds:kde))
      allocate(src_ref(nVar,ids:ide,kds:kde))
      
      allocate(rho_p(     ids:ide,kds:kde))
      
      allocate(P_ref (    ids:ide,kds:kde))
      allocate(PL_ref(    ids:ide,kds:kde))
      allocate(PR_ref(    ids:ide,kds:kde))
      allocate(PB_ref(    ids:ide,kds:kde))
      allocate(PT_ref(    ids:ide,kds:kde))
      
      allocate(sqrtGP    (    ids:ide,kds:kde))
      allocate(sqrtGG13P (    ids:ide,kds:kde))
      
      ! Set reference pressure
      q_ext = ref%q
      call bdy_condition(q_ext,q_ext,ref%q,src)
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
          
          P_ref    (i,k) = calc_pressure(sqrtG(i,k),q_ext(:,i,k))
          sqrtGP   (i,k) = sqrtG(i,k)            * P_ref(i,k)
          sqrtGG13P(i,k) = sqrtG(i,k) * G13(i,k) * P_ref(i,k)
          
          PL_ref(i,k) = calc_pressure(sqrtGL(i,k),qL(:,i,k))
          PR_ref(i,k) = calc_pressure(sqrtGR(i,k),qR(:,i,k))
          PB_ref(i,k) = calc_pressure(sqrtGB(i,k),qB(:,i,k))
          PT_ref(i,k) = calc_pressure(sqrtGT(i,k),qT(:,i,k))
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! Set referenc source term
      src_ref = 0
      do k = kds,kde
        do i = ids,ide
          if(i==ids  )then
            dpdx =-( 3. * sqrtGP(i,k) - 4. * sqrtGP(i+1,k) + sqrtGP(i+2,k) ) / dx
          elseif(i==ide  )then
            dpdx = ( 3. * sqrtGP(i,k) - 4. * sqrtGP(i-1,k) + sqrtGP(i-2,k) ) / dx
          elseif(i==ids+1)then
            dpdx =-( 2. * sqrtGP(i-1,k) + 3. * sqrtGP(i,k) - 6. * sqrtGP(i+1,k) + sqrtGP(i+2,k) ) / ( 6. * dx )
          elseif(i==ide-1)then
            dpdx = ( 2. * sqrtGP(i+1,k) + 3. * sqrtGP(i,k) - 6. * sqrtGP(i-1,k) + sqrtGP(i-2,k) ) / ( 6. * dx )
          else
            dpdx = ( sqrtGP(i-2,k) - 8. * sqrtGP(i-1,k) + 8. * sqrtGP(i+1,k) - sqrtGP(i+2,k) ) / ( 12. * dx )
          endif
          
          if(k==kds  )then
            dpdeta =-( 3. * sqrtGG13P(i,k) - 4. * sqrtGG13P(i,k+1) + sqrtGG13P(i,k+2) ) / deta
          elseif(k==kde  )then
            dpdeta = ( 3. * sqrtGG13P(i,k) - 4. * sqrtGG13P(i,k-1) + sqrtGG13P(i,k-2) ) / deta
          elseif(k==kds+1)then
            dpdeta =-( 2. * sqrtGG13P(i,k-1) + 3. * sqrtGG13P(i,k) - 6. * sqrtGG13P(i,k+1) + sqrtGG13P(i,k+2) ) / ( 6. * deta )
          elseif(k==kde-1)then
            dpdeta = ( 2. * sqrtGG13P(i,k+1) + 3. * sqrtGG13P(i,k) - 6. * sqrtGG13P(i,k-1) + sqrtGG13P(i,k-2) ) / ( 6. * deta )
          else
            dpdeta = ( sqrtGG13P(i,k-2) - 8. * sqrtGG13P(i,k-1) + 8. * sqrtGG13P(i,k+1) - sqrtGG13P(i,k+2) ) / ( 12. * deta )
          endif
          
          !if(i>ids.and.i<ide)then
          !  print*,k,i,dpdx,(sqrtGP(i+1,k)-sqrtGP(i-1,k))/dx/2.
          !elseif(i==ids)then
          !  print*,k,i,dpdx,(sqrtGP(i+1,k)-sqrtGP(i,k))/dx
          !elseif(i==ide)then
          !  print*,k,i,dpdx,(sqrtGP(i,k)-sqrtGP(i-1,k))/dx
          !endif
          
          !if(k>kds.and.k<kde)then
          !  print*,k,i,dpdeta,(sqrtGG13P(i,k+1)-sqrtGG13P(i,k-1))/deta/2.
          !elseif(k==kds)then
          !  print*,k,i,dpdeta,(sqrtGG13P(i,k+1)-sqrtGG13P(i,k))/deta
          !elseif(k==kde)then
          !  print*,k,i,dpdeta,(sqrtGG13P(i,k)-sqrtGG13P(i,k-1))/deta
          !endif
          
          src_ref(2,i,k) = - dpdx - dpdeta
        enddo
      enddo
      
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
          
          if(case_num==1)then
            if(i==ids)qL(2,i,k) = 0
            if(i==ide)qR(2,i,k) = 0
            if(k==kds)qB(3,i,k) = 0
            if(k==kde)qT(3,i,k) = 0
          endif
          
          PL(i,k) = calc_pressure(sqrtGL(i,k),qL(:,i,k))
          PR(i,k) = calc_pressure(sqrtGR(i,k),qR(:,i,k))
          PB(i,k) = calc_pressure(sqrtGB(i,k),qB(:,i,k))
          PT(i,k) = calc_pressure(sqrtGT(i,k),qT(:,i,k))
          
          PL(i,k) = merge( PL(i,k) - PL_ref(i,k), 0., abs( PL(i,k) - PL_ref(i,k) ) / PL_ref(i,k) > 1.e-13 )
          PR(i,k) = merge( PR(i,k) - PR_ref(i,k), 0., abs( PR(i,k) - PR_ref(i,k) ) / PR_ref(i,k) > 1.e-13 )
          PB(i,k) = merge( PB(i,k) - PB_ref(i,k), 0., abs( PB(i,k) - PB_ref(i,k) ) / PB_ref(i,k) > 1.e-13 )
          PT(i,k) = merge( PT(i,k) - PT_ref(i,k), 0., abs( PT(i,k) - PT_ref(i,k) ) / PT_ref(i,k) > 1.e-13 )
          
          FL(:,i,k) = calc_F(sqrtGL(i,k),qL(:,i,k),PL(i,k))
          FR(:,i,k) = calc_F(sqrtGR(i,k),qR(:,i,k),PR(i,k))
          
          HB(:,i,k) = calc_H(sqrtGB(i,k),G13B(i,k),qB(:,i,k),PB(i,k))
          HT(:,i,k) = calc_H(sqrtGT(i,k),G13T(i,k),qT(:,i,k),PT(i,k))
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! Fill boundary
      if(case_num==1)then
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
      elseif(case_num==2)then
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
      endif
      
      rho_p = ( stat%q(1,ids:ide,kds:kde) + stat%q(5,ids:ide,kds:kde) &
              - ref%q (1,ids:ide,kds:kde) - ref%q (5,ids:ide,kds:kde) ) / sqrtG(ids:ide,kds:kde)
      
      where(abs(rho_p)<=1.E-13)rho_p=0.
      
      src(3,ids:ide,kds:kde) = src(3,ids:ide,kds:kde) - sqrtG(ids:ide,kds:kde) * rho_p(ids:ide,kds:kde) * gravity
      
      ! calc x flux
      !$OMP PARALLEL DO PRIVATE(i,im1,eigenvalue_x,maxeigen_x,iVar)
      do k = kds,kde
        do i = ids+1,ide
          im1 = i - 1
          eigenvalue_x(:,1) = calc_eigenvalue_x(sqrtG(im1,k),q_ext(:,im1,k))
          eigenvalue_x(:,2) = calc_eigenvalue_x(sqrtG(i  ,k),q_ext(:,i  ,k))
          
          maxeigen_x = maxval(abs(eigenvalue_x))
          
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
      !$OMP PARALLEL DO PRIVATE(k,km1,eigenvalue_z,maxeigen_z,iVar)
      do i = ids,ide
        do k = kds,kde+1
          km1 = k - 1
          eigenvalue_z(:,1) = calc_eigenvalue_z(sqrtG(i,km1),G13(i,km1),q_ext(:,i,km1))
          eigenvalue_z(:,2) = calc_eigenvalue_z(sqrtG(i  ,k),G13(i  ,k),q_ext(:,i  ,k))
          
          maxeigen_z = maxval(abs(eigenvalue_z))
          
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
      
      ! Fill boundary
      if(case_num==1)then
        ! Fill no-flux bdy
        Fe(1,ids  ,kds:kde) = 0
        Fe(1,ide+1,kds:kde) = 0
        Fe(2,ids  ,kds:kde) = sqrtGL(ids,kds:kde) * PL(ids,kds:kde)
        Fe(2,ide+1,kds:kde) = sqrtGR(ide,kds:kde) * PR(ide,kds:kde)
        Fe(3,ids  ,kds:kde) = 0
        Fe(3,ide+1,kds:kde) = 0
        Fe(4,ids  ,kds:kde) = 0
        Fe(4,ide+1,kds:kde) = 0
        Fe(5,ids  ,kds:kde) = 0
        Fe(5,ide+1,kds:kde) = 0
        
        He(1,ids:ide,kds  ) = 0
        He(1,ids:ide,kde+1) = 0
        He(2,ids:ide,kds  ) = sqrtGB(ids:ide,kds) * G13B(ids:ide,kds) * PB(ids:ide,kds)
        He(2,ids:ide,kde+1) = sqrtGT(ids:ide,kde) * G13T(ids:ide,kde) * PT(ids:ide,kde)
        He(3,ids:ide,kds  ) = PB(ids:ide,kds)
        He(3,ids:ide,kde+1) = PT(ids:ide,kde)
        He(4,ids:ide,kds  ) = 0
        He(4,ids:ide,kde+1) = 0
        He(5,ids:ide,kds  ) = 0
        He(5,ids:ide,kde+1) = 0
      elseif(case_num==2)then
        Fe(1,ids  ,kds:kde) = ref%q(2,ids-1,kds:kde)
        Fe(1,ide+1,kds:kde) = ref%q(2,ide+1,kds:kde)
        Fe(2,ids  ,kds:kde) = ref%q(2,ids-1,kds:kde) * ref%q(2,ids-1,kds:kde) / ref%q(1,ids-1,kds:kde) + sqrtGL(ids,kds:kde) * PL(ids,kds:kde)
        Fe(2,ide+1,kds:kde) = ref%q(2,ide+1,kds:kde) * ref%q(2,ide+1,kds:kde) / ref%q(1,ide+1,kds:kde) + sqrtGR(ide,kds:kde) * PR(ide,kds:kde)
        Fe(3,ids  ,kds:kde) = ref%q(3,ids-1,kds:kde) * ref%q(2,ids-1,kds:kde) / ref%q(1,ids-1,kds:kde)
        Fe(3,ide+1,kds:kde) = ref%q(3,ide+1,kds:kde) * ref%q(2,ide+1,kds:kde) / ref%q(1,ide+1,kds:kde)
        Fe(4,ids  ,kds:kde) = ref%q(4,ids-1,kds:kde) * ref%q(2,ids-1,kds:kde) / ref%q(1,ids-1,kds:kde)
        Fe(4,ide+1,kds:kde) = ref%q(4,ide+1,kds:kde) * ref%q(2,ide+1,kds:kde) / ref%q(1,ide+1,kds:kde)
        Fe(5,ids  ,kds:kde) = 0
        Fe(5,ide+1,kds:kde) = 0
        
        He(1,ids:ide,kds  ) = 0
        He(1,ids:ide,kde+1) = 0
        He(2,ids:ide,kds  ) = G13B(ids:ide,kds) * qB(2,ids:ide,kds)**2 / qB(1,ids:ide,kds) + sqrtGB(ids:ide,kds) * G13B(ids:ide,kds) * PB(ids:ide,kds)
        He(2,ids:ide,kde+1) = G13T(ids:ide,kde) * qT(2,ids:ide,kde)**2 / qT(1,ids:ide,kde) + sqrtGT(ids:ide,kde) * G13T(ids:ide,kde) * PT(ids:ide,kde)
        He(3,ids:ide,kds  ) = PB(ids:ide,kds)
        He(3,ids:ide,kde+1) = PT(ids:ide,kde)
        He(4,ids:ide,kds  ) = 0
        He(4,ids:ide,kde+1) = 0
        He(5,ids:ide,kds  ) = 0
        He(5,ids:ide,kde+1) = 0
      endif
      
      !$OMP PARALLEL DO PRIVATE(i,ip1,kp1)
      do k = kds,kde
        do i = ids,ide
          ip1 = i + 1
          kp1 = k + 1
          tend%q(:,i,k) = - ( ( Fe(:,ip1,k) - Fe(:,i,k) ) / dx + ( He(:,i,kp1) - He(:,i,k) ) / deta ) + src(:,i,k)! + src_ref(:,i,k)
        enddo
      enddo
      !$OMP END PARALLEL DO
      i = 94
      k = kds
      ip1 = i + 1
      kp1 = k + 1
      !print*,-( Fe(2,ip1,k) - Fe(2,i,k) ) / dx,-( He(2,i,kp1) - He(2,i,k) ) / deta, src_ref(2,i,k)
      print*, Fe(2,ip1,k), Fe(2,i,k), He(2,i,kp1), He(2,i,k)
      
    end subroutine spatial_operator
    
    subroutine bdy_condition(q_ext,q,q_ref,src)
      real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(out  ) :: q_ext
      real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(in   ) :: q
      real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(in   ) :: q_ref
      real(r_kind), dimension(nVar,ids:ide,kds:kde), intent(inout) :: src
      
      real(r_kind), dimension(nVar,kds:kde) :: dqx
      real(r_kind), dimension(nVar,ids:ide) :: dqz
      
      integer(i_kind), parameter :: vs = 1
      integer(i_kind), parameter :: ve = 5
      integer(i_kind), parameter :: bdy_width = 20
      real   (r_kind), parameter :: exp_ceof  = 2
      
      integer(i_kind) dir
      integer(i_kind) i,k,iVar
      
      integer(i_kind) il,ir
      integer(i_kind) kt
      
      integer(i_kind) kls,kle ! pure lateral boundary layer indices
      integer(i_kind) its,ite ! pure top boundary layer indices
      
      real(r_kind) :: relax_coef(bdy_width)
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
        q_ext(:,ics:ids-1,:) = q_ref(:,ics:ids-1,:)
        
        ! right
        q_ext(:,ide+1:ice,:) = q_ref(:,ide+1:ice,:)
        
        ! bottom
        q_ext(:,:,kcs:kds-1) = FillValue
        
        ! top
        q_ext(:,:,kde+1:kce) = FillValue
      
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
        
        ! pure zone
        do i = 1,bdy_width
          il = i
          ir = ide-i+1
          kt = kde-i+1
          do iVar = vs,ve
            src(iVar,il     ,kls:kle) = - relax_coef(i) * ( q(iVar,il,kls:kle) - q_ref(iVar,il,kls:kle) )
            src(iVar,ir     ,kls:kle) = - relax_coef(i) * ( q(iVar,ir,kls:kle) - q_ref(iVar,ir,kls:kle) )
            src(iVar,its:ite,kt     ) = - relax_coef(i) * ( q(iVar,its:ite,kt) - q_ref(iVar,its:ite,kt) )
          enddo
        enddo
        
        !overlap zone
        do k = 1,bdy_width
          do i = 1,bdy_width
            il = i
            ir = ide-i+1
            kt = kde-i+1
            do iVar = vs,ve
              src(iVar,il,kt) = - max( relax_coef(i), relax_coef(k) ) * ( q(iVar,il,kt) - q_ref(iVar,il,kt) )
              src(iVar,ir,kt) = - max( relax_coef(i), relax_coef(k) ) * ( q(iVar,ir,kt) - q_ref(iVar,ir,kt) )
            enddo
          enddo
        enddo
      endif
      
      !! Nonreflecting condition
      !kls = kds
      !kle = kde-bdy_width
      !its = ids+bdy_width
      !ite = ide-bdy_width
      !! calculate relax coefficients
      !!max_exp = exp( ( real( bdy_width - 1 ) / real(bdy_width) )**exp_ceof ) - 1.
      !do i = 1,bdy_width
      !  relax_coef(i) = ( real( bdy_width - i + 1 ) / real( bdy_width ) )**4 / dt
      !  !relax_coef(i) = ( exp( ( real( bdy_width - i ) / real(bdy_width) )**exp_ceof ) - 1. ) / ( max_exp * dt )
      !enddo
      !
      !! top only
      !do i = 1,bdy_width
      !  il = i
      !  ir = ide-i+1
      !  kt = kde-i+1
      !  src(vs:ve,ids:ide,kt) = - relax_coef(i) * ( q(vs:ve,ids:ide,kt) - q_ref(vs:ve,ids:ide,kt) )
      !enddo
      !
      !! pure zone
      !do i = 1,bdy_width
      !  il = i
      !  ir = ide-i+1
      !  kt = kde-i+1
      !  do iVar = vs,ve
      !    src(iVar,il     ,kls:kle) = - relax_coef(i) * ( q(iVar,il,kls:kle) - q_ref(iVar,il,kls:kle) )
      !    src(iVar,ir     ,kls:kle) = - relax_coef(i) * ( q(iVar,ir,kls:kle) - q_ref(iVar,ir,kls:kle) )
      !    src(iVar,its:ite,kt     ) = - relax_coef(i) * ( q(iVar,its:ite,kt) - q_ref(iVar,its:ite,kt) )
      !  enddo
      !enddo
      !
      !!overlap zone
      !do k = 1,bdy_width
      !  do i = 1,bdy_width
      !    il = i
      !    ir = ide-i+1
      !    kt = kde-i+1
      !    do iVar = vs,ve
      !      src(iVar,il,kt) = - max( relax_coef(i), relax_coef(k) ) * ( q(iVar,il,kt) - q_ref(iVar,il,kt) )
      !      src(iVar,ir,kt) = - max( relax_coef(i), relax_coef(k) ) * ( q(iVar,ir,kt) - q_ref(iVar,ir,kt) )
      !    enddo
      !  enddo
      !enddo
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
    
    function calc_F(sqrtG,q,pp)
      real(r_kind),dimension(5) :: calc_F
      real(r_kind)              :: sqrtG
      real(r_kind),dimension(5) :: q(5)
      real(r_kind)              :: pp   ! pressure  perturbation
      
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
      calc_F(2) = w2 * u + sqrtG * pp
      calc_F(3) = w3 * u
      calc_F(4) = w4 * u
      calc_F(5) = w5 * u
      
    end function calc_F
    
    function calc_H(sqrtG,G13,q,pp)
      real(r_kind),dimension(5) :: calc_H
      real(r_kind)              :: sqrtG
      real(r_kind)              :: G13
      real(r_kind),dimension(5) :: q(5)
      real(r_kind)              :: pp   ! pressure perturbation
      
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
      
      ww = ( w + sqrtG * G13 * u ) / sqrtG
      
      calc_H(1) = w1 * ww
      calc_H(2) = w2 * ww + sqrtG * G13 * pp
      calc_H(3) = w3 * ww + pp
      calc_H(4) = w4 * ww
      calc_H(5) = w5 * ww
      
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
      
      coef1 = cvd**2*w1**3*w2*w4*(w1 + w5)*(w1 + w5 + eq*w5) + &
            2.*cvd*cvv*w1**2*w2*w4*w5*(w1 + w5)*(w1 + w5 + eq*w5) + &
            cvv**2*w1*w2*w4*w5**2*(w1 + w5)*(w1 + w5 + eq*w5)
      
      coef2 = sqrt( p0*sqrtG*w1**2*w4**2*(w1 + w5)**3*(cpd*w1 + cpv*w5)*(cvd*w1 + cvv*w5)**3* &
      (w1 + w5 + eq*w5)**2*((Rd*w4*(w1 + w5 + eq*w5))/(p0*sqrtG*w1))**((cpd*w1 + cpv*w5)/(cvd*w1 + cvv*w5)) )
      
      coef3 = w1*w4*(w1 + w5)**2*(cvd*w1 + cvv*w5)**2*(w1 + w5 + eq*w5)
      
      calc_eigenvalue_x(1) = w2 / ( w1 + w5 )
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
      
      real(r_kind) sqrtGrho
      real(r_kind) u
      real(r_kind) w
      
      real(r_kind) drhoetadt
      
      real(r_kind) coef1,coef2,coef3
      
      w1 = q(1)
      w2 = q(2)
      w3 = q(3)
      w4 = q(4)
      w5 = q(5)
      
      sqrtGrho = w1 + w5
      u        = w2 / sqrtGrho
      w        = w3 / sqrtGrho
      
      coef1 = cvd**2*w1**3*(G13*sqrtG*w2 + w3)*w4*(w1 + w5)*(w1 + w5 + eq*w5) +      &
            2.*cvd*cvv*w1**2*(G13*sqrtG*w2 + w3)*w4*w5*(w1 + w5)*(w1 + w5 + eq*w5) + &
            cvv**2*w1*(G13*sqrtG*w2 + w3)*w4*w5**2*(w1 + w5)*(w1 + w5 + eq*w5)
      
      coef2 = sqrt( p0*sqrtG*(1 + G13**2*sqrtG**2)*w1**2*                             &
            w4**2*(w1 + w5)**3*(cpd*w1 + cpv*w5)*(cvd*w1 + cvv*w5)**3*                &
           (w1 + w5 + eq*w5)**2*((Rd*w4*(w1 + w5 + eq*w5))/(p0*sqrtG*w1))**((cpd*w1 + &
           cpv*w5)/(cvd*w1 + cvv*w5)) )
      
      coef3 = sqrtG*w1*w4*(w1 + w5)**2*(cvd*w1 + cvv*w5)**2*(w1 + w5 + eq*w5)
      
      drhoetadt = (w + sqrtG*G13*u)/sqrtG
      
      calc_eigenvalue_z(1) = drhoetadt
      calc_eigenvalue_z(2) = drhoetadt
      calc_eigenvalue_z(3) = drhoetadt
      calc_eigenvalue_z(4) = ( coef1 - coef2 ) / coef3
      calc_eigenvalue_z(5) = ( coef1 + coef2 ) / coef3
      
    end function calc_eigenvalue_z
    
END MODULE spatial_operators_mod

