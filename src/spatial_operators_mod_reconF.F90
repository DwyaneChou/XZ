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
  
  real(r_kind), dimension(:,:  ), allocatable :: PF    ! sqrtG * P
  real(r_kind), dimension(:,:  ), allocatable :: PH    ! sqrtG * G13 * P
  real(r_kind), dimension(:,:  ), allocatable :: PFL   ! Reconstructed PF_(i-1/2,k)
  real(r_kind), dimension(:,:  ), allocatable :: PFR   ! Reconstructed PF_(i+1/2,k)
  real(r_kind), dimension(:,:  ), allocatable :: PHB   ! Reconstructed PH_(i,k-1/2)
  real(r_kind), dimension(:,:  ), allocatable :: PHT   ! Reconstructed PH_(i,k+1/2)
  
  real(r_kind), dimension(:,:,:), allocatable :: Fe    ! F on edges of each cell
  real(r_kind), dimension(:,:,:), allocatable :: He    ! H on edges of each cell
  
  real(r_kind), dimension(:,:,:), allocatable :: src   ! source term
  
  real(r_kind), dimension(  :,:), allocatable :: rho_p ! density perturbation
  
  real(r_kind), dimension(:,:  ), allocatable :: P_ref ! reference pressure
  real(r_kind), dimension(:,:  ), allocatable :: PL_ref! Reconstructed P_ref_(i-1/2,k)
  real(r_kind), dimension(:,:  ), allocatable :: PR_ref! Reconstructed P_ref_(i+1/2,k)
  real(r_kind), dimension(:,:  ), allocatable :: PB_ref! Reconstructed P_ref_(i,k-1/2)
  real(r_kind), dimension(:,:  ), allocatable :: PT_ref! Reconstructed P_ref_(i,k+1/2)
  
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
      
      allocate(F    (nVar,ics:ice,kcs:kce))
      allocate(H    (nVar,ics:ice,kcs:kce))
      allocate(P    (     ics:ice,kcs:kce))
      
      allocate(FL   (nVar,ils:ile,kls:kle))
      allocate(FR   (nVar,irs:ire,krs:kre))
      
      allocate(HB   (nVar,ibs:ibe,kbs:kbe))
      allocate(HT   (nVar,its:ite,kts:kte))
      
      allocate(PL   (     ils:ile,kls:kle))
      allocate(PR   (     irs:ire,krs:kre))
      allocate(PB   (     ibs:ibe,kbs:kbe))
      allocate(PT   (     its:ite,kts:kte))
      
      allocate(PF   (     ids:ide,kds:kde))
      allocate(PH   (     ids:ide,kds:kde))
      allocate(PFL  (     ids:ide,kds:kde))
      allocate(PFR  (     ids:ide,kds:kde))
      allocate(PHB  (     ids:ide,kds:kde))
      allocate(PHT  (     ids:ide,kds:kde))
      
      allocate(Fe   (nVar,ids:ide+1,kds:kde  ))
      allocate(He   (nVar,ids:ide  ,kds:kde+1))
      
      allocate(src  (nVar,ids:ide,kds:kde))
      
      allocate(rho_p (    ids:ide,kds:kde))
      
      allocate(P_ref (    ics:ice,kcs:kce))
      allocate(PL_ref(    ids:ide,kds:kde))
      allocate(PR_ref(    ids:ide,kds:kde))
      allocate(PB_ref(    ids:ide,kds:kde))
      allocate(PT_ref(    ids:ide,kds:kde))
      
      ! Set reference pressure
      q_ext = ref%q
      call bdy_condition(q_ext,q_ext,ref%q,src)
      
      P     = FillValue
      P_ref = FillValue
      PF    = FillValue
      PH    = FillValue
      F     = FillValue
      H     = FillValue
      
      do k = kds,kde
        do i = ids,ide
          P_ref(i,k) = calc_pressure(sqrtG(i,k),q_ext(:,i,k))
        enddo
      enddo
      
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
          ! x-dir
          q_weno = P_ref(im2:ip2,k)
          
          dir = -1
          call WENO_limiter(PL_ref(i,k),q_weno,dir)
          dir = 1
          call WENO_limiter(PR_ref(i,k),q_weno,dir)
          
          ! z-dir
          q_weno = P_ref(i,km2:kp2)
          
          dir = -1
          call WENO_limiter(PB_ref(i,k),q_weno,dir)
          dir = 1
          call WENO_limiter(PT_ref(i,k),q_weno,dir)
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
      
      real(r_kind) dFe
      real(r_kind) dHe
      
      ! copy stat
      q_ext = stat%q
      
      ! initialize source terms
      src = 0.
      
      ! Set no flux and nonreflecting boundary condition
      call bdy_condition(q_ext,stat%q,ref%q,src)
      
      do k = kds,kde
        do i = ids,ide
          P(  i,k) = calc_pressure(sqrtG(i,k),q_ext(:,i,k))
          F(:,i,k) = calc_F(sqrtG(i,k),q_ext(:,i,k),P(i,k))
          H(:,i,k) = calc_H(sqrtG(i,k),G13(i,k),q_ext(:,i,k),P(i,k),P_ref(i,k))
          
          PF(i,k) = sqrtG(i,k) * P(i,k)
          PH(i,k) = sqrtG(i,k) * G13(i,k) * P(i,k)
        enddo
      enddo
      
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
            ! Reconstruct q
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
          
            ! Reconstruct F and H
            ! x-dir
            q_weno = F(iVar,im2:ip2,k)
            
            dir = -1
            call WENO_limiter(FL(iVar,i,k),q_weno,dir)
            dir = 1
            call WENO_limiter(FR(iVar,i,k),q_weno,dir)
            
            ! z-dir
            q_weno = H(iVar,i,km2:kp2)
            
            dir = -1
            call WENO_limiter(HB(iVar,i,k),q_weno,dir)
            dir = 1
            call WENO_limiter(HT(iVar,i,k),q_weno,dir)
          enddo
          
          ! Reconstruct P
          ! x-dir
          q_weno = P(im2:ip2,k)
          
          dir = -1
          call WENO_limiter(PL(i,k),q_weno,dir)
          dir = 1
          call WENO_limiter(PR(i,k),q_weno,dir)
          
          ! z-dir
          q_weno = P(i,km2:kp2)
          
          dir = -1
          call WENO_limiter(PB(i,k),q_weno,dir)
          dir = 1
          call WENO_limiter(PT(i,k),q_weno,dir)
          
          
          ! Reconstruct PF and PH
          ! x-dir
          q_weno = PF(im2:ip2,k)
          
          dir = -1
          call WENO_limiter(PFL(i,k),q_weno,dir)
          dir = 1
          call WENO_limiter(PFR(i,k),q_weno,dir)
          
          ! z-dir
          q_weno = PH(i,km2:kp2)
          
          dir = -1
          call WENO_limiter(PHB(i,k),q_weno,dir)
          dir = 1
          call WENO_limiter(PHT(i,k),q_weno,dir)
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      if(case_num==1)then
        qL(2,ids,:) = 0
        qR(2,ide,:) = 0
        qB(3,:,kds) = 0
        qT(3,:,kde) = 0
        
        do k = kds,kde
          FL(1,ids,k) = 0
          FL(2,ids,k) = PFL(ids,k)
          FL(3,ids,k) = 0
          FL(4,ids,k) = 0
          FL(5,ids,k) = 0
          FR(1,ide,k) = 0
          FR(2,ide,k) = PFR(ide,k)
          FR(3,ide,k) = 0
          FR(4,ide,k) = 0
          FR(5,ide,k) = 0
        enddo
        
        do i = ids,ide
          HB(1,i,kds) = 0
          HB(2,i,kds) = 0
          if(abs(PB(i,kds) - PB_ref(i,kds))/PB_ref(i,kds)<1.e-13)then
            HB(3,i,kds) = 0
          else
            HB(3,i,kds) = PB(i,kds) - PB_ref(i,kds)
          endif
          HB(4,i,kds) = 0
          HB(5,i,kds) = 0
          
          HT(1,i,kde) = 0
          HT(2,i,kde) = 0
          if(abs(PT(i,kde) - PT_ref(i,kde))/PT_ref(i,kde)<1.e-13)then
            HT(3,i,kde) = 0
          else
            HT(3,i,kde) = PT(i,kde) - PT_ref(i,kde)
          endif
          HT(4,i,kde) = 0
          HT(5,i,kde) = 0
        enddo
      elseif(case_num==2)then
        qL(2,ids,:) = ref%q(2,ids,kds:kde)
        qR(2,ide,:) = ref%q(2,ide,kds:kde)
        qB(3,:,kds) = 0
        qT(3,:,kde) = 0
          
        do k = kds,kde
          FL(1,ids,k) = ref%q(2,ids,k)
          FL(2,ids,k) = ref%q(2,ids,k)**2 / ref%q(1,ids,k) + PFL(ids,k)
          FL(3,ids,k) = 0
          FL(4,ids,k) = ref%q(4,ids,k) * ref%q(2,ids,k) / ref%q(1,ids,k)
          FL(5,ids,k) = 0
          FR(1,ide,k) = ref%q(2,ide,k)
          FR(2,ide,k) = ref%q(2,ide,k)**2 / ref%q(1,ide,k) + PFR(ide,k)
          FR(3,ide,k) = 0
          FR(4,ide,k) = ref%q(4,ide,k) * ref%q(2,ide,k) / ref%q(1,ide,k)
          FR(5,ide,k) = 0
        enddo
        
        do i = ids,ide
          HB(:,i,kds) = calc_H(sqrtGB(i,kds),G13B(i,kds),qB(:,i,kds),PB(i,kds),PB_ref(i,kds))
          HB(3,i,kds) = PB(i,kds) - PB_ref(i,kds)
          HT(:,i,kde) = calc_H(sqrtGT(i,kde),G13T(i,kde),qT(:,i,kde),PT(i,kde),PT_ref(i,kde))
          !HT(1,i,kde) = 0
          !HT(2,i,kde) = 0
          !HT(3,i,kde) = PT(i,kde) - PT_ref(i,kde)
          !HT(4,i,kde) = 0
          !HT(5,i,kde) = 0
        enddo
      endif
      
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
      
      ! calc x flux
      !$OMP PARALLEL DO PRIVATE(i,im1,eigenvalue_x,maxeigen_x,iVar)
      do k = kds,kde
        do i = ids,ide+1
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
            !if(iVar==2.and.i==109.and.k==1)print*,Fe(iVar,i,k),FL(iVar,i,k),FR(iVar,im1,k),maxeigen_x,qL(iVar,i,k),qR(iVar,im1,k)
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
            !if(iVar==3.and.i==109.and.k==2)print*,He(iVar,i,k),HB(iVar,i,k),HT(iVar,i,km1),maxeigen_z,qB(iVar,i,k),qT(iVar,i,km1)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      !$OMP PARALLEL DO PRIVATE(i,iVar,ip1,kp1,dFe,dHe)
      do k = kds,kde
        do i = ids,ide
          ip1 = i + 1
          kp1 = k + 1
          do iVar = 1,nVar
            dFe = Fe(iVar,ip1,k) - Fe(iVar,i,k)
            dHe = He(iVar,i,kp1) - He(iVar,i,k)
            
            if(abs(dFe/Fe(iVar,i,k))<1.e-15)dFe=0
            if(abs(dHe/He(iVar,i,k))<1.e-15)dHe=0
            
            tend%q(iVar,i,k) = - ( dFe / dx + dHe / deta ) + src(iVar,i,k)
          enddo
          !if(k==kds)then
          !  iVar = 2
          !  print*,i,tend%q(iVar,i,k),( Fe(iVar,ip1,k) - Fe(iVar,i,k) ) / dx, ( He(iVar,i,kp1) - He(iVar,i,k) ) / deta, src(iVar,i,k), Fe(iVar,ip1,k), Fe(iVar,i,k)
          !endif
        enddo
      enddo
      !$OMP END PARALLEL DO
      !stop 'Check hydrostatic'
    end subroutine spatial_operator
    
    subroutine bdy_condition(q_ext,q,q_ref,src)
      real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(out  ) :: q_ext
      real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(in   ) :: q
      real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(in   ) :: q_ref
      real(r_kind), dimension(nVar,ids:ide,kds:kde), intent(inout) :: src
      
      real(r_kind), dimension(nVar,kds:kde) :: dqx
      real(r_kind), dimension(nVar,ids:ide) :: dqz
      
      integer(i_kind), parameter :: vs = 2
      integer(i_kind), parameter :: ve = 3
      integer(i_kind), parameter :: bdy_width = 30
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
        
        !! top only
        !do i = 1,bdy_width
        !  il = i
        !  ir = ide-i+1
        !  kt = kde-i+1
        !  do iVar = vs,ve
        !    src(iVar,ids:ide,kt) = - relax_coef(i) * ( q(iVar,ids:ide,kt) - q_ref(iVar,ids:ide,kt) )
        !  enddo
        !enddo
        
        !! lateral only
        !do i = 1,bdy_width
        !  il = i
        !  ir = ide-i+1
        !  kt = kde-i+1
        !  do iVar = vs,ve
        !    src(iVar,il,kds:kde) = - relax_coef(i) * ( q(iVar,il,kds:kde) - q_ref(iVar,il,kds:kde) )
        !    src(iVar,ir,kds:kde) = - relax_coef(i) * ( q(iVar,ir,kds:kde) - q_ref(iVar,ir,kds:kde) )
        !  enddo
        !enddo
        
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
    
    function calc_F(sqrtG,q,P)
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
      !p        = p0*((Rd*w4*(w1 + w5 + eq*w5))/(p0*sqrtG*w1))**((cpd*w1 + cpv*w5)/(cvd*w1 + cvv*w5))
      
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
      !p        = p0*((Rd*w4*(w1 + w5 + eq*w5))/(p0*sqrtG*w1))**((cpd*w1 + cpv*w5)/(cvd*w1 + cvv*w5))
      p_pert   = p - p_ref
      if(abs(p_pert)/p_ref<1.e-14)p_pert=0
      
      ww = w / sqrtG + G13 * u
      
      calc_H(1) = w1 * ww
      calc_H(2) = w2 * ww + sqrtG * G13 * p
      calc_H(3) = w3 * ww + p_pert
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
      
      if(any(q==FillValue))then
        calc_eigenvalue_x = 0
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
        
        calc_eigenvalue_x(1) = w2 / ( w1 + w5 )
        calc_eigenvalue_x(2) = calc_eigenvalue_x(1)
        calc_eigenvalue_x(3) = calc_eigenvalue_x(1)
        calc_eigenvalue_x(4) = ( coef1 - coef2 ) / coef3
        calc_eigenvalue_x(5) = ( coef1 + coef2 ) / coef3
      endif
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
      
      if(any(q==FillValue))then
        calc_eigenvalue_z = 0
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
        
        calc_eigenvalue_z(1) = drhoetadt
        calc_eigenvalue_z(2) = drhoetadt
        calc_eigenvalue_z(3) = drhoetadt
        calc_eigenvalue_z(4) = ( coef1 - coef2 ) / coef3
        calc_eigenvalue_z(5) = ( coef1 + coef2 ) / coef3
      endif
    end function calc_eigenvalue_z
    
END MODULE spatial_operators_mod

