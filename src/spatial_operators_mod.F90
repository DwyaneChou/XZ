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
  
  real(r_kind), dimension(:,:,:), allocatable :: q_ext ! Extended forecast variables
  
  real(r_kind), dimension(:,:,:), allocatable :: qL    ! Reconstructed q_(i-1/2,k)
  real(r_kind), dimension(:,:,:), allocatable :: qR    ! Reconstructed q_(i+1/2,k)
  real(r_kind), dimension(:,:,:), allocatable :: qB    ! Reconstructed q_(i,k-1/2)
  real(r_kind), dimension(:,:,:), allocatable :: qT    ! Reconstructed q_(i,k+1/2)
  
  real(r_kind), dimension(:,:,:), allocatable :: F
  real(r_kind), dimension(:,:,:), allocatable :: H
  
  real(r_kind), dimension(:,:,:), allocatable :: FL    ! Reconstructed F_(i-1/2,k)
  real(r_kind), dimension(:,:,:), allocatable :: FR    ! Reconstructed F_(i+1/2,k)
  real(r_kind), dimension(:,:,:), allocatable :: HB    ! Reconstructed H_(i,k-1/2)
  real(r_kind), dimension(:,:,:), allocatable :: HT    ! Reconstructed H_(i,k+1/2)
  
  real(r_kind), dimension(  :,:), allocatable :: P
  real(r_kind), dimension(  :,:), allocatable :: P_ref ! reference pressure
  
  real(r_kind), dimension(:,:  ), allocatable :: PL    ! Reconstructed P_(i-1/2,k)
  real(r_kind), dimension(:,:  ), allocatable :: PR    ! Reconstructed P_(i+1/2,k)
  real(r_kind), dimension(:,:  ), allocatable :: PB    ! Reconstructed P_(i,k-1/2)
  real(r_kind), dimension(:,:  ), allocatable :: PT    ! Reconstructed P_(i,k+1/2)
  
  real(r_kind), dimension(:,:  ), allocatable :: PL_ref! Reconstructed P_ref_(i-1/2,k)
  real(r_kind), dimension(:,:  ), allocatable :: PR_ref! Reconstructed P_ref_(i+1/2,k)
  real(r_kind), dimension(:,:  ), allocatable :: PB_ref! Reconstructed P_ref_(i,k-1/2)
  real(r_kind), dimension(:,:  ), allocatable :: PT_ref! Reconstructed P_ref_(i,k+1/2)
  
  real(r_kind), dimension(  :,:), allocatable :: w_eta ! deta/dt
  
  real(r_kind), dimension(:,:,:), allocatable :: Fe    ! F on edges of each cell
  real(r_kind), dimension(:,:,:), allocatable :: He    ! H on edges of each cell
  
  real(r_kind), dimension(:,:,:), allocatable :: src   ! source term
  
  real(r_kind), dimension(  :,:), allocatable :: rho_p ! density perturbation
      
  real(r_kind), dimension(  :,:), allocatable :: eig_x
  real(r_kind), dimension(  :,:), allocatable :: eig_z
      
  real(r_kind), dimension(:,:), allocatable :: relax_coef ! Relax coefficient of Rayleigh damping
  
  real(r_kind) maxeigen_x
  real(r_kind) maxeigen_z
  
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
      
      allocate(F    (nVar,ics:ice,kcs:kce))
      allocate(H    (nVar,ics:ice,kcs:kce))
      
      allocate(FL   (nVar,ics:ice,kcs:kce))
      allocate(FR   (nVar,ics:ice,kcs:kce))
      allocate(HB   (nVar,ics:ice,kcs:kce))
      allocate(HT   (nVar,ics:ice,kcs:kce))
      
      allocate(P    (     ics:ice,kcs:kce))
      allocate(PL   (     ics:ice,kcs:kce))
      allocate(PR   (     ics:ice,kcs:kce))
      allocate(PB   (     ics:ice,kcs:kce))
      allocate(PT   (     ics:ice,kcs:kce))
      
      allocate(P_ref (    ics:ice,kcs:kce))
      allocate(PL_ref(    ics:ice,kcs:kce))
      allocate(PR_ref(    ics:ice,kcs:kce))
      allocate(PB_ref(    ics:ice,kcs:kce))
      allocate(PT_ref(    ics:ice,kcs:kce))
      
      allocate(w_eta(     ics:ice,kcs:kce))
      
      allocate(Fe   (nVar,ids:ide+1,kds:kde  ))
      allocate(He   (nVar,ids:ide  ,kds:kde+1))
      
      allocate(src  (nVar,ids:ide,kds:kde))
      
      allocate(rho_p(     ids:ide,kds:kde))
      
      allocate(eig_x (ics:ice,kcs:kce))
      allocate(eig_z (ics:ice,kcs:kce))
      
      allocate(relax_coef(ics:ice,kcs:kce))
      
      P     = FillValue
      P_ref = FillValue
      F     = FillValue
      H     = FillValue

      eig_x = 0
      eig_z = 0
      
      ! Set reference pressure
      q_ext = ref%q
      
      do k = kds,kde
        do i = ids,ide
          P_ref(i,k) = calc_pressure(sqrtG(i,k),q_ext(:,i,k))
        enddo
      enddo
      
      !$OMP PARALLEL DO PRIVATE(kp1,km1,kp2,km2,i,ip1,im1,ip2,im2,q_weno,dir)
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
      
      ! Calculate Rayleigh damping coef
      call Rayleigh_coef(relax_coef)
      
    end subroutine init_spatial_operator
    
    subroutine spatial_operator(stat,tend)
      type(stat_field), target, intent(inout) :: stat
      type(tend_field), target, intent(inout) :: tend
      
      real(r_kind), dimension(5) :: q_weno
      
      integer(i_kind) dir
      
      integer(i_kind) i,k,iVar
      
      integer(i_kind) ip1,im1
      integer(i_kind) kp1,km1
      integer(i_kind) ip2,im2
      integer(i_kind) kp2,km2
      
      real(r_kind) dFe
      real(r_kind) dHe
      
      ! Attension stat is changed here!
      call Rayleigh_damping(stat%q,ref%q)
      
      ! copy stat
      q_ext = FillValue
      q_ext(:,ids:ide,kds:kde) = stat%q(:,ids:ide,kds:kde)
      
      if(case_num==2)then
        call fill_constant_inflow(q_ext,ref%q)
      endif
      
      ! Calculate P, F, H and eigenvalues
      !$OMP PARALLEL DO PRIVATE(i)
      do k = kds,kce
        do i = ics,ice
          P    (  i,k) = calc_pressure(sqrtG(i,k),q_ext(:,i,k))
          F    (:,i,k) = calc_F(sqrtG(i,k),q_ext(:,i,k),P(i,k))
          H    (:,i,k) = calc_H(sqrtG(i,k),G13(i,k),q_ext(:,i,k),P(i,k),P_ref(i,k))
          
          eig_x(  i,k) = calc_eigenvalue_x(sqrtG(i,k)         ,q_ext(:,i,k))
          eig_z(  i,k) = calc_eigenvalue_z(sqrtG(i,k),G13(i,k),q_ext(:,i,k))
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! Reconstruction X
      !$OMP PARALLEL DO PRIVATE(i,ip1,im1,ip2,im2,iVar,q_weno,dir)
      do k = kds,kde
        do i = ids-1,ide+1
          ip1 = i + 1
          im1 = i - 1
          ip2 = i + 2
          im2 = i - 2
          do iVar = 1,nVar
            q_weno = F(iVar,im2:ip2,k)
            dir = -1
            call WENO_limiter(FL(iVar,i,k),q_weno,dir)
            dir = 1
            call WENO_limiter(FR(iVar,i,k),q_weno,dir)
            
            q_weno = q_ext(iVar,im2:ip2,k)
            dir = -1
            call WENO_limiter(qL(iVar,i,k),q_weno,dir)
            dir = 1
            call WENO_limiter(qR(iVar,i,k),q_weno,dir)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! Reconstruction Z
      !$OMP PARALLEL DO PRIVATE(kp1,km1,kp2,km2,i,iVar,q_weno,dir)
      do k = kds,kde
        kp1 = k + 1
        km1 = k - 1
        kp2 = k + 2
        km2 = k - 2
        do i = ids,ide
          do iVar = 1,nVar
            q_weno = H(iVar,i,km2:kp2)
            dir = -1
            call WENO_limiter(HB(iVar,i,k),q_weno,dir)
            dir = 1
            call WENO_limiter(HT(iVar,i,k),q_weno,dir)
            
            q_weno = q_ext(iVar,i,km2:kp2)
            dir = -1
            call WENO_limiter(qB(iVar,i,k),q_weno,dir)
            dir = 1
            call WENO_limiter(qT(iVar,i,k),q_weno,dir)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      if(case_num==2)then
        qL(2,ids,kds:kde) = ref%q(2,ids,kds:kde)
        qR(2,ide,kds:kde) = ref%q(2,ide,kds:kde)
        qB(3,ids:ide,kds) = -sqrtGB(:,kds) * G13B(:,kds) * qB(2,:,kds)
        qT(3,ids:ide,kde) = -sqrtGT(:,kde) * G13T(:,kde) * qT(2,:,kde)
      endif
      
      ! initialize source terms
      src = 0.
      
      call bdy_condition(P,stat%q,ref%q,src)
      
      rho_p = ( stat%q(1,ids:ide,kds:kde) + stat%q(5,ids:ide,kds:kde) &
              - ref%q (1,ids:ide,kds:kde) - ref%q (5,ids:ide,kds:kde) ) / sqrtG(ids:ide,kds:kde)
      
      where(abs(rho_p)<=1.E-13)rho_p=0.
      
      src(3,ids:ide,kds:kde) = src(3,ids:ide,kds:kde) - sqrtG(ids:ide,kds:kde) * rho_p(ids:ide,kds:kde) * gravity
      
      ! calc x flux
      !$OMP PARALLEL DO PRIVATE(i,im1,maxeigen_x,iVar)
      do k = kds,kde
        do i = ids,ide+1
        !do i = ids+1,ide
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
        !do k = kds+1,kde
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
      
      !$OMP PARALLEL DO PRIVATE(i,ip1,kp1,iVar,dFe,dHe)
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
        enddo
      enddo
      !$OMP END PARALLEL DO
      
    end subroutine spatial_operator
    
    subroutine bdy_condition(P,q,q_ref,src)
      real(r_kind), dimension(     ics:ice  ,kcs:kce  ), intent(in   ) :: P
      real(r_kind), dimension(nVar,ics:ice  ,kcs:kce  ), intent(in   ) :: q
      real(r_kind), dimension(nVar,ics:ice  ,kcs:kce  ), intent(in   ) :: q_ref
      real(r_kind), dimension(nVar,ids:ide  ,kds:kde  ), intent(inout) :: src
      
      integer(i_kind) dir,sign
      integer(i_kind) i,k,iVar
      
      real(r_kind) :: qrec(3)
      
      real(r_kind), dimension(kds:kde) :: PpL
      real(r_kind), dimension(kds:kde) :: PpR
      real(r_kind), dimension(ids:ide) :: PpB
      real(r_kind), dimension(ids:ide) :: PpT
      
      real(r_kind), dimension(kds:kde) :: sqrtG_P_L
      real(r_kind), dimension(kds:kde) :: sqrtG_P_R
      real(r_kind), dimension(ids:ide) :: sqrtG_G13_P_B
      real(r_kind), dimension(ids:ide) :: sqrtG_G13_P_T
      
      !$OMP PARALLEL DO PRIVATE(qrec)
      do k = kds,kde
        !qrec   = P(ids:ids+2,k) - P_ref(ids:ids+2,k)
        !PpL(k) = left_side_recon3(qrec)
        !qrec   = P(ide-2:ide,k) - P_ref(ide-2:ide,k)
        !PpR(k) = right_side_recon3(qrec)
        
        qrec         = sqrtG(ids:ids+2,k) * P(ids:ids+2,k)
        sqrtG_P_L(k) = left_side_recon3(qrec)
        qrec         = sqrtG(ide-2:ide,k) * P(ide-2:ide,k)
        sqrtG_P_R(k) = right_side_recon3(qrec)
      enddo
      !$OMP END PARALLEL DO
      
      !$OMP PARALLEL DO PRIVATE(qrec)
      do i = ids,ide
        qrec   = P(i,kds:kds+2) - P_ref(i,kds:kds+2)
        where(abs(qrec)/P_ref(i,kds:kds+2)<1.e-13)qrec=0
        PpB(i) = left_side_recon3(qrec)
        
        qrec   = P(i,kde-2:kde) - P_ref(i,kde-2:kde)
        where(abs(qrec)/P_ref(i,kde-2:kde)<1.e-13)qrec=0
        PpT(i) = right_side_recon3(qrec)
        
        qrec             = sqrtG(i,kds:kds+2) * G13(i,kds:kds+2) * P(i,kds:kds+2)
        sqrtG_G13_P_B(i) = left_side_recon3(qrec)
        qrec             = sqrtG(i,kde-2:kde) * G13(i,kde-2:kde) * P(i,kde-2:kde)
        sqrtG_G13_P_T(i) = right_side_recon3(qrec)
      enddo
      !$OMP END PARALLEL DO
      
      if(case_num==1)then
        FL(1,ids,kds:kde) = 0
        FL(2,ids,kds:kde) = sqrtG_P_L!sqrtGL(ids,kds:kde) * PL(ids,kds:kde)
        FL(3,ids,kds:kde) = 0
        FL(4,ids,kds:kde) = 0
        FL(5,ids,kds:kde) = 0
        
        FR(1,ide,kds:kde) = 0
        FR(2,ide,kds:kde) = sqrtG_P_R!sqrtGR(ide,kds:kde) * PR(ide,kds:kde)
        FR(3,ide,kds:kde) = 0
        FR(4,ide,kds:kde) = 0
        FR(5,ide,kds:kde) = 0
        
        HB(1,ids:ide,kds) = 0
        HB(2,ids:ide,kds) = 0
        HB(3,ids:ide,kds) = PpB!PB(ids:ide,kds) - PB_ref(ids:ide,kds)
        HB(4,ids:ide,kds) = 0
        HB(5,ids:ide,kds) = 0
        
        HT(1,ids:ide,kde) = 0
        HT(2,ids:ide,kde) = 0
        HT(3,ids:ide,kde) = PpT!PT(ids:ide,kde) - PT_ref(ids:ide,kde)
        HT(4,ids:ide,kde) = 0
        HT(5,ids:ide,kde) = 0
      
        FR(:,ids-1,kds:kde) = FL(:,ids,kds:kde)
        FL(:,ide+1,kds:kde) = FR(:,ide,kds:kde)
        HT(:,ids:ide,kds-1) = HB(:,ids:ide,kds)
        HB(:,ids:ide,kde+1) = HT(:,ids:ide,kde)
        
        qR(:,ids-1,kds:kde) = qL(:,ids,kds:kde)
        qL(:,ide+1,kds:kde) = qR(:,ide,kds:kde)
        qT(:,ids:ide,kds-1) = qB(:,ids:ide,kds)
        qB(:,ids:ide,kde+1) = qT(:,ids:ide,kde)
      elseif(case_num==2)then
        !FL(1,ids,kds:kde) = q_ref(2,ids,kds:kde)
        !FL(2,ids,kds:kde) = q_ref(2,ids,kds:kde) * q_ref(2,ids,kds:kde) / q_ref(1,ids,kds:kde) + sqrtG_P_L!sqrtGL(ids,kds:kde) * PL(ids,kds:kde)
        !FL(3,ids,kds:kde) = q_ref(3,ids,kds:kde) * q_ref(2,ids,kds:kde) / q_ref(1,ids,kds:kde)
        !FL(4,ids,kds:kde) = q_ref(4,ids,kds:kde) * q_ref(2,ids,kds:kde) / q_ref(1,ids,kds:kde)
        !FL(5,ids,kds:kde) = 0
        
        !FR(1,ide,kds:kde) = q_ref(2,ide,kds:kde)
        !FR(2,ide,kds:kde) = q_ref(2,ide,kds:kde) * q_ref(2,ide,kds:kde) / q_ref(1,ide,kds:kde) + sqrtG_P_R!sqrtGR(ide,kds:kde) * PR(ide,kds:kde)
        !FR(3,ide,kds:kde) = q_ref(3,ide,kds:kde) * q_ref(2,ide,kds:kde) / q_ref(1,ide,kds:kde)
        !FR(4,ide,kds:kde) = q_ref(4,ide,kds:kde) * q_ref(2,ide,kds:kde) / q_ref(1,ide,kds:kde)
        !FR(5,ide,kds:kde) = 0
        
        HB(1,ids:ide,kds) = 0
        HB(2,ids:ide,kds) = sqrtG_G13_P_B!sqrtGB(ids:ide,kds) * G13B(ids:ide,kds) * PB(ids:ide,kds)
        HB(3,ids:ide,kds) = PpB!PB(ids:ide,kds) - PB_ref(ids:ide,kds)
        HB(4,ids:ide,kds) = 0
        HB(5,ids:ide,kds) = 0
        
        HT(1,ids:ide,kde) = 0
        HT(2,ids:ide,kde) = sqrtG_G13_P_T!sqrtGT(ids:ide,kde) * G13T(ids:ide,kde) * PT(ids:ide,kde)
        HT(3,ids:ide,kde) = PpT!PT(ids:ide,kde) - PT_ref(ids:ide,kde)
        HT(4,ids:ide,kde) = 0
        HT(5,ids:ide,kde) = 0
        
        !FR(:,ids-1,kds:kde) = FL(:,ids,kds:kde)
        !FL(:,ide+1,kds:kde) = FR(:,ide,kds:kde)
        HT(:,ids:ide,kds-1) = HB(:,ids:ide,kds)
        HB(:,ids:ide,kde+1) = HT(:,ids:ide,kde)
        
        !qR(:,ids-1,kds:kde) = qL(:,ids,kds:kde)
        !qL(:,ide+1,kds:kde) = qR(:,ide,kds:kde)
        qT(:,ids:ide,kds-1) = qB(:,ids:ide,kds)
        qB(:,ids:ide,kde+1) = qT(:,ids:ide,kde)
      endif
      
    end subroutine bdy_condition
    
    subroutine fill_ghost(q_ext,q,dir,sign)
      real   (r_kind), dimension(ics:ice,kcs:kce), intent(out) :: q_ext
      real   (r_kind), dimension(ics:ice,kcs:kce), intent(in ) :: q
      integer(i_kind)                            , intent(in ) :: dir
      integer(i_kind)                            , intent(in ) :: sign
      
      integer(i_kind) i,k
      
      q_ext(ids:ide,kds:kde) = q(ids:ide,kds:kde)
      
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
    
    subroutine fill_periodic_x(q)
      real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(inout) :: q
      
      integer(i_kind) i,k
      
      do i = 1,extPts
        ! left side
        q(:,ids-i,kds:kde) = q(:,ide-i+1,kds:kde)
        ! right side
        q(:,ide+i,kds:kde) = q(:,ids+i-1,kds:kde)
      enddo
      
    end subroutine fill_periodic_x
    
    subroutine fill_constant_inflow(q,q_ref)
      real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(inout) :: q
      real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(in   ) :: q_ref
      
      integer(i_kind) i,k
      
      ! left side
      q(:,ics:ids-1,kds:kde) = q_ref(:,ics:ids-1,kds:kde)
      ! right side
      q(:,ide+1:ice,kds:kde) = q_ref(:,ide+1:ice,kds:kde)
      !! top
      !q(:,ids:ide,kde+1:kce) = q_ref(:,ids:ide,kde+1:kce)
      
    end subroutine fill_constant_inflow
    
    subroutine Rayleigh_damping(q,q_ref)
      real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(inout) :: q
      real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(in   ) :: q_ref
      
      integer(i_kind), parameter :: vs = 3
      integer(i_kind), parameter :: ve = 3
      
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
      
      real(r_kind), parameter :: topSpongeThickness   = 10000
      real(r_kind), parameter :: leftSpongeThickness  = 10000
      real(r_kind), parameter :: rightSpongeThickness = 10000
      
      real(r_kind), parameter :: mu_max_top = 0.02
      real(r_kind), parameter :: mu_max_lat = 0.02
      
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
    
    function calc_w_eta(sqrtG,G13,q)
      real(r_kind)              :: calc_w_eta
      real(r_kind)              :: sqrtG
      real(r_kind)              :: G13
      real(r_kind),dimension(5) :: q(5)
      
      real(r_kind) :: u
      real(r_kind) :: w

      u = q(2) / ( q(1) + q(5) )
      w = q(3) / ( q(1) + q(5) )
      
      calc_w_eta = w / sqrtG + G13 * u
    
    end function calc_w_eta
    
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
      if(abs(p_pert)/p_ref<1.e-13)p_pert=0
      
      ww = w / sqrtG + G13 * u
      
      calc_H(1) = w1 * ww
      calc_H(2) = w2 * ww + sqrtG * G13 * p
      calc_H(3) = w3 * ww + p_pert
      calc_H(4) = w4 * ww
      calc_H(5) = w5 * ww
      
    end function calc_H
    
    function calc_H_w_eta(sqrtG,G13,q,p,p_ref,w_eta)
      real(r_kind),dimension(5) :: calc_H_w_eta
      real(r_kind)              :: sqrtG
      real(r_kind)              :: G13
      real(r_kind),dimension(5) :: q(5)
      real(r_kind)              :: p      ! pressure
      real(r_kind)              :: p_ref  ! reference pressure
      real(r_kind)              :: w_eta  ! deta/dt
      
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
      if(abs(p_pert)/p_ref<1.e-13)p_pert=0
      
      ww = w_eta
      
      calc_H_w_eta(1) = w1 * ww
      calc_H_w_eta(2) = w2 * ww + sqrtG * G13 * p
      calc_H_w_eta(3) = w3 * ww + p_pert
      calc_H_w_eta(4) = w4 * ww
      calc_H_w_eta(5) = w5 * ww
      
    end function calc_H_w_eta
    
    function calc_eigenvalue_x(sqrtG,q)
      real(r_kind)              :: calc_eigenvalue_x
      real(r_kind)              :: sqrtG
      real(r_kind),dimension(5) :: q(5)
      
      real(r_kind) w1
      real(r_kind) w2
      real(r_kind) w3
      real(r_kind) w4
      real(r_kind) w5
      
      real(r_kind) eig(5)
      
      real(r_kind) coef1,coef2,coef3
      
      if(any(q==FillValue))then
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
      real(r_kind)              :: calc_eigenvalue_z
      real(r_kind)              :: sqrtG
      real(r_kind)              :: G13
      real(r_kind),dimension(5) :: q(5)
      
      real(r_kind) w1
      real(r_kind) w2
      real(r_kind) w3
      real(r_kind) w4
      real(r_kind) w5
      
      real(r_kind) eig(5)
      
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

