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
  real(r_kind), dimension(:,:,:), allocatable :: q_diff
  
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
      
  real(r_kind), dimension(:,:), allocatable :: relax_coef ! Relax coefficient of Rayleigh damping
  
  real(r_kind), dimension(nVar,nVar) :: eigen_mtx_x
  real(r_kind), dimension(nVar,nVar) :: eigen_mtx_z
  
  real(r_kind)                  :: sqrtGe
  real(r_kind)                  :: G13e
  real(r_kind), dimension(nVar) :: qe
  
    contains
    subroutine init_spatial_operator
      integer(i_kind) dir
      integer(i_kind) i,k,iVar
      
      real(r_kind), dimension(5) :: q_weno
      
      integer(i_kind) ip1,im1
      integer(i_kind) kp1,km1
      integer(i_kind) ip2,im2
      integer(i_kind) kp2,km2
      
      allocate(qC    (nVar,ics:ice,kcs:kce))
      allocate(q_diff(nVar,ics:ice,kcs:kce))
                     
      allocate(qL    (nVar,ics:ice,kcs:kce))
      allocate(qR    (nVar,ics:ice,kcs:kce))
      allocate(qB    (nVar,ics:ice,kcs:kce))
      allocate(qT    (nVar,ics:ice,kcs:kce))
                     
      allocate(F     (nVar,ics:ice,kcs:kce))
      allocate(H     (nVar,ics:ice,kcs:kce))
                     
      allocate(FL    (nVar,ics:ice,kcs:kce))
      allocate(FR    (nVar,ics:ice,kcs:kce))
      allocate(HB    (nVar,ics:ice,kcs:kce))
      allocate(HT    (nVar,ics:ice,kcs:kce))
                     
      allocate(P     (     ics:ice,kcs:kce))
      allocate(PL    (     ics:ice,kcs:kce))
      allocate(PR    (     ics:ice,kcs:kce))
      allocate(PB    (     ics:ice,kcs:kce))
      allocate(PT    (     ics:ice,kcs:kce))
      
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
      
      allocate(relax_coef(ics:ice,kcs:kce))
      
      P     = FillValue
      P_ref = FillValue
      F     = FillValue
      H     = FillValue
      
      ! Set reference pressure
      qC = ref%q
      
      do k = kds,kde
        do i = ids,ide
          P_ref(i,k) = calc_pressure(sqrtG(i,k),qC(:,i,k))
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
      
      ! Fill out qC
      qC = FillValue
      
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
      if(case_num==2)call Rayleigh_damping(stat%q,ref%q)
      
      ! copy stat
      qC(:,ids:ide,kds:kde) = stat%q(:,ids:ide,kds:kde)
      
      if(case_num==2)then
        call fill_constant_inflow(qC,ref%q)
      endif
      
      ! Calculate P, F, H and eigenvalues
      !$OMP PARALLEL DO PRIVATE(i)
      do k = kds,kde
        do i = ids,ide
          P    (  i,k) = calc_pressure(sqrtG(i,k),qC(:,i,k))
          F    (:,i,k) = calc_F(sqrtG(i,k),qC(:,i,k),P(i,k))
          H    (:,i,k) = calc_H(sqrtG(i,k),G13(i,k),qC(:,i,k),P(i,k),P_ref(i,k))
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! Reconstruction X
      !$OMP PARALLEL DO PRIVATE(i,ip1,im1,ip2,im2,iVar,q_weno,dir)
      do k = kds,kde
        !do i = ids-1,ide+1
        do i = ids,ide
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
            
            q_weno = qC(iVar,im2:ip2,k)
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
            
            q_weno = qC(iVar,i,km2:kp2)
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
      
      rho_p = ( stat%q(1,ids:ide,kds:kde) - ref%q (1,ids:ide,kds:kde) ) / sqrtG(ids:ide,kds:kde)
      
      where(abs(rho_p)<=1.E-13)rho_p=0.
      
      src(3,ids:ide,kds:kde) = src(3,ids:ide,kds:kde) - sqrtG(ids:ide,kds:kde) * rho_p(ids:ide,kds:kde) * gravity
      
      ! calc x flux
      !$OMP PARALLEL DO PRIVATE(i,im1,eigen_mtx_x,qe,sqrtGe)
      do k = kds,kde
        do i = ids,ide+1
          im1 = i - 1
          
          qe          = 0.5 * ( qL(:,i,k) + qR(:,im1,k) )
          sqrtGe      = 0.5 * ( sqrtGL(i,k) + sqrtGR(im1,k) )
          eigen_mtx_x = calc_eigen_matrix_x( qe, sqrtGe )
          
          Fe(:,i,k) = 0.5 * ( FL(:,i,k) + FR(:,im1,k) - matmul( eigen_mtx_x, ( FL(:,i,k) - FR(:,im1,k) ) ) )
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! calc z flux
      !$OMP PARALLEL DO PRIVATE(k,km1,eigen_mtx_z,qe,sqrtGe,G13e)
      do i = ids,ide
        do k = kds,kde+1
          km1 = k - 1
          
          qe          = 0.5 * ( qB(:,i,k) + qT(:,i,km1) )
          sqrtGe      = 0.5 * ( sqrtGB(i,k) + sqrtGT(i,km1) )
          G13e        = 0.5 * ( G13B(i,k) + G13T(i,km1) )
          eigen_mtx_z = calc_eigen_matrix_z( qe, sqrtGe, G13e )
          
          He(:,i,k) = 0.5 * ( HB(:,i,k) + HT(:,i,km1) - matmul( eigen_mtx_z, ( HB(:,i,k) - HT(:,i,km1) ) ) )
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! Viscosity terms for Density Current case only
      if(case_num==3)then
        do iVar = 2,4
          q_diff(iVar,ids:ide,kds:kde) = qC(iVar,ids:ide,kds:kde) / qC(1,ids:ide,kds:kde)
        enddo
        
        ! Scheme 1, 1st derivative flux
        do iVar = 2,4
          do k = kde,kde
            ! Left bdy
            i = ids
            Fe(iVar,i,k) = Fe(iVar,i,k) - viscosity_coef *  qL(1,i,k) * dqdxL(q_diff(iVar,i:i+2,k),dx)
            i = ids + 1
            Fe(iVar,i,k) = Fe(iVar,i,k) - viscosity_coef * ( qL(1,i,k) + qR(1,i-1,k) ) / 2. * dqdxC(q_diff(iVar,i-1:i,k),dx)
            
            ! Right bdy
            i = ide
            Fe(iVar,i,k) = Fe(iVar,i,k) - viscosity_coef * qR(1,i,k) * dqdxR(q_diff(iVar,i-2:i,k),dx)
            i = ide - 1
            Fe(iVar,i,k) = Fe(iVar,i,k) - viscosity_coef * ( qR(1,i,k) + qL(1,i+1,k) ) / 2. * dqdxC(q_diff(iVar,i-1:i,k),dx)
          enddo
        
          do i = ids,ide
            ! Bottom bdy
            k = kds
            He(iVar,i,k) = He(iVar,i,k) - viscosity_coef * qB(1,i,k) * dqdxL(q_diff(iVar,i,k:k+2),deta) / sqrtGB(i,k)**2
            k = kds + 1
            He(iVar,i,k) = He(iVar,i,k) - viscosity_coef * ( qB(1,i,k) + qT(1,i,k-1) ) / 2. * dqdxC(q_diff(iVar,i,k-1:k),deta) / ( ( sqrtGB(i,k) + sqrtGT(i,k-1) ) / 2. )**2
            
            ! Top bdy
            k = kde
            He(iVar,i,k+1) = He(iVar,i,k+1) - viscosity_coef * qT(1,i,k) * dqdxR(q_diff(iVar,i,k-2:k),deta) / sqrtGT(i,k)**2
            k = kde - 1
            He(iVar,i,k+1) = He(iVar,i,k+1) - viscosity_coef * ( qB(1,i,k+1) + qT(1,i,k) ) / 2. * dqdxC(q_diff(iVar,i,k:k+1),deta) / ( ( sqrtGB(i,k+1) + sqrtGT(i,k) ) / 2. )**2
          enddo
          
          ! Center domain
          do k = kds,kde
            do i = ids+2,ide-2
              Fe(iVar,i,k) = Fe(iVar,i,k) - viscosity_coef * ( qR(1,i,k) + qL(1,i+1,k) ) / 2. * dqdx(q_diff(iVar,i-2:i+1,k),dx)
            enddo
          enddo
          
          do k = kds+1,kde-2
            do i = ids,ide
              He(iVar,i,k+1) = He(iVar,i,k+1) - viscosity_coef * ( qT(1,i,k) + qB(1,i,k+1) ) / 2. * dqdx(q_diff(iVar,i,k-1:k+2),deta) / ( ( sqrtGB(i,k+1) + sqrtGT(i,k) ) / 2. )**2
            enddo
          enddo
        enddo
        
        !! Scheme 2, directionly calculate 2nd derivative
        !!$OMP PARALLEL DO PRIVATE(i,iVar)
        !do k = kds+1,kde-1
        !  do i = ids+1,ide-1
        !    do iVar = 2,4
        !      src(iVar,i,k) = src(iVar,i,k) + viscosity_coef * qC(1,i,k) * ( ( q_diff(iVar,i+1,k) - 2. * q_diff(iVar,i,k) + q_diff(iVar,i-1,k) ) / dx  **2 &
        !                                                                   + ( q_diff(iVar,i,k+1) - 2. * q_diff(iVar,i,k) + q_diff(iVar,i,k-1) ) / deta**2 / sqrtG(i,k)**2 )
        !    enddo
        !  enddo
        !enddo
        !!$OMP END PARALLEL DO
      endif
      
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
        qrec   = P(ids:ids+2,k) - P_ref(ids:ids+2,k)
        PpL(k) = left_side_recon3(qrec)
        qrec   = P(ide-2:ide,k) - P_ref(ide-2:ide,k)
        PpR(k) = right_side_recon3(qrec)
        
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
      
      if(case_num==1.or.case_num==3)then
        FL(1,ids,kds:kde) = 0
        FL(2,ids,kds:kde) = sqrtG_P_L!sqrtGL(ids,kds:kde) * PL(ids,kds:kde)
        FL(3,ids,kds:kde) = 0
        FL(4,ids,kds:kde) = 0
        
        FR(1,ide,kds:kde) = 0
        FR(2,ide,kds:kde) = sqrtG_P_R!sqrtGR(ide,kds:kde) * PR(ide,kds:kde)
        FR(3,ide,kds:kde) = 0
        FR(4,ide,kds:kde) = 0
        
        HB(1,ids:ide,kds) = 0
        HB(2,ids:ide,kds) = 0
        HB(3,ids:ide,kds) = PpB!PB(ids:ide,kds) - PB_ref(ids:ide,kds)
        HB(4,ids:ide,kds) = 0
        
        HT(1,ids:ide,kde) = 0
        HT(2,ids:ide,kde) = 0
        HT(3,ids:ide,kde) = PpT!PT(ids:ide,kde) - PT_ref(ids:ide,kde)
        HT(4,ids:ide,kde) = 0
      
        FR(:,ids-1,kds:kde) = FL(:,ids,kds:kde)
        FL(:,ide+1,kds:kde) = FR(:,ide,kds:kde)
        HT(:,ids:ide,kds-1) = HB(:,ids:ide,kds)
        HB(:,ids:ide,kde+1) = HT(:,ids:ide,kde)
        
        qR(:,ids-1,kds:kde) = qL(:,ids,kds:kde)
        qL(:,ide+1,kds:kde) = qR(:,ide,kds:kde)
        qT(:,ids:ide,kds-1) = qB(:,ids:ide,kds)
        qB(:,ids:ide,kde+1) = qT(:,ids:ide,kde)
      elseif(case_num==2)then
        FL(1,ids,kds:kde) = q_ref(2,ids,kds:kde)
        FL(2,ids,kds:kde) = q_ref(2,ids,kds:kde) * q_ref(2,ids,kds:kde) / q_ref(1,ids,kds:kde) + sqrtG_P_L!sqrtGL(ids,kds:kde) * PL(ids,kds:kde)
        FL(3,ids,kds:kde) = q_ref(3,ids,kds:kde) * q_ref(2,ids,kds:kde) / q_ref(1,ids,kds:kde)
        FL(4,ids,kds:kde) = q_ref(4,ids,kds:kde) * q_ref(2,ids,kds:kde) / q_ref(1,ids,kds:kde)
        
        FR(1,ide,kds:kde) = q_ref(2,ide,kds:kde)
        FR(2,ide,kds:kde) = q_ref(2,ide,kds:kde) * q_ref(2,ide,kds:kde) / q_ref(1,ide,kds:kde) + sqrtG_P_R!sqrtGR(ide,kds:kde) * PR(ide,kds:kde)
        FR(3,ide,kds:kde) = q_ref(3,ide,kds:kde) * q_ref(2,ide,kds:kde) / q_ref(1,ide,kds:kde)
        FR(4,ide,kds:kde) = q_ref(4,ide,kds:kde) * q_ref(2,ide,kds:kde) / q_ref(1,ide,kds:kde)
        
        HB(1,ids:ide,kds) = 0
        HB(2,ids:ide,kds) = sqrtG_G13_P_B!sqrtGB(ids:ide,kds) * G13B(ids:ide,kds) * PB(ids:ide,kds)
        HB(3,ids:ide,kds) = PpB!PB(ids:ide,kds) - PB_ref(ids:ide,kds)
        HB(4,ids:ide,kds) = 0
        
        HT(1,ids:ide,kde) = 0
        HT(2,ids:ide,kde) = sqrtG_G13_P_T!sqrtGT(ids:ide,kde) * G13T(ids:ide,kde) * PT(ids:ide,kde)
        HT(3,ids:ide,kde) = PpT!PT(ids:ide,kde) - PT_ref(ids:ide,kde)
        HT(4,ids:ide,kde) = 0
        
        FR(:,ids-1,kds:kde) = FL(:,ids,kds:kde)
        FL(:,ide+1,kds:kde) = FR(:,ide,kds:kde)
        HT(:,ids:ide,kds-1) = HB(:,ids:ide,kds)
        HB(:,ids:ide,kde+1) = HT(:,ids:ide,kde)
        
        qR(:,ids-1,kds:kde) = qL(:,ids,kds:kde)
        qL(:,ide+1,kds:kde) = qR(:,ide,kds:kde)
        qT(:,ids:ide,kds-1) = qB(:,ids:ide,kds)
        qB(:,ids:ide,kde+1) = qT(:,ids:ide,kde)
      endif
      
    end subroutine bdy_condition
    
    subroutine fill_ghost(qC,q,dir,sign)
      real   (r_kind), dimension(ics:ice,kcs:kce), intent(out) :: qC
      real   (r_kind), dimension(ics:ice,kcs:kce), intent(in ) :: q
      integer(i_kind)                            , intent(in ) :: dir
      integer(i_kind)                            , intent(in ) :: sign
      
      integer(i_kind) i,k
      
      qC(ids:ide,kds:kde) = q(ids:ide,kds:kde)
      
      if(dir == 1)then
        ! x-dir
        do i = 1,extPts
          qC(ids-i,kds:kde) = sign * q(ids+i-1,kds:kde)
          qC(ide+i,kds:kde) = sign * q(ide-i+1,kds:kde)
        enddo
      elseif(dir == 2)then
        ! z-dir
        do k = 1,extPts
          qC(ids:ide,kds-k) = sign * q(ids:ide,kds+k-1)
          qC(ids:ide,kde+k) = sign * q(ids:ide,kde-k+1)
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
      
      integer(i_kind), parameter :: vs = 1
      integer(i_kind), parameter :: ve = 4
      
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
      
      real(r_kind), parameter :: mu_max_top = 0.15
      real(r_kind), parameter :: mu_max_lat = 0.15
      
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
      real(r_kind)                 :: calc_pressure
      real(r_kind)                 :: sqrtG
      real(r_kind),dimension(nVar) :: q
      
      real(r_kind) w1
      real(r_kind) w4
      
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
      
      rho      = w1 / sqrtG
      gamma    = 0
      sh       = gamma / ( 1. + gamma )
      R        = ( 1. + eq * sh ) * Rd
      theta    = w4 / w1
      cp       = cpd + ( cpv - cpd ) * sh
      cv       = cvd + ( cvv - cvd ) * sh
      kappa    = cp / cv
      
      calc_pressure = p0 * ( rho * theta * R / p0 )**kappa
      
    end function calc_pressure
    
    function calc_w_eta(sqrtG,G13,q)
      real(r_kind)                 :: calc_w_eta
      real(r_kind)                 :: sqrtG
      real(r_kind)                 :: G13
      real(r_kind),dimension(nVar) :: q
      
      real(r_kind) :: u
      real(r_kind) :: w

      u = q(2) / q(1)
      w = q(3) / q(1)
      
      calc_w_eta = w / sqrtG + G13 * u
    
    end function calc_w_eta
    
    function calc_F(sqrtG,q,P)
      real(r_kind),dimension(nVar) :: calc_F
      real(r_kind)                 :: sqrtG
      real(r_kind),dimension(nVar) :: q
      real(r_kind)                 :: p      ! pressure
      
      real(r_kind) w1
      real(r_kind) w2
      real(r_kind) w3
      real(r_kind) w4
      
      real(r_kind) sqrtGrho
      real(r_kind) u
      
      w1 = q(1)
      w2 = q(2)
      w3 = q(3)
      w4 = q(4)
      
      sqrtGrho = w1
      u        = w2 / sqrtGrho
      
      calc_F(1) = w1 * u
      calc_F(2) = w2 * u + sqrtG * p
      calc_F(3) = w3 * u
      calc_F(4) = w4 * u
      
    end function calc_F
    
    function calc_H(sqrtG,G13,q,p,p_ref)
      real(r_kind),dimension(nVar) :: calc_H
      real(r_kind)                 :: sqrtG
      real(r_kind)                 :: G13
      real(r_kind),dimension(nVar) :: q
      real(r_kind)                 :: p      ! pressure
      real(r_kind)                 :: p_ref  ! reference pressure
      
      real(r_kind)                 :: p_pert ! pressure  perturbation
      
      real(r_kind) w1
      real(r_kind) w2
      real(r_kind) w3
      real(r_kind) w4
      real(r_kind) ww
      
      real(r_kind) sqrtGrho
      real(r_kind) u
      real(r_kind) w
      
      w1 = q(1)
      w2 = q(2)
      w3 = q(3)
      w4 = q(4)
      
      sqrtGrho = w1
      u        = w2 / sqrtGrho
      w        = w3 / sqrtGrho
      p_pert   = p - p_ref
      if(abs(p_pert)/p_ref<1.e-13)p_pert=0
      
      ww = w / sqrtG + G13 * u
      
      calc_H(1) = w1 * ww
      calc_H(2) = w2 * ww + sqrtG * G13 * p
      calc_H(3) = w3 * ww + p_pert
      calc_H(4) = w4 * ww
      
    end function calc_H
    
    function calc_H_w_eta(sqrtG,G13,q,p,p_ref,w_eta)
      real(r_kind),dimension(nVar) :: calc_H_w_eta
      real(r_kind)                 :: sqrtG
      real(r_kind)                 :: G13
      real(r_kind),dimension(nVar) :: q
      real(r_kind)                 :: p      ! pressure
      real(r_kind)                 :: p_ref  ! reference pressure
      real(r_kind)                 :: w_eta  ! deta/dt
      
      real(r_kind)                 :: p_pert ! pressure  perturbation
      
      real(r_kind) w1
      real(r_kind) w2
      real(r_kind) w3
      real(r_kind) w4
      real(r_kind) ww
      
      real(r_kind) sqrtGrho
      real(r_kind) u
      real(r_kind) w
      
      w1 = q(1)
      w2 = q(2)
      w3 = q(3)
      w4 = q(4)
      
      sqrtGrho = w1
      u        = w2 / sqrtGrho
      w        = w3 / sqrtGrho
      p_pert   = p - p_ref
      if(abs(p_pert)/p_ref<1.e-13)p_pert=0
      
      ww = w_eta
      
      calc_H_w_eta(1) = w1 * ww
      calc_H_w_eta(2) = w2 * ww + sqrtG * G13 * p
      calc_H_w_eta(3) = w3 * ww + p_pert
      calc_H_w_eta(4) = w4 * ww
      
    end function calc_H_w_eta
    
    function calc_eigenvalue_x(sqrtG,q)
      real(r_kind),dimension(nVar) :: calc_eigenvalue_x
      real(r_kind)                 :: sqrtG
      real(r_kind),dimension(nVar) :: q
      
      real(r_kind) w1
      real(r_kind) w2
      real(r_kind) w3
      real(r_kind) w4
      
      real(r_kind) eig(nVar)
      
      real(r_kind) coef1,coef2,coef3
      
      if(any(q==FillValue))then
        eig = 0
      else
        w1 = q(1)
        w2 = q(2)
        w3 = q(3)
        w4 = q(4)
        
        coef1 = w2
        
        coef2 = sqrt( cpd * p0 * sqrtG * w1 / cvd ) * ((Rd*w4)/(p0*sqrtG))**(cpd/(2*cvd))
        
        coef3 = w1
        
        eig(1) = w2 / w1
        eig(2) = eig(1)
        eig(3) = ( coef1 - coef2 ) / coef3
        eig(4) = ( coef1 + coef2 ) / coef3
      endif
      
      calc_eigenvalue_x = eig
      
    end function calc_eigenvalue_x
    
    function calc_eigenvalue_z(sqrtG,G13,q)
      real(r_kind),dimension(nVar) :: calc_eigenvalue_z
      real(r_kind)                 :: sqrtG
      real(r_kind)                 :: G13
      real(r_kind),dimension(nVar) :: q
      
      real(r_kind) w1
      real(r_kind) w2
      real(r_kind) w3
      real(r_kind) w4
      
      real(r_kind) eig(nVar)
      
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
        
        coef1 = cvd*sqrtG**2*w1**3*(G13*sqrtG*w2 + w3)*w4
        
        coef2 = sqrt( cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7*w4**2&
                    *((Rd*w4)/(p0*sqrtG))**(cpd/cvd) )
        
        coef3 = cvd*sqrtG**3*w1**4*w4
        
        drhoetadt = (G13*sqrtG*w2 + w3)/ (sqrtG*w1)
        
        eig(1) = drhoetadt
        eig(2) = drhoetadt
        eig(3) = ( coef1 - coef2 ) / coef3
        eig(4) = ( coef1 + coef2 ) / coef3
        
        !eig(3) = (sqrtG**2*w1**3*(G13*sqrtG*w2 + w3) -                           &
        !            Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7*w4**2*  &
        !                  ((Rd*w4)/(p0*sqrtG))**(cpd/cvd))/(cvd*w4))/            &
        !         (sqrtG**3*w1**4)
        !eig(4) = (cvd*sqrtG**2*w1**3*(G13*sqrtG*w2 + w3)*w4 +                    &
        !             Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7*w4**2* &
        !                 ((Rd*w4)/(p0*sqrtG))**(cpd/cvd)))/                      &
        !          (cvd*sqrtG**3*w1**4*w4)
      endif
      
      calc_eigenvalue_z = eig
      
    end function calc_eigenvalue_z
    
    ! A_x = R_x \lambda_x L_x, L_x = R_x^-1
    function calc_eigen_matrix_x(q,sqrtG)
      real(r_kind), dimension(nVar,nVar)              :: calc_eigen_matrix_x
      real(r_kind), dimension(nVar     ), intent(in ) :: q
      real(r_kind),                       intent(in ) :: sqrtG
      
      real(r_kind), dimension(nVar,nVar)              :: mtx
      
      real(r_kind) w1
      real(r_kind) w2
      real(r_kind) w3
      real(r_kind) w4
      
      real(r_kind) a,b,c
      
      real(r_kind), dimension(nVar) :: eigen_value
      
      w1 = q(1)
      w2 = q(2)
      w3 = q(3)
      w4 = q(4)
      
      eigen_value = calc_eigenvalue_x(sqrtG,q)
      a = sign( 1._r_kind, eigen_value(1) )
      b = sign( 1._r_kind, eigen_value(3) )
      c = sign( 1._r_kind, eigen_value(4) )
      
      mtx(1,1) = (b*Sqrt(cvd)*w2 - c*Sqrt(cvd)*w2 +    &
                  2*a*Sqrt(cpd*p0*sqrtG*w1)            &
                   *((Rd*w4)/(p0*sqrtG))**             &
                     (cpd/(2*cvd)))/                   &
                 (((Rd*w4)/(p0*sqrtG))**(cpd/(2*cvd))* &
                 (2*Sqrt(cpd*p0*sqrtG*w1)))
      
      mtx(1,2) = ((-b + c)*Sqrt(cvd*w1))/               &
                 (((Rd*w4)/(p0*sqrtG))**(cpd/(2*cvd))*  &
                  (2*Sqrt(cpd*p0*sqrtG)))
      
      mtx(1,3) = 0
      
      mtx(1,4) = ((-2*a + b + c)*w1)/(2*w4)
      
      mtx(2,1) = (w2*(2*a*Sqrt(w1) - b*Sqrt(w1) -        &
               c*Sqrt(w1) + (b*Sqrt(cvd)*w2)/            &
                 (((Rd*w4)/(p0*sqrtG))**(cpd/(2*cvd))*   &
                    (Sqrt(cpd*p0*sqrtG))) -              &
               (c*Sqrt(cvd)*w2)/(((Rd*w4)/(p0*sqrtG))**  &
               (cpd/(2*cvd))*(Sqrt(cpd*p0*sqrtG)))))/(2*w1**(3./2.))
      
      mtx(2,2) = 0.5 * (b + c - (b*Sqrt(cvd)*w2)/   &
             (((Rd*w4)/(p0*sqrtG))**(cpd/(2*cvd))*  &
                (Sqrt(cpd*p0*sqrtG*w1)))            &
                + (c*Sqrt(cvd)*w2)/                 &
             (((Rd*w4)/(p0*sqrtG))**(cpd/(2*cvd))*  &
                (Sqrt(cpd*p0*sqrtG*w1))))
      
      mtx(2,3) = 0
      
      mtx(2,4) = -((2*a*w2 - b*w2 - c*w2 +           &
             (b*Sqrt(cpd*p0*sqrtG*w1)                &
              *((Rd*w4)/(p0*sqrtG))**                &
                    (cpd/(2*cvd)))/Sqrt(cvd) -       &
             (c*Sqrt(cpd*p0*sqrtG*w1)                &
             *((Rd*w4)/(p0*sqrtG))**                 &
                    (cpd/(2*cvd)))/Sqrt(cvd))/(2*w4))
      
      mtx(3,1) = ((b - c)*Sqrt(cvd)*w2*w3)/          &
              (((Rd*w4)/(p0*sqrtG))**(cpd/(2*cvd))*  &
                 (2*Sqrt(cpd*p0*sqrtG)*              &
                    w1**(3./2.)))
      
      mtx(3,2) = -(((b - c)*Sqrt(cvd)*w3)/           &
             (((Rd*w4)/(p0*sqrtG))**(cpd/(2*cvd))*   &
                (2*Sqrt(cpd*p0*sqrtG*w1))))
      
      mtx(3,3) = a
      
      mtx(3,4) = ((-2*a + b + c)*w3)/(2*w4)
      
      mtx(4,1) = ((b - c)*Sqrt(cvd)*w2*w4)/           &
                (((Rd*w4)/(p0*sqrtG))**(cpd/(2*cvd))* &
                   (2*Sqrt(cpd*p0*sqrtG)* &
                      w1**(3./2.)))
      
      mtx(4,2) = -(((b - c)*Sqrt(cvd)*w4)/  &
      (((Rd*w4)/(p0*sqrtG))**(cpd/(2*cvd))* &
         (2*Sqrt(cpd*p0*sqrtG*w1))))
            
      mtx(4,3) = 0
      
      mtx(4,4) = ( b + c ) / 2.
      
      calc_eigen_matrix_x = mtx
    end function calc_eigen_matrix_x
    
    ! A_z = R_z \lambda_z L_z, L_z = R_z^-1
    function calc_eigen_matrix_z(q,sqrtG,G13)
      real(r_kind), dimension(nVar,nVar)              :: calc_eigen_matrix_z
      real(r_kind), dimension(nVar     ), intent(in ) :: q
      real(r_kind),                       intent(in ) :: sqrtG
      real(r_kind),                       intent(in ) :: G13
      
      real(r_kind), dimension(nVar,nVar)              :: mtx
      
      real(r_kind) w1
      real(r_kind) w2
      real(r_kind) w3
      real(r_kind) w4
      
      real(r_kind) a,b,c
      
      real(r_kind), dimension(nVar) :: eigen_value
      
      w1 = q(1)
      w2 = q(2)
      w3 = q(3)
      w4 = q(4)
      
      eigen_value = calc_eigenvalue_z(sqrtG,G13,q)
      a = sign( 1._r_kind, eigen_value(1) )
      b = sign( 1._r_kind, eigen_value(3) )
      c = sign( 1._r_kind, eigen_value(4) )
      
      mtx(1,1) = (b*cvd*sqrtG**2*w1**3*(G13*sqrtG*w2 + w3)*w4 -               &
                    c*cvd*sqrtG**2*w1**3*(G13*sqrtG*w2 + w3)*w4 +             &
                    2*a*Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7* &
                          w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)))/            &
                 (2*Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7*     &
                        w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)))
          
      mtx(1,2) = -(((b - c)*cvd*G13*sqrtG**3*w1**4*w4)/                      &
                 (2*Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7*    &
                        w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd))))
      
      mtx(1,3) = -(((b - c)*cvd*sqrtG**2*w1**4*w4)/                      &
                (2*Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7* &
                       w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd))))
      
      mtx(1,4) = ((-2*a + b + c)*w1)/(2*w4)
      
      mtx(2,1) = -((sqrtG*(G13*sqrtG*w2 + w3)*                                 &
                   ((-b)*cvd*sqrtG*(1 + G13**2*sqrtG**2)*w1**3*w2*w4 +         &
                      c*cvd*sqrtG*(1 + G13**2*sqrtG**2)*w1**3*w2*w4 -          &
                      2*a*G13*Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*  &
                            w1**7*w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)) +     &
                      b*G13*Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*    &
                            w1**7*w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)) +     &
                      c*G13*Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*    &
                            w1**7*w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd))))/    &
                (2*(w1 + G13**2*sqrtG**2*w1)*                                  &
                   Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7*w4**2* &
                       ((Rd*w4)/(p0*sqrtG))**(cpd/cvd))))
      
      mtx(2,2) = (2*a*Sqrt(                                                      &
                  cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7*               &
                   w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)) -                      &
                   b*G13*sqrtG**2*(cvd*sqrtG*(1 + G13**2*sqrtG**2)*w1**3*w2*w4 - &
                   G13*Sqrt(                                                     &
                     cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7*            &
                      w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd))) +                  &
                   c*G13*sqrtG**2*(cvd*sqrtG*(1 + G13**2*sqrtG**2)*w1**3*w2*w4 + &
                   G13*Sqrt(                                                     &
                     cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7*            &
                      w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd))))/                  &
                (2*(1 + G13**2*sqrtG**2)*                                        &
                Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7*            &
                  w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)))
      
      mtx(2,3) = 0.5*sqrtG*(-((2*a*G13)/(1 + G13**2*sqrtG**2)) + (b*G13)/(1 + &
                  G13**2*sqrtG**2) + (c*G13)/(1 + G13**2*sqrtG**2) -          &
                  (b*cvd*sqrtG*w1**3*w2*w4)/                                  &
                Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7*         &
                  w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)) +                    &
                  (c*cvd*sqrtG*w1**3*w2*w4)/                                  &
                Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7*         &
                  w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)))
      
      mtx(2,4) = (-2*a*cvd*sqrtG*(1 + G13**2*sqrtG**2)*w1**3*w2*w4 +     &
                  b*cvd*sqrtG*(1 + G13**2*sqrtG**2)*w1**3*w2*w4 +        &
                     c*cvd*sqrtG*(1 + G13**2*sqrtG**2)*w1**3*w2*w4 -     &
                     b*G13*                                              &
                   Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7* &
                     w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)) +            &
                     c*G13*                                              &
                   Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7* &
                     w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)))/            &
                  (2*cvd*sqrtG*(1 + G13**2*sqrtG**2)*w1**3*w4**2)
      
      mtx(3,1) = -(((G13*sqrtG*w2 + w3)*                                         &
                   (-2*a*Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*         &
                            w1**7*w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)) +       &
                      b*((-cvd)*sqrtG**2*(1 + G13**2*sqrtG**2)*w1**3*w3*w4 +     &
                           Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7* &
                               w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd))) +         &
                      c*(cvd*sqrtG**2*(1 + G13**2*sqrtG**2)*w1**3*w3*w4 +        &
                           Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7* &
                               w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)))))/        &
                (2*(w1 + G13**2*sqrtG**2*w1)*                                    &
                   Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7*w4**2*   &
                       ((Rd*w4)/(p0*sqrtG))**(cpd/cvd))))
      
      mtx(3,2) = (G13*sqrtG*(-2*a*Sqrt(cpd*cvd*p0*sqrtG**5*                         &
                           (1 + G13**2*sqrtG**2)*w1**7*w4**2*((Rd*w4)/(p0*sqrtG))** &
                             (cpd/cvd)) +                                           &
                     b*((-cvd)*sqrtG**2*(1 + G13**2*sqrtG**2)*w1**3*w3*w4 +         &
                          Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7*     &
                              w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd))) +             &
                     c*(cvd*sqrtG**2*(1 + G13**2*sqrtG**2)*w1**3*w3*w4 +            &
                          Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7*     &
                              w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)))))/            &
               (2*(1 + G13**2*sqrtG**2)*Sqrt(cpd*cvd*p0*sqrtG**5*                   &
                      (1 + G13**2*sqrtG**2)*w1**7*w4**2*((Rd*w4)/(p0*sqrtG))**      &
                        (cpd/cvd)))
      
      mtx(3,3) = (2*a*G13**2*sqrtG**2*Sqrt(cpd*cvd*p0*sqrtG**5*                  &
                     (1 + G13**2*sqrtG**2)*w1**7*w4**2*((Rd*w4)/(p0*sqrtG))**    &
                       (cpd/cvd)) + b*((-cvd)*sqrtG**2*(1 + G13**2*sqrtG**2)*    &
                      w1**3*w3*w4 + Sqrt(cpd*cvd*p0*sqrtG**5*                    &
                        (1 + G13**2*sqrtG**2)*w1**7*w4**2*((Rd*w4)/(p0*sqrtG))** &
                          (cpd/cvd))) + c*(cvd*sqrtG**2*(1 + G13**2*sqrtG**2)*   &
                      w1**3*w3*w4 + Sqrt(cpd*cvd*p0*sqrtG**5*                    &
                        (1 + G13**2*sqrtG**2)*w1**7*w4**2*((Rd*w4)/(p0*sqrtG))** &
                          (cpd/cvd))))/(2*(1 + G13**2*sqrtG**2)*                 &
               Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7*w4**2*       &
                   ((Rd*w4)/(p0*sqrtG))**(cpd/cvd)))
      
      mtx(3,4) = (-2*a*w3*w4 + b*w3*w4 + c*w3*w4 -                         &
                 (b*Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7*  &
                          w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)))/         &
                   (sqrtG**2*(cvd + cvd*G13**2*sqrtG**2)*w1**3) +          &
                 (c*Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7*  &
                          w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)))/         &
                   (sqrtG**2*(cvd + cvd*G13**2*sqrtG**2)*w1**3))/(2*w4**2)
      
      mtx(4,1) = ((b - c)*cvd*sqrtG**2*w1**2*(G13*sqrtG*w2 + w3)*w4**2)/   &
                 (2*Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7*  &
                        w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd)))
      
      mtx(4,2) = -(((b - c)*cvd*G13*sqrtG**3*w1**3*w4**2)/                 &
                 (2*Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7*  &
                        w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd))))
            
      mtx(4,3) = -(((b - c)*cvd*sqrtG**2*w1**3*w4**2)/                   &
                (2*Sqrt(cpd*cvd*p0*sqrtG**5*(1 + G13**2*sqrtG**2)*w1**7* &
                       w4**2*((Rd*w4)/(p0*sqrtG))**(cpd/cvd))))
      
      mtx(4,4) = ( b + c ) / 2.
      
      calc_eigen_matrix_z = mtx
    end function calc_eigen_matrix_z
    
END MODULE spatial_operators_mod

