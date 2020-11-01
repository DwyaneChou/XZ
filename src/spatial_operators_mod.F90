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
  real(r_kind), dimension(:,:,:), allocatable :: q_p   ! q perturbation
  
  real(r_kind), dimension(:,:,:), allocatable :: F
  real(r_kind), dimension(:,:,:), allocatable :: H
  
  real(r_kind), dimension(  :,:), allocatable :: P
  real(r_kind), dimension(  :,:), allocatable :: P_ref ! reference pressure
  
  real(r_kind), dimension(:,:,:), allocatable :: Fn     ! Reconstructed F^-
  real(r_kind), dimension(:,:,:), allocatable :: FnL    ! Reconstructed F^-_(i-1/2,k)
  real(r_kind), dimension(:,:,:), allocatable :: FnR    ! Reconstructed F^-_(i+1/2,k)
  real(r_kind), dimension(:,:,:), allocatable :: Fp     ! Reconstructed F^+
  real(r_kind), dimension(:,:,:), allocatable :: FpL    ! Reconstructed F^+_(i-1/2,k)
  real(r_kind), dimension(:,:,:), allocatable :: FpR    ! Reconstructed F^+_(i+1/2,k)
      
  real(r_kind), dimension(:,:,:), allocatable :: Hn     ! Reconstructed H^-
  real(r_kind), dimension(:,:,:), allocatable :: HnB    ! Reconstructed H^-_(i,k-1/2)
  real(r_kind), dimension(:,:,:), allocatable :: HnT    ! Reconstructed H^-_(i,k+1/2)
  real(r_kind), dimension(:,:,:), allocatable :: Hp     ! Reconstructed H^+
  real(r_kind), dimension(:,:,:), allocatable :: HpB    ! Reconstructed H^+_(i,k-1/2)
  real(r_kind), dimension(:,:,:), allocatable :: HpT    ! Reconstructed H^+_(i,k+1/2)
  
  real(r_kind), dimension(:,:,:), allocatable :: Fe    ! F on edges of each cell
  real(r_kind), dimension(:,:,:), allocatable :: He    ! H on edges of each cell
  
  real(r_kind), dimension(:,:,:), allocatable :: src   ! source term
  
  real(r_kind), dimension(  :,:), allocatable :: rho_p ! density perturbation
      
  real(r_kind), dimension(  :,:), allocatable :: eig_x
  real(r_kind), dimension(  :,:), allocatable :: eig_z
  
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
      
      allocate(q_ext(nVar,ics:ice,kcs:kce))
      allocate(q_p  (nVar,ics:ice,kcs:kce))
      
      allocate(F    (nVar,ics:ice,kcs:kce))
      allocate(H    (nVar,ics:ice,kcs:kce))
      
      allocate(P    (     ics:ice,kcs:kce))
      allocate(P_ref(     ics:ice,kcs:kce))
      
      allocate(Fn   (nVar,ics:ice,kcs:kce))
      allocate(FnL  (nVar,ics:ice,kcs:kce))
      allocate(FnR  (nVar,ics:ice,kcs:kce))
      allocate(Fp   (nVar,ics:ice,kcs:kce))
      allocate(FpL  (nVar,ics:ice,kcs:kce))
      allocate(FpR  (nVar,ics:ice,kcs:kce))
      
      allocate(Hn   (nVar,ics:ice,kcs:kce))
      allocate(HnB  (nVar,ics:ice,kcs:kce))
      allocate(HnT  (nVar,ics:ice,kcs:kce))
      allocate(Hp   (nVar,ics:ice,kcs:kce))
      allocate(HpB  (nVar,ics:ice,kcs:kce))
      allocate(HpT  (nVar,ics:ice,kcs:kce))
      
      allocate(Fe   (nVar,ids:ide+1,kds:kde  ))
      allocate(He   (nVar,ids:ide  ,kds:kde+1))
      
      allocate(src  (nVar,ids:ide,kds:kde))
      
      allocate(rho_p(     ids:ide,kds:kde))
      
      allocate(eig_x (ics:ice,kcs:kce))
      allocate(eig_z (ics:ice,kcs:kce))
      
      ! Set reference pressure
      q_ext = ref%q
      call bdy_condition(q_ext,q_ext,ref%q,src)
      
      P     = FillValue
      P_ref = FillValue
      F     = FillValue
      H     = FillValue
      
      do k = kds,kde
        do i = ids,ide
          P_ref(i,k) = calc_pressure(sqrtG(i,k),q_ext(:,i,k))
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
      
      real(r_kind) dFe
      real(r_kind) dHe
      
      ! copy stat
      q_ext = stat%q
      
      q_p(1,:,:) = stat%q(1,:,:) - ref%q(1,:,:)
      q_p(2,:,:) = stat%q(2,:,:)
      q_p(3,:,:) = stat%q(3,:,:)
      q_p(4,:,:) = stat%q(4,:,:) - ref%q(4,:,:)
      q_p(5,:,:) = stat%q(5,:,:)
      
      ! initialize source terms
      src = 0.
      
      ! Set no flux and nonreflecting boundary condition
      call bdy_condition(q_ext,stat%q,ref%q,src)
      
      do k = kcs,kce
        do i = ics,ice
          P    (  i,k) = calc_pressure(sqrtG(i,k),q_ext(:,i,k))
          F    (:,i,k) = calc_F(sqrtG(i,k),q_ext(:,i,k),P(i,k))
          H    (:,i,k) = calc_H(sqrtG(i,k),G13(i,k),q_ext(:,i,k),P(i,k),P_ref(i,k))
          
          eig_x(  i,k) = calc_eigenvalue_x(sqrtG(i,k)         ,q_ext(:,i,k))
          eig_z(  i,k) = calc_eigenvalue_z(sqrtG(i,k),G13(i,k),q_ext(:,i,k))
        enddo
      enddo
      
      ! calc Fn and Fp
      do i = ids-1,ide+1
        ip2 = i + 2
        im2 = i - 2
        do k = kds,kde
          Fn(:,i,k) = F(:,i,k) - maxval(eig_x(im2:ip2,k)) * q_ext(:,i,k)
          Fp(:,i,k) = F(:,i,k) + maxval(eig_x(im2:ip2,k)) * q_ext(:,i,k)
        enddo
      enddo
      
      ! calc Hn and Hp
      do k = kds-1,kde+1
        kp2 = k + 2
        km2 = k - 2
        do i = ids,ide
          Hn(:,i,k) = H(:,i,k) - maxval(eig_z(i,km2:kp2)) * q_ext(:,i,k)
          Hp(:,i,k) = H(:,i,k) + maxval(eig_z(i,km2:kp2)) * q_ext(:,i,k)
        enddo
      enddo
      
      ! Reconstruct FnL and FpR
      !$OMP PARALLEL DO PRIVATE(i,iVar,ip1,im1,ip2,im2,q_weno,dir,kp1,km1,kp2,km2)
      do k = kds,kde
        kp1 = k + 1
        km1 = k - 1
        kp2 = k + 2
        km2 = k - 2
        do i = ids-1,ide+1
          ip1 = i + 1
          im1 = i - 1
          ip2 = i + 2
          im2 = i - 2
          do iVar = 1,nVar
            ! F^-
            q_weno = Fn(iVar,im2:ip2,k)
            dir = -1
            call WENO_limiter(FnL(iVar,i,k),q_weno,dir)
            ! F^+
            q_weno = Fp(iVar,im2:ip2,k)
            dir = 1
            call WENO_limiter(FpR(iVar,i,k),q_weno,dir)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! Reconstruct HnB and HpT
      !$OMP PARALLEL DO PRIVATE(i,iVar,ip1,im1,ip2,im2,q_weno,dir,kp1,km1,kp2,km2)
      do k = kds-1,kde+1
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
            ! H^-
            q_weno = Hn(iVar,i,km2:kp2)
            dir = -1
            call WENO_limiter(HnB(iVar,i,k),q_weno,dir)
            ! H^+
            q_weno = Hp(iVar,i,km2:kp2)
            dir = 1
            call WENO_limiter(HpT(iVar,i,k),q_weno,dir)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      rho_p = ( stat%q(1,ids:ide,kds:kde) + stat%q(5,ids:ide,kds:kde) &
              - ref%q (1,ids:ide,kds:kde) - ref%q (5,ids:ide,kds:kde) ) / sqrtG(ids:ide,kds:kde)
      
      where(abs(rho_p)<=1.E-13)rho_p=0.
      
      src(3,ids:ide,kds:kde) = src(3,ids:ide,kds:kde) - sqrtG(ids:ide,kds:kde) * rho_p(ids:ide,kds:kde) * gravity
      
      ! calc x flux
      !$OMP PARALLEL DO PRIVATE(i,im1)
      do k = kds,kde
        do i = ids,ide+1
          im1 = i - 1
          Fe(:,i,k) = 0.5 * ( FnL(:,i,k) + FpR(:,im1,k) )
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! calc z flux
      !$OMP PARALLEL DO PRIVATE(k,km1)
      do i = ids,ide
        do k = kds,kde+1
          km1 = k - 1
          He(:,i,k) = 0.5 * ( HnB(:,i,k) + HpT(:,i,km1) )
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
      
      integer(i_kind), parameter :: vs = 2
      integer(i_kind), parameter :: ve = 3
      integer(i_kind), parameter :: bdy_width = 30
      real   (r_kind), parameter :: exp_ceof  = 2
      
      integer(i_kind) dir,sign
      integer(i_kind) i,k,iVar
      
      integer(i_kind) il,ir
      integer(i_kind) kt
      
      integer(i_kind) kls,kle ! pure lateral boundary layer indices
      integer(i_kind) its,ite ! pure top boundary layer indices
      
      real(r_kind) :: relax_coef(bdy_width)
      real(r_kind) :: max_exp
      
      if(case_num==1)then
        ! x-dir
        call fill_ghost(q_ext(1,:,:),q_p(1,:,:),dir=1,sign= 1)
        call fill_ghost(q_ext(2,:,:),q_p(2,:,:),dir=1,sign=-1)
        call fill_ghost(q_ext(3,:,:),q_p(3,:,:),dir=1,sign= 1)
        call fill_ghost(q_ext(4,:,:),q_p(4,:,:),dir=1,sign= 1)
        call fill_ghost(q_ext(5,:,:),q_p(5,:,:),dir=1,sign= 1)
        ! z-dir
        call fill_ghost(q_ext(1,:,:),q_p(1,:,:),dir=2,sign= 1)
        call fill_ghost(q_ext(2,:,:),q_p(2,:,:),dir=2,sign= 1)
        call fill_ghost(q_ext(3,:,:),q_p(3,:,:),dir=2,sign=-1)
        call fill_ghost(q_ext(4,:,:),q_p(4,:,:),dir=2,sign= 1)
        call fill_ghost(q_ext(5,:,:),q_p(5,:,:),dir=2,sign= 1)
        
        q_ext(1,:,:) = q_ext(1,:,:) + q_ref(1,:,:)
        q_ext(4,:,:) = q_ext(4,:,:) + q_ref(4,:,:)
      elseif(case_num==2)then
        ! z-dir
        call fill_ghost(q_ext(1,:,:),q_p(1,:,:),dir=2,sign= 1)
        call fill_ghost(q_ext(2,:,:),q_p(2,:,:),dir=2,sign= 1)
        call fill_ghost(q_ext(3,:,:),q_p(3,:,:),dir=2,sign=-1)
        call fill_ghost(q_ext(4,:,:),q_p(4,:,:),dir=2,sign= 1)
        call fill_ghost(q_ext(5,:,:),q_p(5,:,:),dir=2,sign= 1)
        
        q_ext(1,:,:) = q_ext(1,:,:) + q_ref(1,:,:)
        q_ext(4,:,:) = q_ext(4,:,:) + q_ref(4,:,:)
        
        ! left
        q_ext(:,ics:ids-1,:) = q_ref(:,ics:ids-1,:)
        
        ! right
        q_ext(:,ide+1:ice,:) = q_ref(:,ide+1:ice,:)
        
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

