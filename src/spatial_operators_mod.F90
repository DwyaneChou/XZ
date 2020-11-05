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
      
      real(r_kind) dpdx
      real(r_kind) dpdeta
      
      allocate(q_ext(nVar,ics:ice,kcs:kce))
      allocate(q_p  (nVar,ics:ice,kcs:kce))
      
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
      
      allocate(P_ref (     ics:ice,kcs:kce))
      allocate(PL_ref(     ics:ice,kcs:kce))
      allocate(PR_ref(     ics:ice,kcs:kce))
      allocate(PB_ref(     ics:ice,kcs:kce))
      allocate(PT_ref(     ics:ice,kcs:kce))
      
      allocate(w_eta(     ics:ice,kcs:kce))
      
      allocate(Fe   (nVar,ids:ide+1,kds:kde  ))
      allocate(He   (nVar,ids:ide  ,kds:kde+1))
      
      allocate(src  (nVar,ids:ide,kds:kde))
      
      allocate(rho_p(     ids:ide,kds:kde))
      
      allocate(eig_x (ics:ice,kcs:kce))
      allocate(eig_z (ics:ice,kcs:kce))
      
      P     = FillValue
      P_ref = FillValue
      F     = FillValue
      H     = FillValue
      
      ! Set reference pressure
      q_ext = ref%q
      
      do k = kcs,kce
        do i = ics,ice
          P_ref(i,k) = calc_pressure(sqrtG(i,k),q_ext(:,i,k))
        enddo
      enddo
      
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
      
      ! Calculate P, F, H and eigenvalues
      do k = kds,kde
        do i = ids,ide
          P    (  i,k) = calc_pressure(sqrtG(i,k),q_ext(:,i,k))
          F    (:,i,k) = calc_F(sqrtG(i,k),q_ext(:,i,k),P(i,k))
          H    (:,i,k) = calc_H(sqrtG(i,k),G13(i,k),q_ext(:,i,k),P(i,k),P_ref(i,k))
          !H    (:,i,k) = calc_H_w_eta(sqrtG(i,k),G13(i,k),q_ext(:,i,k),P(i,k),P_ref(i,k),w_eta(i,k))
          
          eig_x(  i,k) = calc_eigenvalue_x(sqrtG(i,k)         ,q_ext(:,i,k))
          eig_z(  i,k) = calc_eigenvalue_z(sqrtG(i,k),G13(i,k),q_ext(:,i,k))
        enddo
      enddo
      
      ! Reconstruct FL, FR, HB and HT
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
            q_weno = F(iVar,im2:ip2,k)
            dir = -1
            call WENO_limiter(FL(iVar,i,k),q_weno,dir)
            dir = 1
            call WENO_limiter(FR(iVar,i,k),q_weno,dir)
            
            q_weno = H(iVar,i,km2:kp2)
            dir = -1
            call WENO_limiter(HB(iVar,i,k),q_weno,dir)
            dir = 1
            call WENO_limiter(HT(iVar,i,k),q_weno,dir)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! Reconstruct qL, qR, qB and qT
      !$OMP PARALLEL DO PRIVATE(kp1,km1,kp2,km2,i,ip1,im1,ip2,im2,iVar,q_weno,dir)
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
            q_weno = q_ext(iVar,im2:ip2,k)
            dir = -1
            call WENO_limiter(qL(iVar,i,k),q_weno,dir)
            dir = 1
            call WENO_limiter(qR(iVar,i,k),q_weno,dir)
            q_weno = q_ext(iVar,i,km2:kp2)
            dir = -1
            call WENO_limiter(qB(iVar,i,k),q_weno,dir)
            dir = 1
            call WENO_limiter(qT(iVar,i,k),q_weno,dir)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      rho_p = ( stat%q(1,ids:ide,kds:kde) + stat%q(5,ids:ide,kds:kde) &
              - ref%q (1,ids:ide,kds:kde) - ref%q (5,ids:ide,kds:kde) ) / sqrtG(ids:ide,kds:kde)
      
      where(abs(rho_p)<=1.E-13)rho_p=0.
      
      src(3,ids:ide,kds:kde) = src(3,ids:ide,kds:kde) - sqrtG(ids:ide,kds:kde) * rho_p(ids:ide,kds:kde) * gravity
      
      ! calc x flux
      !$OMP PARALLEL DO PRIVATE(i,im1,maxeigen_x,iVar)
      do k = kds,kde
        !do i = ids,ide+1
        do i = ids+1,ide
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
        !do k = kds,kde+1
        do k = kds+1,kde
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
      
      ! initialize source terms
      src = 0.
      
      ! Set no flux and nonreflecting boundary condition
      !$OMP PARALLEL DO PRIVATE(i)
      do k = kds,kde
        do i = ids,ide
          PL(i,k) = calc_pressure(sqrtGL(i,k),qL(:,i,k))
          PR(i,k) = calc_pressure(sqrtGR(i,k),qR(:,i,k))
          PB(i,k) = calc_pressure(sqrtGB(i,k),qB(:,i,k))
          PT(i,k) = calc_pressure(sqrtGT(i,k),qT(:,i,k))
        enddo
      enddo
      !$OMP END PARALLEL DO
      call bdy_condition(Fe,He,P,stat%q,ref%q,src)
      
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
    
    subroutine bdy_condition(Fe,He,P,q,q_ref,src)
      real(r_kind), dimension(nVar,ids:ide+1,kds:kde  ), intent(inout) :: Fe
      real(r_kind), dimension(nVar,ids:ide  ,kds:kde+1), intent(inout) :: He
      real(r_kind), dimension(nVar,ics:ice  ,kcs:kce  ), intent(in   ) :: P
      real(r_kind), dimension(nVar,ics:ice  ,kcs:kce  ), intent(in   ) :: q
      real(r_kind), dimension(nVar,ics:ice  ,kcs:kce  ), intent(in   ) :: q_ref
      real(r_kind), dimension(nVar,ids:ide  ,kds:kde  ), intent(inout) :: src
      
      integer(i_kind), parameter :: vs = 1
      integer(i_kind), parameter :: ve = 5
      integer(i_kind), parameter :: bdy_width = 40
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
        Fe(1,ids  ,kds:kde) = 0
        Fe(2,ids  ,kds:kde) = sqrtGL(ids,kds:kde) * PL(ids,kds:kde)
        Fe(3,ids  ,kds:kde) = 0
        Fe(4,ids  ,kds:kde) = 0
        Fe(5,ids  ,kds:kde) = 0
        
        Fe(1,ide+1,kds:kde) = 0
        Fe(2,ide+1,kds:kde) = sqrtGR(ide,kds:kde) * PR(ide,kds:kde)
        Fe(3,ide+1,kds:kde) = 0
        Fe(4,ide+1,kds:kde) = 0
        Fe(5,ide+1,kds:kde) = 0
        
        He(1,ids:ide,kds  ) = 0
        He(2,ids:ide,kds  ) = sqrtGB(ids:ide,kds) * G13B(ids:ide,kds) * PB(ids:ide,kds)
        He(3,ids:ide,kds  ) = PB(ids:ide,kds) - PB_ref(ids:ide,kds)
        He(4,ids:ide,kds  ) = 0
        He(5,ids:ide,kds  ) = 0
        
        He(1,ids:ide,kde+1) = 0
        He(2,ids:ide,kde+1) = sqrtGT(ids:ide,kde) * G13T(ids:ide,kde) * PT(ids:ide,kde)
        He(3,ids:ide,kde+1) = PT(ids:ide,kde) - PT_ref(ids:ide,kde)
        He(4,ids:ide,kde+1) = 0
        He(5,ids:ide,kde+1) = 0
      elseif(case_num==2)then
        Fe(1,ids,kds:kde) = 0
        Fe(2,ids,kds:kde) = q_ref(2,ids,kds:kde)**2 / q_ref(1,ids,kds:kde) + sqrtGL(ids,kds:kde) * PL(ids,kds:kde)
        Fe(3,ids,kds:kde) = 0
        Fe(4,ids,kds:kde) = q_ref(4,ids,kds:kde) * q_ref(2,ids,kds:kde) /  q_ref(1,ids,kds:kde)
        Fe(5,ids,kds:kde) = 0
        
        Fe(1,ide+1,kds:kde) = 0
        Fe(2,ide+1,kds:kde) = q_ref(2,ide,kds:kde)**2 / q_ref(1,ide,kds:kde) + sqrtGR(ide,kds:kde) * PR(ide,kds:kde)
        Fe(3,ide+1,kds:kde) = 0
        Fe(4,ide+1,kds:kde) = q_ref(4,ide,kds:kde) * q_ref(2,ide,kds:kde) /  q_ref(1,ide,kds:kde)
        Fe(5,ide+1,kds:kde) = 0
        
        He(1,ids:ide,kds) = 0
        He(2,ids:ide,kds) = sqrtGB(ids:ide,kds) * G13B(ids:ide,kds) * PB(ids:ide,kds)
        He(3,ids:ide,kds) = PB(ids:ide,kds) - PB_ref(ids:ide,kds)
        He(4,ids:ide,kds) = 0
        He(5,ids:ide,kds) = 0
        
        He(1,ids:ide,kde+1) = 0
        He(2,ids:ide,kde+1) = sqrtGT(ids:ide,kde) * G13T(ids:ide,kde) * PT(ids:ide,kde)
        He(3,ids:ide,kde+1) = PT(ids:ide,kde) - PT_ref(ids:ide,kde)
        He(4,ids:ide,kde+1) = 0
        He(5,ids:ide,kde+1) = 0
        
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
      if(abs(p_pert)/p_ref<1.e-14)p_pert=0
      
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
      if(abs(p_pert)/p_ref<1.e-14)p_pert=0
      
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

