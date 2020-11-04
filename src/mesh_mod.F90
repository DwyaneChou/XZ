module mesh_mod
  use constants_mod
  use parameters_mod
  use reconstruction_mod, only: WENO_limiter
  use math_mod
  implicit none
  
  real(r_kind), dimension(:,:), allocatable :: x
  
  real(r_kind), dimension(:,:), allocatable :: eta
  real(r_kind), dimension(:,:), allocatable :: xi
  real(r_kind), dimension(:,:), allocatable :: z
  
  ! For dzdeta ( sqrt(G) )
  real(r_kind), dimension(:,:), allocatable :: dzdxi
  real(r_kind), dimension(:,:), allocatable :: dxideta
  real(r_kind), dimension(:,:), allocatable :: dzdeta
  
  ! For detadx ( G13 )
  real(r_kind), dimension(:,:), allocatable :: detadxi
  real(r_kind), dimension(:,:), allocatable :: dxidx
  real(r_kind), dimension(:,:), allocatable :: dzsdx
  real(r_kind), dimension(:,:), allocatable :: detadx
  
  ! sqrt(G)
  real(r_kind), dimension(:,:), allocatable :: sqrtG
  real(r_kind), dimension(:,:), allocatable :: sqrtGL
  real(r_kind), dimension(:,:), allocatable :: sqrtGR
  real(r_kind), dimension(:,:), allocatable :: sqrtGB
  real(r_kind), dimension(:,:), allocatable :: sqrtGT
  
  ! G13
  real(r_kind), dimension(:,:), allocatable :: G13
  real(r_kind), dimension(:,:), allocatable :: G13L
  real(r_kind), dimension(:,:), allocatable :: G13R
  real(r_kind), dimension(:,:), allocatable :: G13B
  real(r_kind), dimension(:,:), allocatable :: G13T
  
  real(r_kind), dimension(:,:), allocatable :: zs
  
  real(r_kind) deta
  
  ! Extened fields
  real(r_kind), dimension(:,:), allocatable :: x_ext
  
  real(r_kind), dimension(:,:), allocatable :: eta_ext
  real(r_kind), dimension(:,:), allocatable :: xi_ext
  real(r_kind), dimension(:,:), allocatable :: z_ext
  
  ! For dzdeta ( sqrt(G) )
  real(r_kind), dimension(:,:), allocatable :: dzdxi_ext
  real(r_kind), dimension(:,:), allocatable :: dxideta_ext
  real(r_kind), dimension(:,:), allocatable :: dzdeta_ext
  
  ! For detadx ( G13 )
  real(r_kind), dimension(:,:), allocatable :: detadxi_ext
  real(r_kind), dimension(:,:), allocatable :: dxidx_ext
  real(r_kind), dimension(:,:), allocatable :: dzsdx_ext
  real(r_kind), dimension(:,:), allocatable :: detadx_ext
  
  real(r_kind), dimension(:,:), allocatable :: zs_ext
  
  
    contains
    subroutine init_mesh
      
      allocate( x      (ics:ice,kcs:kce) )
      
      allocate( eta    (ics:ice,kcs:kce) )
      allocate( xi     (ics:ice,kcs:kce) )
      allocate( z      (ics:ice,kcs:kce) )
      
      allocate( dzdxi  (ics:ice,kcs:kce) )
      allocate( dxideta(ics:ice,kcs:kce) )
      allocate( dzdeta (ics:ice,kcs:kce) )
      
      allocate( detadxi(ics:ice,kcs:kce) )
      allocate( dxidx  (ics:ice,kcs:kce) )
      allocate( dzsdx  (ics:ice,kcs:kce) )
      allocate( detadx (ics:ice,kcs:kce) )
      
      allocate( sqrtG  (ics:ice,kcs:kce) )
      allocate( sqrtGL (ics:ice,kcs:kce) )
      allocate( sqrtGR (ics:ice,kcs:kce) )
      allocate( sqrtGB (ics:ice,kcs:kce) )
      allocate( sqrtGT (ics:ice,kcs:kce) )
      
      allocate( G13    (ics:ice,kcs:kce) )
      allocate( G13L   (ics:ice,kcs:kce) )
      allocate( G13R   (ics:ice,kcs:kce) )
      allocate( G13B   (ics:ice,kcs:kce) )
      allocate( G13T   (ics:ice,kcs:kce) )
      
      allocate( zs     (ics:ice,kcs:kce) )
      
      ! Extended fields
      allocate( x_ext      (ies:iee,kes:kee) )
      allocate( eta_ext    (ies:iee,kes:kee) )
      allocate( xi_ext     (ies:iee,kes:kee) )
      allocate( z_ext      (ies:iee,kes:kee) )
      allocate( dzdxi_ext  (ies:iee,kes:kee) )
      allocate( dxideta_ext(ies:iee,kes:kee) )
      allocate( dzdeta_ext (ies:iee,kes:kee) )
      allocate( detadxi_ext(ies:iee,kes:kee) )
      allocate( dxidx_ext  (ies:iee,kes:kee) )
      allocate( dzsdx_ext  (ies:iee,kes:kee) )
      allocate( detadx_ext (ies:iee,kes:kee) )
      allocate( zs_ext     (ies:iee,kes:kee) )
      
      x       = FillValue
      
      eta     = FillValue
      xi      = FillValue
      z       = FillValue
      
      dzdxi   = FillValue
      dxideta = FillValue
      dzdeta  = FillValue
      
      detadxi = FillValue
      dxidx   = FillValue
      dzsdx   = FillValue
      detadx  = FillValue
      
      sqrtG   = FillValue
      G13     = FillValue
      
      zs      = FillValue
      
      call init_horizontal_mesh
      call init_vertical_distribiution
      
    end subroutine init_mesh
    
    subroutine init_horizontal_mesh
      integer i,k
      do k = kcs,kce
        do i = ics,ice
          x (i,k) = x_min + ( real(i) - 0.5 ) * dx
        enddo
      enddo
      
      do k = kes,kee
        do i =ies,iee
          x_Ext (i,k) = x_min + ( real(i) - 1. ) * dx / 2.
        enddo
      enddo
      
    end subroutine init_horizontal_mesh
    
    subroutine init_vertical_distribiution
      integer(i_kind) i,k
      
      real(r_kind) dxi
      real(r_kind) xx ! x parameter in atan vertical distribiution
      real(r_kind) xx_min
      real(r_kind) xx_max
    
      deta = 1. / nz
      
      if( vertical_distribution == 1 )then
        dxi = ( z_max - z_min ) / nz
        
        do k = kcs,kce
          eta(:,k) = ( real(k) - 0.5 ) * deta
          xi (:,k) = ( real(k) - 0.5 ) * dxi
        enddo
        
        do k = kes,kee
          eta_ext(:,k) = ( real(k) - 1. ) * deta / 2.
          xi_ext (:,k) = ( real(k) - 1. ) * dxi  / 2.
        enddo
        
        dxideta = z_max - z_min
        detadxi = 1. / dxideta
        
        dxideta_ext = z_max - z_min
        detadxi_ext = 1. / dxideta_Ext
      elseif( vertical_distribution == 2 )then
        do k = kcs,kce
          eta(:,k) = ( real(i) - 0.5 ) * deta
          xx       = 2. * pi * m_coef * ( eta(1,k) - 0.5 )
          xx_min   = -pi * m_coef
          xx_max   = xx
          xi     (:,k) = z_max / ( pi**2 * m_coef ) * ( xx_max * atan(xx_max) - xx_min * atan(xx_min) + log( sqrt( (1.+xx_min**2) / (1.+xx_max**2) ) ) ) + z_max * eta(1,k)
          dxideta(:,k) = 2. * z_max / pi * ( atan(xx) + pi / 2. )
          detadxi(:,k) = 1. / dxideta(:,k)
        enddo
        
        do k = kes,kee
          eta_ext    (:,k) = ( real(i) - 1. ) * deta / 2.
          xx               = 2. * pi * m_coef * ( eta_ext(1,k) - 0.5 )
          xx_min           = -pi * m_coef
          xx_max           = xx
          xi_ext     (:,k) = z_max / ( pi**2 * m_coef ) * ( xx_max * atan(xx_max) - xx_min * atan(xx_min) + log( sqrt( (1.+xx_min**2) / (1.+xx_max**2) ) ) ) + z_max * eta_ext(1,k)
          dxideta_ext(:,k) = 2. * z_max / pi * ( atan(xx) + pi / 2. )
          detadxi_ext(:,k) = 1. / dxideta_ext(:,k)
        enddo
      endif
      
      ! Check necessary fields
      if(any(eta     == FillValue))stop 'eta     is not fully filled'
      if(any(xi      == FillValue))stop 'xi      is not fully filled'
      if(any(dxideta == FillValue))stop 'dxideta is not fully filled'
      if(any(detadxi == FillValue))stop 'detadxi is not fully filled'
      
    end subroutine init_vertical_distribiution
    
    subroutine init_vertical_coordinate
      ! For reconstruction
      integer(i_kind), parameter :: nMertric = 2 ! number of metric tensor
      
      real(r_kind), dimension(:,:,:), allocatable :: q_ext ! Extended forecast variables
      real(r_kind), dimension(:,:,:), allocatable :: qL    ! Reconstructed q_(i-1/2,k)
      real(r_kind), dimension(:,:,:), allocatable :: qR    ! Reconstructed q_(i+1/2,k)
      real(r_kind), dimension(:,:,:), allocatable :: qB    ! Reconstructed q_(i,k-1/2)
      real(r_kind), dimension(:,:,:), allocatable :: qT    ! Reconstructed q_(i,k+1/2)
      
      real(r_kind), dimension(5) :: q_weno
      
      integer(i_kind) dir
      integer(i_kind) ip1,im1
      integer(i_kind) kp1,km1
      integer(i_kind) ip2,im2
      integer(i_kind) kp2,km2
      
      real(r_kind) angle
      
      integer i,k,iVar
      
      ! For Schar, 2002
      real(r_kind) :: H
      real(r_kind) :: s
      
      allocate(q_ext(nMertric,ics:ice,kcs:kce))
      
      allocate(qL   (nMertric,ids:ide,kds:kde))
      allocate(qR   (nMertric,ids:ide,kds:kde))
      allocate(qB   (nMertric,ids:ide,kds:kde))
      allocate(qT   (nMertric,ids:ide,kds:kde))
      
      ! Check necessary fields
      if(any(zs    == FillValue)) stop 'zs    is not fully filled'
      if(any(dzsdx == FillValue)) stop 'dzsdx is not fully filled'
      
      ! Calculate sqrtG and G13 by specific coordinate scheme
      if(vertical_coordinate==1)then
        ! Gal-Chen, 1975
        do k = kcs,kce
          do i = ics,ice
            dzdxi(i,k) = ( z_max - zs(i,k) ) / z_max
            dxidx(i,k) = ( xi(i,k) - z_max ) / ( z_max - zs(i,k) ) * dzsdx(i,k)
            
            dzdeta(i,k) = dzdxi(i,k) * dxideta(i,k)
            detadx(i,k) = detadxi(i,k) * dxidx(i,k)
            
            z(i,k) = ( z_max - zs(i,k) ) / z_max * xi(i,k)  + zs(i,k)
          enddo
        enddo
        
        ! Gal-Chen, 1975
        do k = kes,kee
          do i = ies,iee
            dzdxi_ext(i,k) = ( z_max - zs_ext(i,k) ) / z_max
            dxidx_ext(i,k) = ( xi_ext(i,k) - z_max ) / ( z_max - zs_ext(i,k) ) * dzsdx_ext(i,k)
            
            dzdeta_ext(i,k) = dzdxi_ext(i,k) * dxideta_ext(i,k)
            detadx_ext(i,k) = detadxi_ext(i,k) * dxidx_ext(i,k)
            
            z_ext(i,k) = ( z_max - zs_ext(i,k) ) / z_max * xi_ext(i,k)  + zs_ext(i,k)
          enddo
        enddo
      elseif(vertical_coordinate==2)then
        ! Schar, 2002
        H = z_max - z_min
        if(case_num==2)then
          s = 3000.
        else
          s = H / exp(1.)
        endif
        
        do k = kcs,kce
          do i = ics,ice
            dzdxi(i,k) = 1. - zs(i,k) * cosh( ( H - xi(i,k) ) / s ) / ( s * sinh( H / s ) )
            dxidx(i,k) = -sinh( ( z_max - xi(i,k) ) / s ) / sinh( z_max / s ) &
                       / ( 1. - zs(i,k) / s * cosh( ( z_max - xi(i,k) ) / s ) / sinh( z_max / s ) ) * dzsdx(i,k)
            
            dzdeta(i,k) = dzdxi(i,k) * dxideta(i,k)
            detadx(i,k) = detadxi(i,k) * dxidx(i,k)
            
            z(i,k) = xi(i,k) + zs(i,k) * sinh( ( H - xi(i,k) ) / s ) / sinh( H / s )
          enddo
        enddo
        
        do k = kes,kee
          do i = ies,iee
            dzdxi_ext(i,k) = 1. - zs_ext(i,k) * cosh( ( H - xi_ext(i,k) ) / s ) / ( s * sinh( H / s ) )
            dxidx_ext(i,k) = -sinh( ( z_max - xi_ext(i,k) ) / s ) / sinh( z_max / s ) &
                       / ( 1. - zs_ext(i,k) / s * cosh( ( z_max - xi_ext(i,k) ) / s ) / sinh( z_max / s ) ) * dzsdx_ext(i,k)
            
            dzdeta_ext(i,k) = dzdxi_ext(i,k) * dxideta_ext(i,k)
            detadx_ext(i,k) = detadxi_ext(i,k) * dxidx_ext(i,k)
            
            z_ext(i,k) = xi_ext(i,k) + zs_ext(i,k) * sinh( ( H - xi_ext(i,k) ) / s ) / sinh( H / s )
          enddo
        enddo
      else
        stop 'Unknown vertical_coordinate'
      endif
        
      !sqrtG = dzdeta
      !G13   = detadx
        
      ! Reconstruct mertric tensor by analytical value
      do k = kds,kde
        do i = ids,ide
      !do k = kcs,kce
      !  do i = ics,ice
          sqrtG (i,k) = dzdeta_ext(i*2  ,k*2  )
          sqrtGL(i,k) = dzdeta_ext(i*2-1,k*2  )
          sqrtGR(i,k) = dzdeta_ext(i*2+1,k*2  )
          sqrtGB(i,k) = dzdeta_ext(i*2  ,k*2-1)
          sqrtGT(i,k) = dzdeta_ext(i*2  ,k*2+1)
          
          G13   (i,k) = detadx_ext(i*2  ,k*2  )
          G13L  (i,k) = detadx_ext(i*2-1,k*2  )
          G13R  (i,k) = detadx_ext(i*2+1,k*2  )
          G13B  (i,k) = detadx_ext(i*2  ,k*2-1)
          G13T  (i,k) = detadx_ext(i*2  ,k*2+1)
        enddo
      enddo
      
      !! Reconstruct mertric tensor by  Mid-point average
      !do k = kds,kde
      !  do i = ids,ide
      !    sqrtGL(i,k) = ( sqrtG(i-1,k  ) + sqrtG(i  ,k  ) ) / 2.
      !    sqrtGR(i,k) = ( sqrtG(i  ,k  ) + sqrtG(i+1,k  ) ) / 2.
      !    sqrtGB(i,k) = ( sqrtG(i  ,k-1) + sqrtG(i  ,k  ) ) / 2.
      !    sqrtGT(i,k) = ( sqrtG(i  ,k  ) + sqrtG(i  ,k+1) ) / 2.
      !    
      !    G13L  (i,k) = ( G13(i-1,k  ) + G13(i  ,k  ) ) / 2.
      !    G13R  (i,k) = ( G13(i  ,k  ) + G13(i+1,k  ) ) / 2.
      !    G13B  (i,k) = ( G13(i  ,k-1) + G13(i  ,k  ) ) / 2.
      !    G13T  (i,k) = ( G13(i  ,k  ) + G13(i  ,k+1) ) / 2.
      !  enddo
      !enddo
      
      ! Reconstruct mertric tensor by WENO
      q_ext(1,:,:) = sqrtG
      q_ext(2,:,:) = G13
      !$OMP PARALLEL DO PRIVATE(i,ip1,im1,ip2,im2,iVar,q_weno,dir,kp1,km1,kp2,km2)
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
          do iVar = 1,nMertric
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
        enddo
      enddo
      !$OMP END PARALLEL DO
      sqrtGL(ids:ide,kds:kde) = qL(1,:,:)
      sqrtGR(ids:ide,kds:kde) = qR(1,:,:)
      sqrtGB(ids:ide,kds:kde) = qB(1,:,:)
      sqrtGT(ids:ide,kds:kde) = qT(1,:,:)
      
      G13L  (ids:ide,kds:kde) = qL(2,:,:)
      G13R  (ids:ide,kds:kde) = qR(2,:,:)
      G13B  (ids:ide,kds:kde) = qB(2,:,:)
      G13T  (ids:ide,kds:kde) = qT(2,:,:)
      
    end subroutine init_vertical_coordinate
end module mesh_mod
    