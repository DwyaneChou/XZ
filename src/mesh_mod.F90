module mesh_mod
  use constants_mod
  use parameters_mod
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
  real(r_kind), dimension(:,:), allocatable :: dxidz
  real(r_kind), dimension(:,:), allocatable :: dzdx
  real(r_kind), dimension(:,:), allocatable :: dzsdx
  real(r_kind), dimension(:,:), allocatable :: detadx
  
  real(r_kind), dimension(:,:), allocatable :: sqrtG
  real(r_kind), dimension(:,:), allocatable :: G13
  
  real(r_kind), dimension(:,:), allocatable :: zs
  
  real(r_kind) deta
  
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
      allocate( dxidz  (ics:ice,kcs:kce) )
      allocate( dzdx   (ics:ice,kcs:kce) )
      allocate( dzsdx  (ics:ice,kcs:kce) )
      allocate( detadx (ics:ice,kcs:kce) )
      
      allocate( sqrtG  (ics:ice,kcs:kce) )
      allocate( G13    (ics:ice,kcs:kce) )
      allocate( zs     (ics:ice,kcs:kce) )
      
      x       = FillValue
      
      eta     = FillValue
      xi      = FillValue
      z       = FillValue
      
      dzdxi   = FillValue
      dxideta = FillValue
      dzdeta  = FillValue
      
      detadxi = FillValue
      dxidz   = FillValue
      dzsdx   = FillValue
      dzdx    = FillValue
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
          x (i,k) = ( real(i) - 0.5 ) * dx
        enddo
      enddo
    end subroutine init_horizontal_mesh
    
    subroutine init_vertical_distribiution
      integer i,k
      
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
        dxideta = z_max - z_min
        detadxi = 1. / dxideta
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
      endif
      
      ! Check necessary fields
      if(any(eta     == FillValue))stop 'eta     is not fully filled'
      if(any(xi      == FillValue))stop 'xi      is not fully filled'
      if(any(dxideta == FillValue))stop 'dxideta is not fully filled'
      if(any(detadxi == FillValue))stop 'detadxi is not fully filled'
      
    end subroutine init_vertical_distribiution
    
    subroutine init_vertical_coordinate
      integer i,k
      
      ! Check necessary fields
      if(any(zs    == FillValue)) stop 'zs    is not fully filled'
      if(any(dzsdx == FillValue)) stop 'dzsdx is not fully filled'
      
      if(vertical_coordinate==1)then
        ! Gal-Chen coordinate
        do k = kcs,kce
          do i = ics,ice
            dzdxi(i,k) = ( z_max + zs(i,k) ) / z_max
            dxidz(i,k) = z_max / ( z_max + zs(i,k) )
            dzdx (i,k) = dzsdx(i,k)
            
            dzdeta(i,k) = dzdxi(i,k) * dxideta(i,k)
            detadx(i,k) = detadxi(i,k) * dxidz(i,k) * dzdx(i,k)
            
            z(i,k) = ( z_max + zs(i,k) ) / z_max * xi(i,k)  + zs(i,k)
          enddo
        enddo
        
        sqrtG = dzdeta
        G13   = detadx
      else
        stop 'Unknown vertical_coordinate'
      endif
      
    end subroutine init_vertical_coordinate
end module mesh_mod
    