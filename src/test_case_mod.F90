module test_case_mod
  use constants_mod
  use parameters_mod
  use mesh_mod
  use stat_mod
  use tend_mod
  implicit none
  
  real(r_kind), dimension(:,:,:), allocatable :: q_ref      ! reference q
  
    contains
    subroutine init_test_case
      integer(i_kind) iT
      integer(i_kind) i,k
      
      allocate(q_ref(nVar,ics:ice,kcs:kce))
      
      iT = 0
      
      if(case_num==1)then
        call thermal_bubble(stat(iT))
      elseif(case_num==2)then
        call schar_mountain(stat(iT))
      elseif(case_num==3)then
        call density_current(stat(iT))
      else
        stop 'Unknown case_num'
      endif
      
    end subroutine init_test_case
    
    ! thermal_bubble according to Li 2013
    subroutine thermal_bubble(stat)
      type(stat_field),intent(inout) :: stat
      
      real(r_kind), dimension(:,:,:), allocatable :: rho
      real(r_kind), dimension(:,:,:), allocatable :: theta
      real(r_kind), dimension(:,:,:), allocatable :: u
      real(r_kind), dimension(:,:,:), allocatable :: w
      real(r_kind), dimension(:,:,:), allocatable :: q
      
      real(r_kind), dimension(:,:,:), allocatable :: r
      real(r_kind), dimension(:,:,:), allocatable :: exner
      real(r_kind), dimension(:,:,:), allocatable :: p
      real(r_kind), dimension(:,:,:), allocatable :: T
      
      real(r_kind), dimension(:,:,:), allocatable :: dexner
      
      real(r_kind) theta_bar
      real(r_kind) dtheta
      real(r_kind) R_bubble
      real(r_kind) x0
      real(r_kind) z0
      
      integer i,j,k,iVar
      
      allocate(rho  (nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(theta(nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(u    (nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(w    (nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(q    (nQuadPointsOnCell,ics:ice,kcs:kce))
      
      allocate(r    (nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(exner(nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(p    (nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(T    (nQuadPointsOnCell,ics:ice,kcs:kce))
      
      allocate(dexner(nQuadPointsOnCell,ics:ice,kcs:kce))
      
      zs    = 0.
      dzsdx = 0.
      
      call init_vertical_coordinate
      
      theta_bar = 300.
      dtheta    = 2.
      R_bubble  = 2000.
      x0        = 10000.
      z0        = 2000.
      
      do k = kcs,kce
        do i = ics,ice
          do j = 1,nQuadPointsOncell
            r     (j,i,k) = sqrt( ( x(j,i,k) - x0 )**2 + ( z(j,i,k) - z0 )**2 )
            theta (j,i,k) = theta_bar
            dexner(j,i,k) = -gravity / ( Cpd * theta(j,i,k) )
          enddo
        enddo
      enddo
      
      do i = ics,ice
        do k = kcs,kce
          do j = 1,nQuadPointsOnCell
            exner(j,i,k) = 1. + dexner(j,i,k) * z(j,i,k)
            p    (j,i,k) = p0 * exner(j,i,k)**(Cpd/Rd)
            T    (j,i,k) = exner(j,i,k) * theta(j,i,k)
            rho  (j,i,k) = p(j,i,k) / ( Rd * T(j,i,k) )
            u    (j,i,k) = 0.
            w    (j,i,k) = 0.
            q    (j,i,k) = 0.
          enddo
        enddo
      enddo
      
      print*,'Set reference fields'
      do k = kcs,kce
        do i = ics,ice
          q_ref(1,i,k) = Gaussian_quadrature_2d(sqrtG(:,i,k)*rho(:,i,k))
          q_ref(2,i,j) = 0.
          q_ref(3,i,j) = 0.
          q_ref(4,i,k) = q_ref(1,i,k) * theta_bar
          q_ref(5,i,j) = 0.
        enddo
      enddo
      
      ref%q = q_ref
      
      ! Set theta perturbation
      do k = kcs,kce
        do i = ics,ice
          do j = 1,nQuadPointsOnCell
            theta(j,i,k) = theta_bar + dtheta * max( 0., 1. - r(j,i,k) / R_bubble )
          enddo
        enddo
      enddo
      
      do k = kcs,kce
        do i = ics,ice
          stat%q(1,i,k) = Gaussian_quadrature_2d( sqrtG(:,i,k) * rho(:,i,k)                )
          stat%q(2,i,k) = Gaussian_quadrature_2d( sqrtG(:,i,k) * rho(:,i,k) * u    (:,i,k) )
          stat%q(3,i,k) = Gaussian_quadrature_2d( sqrtG(:,i,k) * rho(:,i,k) * w    (:,i,k) )
          stat%q(4,i,k) = Gaussian_quadrature_2d( sqrtG(:,i,k) * rho(:,i,k) * theta(:,i,k) )
          stat%q(5,i,k) = Gaussian_quadrature_2d( sqrtG(:,i,k) * rho(:,i,k) * q    (:,i,k) )
        enddo
      enddo
      
      print*,'max/min value of sqrtG           ',maxval(sqrtG        ),minval(sqrtG        )
      print*,'max/min value of G13             ',maxval(G13          ),minval(G13          )
      print*,'max/min value of sqrtG*rho       ',maxval(stat%q(1,:,:)),minval(stat%q(1,:,:))
      print*,'max/min value of sqrtG*rho*u     ',maxval(stat%q(2,:,:)),minval(stat%q(2,:,:))
      print*,'max/min value of sqrtG*rho*v     ',maxval(stat%q(3,:,:)),minval(stat%q(3,:,:))
      print*,'max/min value of sqrtG*rho*theta ',maxval(stat%q(4,:,:)),minval(stat%q(4,:,:))
      print*,'max/min value of sqrtG*rho*q     ',maxval(stat%q(5,:,:)),minval(stat%q(5,:,:))
      
    end subroutine thermal_bubble
    
    subroutine schar_mountain(stat)
      type(stat_field),intent(inout) :: stat
      
      real(r_kind), dimension(:,:,:), allocatable :: rho
      real(r_kind), dimension(:,:,:), allocatable :: theta
      real(r_kind), dimension(:,:,:), allocatable :: u
      real(r_kind), dimension(:,:,:), allocatable :: w
      real(r_kind), dimension(:,:,:), allocatable :: q
      
      real(r_kind), dimension(:,:,:), allocatable :: r
      real(r_kind), dimension(:,:,:), allocatable :: exner
      real(r_kind), dimension(:,:,:), allocatable :: p
      real(r_kind), dimension(:,:,:), allocatable :: T
      
      real(r_kind), dimension(:,:,:), allocatable :: dexner
      
      real(r_kind), dimension(:,:,:), allocatable :: theta_bar
      
      real(r_kind) dtheta
      real(r_kind) a0
      real(r_kind) h0
      real(r_kind) lambda0
      real(r_kind) theta0
      real(r_kind) N0
      
      integer i,j,k,iVar
      
      allocate(rho  (nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(theta(nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(u    (nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(w    (nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(q    (nQuadPointsOnCell,ics:ice,kcs:kce))
      
      allocate(r    (nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(exner(nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(p    (nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(T    (nQuadPointsOnCell,ics:ice,kcs:kce))
      
      allocate(dexner(nQuadPointsOnCell,ics:ice,kcs:kce))
      
      allocate(theta_bar(nQuadPointsOnCell,ics:ice,kcs:kce))
      
      h0      = 250.
      a0      = 5000.
      lambda0 = 4000.
      theta0  = 280.
      N0      = 0.01
      
      zs    = h0 * exp( -( x / a0 )**2 ) * cos( pi * x / lambda0 )**2
      dzsdx = -((2.*h0*Cos((pi*x)/lambda0)*(lambda0*x*Cos((pi*x)/lambda0) + a0**2*pi*Sin((pi*x)/lambda0)))/(Exp(x**2/a0**2)*(a0**2*lambda0)))
      
      call init_vertical_coordinate
      
      do k = kcs,kce
        do i = ics,ice
          do j = 1,nQuadPointsOncell
            theta_bar(j,i,k) = theta0 * exp( N0**2 /gravity * z(j,i,k) )
            theta    (j,i,k) = theta_bar(j,i,k)
            dexner   (j,i,k) = -gravity / ( Cpd * theta(j,i,k) )
          enddo
        enddo
      enddo
      
      ! Initialize fields
      do k = kcs,kce
        do i = ics,ice
          do j = 1,nQuadPointsOnCell
            exner(j,i,k) = 1. + gravity**2 / ( cpd * theta0 * N0**2 ) * ( exp( -N0**2 / gravity * z(j,i,k) ) - 1. )
            p    (j,i,k) = p0 * exner(j,i,k)**(Cpd/Rd)
            T    (j,i,k) = exner(j,i,k) * theta(j,i,k)
            rho  (j,i,k) = p(j,i,k) / ( Rd * T(j,i,k) )
            u    (j,i,k) = 10.
            w    (j,i,k) = 0.
            q    (j,i,k) = 0.
          enddo
        enddo
      enddo
      
      print*,'Set reference fields'
      do k = kcs,kce
        do i = ics,ice
          q_ref(1,i,k) = Gaussian_quadrature_2d( sqrtG(:,i,k) * rho(:,i,k)                )
          q_ref(2,i,k) = Gaussian_quadrature_2d( sqrtG(:,i,k) * rho(:,i,k) * u    (:,i,k) )
          q_ref(3,i,k) = Gaussian_quadrature_2d( sqrtG(:,i,k) * rho(:,i,k) * w    (:,i,k) )
          q_ref(4,i,k) = Gaussian_quadrature_2d( sqrtG(:,i,k) * rho(:,i,k) * theta(:,i,k) )
          q_ref(5,i,k) = Gaussian_quadrature_2d( sqrtG(:,i,k) * rho(:,i,k) * q    (:,i,k) )
        enddo
      enddo
      
      ref%q = q_ref
      
      do k = kcs,kce
        do i = ics,ice
          stat%q(1,i,k) = Gaussian_quadrature_2d( sqrtG(:,i,k) * rho(:,i,k)                )
          stat%q(2,i,k) = Gaussian_quadrature_2d( sqrtG(:,i,k) * rho(:,i,k) * u    (:,i,k) )
          stat%q(3,i,k) = Gaussian_quadrature_2d( sqrtG(:,i,k) * rho(:,i,k) * w    (:,i,k) )
          stat%q(4,i,k) = Gaussian_quadrature_2d( sqrtG(:,i,k) * rho(:,i,k) * theta(:,i,k) )
          stat%q(5,i,k) = Gaussian_quadrature_2d( sqrtG(:,i,k) * rho(:,i,k) * q    (:,i,k) )
        enddo
      enddo
      
      print*,'max/min value of sqrtG           ',maxval(sqrtG        ),minval(sqrtG        )
      print*,'max/min value of G13             ',maxval(G13          ),minval(G13          )
      print*,'max/min value of sqrtG*rho       ',maxval(stat%q(1,:,:)),minval(stat%q(1,:,:))
      print*,'max/min value of sqrtG*rho*u     ',maxval(stat%q(2,:,:)),minval(stat%q(2,:,:))
      print*,'max/min value of sqrtG*rho*w     ',maxval(stat%q(3,:,:)),minval(stat%q(3,:,:))
      print*,'max/min value of sqrtG*rho*theta ',maxval(stat%q(4,:,:)),minval(stat%q(4,:,:))
      print*,'max/min value of sqrtG*rho*q     ',maxval(stat%q(5,:,:)),minval(stat%q(5,:,:))
    end subroutine schar_mountain
    
    ! Density current according to Li 2013
    subroutine density_current(stat)
      type(stat_field),intent(inout) :: stat
      
      real(r_kind), dimension(:,:,:), allocatable :: rho
      real(r_kind), dimension(:,:,:), allocatable :: theta
      real(r_kind), dimension(:,:,:), allocatable :: u
      real(r_kind), dimension(:,:,:), allocatable :: w
      real(r_kind), dimension(:,:,:), allocatable :: q
      
      real(r_kind), dimension(:,:,:), allocatable :: r
      real(r_kind), dimension(:,:,:), allocatable :: exner
      real(r_kind), dimension(:,:,:), allocatable :: p
      real(r_kind), dimension(:,:,:), allocatable :: T
      
      real(r_kind), dimension(:,:,:), allocatable :: dexner
      
      real(r_kind) theta_bar
      real(r_kind) dtheta
      real(r_kind) xr
      real(r_kind) zr
      real(r_kind) x0
      real(r_kind) z0
      
      integer i,j,k,iVar
      
      allocate(rho  (nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(theta(nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(u    (nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(w    (nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(q    (nQuadPointsOnCell,ics:ice,kcs:kce))
      
      allocate(r    (nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(exner(nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(p    (nQuadPointsOnCell,ics:ice,kcs:kce))
      allocate(T    (nQuadPointsOnCell,ics:ice,kcs:kce))
      
      allocate(dexner(nQuadPointsOnCell,ics:ice,kcs:kce))
      
      zs    = 0.
      dzsdx = 0.
      
      call init_vertical_coordinate
      
      theta_bar = 300.
      dtheta    = -15.
      xr        = 4000.
      zr        = 2000.
      x0        = 0.
      z0        = 3000.
      
      do k = kcs,kce
        do i = ics,ice
          do j = 1,nQuadPointsOncell
            r     (j,i,k) = sqrt( ( x(j,i,k) - x0 )**2 + ( z(j,i,k) - z0 )**2 )
            theta (j,i,k) = theta_bar
            dexner(j,i,k) = -gravity / ( Cpd * theta(j,i,k) )
          enddo
        enddo
      enddo
      
      do i = ics,ice
        do k = kcs,kce
          do j = 1,nQuadPointsOnCell
            exner(j,i,k) = 1. + dexner(j,i,k) * z(j,i,k)
            p    (j,i,k) = p0 * exner(j,i,k)**(Cpd/Rd)
            T    (j,i,k) = exner(j,i,k) * theta(j,i,k)
            rho  (j,i,k) = p(j,i,k) / ( Rd * T(j,i,k) )
            u    (j,i,k) = 0.
            w    (j,i,k) = 0.
            q    (j,i,k) = 0.
          enddo
        enddo
      enddo
      
      print*,'Set reference fields'
      do k = kcs,kce
        do i = ics,ice
          q_ref(1,i,k) = Gaussian_quadrature_2d(sqrtG(:,i,k)*rho(:,i,k))
          q_ref(2,i,j) = 0.
          q_ref(3,i,j) = 0.
          q_ref(4,i,k) = q_ref(1,i,k) * theta_bar
          q_ref(5,i,j) = 0.
        enddo
      enddo
      
      ref%q = q_ref
      
      ! Set theta perturbation
      where(r<=1.)theta = theta_bar + dtheta * ( cos(pi*r) + 1. )/2.
      !where(r<=1.)theta = theta_bar + dtheta * cos(pi*r/2)**2
      
      viscosity_coef = 75
      
      do k = kcs,kce
        do i = ics,ice
          stat%q(1,i,k) = Gaussian_quadrature_2d( sqrtG(:,i,k) * rho(:,i,k)                )
          stat%q(2,i,k) = Gaussian_quadrature_2d( sqrtG(:,i,k) * rho(:,i,k) * u    (:,i,k) )
          stat%q(3,i,k) = Gaussian_quadrature_2d( sqrtG(:,i,k) * rho(:,i,k) * w    (:,i,k) )
          stat%q(4,i,k) = Gaussian_quadrature_2d( sqrtG(:,i,k) * rho(:,i,k) * theta(:,i,k) )
          stat%q(5,i,k) = Gaussian_quadrature_2d( sqrtG(:,i,k) * rho(:,i,k) * q    (:,i,k) )
        enddo
      enddo
      
      print*,'max/min value of sqrtG           ',maxval(sqrtG        ),minval(sqrtG        )
      print*,'max/min value of G13             ',maxval(G13          ),minval(G13          )
      print*,'max/min value of sqrtG*rho       ',maxval(stat%q(1,:,:)),minval(stat%q(1,:,:))
      print*,'max/min value of sqrtG*rho*u     ',maxval(stat%q(2,:,:)),minval(stat%q(2,:,:))
      print*,'max/min value of sqrtG*rho*v     ',maxval(stat%q(3,:,:)),minval(stat%q(3,:,:))
      print*,'max/min value of sqrtG*rho*theta ',maxval(stat%q(4,:,:)),minval(stat%q(4,:,:))
      print*,'max/min value of sqrtG*rho*q     ',maxval(stat%q(5,:,:)),minval(stat%q(5,:,:))
      
    end subroutine density_current
    
end module test_case_mod
    