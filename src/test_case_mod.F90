module test_case_mod
  use constants_mod
  use parameters_mod
  use mesh_mod
  use stat_mod
  use tend_mod
  use linear_integration_mod
  use spatial_operators_mod
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
      
      real(r_kind), dimension(:,:), allocatable :: rho
      real(r_kind), dimension(:,:), allocatable :: theta
      real(r_kind), dimension(:,:), allocatable :: u
      real(r_kind), dimension(:,:), allocatable :: w
      real(r_kind), dimension(:,:), allocatable :: q
      
      real(r_kind), dimension(:,:), allocatable :: r
      real(r_kind), dimension(:,:), allocatable :: exner
      real(r_kind), dimension(:,:), allocatable :: p
      real(r_kind), dimension(:,:), allocatable :: T
      
      real(r_kind), dimension(:,:), allocatable :: dexner
      
      real(r_kind) theta_bar
      real(r_kind) dtheta
      real(r_kind) R_bubble
      real(r_kind) x0
      real(r_kind) z0
      
      integer i,k,iVar
      
      allocate(rho  (ics:ice,kcs:kce))
      allocate(theta(ics:ice,kcs:kce))
      allocate(u    (ics:ice,kcs:kce))
      allocate(w    (ics:ice,kcs:kce))
      allocate(q    (ics:ice,kcs:kce))
      
      allocate(r    (ics:ice,kcs:kce))
      allocate(exner(ics:ice,kcs:kce))
      allocate(p    (ics:ice,kcs:kce))
      allocate(T    (ics:ice,kcs:kce))
      
      allocate(dexner(ics:ice,kcs:kce))
      
      zs    = 0.
      dzsdx = 0.
      
      zs_ext    = 0.
      dzsdx_ext = 0.
      
      call init_vertical_coordinate
      
      theta_bar = 300.
      dtheta    = 2.
      R_bubble  = 2000.
      x0        = 10000.
      z0        = 2000.
      
      do k = kcs,kce
        do i = ics,ice
          r     (i,k) = sqrt( ( x(i,k) - x0 )**2 + ( z(i,k) - z0 )**2 )
          theta (i,k) = theta_bar
          dexner(i,k) = -gravity / ( Cpd * theta(i,k) )
        enddo
      enddo
      
      !do i = ics,ice
      !  call spline2_integration(nz_ext-1,z(i,:),dexner(i,:),0.,0.,nz_ext,z(i,:),dexner(i,:))
      !  !call spline4_integration(nz_ext-1,z(i,:),dexner(i,:)      ,nz_ext,z(i,:),dexner(i,:))
      !  
      !  exner(i,1) = 1.
      !  
      !  do k = kds+1,kce
      !    exner(i,k) = exner(i,1) + sum(dexner(i,kds+1:k))
      !  enddo
      !  
      !  do k = kcs,kds-1
      !    exner(i,k) = 1. - sum(dexner(i,kds:k+1:-1))
      !  enddo
      !  
      !  do k = kcs,kce
      !    p    (i,k) = p0 * exner(i,k)**(Cpd/Rd)
      !    T    (i,k) = exner(i,k) * theta(i,k)
      !    rho  (i,k) = p(i,k) / ( Rd * T(i,k) )
      !    u    (i,k) = 0.
      !    w    (i,k) = 0.
      !    q    (i,k) = 0.
      !  enddo
      !enddo
      
      do i = ics,ice
        do k = kcs,kce
          exner(i,k) = 1. + dexner(i,k) * z(i,k)
          p    (i,k) = p0 * exner(i,k)**(Cpd/Rd)
          T    (i,k) = exner(i,k) * theta(i,k)
          rho  (i,k) = p(i,k) / ( Rd * T(i,k) )
          u    (i,k) = 0.
          w    (i,k) = 0.
          q    (i,k) = 0.
        enddo
      enddo
      
      print*,'Set reference fields'
      do k = kcs,kce
        q_ref(1,ids:ide,k) = sum(sqrtG(ids:ide,k)*rho(ids:ide,k)) / nx
      enddo
      
      q_ref(2,:,:) = 0.
      q_ref(3,:,:) = 0.
      q_ref(4,:,:) = 300.
      q_ref(5,:,:) = 0.
      
      do iVar = 2,nVar
        q_ref(iVar,:,:) = q_ref(iVar,:,:) * q_ref(1,:,:)
      enddo
      
      ref%q = q_ref
      
      ! Set theta perturbation
      do k = kcs,kce
        do i = ics,ice
          theta(i,k) = theta_bar + dtheta * max( 0., 1. - r(i,k) / R_bubble )
        enddo
      enddo
      
      stat%q(1,:,:) = sqrtG * rho
      stat%q(2,:,:) = sqrtG * rho * u
      stat%q(3,:,:) = sqrtG * rho * w
      stat%q(4,:,:) = sqrtG * rho * theta
      stat%q(5,:,:) = sqrtG * rho * q
      
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
      
      real(r_kind), dimension(:,:), allocatable :: rho
      real(r_kind), dimension(:,:), allocatable :: theta
      real(r_kind), dimension(:,:), allocatable :: u
      real(r_kind), dimension(:,:), allocatable :: w
      real(r_kind), dimension(:,:), allocatable :: q
      
      real(r_kind), dimension(:,:), allocatable :: r
      real(r_kind), dimension(:,:), allocatable :: exner
      real(r_kind), dimension(:,:), allocatable :: p
      real(r_kind), dimension(:,:), allocatable :: T
      
      real(r_kind), dimension(:,:), allocatable :: dexner
      
      real(r_kind), dimension(:,:), allocatable :: theta_bar
      
      real(r_kind) dtheta
      real(r_kind) a0
      real(r_kind) h0
      real(r_kind) lambda0
      real(r_kind) theta0
      real(r_kind) N0
      
      integer i,k,iVar
      
      allocate(rho  (ics:ice,kcs:kce))
      allocate(theta(ics:ice,kcs:kce))
      allocate(u    (ics:ice,kcs:kce))
      allocate(w    (ics:ice,kcs:kce))
      allocate(q    (ics:ice,kcs:kce))
      
      allocate(r    (ics:ice,kcs:kce))
      allocate(exner(ics:ice,kcs:kce))
      allocate(p    (ics:ice,kcs:kce))
      allocate(T    (ics:ice,kcs:kce))
      
      allocate(dexner(ics:ice,kcs:kce))
      
      allocate(theta_bar(ics:ice,kcs:kce))
      
      h0      = 250.
      a0      = 5000.
      lambda0 = 4000.
      theta0  = 280.
      N0      = 0.01
      
      zs    = h0 * exp( -( x / a0 )**2 ) * cos( pi * x / lambda0 )**2
      dzsdx = -((2.*h0*Cos((pi*x)/lambda0)*(lambda0*x*Cos((pi*x)/lambda0) + a0**2*pi*Sin((pi*x)/lambda0)))/(Exp(x**2/a0**2)*(a0**2*lambda0)))
      
      zs_ext    = h0 * exp( -( x_ext / a0 )**2 ) * cos( pi * x_ext / lambda0 )**2
      dzsdx_ext = -((2.*h0*Cos((pi*x_ext)/lambda0)*(lambda0*x_ext*Cos((pi*x_ext)/lambda0) + a0**2*pi*Sin((pi*x_ext)/lambda0)))/(Exp(x_ext**2/a0**2)*(a0**2*lambda0)))
      
      call init_vertical_coordinate
      
      !do k = kcs,kce
      !  do i = ics,ice
      !    theta_bar(i,k) = theta0 * exp( N0**2 /gravity * z(i,k) )
      !    theta    (i,k) = theta_bar(i,k)
      !    dexner   (i,k) = -gravity / ( Cpd * theta(i,k) )
      !  enddo
      !enddo
      !
      !! Set exner on top
      !do i = ics,ice
      !  exner(i,1) = 1.
      !  
      !  call spline2_integration(nz_ext-1,xi(i,:),dexner(i,:),0.,0.,nz_ext,xi(i,:),dexner(i,:))
      !  !call spline4_integration(nz_ext-1,z(i,:),dexner(i,:)      ,nz_ext,z(i,:),dexner(i,:))
      !  
      !  do k = kds+1,kce
      !    exner(i,k) = exner(i,1) + sum(dexner(i,kds+1:k))
      !  enddo
      !  
      !enddo
      !exner(:,kce) = sum( exner(ids:ide,kce) ) / nx
      !
      !do k = kcs,kce
      !  do i = ics,ice
      !    theta_bar(i,k) = theta0 * exp( N0**2 /gravity * z(i,k) )
      !    theta    (i,k) = theta_bar(i,k)
      !    dexner   (i,k) = -gravity / ( Cpd * theta(i,k) )
      !  enddo
      !enddo
      !
      !! Initialize fields
      !do i = ics,ice
      !  call spline2_integration(nz_ext-1,z(i,:),dexner(i,:),0.,0.,nz_ext,z(i,:),dexner(i,:))
      !  !call spline4_integration(nz_ext-1,z(i,:),dexner(i,:)      ,nz_ext,z(i,:),dexner(i,:))
      !  
      !  do k = kce,kcs,-1
      !    exner(i,k) = exner(i,kce) - sum(dexner(i,k:kce))
      !  enddo
      !  
      !  do k = kcs,kce
      !    p    (i,k) = p0 * exner(i,k)**(Cpd/Rd)
      !    T    (i,k) = exner(i,k) * theta(i,k)
      !    rho  (i,k) = p(i,k) / ( Rd * T(i,k) )
      !    u    (i,k) = 10.
      !    w    (i,k) = 0.
      !    q    (i,k) = 0.
      !  enddo
      !enddo
      
      do k = kcs,kce
        do i = ics,ice
          theta_bar(i,k) = theta0 * exp( N0**2 /gravity * z(i,k) )
          theta    (i,k) = theta_bar(i,k)
          dexner   (i,k) = -gravity / ( Cpd * theta(i,k) )
        enddo
      enddo
      
      ! Initialize fields
      do i = ics,ice
        do k = kcs,kce
          exner(i,k) = 1. + gravity**2 / ( cpd * theta0 * N0**2 ) * ( exp( -N0**2 / gravity * z(i,k) ) - 1. )
          p    (i,k) = p0 * exner(i,k)**(Cpd/Rd)
          T    (i,k) = exner(i,k) * theta(i,k)
          rho  (i,k) = p(i,k) / ( Rd * T(i,k) )
          u    (i,k) = 10.
          w    (i,k) = 0.
          q    (i,k) = 0.
        enddo
      enddo
      
      print*,'Set reference fields'
      q_ref(1,:,:) = sqrtG * rho
      q_ref(2,:,:) = sqrtG * rho * u
      q_ref(3,:,:) = sqrtG * rho * w
      q_ref(4,:,:) = sqrtG * rho * theta
      q_ref(5,:,:) = sqrtG * rho * q
      
      ref%q = q_ref
      
      stat%q(1,:,:) = sqrtG * rho
      stat%q(2,:,:) = sqrtG * rho * u
      stat%q(3,:,:) = sqrtG * rho * w
      stat%q(4,:,:) = sqrtG * rho * theta
      stat%q(5,:,:) = sqrtG * rho * q
      
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
      
      real(r_kind), dimension(:,:), allocatable :: rho
      real(r_kind), dimension(:,:), allocatable :: theta
      real(r_kind), dimension(:,:), allocatable :: u
      real(r_kind), dimension(:,:), allocatable :: w
      real(r_kind), dimension(:,:), allocatable :: q
      
      real(r_kind), dimension(:,:), allocatable :: r
      real(r_kind), dimension(:,:), allocatable :: exner
      real(r_kind), dimension(:,:), allocatable :: p
      real(r_kind), dimension(:,:), allocatable :: T
      
      real(r_kind), dimension(:,:), allocatable :: dexner
      
      real(r_kind) theta_bar
      real(r_kind) dtheta
      real(r_kind) xr
      real(r_kind) zr
      real(r_kind) x0
      real(r_kind) z0
      
      integer i,k,iVar
      
      allocate(rho  (ics:ice,kcs:kce))
      allocate(theta(ics:ice,kcs:kce))
      allocate(u    (ics:ice,kcs:kce))
      allocate(w    (ics:ice,kcs:kce))
      allocate(q    (ics:ice,kcs:kce))
      
      allocate(r    (ics:ice,kcs:kce))
      allocate(exner(ics:ice,kcs:kce))
      allocate(p    (ics:ice,kcs:kce))
      allocate(T    (ics:ice,kcs:kce))
      
      allocate(dexner(ics:ice,kcs:kce))
      
      zs    = 0.
      dzsdx = 0.
      
      zs_ext    = 0.
      dzsdx_ext = 0.
      
      call init_vertical_coordinate
      
      theta_bar = 300.
      dtheta    = -15.
      xr        = 4000.
      zr        = 2000.
      x0        = 0.
      z0        = 3000.
           
      do k = kcs,kce
        do i = ics,ice
          r     (i,k) = sqrt( ( ( x(i,k) - x0 ) / xr )**2 + ( ( z(i,k) - z0 ) / zr )**2 )
          theta (i,k) = theta_bar
          dexner(i,k) = -gravity / ( Cpd * theta(i,k) )
        enddo
      enddo
      
      do i = ics,ice
        do k = kcs,kce
          exner(i,k) = 1. + dexner(i,k) * z(i,k)
          p    (i,k) = p0 * exner(i,k)**(Cpd/Rd)
          T    (i,k) = exner(i,k) * theta(i,k)
          rho  (i,k) = p(i,k) / ( Rd * T(i,k) )
          u    (i,k) = 0.
          w    (i,k) = 0.
          q    (i,k) = 0.
        enddo
      enddo
      
      print*,'Set reference fields'
      do k = kcs,kce
        q_ref(1,ids:ide,k) = sum(sqrtG(ids:ide,k)*rho(ids:ide,k)) / nx
      enddo
      
      q_ref(2,:,:) = 0.
      q_ref(3,:,:) = 0.
      q_ref(4,:,:) = 300.
      q_ref(5,:,:) = 0.
      
      do iVar = 2,nVar
        q_ref(iVar,:,:) = q_ref(iVar,:,:) * q_ref(1,:,:)
      enddo
      
      ref%q = q_ref
      
      ! Set theta perturbation
      where(r<=1.)theta = theta_bar + dtheta * ( cos(pi*r) + 1. )/2.
      !where(r<=1.)theta = theta_bar + dtheta * cos(pi*r/2)**2
      
      viscosity_coef = 75
      
      stat%q(1,:,:) = sqrtG * rho
      stat%q(2,:,:) = sqrtG * rho * u
      stat%q(3,:,:) = sqrtG * rho * w
      stat%q(4,:,:) = sqrtG * rho * theta
      stat%q(5,:,:) = sqrtG * rho * q
      
      print*,'max/min value of sqrtG           ',maxval(sqrtG        ),minval(sqrtG        )
      print*,'max/min value of G13             ',maxval(G13          ),minval(G13          )
      print*,'max/min value of sqrtG*rho       ',maxval(stat%q(1,:,:)),minval(stat%q(1,:,:))
      print*,'max/min value of sqrtG*rho*u     ',maxval(stat%q(2,:,:)),minval(stat%q(2,:,:))
      print*,'max/min value of sqrtG*rho*v     ',maxval(stat%q(3,:,:)),minval(stat%q(3,:,:))
      print*,'max/min value of sqrtG*rho*theta ',maxval(stat%q(4,:,:)),minval(stat%q(4,:,:))
      print*,'max/min value of sqrtG*rho*q     ',maxval(stat%q(5,:,:)),minval(stat%q(5,:,:))
      
    end subroutine density_current
    
end module test_case_mod
    