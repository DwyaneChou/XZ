module test_case_mod
  use constants_mod
  use parameters_mod
  use mesh_mod
  use stat_mod
  use tend_mod
  use linear_integration_mod
  implicit none
  
  real(r_kind), dimension(:,:,:), allocatable :: q_ref      ! reference q
  real(r_kind), dimension(:,:,:), allocatable :: q_init_ext ! extended intial q
  
    contains
    subroutine init_test_case
      integer iT
      
      allocate(q_ref     (nVar,ics:ice,kcs:kce))
      allocate(q_init_ext(nVar,ics:ice,kcs:kce))
      
      iT = 0
      
      if(case_num==1)then
        call thermal_bubble(stat(iT))
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
      
      do i = ics,ice
        call spline2_integration(nz_ext-1,z(i,:),dexner(i,:),0.,0.,nz_ext,z(i,:),dexner(i,:))
        !call spline4_integration(nz_ext-1,z(i,:),dexner(i,:)      ,nz_ext,z(i,:),dexner(i,:))
        
        do k = kds,kce
          exner(i,k) = 1. + sum(dexner(i,kds+1:k))
        enddo
        
        do k = kcs,kds-1
          exner(i,k) = 1. - sum(dexner(i,kds:k+1:-1))
        enddo
        
        do k = kcs,kce
          p    (i,k) = p0 * exner(i,k)**(Cpd/Rd)
          T    (i,k) = exner(i,k) * theta(i,k)
          rho  (i,k) = p(i,k) / Rd / T(i,k)
          u    (i,k) = 0.
          w    (i,k) = 0.
          q    (i,k) = 0.
        enddo
      enddo
      
      !do k = kcs,kce
      !  do i = ics,ice
      !    theta(i,k) = theta_bar + dtheta * max( 0., 1. - r(i,k) / R_bubble )
      !  enddo
      !enddo
      
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
      
      print*,'Reset reference fields'
      do k = kcs,kce
        !q_ref(1,:,k) = stat%q(1,:,k)
        q_ref(1,:,k) = sum(sqrtG(:,k)*rho(:,k)) / nx_ext
      enddo
      
      q_ref(2,:,:) = 0.
      q_ref(3,:,:) = 0.
      q_ref(4,:,:) = 300.
      q_ref(5,:,:) = 0.
      
      do iVar = 2,nVar
        q_ref(iVar,:,:) = q_ref(iVar,:,:) * q_ref(1,:,:)
      enddo
      
    end subroutine thermal_bubble
    
end module test_case_mod
    