module test_case_mod
  use constants_mod
  use parameters_mod
  use mesh_mod
  use stat_mod
  use tend_mod
  use linear_integration_mod
  use spatial_operators_mod, only: WENO_limiter, calc_H, calc_pressure, calc_eigenvalue_z
  implicit none
  
  real(r_kind), dimension(:,:,:), allocatable :: q_ref      ! reference q
  
    contains
    subroutine init_test_case
      integer iT
      
      allocate(q_ref(nVar,ics:ice,kcs:kce))
      
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
        !call spline2_integration(nz_ext-1,z(i,:),dexner(i,:),0.,0.,nz_ext,z(i,:),dexner(i,:))
        call spline4_integration(nz_ext-1,z(i,:),dexner(i,:)      ,nz_ext,z(i,:),dexner(i,:))
        
        exner(i,1) = 1.
        
        do k = kds+1,kce
          exner(i,k) = exner(i,1) + sum(dexner(i,kds+1:k))
        enddo
        
        do k = kcs,kds-1
          exner(i,k) = 1. - sum(dexner(i,kds:k+1:-1))
        enddo
        
        do k = kcs,kce
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
        q_ref(1,:,k) = sum(sqrtG(:,k)*rho(:,k)) / nx_ext
      enddo
      
      q_ref(2,:,:) = 0.
      q_ref(3,:,:) = 0.
      q_ref(4,:,:) = 300.
      q_ref(5,:,:) = 0.
      
      do iVar = 2,nVar
        q_ref(iVar,:,:) = q_ref(iVar,:,:) * q_ref(1,:,:)
      enddo
      
      ref%q = q_ref
      
      call hydrostatic(ref%q)
      
      !! Set theta perturbation
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
      
    end subroutine thermal_bubble
    
    subroutine hydrostatic(q)
      real(r_kind), dimension(nVar,ics:ice,kcs:kce), intent(inout) :: q
      
      real(r_kind), dimension(:,:), allocatable :: rho
      real(r_kind), dimension(:,:), allocatable :: rho_theta
      real(r_kind), dimension(:,:), allocatable :: theta
      real(r_kind), dimension(:,:), allocatable :: p
      
      real(r_kind), dimension(:,:,:), allocatable :: q_ext
      
      real(r_kind), dimension(:,:,:), allocatable :: qB    ! Reconstructed q_(i,k-1/2)
      real(r_kind), dimension(:,:,:), allocatable :: qT    ! Reconstructed q_(i,k+1/2)
      
      real(r_kind), dimension(:,:,:), allocatable :: HB    ! Reconstructed H_(i,k-1/2)
      real(r_kind), dimension(:,:,:), allocatable :: HT    ! Reconstructed H_(i,k+1/2)
      
      real(r_kind), dimension(:,:  ), allocatable :: PB    ! Reconstructed P_(i,k-1/2)
      real(r_kind), dimension(:,:  ), allocatable :: PT    ! Reconstructed P_(i,k+1/2)
      
      real(r_kind), dimension(:,:,:), allocatable :: He    ! H on edges of each cell
      
      real(r_kind), dimension(5) :: q_weno
      
      real(r_kind) eigenvalue_z(5,2)
      real(r_kind) maxeigen_z
      
      integer(i_kind) dir
      
      integer(i_kind) ip1,im1
      integer(i_kind) kp1,km1
      integer(i_kind) ip2,im2
      integer(i_kind) kp2,km2
      
      integer iter
      
      integer(i_kind),parameter :: max_iter = 100
      real   (r_kind)           :: residual_error
      
      integer i,k,iVar
      
      
      allocate(rho      (ics:ice,kcs:kce))
      allocate(rho_theta(ics:ice,kcs:kce))
      allocate(theta    (ics:ice,kcs:kce))
      allocate(p        (ics:ice,kcs:kce))
      
      allocate(q_ext(nVar,ics:ice,kcs:kce))
      allocate(qB   (nVar,ids:ide,kds:kde+1))
      allocate(qT   (nVar,ids:ide,kds-1:kde))
      
      allocate(HB   (nVar,ids:ide,kds:kde+1))
      allocate(HT   (nVar,ids:ide,kds-1:kde))
      
      allocate(PB(ids:ide,kds:kde+1))
      allocate(PT(ids:ide,kds-1:kde))
      
      allocate(He(nVar,ids:ide,kds:kde+1))
      
      theta = q(4,:,:) / q(1,:,:)
      
      iter           = 0
      residual_error = 10
      do while(residual_error>1.e-10 .and. iter<=max_iter)
        q_ext = q
        !$OMP PARALLEL DO PRIVATE(i,iVar,q_weno,dir,kp1,km1,kp2,km2)
        do k = kds,kde
          kp1 = k + 1
          km1 = k - 1
          kp2 = k + 2
          km2 = k - 2
          do i = ids,ide
            do iVar = 1,nVar
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
            if(k==kds)qB(3,i,k) = 0
            if(k==kde)qT(3,i,k) = 0
            
            PB(i,k) = calc_pressure(sqrtG(i,k),qB(:,i,k))
            PT(i,k) = calc_pressure(sqrtG(i,k),qT(:,i,k))
            
            HB(:,i,k) = calc_H(sqrtG(i,k),G13(i,k),qB(:,i,k),PB(i,k))
            HT(:,i,k) = calc_H(sqrtG(i,k),G13(i,k),qT(:,i,k),PT(i,k))
          enddo
        enddo
        !$OMP END PARALLEL DO
        HB(3,:,kde+1) = PT(:,kde)
        HT(3,:,kds-1) = PB(:,kds)
        
        ! calc z flux
        !$OMP PARALLEL DO PRIVATE(k,km1,eigenvalue_z,maxeigen_z,iVar)
        do i = ids,ide
          do k = kds,kde+1
            km1 = k - 1
            eigenvalue_z(:,1) = calc_eigenvalue_z(sqrtG(i,km1),G13(i,km1),q_ext(:,i,km1))
            eigenvalue_z(:,2) = calc_eigenvalue_z(sqrtG(i  ,k),G13(i  ,k),q_ext(:,i  ,k))
            
            maxeigen_z = maxval(abs(eigenvalue_z))
            
            He(:,i,k) = 0.5 * ( HB(:,i,k) + HT(:,i,km1) - maxeigen_z * ( qB(:,i,k) - qT(:,i,km1) ) )
            !do iVar = 1,nVar
            !  if(abs(HB(iVar,i,k) + HT(iVar,i,km1))<=1.E-15)then
            !    He(iVar,i,k) = 0
            !  else
            !    He(iVar,i,k) = 0.5 * ( HB(iVar,i,k) + HT(iVar,i,km1) - maxeigen_z * ( qB(iVar,i,k) - qT(iVar,i,km1) ) )
            !  endif
            !enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
        !He(1,ids:ide,kds  ) = 0
        !He(1,ids:ide,kde+1) = 0
        !He(2,ids:ide,kds  ) = sqrtG(ids:ide,kds) * G13(ids:ide,kds) * PB(ids:ide,kds)
        !He(2,ids:ide,kde+1) = sqrtG(ids:ide,kde) * G13(ids:ide,kde) * PT(ids:ide,kde)
        !He(3,ids:ide,kds  ) = PB(ids:ide,kds)
        !He(3,ids:ide,kde+1) = PT(ids:ide,kde)
        !He(4,ids:ide,kds  ) = 0
        !He(4,ids:ide,kde+1) = 0
        !He(5,ids:ide,kds  ) = 0
        !He(5,ids:ide,kde+1) = 0
        
        do i = ids,ide
          do k = kds,kde
            rho      (i,k) = -( He(3,i,k+1) - He(3,i,k) ) / ( deta * sqrtG(i,k) * gravity )
            rho_theta(i,k) = rho(i,k) * theta(i,k)
          enddo
        enddo
        
        !i = ids
        !do k = kds,kde
        !  print*,iter, k, q(1,i,k)/sqrtG(i,k), rho(i,k), rho_theta(i,k), He(3,i,k+1), He(3,i,k)
        !enddo
        
        q(1,ids:ide,kds:kde) = sqrtG(ids:ide,kds:kde) * rho      (ids:ide,kds:kde)
        q(4,ids:ide,kds:kde) = sqrtG(ids:ide,kds:kde) * rho_theta(ids:ide,kds:kde)
        
        residual_error = max( maxval( abs( q(1,ids:ide,kds:kde) - q_ext(1,ids:ide,kds:kde) ) / q_ext(1,ids:ide,kds:kde) ),&
                              maxval( abs( q(4,ids:ide,kds:kde) - q_ext(4,ids:ide,kds:kde) ) / q_ext(4,ids:ide,kds:kde) ) )
        
        iter = iter + 1
        
        print*,'iter, residual_error:',iter,residual_error
      enddo
    end subroutine hydrostatic
    
end module test_case_mod
    