module test_case_mod
  use constants_mod
  use parameters_mod
  use mesh_mod
  use stat_mod
  use tend_mod
  use linear_integration_mod
  use spatial_operators_mod
  implicit none
  
    contains
    subroutine init_test_case
      integer(i_kind) iT
      integer(i_kind) i,k
      
      iT = 0
      
      if(case_num==1)then
        call case1(stat(iT))
      else
        stop 'Unknown case_num'
      endif
      
    end subroutine init_test_case
    
    subroutine case1(stat)
      type(stat_field),intent(inout) :: stat
      
      real(r_kind), dimension(:,:), allocatable :: rho
      real(r_kind), dimension(:,:), allocatable :: u
      real(r_kind), dimension(:,:), allocatable :: w
      
      real(r_kind), parameter :: x0 = -50000.
      real(r_kind), parameter :: z0 = 9000.
      
      real(r_kind), parameter :: z1 = 4000.
      real(r_kind), parameter :: z2 = 5000.
      
      real(r_kind), parameter :: u0 = 10.
      
      real(r_kind), parameter :: a0      = 25000.
      real(r_kind), parameter :: h0      = 3000.
      real(r_kind), parameter :: lambda0 = 8000
      real(r_kind), parameter :: rho0    = 1.
      real(r_kind), parameter :: Ax      = 25000.
      real(r_kind), parameter :: Az      = 3000.
      
      real(r_kind) r
      
      integer i,k,iVar
      
      allocate(rho  (ics:ice,kcs:kce))
      allocate(u    (ics:ice,kcs:kce))
      allocate(w    (ics:ice,kcs:kce))
      
      !! flat
      !zs    = 0.
      !dzsdx = 0.
      !
      !zs_ext    = 0.
      !dzsdx_ext = 0.
      
      !mountain
      zs    = h0 * cos( pi * x / lambda0 )**2 * cos( pi * x / a0 / 2. )**2
      dzsdx = -((h0*pi*Cos((pi*x)/(2.*a0))*Cos((pi*x)/lambda0)*(lambda0*Cos((pi*x)/lambda0)*Sin((pi*x)/(2.*a0)) + 2.*a0*Cos((pi*x)/(2.*a0))*Sin((pi*x)/lambda0)))/(a0*lambda0))
      
      zs_ext    = h0 * cos( pi * x_ext / lambda0 )**2 * cos( pi * x_ext / a0 / 2. )**2
      dzsdx_ext = -((h0*pi*Cos((pi*x_ext)/(2.*a0))*Cos((pi*x_ext)/lambda0)*(lambda0*Cos((pi*x_ext)/lambda0)*Sin((pi*x_ext)/(2.*a0)) + 2.*a0*Cos((pi*x_ext)/(2.*a0))*Sin((pi*x_ext)/lambda0)))/(a0*lambda0))
      
      where( abs(x) >= a0 )
        zs    = 0
        dzsdx = 0
      endwhere
      
      where( abs(x_ext) >= a0 )
        zs_ext    = 0
        dzsdx_ext = 0
      endwhere
      
      call init_vertical_coordinate
      
      do k = kcs,kce
        do i = ics,ice
          r = sqrt( ( ( X(i,k) - X0 ) / Ax )**2 + ( ( Z(i,k) - Z0 ) / Az )**2 )
          
          if(r<=1)then
            rho(i,k) = rho0 * cos( pi * r / 2. )**2
          else
            rho(i,k) = 0
          endif
          
          if(z(i,k)>z2)then
            u(i,k) = 1
          elseif(z(i,k)<z1)then
            u(i,k) = 0
          else
            u(i,k) = sin( pi/2. * ( z(i,k) - z1 ) / ( z2 - z1 ) )**2
          endif
          u(i,k) = u0 * u(i,k)
          
          !rho(i,k) = sin(x(i,k)/(x_max-x_min)*40.*pi)*sin(z(i,k)/(z_max-z_min)*30*pi) + 1.
          !u(i,k) = u0
        
          w(i,k) = 0
        enddo
      enddo
      
      stat%q(1,:,:) = sqrtG * rho
      stat%q(2,:,:) = u
      stat%q(3,:,:) = w
      
      print*,'max/min value of sqrtG           ',maxval(sqrtG        ),minval(sqrtG        )
      print*,'max/min value of G13             ',maxval(G13          ),minval(G13          )
      print*,'max/min value of sqrtG*rho       ',maxval(stat%q(1,:,:)),minval(stat%q(1,:,:))
      print*,'max/min value of u               ',maxval(stat%q(2,:,:)),minval(stat%q(2,:,:))
      print*,'max/min value of v               ',maxval(stat%q(3,:,:)),minval(stat%q(3,:,:))
      
    end subroutine case1
    
end module test_case_mod
    