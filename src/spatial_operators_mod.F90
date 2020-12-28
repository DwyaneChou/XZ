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
  
  real(r_kind), dimension(:,:,:), allocatable :: Fe    ! F on edges of each cell
  real(r_kind), dimension(:,:,:), allocatable :: He    ! H on edges of each cell
  
  real(r_kind), dimension(:,:,:), allocatable :: src   ! source term
  
  real(r_kind), dimension(  :,:), allocatable :: u
  real(r_kind), dimension(  :,:), allocatable :: uL
  real(r_kind), dimension(  :,:), allocatable :: uR
  real(r_kind), dimension(  :,:), allocatable :: uB
  real(r_kind), dimension(  :,:), allocatable :: uT
  
  real(r_kind), dimension(  :,:), allocatable :: detadt
  real(r_kind), dimension(  :,:), allocatable :: detadtB
  real(r_kind), dimension(  :,:), allocatable :: detadtT
  
    contains
    subroutine init_spatial_operator
      integer(i_kind) dir
      integer(i_kind) i,k,iVar
      
      real(r_kind), dimension(5) :: q_weno
      
      integer(i_kind) ip1,im1
      integer(i_kind) kp1,km1
      integer(i_kind) ip2,im2
      integer(i_kind) kp2,km2
      
      allocate(qC(nVar,ics:ice,kcs:kce))
      
      allocate(qL   (nVar,ics:ice,kcs:kce))
      allocate(qR   (nVar,ics:ice,kcs:kce))
      allocate(qB   (nVar,ics:ice,kcs:kce))
      allocate(qT   (nVar,ics:ice,kcs:kce))
      
      allocate(F    (nVar,ids:ide,kds:kce))
      allocate(H    (nVar,ids:ide,kds:kce))
      
      allocate(FL   (nVar,ics:ice,kcs:kce))
      allocate(FR   (nVar,ics:ice,kcs:kce))
      
      allocate(HB   (nVar,ics:ice,kcs:kce))
      allocate(HT   (nVar,ics:ice,kcs:kce))
      
      allocate(Fe   (nVar,ids:ide+1,kds:kde  ))
      allocate(He   (nVar,ids:ide  ,kds:kde+1))
      
      allocate(src  (nVar,ids:ide,kds:kde))
      
      allocate(u (ics:ice,kcs:kce))
      allocate(uL(ics:ice,kcs:kce))
      allocate(uR(ics:ice,kcs:kce))
      allocate(uB(ics:ice,kcs:kce))
      allocate(uT(ics:ice,kcs:kce))
      
      allocate(detadt (ics:ice,kcs:kce))
      allocate(detadtB(ics:ice,kcs:kce))
      allocate(detadtT(ics:ice,kcs:kce))
      
      ! Attension, error !
      do k = kcs,kce
        do i = ics,ice
          u (i,k) = stat(0)%q(2,i,k)
          uL(i,k) = stat(0)%q(2,i,k)
          uR(i,k) = stat(0)%q(2,i,k)
          uB(i,k) = stat(0)%q(2,i,k)
          uT(i,k) = stat(0)%q(2,i,k)
          
          detadt (i,k) = ( stat(0)%q(3,i,k) + sqrtG (i,k) * G13 (i,k) * u (i,k) ) / sqrtG (i,k)
          detadtB(i,k) = ( stat(0)%q(3,i,k) + sqrtGB(i,k) * G13B(i,k) * uB(i,k) ) / sqrtGB(i,k)
          detadtT(i,k) = ( stat(0)%q(3,i,k) + sqrtGT(i,k) * G13T(i,k) * uT(i,k) ) / sqrtGT(i,k)
        enddo
      enddo
      
      src   = 0
      
      ! Fill out qC
      qC = FillValue
      
    end subroutine init_spatial_operator
    
    subroutine spatial_operator(stat,tend)
      type(stat_field), target, intent(inout) :: stat
      type(tend_field), target, intent(inout) :: tend
      
      real(r_kind), dimension(5) :: q_weno
      
      real(r_kind) :: detadt
      
      real(r_kind) maxeigen_x
      real(r_kind) maxeigen_z
      
      integer(i_kind) dir
      
      integer(i_kind) i,k,iVar
      
      integer(i_kind) ip1,im1
      integer(i_kind) kp1,km1
      integer(i_kind) ip2,im2
      integer(i_kind) kp2,km2
      
      ! Fill periodic boundary
      stat%q(1,ics:ids-1,kds:kde) = stat%q(1,ide-extPts+1:ide,kds:kde)
      stat%q(1,ide+1:ice,kds:kde) = stat%q(1,ids:ids+extPts-1,kds:kde)
      
      ! copy stat
      qC = stat%q
      
      ! Reconstruct X
      !$OMP PARALLEL DO PRIVATE(i,iVar,ip1,im1,ip2,im2,q_weno,dir)
      do k = kds,kde
        do i = ids-1,ide+1
          ip1 = i + 1
          im1 = i - 1
          ip2 = i + 2
          im2 = i - 2
          iVar = 1
          ! x-dir
          q_weno = qC(iVar,im2:ip2,k)
          
          dir = -1
          call WENO_limiter(qL(iVar,i,k),q_weno,dir)
          dir = 1
          call WENO_limiter(qR(iVar,i,k),q_weno,dir)
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! Reconstruct Z
      !$OMP PARALLEL DO PRIVATE(i,iVar,q_weno,dir,kp1,km1,kp2,km2)
      do k = kds,kde
        kp1 = k + 1
        km1 = k - 1
        kp2 = k + 2
        km2 = k - 2
        do i = ids,ide
          iVar = 1
          ! z-dir
          q_weno = qC(iVar,i,km2:kp2)
          
          dir = -1
          call WENO_limiter(qB(iVar,i,k),q_weno,dir)
          dir = 1
          call WENO_limiter(qT(iVar,i,k),q_weno,dir)
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! Calculate functions X
      !$OMP PARALLEL DO PRIVATE(i)
      do k = kds,kde
        do i = ids-1,ide+1
          FL(:,i,k) = qL(1,i,k) * uL(i,k)
          FR(:,i,k) = qR(1,i,k) * uR(i,k)
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! Calculate functions Z
      !$OMP PARALLEL DO PRIVATE(i)
      do k = kds,kde
        do i = ids,ide
          HB(:,i,k) = qB(1,i,k) * detadtB(i,k)
          HT(:,i,k) = qT(1,i,k) * detadtT(i,k)
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      !src = 0
      
      ! calc x flux
      iVar = 1
      !$OMP PARALLEL DO PRIVATE(i,im1,maxeigen_x)
      do k = kds,kde
        do i = ids,ide+1
          im1 = i - 1
          maxeigen_x = sign(1.,uL(i,k))
          Fe(iVar,i,k) = 0.5 * ( FL(iVar,i,k) + FR(iVar,im1,k) - maxeigen_x * ( FL(iVar,i,k) - FR(iVar,im1,k) ) )
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! calc z flux
      iVar = 1
      !$OMP PARALLEL DO PRIVATE(k,km1,maxeigen_z)
      do i = ids,ide
        do k = kds,kde+1
          km1 = k - 1
          maxeigen_z = sign(1.,detadtB(i,k))
          He(iVar,i,k) = 0.5 * ( HB(iVar,i,k) + HT(iVar,i,km1) - maxeigen_z * ( HB(iVar,i,k) - HT(iVar,i,km1) ) )
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! Calculate tend
      iVar = 1
      !$OMP PARALLEL DO PRIVATE(i,ip1,kp1)
      do k = kds,kde
        do i = ids,ide
          ip1 = i + 1
          kp1 = k + 1
          tend%q(iVar,i,k) = - ( ( Fe(iVar,ip1,k) - Fe(iVar,i,k) ) / dx + ( He(iVar,i,kp1) - He(iVar,i,k) ) / deta ) + src(iVar,i,k)
        enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine spatial_operator
    
END MODULE spatial_operators_mod

