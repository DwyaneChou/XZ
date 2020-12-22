module temporal_mod
  use constants_mod
  use parameters_mod
  use mesh_mod
  use stat_mod
  use tend_mod
  use spatial_operators_mod
  use io_mod
  implicit none
  
  private
  
  public RK4, RK3_TVD, IRK2, PC2, SSPRK
  
  integer, parameter :: k1 = -1
  integer, parameter :: k2 = -2
  integer, parameter :: k3 = -3
  integer, parameter :: k4 = -4
  integer, parameter :: new = 1
  integer, parameter :: old = 0
  
    contains

    subroutine RK4(stat_new,stat_old)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(inout) :: stat_old
      
      integer :: iPoint
      
      real,dimension(4),parameter :: RK4_weight = [1./6.,1./3.,1./3.,1./6.]
      
      !call copyStat(stat(k1),stat_old)
      
      call spatial_operator (stat_old, tend(k1))
      call update_stat      (stat(k2), stat_old, tend(k1), 0.5 * dt)
      
      call spatial_operator (stat(k2), tend(k2))
      call update_stat      (stat(k3), stat_old, tend(k2), 0.5 * dt)
      
      call spatial_operator (stat(k3), tend(k3))
      call update_stat      (stat(k4), stat_old, tend(k3), dt)
      
      call spatial_operator(stat(k4), tend(k4))
            
      tend(new)%q = RK4_weight(1) * tend(k1)%q + RK4_weight(2) * tend(k2)%q + RK4_weight(3) * tend(k3)%q + RK4_weight(4) * tend(k4)%q
      
      call update_stat (stat_new, stat_old, tend(new), dt)

    end subroutine RK4
    
    subroutine RK3_TVD(stat_new,stat_old)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(inout) :: stat_old
      
      !call copyStat(stat(k1),stat_old)
      
      call spatial_operator (stat_old, tend(k1))
      call update_stat      (stat(k2), stat_old, tend(k1), dt)
      
      call spatial_operator (stat(k2), tend(k2))
      call update_stat_RK3_TVD_1(stat(k3), stat_old, stat(k2), tend(k2))
      
      call spatial_operator (stat(k3), tend(k3))
      call update_stat_RK3_TVD_2(stat_new, stat_old, stat(k3), tend(k3))

    end subroutine RK3_TVD
    
    subroutine IRK2(stat_new,stat_old)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(inout) :: stat_old
      
      real,parameter :: b1 = 0.5
      real,parameter :: b2 = 0.5
      
      integer            iter
      real :: residual_old
      real :: residual_new
      real :: dresidual
      
      call spatial_operator (stat_old, tend(0))
      
      call update_stat_IRK2_1(stat(k1), stat_old, tend(0), tend(0))
      
      call spatial_operator (stat(k1), tend(k1))
      
      call update_stat_IRK2_2(stat(k2), stat_old, tend(k1), tend(0))
      
      call spatial_operator (stat(k2), tend(k2))
      
      iter         = 0
      residual_new = Inf
      residual_old = Inf
      dresidual    = Inf
      
      do while( residual_old>=IRK_residual) !.and. dresidual>=1.e-15 )
        iter = iter + 1
          
        call copyTend(tend(0), tend(k2))
        
        call update_stat_IRK2_1(stat(k1), stat_old, tend(k1), tend(k2))
        
        call spatial_operator (stat(k1), tend(k1))
        
        call update_stat_IRK2_2(stat(k2), stat_old, tend(k1), tend(k2))
        
        call spatial_operator (stat(k2), tend(k2))
        
        call calc_residual(residual_new, tend(0), tend(k2))
        
        dresidual = abs( ( residual_new - residual_old ) / residual_old )
        
        print*,'residual error in IRK2 = ',residual_new,'after ',iter,'loop'!,', dresidual =',dresidual
        
        residual_old = residual_new
      enddo
      
      tend(new)%q  = b1 * tend(k1)%q  + b2 * tend(k2)%q
      
      call update_stat (stat_new, stat_old, tend(new), dt)
      
    end subroutine IRK2
    
    subroutine SSPRK(stat_new,stat_old)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(inout) :: stat_old
    
      call spatial_operator (stat_old, tend(k1))
      call update_stat      (stat(k1), stat_old, tend(k1), 0.5 * dt)
      
      call spatial_operator (stat(k1), tend(k2))
      call update_stat      (stat(k2), stat_old, tend(k2), 0.5 * dt)
      
      call spatial_operator (stat(k2), tend(k3))
      call update_stat_SSPRK(stat(k3), stat_old, stat(k2), tend(k3))
      
      call spatial_operator (stat(k3), tend(k4))
      call update_stat      (stat_new, stat(k3), tend(k4), 0.5 * dt)
      
    end subroutine SSPRK
    
    subroutine PC2(stat_new,stat_old)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(inout) :: stat_old
      
      !call copyStat(stat(k1),stat_old)
      
      call spatial_operator (stat_old, tend(k1))
      call update_stat      (stat(k2), stat_old, tend(k1), 0.5 * dt)
      
      call spatial_operator (stat(k2), tend(k2))
      call update_stat      (stat(k3), stat_old, tend(k2), 0.5 * dt)
      
      call spatial_operator(stat(k3), tend(k3))
      
      call update_stat (stat_new, stat_old, tend(k3), dt)
    end subroutine PC2
    
    subroutine update_stat(stat_new, stat_old, tend, inc_t)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old
      type(tend_field), intent(in   ) :: tend
      real            , intent(in   ) :: inc_t
      
      integer iVar
      
      print*,'start'
      iVar = 2
      !print*,maxval(abs(stat_old%q(:,1:nx/2,kds:kde)) - abs(stat_old%q(:,nx:nx/2+1:-1,kds:kde))),maxval(stat_old%q(:,ids:ide,kds:kde))
      !print*,maxval(abs(tend%q(:,1:nx/2,kds:kde)) - abs(tend%q(:,nx:nx/2+1:-1,kds:kde))),maxval(tend%q(:,ids:ide,kds:kde))
      print*,maxval(abs(stat_old%q(iVar,1:nx/2,kds:kde) + stat_old%q(iVar,nx:nx/2+1:-1,kds:kde)))
      print*,maxval(abs(tend%q(iVar,1:nx/2,kds:kde) + tend%q(iVar,nx:nx/2+1:-1,kds:kde)))
      stat_new%q = stat_old%q + inc_t * tend%q
      !print*,maxval(abs(stat_new%q(:,1:nx/2,kds:kde)) - abs(stat_new%q(:,nx:nx/2+1:-1,kds:kde))),maxval(stat_new%q(:,ids:ide,kds:kde))
      print*,maxval(abs(stat_new%q(iVar,1:nx/2,kds:kde) + stat_new%q(iVar,nx:nx/2+1:-1,kds:kde)))
      
    end subroutine update_stat
    
    subroutine update_stat_RK3_TVD_1(stat_new, stat_old,stat1, tend)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old
      type(stat_field), intent(in   ) :: stat1
      type(tend_field), intent(in   ) :: tend
      
      stat_new%q = 0.75 * stat_old%q + 0.25 * stat1%q + 0.25 * dt * tend%q
      
    end subroutine update_stat_RK3_TVD_1
    
    subroutine update_stat_RK3_TVD_2(stat_new, stat_old,stat2, tend)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old
      type(stat_field), intent(in   ) :: stat2
      type(tend_field), intent(in   ) :: tend
      
      real,dimension(3),parameter :: weight = [1./3., 2./3., 2./3.]
      
      stat_new%q = weight(1) * stat_old%q + weight(2) * stat2%q + weight(3) * dt * tend%q
    end subroutine update_stat_RK3_TVD_2
    
    subroutine update_stat_IRK2_1(stat_new, stat_old, tend1, tend2)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old
      type(tend_field), intent(in   ) :: tend1
      type(tend_field), intent(in   ) :: tend2
      
      real,parameter :: a11 = 0.25
      real,parameter :: a12 = 0.25 - sqrt(3.)/6.
      
      stat_new%q  = stat_old%q  + dt * ( a11 * tend1%q  + a12 * tend2%q  )
      
    end subroutine update_stat_IRK2_1
    
    subroutine update_stat_IRK2_2(stat_new, stat_old, tend1, tend2)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old
      type(tend_field), intent(in   ) :: tend1
      type(tend_field), intent(in   ) :: tend2
      
      real,parameter :: a21 = 0.25 + sqrt(3.)/6.
      real,parameter :: a22 = 0.25
      
      stat_new%q  = stat_old%q  + dt * ( a21 * tend1%q  + a22 * tend2%q  )
      
    end subroutine update_stat_IRK2_2
    
    subroutine update_stat_SSPRK(stat_new, stat_old, stat2, tend3)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old
      type(stat_field), intent(in   ) :: stat2
      type(tend_field), intent(in   ) :: tend3
      
      real,dimension(3),parameter :: coef = [2./3.,1./3.,1./6.]
      
      stat_new%q = coef(1) * stat_old%q + coef(2) * stat2%q + coef(3) * dt * tend3%q
      
    end subroutine update_stat_SSPRK
    
    subroutine calc_residual(residual, tend_old, tend_new)
      real            , intent(out) :: residual
      type(tend_field), intent(in ) :: tend_old
      type(tend_field), intent(in ) :: tend_new
    
      residual = abs(maxval(tend_new%q - tend_old%q))
    end subroutine calc_residual
    
end module temporal_mod
