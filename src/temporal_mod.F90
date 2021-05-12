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
      
      real(r_kind),dimension(4),parameter :: RK4_weight = [1./6.,1./3.,1./3.,1./6.]
      
      !call copyStat(stat(k1),stat_old)
      
      call spatial_operator (stat_old, tend(k1))
      call update_stat      (stat(k2), stat_old, tend(k1), 0.5_r_kind * dt)
      
      call spatial_operator (stat(k2), tend(k2))
      call update_stat      (stat(k3), stat_old, tend(k2), 0.5_r_kind * dt)
      
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
      
      !call history_write_stat(stat(k2),2)
      !stop 'RK3_TVD'
      
      call spatial_operator (stat(k2), tend(k2))
      call update_stat_RK3_TVD_1(stat(k3), stat_old, stat(k2), tend(k2))
      
      !call history_write_stat(stat(k3),2)
      !stop 'RK3_TVD'
      
      call spatial_operator (stat(k3), tend(k3))
      call update_stat_RK3_TVD_2(stat_new, stat_old, stat(k3), tend(k3))
      
      !call history_write_stat(stat_new,2)
      !stop 'RK3_TVD'

    end subroutine RK3_TVD
    
    subroutine IRK2(stat_new,stat_old)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(inout) :: stat_old
      
      integer            iter
      real(r_kind) :: residual_old
      real(r_kind) :: residual_new
      real(r_kind) :: dresidual
      
      call spatial_operator(stat_old, tend(k1))
      
      call update_stat(stat(k1), stat_old, tend(k1), 0.5 * dt)
      
      iter         = 0
      residual_new = Inf
      residual_old = Inf
      dresidual    = Inf
      
      do while( residual_old>=IRK_residual .and. dresidual>=1.e-15 )
        iter = iter + 1
          
        call copyTend(tend(0), tend(k1))
        
        call spatial_operator(stat(k1), tend(k1))
        
        call update_stat(stat(k1), stat_old, tend(k1), 0.5 * dt)
        
        call calc_residual(residual_new, tend(0), tend(k1))
        
        dresidual = abs( ( residual_new - residual_old ) / residual_old )
        
        print*,'residual error in IRK2 = ',residual_new,'after ',iter,'loop'!,', dresidual =',dresidual
        
        residual_old = residual_new
        
        if(any(isnan(stat(k1)%q(1,ids:ide,kds:kde))))then
          call history_write_stat(stat(k1),2)
          
          print*,'Nan appears at i k: ',maxloc(stat(k1)%q(1,ids:ide,kds:kde),isnan(stat(k1)%q(1,ids:ide,kds:kde)))
          print*,'after ',iter,' iterations'
          stop 'Model blow up, maybe smaller dt would help'
        endif
      enddo
      
      call update_stat (stat_new, stat_old, tend(k1), dt)
      
    end subroutine IRK2
    
    !subroutine IRK2(stat_new,stat_old)
    !  type(stat_field), intent(inout) :: stat_new
    !  type(stat_field), intent(inout) :: stat_old
    !  
    !  real(r_kind),parameter :: b1 = 0.5
    !  real(r_kind),parameter :: b2 = 0.5
    !  
    !  integer            iter
    !  real(r_kind) :: residual_old
    !  real(r_kind) :: residual_new
    !  real(r_kind) :: dresidual
    !  
    !  call spatial_operator (stat_old, tend(0))
    !  
    !  call update_stat_IRK2_1(stat(k1), stat_old, tend(0), tend(0))
    !  
    !  call spatial_operator (stat(k1), tend(k1))
    !  
    !  call update_stat_IRK2_2(stat(k2), stat_old, tend(k1), tend(0))
    !  
    !  call spatial_operator (stat(k2), tend(k2))
    !  
    !  iter         = 0
    !  residual_new = Inf
    !  residual_old = Inf
    !  dresidual    = Inf
    !  
    !  do while( residual_old>=IRK_residual) !.and. dresidual>=1.e-15 )
    !    iter = iter + 1
    !      
    !    call copyTend(tend(0), tend(k2))
    !    
    !    call update_stat_IRK2_1(stat(k1), stat_old, tend(k1), tend(k2))
    !    
    !    call spatial_operator (stat(k1), tend(k1))
    !    
    !    call update_stat_IRK2_2(stat(k2), stat_old, tend(k1), tend(k2))
    !    
    !    call spatial_operator (stat(k2), tend(k2))
    !    
    !    call calc_residual(residual_new, tend(0), tend(k2))
    !    
    !    dresidual = abs( ( residual_new - residual_old ) / residual_old )
    !    
    !    print*,'residual error in IRK2 = ',residual_new,'after ',iter,'loop'!,', dresidual =',dresidual
    !    
    !    residual_old = residual_new
    !  enddo
    !  
    !  tend(new)%q  = b1 * tend(k1)%q  + b2 * tend(k2)%q
    !  
    !  call update_stat (stat_new, stat_old, tend(new), dt)
    !  
    !end subroutine IRK2
    
    subroutine SSPRK(stat_new,stat_old)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(inout) :: stat_old
    
      call spatial_operator (stat_old, tend(k1))
      call update_stat      (stat(k1), stat_old, tend(k1), 0.5_r_kind * dt)
      
      call spatial_operator (stat(k1), tend(k2))
      call update_stat      (stat(k2), stat_old, tend(k2), 0.5_r_kind * dt)
      
      call spatial_operator (stat(k2), tend(k3))
      call update_stat_SSPRK(stat(k3), stat_old, stat(k2), tend(k3))
      
      call spatial_operator (stat(k3), tend(k4))
      call update_stat      (stat_new, stat(k3), tend(k4), 0.5_r_kind * dt)
      
    end subroutine SSPRK
    
    subroutine PC2(stat_new,stat_old)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(inout) :: stat_old
      
      call spatial_operator (stat_old, tend(k1))
      call update_stat      (stat(k2), stat_old, tend(k1), 0.5_r_kind * dt)
      
      call spatial_operator (stat(k2), tend(k2))
      call update_stat      (stat(k3), stat_old, tend(k2), 0.5_r_kind * dt)
      
      call spatial_operator(stat(k3), tend(k3))
      
      call update_stat (stat_new, stat_old, tend(k3), dt)
    end subroutine PC2
    
    subroutine update_stat(stat_new, stat_old, tend, inc_t)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old
      type(tend_field), intent(in   ) :: tend
      real(r_kind)    , intent(in   ) :: inc_t
      
      integer(i_kind) :: iVar,i,k

      !$OMP PARALLEL DO PRIVATE(i,iVar) COLLAPSE(3)
      do k = kds,kde
        do i = ids,ide
          do iVar = 1,nVar
            stat_new%q(iVar,i,k) = stat_old%q(iVar,i,k) + inc_t * tend%q(iVar,i,k)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO

    end subroutine update_stat
    
    subroutine update_stat_RK3_TVD_1(stat_new, stat_old,stat1, tend)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old
      type(stat_field), intent(in   ) :: stat1
      type(tend_field), intent(in   ) :: tend

      integer(i_kind) :: iVar,i,k
      
      !$OMP PARALLEL DO PRIVATE(i,iVar) COLLAPSE(3)
      do k = kds,kde
        do i = ids,ide
          do iVar = 1,nVar
            stat_new%q(iVar,i,k) = 0.75 * stat_old%q(iVar,i,k) + 0.25 * stat1%q(iVar,i,k) + 0.25 * dt * tend%q(iVar,i,k)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine update_stat_RK3_TVD_1
    
    subroutine update_stat_RK3_TVD_2(stat_new, stat_old,stat2, tend)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old
      type(stat_field), intent(in   ) :: stat2
      type(tend_field), intent(in   ) :: tend
      
      real(r_kind),dimension(3),parameter :: weight = [1./3., 2./3., 2./3.]
      
      integer(i_kind) :: iVar,i,k
      
      !$OMP PARALLEL DO PRIVATE(i,iVar) COLLAPSE(3)
      do k = kds,kde
        do i = ids,ide
          do iVar = 1,nVar
            stat_new%q(iVar,i,k) = weight(1) * stat_old%q(iVar,i,k) + weight(2) * stat2%q(iVar,i,k) + weight(3) * dt * tend%q(iVar,i,k)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine update_stat_RK3_TVD_2
    
    subroutine update_stat_IRK2_1(stat_new, stat_old, tend1, tend2)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old
      type(tend_field), intent(in   ) :: tend1
      type(tend_field), intent(in   ) :: tend2
      
      real(r_kind),parameter :: a11 = 0.25
      real(r_kind),parameter :: a12 = 0.25 - sqrt(3.)/6.
      
      stat_new%q  = stat_old%q  + dt * ( a11 * tend1%q  + a12 * tend2%q  )
      
    end subroutine update_stat_IRK2_1
    
    subroutine update_stat_IRK2_2(stat_new, stat_old, tend1, tend2)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old
      type(tend_field), intent(in   ) :: tend1
      type(tend_field), intent(in   ) :: tend2
      
      real(r_kind),parameter :: a21 = 0.25 + sqrt(3.)/6.
      real(r_kind),parameter :: a22 = 0.25
      
      stat_new%q  = stat_old%q  + dt * ( a21 * tend1%q  + a22 * tend2%q  )
      
    end subroutine update_stat_IRK2_2
    
    subroutine update_stat_SSPRK(stat_new, stat_old, stat2, tend3)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old
      type(stat_field), intent(in   ) :: stat2
      type(tend_field), intent(in   ) :: tend3
      
      real(r_kind),dimension(3),parameter :: coef = [2./3.,1./3.,1./6.]
      
      stat_new%q = coef(1) * stat_old%q + coef(2) * stat2%q + coef(3) * dt * tend3%q
      
    end subroutine update_stat_SSPRK
    
    subroutine calc_residual(residual, tend_old, tend_new)
      real(r_kind)    , intent(out) :: residual
      type(tend_field), intent(in ) :: tend_old
      type(tend_field), intent(in ) :: tend_new
    
      residual = maxval(abs(tend_new%q(:,ids:ide,kds:kde) - tend_old%q(:,ids:ide,kds:kde)))
    end subroutine calc_residual
    
end module temporal_mod
