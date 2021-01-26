    program xz
      use constants_mod
      use parameters_mod
      use mesh_mod
      use stat_mod
      use tend_mod
      use test_case_mod
      use spatial_operators_mod
      use temporal_mod
      use io_mod
      implicit none
      
      integer :: it
      
      integer :: old = 0
      integer :: new = 1
      
      real(r_kind) :: total_mass0
      real(r_kind) :: total_mass
      real(r_kind) :: MCR
      
      integer :: output_idx, total_output_num
      
      real(r_kind), parameter :: timer_coef = 1.e9
      
      integer(i4) :: timeStart,timeEnd
      
      character(20) :: procedure_name
      integer  (i8) :: temporal_timer1,temporal_timer2
      integer  (i8) :: nan_finder_timer1,nan_finder_timer2
      integer  (i8) :: io_timer1,io_timer2
      
      integer :: output_interval
      
      !Timing start
      call SYSTEM_CLOCK(timeStart)
      
      call initParameters
      call init_reconstruction
      call init_mesh
      call init_stat
      call init_tend
      call init_test_case
      call init_spatial_operator
      
      output_idx       = 1
      total_output_num = nsteps * dt / history_interval
      output_interval  = nint( history_interval / dt )
      
      if( abs( history_interval / dt - int(history_interval / dt) ) > 1.e-14 )then
        print*,'history_interval divides dt must be integer'
        print*,abs( history_interval / dt - int(history_interval / dt) )
        stop
      endif
      
      call history_init
      call history_write_stat(stat(old),output_idx)
      print*,'output index/total',output_idx-1,'/',total_output_num
            
      total_mass0 = sum(stat(old)%q(1,ids:ide,kds:kde))
      
      do it = 1,nsteps
        call system_clock(temporal_timer1)
        if(trim(adjustl(integral_scheme)) == 'RK3_TVD')then
          call RK3_TVD(stat(new),stat(old))
        elseif(trim(adjustl(integral_scheme)) == 'RK4')then
          call RK4(stat(new),stat(old))
        elseif(trim(adjustl(integral_scheme)) == 'IRK2')then
          call IRK2(stat(new),stat(old))
        elseif(trim(adjustl(integral_scheme)) == 'PC2')then
          call PC2(stat(new),stat(old))
        elseif(trim(adjustl(integral_scheme)) == 'SSPRK')then
          call SSPRK(stat(new),stat(old))
        endif
        call system_clock(temporal_timer2)
        
        call system_clock(nan_finder_timer1)
        if(any(isnan(stat(new)%q(1,ids:ide,kds:kde))))then
          print*,'Nan appears at i k: ',maxloc(stat(new)%q(1,ids:ide,kds:kde),isnan(stat(new)%q(1,ids:ide,kds:kde)))
          print*,'after ',it,' steps'
          stop 'Model blow up, maybe smaller dt would help'
        endif
        call system_clock(nan_finder_timer2)
        
        if( mod( it, output_interval )==0 .and. ( it >= output_interval ) )then
          total_mass = sum(stat(new)%q(1,ids:ide,kds:kde))!calc_total_mass     (stat(new))
          
          MCR = ( total_mass - total_mass0 ) / total_mass0
          
          print*,'MCR = ',MCR
            
          output_idx = output_idx + 1
          call system_clock(io_timer1)
          call history_write_stat(stat(new),output_idx)
          call system_clock(io_timer2)
          print*,'output index/total',output_idx-1,'/',total_output_num
          procedure_name = 'IO'
          call procedure_timer(procedure_name,io_timer1,io_timer2,it)
        endif
        
        procedure_name = 'Temporal integration'
        call procedure_timer(procedure_name,temporal_timer1,temporal_timer2,it)
        procedure_name = 'NaN finder'
        call procedure_timer(procedure_name,nan_finder_timer1,nan_finder_timer2,it)
        print*,''
        
        call copyStat(stat(old),stat(new))
      enddo
      
      ! Timing end
      call SYSTEM_CLOCK(timeEnd)
      
      print*,'It took ',real(timeEnd-timeStart,r_kind)/timer_coef,' seconds to run this program'
      
    contains
      subroutine procedure_timer(procedure_name,t1,t2,step_num)
        character(20) :: procedure_name
        integer  (i8) :: t1
        integer  (i8) :: t2
        integer  (i4) :: step_num
        
        print*,procedure_name//' elapse time', real( t2 - t1, r_kind ) / timer_coef, 'on step ', step_num
        
      end subroutine procedure_timer
    
    end program xz
