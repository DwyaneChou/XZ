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
      
      character(8 ) :: date_timer
      character(10) :: timeStart,timeEnd
      
      character(20) :: procedure_name
      character(10) :: temporal_timer1,temporal_timer2
      character(10) :: nan_finder_timer1,nan_finder_timer2
      character(10) :: io_timer1,io_timer2
      
      real(r_kind) :: time1, time2
      
      integer :: output_interval
      
      !Timing start
      call DATE_AND_TIME(date_timer,timeStart)
      
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
      
      if( abs( history_interval / dt - output_interval ) > tolerance )then
        print*,'history_interval divides dt must be integer'
        print*,abs( history_interval / dt - output_interval )
        stop
      endif
      
      call history_init
      call history_write_stat(stat(old),output_idx)
      print*,'output index/total',output_idx-1,'/',total_output_num
            
      total_mass0 = sum(stat(old)%q(1,ids:ide,kds:kde))
      
      do it = 1,nsteps
        call DATE_AND_TIME(date_timer,temporal_timer1)
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
        call DATE_AND_TIME(date_timer,temporal_timer2)
        
        call DATE_AND_TIME(date_timer,nan_finder_timer1)
        if(any(isnan(stat(new)%q(1,ids:ide,kds:kde))))then
          print*,'Nan appears at i k: ',maxloc(stat(new)%q(1,ids:ide,kds:kde),isnan(stat(new)%q(1,ids:ide,kds:kde)))
          print*,'after ',it,' steps'
          stop 'Model blow up, maybe smaller dt would help'
        endif
        call DATE_AND_TIME(date_timer,nan_finder_timer2)
        
        if( mod( it, output_interval )==0 .and. ( it >= output_interval ) )then
          total_mass = sum(stat(new)%q(1,ids:ide,kds:kde))!calc_total_mass     (stat(new))
          
          MCR = ( total_mass - total_mass0 ) / total_mass0
          
          print*,'MCR = ',MCR
            
          output_idx = output_idx + 1
          call DATE_AND_TIME(date_timer,io_timer1)
          call history_write_stat(stat(new),output_idx)
          call DATE_AND_TIME(date_timer,io_timer2)
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
      call DATE_AND_TIME(date_timer,timeEnd)
      
      read(timeStart,*)time1
      read(timeEnd  ,*)time2
      print*,'It took ',time2-time1,' seconds to run this program'
      
    contains
      subroutine procedure_timer(procedure_name,t1,t2,step_num)
        character(20) :: procedure_name
        character(10) :: t1
        character(10) :: t2
        integer  (i4) :: step_num
        
        real(r_kind) :: time1, time2
        
        read(t1,*)time1
        read(t2,*)time2
      
        print*,procedure_name//' elapse time', time2 - time1, 'on step ', step_num
        
      end subroutine procedure_timer
    
    end program xz
