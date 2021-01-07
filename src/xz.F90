    program xz
      use constants_mod
      use parameters_mod
      use mesh_mod
      use stat_mod
      use tend_mod
      use test_case_mod
      use io_mod
      use spatial_operators_mod
      implicit none
      
      integer :: it
      
      integer :: old = 0
      integer :: new = 1
      
      real(r_kind) :: total_mass0
      real(r_kind) :: total_mass
      real(r_kind) :: MCR
      
      integer :: output_idx, total_output_num
      
      integer :: timeStart,timeEnd
      
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
      
      if( abs( history_interval / dt - int(history_interval / dt) ) > tolerance )then
        print*,'history_interval divides dt must be integer'
        stop
      endif
      
      call history_init
      call history_write_stat(stat(old),output_idx)
      print*,'output index/total',output_idx-1,'/',total_output_num
      
      total_mass0 = sum(stat(old)%q(1,ids:ide,kds:kde))
      
      ! Timing end
      call SYSTEM_CLOCK(timeEnd)
      
      print*,'It took ',dble(timeEnd-timeStart)/1000.0,' seconds to run this program'
    end program xz
