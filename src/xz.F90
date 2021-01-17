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
      
      do it = 1,nsteps
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
        
        if( mod( it, output_interval )==0 .and. ( it >= output_interval ) )then
          total_mass = sum(stat(new)%q(1,ids:ide,kds:kde))!calc_total_mass     (stat(new))
          
          MCR = ( total_mass - total_mass0 ) / total_mass0
          
          print*,'MCR = ',MCR
            
          output_idx = output_idx + 1
          call history_write_stat(stat(new),output_idx)
          print*,'output index/total',output_idx-1,'/',total_output_num
        endif
        
        call copyStat(stat(old),stat(new))
      enddo
      
      ! Timing end
      call SYSTEM_CLOCK(timeEnd)
      
      print*,'It took ',dble(timeEnd-timeStart)/1000.0,' seconds to run this program'
    end program xz
