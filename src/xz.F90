    program xz
      use constants_mod
      use parameters_mod
      use mesh_mod
      use stat_mod
      use tend_mod
      use test_case_mod
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
      
      ! Timing end
      call SYSTEM_CLOCK(timeEnd)
      
      print*,'It took ',dble(timeEnd-timeStart)/1000.0,' seconds to run this program'
    end program xz
