module stat_mod
  use constants_mod
  use parameters_mod
  use mesh_mod
  implicit none
  type stat_field
    real(r_kind), dimension(:,:,:), allocatable :: q
    
    real(r_kind), dimension(:,:,:), allocatable :: vc ! vectors in computatinal space
    real(r_kind), dimension(:,:,:), allocatable :: vp ! vectors in physical space
  end type stat_field
  
  type(stat_field), dimension(:), allocatable :: stat  ! allocated by n time points, which is used by temporal integration schemes
  type(stat_field)                            :: ref   ! reference fields
  
    contains
    subroutine init_stat
      integer :: iT
      
      allocate(stat(-nIntegralSubSteps:1))
      
      do iT = -nIntegralSubSteps, 1
        allocate(stat(iT)%q (nVar,ics:ice,kcs:kce))
        allocate(stat(iT)%vc(2   ,ics:ice,kcs:kce))
        allocate(stat(iT)%vp(2   ,ics:ice,kcs:kce))
      enddo
      
      allocate(ref%q (nVar,ics:ice,kcs:kce))
      allocate(ref%vc(2   ,ics:ice,kcs:kce))
      allocate(ref%vp(2   ,ics:ice,kcs:kce))
      
    end subroutine init_stat

    subroutine copyStat(stat_out,stat_in)
      type(stat_field),intent(out) :: stat_out
      type(stat_field),intent(in ) :: stat_in
    
      stat_out%q  = stat_in%q
      !stat_out%vc = stat_in%vc
      !stat_out%vp = stat_in%vp
    end subroutine copyStat
end module stat_mod
    