module tend_mod
  use constants_mod
  use parameters_mod
  use mesh_mod
  implicit none
  
  type tend_field
    real, dimension(:,:,:), allocatable :: q
  end type tend_field
  
  type(tend_field), dimension(:), allocatable :: tend ! allocated by n time points, which is used by temporal integration schemes
  
  contains
  
  subroutine init_tend
    integer :: iT
    
    allocate( tend(-nIntegralSubSteps:1) )
    
    do iT = -nIntegralSubSteps, 1
      allocate(tend(iT)%q(nVar,ids:ide,kds:kde))
    enddo
    
  end subroutine init_tend

  subroutine copyTend(tend_out,tend_in)
    type(tend_field),intent(inout) :: tend_out
    type(tend_field),intent(in   ) :: tend_in
  
    tend_out%q = tend_in%q
  end subroutine copyTend
end module tend_mod
    