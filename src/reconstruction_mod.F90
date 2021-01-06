    module reconstruction_mod
      use constants_mod
      use parameters_mod
      implicit none
      
      integer(i_kind) :: nRecCells
      
      integer(i_kind) :: nRecTerms
      
      real   (r_kind), dimension(:,:,:), allocatable :: polyCoordCoef
      
    contains
    
    subroutine init_reconstruction
      
      nRecCells = stencil_width**2
      nRecTerms = (recPolyDegree+1) * (recPolyDegree+2) / 2
      
      allocate( polyCoordCoef(nRecTerms,ics:ice,kcs:kce) )
      
    end subroutine init_reconstruction
      
    end module reconstruction_mod
    