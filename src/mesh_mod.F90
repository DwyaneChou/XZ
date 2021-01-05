module mesh_mod
  use constants_mod
  use parameters_mod
  use math_mod
  implicit none
  
  integer(i_kind) :: nCells
  integer(i_kind) :: nEdges
  integer(i_kind) :: nEdgesOnCell
  integer(i_kind) :: nPointsOnEdge
  integer(i_kind) :: nQuadraturePointsOnCell
  
  real(r_kind), dimension(:,:), allocatable :: xLL
  real(r_kind), dimension(:,:), allocatable :: xLR
  real(r_kind), dimension(:,:), allocatable :: xUL
  real(r_kind), dimension(:,:), allocatable :: xUR
  
  real(r_kind), dimension(:,:), allocatable :: zLL
  real(r_kind), dimension(:,:), allocatable :: zLR
  real(r_kind), dimension(:,:), allocatable :: zUL
  real(r_kind), dimension(:,:), allocatable :: zUR
  
  real(r_kind), dimension(:,:,:), allocatable :: xQP ! x coordinate of Quadrature Points
  real(r_kind), dimension(:,:,:), allocatable :: zQP ! z coordinate of Quadrature Points
  
  real(r_kind), dimension(:,:,:), allocatable :: zsQP ! topography on Quadrature Points
end module mesh_mod
