module parameters_mod
  use constants_mod
  implicit none
  
  ! Namelist parameters
  ! time_settings
  integer         :: run_days
  integer         :: run_hours
  integer         :: run_minutes
  integer         :: run_seconds
  real   (r_kind) :: dt               ! time step
  real   (r_kind) :: history_interval ! output interval in seconds
  real            :: IRK_residual
  
  character*200 :: integral_scheme
  
  ! Case select
  integer :: case_num
  
  ! Domain
  real   (r_kind) :: dx        !  grid space in x-direction
  integer(i_kind) :: nx        !  grid number in x-direction
  integer(i_kind) :: nz        !  grid number in z-direction
  real   (r_kind) :: x_max     !  max x coordinate
  real   (r_kind) :: x_min     !  min x coordinate
  real   (r_kind) :: z_max     !  max z coordinate
  real   (r_kind) :: z_min     !  min z coordinate
  real   (r_kind) :: m_coef    !  m coefficient in atan vertical distribution
  
  integer(i_kind) :: vertical_distribution ! 1 for even, 2 for atan function
  integer(i_kind) :: vertical_coordinate   ! 1 for Gal-Chen, 2 for Klemp 2011
  
  integer(i_kind) ids,ide
  integer(i_kind) kds,kde
  
  integer(i_kind) :: nIntegralSubSteps ! number of integral substeps in temporal integration scheme
  integer(i_kind) :: nsteps            ! total integral steps
  
  ! Model run time control variables
  integer(i_kind) :: total_run_time   ! total run time for this model in seconds, this variable is determined by run_days, run_hours ...
  integer(i_kind) :: total_run_steps  ! total run steps for this model in seconds, this variable is determined by total_run_time and dt
  
  ! dynamic_opt
  real(r_kind) :: ref_u     ! reference u wind
  real(r_kind) :: ref_w     ! reference w wind
  real(r_kind) :: ref_theta ! reference potential temperature
  real(r_kind) :: ref_q     ! reference specific humidity
  
  namelist /time_settings/ dt               ,&
                           run_days         ,&
                           run_hours        ,&
                           run_minutes      ,&
                           run_seconds      ,&
                           history_interval ,&
                           integral_scheme  ,&
                           IRK_residual
  
  namelist /case_select/   case_num
  namelist /domain/ dx                   ,&
                    nz                   ,&
                    x_max                ,&
                    x_min                ,&
                    z_max                ,&
                    z_min                ,&
                    m_coef               ,&
                    vertical_distribution,&
                    vertical_coordinate
  
  namelist /dynamic_opt/ ref_u    ,&
                         ref_w    ,&
                         ref_theta,&
                         ref_q    
  
  contains
  
  subroutine readNamelist
    
    open(1, file = 'namelist.input',status='old')
    read(1, nml  = time_settings)
    read(1, nml  = case_select  )
    read(1, nml  = domain       )
    close(1)
    
  end subroutine readNamelist
  
  subroutine initParameters
    
    ! Setting default values
    dt    = 300.
    dx    = 2000.
    nz    = 60
    x_max = 100000.
    x_min = -100000.
    z_max = 60000.
    z_min = 0.
    
    vertical_distribution = 1
    vertical_coordinate   = 1
    
    run_days         = 1
    run_hours        = 0
    run_minutes      = 0
    run_seconds      = 0
    history_interval = 360
    integral_scheme  = 'RK4'
    
    ref_u     = 0.
    ref_w     = 0.
    ref_theta = 300.
    ref_q     = 0.
    
    ! Read namelist
    call readNamelist
    
    ! Calculate total run time in seconds
    total_run_time  = run_days * 86400 + run_hours * 3600 + run_minutes * 60 + run_seconds
    
    ! Calculate total run steps
    total_run_steps = ceiling(total_run_time/dt)
    
    ! Setting the number of substeps in temporal integration scheme
    if(trim(adjustl(integral_scheme)) == 'RK3')then
      nIntegralSubSteps = 3
    elseif(trim(adjustl(integral_scheme)) == 'RK4')then
      nIntegralSubSteps = 4
    else
      stop 'Unknown integral scheme, please select from RK3 or RK4 ...'
    endif
    
    if(case_num==1)then
      print*,'Thermal Bubble case is selected'
      print*,'Reset x_min to 0  km'
      print*,'Reset x_max to 20 km'
      print*,'Reset z_min to 0  km'
      print*,'Reset z_max to 10 km'
      x_min = 0.
      x_max = 20. * 1000.
      z_min = 0.
      z_max = 10. * 1000.
    else
      stop 'Unknow case_num'
    endif
    
    ids = 1
    ide = ( x_max - x_min ) / dx + 1
    kds = 1
    kde = nz
    
    nx = ( x_max - x_min ) / dx + 1
    
    nsteps = total_run_time / dt
    
  end subroutine initParameters
  
end module parameters_mod
    