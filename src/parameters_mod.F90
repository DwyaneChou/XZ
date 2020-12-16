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
  
  integer(i_kind) :: nIntegralSubSteps ! number of integral substeps in temporal integration scheme
  integer(i_kind) :: nsteps            ! total integral steps
  
  ! Model run time control variables
  integer(i_kind) :: total_run_time   ! total run time for this model in seconds, this variable is determined by run_days, run_hours ...
  integer(i_kind) :: total_run_steps  ! total run steps for this model in seconds, this variable is determined by total_run_time and dt
  
  integer(i_kind), parameter :: extPts = 3
  
  integer(i_kind) :: ids, ide
  integer(i_kind) :: kds, kde
  integer(i_kind) :: ics, ice
  integer(i_kind) :: kcs, kce
  integer(i_kind) :: ies, iee
  integer(i_kind) :: kes, kee
  
  integer(i_kind) :: nx_ext
  integer(i_kind) :: nz_ext
  
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
    
    ! Read namelist
    call readNamelist
    
    ! Calculate total run time in seconds
    total_run_time  = run_days * 86400 + run_hours * 3600 + run_minutes * 60 + run_seconds
    
    ! Calculate total run steps
    total_run_steps = ceiling(total_run_time/dt)
    
    ! Setting the number of substeps in temporal integration scheme
    if(trim(adjustl(integral_scheme)) == 'RK3')then
      nIntegralSubSteps = 3
    elseif(trim(adjustl(integral_scheme)) == 'RK3_TVD')then
      nIntegralSubSteps = 3
    elseif(trim(adjustl(integral_scheme)) == 'RK4')then
      nIntegralSubSteps = 4
    elseif(trim(adjustl(integral_scheme)) == 'SSPRK')then
      nIntegralSubSteps = 4
    elseif(trim(adjustl(integral_scheme)) == 'PC2')then
      nIntegralSubSteps = 3
    elseif(trim(adjustl(integral_scheme)) == 'IRK2')then
      nIntegralSubSteps = 2
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
    elseif(case_num==2)then
      print*,'Schar Mountain case is selected'
      print*,'Reset x_min to -25 km'
      print*,'Reset x_max to  25 km'
      print*,'Reset z_min to  0  km'
      print*,'Reset z_max to  21 km'
      x_min = -25000.
      x_max = 25000.
      z_min = 0.
      z_max = 21000.
    elseif(case_num==3)then
      print*,'Density Current case is selected'
      print*,'Reset x_min to -26.5 km'
      print*,'Reset x_max to  26.5 km'
      print*,'Reset z_min to  0  km'
      print*,'Reset z_max to  6.4 km'
      x_min = -26500.
      x_max = 26500.
      z_min = 0.
      z_max = 6400.
    else
      stop 'Unknow case_num'
    endif
    
    ids = 1
    ide = ( x_max - x_min ) / dx
    kds = 1
    kde = nz
    
    ics = ids - extPts
    ice = ide + extPts
    kcs = kds - extPts
    kce = kde + extPts
      
    nx = ide - ids + 1
    
    nx_ext = ice - ics + 1
    nz_ext = kce - kcs + 1
    
    ies = ics * 2 - 1
    iee = ice * 2 + 1
    kes = kcs * 2 - 1
    kee = kce * 2 + 1
    
    nsteps = total_run_time / dt
    
  end subroutine initParameters
  
end module parameters_mod
    