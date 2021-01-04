module netcdf_mod
  use netcdf
  implicit none
  
  include 'netcdf.inc'

  contains
  
  !************************************************************************
  !!handle_err
  !************************************************************************
  !
  !!ROUTINE:      handle_err
  !!DESCRIPTION:  error handler
  !--------------------------------------------------------------------------
  
  subroutine handle_err(status)
    
    implicit         none
    
    integer          status
    
    if (status .ne. nf_noerr) then
      print *, nf90_strerror(status)
      stop 'Stopped'
    endif
    
  end subroutine handle_err

end module netcdf_mod

