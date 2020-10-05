module io_mod
  use netcdf_mod
  use parameters_mod
  use mesh_mod
  use stat_mod
  use tend_mod
  implicit none
  
  character(500) :: stat_output_file = 'output.nc'
    
    contains
    subroutine history_init
      integer ncid
      integer status
      
      integer nx_dimid
      integer nz_dimid
      
      integer x_id
      integer z_id

      integer sqrtG_id
      integer G13_id
      
      integer rho_id  
      integer u_id    
      integer w_id    
      integer theta_id
      integer q_id    
      integer p_id    
      integer T_id    
      
      character*1 case_num_c
      
      write(case_num_c ,'(i1)')case_num
      
      stat_output_file = 'output_xz_'//case_num_c//'.nc'
      
      status = nf90_create(trim(stat_output_file), NF90_CLOBBER + NF90_64BIT_OFFSET , ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_def_dim(ncid,'nx', nx, nx_dimid)
      status = nf90_def_dim(ncid,'nz', nz, nz_dimid)
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_put_att(ncid,nf90_global, 'case_num' , case_num   )
      status = nf90_put_att(ncid,nf90_global, 'Rd'       , real(Rd,r8))
      status = nf90_put_att(ncid,nf90_global, 'Rv'       , real(Rv,r8))
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_def_var(ncid,'x'    ,NF90_DOUBLE,(/nx_dimid,nz_dimid/),x_id     )
      status = nf90_def_var(ncid,'z'    ,NF90_DOUBLE,(/nx_dimid,nz_dimid/),z_id     )
      
      status = nf90_def_var(ncid,'sqrtG',NF90_DOUBLE,(/nx_dimid,nz_dimid/),sqrtG_id )
      status = nf90_def_var(ncid,'G13'  ,NF90_DOUBLE,(/nx_dimid,nz_dimid/),G13_id   )
      
      status = nf90_def_var(ncid,'rho'  ,NF90_DOUBLE,(/nx_dimid,nz_dimid/),rho_id   )
      status = nf90_def_var(ncid,'u'    ,NF90_DOUBLE,(/nx_dimid,nz_dimid/),u_id     )
      status = nf90_def_var(ncid,'w'    ,NF90_DOUBLE,(/nx_dimid,nz_dimid/),w_id     )
      status = nf90_def_var(ncid,'theta',NF90_DOUBLE,(/nx_dimid,nz_dimid/),theta_id )
      status = nf90_def_var(ncid,'q'    ,NF90_DOUBLE,(/nx_dimid,nz_dimid/),q_id     )
      status = nf90_def_var(ncid,'p'    ,NF90_DOUBLE,(/nx_dimid,nz_dimid/),p_id     )
      status = nf90_def_var(ncid,'T'    ,NF90_DOUBLE,(/nx_dimid,nz_dimid/),T_id     )
      
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_enddef(ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_put_var(ncid,x_id     , real(x    ,r8))
      status = nf90_put_var(ncid,z_id     , real(z    ,r8))
      status = nf90_put_var(ncid,sqrtG_id , real(sqrtG,r8))
      status = nf90_put_var(ncid,G13_id   , real(G13  ,r8))
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_close(ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
    end subroutine history_init
    
    subroutine history_write_stat(stat,time_slot_num)
      type(stat_field), intent(inout) :: stat
      integer         , intent(in   ) :: time_slot_num
      
      integer status
      integer ncid
      
      integer rho_id  
      integer u_id    
      integer w_id    
      integer theta_id
      integer q_id    
      integer p_id    
      integer T_id    
      
      real(r_kind), dimension(nx,nz) :: rho
      real(r_kind), dimension(nx,nz) :: u
      real(r_kind), dimension(nx,nz) :: w
      real(r_kind), dimension(nx,nz) :: theta
      real(r_kind), dimension(nx,nz) :: q
      real(r_kind), dimension(nx,nz) :: T
      real(r_kind), dimension(nx,nz) :: p
      real(r_kind), dimension(nx,nz) :: kappa
      real(r_kind), dimension(nx,nz) :: Ra ! 比气体常数
      real(r_kind), dimension(nx,nz) :: Cp
      real(r_kind), dimension(nx,nz) :: Cv
      
      status = nf90_open(stat_output_file,NF90_WRITE,ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
      !print*,'nf90_inq_varid'
      status = nf90_inq_varid(ncid,'rho'  , rho_id  )
      status = nf90_inq_varid(ncid,'u'    , u_id    )
      status = nf90_inq_varid(ncid,'w'    , w_id    )
      status = nf90_inq_varid(ncid,'theta', theta_id)
      status = nf90_inq_varid(ncid,'q'    , q_id    )
      status = nf90_inq_varid(ncid,'p'    , p_id    )
      status = nf90_inq_varid(ncid,'T'    , T_id    )
      if(status/=nf90_noerr) call handle_err(status)
      
      ! diag air status
      rho   = stat%q(1,:,:)/sqrtG
      u     = stat%q(2,:,:)/stat%q(1,:,:)
      w     = stat%q(3,:,:)/stat%q(1,:,:)
      theta = stat%q(4,:,:)/stat%q(1,:,:)
      q     = stat%q(5,:,:)/stat%q(1,:,:)
      Cp    = Cpd + ( Cpv - Cpd ) * q
      Cv    = Cvd + ( Cvv - Cvd ) * q
      kappa = Cp / Cv
      Ra    = Rd * ( 1. + eq * q )
      p     = p0 * ( rho * Ra * theta / p0 )** kappa
      T     = p / ( rho * Ra )
      
      where( rho == 0. )
        theta = 300.
        p     = 0.
        T     = 0.
      endwhere
      
      !if(time_slot_num==1.and.case_num==1)theta(:,kde)=300.
      
      ! diag air status
      
      !print*,'nf90_put_var'
      status = nf90_put_var(ncid, rho_id  , rho  , start=(/1,1,time_slot_num/),count=(/nx,nz,1/))
      status = nf90_put_var(ncid, u_id    , u    , start=(/1,1,time_slot_num/),count=(/nx,nz,1/))
      status = nf90_put_var(ncid, w_id    , w    , start=(/1,1,time_slot_num/),count=(/nx,nz,1/))
      status = nf90_put_var(ncid, theta_id, theta, start=(/1,1,time_slot_num/),count=(/nx,nz,1/))
      status = nf90_put_var(ncid, q_id    , q    , start=(/1,1,time_slot_num/),count=(/nx,nz,1/))
      status = nf90_put_var(ncid, p_id    , p    , start=(/1,1,time_slot_num/),count=(/nx,nz,1/))
      status = nf90_put_var(ncid, T_id    , T    , start=(/1,1,time_slot_num/),count=(/nx,nz,1/))
      if(status/=nf90_noerr) call handle_err(status)
      
      !print*,'nf90_close'
      status = nf90_close(ncid)
      if(status/=nf90_noerr) call handle_err(status)
    
    end subroutine history_write_stat
end module io_mod
    