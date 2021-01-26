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
      integer time_dimid
      
      integer x_id
      integer z_id
      
      integer zs_id

      integer sqrtGC_id
      integer sqrtGL_id
      integer sqrtGR_id
      integer sqrtGB_id
      integer sqrtGT_id
      
      integer G13C_id
      integer G13L_id
      integer G13R_id
      integer G13B_id
      integer G13T_id
      
      integer rho_id  
      integer u_id    
      integer w_id    
      integer theta_id
      integer q_id    
      integer p_id    
      integer T_id    
      
      real(r_kind), dimension(ids:ide,kds:kde) :: sqrtGLe
      real(r_kind), dimension(ids:ide,kds:kde) :: sqrtGRe
      real(r_kind), dimension(ids:ide,kds:kde) :: sqrtGBe
      real(r_kind), dimension(ids:ide,kds:kde) :: sqrtGTe
      real(r_kind), dimension(ids:ide,kds:kde) :: G13Le
      real(r_kind), dimension(ids:ide,kds:kde) :: G13Re
      real(r_kind), dimension(ids:ide,kds:kde) :: G13Be
      real(r_kind), dimension(ids:ide,kds:kde) :: G13Te
      
      integer i,k
      
      character*1 case_num_c
      
      do k = kds,kde
        do i = ids,ide
          sqrtGLe(i,k) = Gaussian_quadrature_1d(sqrtGL(:,i,k))
          sqrtGRe(i,k) = Gaussian_quadrature_1d(sqrtGR(:,i,k))
          sqrtGBe(i,k) = Gaussian_quadrature_1d(sqrtGB(:,i,k))
          sqrtGTe(i,k) = Gaussian_quadrature_1d(sqrtGT(:,i,k))
          G13Le  (i,k) = Gaussian_quadrature_1d(G13L  (:,i,k))
          G13Re  (i,k) = Gaussian_quadrature_1d(G13R  (:,i,k))
          G13Be  (i,k) = Gaussian_quadrature_1d(G13B  (:,i,k))
          G13Te  (i,k) = Gaussian_quadrature_1d(G13T  (:,i,k))
        enddo
      enddo
      
      write(case_num_c ,'(i1)')case_num
      
      stat_output_file = 'output_xz_'//case_num_c//'.nc'
      
      status = nf90_create(trim(stat_output_file), NF90_CLOBBER + NF90_64BIT_OFFSET , ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_def_dim(ncid,'nx'  , nx            , nx_dimid  )
      status = nf90_def_dim(ncid,'nz'  , nz            , nz_dimid  )
      status = nf90_def_dim(ncid,'time', NF90_UNLIMITED, time_dimid)
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_put_att(ncid,nf90_global, 'case_num' , case_num   )
      status = nf90_put_att(ncid,nf90_global, 'Rd'       , real(Rd,r8))
      status = nf90_put_att(ncid,nf90_global, 'Rv'       , real(Rv,r8))
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_def_var(ncid,'x'    ,NF90_DOUBLE,(/nx_dimid,nz_dimid/),x_id     )
      status = nf90_def_var(ncid,'z'    ,NF90_DOUBLE,(/nx_dimid,nz_dimid/),z_id     )
      status = nf90_def_var(ncid,'zs'   ,NF90_DOUBLE,(/nx_dimid,nz_dimid/),zs_id    )
      
      status = nf90_def_var(ncid,'sqrtGC',NF90_DOUBLE,(/nx_dimid,nz_dimid/),sqrtGC_id )
      status = nf90_def_var(ncid,'sqrtGL',NF90_DOUBLE,(/nx_dimid,nz_dimid/),sqrtGL_id )
      status = nf90_def_var(ncid,'sqrtGR',NF90_DOUBLE,(/nx_dimid,nz_dimid/),sqrtGR_id )
      status = nf90_def_var(ncid,'sqrtGB',NF90_DOUBLE,(/nx_dimid,nz_dimid/),sqrtGB_id )
      status = nf90_def_var(ncid,'sqrtGT',NF90_DOUBLE,(/nx_dimid,nz_dimid/),sqrtGT_id )
      
      status = nf90_def_var(ncid,'G13C'  ,NF90_DOUBLE,(/nx_dimid,nz_dimid/),G13C_id   )
      status = nf90_def_var(ncid,'G13L'  ,NF90_DOUBLE,(/nx_dimid,nz_dimid/),G13L_id   )
      status = nf90_def_var(ncid,'G13R'  ,NF90_DOUBLE,(/nx_dimid,nz_dimid/),G13R_id   )
      status = nf90_def_var(ncid,'G13B'  ,NF90_DOUBLE,(/nx_dimid,nz_dimid/),G13B_id   )
      status = nf90_def_var(ncid,'G13T'  ,NF90_DOUBLE,(/nx_dimid,nz_dimid/),G13T_id   )
      
      status = nf90_def_var(ncid,'rho'  ,NF90_DOUBLE,(/nx_dimid,nz_dimid,time_dimid/),rho_id   )
      status = nf90_def_var(ncid,'u'    ,NF90_DOUBLE,(/nx_dimid,nz_dimid,time_dimid/),u_id     )
      status = nf90_def_var(ncid,'w'    ,NF90_DOUBLE,(/nx_dimid,nz_dimid,time_dimid/),w_id     )
      status = nf90_def_var(ncid,'theta',NF90_DOUBLE,(/nx_dimid,nz_dimid,time_dimid/),theta_id )
      status = nf90_def_var(ncid,'q'    ,NF90_DOUBLE,(/nx_dimid,nz_dimid,time_dimid/),q_id     )
      status = nf90_def_var(ncid,'p'    ,NF90_DOUBLE,(/nx_dimid,nz_dimid,time_dimid/),p_id     )
      status = nf90_def_var(ncid,'T'    ,NF90_DOUBLE,(/nx_dimid,nz_dimid,time_dimid/),T_id     )
      
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_enddef(ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_put_var(ncid,x_id      , real(xC    (ids:ide,kds:kde),r8))
      status = nf90_put_var(ncid,z_id      , real(zC    (ids:ide,kds:kde),r8))
      status = nf90_put_var(ncid,zs_id     , real(zsC   (ids:ide,kds:kde),r8))
      
      status = nf90_put_var(ncid,sqrtGC_id , real(sqrtGC (ids:ide,kds:kde),r8))
      status = nf90_put_var(ncid,sqrtGL_id , real(sqrtGLe(ids:ide,kds:kde),r8))
      status = nf90_put_var(ncid,sqrtGR_id , real(sqrtGRe(ids:ide,kds:kde),r8))
      status = nf90_put_var(ncid,sqrtGB_id , real(sqrtGBe(ids:ide,kds:kde),r8))
      status = nf90_put_var(ncid,sqrtGT_id , real(sqrtGTe(ids:ide,kds:kde),r8))
      status = nf90_put_var(ncid,G13C_id   , real(G13C   (ids:ide,kds:kde),r8))
      status = nf90_put_var(ncid,G13L_id   , real(G13Le  (ids:ide,kds:kde),r8))
      status = nf90_put_var(ncid,G13R_id   , real(G13Re  (ids:ide,kds:kde),r8))
      status = nf90_put_var(ncid,G13B_id   , real(G13Be  (ids:ide,kds:kde),r8))
      status = nf90_put_var(ncid,G13T_id   , real(G13Te  (ids:ide,kds:kde),r8))
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_close(ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
    end subroutine history_init
    
    subroutine history_write_stat(stat,time_slot_num)
      type(stat_field), intent(in) :: stat
      integer         , intent(in) :: time_slot_num
      
      integer status
      integer ncid
      
      integer rho_id  
      integer u_id    
      integer w_id    
      integer theta_id
      integer q_id    
      integer p_id    
      integer T_id    
      
      real(r_kind), dimension(nx,nz) :: w1
      real(r_kind), dimension(nx,nz) :: w2
      real(r_kind), dimension(nx,nz) :: w3
      real(r_kind), dimension(nx,nz) :: w4
      real(r_kind), dimension(nx,nz) :: w5
      
      real(r_kind), dimension(nx,nz) :: rho
      real(r_kind), dimension(nx,nz) :: sqrtGrho
      real(r_kind), dimension(nx,nz) :: u
      real(r_kind), dimension(nx,nz) :: w
      real(r_kind), dimension(nx,nz) :: theta
      real(r_kind), dimension(nx,nz) :: q     !��ʪ
      real(r_kind), dimension(nx,nz) :: gamma !��ϱ�
      real(r_kind), dimension(nx,nz) :: T
      real(r_kind), dimension(nx,nz) :: p
      real(r_kind), dimension(nx,nz) :: kappa
      real(r_kind), dimension(nx,nz) :: Ra ! �����峣��
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
      w1 = stat%q(1,ids:ide,kds:kde)
      w2 = stat%q(2,ids:ide,kds:kde)
      w3 = stat%q(3,ids:ide,kds:kde)
      w4 = stat%q(4,ids:ide,kds:kde)
      w5 = stat%q(5,ids:ide,kds:kde)
      
      sqrtGrho = w1 + w5
      rho      = sqrtGrho / sqrtGC(ids:ide,kds:kde)
      u        = w2 / sqrtGrho
      w        = w3 / sqrtGrho
      theta    = w4 / w1
      gamma    = w5 / w1
      q        = gamma / ( 1. + gamma )
      Cp       = Cpd + ( Cpv - Cpd ) * q
      Cv       = Cvd + ( Cvv - Cvd ) * q
      kappa    = Cp / Cv
      Ra       = Rd * ( 1. + eq * q )
      p        = p0 * ( rho * Ra * theta / p0 )** kappa
      T        = p / ( rho * Ra )
      
      !print*,'nf90_put_var'
      status = nf90_put_var(ncid, rho_id  , real(rho  ,r8), start=(/1,1,time_slot_num/),count=(/nx,nz,1/))
      status = nf90_put_var(ncid, u_id    , real(u    ,r8), start=(/1,1,time_slot_num/),count=(/nx,nz,1/))
      status = nf90_put_var(ncid, w_id    , real(w    ,r8), start=(/1,1,time_slot_num/),count=(/nx,nz,1/))
      status = nf90_put_var(ncid, theta_id, real(theta,r8), start=(/1,1,time_slot_num/),count=(/nx,nz,1/))
      status = nf90_put_var(ncid, q_id    , real(q    ,r8), start=(/1,1,time_slot_num/),count=(/nx,nz,1/))
      status = nf90_put_var(ncid, p_id    , real(p    ,r8), start=(/1,1,time_slot_num/),count=(/nx,nz,1/))
      status = nf90_put_var(ncid, T_id    , real(T    ,r8), start=(/1,1,time_slot_num/),count=(/nx,nz,1/))
      if(status/=nf90_noerr) call handle_err(status)
      
      !print*,'nf90_close'
      status = nf90_close(ncid)
      if(status/=nf90_noerr) call handle_err(status)
    
    end subroutine history_write_stat
end module io_mod
    