module linear_integration_mod
    implicit none

    contains
    subroutine spline2_integration(n,x,y,d2fa,d2fb,nt,t,ty)
    !---------------------------------subroutine  comment
    !
    !  Purpose   :   三次样条之第二类边界条件
    !  该函数不返回插值结果，而是返回积分结果，用于静力平衡的计算
    !-----------------------------------------------------
    !  Input  parameters  :
    !       1.   n-----插值节点个数减1，如有九个节点 则N=8
    !       2.   x ---节点自变量  为（0：N）维向量
    !       3.   y----节点因变量  （0：N）维向量
    !       4.   nt 要计算向量的维数
    !       5.   t 要计算的向量  (1:nt) 维向量
    !       6.   
    !       7.    d2fa,d2fb  起点于终点处的二阶导数条件，
    !              如果是自然边界条件，则二者为皆为0
    !  Output parameters  :
    !       1.   ty --积分结果向量，(1：nt)维向量
    !
    !  Common parameters  :
    !
    !----------------------------------------------------
    !  Post Script :
    !          
    !       1.   自然边界条件是第一类边界条件的特殊情况
    !            本程序可以直接处理第一类边界条件
    !       2.     
    !            对于要插值的向量分量，可以不必按大小排列
    !----------------------------------------------------
    
    implicit real(a-z)
    
    !n为插值节点数减1，即如果有9个节点，则n=8
    !nt 要插值点的个数，及向量t,ty的维数
    integer, intent(in ) :: n, nt
    real   , intent(in ) :: x(0:n),y(0:n)
    real   , intent(in ) :: t(nt)
    real   , intent(in ) :: d2fa
    real   , intent(in ) :: d2fb
    real   , intent(out) :: ty(nt)
    
    integer::i,j,k
    
    real::h(0:n-1)
    
    real::f1(0:n-1),f2(1:n-1)
    
    real::u(1:n-1),lambda(1:n-1),d(1:n-1)
    
    real::M(0:n),v(1:n-1)
    
    real::A(1:n-1,1:n-1)
    
    integer :: node(nt)  ! 记录每个被插点坐落于哪个节点范围内，如节点1到节点2之间都是节点1的控制范围
    integer :: dnode(nt) ! 记录每个被插点之间隔了多少个节点
    
    
    M(0) = d2fa
    M(n) = d2fb
    
    do i=0,n-1
      h (i) = x(i+1) - x(i)
      f1(i) = (y(i+1) - y(i)) / h(i)
    end do
    
    
    !求得 u, lambda, d
    do i=1,n-1
      u     (i) = h(i-1) / ( h(i-1) + h(i) )
      lambda(i) = 1 - u(i)
      
      f2(i) = ( f1(i-1) - f1(i) ) / ( x(i-1) - x(i+1) )
      
      d(i) = 6. * f2(i) 
    end do
    
    !设置A矩阵值
    A = 0
    do i = 1, n-1
      a(i,i) = 2.
    end do
    
    do i = 2,n-1
      a(i,i-1) = u(i)
    end do
    
    do i = 1, n-2
      a(i,i+1) = lambda(i)
    end do
    
    ! 设置右向量值
    d(1  ) = d(1  ) - u     (1  ) * M(0)
    d(n-1) = d(n-1) - lambda(n-1) * M(n)
    
    call chase(a,d,v,N-1)
    
    do i=1,n-1
      M(i) = v(i)
    end do
    
    !--------以上以及求得系数
    !已经完成插值多项式的建立
    
    !------------------------------------------------
    ! 以下开始计算具体值
    do k=1,nt
      !------------
      !  对要插值向量每个分量而言，先找到其在数据中的位置
      do i=1,n-1
        if (t(k)<x(i+1)) exit
      end do
      
      if(i==n)then
        ! 外插
        i = n - 1
      endif
      
      node(k) = i
    enddo
    dnode(2:nt) = node(2:nt) - node(1:nt-1)
    dnode(1   ) = 0
    
    !分段积分
    ty    = 0.
    k     = 1
    ty(k) = 0.
    do k = 2, nt
      if(dnode(k)==0)then
        ty(k) = ty(k) + intfunc( t(k),node(k) ) - intfunc( t(k-1),node(k) )
      else
        ty(k) = ty(k) + intfunc( x(node(k-1)+1),node(k-1) ) - intfunc( t(k-1),node(k-1) )
        
        if(dnode(k)>1)then
          do i = 1, dnode(k)
            ty(k) = ty(k) + intfunc( x(node(k-1)+i+1),node(k-1)+i ) - intfunc( x(node(k-1)+i),node(k-1)+i )
          enddo
        endif
        
        ty(k) = ty(k) + intfunc( t(k),node(k) ) - intfunc( x(node(k)),node(k) )
      endif
    enddo
    
    contains

      real function intfunc(tcoord,i) result(res)
        real   , intent(in ) :: tcoord ! target coordinate
        integer, intent(in ) :: i
        
        res = - M(i) * ( x(i+1) - tcoord )**4 / 24. / h(i) + M(i+1) * ( tcoord - x(i) )**4 / 24. / h(i) &
              - ( y(i  ) - M(i  ) * h(i)**2 / 6. ) / h(i) * 0.5 * ( x(i+1) - tcoord )**2                &
              + ( y(i+1) - M(i+1) * h(i)**2 / 6. ) / h(i) * 0.5 * ( tcoord - x(i  ) )**2
      
      end function intfunc
    
    end subroutine spline2_integration
    
    subroutine spline4_integration(n,x,y,nt,t,ty)
    !---------------------------------subroutine  comment
    !  Purpose   :   三次样条之第二类边界条件，端点二阶导数不变
    !  即M(1)=M(2); M(end-1)=SM(end)
    !  该函数不返回插值结果，而是返回积分结果，用于静力平衡的计算
    !-----------------------------------------------------
    !  Input  parameters  :
    !       1.   n-----插值节点个数减1，如有九个节点 则N=8
    !       2.   x ---节点自变量  为（0：N）维向量
    !       3.   y----节点因变量 （0：N）维向量
    !       4.   nt---要计算向量的维数
    !       5.   t----被插点坐标 (1:nt)维向量
    !  Output parameters  :
    !       1.   ty --积分结果向量，(1：nt)维向量
    !
    !  Common parameters  :
    !
    !----------------------------------------------------
    !  Post Script :
    !          
    !       1.   自然边界条件是第一类边界条件的特殊情况
    !            本程序可以直接处理第一类边界条件
    !       2.     
    !            对于要插值的向量分量，可以不必按大小排列
    !----------------------------------------------------
    
    implicit real(a-z)
    
    !n为插值节点数减1，即如果有9个节点，则n=8
    !nt 要插值点的个数，及向量t,ty的维数
    integer, intent(in ) :: n, nt
    real   , intent(in ) :: x(0:n),y(0:n)
    real   , intent(in ) :: t(nt)
    real   , intent(out) :: ty(nt)
    
    integer::i,j,k
    
    real::h(0:n-1)
    
    real::f1(0:n-1),f2(1:n-1)
    
    real::u(1:n-1),lambda(1:n-1),d(1:n-1)
    
    real::M(0:n),v(1:n-1)
    
    real::A(1:n-1,1:n-1)
    
    integer :: node(nt)  ! 记录每个被插点坐落于哪个节点范围内，如节点1到节点2之间都是节点1的控制范围
    integer :: dnode(nt) ! 记录每个被插点之间隔了多少个节点
    
    do i=0,n-1
      h (i) = x(i+1) - x(i)
      f1(i) = (y(i+1) - y(i)) / h(i)
    end do
    
    
    !求得 u, lambda, d
    do i=1,n-1
      u     (i) = h(i-1) / ( h(i-1) + h(i) )
      lambda(i) = 1 - u(i)
      
      f2(i) = ( f1(i-1) - f1(i) ) / ( x(i-1) - x(i+1) )
      
      d(i) = 6. * f2(i) 
    end do
    
    !设置A矩阵值
    A = 0.
    
    A(1  ,1  ) = u(1) + 2.
    A(n-1,n-1) = 2. + lambda(n-1)
    do i = 2, n-2
      a(i,i) = 2.
    end do
    
    do i = 2,n-1
      a(i,i-1) = u(i)
    end do
    
    do i = 1, n-2
      a(i,i+1) = lambda(i)
    end do
    
    call chase(a,d,v,N-1)
    
    do i=1,n-1
      M(i) = v(i)
    end do
    M(0) = M(1)
    M(n) = M(n-1)
    
    !--------以上以及求得系数
    !已经完成插值多项式的建立
    
    !------------------------------------------------
    ! 以下开始计算具体值
    do k=1,nt
      !------------
      !  对要插值向量每个分量而言，先找到其在数据中的位置
      do i=1,n-1
        if (t(k)<x(i+1)) exit
      end do
      
      if(i==n)then
        ! 外插
        i = n - 1
      endif
      
      node(k) = i
    enddo
    dnode(2:nt) = node(2:nt) - node(1:nt-1)
    dnode(1   ) = 0
    
    !分段积分
    ty    = 0.
    k     = 1
    ty(k) = 0.
    do k = 2, nt
      if(dnode(k)==0)then
        ty(k) = ty(k) + intfunc( t(k),node(k) ) - intfunc( t(k-1),node(k) )
      else
        ty(k) = ty(k) + intfunc( x(node(k-1)+1),node(k-1) ) - intfunc( t(k-1),node(k-1) )
        
        if(dnode(k)>1)then
          do i = 1, dnode(k)
            ty(k) = ty(k) + intfunc( x(node(k-1)+i+1),node(k-1)+i ) - intfunc( x(node(k-1)+i),node(k-1)+i )
          enddo
        endif
        
        ty(k) = ty(k) + intfunc( t(k),node(k) ) - intfunc( x(node(k)),node(k) )
      endif
    enddo
    
    contains

      real function intfunc(tcoord,i) result(res)
        real   , intent(in ) :: tcoord ! target coordinate
        integer, intent(in ) :: i
        
        res = - M(i) * ( x(i+1) - tcoord )**4 / 24. / h(i) + M(i+1) * ( tcoord - x(i) )**4 / 24. / h(i) &
              - ( y(i  ) - M(i  ) * h(i)**2 / 6. ) / h(i) * 0.5 * ( x(i+1) - tcoord )**2                &
              + ( y(i+1) - M(i+1) * h(i)**2 / 6. ) / h(i) * 0.5 * ( tcoord - x(i  ) )**2
      
      end function intfunc
    
    end subroutine spline4_integration
    
    subroutine chase(A,f,x,N)
    !---------------------------------subroutine  comment
    !  Version   :  V1.0    
    !  Coded by  :  syz 
    !  Date      :  2010-4-9
    !-----------------------------------------------------
    !  Purpose   :  追赶法计算三对角方程组
    !              Ax=f
    !-----------------------------------------------------
    !  Input  parameters  :
    !       1.  A系数矩阵
    !       2. f 右向量
    !  Output parameters  :
    !       1.  x方程的解
    !       2.  N维数
    !  Common parameters  :
    !
    !----------------------------------------------------
    !  Post Script :
    !       1.   注意：该方法仅适用于三对角方程组
    !       2.
    !---------------------------------------------------
     
     implicit real(a-z)
     integer::N
     real::A(N,N),f(N),x(N)
     real::L(2:N),u(N),d(1:N-1)
     
     real::c(1:N-1),b(N),e(2:N)
    
     integer::i
     
     real::y(N)
    
    !---------------把A矩阵复制给向量e,b,c 
    do i=1,N
      b(i)=a(i,i)
    end do 
    
    do i=1,N-1
      c(i)=a(i,i+1)
    end do
    
    do i=2,N
      e(i)=a(i,i-1)
    end do
    !------------------------
    
    do i=1,N-1
      d(i)=c(i)
    end do 
    
    u(1)=b(1)
    
    do i=2,N
      L(i)=e(i)/u(i-1)
      u(i)=b(i)-L(i)*c(i-1)
    end do
    
    !------开始回带,求得y
    
    y(1)=f(1)
    do i=2,N
      y(i)=f(i)-L(i)*y(i-1)
    end do
    
    !-----开始回带，求得x
    
    x(n)=y(n)/u(n)
    
    do i=n-1,1,-1
      x(i)=(y(i)-c(i)*x(i+1))/u(i)
    end do
    
    end subroutine chase
end module linear_integration_mod
