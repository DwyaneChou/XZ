module linear_integration_mod
    implicit none

    contains
    subroutine spline2_integration(n,x,y,d2fa,d2fb,nt,t,ty)
    !---------------------------------subroutine  comment
    !
    !  Purpose   :   ��������֮�ڶ���߽�����
    !  �ú��������ز�ֵ��������Ƿ��ػ��ֽ�������ھ���ƽ��ļ���
    !-----------------------------------------------------
    !  Input  parameters  :
    !       1.   n-----��ֵ�ڵ������1�����оŸ��ڵ� ��N=8
    !       2.   x ---�ڵ��Ա���  Ϊ��0��N��ά����
    !       3.   y----�ڵ������  ��0��N��ά����
    !       4.   nt Ҫ����������ά��
    !       5.   t Ҫ���������  (1:nt) ά����
    !       6.   
    !       7.    d2fa,d2fb  ������յ㴦�Ķ��׵���������
    !              �������Ȼ�߽������������Ϊ��Ϊ0
    !  Output parameters  :
    !       1.   ty --���ֽ��������(1��nt)ά����
    !
    !  Common parameters  :
    !
    !----------------------------------------------------
    !  Post Script :
    !          
    !       1.   ��Ȼ�߽������ǵ�һ��߽��������������
    !            ���������ֱ�Ӵ����һ��߽�����
    !       2.     
    !            ����Ҫ��ֵ���������������Բ��ذ���С����
    !----------------------------------------------------
    
    implicit real(a-z)
    
    !nΪ��ֵ�ڵ�����1���������9���ڵ㣬��n=8
    !nt Ҫ��ֵ��ĸ�����������t,ty��ά��
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
    
    integer :: node(nt)  ! ��¼ÿ��������������ĸ��ڵ㷶Χ�ڣ���ڵ�1���ڵ�2֮�䶼�ǽڵ�1�Ŀ��Ʒ�Χ
    integer :: dnode(nt) ! ��¼ÿ�������֮����˶��ٸ��ڵ�
    
    
    M(0) = d2fa
    M(n) = d2fb
    
    do i=0,n-1
      h (i) = x(i+1) - x(i)
      f1(i) = (y(i+1) - y(i)) / h(i)
    end do
    
    
    !��� u, lambda, d
    do i=1,n-1
      u     (i) = h(i-1) / ( h(i-1) + h(i) )
      lambda(i) = 1 - u(i)
      
      f2(i) = ( f1(i-1) - f1(i) ) / ( x(i-1) - x(i+1) )
      
      d(i) = 6. * f2(i) 
    end do
    
    !����A����ֵ
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
    
    ! ����������ֵ
    d(1  ) = d(1  ) - u     (1  ) * M(0)
    d(n-1) = d(n-1) - lambda(n-1) * M(n)
    
    call chase(a,d,v,N-1)
    
    do i=1,n-1
      M(i) = v(i)
    end do
    
    !--------�����Լ����ϵ��
    !�Ѿ���ɲ�ֵ����ʽ�Ľ���
    
    !------------------------------------------------
    ! ���¿�ʼ�������ֵ
    do k=1,nt
      !------------
      !  ��Ҫ��ֵ����ÿ���������ԣ����ҵ����������е�λ��
      do i=1,n-1
        if (t(k)<x(i+1)) exit
      end do
      
      if(i==n)then
        ! ���
        i = n - 1
      endif
      
      node(k) = i
    enddo
    dnode(2:nt) = node(2:nt) - node(1:nt-1)
    dnode(1   ) = 0
    
    !�ֶλ���
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
    !  Purpose   :   ��������֮�ڶ���߽��������˵���׵�������
    !  ��M(1)=M(2); M(end-1)=SM(end)
    !  �ú��������ز�ֵ��������Ƿ��ػ��ֽ�������ھ���ƽ��ļ���
    !-----------------------------------------------------
    !  Input  parameters  :
    !       1.   n-----��ֵ�ڵ������1�����оŸ��ڵ� ��N=8
    !       2.   x ---�ڵ��Ա���  Ϊ��0��N��ά����
    !       3.   y----�ڵ������ ��0��N��ά����
    !       4.   nt---Ҫ����������ά��
    !       5.   t----��������� (1:nt)ά����
    !  Output parameters  :
    !       1.   ty --���ֽ��������(1��nt)ά����
    !
    !  Common parameters  :
    !
    !----------------------------------------------------
    !  Post Script :
    !          
    !       1.   ��Ȼ�߽������ǵ�һ��߽��������������
    !            ���������ֱ�Ӵ����һ��߽�����
    !       2.     
    !            ����Ҫ��ֵ���������������Բ��ذ���С����
    !----------------------------------------------------
    
    implicit real(a-z)
    
    !nΪ��ֵ�ڵ�����1���������9���ڵ㣬��n=8
    !nt Ҫ��ֵ��ĸ�����������t,ty��ά��
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
    
    integer :: node(nt)  ! ��¼ÿ��������������ĸ��ڵ㷶Χ�ڣ���ڵ�1���ڵ�2֮�䶼�ǽڵ�1�Ŀ��Ʒ�Χ
    integer :: dnode(nt) ! ��¼ÿ�������֮����˶��ٸ��ڵ�
    
    do i=0,n-1
      h (i) = x(i+1) - x(i)
      f1(i) = (y(i+1) - y(i)) / h(i)
    end do
    
    
    !��� u, lambda, d
    do i=1,n-1
      u     (i) = h(i-1) / ( h(i-1) + h(i) )
      lambda(i) = 1 - u(i)
      
      f2(i) = ( f1(i-1) - f1(i) ) / ( x(i-1) - x(i+1) )
      
      d(i) = 6. * f2(i) 
    end do
    
    !����A����ֵ
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
    
    !--------�����Լ����ϵ��
    !�Ѿ���ɲ�ֵ����ʽ�Ľ���
    
    !------------------------------------------------
    ! ���¿�ʼ�������ֵ
    do k=1,nt
      !------------
      !  ��Ҫ��ֵ����ÿ���������ԣ����ҵ����������е�λ��
      do i=1,n-1
        if (t(k)<x(i+1)) exit
      end do
      
      if(i==n)then
        ! ���
        i = n - 1
      endif
      
      node(k) = i
    enddo
    dnode(2:nt) = node(2:nt) - node(1:nt-1)
    dnode(1   ) = 0
    
    !�ֶλ���
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
    !  Purpose   :  ׷�Ϸ��������ԽǷ�����
    !              Ax=f
    !-----------------------------------------------------
    !  Input  parameters  :
    !       1.  Aϵ������
    !       2. f ������
    !  Output parameters  :
    !       1.  x���̵Ľ�
    !       2.  Nά��
    !  Common parameters  :
    !
    !----------------------------------------------------
    !  Post Script :
    !       1.   ע�⣺�÷��������������ԽǷ�����
    !       2.
    !---------------------------------------------------
     
     implicit real(a-z)
     integer::N
     real::A(N,N),f(N),x(N)
     real::L(2:N),u(N),d(1:N-1)
     
     real::c(1:N-1),b(N),e(2:N)
    
     integer::i
     
     real::y(N)
    
    !---------------��A�����Ƹ�����e,b,c 
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
    
    !------��ʼ�ش�,���y
    
    y(1)=f(1)
    do i=2,N
      y(i)=f(i)-L(i)*y(i-1)
    end do
    
    !-----��ʼ�ش������x
    
    x(n)=y(n)/u(n)
    
    do i=n-1,1,-1
      x(i)=(y(i)-c(i)*x(i+1))/u(i)
    end do
    
    end subroutine chase
end module linear_integration_mod
