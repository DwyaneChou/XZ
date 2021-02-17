    module math_mod
    use constants_mod
    contains
    subroutine calc_polynomial_matrix(d,m,n,xi,eta,A)
      integer(i_kind), intent(in ) :: d ! polynomial degree
      integer(i_kind), intent(in ) :: m ! number of points
      integer(i_kind), intent(in ) :: n ! number of terms of polynomial, n = (d+1)*(d+2)/2
      real   (r_kind), intent(in ) :: xi (m)
      real   (r_kind), intent(in ) :: eta(m)
      real   (r_kind), intent(out) :: A  (m,n)
      
      real   (r_kind) :: x
      real   (r_kind) :: y
      integer(i_kind) :: iPOC
      integer(i_kind)  :: i,j,k
      
      do iPOC = 1,m
        x = xi (iPOC)
        y = eta(iPOC)
        
        k = 0
        do j = 0,d
          do i = 0,j
            k = k + 1
            A(iPOC,k) = x**real(j-i,r_kind)*y**real(i,r_kind)
          enddo
        enddo
      enddo
  
    end subroutine calc_polynomial_matrix
    
    subroutine calc_polynomial_deriv_matrix(d,m,n,xi,eta,px,py)
      integer(i_kind), intent(in ) :: d ! polynomial degree
      integer(i_kind), intent(in ) :: m ! number of points
      integer(i_kind), intent(in ) :: n ! number of points on cell, n = (d+1)*(d+2)/2
      real   (r_kind), intent(in ) :: xi (m)
      real   (r_kind), intent(in ) :: eta(m)
      real   (r_kind), intent(out) :: px (m,n)
      real   (r_kind), intent(out) :: py (m,n)
      
      real   (r16) :: x
      real   (r16) :: y
      integer(i4 ) :: iPOC
      integer(i4)  :: i,j,k
      real   (r16) :: powx1,powx2
      real   (r16) :: powy1,powy2
      real   (r16) :: coefx,coefy
      
      do iPOC = 1,m
        x = xi (iPOC)
        y = eta(iPOC)
        
        k = 0
        do j = 0,d
          do i = 0,j
            k = k + 1
            powx1 = merge( 0._r_kind, real(j-i-1,r_kind), real(j-i-1,r_kind)<0._r_kind )
            powx2 = real(i  ,r_kind)
            powy1 = real(j-i,r_kind)
            powy2 = merge( 0._r_kind, real(i-1  ,r_kind), real(i-1  ,r_kind)<0._r_kind )
            coefx = merge( 0._r_kind, real(j-i  ,r_kind), real(j-i-1,r_kind)<0._r_kind )
            coefy = merge( 0._r_kind, real(i    ,r_kind), real(i-1  ,r_kind)<0._r_kind )
            
            px(iPOC,k) = coefx * x**powx1 * y**powx2
            py(iPOC,k) = coefy * x**powy1 * y**powy2
          enddo
        enddo
      enddo
    
    end subroutine calc_polynomial_deriv_matrix
  
    subroutine  calc_polynomial_triangle_integration(d,c)
      integer(i_kind), intent(in ) :: d ! degree of polynomial
      real   (r_kind), intent(out) :: c(:)
      
      integer(i_kind) :: i,j,k,r
      
      k = 0
      c = 0
      do j = 0,d
        do i = 0,j
          k = k + 1
          do r = 0,i+1
            c(k) = c(k) + (-1)**(r) * nchoosek(i+1,r) / real( j - i + r + 1, r_kind )
          enddo
          c(k) = c(k) / real( i + 1, r_kind )
        enddo
      enddo
      
    end subroutine  calc_polynomial_triangle_integration
    
    subroutine  calc_polynomial_square_integration(d,x_min,x_max,y_min,y_max,c)
      integer(i_kind), intent(in ) :: d ! degree of polynomial
      real   (r_kind), intent(in ) :: x_min
      real   (r_kind), intent(in ) :: x_max
      real   (r_kind), intent(in ) :: y_min
      real   (r_kind), intent(in ) :: y_max
      real   (r_kind), intent(out) :: c(:)
      
      integer(i_kind) :: i,j,k
      
      k = 0
      c = 0
      do j = 0,d
        do i = 0,j
          k = k + 1
          c(k) = ( x_max**(j-i+1) - x_min**(j-i+1) ) * ( y_max**(i+1) - y_min**(i+1) ) / real( ( i + 1 ) * ( j - i + 1 ), r_kind )
        enddo
      enddo
      
    end subroutine  calc_polynomial_square_integration
    
    subroutine calc_rectangle_poly_matrix(nx,ny,m,xi,eta,A)
      integer(i_kind), intent(in ) :: nx ! number of points on x direction for reconstruction
      integer(i_kind), intent(in ) :: ny ! number of points on y direction for reconstruction
      integer(i_kind), intent(in ) :: m  ! number of unkonwn point values
      real   (r_kind), intent(in ) :: xi (m)
      real   (r_kind), intent(in ) :: eta(m)
      real   (r_kind), intent(out) :: A  (m,nx*ny)
      
      real   (r_kind) :: x
      real   (r_kind) :: y
      integer(i_kind) :: iPOC
      integer(i_kind)  :: i,j,k
      
      do iPOC = 1,m
        x = xi (iPOC)
        y = eta(iPOC)
        
        k = 0
        do j = 0,ny-1
          do i = 0,nx-1
            k = k + 1
            A(iPOC,k) = x**real(i,r_kind) * y**real(j,r_kind)
          enddo
        enddo
      enddo
  
    end subroutine calc_rectangle_poly_matrix
    
    subroutine calc_rectangle_poly_integration(nx,ny,x_min,x_max,y_min,y_max,c)
      integer(i_kind), intent(in ) :: nx ! number of points on x direction
      integer(i_kind), intent(in ) :: ny ! number of points on y direction
      real   (r_kind), intent(in ) :: x_min
      real   (r_kind), intent(in ) :: x_max
      real   (r_kind), intent(in ) :: y_min
      real   (r_kind), intent(in ) :: y_max
      real   (r_kind), intent(out) :: c(nx*ny)
      
      integer(i_kind) :: i,j,k
      
      k = 0
      c = 0
      do j = 0,ny-1
        do i = 0,nx-1
          k = k + 1
          c(k) = ( x_max**(i+1) - x_min**(i+1) ) * ( y_max**(j+1) - y_min**(j+1) ) / real( ( i + 1 ) * ( j + 1 ), r_kind )
        enddo
      enddo
      
    end subroutine  calc_rectangle_poly_integration
    
    function nchoosek(n,k) ! same as nchoosek in matlab
      real   (r_kind) :: nchoosek
      integer(i_kind) :: n
      integer(i_kind) :: k
      
      nchoosek = factorial(n) / ( factorial(n-k) * factorial(k) )
      
    end function nchoosek
    
    function factorial(n)
      real   (r_kind) :: factorial
      integer(i_kind) :: n
      
      factorial = gamma(real(n+1,r_kind))
    
    end function factorial
    
    ! calculate inverse matrix of A_input
    ! N is the order of matrix A_input and A
    ! A is inverse A_input
    ! L is status info
    SUBROUTINE BRINV(N,A_input,A,L)
    implicit none
    integer(i_kind),intent(in )           :: N
    real   (r_kind),intent(in )           :: A_input(N,N)
    real   (r_kind),intent(out)           :: A      (N,N)
    integer(i_kind),intent(out), optional :: L
    
    real    :: T,D
    integer :: IS(N),JS(N)
    integer :: i,j,k
    
    A = A_input
    
    if(present(L))L=1
    do K=1,N
      D=0.
      do I=K,N
        do J=K,N
          IF (ABS(A(I,J)).GT.D) THEN
            D=ABS(A(I,J))
            IS(K)=I
            JS(K)=J
          END IF
        enddo
      enddo
    
      IF (D+1.0.EQ.1.0) THEN
        if(present(L))L=0
        WRITE(*,*)'ERR**NOT INV'
        RETURN
      END IF
      
      do J=1,N
        T=A(K,J)
        A(K,J)=A(IS(K),J)
        A(IS(K),J)=T
      enddo
      
      do I=1,N
        T=A(I,K)
        A(I,K)=A(I,JS(K))
        A(I,JS(K))=T
      enddo
      
      A(K,K)=1/A(K,K)
      do J=1,N
        IF (J.NE.K) THEN
          A(K,J)=A(K,J)*A(K,K)
        END IF
      enddo
      
      do I=1,N
        IF (I.NE.K) THEN
          do J=1,N
            IF (J.NE.K) THEN
              A(I,J)=A(I,J)-A(I,K)*A(K,J)
            END IF
          enddo
        END IF
      enddo
      
      do I=1,N
        IF (I.NE.K) THEN
          A(I,K)=-A(I,K)*A(K,K)
        END IF
      enddo
    enddo
    
    do K=N,1,-1
      do J=1,N
        T=A(K,J)
        A(K,J)=A(JS(K),J)
        A(JS(K),J)=T
      enddo
      do I=1,N
        T=A(I,K)
        A(I,K)=A(I,IS(K))
        A(I,IS(K))=T
      enddo
    enddo
    RETURN
    END SUBROUTINE BRINV
    
    ! Gaussian quadrature from http://bbs.fcode.cn/thread-219-1-1.html, http://fcode.cn/algorithm-73-1.html
    Subroutine Gaussian_Legendre(n, p, w)
      Integer(i_kind),intent(in ) :: n
      Real   (r_kind),intent(out) :: p(n), w(n)
      
      Real   (r16   ) :: fn(n), ak(n)
      Real   (r16   ) :: m
      !定义数组,大小n由module开始声明。    
      Integer(i_kind) :: i, j
      j = 0 !赋值控制循环变量的初值           
      m = -1.000001 !设置计算域[-1，1] 的下限，即代替-1 
      Do i = 1, 200000 !这个循环次数应该是由步长0.00001决 定,计算方法：200000=2/0.000001     
        If (legendreP(m,n)*legendreP(m+0.00001,n)<0) Then !从下限处开始往上逐步累加，
          !由步长0.00001说明最多求解10^5个解
          j = j + 1 !记录这是第几个解
          fn(j) = bis(m, m+0.00001, n)
          !调用二分法求解程序在分好的一小段上求解，
          !将解存储在fn（j）
          ak(j) = 2.0_r_kind/(n*dLegendreP1(fn(j),n)*dLegendrePn(fn(j),n)) !高斯点的权重
          !Write (*, *) '高斯点序号', j
          !Write (*, *) '高斯点', fn(j)
          !Write (*, *) '高斯点权重', ak(j)
          p(j) = fn(j)
          w(j) = ak(j)
        End If
        m = m + 0.00001 !执行完一次判断m向前推进一步
      End Do
    End Subroutine Gaussian_Legendre
    
    Real(r16) Function legendreP(x,n) !定义Legendre函数
      Real   (r16   ),intent(in ) :: x
      Integer(i_kind),intent(in ) :: n
      
      Real   (r16   ) :: a(n) !a(n)代表n阶勒让德多项式
      Integer(i_kind) :: i
      a(1) = x !1阶勒让德多项式
      a(2) = 1.5_r_kind*(x**2) - 0.5_r_kind !2阶勒让德多项式
      Do i = 3, n
        a(i) = (2*i-1)*x*a(i-1)/i - (i-1)*a(i-2)/i
        !利用递推关系产生n阶勒让德多项式
      End Do
      legendreP = a(n) !生成的n阶勒让德多项式   
    End Function legendreP
    
    Real(r16) Function dLegendreP1(x,n)!生成的（n-1）阶勒让德多项式  
      Real   (r16   ),intent(in ) :: x
      Integer(i_kind),intent(in ) :: n   
      
      Real   (r16   ) :: a(n) !a(n-1)代表（n-1）阶勒让德多项式
      Integer :: i
      a(1) = x
      a(2) = 1.5_r_kind*x**2 - 0.5_r_kind
      Do i = 3, n - 1
        a(i) = (2*i-1)*x*a(i-1)/i - (i-1)*a(i-2)/i
      End Do
      dLegendreP1 = a(n-1) 
    End Function dLegendreP1
    
    Real (r16) Function dLegendrePn(x,n)!生成n阶勒让德多项式的导数表达式
      Real   (r16   ),intent(in ) :: x
      Integer(i_kind),intent(in ) :: n
      
      Real   (r16   ) :: a(n)
      Integer :: i
      a(1) = x
      a(2) = 1.5_r_kind*x**2 - 0.5_r_kind
      Do i = 3, n
        a(i) = (2*i-1)*x*a(i-1)/i - (i-1)*a(i-2)/i
      End Do
      dLegendrePn = n*a(n-1)/(1-x**2) - n*x*a(n)/(1-x**2)
      
    End Function dLegendrePn
    
    Real (r16) Function bis(a_in, b_in,n) !二分法求解函数的解
      Real   (r16   ), intent(in ) :: a_in, b_in
      integer(i_kind), intent(in ) :: n
      
      Real   (r16   )  :: a, b
      Real   (r16   ) :: c
      !a,b是传递进来的划分好的有一个解存在的区间
      a = a_in
      b = b_in
      
      Do
        c = (a+b)/2.0_r_kind
        If (legendreP(c,n)*legendreP(a,n)<0) Then
          b = c
        Else
          a = c
        End If
        If ((b-a)<1.e-16) exit 
      End Do
      bis = c!bis即是利用二分法求得的解
    End Function bis
    
    end module math_mod
    