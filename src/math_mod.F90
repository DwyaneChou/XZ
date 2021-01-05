    module math_mod
    use constants_mod
    contains
    subroutine calc_polynomial_matrix(d,n,xi,eta,A)
      integer(i_kind), intent(in ) :: d ! polynomial degree
      integer(i_kind), intent(in ) :: n ! number of points on cell, n = (d+1)*(d+2)/2
      real   (r_kind), intent(in ) :: xi (n)
      real   (r_kind), intent(in ) :: eta(n)
      real   (r_kind), intent(out) :: A  (n,n)
      
      real   (r_kind) :: x
      real   (r_kind) :: y
      integer(i_kind) :: iPOC
      integer(i_kind)  :: i,j,k
      
      do iPOC = 1,n
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
    
    subroutine calc_polynomial_deriv_matrix(d,n,xi,eta,px,py)
      integer(i_kind), intent(in ) :: d ! polynomial degree
      integer(i_kind), intent(in ) :: n ! number of points on cell, n = (d+1)*(d+2)/2
      real   (r_kind), intent(in ) :: xi (n)
      real   (r_kind), intent(in ) :: eta(n)
      real   (r_kind), intent(out) :: px (n,n)
      real   (r_kind), intent(out) :: py (n,n)
      
      real   (r16) :: x
      real   (r16) :: y
      integer(i4 ) :: iPOC
      integer(i4)  :: i,j,k
      real   (r16) :: powx1,powx2
      real   (r16) :: powy1,powy2
      real   (r16) :: coefx,coefy
      
      do iPOC = 1,n
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
      integer(i_kind),               intent(in ) :: d ! degree of polynomial
      real   (r_kind), dimension(:), intent(out) :: c
      
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
      integer(i_kind),               intent(in ) :: d ! degree of polynomial
      real   (r_kind)              , intent(in ) :: x_min
      real   (r_kind)              , intent(in ) :: x_max
      real   (r_kind)              , intent(in ) :: y_min
      real   (r_kind)              , intent(in ) :: y_max
      real   (r_kind), dimension(:), intent(out) :: c
      
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
    
    end module math_mod
    