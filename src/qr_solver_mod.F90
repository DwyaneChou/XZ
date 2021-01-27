    module qr_solver_mod
    use constants_mod
    implicit none
    !----------------------------------------module coment
    !  Version     :  V1.0    
    !  Coded by    :  syz 
    !  Date        :  
    !-----------------------------------------------------
    !  Description :   修正的Gram-Schimdt正交化求最小二
    !                    乘问题模块
    !    
    !-----------------------------------------------------
    !  Parameters  :
    !      1.    solve     解超定方程 方法函数
    !      2.    gram_dec  G-S ,QR分解
    !      3.    uptri     上三角方程回带函数
    !      4.
    !-----------------------------------------------------
    !  Post Script :
    !      1.      即可以单独调用 QR分解函数
    !      2.      也可以调用解方程函数
    !-----------------------------------------------------
    
    contains
    
      subroutine qr_solver(A,b,x,M,N)
        !---------------------------------subroutine  comment
        !  Version   :  V1.0    
        !  Coded by  :  syz 
        !  Date      :  
        !-----------------------------------------------------
        !  Purpose   :  通过修正的Gram-Schmidt正交化求最小二问题
        !               方法函数
        !-----------------------------------------------------
        !  Method    :
        !               对超定方程  A进行QR分解后   方程变为
        !                   QR x=b
        !                => Rx=Q'b   R为上三角阵
        !                => 回带，可以求得最小二乘意义下的解
        !-----------------------------------------------------
        !  Post Script :
        !       1.       即求解 超定方程组  Ax=b    其中 A(M,N)  M>N    
        !       2.
        !----------------------------------------------------
        implicit real(r_kind)(a-z)
        
        integer(i_kind)::M,N
        real   (r_kind)::A(M,N),Q(M,N),R(N,N)
        real   (r_kind)::b(M)
        real   (r_kind)::QT(N,M)  !Q的转置矩阵
        real   (r_kind)::QTb(N)   !Q'b
        real   (r_kind)::x(N)
        real   (r_kind)::AT (N,M)
        real   (r_kind)::ATA(N,N)
        real   (r_kind)::ATb(N)
        
        !  For LAPACK only
        real   (r_kind), dimension(m+n) :: work
        real   (r_kind), dimension(n  ) :: tau
        integer(i_kind) :: INFO
        
        integer i,j
        
        ! Modified Gram-Schmidt
        call gram_dec(A,Q,R,M,N)
        
        QT=transpose(Q)
        QTb=matmul(QT,b)  !  Rx=Q'b
        
        call uptri(R,QTb,x,N) !回带
        
        !! Cholesky QR(factorization)
        !AT  = transpose(A)
        !ATA = matmul( AT, A )
        !ATb = matmul( AT, b )
        !call chol_eq(ATA,ATb,x,N)
        
      end subroutine qr_solver
      
      subroutine gram_dec(A,Q,R,M,N)
        !---------------------------------subroutine  comment
        !  Version   :  V1.0    
        !  Coded by  :  syz 
        !  Date      :  
        !-----------------------------------------------------
        !  Purpose   :   采用修正的 Gram-Schmidt分解求矩阵的QR分解
        !    
        !-----------------------------------------------------
        !  Input  parameters  :
        !       1.    A原始矩阵
        !       2.    A(M,N)
        !  Output parameters  :
        !       1.    分解结果为   Q(M,N):注意 Q不是方阵，Q列向量为标准正交基
        !       2.                 R(N,N)：R是方阵
        !       3.   
        !----------------------------------------------------
        !  Post Script :
        !       1.  注意矩阵的维数，分解后Q列向量是正交的
        !       2.  关于编程方法可以参看《矩阵分析与应用》张贤达编著
        !       3.  详细的数学解释，可以参看 麻省理工学院的
        !           线性代数教材《Linear Algebra with Application》
        !----------------------------------------------------
        implicit real(r_kind)(a-z)
        
        integer(i_kind)::M,N
        integer(i_kind)::i,j,k
        
        real(r_kind)::A(M,N),Q(M,N),R(N,N)
        
        real(r_kind)::vec_temp(M)
        
        R(1,1)=sqrt(dot_product(A(:,1),A(:,1)))
        Q(:,1)=A(:,1)/R(1,1)
        
        do k=2,N
          do j=1,k-1
            R(j,k)=dot_product(Q(:,j),A(:,k))
          end do
          
          vec_temp=A(:,k)
          
          do j=1,k-1
            vec_temp=vec_temp-Q(:,j)*R(j,k)
          end do
        
          R(k,k)=sqrt(dot_product(vec_temp,vec_temp))
          Q(:,k)=vec_temp/R(k,k)
        end do
       
      end subroutine gram_dec
      
      subroutine chol_eq(A,b,x,N)
        !---------------------------------subroutine  comment
        !  Version   :  V1.0    
        !  Coded by  :  syz 
        !  Date      :  
        !-----------------------------------------------------
        !  Purpose   :  用cholesky分解方法解方程
        !    
        !-----------------------------------------------------
        integer(i_kind)N
        real   (r_kind)::A(N,N),b(N),x(N)
        real   (r_kind)::L(N,N),y(N),LT(N,N)!LT 为L的转置矩阵
        
        integer(i_kind)i,j
        
        call chol(A,L,N)
        
        call downtri(L,b,y,N)
        
        LT = transpose(L)
        
        call uptri(LT,y,x,N)  !这一步已经算出了x
      
      end subroutine chol_eq
      
      subroutine  chol(A,L,N)
        !---------------------------------subroutine  comment
        !  Version   :  V1.0    
        !  Coded by  :  syz 
        !  Date      :  
        !-----------------------------------------------------
        !  Purpose   :  Cholesky分解子程序
        !-----------------------------------------------------  
        integer(i_kind):: N
        real   (r_kind):: A(N,N),L(N,N)
        real   (r_kind):: s
        integer(i_kind):: i,j,k
        
        L = 0
        
        L(1,1) = sqrt(A(1,1))
        
        L(2:,1) = A(2:,1) / L(1,1)
        
        do j=2,N
          s = dot_product(L(j,1:j-1),L(j,1:j-1))
          
          L(j,j) = sqrt( A(j,j) - s )
          
          !注意i范围
          do i=j+1,N
            s = dot_product(L(i,1:j-1),L(j,1:j-1))
            
            L(i,j) = ( A(i,j) - s ) / L(j,j)
          end do
        end do 
      
      end subroutine chol
      
      subroutine uptri(A,b,x,N)
        !---------------------------------subroutine  comment
        !  Version   :  V1.0    
        !  Coded by  :  syz 
        !  Date      :  2010-4-8
        !-----------------------------------------------------
        !  Purpose   :  上三角方程组的回带方法
        !                 Ax=b
        !-----------------------------------------------------
        !  Input  parameters  :
        !       1.   A(N,N)系数矩阵
        !       2.   b(N)右向量
        !       3.   N方程维数
        !  Output parameters  :
        !       1.  x  方程的根
        !       2.
        !  Common parameters  :
        !
        !----------------------------------------------------
        integer(i_kind):: N
        real   (r_kind):: A(N,N),b(N),x(N)
        integer(i_kind):: i,j,k
        
        x(N) = b(N) / A(N,N)
        
        !回带部分
        do i=n-1,1,-1
          x(i) = ( b(i) - dot_product( A(i,i+1:N), x(i+1:N)) ) / A(i,i)
        end do
      end subroutine uptri
      
      subroutine downtri(A,b,x,N)
        !---------------------------------subroutine  comment
        !  Version   :  V1.0    
        !  Coded by  :  syz 
        !  Date      :  2010-4-9
        !-----------------------------------------------------
        !  Purpose   :  下三角方程组的回带方法
        !                 Ax=b
        !-----------------------------------------------------
        integer(i_kind):: N
        real   (r_kind):: A(N,N),b(N),x(N)
        integer(i_kind):: i,j,k
        
        x(1) = b(1) / a(1,1)
        
        do k=2,N
          x(k) = ( b(k) - dot_product( A(k,1:k-1), x(1:k-1)) ) / A(k,k)
        end do
      
      end subroutine downtri
      
    end module qr_solver_mod
    
