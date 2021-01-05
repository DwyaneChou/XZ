module gram_sch
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

subroutine  solve(A,b,x,M,N)
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
implicit real*8(a-z)

integer::M,N
real*8::A(M,N),Q(M,N),R(N,N)
real*8::b(M)
real*8::QT(N,M)  !Q的转置矩阵
real*8::QTb(N)   !Q'b
real*8::x(N)
 
call gram_dec(A,Q,R,M,N)

QT=transpose(Q)
QTb=matmul(QT,b)  !  Rx=Q'b

call uptri(R,QTb,x,N) !回带

end subroutine

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
implicit real*8(a-z)

integer::M,N
integer::i,j,k

real*8::A(M,N),Q(M,N),R(N,N)

real*8::vec_temp(M)



R(1,1)=dsqrt(dot_product(a(:,1),a(:,1)))

Q(:,1)=a(:,1)/R(1,1)


do k=2,N

      do j=1,k-1
        R(j,k)=dot_product(Q(:,j),A(:,k))   
      end do
   
      vec_temp=A(:,k)
   
      do j=1,k-1
   
        vec_temp=vec_temp-Q(:,j)*R(j,k)
   
      end do


    R(k,k)=dsqrt(dot_product(vec_temp,vec_temp))

     
    Q(:,k)=vec_temp/R(k,k)

end do
 
end subroutine gram_dec

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

implicit real*8(a-z)

integer::i,j,k,N

real*8::A(N,N),b(N),x(N)

x(N)=b(N)/A(N,N)

!回带部分
do i=n-1,1,-1
   
    x(i)=b(i)
   do j=i+1,N
    x(i)=x(i)-a(i,j)*x(j)
   end do
    x(i)=x(i)/A(i,i)

end do

end subroutine uptri



end module gram_sch



module driver
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  
!-----------------------------------------------------
!  Description : 
!                驱动程序模块
!                分别调用QR分解 ，以及用之解决最小二乘问题
!-----------------------------------------------------

contains

subroutine dri_main(M,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  
!-----------------------------------------------------
!  Purpose   :  驱动模块主函数
!    
!-----------------------------------------------------


use gram_sch
integer::M,N
real*8::A(m,n),Q(m,n),R(n,n)
real*8::b(M),x(N)

!读入矩阵
read(11,*)((A(i,j),j=1,n),i=1,m)

!读入b向量
read(11,*)b


call gram_dec(A,Q,R,m,n)

!-------------------这段程序用于输出QR分解
write(12,101)
101 format(T10,'修正的Gram-Schmidt方法QR分解计算最小二乘问题',/,T3,'Q=',/)

write(12,102)((Q(i,j),j=1,N),i=1,M)

!变量输出格式只针对 IVF编译器，在CVF中不支持
102 format(<M>(<N>F10.5/))

write(12,103)
103 format(T3,'R=',/)

write(12,104)((R(i,j),j=1,n),i=1,n)

!变量输出格式只针对 IVF编译器，在CVF中不支持
104 format(<N>(<N>F10.5/))

!------------------------------


!调用求解函数
call solve(A,b,x,M,N)

write(12,105)x

105 format(T3,'最小二乘意义下的解为',//,(<N>F10.5))

end subroutine dri_main

end module driver



program main

!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-11
!-----------------------------------------------------
!  Purpose   :  1. 采用修正的Gram-Schimidt方法进行QR分解
!               2. 分解后进行超定方程的最小二乘求解 
!-----------------------------------------------------
!  In put data  files :
!       1.         fin.txt 输入文件
!       2.   
!  Output data files  :
!       1.         fou.txt 输出文件
!       2.
!-----------------------------------------------------
!  Post Script :
!       1.         由主函数引导驱动程序
!       2.         该程序可以处以一般的最小二乘问题
!                  只要按照输入卡片输入数据，即可以计算
!-----------------------------------------------------
use driver
integer::M,N

open(unit=11,file='fin.txt')
open(unit=12,file='fout.txt')


read(11,*)
!读入A矩阵

!读入方程维数系数
read(11,*)M,N

!调用驱动函数
call dri_main(M,N)

end program main 