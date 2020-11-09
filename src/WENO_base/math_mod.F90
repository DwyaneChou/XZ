module math_mod
  
  contains
  ! calculate inverse matrix of A_input
  ! N is the order of matrix A_input and A
  ! A is inverse A_input
  ! L is status info
  SUBROUTINE BRINV(N,A_input,A,L)
  implicit none
  integer,intent(in )           :: N
  real   ,intent(in )           :: A_input(N,N)
  real   ,intent(out)           :: A      (N,N)
  integer,intent(out), optional :: L
  
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
    