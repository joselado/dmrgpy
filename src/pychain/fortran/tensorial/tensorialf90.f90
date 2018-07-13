subroutine tensorial_operator(op1,op2,op3,n1,n2)
implicit none
integer, intent(in) :: n1,n2 ! dimensions
complex (kind=8), intent(in) :: op1(n1,n1),op2(n2,n2)  ! input matrices
complex (kind=8), intent(out) :: op3(n1*n2,n1*n2)  ! output matrix
integer :: i1,j1,i2,j2,i3,j3

! four loops
do j1=1,n1 ! loop over states j op1
  do j2=1,n2  ! loop over states j op2
    do i1=1,n1 ! loop over states i op1
      do i2=1,n2  ! loop over states i op2
        j3 = n2*(j1-1) + j2
        i3 = n2*(i1-1) + i2
        op3(i3,j3) = op1(i1,j1)*op2(i2,j2)
      enddo
    enddo
  enddo
enddo




return
end subroutine


subroutine sparse_tensorial_operator(data1,row1,col1,data2,row2,col2, &
                                      data3,row3,col3,n1,n2,dim1,dim2)
implicit none
integer, intent(in) :: n1,n2 ! number of non-zero entries 
integer, intent(in) :: dim1,dim2 ! dimension of the matrices 
complex (kind=8), intent(in) :: data1(n1),data2(n2)  ! input matrices
integer, intent(in) :: row1(n1),row2(n2)  ! input matrices
integer, intent(in) :: col1(n1),col2(n2)  ! input matrices
! outputs
complex (kind=8), intent(out) :: data3(n1*n2) ! out matrix
integer, intent(out) :: col3(n1*n2),row3(n1*n2)  ! output indexes matrices
integer :: i1,j1,i2,j2,i3,j3,ii,jj,kk

! two loops
kk = 1 ! index for output
do ii=1,n1 ! loop over non zero entries in op1
  do jj=1,n2  ! loop over non zero entries in op2
      j1 = col1(ii) ! index of first matrix
      i1 = row1(ii) ! index of first matrix
      j2 = col2(jj) ! index of first matrix
      i2 = row2(jj) ! index of first matrix
!      write(*,*) i1,j1
      i3 = dim2*(i1-1) + i2 ! index for product
      j3 = dim2*(j1-1) + j2 ! index for product
      data3(kk) = data1(ii)*data2(jj) ! product
      row3(kk) = i3 ! product
      col3(kk) = j3 ! product
      kk = kk + 1
  enddo
enddo


return
end subroutine




subroutine tensorial_abv(data1,row1,col1,data2,row2,col2,vin,vout,n1,n2,dim1,dim2)
implicit none
integer, intent(in) :: n1,n2 ! number of non-zero entries 
integer, intent(in) :: dim1,dim2 ! dimension of the matrices 
complex (kind=8), intent(in) :: data1(n1),data2(n2)  ! input matrices
complex (kind=8), intent(in) :: vin(dim1*dim2) ! input vector
integer, intent(in) :: row1(n1),row2(n2)  ! input matrices
integer, intent(in) :: col1(n1),col2(n2)  ! input matrices
! outputs
complex (kind=8), intent(out) :: vout(dim1*dim2) ! input vector
integer :: i1,j1,i2,j2,i3,j3,ii,jj,kk

vout(:) = (0.d00,0.d00) ! initialize
! two loops
kk = 1 ! index for output
do ii=1,n1 ! loop over non zero entries in op1
  do jj=1,n2  ! loop over non zero entries in op2
      j1 = col1(ii) ! index of first matrix
      i1 = row1(ii) ! index of first matrix
      j2 = col2(jj) ! index of first matrix
      i2 = row2(jj) ! index of first matrix
!      write(*,*) i1,j1
      i3 = dim2*(i1-1) + i2 ! index for product
      j3 = dim2*(j1-1) + j2 ! index for product
      vout(i3) = vout(i3) + data1(ii)*data2(jj)*vin(j3) ! add contribution
  enddo
enddo


return
end subroutine 



subroutine tensorial_idbv(dim1,data2,row2,col2,vin,vout,n2,dim2)
implicit none
integer, intent(in) :: n2 ! number of non-zero entries 
integer, intent(in) :: dim1,dim2 ! dimension of the matrices 
complex (kind=8), intent(in) :: data2(n2)  ! input matrices
complex (kind=8), intent(in) :: vin(dim1*dim2) ! input vector
integer, intent(in) :: row2(n2)  ! input matrices
integer, intent(in) :: col2(n2)  ! input matrices
! outputs
complex (kind=8), intent(out) :: vout(dim1*dim2) ! input vector
integer :: i1,i2,j2,i3,j3,ii,jj,kk

vout(:) = (0.d00,0.d00) ! initialize
! two loops
kk = 1 ! index for output
do ii=1,dim1 ! loop over non zero entries in op1
  kk = dim2*(ii-1)
  do jj=1,n2  ! loop over non zero entries in op2
      j2 = col2(jj) ! index of first matrix
      i2 = row2(jj) ! index of first matrix
      i3 = kk + i2 ! index for product
      j3 = kk + j2 ! index for product
      vout(i3) = vout(i3) + data2(jj)*vin(j3) ! add contribution
  enddo
enddo


return
end subroutine 



subroutine tensorial_aidv(dim2,data1,row1,col1,vin,vout,n1,dim1)
implicit none
integer, intent(in) :: n1 ! number of non-zero entries 
integer, intent(in) :: dim1,dim2 ! dimension of the matrices 
complex (kind=8), intent(in) :: data1(n1)  ! input matrices
complex (kind=8), intent(in) :: vin(dim1*dim2) ! input vector
integer, intent(in) :: row1(n1)  ! input matrices
integer, intent(in) :: col1(n1)  ! input matrices
! outputs
complex (kind=8), intent(out) :: vout(dim1*dim2) ! input vector
integer :: i1,j1,i2,j2,i3,j3,ii,jj,kk1,kk2

vout(:) = (0.d00,0.d00) ! initialize
! two loops
do ii=1,n1 ! loop over non zero entries in op1
  j1 = col1(ii) ! index of first matrix
  i1 = row1(ii) ! index of first matrix
  kk1 = dim2*(i1-1)
  kk2 = dim2*(j1-1)
  do jj=1,dim2  ! loop over non zero entries in op2
      i3 = kk1 + jj ! index for product
      j3 = kk2 + jj ! index for product
      vout(i3) = vout(i3) + data1(ii)*vin(j3) ! add contribution
  enddo
enddo


return
end subroutine 
