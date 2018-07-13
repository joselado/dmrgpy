subroutine traceover(wf,dmat,n1,n2)
implicit none
integer, intent(in) :: n1,n2 ! dimensions of the two spaces
complex (kind=8), intent(in) :: wf(n1*n2) ! wavefunction
complex (kind=8), intent(out) :: dmat(n1,n1) ! density matrix
integer :: i1,j1,i2,j2,i3,j3

dmat(:,:) = (0.d00,0d00) ! initialize density matrix

j3 = 1 ! initialize
do j1=1,n1
  do i1=1,n1
    do i2=1,n2
        j3 = n2*(j1-1) + i2
        i3 = n2*(i1-1) + i2
        dmat(i1,j1) = dmat(i1,j1) + wf(i3)*conjg(wf(j3)) 
    enddo
  enddo
enddo

return
end subroutine
