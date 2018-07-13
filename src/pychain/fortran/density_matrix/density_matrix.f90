subroutine density_matrix(wf,dimwf,dmat)
implicit none
integer, intent(in) :: dimwf ! dimension of the wavefunction
complex (kind=8), intent(in) :: wf(dimwf) ! wavefunction
complex (kind=8), intent(out) :: dmat(dimwf,dimwf) ! wavefunction
integer :: i,j

dmat(:,:) = (0.d00,0.d00) ! initialize

do i=1,dimwf
  do j=1,dimwf
    dmat(i,j) = wf(i)*conjg(wf(j))
  enddo
enddo

return
end subroutine density_matrix



subroutine reduced_density_matrix(wf,basis,nsites,dimwf, &
                                   dsite,dmat,dimdmat)
implicit none
integer, intent(in) :: dimwf ! dimension of the wavefunction
complex (kind=8), intent(in) :: wf(dimwf) ! wavefunction
integer, intent(in) :: nsites ! number of sites that the basis has
! set of integer which label each element
integer, intent(in) :: basis(nsites,dimwf) 
integer, intent(in) :: dsite ! site where we calculate the density matrix
integer, intent(in) :: dimdmat ! dimension of the density matrix
complex (kind=8), intent(out) :: dmat(dimdmat,dimdmat) ! wavefunction

! internal variables
integer :: i,j,ii,jj
logical :: equal

dmat(:,:) = (0.d00,0.d00) ! initialize

do i=1,dimwf
  ii = basis(dsite,i) ! index in the reduced subspace
  do j=1,dimwf
    jj = basis(dsite,j) ! index in the reduced subspace
    ! now check whether the two vectors only differ by dsite 
    call compare_vectors(basis(:,i),basis(:,j),nsites,dsite,equal)
    if (equal) dmat(ii,jj) = dmat(ii,jj) + wf(i)*conjg(wf(j))
  enddo
enddo


end subroutine reduced_density_matrix


! check whether two vectors are the same, skiping an element
subroutine compare_vectors(v1,v2,n,skip,equal)
integer, intent(in) :: n
integer, intent(in) :: v1(n)
integer, intent(in) :: v2(n)
integer, intent(in) :: skip ! index to skip
logical, intent(out) :: equal ! whether they are the same or not

integer :: k ! counter

equal = .true. ! assuma than they are the same

do k=1,n
  if (k.eq.skip) cycle ! skip the special site, next iteration
  if (v1(k).ne.v2(k)) then
    equal = .false. ! say that they are different
    exit ! end the loop
  endif
enddo

return
endsubroutine compare_vectors
