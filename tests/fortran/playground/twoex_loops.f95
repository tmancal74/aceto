
program loops

  integer :: i, j, k, N, M, w
  
  N = 4
  
  M = 100000000
  
  do k = 1, M
  w = 0
  do i = 1, N-1
    do j = i+1, N
       !print *, i, j
       w = w + 1
    end do
  end do
  end do
  
  print *, w
  

end program loops