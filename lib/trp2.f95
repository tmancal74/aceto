!*******************************************************************************
!
!   primitive_trp2, integral_trp2
!
!
!   Accelerated Charge and Energy Transfer Objects (ACETO) library
!
!
!
!*******************************************************************************
!
! Trapeziodal rule integration of a real or complex function
!
! Integration using trapoziodal rule for each interval
! accelerated using openACC standard directives.
!
! Parameters
! ----------
! x : real array
!   x-axis of the function y(x) as an array of values.
!   If the length of the array is the same as the length N of the array y,
!   it is assumed that the distances between points in x may differ.
!   If the lenth is 1, it is understood to be the step of x, which is then
!   understood as having equidistant points. The length of x
!   may only be 1 or N.
!
! y : real or complex array
!   On input, y represents the function to be integrated in the points
!   specified by the array x.
!   On output, y(i) contains the integral from x(1) to x(i). y(N) contains
!   the definite integral over the while interval specified by array x.
!
!*******************************************************************************

!-------------------------------------------------------------------------------
!
!
!   No Acceleration
!
!
!-------------------------------------------------------------------------------
subroutine primitive_trp2_dp_none(x,y, yout)
  use iso_c_binding
  implicit none
  !
  ! Real double precision version of the routine, accelerate by OpenACC directives
  !
  real(c_double), dimension(:), intent(in) :: x
  real(c_double), dimension(:,:), intent(in) :: y
  real(c_double), dimension(:,:), intent(out) :: yout

  ! local variables
  integer :: N, Nd
  integer :: Nx, i, j
  logical :: equidistant
  real(c_double), dimension(size(y,1)) :: yy
  real(c_double) :: dt

  N = size(y,1)
  Nd = size(y,2)
  Nx = size(x,1)

  equidistant = .false.
  if (Nx == 1) then
    equidistant = .true.
    dt = x(1)
  else if (Nx /= N) then
    stop 'Lengths of the arrays x and y have to match or length of x must be 1'
  end if

  if (equidistant) then
    !
    ! implementation with equidistant step
    !
    do j = 1, Nd
        yy(1) = y(1,j)*dt/2.0 ! fake first
        do i = 2, N-1
            yy(i) = yy(i-1) + y(i,j)*dt
        end do
        ! first
        yy(1) = 0.0d0  ! it must be reset to zero
        ! last
        yy(N) = yy(N-1) + y(N,j)*dt/2.0d0
        yout(:,j) = yy
    end do

  else
    !
    ! implementation with variable step
    !
    stop 'Not implemented'

  end if

  print *, "None"

end subroutine primitive_trp2_dp_none

subroutine integral_trp2_dp_none(x,y,yout)
  use iso_c_binding
  implicit none
  !
  ! Real double precision version of the routine, accelerate by OpenACC directives
  !
  real(c_double), dimension(:), intent(in) :: x
  real(c_double), dimension(:,:), intent(in) :: y
  real(c_double), dimension(size(y,2)), intent(out) :: yout

  ! local variables
  integer :: N
  integer :: Nx, Nd, i, j
  logical :: equidistant
  real(c_double) :: dt

  N = size(y,1)
  Nd = size(y,2)
  Nx = size(x,1)

  equidistant = .false.
  if (Nx == 1) then
    equidistant = .true.
    dt = x(1)
  else if (Nx /= N) then
    stop 'Lengths of the arrays x and y have to match or length of x must be 1'
  end if

  if (equidistant) then
    !
    ! implementation with equidistant step
    !
    do j = 1, Nd
        yout(j) = y(1,j)*dt/2.0 ! fake first
        yout(j) = yout(j) + sum(y(2:(N-1),j))*dt
        yout(j) = yout(j) + y(N,j)*dt/2.0d0
    end do


  else
    !
    ! implementation with variable step
    !
    stop 'Not implemented'

  end if

end subroutine integral_trp2_dp_none

!-------------------------------------------------------------------------------
!
!
!   OpenACC
!
!
!-------------------------------------------------------------------------------
subroutine primitive_trp2_dp_openacc(x,y)
  use iso_c_binding
  implicit none
  !
  ! Real double precision version of the routine, accelerate by OpenACC directives
  !
  real(c_double), dimension(:), intent(in) :: x
  real(c_double), dimension(:), intent(inout) :: y

  ! local variables
  integer :: N
  integer :: Nx, i
  logical :: equidistant
  real(c_double), dimension(size(y,1)) :: yy
  real(c_double) :: dt

  N = size(y,1)
  Nx = size(x,1)

  equidistant = .false.
  if (Nx == 1) then
    equidistant = .true.
    dt = x(1)
  else if (Nx /= N) then
    stop 'Lengths of the arrays x and y have to match or length of x must be 1'
  end if
  
  if (equidistant) then
    !
    ! implementation with equidistant step
    !
    yy(1) = y(1)*dt/2.0 ! fake first
    do i = 2, N-1
      yy(i) = yy(i-1) + y(i)*dt
    end do
    ! fistr
    yy(1) = 0.0d0  ! it must be reset to zero
    ! last
    yy(N) = yy(N-1) + y(i)*dt/2.0d0


  else
    !
    ! implementation with variable step
    !
    stop 'Not implemented'

  end if

  y = yy
  
end subroutine primitive_trp2_dp_openacc

!-------------------------------------------------------------------------------
!
!
!   OpenMP
!
!
!-------------------------------------------------------------------------------
subroutine primitive_trp2_dp_openmp(x,y,yout)
  use iso_c_binding
  use omp_lib
  implicit none
  !
  ! Real double precision version of the routine, accelerate by OpenMP directives
  !
  real(c_double), dimension(:), intent(in) :: x
  real(c_double), dimension(:,:), intent(in) :: y
  real(c_double), dimension(:,:), intent(out) :: yout

  ! local variables
  integer :: N, Nd
  integer :: Nx, i, j
  logical :: equidistant
  real(c_double), dimension(size(y,1)) :: yy
  real(c_double), dimension(size(y,1), size(y,2)) :: yy2
  real(c_double) :: yini
  real(c_double) :: dt
  integer :: nthreads, n_perthread, n_remainder, thread_id, st, fi

  N = size(y,1)
  Nd = size(y,2)
  Nx = size(x,1)

  equidistant = .false.
  if (Nx == 1) then
    equidistant = .true.
    dt = x(1)
  else if (Nx /= N) then
    stop 'Lengths of the arrays x and y have to match or length of x must be 1'
  end if
  
  if (equidistant) then
      !
      ! implementation with equidistant step
      !

!      !$omp parallel &
!      !$omp default(none) &
!      !$omp private(i,j,nthreads,thread_id,n_perthread,n_remainder,st,fi,yy) &
!      !$omp shared(N, Nd, y, yout, dt)
!      nthreads = omp_get_num_threads()
!      thread_id = omp_get_thread_num()
!
!      n_perthread = Nd/nthreads
!      n_remainder = Nd - nthreads*n_perthread
!
!      st = thread_id*(n_perthread) + 1
!      fi = (thread_id+1)*(n_perthread)
!      do j = st, fi
      !$omp parallel do
      do j = 1, Nd

          yy(1) = y(1,j)*dt/2.0 ! fake first
          do i = 2, N-1
              yy(i) = yy(i-1) + y(i,j)*dt
          end do
          ! first
          yy(1) = 0.0d0  ! it must be reset to zero
          ! last
          yy(N) = yy(N-1) + y(N,j)*dt/2.0d0

          yout(:,j) = yy
      end do
      !$omp end parallel do

  else
    !
    ! implementation with variable step
    !
    stop 'Not implemented'

  end if

  print *, "OpenMP"
  
end subroutine primitive_trp2_dp_openmp

subroutine integral_trp2_dp_openmp(x,y,yout)
  use iso_c_binding
  use omp_lib
  implicit none
  !
  ! Real double precision version of the routine, accelerate by OpenMP directives
  !
  real(c_double), dimension(:), intent(in) :: x
  real(c_double), dimension(:,:), intent(in) :: y
  real(c_double), dimension(size(y,2)), intent(out) :: yout

  ! local variables
  integer :: N
  integer :: Nx, Nd, i, j
  logical :: equidistant
  real(c_double) :: dt
  integer :: nthreads, n_perthread, n_remainder, thread_id, st, fi

  N = size(y,1)
  Nd = size(y,2)
  Nx = size(x,1)

  equidistant = .false.
  if (Nx == 1) then
    equidistant = .true.
    dt = x(1)
  else if (Nx /= N) then
    stop 'Lengths of the arrays x and y have to match or length of x must be 1'
  end if

  if (equidistant) then
    !
    ! implementation with equidistant step
    !
!    !$omp parallel &
!    !$omp default(none) &
!    !$omp private(i,j,nthreads,thread_id,n_perthread,n_remainder,st,fi) &
!    !$omp shared(N, Nd, y, yout, dt)
!    nthreads = omp_get_num_threads()
!    thread_id = omp_get_thread_num()
!
!    n_perthread = Nd/nthreads
!    n_remainder = Nd - nthreads*n_perthread
!
!    st = thread_id*(n_perthread) + 1
!    fi = (thread_id+1)*(n_perthread)
!    do j = st, fi
    !$omp parallel do
    do j = 1, Nd
        yout(j) = y(1,j)*dt/2.0 ! fake first
        yout(j) = yout(j) + sum(y(2:(N-1),j))*dt
        yout(j) = yout(j) + y(N,j)*dt/2.0d0
    end do
    !$omp end parallel do


  else
    !
    ! implementation with variable step
    !
    stop 'Not implemented'

  end if
end subroutine integral_trp2_dp_openmp
