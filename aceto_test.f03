!******************************************************************
!
! Test driver for: 
! 
! Accelerated Charge and Energy Transfer Objects (ACETO) library 
!
!
! Created by:
! 
!   Tomas Mancal
!   Charles University (2016)
!   
!   email: mancal@karlov.mff.cuni.cz
!
!
!******************************************************************
program aceto_test
  use iso_c_binding
  use acetolib
  implicit none

  integer(c_long) :: N, N1, Nd
  integer :: i
  real(c_double) :: dt  
  real(c_double), dimension(:), allocatable :: x, x1
  real(c_double), dimension(:,:), allocatable :: y, yout

  real(c_double), dimension(:), allocatable :: yI

  type(aceto) :: acc_openacc, acc_openmp, acc_none
  type(aceto_properties) :: prop

  real :: start, finish
  integer :: test, counti, countf, count_rate

  print *, ""
  print *, "************************************************************"
  print *, "*                Testing Aceto library"
  print *, "************************************************************"
  print *, "" 


  ! initialize aceto defaults
  !prop = GPU_OPENACC
  prop = CPU_MULTICORE
  call acc_openmp%init(prop)
  call acc_openacc%init(GPU_OPENACC)
  call acc_none%init(NO_ACCELERATION)

  !*****************************************************************************
  !
  ! numerical quadrature
  !
  !*****************************************************************************
  print *, "Numerical quadrature "

  Nd = 81*81
  N1 = 1
  N = 1000
  dt = 0.1

  print *, N*Nd/(1024*1024), "MB"
  allocate(x1(1:N1), x(1:N), y(1:N,1:Nd))
  allocate(yout(1:N,1:Nd), yI(1:Nd))

  open(unit=10, file="test.in")

  read(10,*) test

  close(10)



  y = 0.0
  do i = 1, N
      x(i) = (i-1)*dt
  end do
  x1(1) = dt

  select case(test)

      case (0)

          do i = 1, Nd
              y(:,i) = sin(100.0*x/real(Nd))
          end do

          call system_clock(counti, count_rate)
          call cpu_time(start)
          call acc_none%primitive_trp2_dp(x1,y,yout)
          call cpu_time(finish)
          call system_clock(countf)
          print '("Time NO_ACCELERATION = ",f12.9," seconds.")',finish-start
          print *, yout(N,Nd)

      case (1)

          print *, "nic"
      !  y = sin(x)
      !
      !  call cpu_time(start)
      !  call acc_openacc%primitive_trp2_dp(x1,y)
      !  call cpu_time(finish)
      !
      !  print '("Time GPU_OPENACC = ",f12.9," seconds.")',finish-start
      !  print *, y(N)

      case(2)

          do i = 1, Nd
              y(:,i) = sin(100.0*x/real(Nd))
          end do

          call system_clock(counti, count_rate)
          call cpu_time(start)
          call acc_openmp%primitive_trp2_dp(x1,y,yout)
          call cpu_time(finish)
          call system_clock(countf)
          print '("Time CPU_MULTICORE = ",f12.9," seconds.")',finish-start
          print *, yout(N,Nd)

      case(3)

          do i = 1, Nd
              y(:,i) = sin(100.0*x/real(Nd))
          end do

          call system_clock(counti, count_rate)
          call cpu_time(start)
          call acc_none%integral_trp2_dp(x1,y,yI)
          call cpu_time(finish)
          call system_clock(countf)
          print '("Time NO_ACCELERATION = ",f12.9," seconds.")',finish-start
          print *, yI(Nd)

      case (4)

          do i = 1, Nd
              y(:,i) = sin(100.0*x/real(Nd))
          end do

          call system_clock(counti, count_rate)
          call cpu_time(start)
          call acc_openmp%integral_trp2_dp(x1,y,yI)
          call cpu_time(finish)
          call system_clock(countf)

          print '("Time CPU_MULTICORE = ",f12.9," seconds.")',finish-start
          print *, yI(Nd)

  end select



  print *, "Time = ", real(countf-counti)/real(count_rate)

  ! close aceto
  call acc_openacc%destroy()
  call acc_openmp%destroy()
  call acc_none%destroy()

  end program aceto_test

