!************************************************************************************
!
! Accelerated Charge and Energy Transfer Objects (ACETO) library
!
! Fortran Implementation
!
! Authors:
!   Tomas Mancal, Charles University
!   
!
!   email: mancal@karlov.mff.cuni.cz
!
!
!************************************************************************************
module acetolib
  use iso_c_binding
  use acetodef
  
  implicit none



  type aceto_properties
    character(64) :: accelerator 
  end type

  type(aceto_properties), parameter :: GPU_CUDA = aceto_properties("gpu_cuda")
  type(aceto_properties), parameter :: GPU_OPENACC = aceto_properties("openacc")
  type(aceto_properties), parameter :: CPU_MULTICORE = aceto_properties("openmp")
  type(aceto_properties), parameter :: NO_ACCELERATION = aceto_properties("no_acceleration")


  type aceto
    character(64) :: name
    procedure(primitive_trp2_dp), pointer, nopass :: primitive_trp2_dp => null()
    procedure(integral_trp2_dp), pointer, nopass :: integral_trp2_dp => null()
    procedure(primitive_trp2_cdp), pointer, nopass :: primitive_trp2_cdp => null()
  contains
    procedure :: init => aceto_init
    procedure :: destroy => aceto_destroy
  end type

  !----------------------------------------------------------------------------------
  !
  ! Abstract interfaces
  !
  !----------------------------------------------------------------------------------
  interface 
    subroutine primitive_trp2_dp(x, y, yout)
      use iso_c_binding
      real(c_double), dimension(:), intent(in) :: x
      real(c_double), dimension(:,:), intent(in) :: y
      real(c_double), dimension(:,:), intent(out) :: yout
    end subroutine primitive_trp2_dp
    subroutine integral_trp2_dp(x, y, yout)
      use iso_c_binding
      real(c_double), dimension(:), intent(in) :: x
      real(c_double), dimension(:,:), intent(in) :: y
      real(c_double), dimension(size(y,2)), intent(out) :: yout
    end subroutine integral_trp2_dp
    subroutine primitive_trp2_cdp(x, y, yout)
      use iso_c_binding
      real(c_double), dimension(:), intent(in) :: x
      complex(c_double_complex), dimension(:), intent(in) :: y
      complex(c_double_complex), dimension(:,:), intent(out) :: yout
    end subroutine primitive_trp2_cdp
  end interface

  !----------------------------------------------------------------------------------
  !
  ! Interfaces to overloaded subroutines
  ! 
  !——--------------------------------------------------------------------------------
  interface primitive_trp2_openacc
    procedure primitive_trp2_dp_openacc
  end interface


  !—---------------------------------------------------------------------------------
  !
  ! Interfaces to external subroutines
  !
  !----------------------------------------------------------------------------------
  interface 

    subroutine primitive_trp2_dp_openacc(x, y, yout)
      use iso_c_binding
      real(c_double), dimension(:), intent(in) :: x
      real(c_double), dimension(:,:), intent(in) :: y
      real(c_double), dimension(:,:), intent(out) :: yout
    end subroutine primitive_trp2_dp_openacc



    subroutine primitive_trp2_dp_openmp(x, y, yout)
      use iso_c_binding
      real(c_double), dimension(:), intent(in) :: x
      real(c_double), dimension(:,:), intent(in) :: y
      real(c_double), dimension(:,:), intent(out) :: yout
    end subroutine primitive_trp2_dp_openmp

    subroutine integral_trp2_dp_openmp(x, y, yout)
      use iso_c_binding
      real(c_double), dimension(:), intent(in) :: x
      real(c_double), dimension(:,:), intent(in) :: y
      real(c_double), dimension(size(y,2)), intent(out) :: yout
    end subroutine integral_trp2_dp_openmp



    subroutine primitive_trp2_dp_none(x, y, yout)
      use iso_c_binding
      real(c_double), dimension(:), intent(in) :: x
      real(c_double), dimension(:,:), intent(in) :: y
      real(c_double), dimension(:,:), intent(out) :: yout
    end subroutine primitive_trp2_dp_none

    subroutine integral_trp2_dp_none(x, y, yout)
      use iso_c_binding
      real(c_double), dimension(:), intent(in) :: x
      real(c_double), dimension(:,:), intent(in) :: y
      real(c_double), dimension(size(y,2)), intent(out) :: yout
    end subroutine integral_trp2_dp_none

  end interface

  ! private names
  private :: aceto_init
  private :: aceto_destroy

contains


  subroutine aceto_init(this, properties)
    !
    ! Initialization of the aceto module
    !
    class(aceto) :: this
    type(aceto_properties) :: properties

    ! choose subroutines to be called
    select case (properties%accelerator)

    case ("openacc")
      !
      ! OpenACC accelerated routines
      !
      this%primitive_trp2_dp => primitive_trp2_dp_openacc
      this%integral_trp2_dp => integral_trp2_dp_none

    case ("openmp")
      !
      ! OpenMP accelerated routines
      !
      this%primitive_trp2_dp => primitive_trp2_dp_openmp
      this%integral_trp2_dp => integral_trp2_dp_openmp

    case ("no_acceleration")
      !
      ! Routines without acceleration
      !
      this%primitive_trp2_dp => primitive_trp2_dp_none
      this%integral_trp2_dp => integral_trp2_dp_none

    case default
      stop "Unknown accelerator "
    end select

  end subroutine aceto_init

  subroutine aceto_destroy(this)
    !
    ! Closing aceto module and releasing resources 
    !
    class(aceto) :: this

  end subroutine aceto_destroy


end module acetolib
