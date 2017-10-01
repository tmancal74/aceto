module acetolab

  use acetodef

  integer, parameter :: FOUR_WAVE_MIXING = 4
      
  type lab_settings
  
    real(dp), dimension(:,:), allocatable :: pol
    real(dp), dimension(3)   :: orient_aver
    integer :: exptype
        
  contains
  
    procedure :: init !=> lab_init
    procedure :: set_laser_polarizations !=> lab_set_laser_polarizations
    procedure :: get_oafactor !=> lab_oafactor
    
  end type

!  private :: lab_init, lab_set_laser_polarizations
!  private :: lab_oafactor

contains

  subroutine init(this, exptype)
!  subroutine lab_init(this, exptype)
    class(lab_settings) :: this
    integer :: exptype
    
    this%exptype = exptype
  
!  end subroutine lab_init
  end subroutine init

  
  subroutine set_laser_polarizations(this, e1, e2, e3, e4)
!  subroutine lab_set_laser_polarizations(this, e1, e2, e3, e4)
    class(lab_settings) :: this
    real(dp), dimension(:) :: e1, e2, e3, e4
    ! local
    real(dp), dimension(3,3) :: M4
    real(dp), dimension(3)   :: F4
    integer :: i
    
    if (.not.allocated(this%pol)) then
       allocate(this%pol(4,3))
    end if
    this%pol(1,:) = e1
    this%pol(2,:) = e2
    this%pol(3,:) = e3
    this%pol(4,:) = e4    
    
    ! initialiying M4
    M4 = -1.0d0
    do i = 1, 3
      M4(i,i) = 4.0d0
    end do
    M4 = M4/30.0d0
    
    F4(1) = dot_product(e4,e3)*dot_product(e2,e1)
    F4(2) = dot_product(e4,e2)*dot_product(e3,e1)
    F4(3) = dot_product(e4,e1)*dot_product(e3,e2)
    
    this%orient_aver = matmul(transpose(M4),F4) 
    
!  end subroutine lab_set_laser_polarizations
  end subroutine set_laser_polarizations
    
  
  function get_oafactor(this, d1, d2, d3, d4) result(ret)
!  function lab_oafactor(this, d1, d2, d3, d4) result(ret)
    class(lab_settings) :: this
    real(dp), dimension(:) :: d1, d2, d3, d4
    real(dp) :: ret
    ! local
    real(dp), dimension(3) :: F4
        
    F4(1) = dot_product(d4,d3)*dot_product(d2,d1)
    F4(2) = dot_product(d4,d2)*dot_product(d3,d1)
    F4(3) = dot_product(d4,d1)*dot_product(d3,d2)
    
    ret = dot_product(this%orient_aver,F4)
                      
!  end function lab_oafactor
  end function get_oafactor
    
end module acetolab