module acetosys
!
!
!
!
!
!

use acetodef

type band_system
  !
  ! System of level bands. Currently, only 2 and 3 band systems can be
  ! correctly handled
  !

  ! Number of bands
  integer :: Nb
  ! total number of states
  integer :: Ne
  ! Number of states in bands
  integer, dimension(:), pointer :: Ns => null()
  ! energies of states
  real(dp), dimension(:), pointer :: en => null()
  ! transition frequencies
  real(dp), dimension(:,:), pointer :: om01 => null()
  ! transition frequencies
  real(dp), dimension(:,:), pointer :: om12 => null()
  ! transition dipole moments between bands
  real(dp), dimension(:,:,:), pointer :: nn01 => null()
  ! transition dipole moments between bands
  real(dp), dimension(:,:,:), pointer :: nn12 => null()
  ! transition dipole moment lengths
  real(dp), dimension(:,:), pointer :: dd01 => null()
  ! transition dipole moment lengths
  real(dp), dimension(:,:), pointer :: dd12 => null()
  ! relaxation rates
  real(dp), dimension(:,:), pointer :: Kr11 => null()
  ! relaxation rates
  real(dp), dimension(:,:), pointer :: Kr22 => null()
  ! dephasing rates
  real(dp), dimension(:,:), pointer :: Kd01 => null()
  ! dephasing rates
  real(dp), dimension(:,:), pointer :: Kd11 => null()
  ! dephasing rates
  real(dp), dimension(:,:), pointer :: Kd12 => null()
      
contains

  procedure :: init => bs_init
  procedure :: set_energies => bs_set_energies
  procedure :: set_dipoles => bs_set_dipoles
  procedure :: set_relaxation_rates => bs_set_relaxation_rates
  procedure :: update_dephasing_rates => bs_update_dephasing_rates
  procedure :: delete_dephasing_rates => bs_delete_dephasing_rates
  
end type band_system

! all procedures below are private
private :: bs_init, bs_set_energies
private :: bs_set_dipoles, bs_set_relaxation_rates
private :: bs_update_dephasing_rates
private :: bs_delete_dephasing_rates

contains

  subroutine bs_init(this, Nb, Ns)
    ! Initializes the band_system object
    !
    ! This accepts all posible values of Nb and Ns
    ! but below, only Nb = 2 and 3 are valid
    !
    class(band_system) :: this
    integer :: Nb
    integer, dimension(:), target :: Ns
    
    if (Nb > 3) stop "Only 2 and 3 band systems are currently supported"
    
    ! set number of bands
    this%Nb = Nb
    ! set number of states in bands
    this%Ns => Ns(1:Nb)
    ! calculate total number of states
    this%Ne = sum(this%Ns)    
    
  end subroutine bs_init


  subroutine bs_set_energies(this, en)
    ! sets energies of the band_system
    !
    !
    !
    class(band_system) :: this
    real(dp), dimension(:), target :: en
          
    this%en => en(1:this%Ne)
              
    if (.not.associated(this%om01)) allocate(this%om01(this%Ns(1),this%Ns(2)))
    if (.not.associated(this%om12)) allocate(this%om12(this%Ns(2),this%Ns(3)))
    
    do i = 1, this%Ns(1)
      do j = 1, this%Ns(2)
        this%om01(i,j) = this%en(i)-this%en(this%Ns(1)+j)
      end do
    end do
    do i = 1, this%Ns(2)
      do j = 1, this%Ns(3)
        this%om12(i,j) = this%en(this%Ns(1)+i)-this%en(this%Ns(1)+this%Ns(2)+j)
      end do
    end do
    
  end subroutine bs_set_energies



  subroutine bs_set_dipoles(this, Nbi, Nbf, dab)
    ! sets transition dipole moments between two bands
    !
    !
    !
    !
    class(band_system) :: this
    integer :: Nbi, Nbf
    real(dp), dimension(:,:,:), target :: dab
    ! local
    integer :: N1, N2
    integer :: i, j
    real(dp) :: dd
    
    if (size(dab,1) /= 3) stop "Wrong size of dimension 1"
    if (size(dab,2) /= this%Ns(Nbi)) stop "Wrong size of dimension 2"
    if (size(dab,3) /= this%Ns(Nbf)) stop "Wrong size of dimension 3"
    
    if ((Nbi == 1) .and. (Nbf == 2)) then
      ! assign dipole matrix
      this%nn01 => dab
      ! calculate dipole lengths
      if (.not.associated(this%dd01)) then
         allocate(this%dd01(this%Ns(Nbi),this%Ns(Nbd)))
      end if
      N1 = this%Ns(Nbi)
      N2 = this%Ns(Nbf)
      do i = 1, N1
        do j = 1, N2 
          dd = sqrt(dot_product(this%nn01(:,i,j),this%nn01(:,i,j)))
          this%dd01(i,j) = dd
          ! normalize the vector
          this%nn01(:,i,j) = this%nn01(:,i,j)/dd          
          end do
        end do

    else if ((Nbi == 2) .and. (Nbf == 3)) then
      ! assign dipole matrix
      this%nn12 => dab
      ! calculate dipole lengths
      if (.not.associated(this%dd12)) then
         allocate(this%dd12(this%Ns(Nbi),this%Ns(Nbd)))
      end if
      N1 = this%Ns(Nbi)
      N2 = this%Ns(Nbf)
      do i = 1, N1
        do j = 1, N2 
          dd = sqrt(dot_product(this%nn12(:,i,j),this%nn12(:,i,j)))
          this%dd12(i,j) = dd
          ! normalize the vector
          this%nn12(:,i,j) = this%nn12(:,i,j)/dd          
          end do
        end do
       
    else
        stop "Attempt to assign unsupported dipole block"
    end if
                   
  end subroutine bs_set_dipoles

  subroutine bs_set_relaxation_rates(this, Nb, RR)
    class(band_system) :: this
    integer :: Nb
    real(dp), dimension(:,:), target :: RR
    
    if (size(RR,1) /= this%Ns(Nb)) stop "Wrong size of dimension 1"
    if (size(RR,2) /= this%Ns(Nb)) stop "Wrong size of dimension 2"
    
    if (Nb == 2) then
       this%Kr11 => RR
    else if (Nb == 3) then
       this%Kr22 => RR
    else
       stop "Attempt to assign unsupported rate block"
    end if
    
    call this%update_dephasing_rates(Nb)
    
  end subroutine bs_set_relaxation_rates
  
  subroutine bs_update_dephasing_rates(this, Nb)
    class(band_system) :: this
    integer :: Nb
    ! local
    integer :: i, j

    ! check existence of dephasing rate matrices
    if (.not. associated(this%Kd01)) then
       allocate(this%Kd01(this%Ns(1),this%Ns(2)))
       this%Kd01 = 0.0d0
    end if
    if (.not. associated(this%Kd12)) then
       allocate(this%Kd12(this%Ns(2),this%Ns(3)))
       this%Kd12 = 0.0d0
    end if
    if (.not. associated(this%Kd11)) then
       allocate(this%Kd11(this%Ns(2),this%Ns(2)))
       this%Kd11 = 0.0d0
    end if
    
    if (Nb == 2) then
       ! dephasing contribution of Kr11 to Kd01, Kd12 and Kd11
       do i = 1, this%Ns(1)
       do j = 1, this%Ns(2)
         this%Kd01(i,j) = this%Kd01(i,j) + this%Kr11(j,j)/2.0d0
       end do
       end do
       do i = 1, this%Ns(2)
       do j = 1, this%Ns(3)
         this%Kd12(i,j) = this%Kd12(i,j) + this%Kr11(i,i)/2.0d0
       end do
       end do
       do i = 1, this%Ns(2)
       do j = 1, this%Ns(2)
         this%Kd11(i,j) = this%Kd11(i,j) + &
                  (this%Kr11(i,i)+this%Kr11(j,j))/2.0d0
       end do
       end do

    else if (Nb == 3) then
       ! dephasing contribution of Kr22 toKd12
       do i = 1, this%Ns(2)
       do j = 1, this%Ns(3)
         this%Kd12(i,j) = this%Kd12(i,j) + this%Kr22(i,i)/2.0d0
       end do
       end do            
    else
       stop "Attempt to update wrong block of dephasing rates"
    end if
    
  end subroutine bs_update_dephasing_rates
  
  subroutine bs_delete_dephasing_rates(this)
    class(band_system) :: this
    ! check existence of dephasing rate matrices
    if (.not. associated(this%Kd01)) then
       allocate(this%Kd01(this%Ns(1),this%Ns(2)))
    end if
    if (.not. associated(this%Kd12)) then
       allocate(this%Kd12(this%Ns(2),this%Ns(3)))
    end if
    if (.not. associated(this%Kd11)) then
       allocate(this%Kd11(this%Ns(2),this%Ns(2)))
    end if
    this%Kd01 = 0.0d0
    this%Kd12 = 0.0d0
    this%Kd11 = 0.0d0
  end subroutine bs_delete_dephasing_rates
    
end module acetosys
