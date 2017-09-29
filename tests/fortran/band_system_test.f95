program band_system_test
  
  use acetodef
  use acetolib
  use acetosys
  use acetolab

  use nr3td
    
  implicit none

  integer, dimension(3) :: Ns
  type(band_system) :: bs
  real(dp), dimension(3,1,2) :: dge    
  real(dp), dimension(3,2,1) :: def
  real(dp), dimension(2,2) :: Kr1
  real(dp), dimension(4) :: en

  real(dp), dimension(3) :: e1    
  type(lab_settings) :: LAB    
  
  complex(dpc), dimension(:,:), allocatable :: resp
  real(dp), dimension(:), allocatable :: t1, t3
  real(dp) :: rwa, rmin
  integer :: i
      
  print *, ""
  print *, "Testing Band System class"

  Ns(1) = 1
  Ns(2) = 2
  Ns(3) = 1
  call bs%init(3,Ns)
  
  en(1) = 0.0
  en(2) = -0.9
  en(3) = 1.1
  en(4) = 2.0
  
  call bs%set_energies(en)
  print *, "Energies"
  print *, bs%en
  print *, "Frequencies"
  print *, bs%om01
  print *, bs%om12

  
  dge = 1.0
  def = 2.0
  
  call bs%set_dipoles(1,2,dge)
  call bs%set_dipoles(2,3,def)
  
  print *, "Transition dipole moments"
  print *, bs%dd01
  print *, bs%dd12
  
  print *, bs%nn01(1,1,2)
  print *, bs%nn12(1,2,1)
  
  Kr1 = 0.1d0
  call bs%set_relaxation_rates(2,Kr1)
  
  print *, bs%Kr11
  print *, bs%Kd01(1,2)
  
  ! Initialization of a laboratory polarizations
  e1 = 0.0d0
  e1(1) = 1.0d0
  call LAB%init(FOUR_WAVE_MIXING)
  call LAB%set_laser_polarizations(e1, e1, e1, e1)

  print *, "Orientation vector"
  print *, LAB%orient_aver
  
  print *, bs%dd01(1,2)    
  
  allocate(resp(10,10))
  allocate(t1(10),t3(10))
  
  resp = 0.0d0 !1.0d0*j1
  
  do i = 1, 10
    t1(i) = (i-1)*0.1
    t3(i) = t1(i)
  end do
  
  rwa = 1.0d0
  rmin = 0.01d0
  
  !print *, "Calculation of response R2g"
  !call nr3_r2g(LAB, bs, 20, t1, t3, rwa, rmin, resp)
  
  !print *, resp
  
end program
