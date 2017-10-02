module nr3td
!
! This module contains interfaces to routines calculating third order
! non-linear response, and it provides some convenience routines 
! with reduced number of arguments
!
!
!
!
use nr3td_int

contains

  subroutine nr3_r2g(LAB, SYS, it2, t1s, t3s, rwa, rmin, resp)
    !
    ! R2g response of an three band multi-level system
    !
    ! Convenience impementation with lab_settings and band_system objects 
    !
    !
    
    use acetodef
    use acetosys
    use acetolab
    
    implicit none
    
    type(lab_settings) :: LAB
    type(band_system)  :: SYS
    integer, intent(in) :: it2
    real(dp), dimension(:), intent(in) :: t1s,t3s
    real(dp), intent(in) :: rwa, rmin
    complex(dpc), dimension(:,:), intent(inout) :: resp
           
           
    call nr3_r2g_fic(LAB%orient_aver, &
                    SYS%Ns, SYS%om01, SYS%nn01, SYS%dd01, SYS%Kd01, SYS%Kd11, &
                    SYS%gofts, SYS%ptn, SYS%SS1, &
                    it2, t1s, t3s, rwa, rmin, resp)
    
  end subroutine nr3_r2g

  subroutine nr3_r2gt10(LAB, SYS, it2, t3s, rwa, rmin, resp)
    !
    ! R2g response of an three band multi-level system, for thw pump-probe
    ! case of t1 = 0
    !
    ! Convenience impementation with lab_settings and band_system objects 
    !
    !
    
    use acetodef
    use acetosys
    use acetolab
    
    implicit none
    
    type(lab_settings) :: LAB
    type(band_system)  :: SYS
    integer, intent(in) :: it2
    real(dp), dimension(:), intent(in) :: t3s
    real(dp), intent(in) :: rwa, rmin
    complex(dpc), dimension(:), intent(inout) :: resp
           
           
    call nr3_r2g_fic(LAB%orient_aver, &
                    SYS%Ns, SYS%om01, SYS%nn01, SYS%dd01, SYS%Kd01, SYS%Kd11, &
                    SYS%gofts, SYS%ptn, SYS%SS1, &
                    it2, t3s, rwa, rmin, resp)
    
  end subroutine nr3_r2gt10


  subroutine nr3_r3g(LAB, SYS, it2, t1s, t3s, rwa, rmin, resp)
    !
    ! R2g response of an three band multi-level system
    !
    ! Convenience impementation with lab_settings and band_system objects 
    !
    !
    
    use acetodef
    use acetosys
    use acetolab
    
    implicit none
    
    type(lab_settings) :: LAB
    type(band_system)  :: SYS
    integer, intent(in) :: it2
    real(dp), dimension(:), intent(in) :: t1s,t3s
    real(dp), intent(in) :: rwa, rmin
    complex(dpc), dimension(:,:), intent(inout) :: resp
           
           
    call nr3_r3g_fic(LAB%orient_aver, &
                    SYS%Ns, SYS%om01, SYS%nn01, SYS%dd01, SYS%Kd01, SYS%Kd11, &
                    SYS%gofts, SYS%ptn, SYS%SS1, &
                    it2, t1s, t3s, rwa, rmin, resp)
    
  end subroutine nr3_r3g


  subroutine nr3_r1g(LAB, SYS, it2, t1s, t3s, rwa, rmin, resp)
    !
    ! R2g response of an three band multi-level system
    !
    ! Convenience impementation with lab_settings and band_system objects 
    !
    !
    
    use acetodef
    use acetosys
    use acetolab
    
    implicit none
    
    type(lab_settings) :: LAB
    type(band_system)  :: SYS
    integer, intent(in) :: it2
    real(dp), dimension(:), intent(in) :: t1s,t3s
    real(dp), intent(in) :: rwa, rmin
    complex(dpc), dimension(:,:), intent(inout) :: resp
           
           
    call nr3_r1g_fic(LAB%orient_aver, &
                    SYS%Ns, SYS%om01, SYS%nn01, SYS%dd01, SYS%Kd01, SYS%Kd11, &
                    SYS%gofts, SYS%ptn, SYS%SS1, &
                    it2, t1s, t3s, rwa, rmin, resp)
    
  end subroutine nr3_r1g

  subroutine nr3_r4g(LAB, SYS, it2, t1s, t3s, rwa, rmin, resp)
    !
    ! R2g response of an three band multi-level system
    !
    ! Convenience impementation with lab_settings and band_system objects 
    !
    !
    
    use acetodef
    use acetosys
    use acetolab
    
    implicit none
    
    type(lab_settings) :: LAB
    type(band_system)  :: SYS
    integer, intent(in) :: it2
    real(dp), dimension(:), intent(in) :: t1s,t3s
    real(dp), intent(in) :: rwa, rmin
    complex(dpc), dimension(:,:), intent(inout) :: resp
           
           
    call nr3_r4g_fic(LAB%orient_aver, &
                    SYS%Ns, SYS%om01, SYS%nn01, SYS%dd01, SYS%Kd01, SYS%Kd11, &
                    SYS%gofts, SYS%ptn, SYS%SS1, &
                    it2, t1s, t3s, rwa, rmin, resp)
    
  end subroutine nr3_r4g

end module nr3td