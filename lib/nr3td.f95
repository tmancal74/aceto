module nr3td
!
! This module contains interfaces to routines calculating third order
! non-linear response, and it provides some convenience routines 
! with reduced number of arguments
!
!
!
!
interface 

  subroutine nr3_r2g_fi(orient_av, Ns, omge, nnge, ddge, Kdge, Kdee, &
                        gofts, ptn, SS1, &
                        it2, t1s, t3s, rwa, rmin, resp)
    !
    ! R2g response of an three band multi-level system
    !
    ! Implementation with full interface
    !
    !
    !
    use acetodef 
    use acetolab
    implicit none
    
    ! arguments representing lab_settings
    real(8), dimension(:) :: orient_av
        
    ! arguments representing band_system
    integer, dimension(:) :: Ns
    real(8), dimension(:,:) :: omge
    real(8), dimension(:,:,:) :: nnge
    real(8), dimension(:,:) :: ddge
    real(8), dimension(:,:) :: Kdge
    real(8), dimension(:,:) :: Kdee
    complex(8), dimension(:,:) :: gofts
    integer, dimension(:,:) :: ptn
    real(8), dimension(:,:) :: SS1
           
    ! original arguments           
    integer, intent(in) :: it2
    real(8), dimension(:), intent(in) :: t1s,t3s
    real(8), intent(in) :: rwa
    real(8), intent(in) :: rmin
    complex(8), dimension(:,:), intent(inout) :: resp  
  end subroutine nr3_r2g_fi
  
  subroutine nr3_r3g_fi(orient_av, Ns, omge, nnge, ddge, Kdge, Kdee, &
                        gofts, ptn, SS1, &
                        it2, t1s, t3s, rwa, rmin, resp)
    !
    ! R2g response of an three band multi-level system
    !
    ! Implementation with full interface
    !
    !
    !
    use acetodef 
    use acetolab
    implicit none
    
    ! arguments representing lab_settings
    real(8), dimension(:) :: orient_av
        
    ! arguments representing band_system
    integer, dimension(:) :: Ns
    real(8), dimension(:,:) :: omge
    real(8), dimension(:,:,:) :: nnge
    real(8), dimension(:,:) :: ddge
    real(8), dimension(:,:) :: Kdge
    real(8), dimension(:,:) :: Kdee
    complex(8), dimension(:,:) :: gofts
    integer, dimension(:,:) :: ptn
    real(8), dimension(:,:) :: SS1
           
    ! original arguments           
    integer, intent(in) :: it2
    real(8), dimension(:), intent(in) :: t1s,t3s
    real(8) :: rwa, rmin
    complex(8), dimension(:,:), intent(inout) :: resp  
  end subroutine nr3_r3g_fi  

  subroutine nr3_r1g_fi(orient_av, Ns, omge, nnge, ddge, Kdge, Kdee, &
                        gofts, ptn, SS1, &
                        it2, t1s, t3s, rwa, rmin, resp)
    !
    ! R2g response of an three band multi-level system
    !
    ! Implementation with full interface
    !
    !
    !
    use acetodef 
    use acetolab
    implicit none
    
    ! arguments representing lab_settings
    real(8), dimension(:) :: orient_av
        
    ! arguments representing band_system
    integer, dimension(:) :: Ns
    real(8), dimension(:,:) :: omge
    real(8), dimension(:,:,:) :: nnge
    real(8), dimension(:,:) :: ddge
    real(8), dimension(:,:) :: Kdge
    real(8), dimension(:,:) :: Kdee
    complex(8), dimension(:,:) :: gofts
    integer, dimension(:,:) :: ptn
    real(8), dimension(:,:) :: SS1
           
    ! original arguments           
    integer, intent(in) :: it2
    real(8), dimension(:), intent(in) :: t1s,t3s
    real(8) :: rwa, rmin
    complex(8), dimension(:,:), intent(inout) :: resp  
  end subroutine nr3_r1g_fi  


  subroutine nr3_r4g_fi(orient_av, Ns, omge, nnge, ddge, Kdge, Kdee, &
                        gofts, ptn, SS1, &
                        it2, t1s, t3s, rwa, rmin, resp)
    !
    ! R2g response of an three band multi-level system
    !
    ! Implementation with full interface
    !
    !
    !
    use acetodef 
    use acetolab
    implicit none
    
    ! arguments representing lab_settings
    real(8), dimension(:) :: orient_av
        
    ! arguments representing band_system
    integer, dimension(:) :: Ns
    real(8), dimension(:,:) :: omge
    real(8), dimension(:,:,:) :: nnge
    real(8), dimension(:,:) :: ddge
    real(8), dimension(:,:) :: Kdge
    real(8), dimension(:,:) :: Kdee
    complex(8), dimension(:,:) :: gofts
    integer, dimension(:,:) :: ptn
    real(8), dimension(:,:) :: SS1
           
    ! original arguments           
    integer, intent(in) :: it2
    real(8), dimension(:), intent(in) :: t1s,t3s
    real(8) :: rwa, rmin
    complex(8), dimension(:,:), intent(inout) :: resp  
  end subroutine nr3_r4g_fi  

end interface

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
           
           
    call nr3_r2g_fi(LAB%orient_aver, &
                    SYS%Ns, SYS%om01, SYS%nn01, SYS%dd01, SYS%Kd01, SYS%Kd11, &
                    SYS%gofts, SYS%ptn, SYS%SS1, &
                    it2, t1s, t3s, rwa, rmin, resp)
    
  end subroutine nr3_r2g

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
           
           
    call nr3_r3g_fi(LAB%orient_aver, &
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
           
           
    call nr3_r1g_fi(LAB%orient_aver, &
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
           
           
    call nr3_r4g_fi(LAB%orient_aver, &
                    SYS%Ns, SYS%om01, SYS%nn01, SYS%dd01, SYS%Kd01, SYS%Kd11, &
                    SYS%gofts, SYS%ptn, SYS%SS1, &
                    it2, t1s, t3s, rwa, rmin, resp)
    
  end subroutine nr3_r4g

end module nr3td