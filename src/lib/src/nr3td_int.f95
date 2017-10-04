module nr3td_int

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
  
  
  subroutine nr3_r2gt10_fi(orient_av, Ns, omge, nnge, ddge, Kdge, Kdee, &
                        gofts, ptn, SS1, &
                        it2, t3s, rwa, rmin, resp)
    !
    ! R2g response of an three band multi-level system, for the case
    ! of t1 = 0 (impulsive pump-probe situation)
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
    real(8), dimension(:), intent(in) :: t3s
    real(8), intent(in) :: rwa
    real(8), intent(in) :: rmin
    complex(8), dimension(:), intent(inout) :: resp  
    
  end subroutine nr3_r2gt10_fi  
  
  
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
  
  
  subroutine nr3_r1f_fi(orient_av, Ns, omge, omef, nnge, ddge, nnef, ddef, &
                        Kdge, Kdee, Kdef, &
                        gofts, ptn, SS1, SS2, &
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
    real(8), dimension(:,:) :: omef        
    real(8), dimension(:,:,:) :: nnge
    real(8), dimension(:,:) :: ddge
    real(8), dimension(:,:,:) :: nnef
    real(8), dimension(:,:) :: ddef        
    real(8), dimension(:,:) :: Kdge
    real(8), dimension(:,:) :: Kdee
    real(8), dimension(:,:) :: Kdef
    complex(8), dimension(:,:) :: gofts
    integer, dimension(:,:) :: ptn
    real(8), dimension(:,:) :: SS1, SS2    
        
    ! original arguments           
    integer, intent(in) :: it2
    real(8), dimension(:), intent(in) :: t1s,t3s
    real(8), intent(in) :: rwa
    real(8), intent(in) :: rmin
    complex(8), dimension(:,:), intent(inout) :: resp  

  end subroutine nr3_r1f_fi


end interface




end module nr3td_int
