!
! This file contains routines calculating response functions of third
! order  
!
!


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
    use acetoutils
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
           
    ! local
    integer :: Ng, Ne
    integer :: Nt1, Nt3
    integer :: it1, it3
    integer :: e1, e2, g1, f1
    real(8) :: t1,t2,t3
    real(8), dimension(:,:), allocatable :: oafac
    complex(8) :: r, prod, exparg
           
    real(8), dimension(:,:,:), allocatable :: ss2
        
    ! lineshape function values
    complex(8), dimension(:), allocatable :: gn_t1, gn_t2, gn_t3
    complex(8), dimension(:), allocatable :: gn_t1t2, gn_t2t3, gn_t1t2t3

    complex(8) :: gg_21_t3, gg_21_t1, gg_21_t2, gg_22_t1t2
    complex(8) :: gg_11_t2t3, gg_21_t1t2t3

    ! minimal value of the dipole factor
    real(8) :: minfac    
           
        
    Nt1 = size(t1s)
    Nt3 = size(t3s)
    
    Ng = Ns(1)
    Ne = Ns(2)
    
    allocate(oafac(Ne,Ne))
    allocate(ss2(Ne,Ne,Ne))
    allocate(gn_t1(Ne), gn_t2(Ne), gn_t3(Ne))
    allocate(gn_t1t2(Ne), gn_t2t3(Ne), gn_t1t2t3(Ne))

    call set_goft_mixing(SS1,ss2)

    ! initial and final states (normally there is a sum over them)
    g1 = 1
    f1 = 1 ! this points to ground state
    
    call set_dipole_factor_g(g1, f1, "R2g", orient_av, Ne, ddge, nnge, oafac, &
                           rmin, minfac)

    t2 = t1s(it2)


    do it1 = 1,Nt1
    do it3 = 1,Nt3
    
      r = 0.0d0
      t1 = t1s(it1)
      t3 = t3s(it3)

      call set_goft_g(gn_t1, it1, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t2, it2, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t3, it3, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2, it1+it2, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t2t3, it2+it3, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2t3, it1+it2+it3, Nt1, gofts, ptn, t1s)
      
!      ! assuming Ng = 1
      do e1 = 1, Ne
      do e2 = 1, Ne
                                 
          if (oafac(e1,e2) > minfac) then

              ! frequencies
              exparg = -j1*((omge(1,e1)+rwa)*t1 - (omge(1,e2)+rwa)*t3) &
                       -j1*(omge(1,e2)-omge(1,e1))*t2

              ! dephasing
              exparg = exparg + & 
                (-Kdge(1,e1)*t1 - Kdge(1,e2)*t3)
            
              ! decay
              if (e1 /= e2) then
                  exparg = exparg - Kdee(e2,e1)*t2
              else
                  ! dephasing rates for e1 == e2 contain -depopulation rates
                  exparg = exparg + Kdee(e2,e1)*t2              
              end if

              gg_21_t1 = dot_product(ss2(:,e2,e1),gn_t1)
              gg_21_t2 = dot_product(ss2(:,e2,e1),gn_t2)
              gg_21_t3 = dot_product(ss2(:,e2,e1),gn_t3)
              gg_22_t1t2 = dot_product(ss2(:,e2,e2),gn_t1t2)
              gg_11_t2t3 = dot_product(ss2(:,e1,e1),gn_t2t3)
              gg_21_t1t2t3 = dot_product(ss2(:,e2,e1),gn_t1t2t3)
        
              ! line-shape functions
              exparg = exparg + &
                (-conjg(gg_21_t1) -conjg(gg_21_t3) + gg_21_t2 &
                 - conjg(gg_22_t1t2) -gg_11_t2t3 + conjg(gg_21_t1t2t3))
            
              ! exp
              prod = oafac(e1,e2)*exp(exparg)
          
              r = r + prod
          
          end if
    
      end do
      end do
    
      resp(it1,it3) = resp(it1,it3) + r
    
    end do
    end do

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
    use acetoutils
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
           
    ! local
    integer :: Ng, Ne
    integer :: Nt1, Nt3
    integer :: it1, it3
    integer :: e1, e2, g1, f1
    real(8) :: t1,t2,t3
    real(8), dimension(:,:), allocatable :: oafac
    complex(8) :: r, prod, exparg
           
    real(8), dimension(:,:,:), allocatable :: ss2
        
    ! lineshape function values
    complex(8), dimension(:), allocatable :: gn_t1, gn_t2, gn_t3
    complex(8), dimension(:), allocatable :: gn_t1t2, gn_t2t3, gn_t1t2t3

    complex(8) :: gg_22_t1, gg_21_t2, gg_21_t1t2, gg_21_t2t3, gg_21_t1t2t3
    complex(8) :: gg_11_t3

    ! minimal value of the dipole factor
    real(8) :: minfac    
           
        
    Nt1 = size(t1s)
    Nt3 = size(t3s)
    
    Ng = Ns(1)
    Ne = Ns(2)
    
    allocate(oafac(Ne,Ne))
    allocate(ss2(Ne,Ne,Ne))
    allocate(gn_t1(Ne), gn_t2(Ne), gn_t3(Ne))
    allocate(gn_t1t2(Ne), gn_t2t3(Ne), gn_t1t2t3(Ne))

    call set_goft_mixing(SS1,ss2)

    ! initial and final states (normally there is a sum over them)
    g1 = 1
    f1 = 1 ! this points to ground state
    
    call set_dipole_factor_g(g1, f1, "R3g", orient_av, Ne, ddge, nnge, oafac, &
                           rmin, minfac)

    t2 = t1s(it2)


    do it1 = 1,Nt1
    do it3 = 1,Nt3
    
      r = 0.0d0
      t1 = t1s(it1)
      t3 = t3s(it3)

      call set_goft_g(gn_t1, it1, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t2, it2, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t3, it3, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2, it1+it2, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t2t3, it2+it3, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2t3, it1+it2+it3, Nt1, gofts, ptn, t1s)
      
!      ! assuming Ng = 1
      do e1 = 1, Ne
      do e2 = 1, Ne
                                 
          if (oafac(e1,e2) > minfac) then

              ! frequencies
              exparg = -j1*((omge(1,e1)+rwa)*t1 - (omge(1,e2)+rwa)*t3) 

              ! dephasing
              exparg = exparg + & 
                (-Kdge(1,e1)*t1 - Kdge(1,e2)*t3)
            
              ! decay
              !if (e1 /= e2) then
              !    exparg = exparg + &
              !      (- Kdee(e2,e1)*t2)
              !end if

              gg_22_t1 = dot_product(ss2(:,e2,e2),gn_t1)
              gg_21_t2 = dot_product(ss2(:,e2,e1),gn_t2)
              gg_21_t1t2 = dot_product(ss2(:,e2,e1),gn_t1t2)
              gg_21_t2t3 = dot_product(ss2(:,e2,e1),gn_t2t3)
              gg_21_t1t2t3 = dot_product(ss2(:,e2,e1),gn_t1t2t3)
              gg_11_t3 = dot_product(ss2(:,e1,e1),gn_t3)
        
              ! line-shape functions
              exparg = exparg + &
                (-conjg(gg_22_t1) +conjg(gg_21_t2) -conjg(gg_21_t1t2) &
                 -conjg(gg_21_t2t3) +conjg(gg_21_t1t2t3) -gg_11_t3)
            
              ! exp
              prod = oafac(e1,e2)*exp(exparg)
          
              r = r + prod
          
          end if
    
      end do
      end do
    
      resp(it1,it3) = resp(it1,it3) + r
    
    end do
    end do

end subroutine nr3_r3g_fi

subroutine nr3_r1g_fi(orient_av, Ns, omge, nnge, ddge, Kdge, Kdee, &
                      gofts, ptn, SS1, &
                      it2, t1s, t3s, rwa, rmin, resp)
    !
    ! R1g response of an three band multi-level system
    !
    ! Implementation with full interface
    !
    !
    !
    use acetodef 
    use acetolab
    use acetoutils
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
           
    ! local
    integer :: Ng, Ne
    integer :: Nt1, Nt3
    integer :: it1, it3
    integer :: e1, e2, g1, f1
    real(8) :: t1,t2,t3
    real(8), dimension(:,:), allocatable :: oafac
    complex(8) :: r, prod, exparg
           
    real(8), dimension(:,:,:), allocatable :: ss2
        
    ! lineshape function values
    complex(8), dimension(:), allocatable :: gn_t1, gn_t2, gn_t3
    complex(8), dimension(:), allocatable :: gn_t1t2, gn_t2t3, gn_t1t2t3

    complex(8) :: gg_12_t3, gg_12_t2t3, gg_22_t2
    complex(8) :: gg_11_t1t2t3, gg_12_t1, gg_12_t1t2   

    ! minimal value of the dipole factor
    real(8) :: minfac    
           
        
    Nt1 = size(t1s)
    Nt3 = size(t3s)
    
    Ng = Ns(1)
    Ne = Ns(2)
    
    allocate(oafac(Ne,Ne))
    allocate(ss2(Ne,Ne,Ne))
    allocate(gn_t1(Ne), gn_t2(Ne), gn_t3(Ne))
    allocate(gn_t1t2(Ne), gn_t2t3(Ne), gn_t1t2t3(Ne))

    call set_goft_mixing(SS1,ss2)

    ! initial and final states (normally there is a sum over them)
    g1 = 1
    f1 = 1 ! this points to ground state
    
    call set_dipole_factor_g(g1, f1, "R1g", orient_av, Ne, ddge, nnge, oafac, &
                           rmin, minfac)

    t2 = t1s(it2)


    do it1 = 1,Nt1
    do it3 = 1,Nt3
    
      r = 0.0d0
      t1 = t1s(it1)
      t3 = t3s(it3)

      call set_goft_g(gn_t1, it1, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t2, it2, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t3, it3, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2, it1+it2, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t2t3, it2+it3, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2t3, it1+it2+it3, Nt1, gofts, ptn, t1s)
      
!      ! assuming Ng = 1
      do e1 = 1, Ne
      do e2 = 1, Ne
                                 
          if (oafac(e1,e2) > minfac) then

              ! frequencies
              exparg = j1*((omge(1,e1)+rwa)*t1 + (omge(1,e2)+rwa)*t3) &
                       +j1*(omge(1,e2)-omge(1,e1))*t2

              ! dephasing
              exparg = exparg + & 
                (-Kdge(1,e1)*t1 - Kdge(1,e2)*t3)
            
              ! decay
              if (e1 /= e2) then
                  exparg = exparg - Kdee(e2,e1)*t2
              else
                  ! dephasing rates for e1 == e2 contain -depopulation rates
                  exparg = exparg + Kdee(e2,e1)*t2
              end if

              gg_12_t3 = dot_product(ss2(:,e1,e2),gn_t3)
              gg_12_t2t3 = dot_product(ss2(:,e1,e2),gn_t2t3)
              gg_22_t2 = dot_product(ss2(:,e2,e2),gn_t2)
              gg_11_t1t2t3 = dot_product(ss2(:,e1,e1),gn_t1t2t3)
              gg_12_t1 = dot_product(ss2(:,e1,e2),gn_t1)
              gg_12_t1t2 = dot_product(ss2(:,e1,e2),gn_t1t2)
        
              ! line-shape functions
              exparg = exparg + &
                (-conjg(gg_12_t3) +conjg(gg_12_t2t3) -conjg(gg_22_t2) &
                 -gg_11_t1t2t3 -gg_12_t1 +gg_12_t1t2)
            
              ! exp
              prod = oafac(e1,e2)*exp(exparg)
          
              r = r + prod
          
          end if
    
      end do
      end do
    
      resp(it1,it3) = resp(it1,it3) + r
    
    end do
    end do

end subroutine nr3_r1g_fi

subroutine nr3_r4g_fi(orient_av, Ns, omge, nnge, ddge, Kdge, Kdee, &
                      gofts, ptn, SS1, &
                      it2, t1s, t3s, rwa, rmin, resp)
    !
    ! R1g response of an three band multi-level system
    !
    ! Implementation with full interface
    !
    !
    !
    use acetodef 
    use acetolab
    use acetoutils
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
           
    ! local
    integer :: Ng, Ne
    integer :: Nt1, Nt3
    integer :: it1, it3
    integer :: e1, e2, g1, f1
    real(8) :: t1,t2,t3
    real(8), dimension(:,:), allocatable :: oafac
    complex(8) :: r, prod, exparg
           
    real(8), dimension(:,:,:), allocatable :: ss2
        
    ! lineshape function values
    complex(8), dimension(:), allocatable :: gn_t1, gn_t2, gn_t3
    complex(8), dimension(:), allocatable :: gn_t1t2, gn_t2t3, gn_t1t2t3

    complex(8) :: gg_11_t1, gg_12_t2, gg_12_t1t2 
    complex(8) :: gg_12_t2t3, gg_12_t1t2t3, gg_22_t3   

    ! minimal value of the dipole factor
    real(8) :: minfac    
           
        
    Nt1 = size(t1s)
    Nt3 = size(t3s)
    
    Ng = Ns(1)
    Ne = Ns(2)
    
    allocate(oafac(Ne,Ne))
    allocate(ss2(Ne,Ne,Ne))
    allocate(gn_t1(Ne), gn_t2(Ne), gn_t3(Ne))
    allocate(gn_t1t2(Ne), gn_t2t3(Ne), gn_t1t2t3(Ne))

    call set_goft_mixing(SS1,ss2)

    ! initial and final states (normally there is a sum over them)
    g1 = 1
    f1 = 1 ! this points to ground state
    
    call set_dipole_factor_g(g1, f1, "R4g", orient_av, Ne, ddge, nnge, oafac, &
                           rmin, minfac)

    t2 = t1s(it2)


    do it1 = 1,Nt1
    do it3 = 1,Nt3
    
      r = 0.0d0
      t1 = t1s(it1)
      t3 = t3s(it3)

      call set_goft_g(gn_t1, it1, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t2, it2, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t3, it3, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2, it1+it2, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t2t3, it2+it3, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2t3, it1+it2+it3, Nt1, gofts, ptn, t1s)
      
!      ! assuming Ng = 1
      do e1 = 1, Ne
      do e2 = 1, Ne
                                 
          if (oafac(e1,e2) > minfac) then

              ! frequencies
              exparg = j1*((omge(1,e1)+rwa)*t1 + (omge(1,e2)+rwa)*t3) 
              
              ! dephasing
              exparg = exparg + & 
                (-Kdge(1,e1)*t1 - Kdge(1,e2)*t3)
            
              ! decay
              !if (e1 /= e2) then
              !    exparg = exparg + &
              !      (- Kdee(e2,e1)*t2)
              !end if

              gg_11_t1 = dot_product(ss2(:,e1,e1),gn_t1)
              gg_12_t2 = dot_product(ss2(:,e1,e2),gn_t2)
              gg_12_t1t2 = dot_product(ss2(:,e1,e2),gn_t1t2)
              gg_12_t2t3 = dot_product(ss2(:,e1,e2),gn_t2t3)
              gg_12_t1t2t3 = dot_product(ss2(:,e1,e2),gn_t1t2t3)
              gg_22_t3 = dot_product(ss2(:,e2,e2),gn_t3)
        
              ! line-shape functions
              exparg = exparg + &
                (-gg_11_t1 -gg_12_t2 +gg_12_t1t2 &
                 +gg_12_t2t3 -gg_12_t1t2t3 -gg_22_t3)
            
              ! exp
              prod = oafac(e1,e2)*exp(exparg)
          
              r = r + prod
          
          end if
    
      end do
      end do
    
      resp(it1,it3) = resp(it1,it3) + r
    
    end do
    end do

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
    use acetoutils
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
           
    ! local
    integer :: Ng, Ne
    integer :: Nt1, Nt3
    integer :: it1, it3
    integer :: e1, e2, g1, f1
    real(8) :: t1,t2,t3
    real(8), dimension(:,:), allocatable :: oafac
    complex(8) :: r, prod, exparg
           
    real(8), dimension(:,:,:), allocatable :: zz2
        
    ! lineshape function values
    complex(8), dimension(:), allocatable :: gn_t1, gn_t2, gn_t3
    complex(8), dimension(:), allocatable :: gn_t1t2, gn_t2t3, gn_t1t2t3

    complex(8) :: gg_21_t3, gg_21_t1, gg_21_t2, gg_22_t1t2
    complex(8) :: gg_11_t2t3, gg_21_t1t2t3

    ! minimal value of the dipole factor
    real(8) :: minfac    
           
        
    Nt1 = size(t1s)
    Nt3 = size(t3s)
    
    Ng = Ns(1)
    Ne = Ns(2)
    
    allocate(oafac(Ne,Ne))
    allocate(zz2(Ne,Ne,Ne))
    allocate(gn_t1(Ne), gn_t2(Ne), gn_t3(Ne))
    allocate(gn_t1t2(Ne), gn_t2t3(Ne), gn_t1t2t3(Ne))

    call set_goft_mixing(SS1,zz2)

    ! initial and final states (normally there is a sum over them)
    g1 = 1
    f1 = 1 ! this points to ground state
    
    call set_dipole_factor_f(g1, f1, "R1f", orient_av, Ne, ddge, nnge, oafac, &
                           rmin, minfac)

    t2 = t1s(it2)


    do it1 = 1,Nt1
    do it3 = 1,Nt3
    
      r = 0.0d0
      t1 = t1s(it1)
      t3 = t3s(it3)

      call set_goft_g(gn_t1, it1, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t2, it2, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t3, it3, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2, it1+it2, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t2t3, it2+it3, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2t3, it1+it2+it3, Nt1, gofts, ptn, t1s)
      
!      ! assuming Ng = 1
      do e1 = 1, Ne
      do e2 = 1, Ne
                                 
          if (oafac(e1,e2) > minfac) then

              ! frequencies
              exparg = -j1*((omge(1,e1)+rwa)*t1 - (omef(e2,1)+rwa)*t3) &
                       -j1*(omge(1,e2)-omge(1,e1))*t2

              ! dephasing
              exparg = exparg + & 
                (-Kdge(1,e1)*t1 - Kdef(e2,1)*t3)
            
              ! decay
              if (e1 /= e2) then
                  exparg = exparg - Kdee(e2,e1)*t2
              else
                  ! dephasing rates for e1 == e2 contain -depopulation rates
                  exparg = exparg + Kdee(e2,e1)*t2              
              end if

              gg_21_t1 = dot_product(zz2(:,e2,e1),gn_t1)
              gg_21_t2 = dot_product(zz2(:,e2,e1),gn_t2)
              gg_21_t3 = dot_product(zz2(:,e2,e1),gn_t3)
              gg_22_t1t2 = dot_product(zz2(:,e2,e2),gn_t1t2)
              gg_11_t2t3 = dot_product(zz2(:,e1,e1),gn_t2t3)
              gg_21_t1t2t3 = dot_product(zz2(:,e2,e1),gn_t1t2t3)
        
              ! line-shape functions
              exparg = exparg + &
                (-conjg(gg_21_t1) -conjg(gg_21_t3) + gg_21_t2 &
                 - conjg(gg_22_t1t2) -gg_11_t2t3 + conjg(gg_21_t1t2t3))
            
              ! exp
              prod = oafac(e1,e2)*exp(exparg)
          
              r = r + prod
          
          end if
    
      end do
      end do
    
      resp(it1,it3) = resp(it1,it3) + r
    
    end do
    end do

end subroutine nr3_r1f_fi


!##############################################################################
!
!  Secular Energy Transfer Containing Pathways
!
!
!##############################################################################


subroutine nr3_r2g_trans_fi(orient_av, Ns, omge, nnge, ddge, Kdge, Kdee, &
                      gofts, ptn, SS1, Ueet2, &
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
    use acetoutils
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
    ! population evolution matrix at time t2
    real(8), dimension(:,:) :: Ueet2
        
    ! original arguments           
    integer, intent(in) :: it2
    real(8), dimension(:), intent(in) :: t1s,t3s
    real(8), intent(in) :: rwa
    real(8), intent(in) :: rmin
    complex(8), dimension(:,:), intent(inout) :: resp
           
    ! local
    integer :: Ng, Ne
    integer :: Nt1, Nt3
    integer :: it1, it3
    integer :: e1, e2, g1, f1
    real(8) :: t1,t2,t3
    real(8), dimension(:,:), allocatable :: oafac
    complex(8) :: r, prod, exparg
           
    real(8), dimension(:,:,:), allocatable :: ss2
        
    ! lineshape function values
    complex(8), dimension(:), allocatable :: gn_t1, gn_t2, gn_t3
    complex(8), dimension(:), allocatable :: gn_t1t2, gn_t2t3, gn_t1t2t3

    complex(8) :: gg_21_t3, gg_21_t1, gg_21_t2, gg_22_t1t2
    complex(8) :: gg_11_t2t3, gg_21_t1t2t3

    ! minimal value of the dipole factor
    real(8) :: minfac    
           
        
    Nt1 = size(t1s)
    Nt3 = size(t3s)
    
    Ng = Ns(1)
    Ne = Ns(2)
    
    allocate(oafac(Ne,Ne))
    allocate(ss2(Ne,Ne,Ne))
    allocate(gn_t1(Ne), gn_t2(Ne), gn_t3(Ne))
    allocate(gn_t1t2(Ne), gn_t2t3(Ne), gn_t1t2t3(Ne))

    call set_goft_mixing(SS1,ss2)

    ! initial and final states (normally there is a sum over them)
    g1 = 1
    f1 = 1 ! this points to ground state
    
    call set_dipole_factor_g(g1, f1, "R2g", orient_av, Ne, ddge, nnge, oafac, &
                           rmin, minfac)

    t2 = t1s(it2)


    do it1 = 1,Nt1
    do it3 = 1,Nt3
    
      r = 0.0d0
      t1 = t1s(it1)
      t3 = t3s(it3)

      call set_goft_g(gn_t1, it1, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t2, it2, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t3, it3, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2, it1+it2, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t2t3, it2+it3, Nt1, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2t3, it1+it2+it3, Nt1, gofts, ptn, t1s)
      
!      ! assuming Ng = 1
      do e1 = 1, Ne
      do e2 = 1, Ne
                                 
          if (oafac(e1,e2) > minfac) then

              ! frequencies
              exparg = -j1*((omge(1,e1)+rwa)*t1 - (omge(1,e2)+rwa)*t3) &
                       -j1*(omge(1,e2)-omge(1,e1))*t2

              ! dephasing
              exparg = exparg + & 
                (-Kdge(1,e1)*t1 - Kdge(1,e2)*t3)
            
              ! decay
              if (e1 /= e2) then
                  exparg = exparg - Kdee(e2,e1)*t2
              else
                  ! dephasing rates for e1 == e2 contain -depopulation rates
                  exparg = exparg + Kdee(e2,e1)*t2              
              end if

              gg_21_t1 = dot_product(ss2(:,e2,e1),gn_t1)
              gg_21_t2 = dot_product(ss2(:,e2,e1),gn_t2)
              gg_21_t3 = dot_product(ss2(:,e2,e1),gn_t3)
              gg_22_t1t2 = dot_product(ss2(:,e2,e2),gn_t1t2)
              gg_11_t2t3 = dot_product(ss2(:,e1,e1),gn_t2t3)
              gg_21_t1t2t3 = dot_product(ss2(:,e2,e1),gn_t1t2t3)
        
              ! line-shape functions
              exparg = exparg + &
                (-conjg(gg_21_t1) -conjg(gg_21_t3) + gg_21_t2 &
                 - conjg(gg_22_t1t2) -gg_11_t2t3 + conjg(gg_21_t1t2t3))
            
              ! exp
              prod = oafac(e1,e2)*exp(exparg)
          
              r = r + prod
          
          end if
    
      end do
      end do
    
      resp(it1,it3) = resp(it1,it3) + r
    
    end do
    end do

end subroutine nr3_r2g_trans_fi