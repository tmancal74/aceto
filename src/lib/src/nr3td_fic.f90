!
! This file contains routines calculating response functions of the third
! order  
!
!


subroutine nr3_r2g_fic(orient_av, Ns, omge, nnge, ddge, Kdge, Kdee, &
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
    use acetoaux
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
    integer :: Nt1, Nt3, Ntsbi
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
           
        
!    print *, "Inside "
        
    Nt1 = size(t1s)
    Nt3 = size(t3s)
    Ntsbi = size(gofts,2)
    
    Ng = Ns(1)
    Ne = Ns(2)
    
!    print *, "Allocation: ", Ne, Ng
    
    allocate(oafac(Ne,Ne))
    allocate(ss2(Ne,Ne,Ne))
    allocate(gn_t1(Ne), gn_t2(Ne), gn_t3(Ne))
    allocate(gn_t1t2(Ne), gn_t2t3(Ne), gn_t1t2t3(Ne))

    call set_goft_mixing(SS1,ss2)

    ! initial and final states (normally there is a sum over them)
    g1 = 1
    f1 = 1 ! this points to ground state
    
    call set_dipole_factor_g(g1, f1, "R2g0", orient_av, Ne, ddge, nnge, oafac, &
                           rmin, minfac)

    t2 = t1s(it2)


    do it1 = 1,Nt1
    do it3 = 1,Nt3
    
      r = 0.0d0
      t1 = t1s(it1)
      t3 = t3s(it3)

      call set_goft_g(gn_t1, it1, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t2, it2, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t3, it3, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2, it1+it2, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t2t3, it2+it3, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2t3, it1+it2+it3, Ntsbi, gofts, ptn, t1s)
      
      ! assuming Ng = 1
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
    
      resp(it3,it1) = resp(it3,it1) + r
    
    end do
    end do

end subroutine nr3_r2g_fic



subroutine nr3_r2gt10_fic(orient_av, Ns, omge, nnge, ddge, Kdge, Kdee, &
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
    use acetoaux
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
           
    ! local
    integer :: Ng, Ne
    integer :: Nt3, Ntsbi
    integer :: it3
    integer :: e1, e2, g1, f1
    real(8) :: t2,t3
    real(8), dimension(:,:), allocatable :: oafac
    complex(8) :: r, prod, exparg
           
    real(8), dimension(:,:,:), allocatable :: ss2
        
    ! lineshape function values
    complex(8), dimension(:), allocatable :: gn_t2, gn_t3
    complex(8), dimension(:), allocatable :: gn_t2t3

    complex(8) :: gg_21_t3, gg_21_t2, gg_22_t2
    complex(8) :: gg_11_t2t3, gg_21_t2t3
    
    ! minimal value of the dipole factor
    real(8) :: minfac  
    
    real(8) :: t1  
           
        
    Nt3 = size(t3s)
    Ntsbi = size(gofts,2)
    
    Ng = Ns(1)
    Ne = Ns(2)
    
    allocate(oafac(Ne,Ne))
    allocate(ss2(Ne,Ne,Ne))
    allocate(gn_t2(Ne), gn_t3(Ne))
    allocate(gn_t2t3(Ne))

    call set_goft_mixing(SS1,ss2)

    ! initial and final states (normally there is a sum over them)
    g1 = 1
    f1 = 1 ! this points to ground state
    
    call set_dipole_factor_g(g1, f1, "R2g0", orient_av, Ne, ddge, nnge, oafac, &
                           rmin, minfac)

    t2 = t3s(it2)


    !do it1 = 1,Nt1
    t1 = 0.0d0
    
    do it3 = 1,Nt3
    
      r = 0.0d0
      !t1 = t1s(it1)
      t3 = t3s(it3)

      !call set_goft_g(gn_t1, it1, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t2, it2, Ntsbi, gofts, ptn, t3s)
      call set_goft_g(gn_t3, it3, Ntsbi, gofts, ptn, t3s)
      !call set_goft_g(gn_t1t2, it1+it2, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t2t3, it2+it3, Ntsbi, gofts, ptn, t3s)
      !call set_goft_g(gn_t1t2t3, it1+it2+it3, Ntsbi, gofts, ptn, t1s)
      
      ! assuming Ng = 1
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

              !gg_21_t1 = dot_product(ss2(:,e2,e1),gn_t1)
              gg_21_t2 = dot_product(ss2(:,e2,e1),gn_t2)
              gg_21_t3 = dot_product(ss2(:,e2,e1),gn_t3)
              gg_22_t2 = dot_product(ss2(:,e2,e2),gn_t2)
              gg_11_t2t3 = dot_product(ss2(:,e1,e1),gn_t2t3)
              gg_21_t2t3 = dot_product(ss2(:,e2,e1),gn_t2t3)
        
              ! line-shape functions
              !exparg = exparg + &
              !  (-conjg(gg_21_t1) -conjg(gg_21_t3) + gg_21_t2 &
              !   - conjg(gg_22_t1t2) -gg_11_t2t3 + conjg(gg_21_t1t2t3))

              exparg = exparg + &
                ( -conjg(gg_21_t3) + gg_21_t2 &
                 - conjg(gg_22_t2) -gg_11_t2t3 + conjg(gg_21_t2t3))

            
              ! exp
              prod = oafac(e1,e2)*exp(exparg)
          
              r = r + prod
          
          end if
    
      end do
      end do
    
      resp(it3) = resp(it3) + r
    
    end do
    !end do

end subroutine nr3_r2gt10_fic


subroutine nr3_r3g_fic(orient_av, Ns, omge, nnge, ddge, Kdge, Kdee, &
                      gofts, ptn, SS1, &
                      it2, t1s, t3s, rwa, rmin, resp)
    !
    ! R3g response of an three band multi-level system
    !
    ! Implementation with full interface
    !
    !
    !
    use acetodef 
    use acetolab
    use acetoaux
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
    integer :: Nt1, Nt3, Ntsbi
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
    Ntsbi = size(gofts,2)
    
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
    
    call set_dipole_factor_g(g1, f1, "R3g0", orient_av, Ne, ddge, nnge, oafac, &
                           rmin, minfac)

    t2 = t1s(it2)


    do it1 = 1,Nt1
    do it3 = 1,Nt3
    
      r = 0.0d0
      t1 = t1s(it1)
      t3 = t3s(it3)

      call set_goft_g(gn_t1, it1, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t2, it2, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t3, it3, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2, it1+it2, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t2t3, it2+it3, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2t3, it1+it2+it3, Ntsbi, gofts, ptn, t1s)
      
      ! assuming Ng = 1
      do e1 = 1, Ne
      do e2 = 1, Ne
                                 
          !if (oafac(e1,e2) > minfac) then

              ! frequencies
              exparg = -j1*((omge(1,e1)+rwa)*t1 - (omge(1,e2)+rwa)*t3) 

              ! dephasing
              exparg = exparg + & 
                (-Kdge(1,e1)*t1 - Kdge(1,e2)*t3)
            
              ! decay only if it is present in the electronic ground state
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
          
          !end if
    
      end do
      end do
    
      resp(it3,it1) = resp(it3,it1) + r
    
    end do
    end do

end subroutine nr3_r3g_fic

subroutine nr3_r1g_fic(orient_av, Ns, omge, nnge, ddge, Kdge, Kdee, &
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
    use acetoaux
    implicit none
    
    ! arguments representing lab_settings
    real(8), dimension(:), intent(in) :: orient_av
        
    ! arguments representing band_system
    integer, dimension(:), intent(in) :: Ns
    real(8), dimension(:,:), intent(in) :: omge
    real(8), dimension(:,:,:), intent(in) :: nnge
    real(8), dimension(:,:), intent(in) :: ddge
    real(8), dimension(:,:), intent(in) :: Kdge
    real(8), dimension(:,:), intent(in) :: Kdee
    complex(8), dimension(:,:), intent(in) :: gofts
    integer, dimension(:,:), intent(in) :: ptn
    real(8), dimension(:,:), intent(in) :: SS1    
        
    ! original arguments           
    integer, intent(in) :: it2
    real(8), dimension(:), intent(in) :: t1s,t3s
    real(8), intent(in) :: rwa
    real(8), intent(in) :: rmin
    complex(8), dimension(:,:), intent(inout) :: resp
           
    ! local
    integer :: Ng, Ne
    integer :: Nt1, Nt3, Ntsbi
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
    Ntsbi = size(gofts,2)
    
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
    
    call set_dipole_factor_g(g1, f1, "R1g0", orient_av, Ne, ddge, nnge, oafac, &
                           rmin, minfac)

    t2 = t1s(it2)


    do it1 = 1,Nt1
    do it3 = 1,Nt3
    
      r = 0.0d0
      t1 = t1s(it1)
      t3 = t3s(it3)

      call set_goft_g(gn_t1, it1, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t2, it2, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t3, it3, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2, it1+it2, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t2t3, it2+it3, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2t3, it1+it2+it3, Ntsbi, gofts, ptn, t1s)
      
      ! assuming Ng = 1
      do e1 = 1, Ne
      do e2 = 1, Ne
                                 
          !if (oafac(e1,e2) > minfac) then

              ! frequencies
              exparg = j1*((omge(1,e1)+rwa)*t1 + (omge(1,e1)+rwa)*t3) &
                       +j1*(omge(1,e1)-omge(1,e2))*t2

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
                (-conjg(gg_12_t3) + conjg(gg_12_t2t3) -conjg(gg_22_t2) &
                 -gg_11_t1t2t3 -gg_12_t1 +gg_12_t1t2)
            
              ! exp
              prod = oafac(e1,e2)*exp(exparg)

          
              r = r + prod
          
          !end if
    
      end do
      end do
    
      resp(it3,it1) = resp(it3,it1) + r
    
    end do
    end do
    
end subroutine nr3_r1g_fic

subroutine nr3_r4g_fic(orient_av, Ns, omge, nnge, ddge, Kdge, Kdee, &
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
    use acetoaux
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
    integer :: Nt1, Nt3, Ntsbi
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
    Ntsbi = size(gofts,2)
    
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
    
    call set_dipole_factor_g(g1, f1, "R4g0", orient_av, Ne, ddge, nnge, oafac, &
                           rmin, minfac)

    t2 = t1s(it2)


    do it1 = 1,Nt1
    do it3 = 1,Nt3
    
      r = 0.0d0
      t1 = t1s(it1)
      t3 = t3s(it3)

      call set_goft_g(gn_t1, it1, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t2, it2, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t3, it3, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2, it1+it2, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t2t3, it2+it3, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2t3, it1+it2+it3, Ntsbi, gofts, ptn, t1s)
      
      ! assuming Ng = 1
      do e1 = 1, Ne
      do e2 = 1, Ne
                                 
          !if (oafac(e1,e2) > minfac) then

              ! frequencies
              exparg = j1*((omge(1,e1)+rwa)*t1 + (omge(1,e2)+rwa)*t3) 
              
              ! dephasing
              exparg = exparg + & 
                (-Kdge(1,e1)*t1 - Kdge(1,e2)*t3)
            
              ! decay only if it is present in the electronic ground state
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
          
          !end if
    
      end do
      end do
    
      resp(it3,it1) = resp(it3,it1) + r
    
    end do
    end do

end subroutine nr3_r4g_fic



!##############################################################################
!
!    ESA containing pathways
!
!##############################################################################


subroutine nr3_r1fs_fic(orient_av, Ns, omge, omef, nnge, ddge, nnef, ddef, &
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
    use acetoaux
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
    complex(8), dimension(:,:), intent(inout) :: gofts
    integer, dimension(:,:) :: ptn
    real(8), dimension(:,:) :: SS1, SS2    
        
    ! original arguments           
    integer, intent(in) :: it2
    real(8), dimension(:), intent(in) :: t1s,t3s
    real(8), intent(in) :: rwa
    real(8), intent(in) :: rmin
    complex(8), dimension(:,:), intent(inout) :: resp
           
    ! local
    integer :: Ng, Ne, Nf
    integer :: Nt1, Nt3, Ntsbi
    integer :: it1, it3
    integer :: ea, eb, g1, f1
    real(8) :: t1,t2,t3
    real(8), dimension(:,:,:), allocatable :: oafac
    complex(8) :: r, prod, exparg, om
           
    real(8), dimension(:,:,:), allocatable :: zz2
    real(8), dimension(:,:,:), allocatable :: aa2
    real(8), dimension(:,:,:), allocatable :: a21
        
    ! lineshape function values
    complex(8), dimension(:), allocatable :: gn_t1, gn_t2, gn_t3
    complex(8), dimension(:), allocatable :: gn_t1t2, gn_t2t3, gn_t1t2t3

    complex(8) :: gg_aa_t1t2t3, gg_ab_t1, gg_ab_t1t2, gg_af_t1t2 
    complex(8) :: gg_af_t1t2t3, gg_af_t3, gg_ab_t3, gg_ab_t2t3
    complex(8) :: gg_bb_t2, gg_fb_t2, gg_fb_t3
    complex(8) :: gg_fb_t2t3, gg_ff_t3
    
    ! minimal value of the dipole factor
    real(8) :: minfac    
           
    
    Nt1 = size(t1s)
    Nt3 = size(t3s)
    Ntsbi = size(gofts,2)
    
    Ng = Ns(1)
    Ne = Ns(2)
    Nf = Ns(3)
    
    allocate(oafac(Ne,Ne,Nf))
    allocate(zz2(Ne,Ne,Ne))
    allocate(aa2(Ne,Nf,Nf))
    allocate(a21(Ne,Nf,Ne))
    allocate(gn_t1(Ne), gn_t2(Ne), gn_t3(Ne))
    allocate(gn_t1t2(Ne), gn_t2t3(Ne), gn_t1t2t3(Ne))

    call set_goft_mixing(SS1,zz2)
    call set_goft_mixing_22(SS2,aa2,Ne)
    call set_goft_mixing_21(SS2,SS1,a21)

    ! initial and final states (normally there is a sum over them)
    g1 = 1
    !f1 = 1 ! this points to ground state

    do f1 = 1, Nf
      call set_dipole_factor_f(g1, f1, "R1f0", orient_av, Ne, ddge, nnge, ddef, &
                             nnef, oafac(:,:,f1), &
                             rmin, minfac)    
    end do

    t2 = t1s(it2)
    
    do it1 = 1,Nt1
    do it3 = 1,Nt3
    
      r = 0.0d0
      t1 = t1s(it1)
      t3 = t3s(it3)

      call set_goft_g(gn_t1, it1, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t2, it2, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t3, it3, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2, it1+it2, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t2t3, it2+it3, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2t3, it1+it2+it3, Ntsbi, gofts, ptn, t1s)
      
      ! assuming Ng = 1

      do f1 = 1, Nf  
    
      do ea = 1, Ne
      do eb = 1, Ne

                                 
          !if (oafac(ea,ea,f1) > minfac) then

              ! frequencies
              exparg = j1*(omef(ea,f1)+rwa)*t3 - j1*(omge(1,ea)+rwa)*t1 &
                       -j1*(omge(1,ea)-omge(1,eb))*t2


                       
              ! dephasing
              exparg = exparg + & 
                (-Kdge(1,ea)*t1 - Kdef(ea,f1)*t3)
            
              ! decay
              if (ea /= eb) then
                  exparg = exparg - Kdee(ea,eb)*t2
              else
                  ! dephasing rates for e1 == e2 contain -depopulation rates
                  exparg = exparg + Kdee(ea,eb)*t2              
              end if

              gg_aa_t1t2t3 = dot_product(zz2(:,ea,ea),gn_t1t2t3)
              gg_ab_t1 = dot_product(zz2(:,ea,eb),gn_t1)
              gg_ab_t1t2 = dot_product(zz2(:,ea,eb),gn_t1t2)
              gg_af_t1t2 = dot_product(a21(:,f1,ea),gn_t1t2)
              gg_af_t1t2t3 = dot_product(a21(:,f1,ea),gn_t1t2t3)
              gg_ab_t3 = dot_product(zz2(:,ea,eb),gn_t3)
              gg_ab_t2t3 = dot_product(zz2(:,ea,eb),gn_t2t3)
              gg_af_t3 = dot_product(a21(:,f1,ea),gn_t3)
              gg_bb_t2 = dot_product(zz2(:,eb,eb),gn_t2)
              gg_fb_t2 = dot_product(a21(:,f1,eb),gn_t2)
              gg_fb_t3 = dot_product(a21(:,f1,eb),gn_t3)
              gg_fb_t2t3 = dot_product(a21(:,f1,eb),gn_t2t3)
              gg_ff_t3 = dot_product(aa2(:,f1,f1),gn_t3)

!
! Quantathei Symbolic expression
!
! -conjg(gg(a, a, t1 + t2 + t3)) - conjg(gg(a, b, t1)) 
! + conjg(gg(a, b, t1 + t2)) - conjg(gg(a, f, t1 + t2))
! + conjg(gg(a, f, t1 + t2 + t3)) - gg(a, b, t3) + gg(a, b, t2 + t3)
! + gg(a, f, t3) - gg(b, b, t2) + gg(f, b, t2) + gg(f, b, t3)
! - gg(f, b, t2 + t3) - gg(f, f, t3)  
!      
              ! line-shape functions
              exparg = exparg + &
                (-conjg(gg_aa_t1t2t3) -conjg(gg_ab_t1) &
                 + conjg(gg_ab_t1t2) - conjg(gg_af_t1t2) &
                 + conjg(gg_af_t1t2t3) - gg_ab_t3 + gg_ab_t2t3 &
                 + gg_af_t3 - gg_bb_t2 + gg_fb_t2 + gg_fb_t3 &
                 - gg_fb_t2t3 &
                 - gg_ff_t3)
            
              ! exp
              prod = -oafac(ea,eb,f1)*exp(exparg)
          
              r = r + prod
          
          !end if
    
      end do
      end do

      end do
    
      resp(it3,it1) = resp(it3,it1) + r
    
    end do
    end do

    !end do  ! f1

end subroutine nr3_r1fs_fic


subroutine nr3_r2fs_fic(orient_av, Ns, omge, omef, nnge, ddge, nnef, ddef, &
                      Kdge, Kdee, Kdef, &
                      gofts, ptn, SS1, SS2, &
                      it2, t1s, t3s, rwa, rmin, resp)
    !
    ! R2f* response of an three band multi-level system
    !
    ! Implementation with full interface
    !
    !
    !
    use acetodef 
    use acetolab
    use acetoaux
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
    complex(8), dimension(:,:), intent(inout) :: gofts
    integer, dimension(:,:) :: ptn
    real(8), dimension(:,:) :: SS1, SS2    
        
    ! original arguments           
    integer, intent(in) :: it2
    real(8), dimension(:), intent(in) :: t1s,t3s
    real(8), intent(in) :: rwa
    real(8), intent(in) :: rmin
    complex(8), dimension(:,:), intent(inout) :: resp
           
    ! local
    integer :: Ng, Ne, Nf
    integer :: Nt1, Nt3, Ntsbi
    integer :: it1, it3
    integer :: ea, eb, g1, f1
    real(8) :: t1,t2,t3
    real(8), dimension(:,:,:), allocatable :: oafac
    complex(8) :: r, prod, exparg, om
           
    real(8), dimension(:,:,:), allocatable :: zz2
    real(8), dimension(:,:,:), allocatable :: aa2
    real(8), dimension(:,:,:), allocatable :: a21
        
    ! lineshape function values
    complex(8), dimension(:), allocatable :: gn_t1, gn_t2, gn_t3
    complex(8), dimension(:), allocatable :: gn_t1t2, gn_t2t3, gn_t1t2t3

    complex(8) :: gg_ab_t2, gg_bb_t2t3, gg_fb_t2
    complex(8) :: gg_fb_t2t3, gg_aa_t1t2, gg_ab_t1
    complex(8) :: gg_ab_t3, gg_ab_t1t2t3, gg_af_t3, gg_af_t1t2 
    complex(8) :: gg_af_t1t2t3, gg_bf_t3, gg_ff_t3 
    
    ! minimal value of the dipole factor
    real(8) :: minfac      
    
    Nt1 = size(t1s)
    Nt3 = size(t3s)
    Ntsbi = size(gofts,2)
    
    Ng = Ns(1)
    Ne = Ns(2)
    Nf = Ns(3)
    
    allocate(oafac(Ne,Ne,Nf))
    allocate(zz2(Ne,Ne,Ne))
    allocate(aa2(Ne,Nf,Nf))
    allocate(a21(Ne,Nf,Ne))
    allocate(gn_t1(Ne), gn_t2(Ne), gn_t3(Ne))
    allocate(gn_t1t2(Ne), gn_t2t3(Ne), gn_t1t2t3(Ne))

    call set_goft_mixing(SS1,zz2)
    call set_goft_mixing_22(SS2,aa2,Ne)
    call set_goft_mixing_21(SS2,SS1,a21)

    ! initial and final states (normally there is a sum over them)
    g1 = 1
    !f1 = 1 ! this points to ground state

    do f1 = 1, Nf
      call set_dipole_factor_f(g1, f1, "R2f0", orient_av, Ne, ddge, nnge, ddef, &
                             nnef, oafac(:,:,f1), &
                             rmin, minfac)    
    end do

    t2 = t1s(it2)
    

    do it1 = 1,Nt1
    do it3 = 1,Nt3
    
      r = 0.0d0
      t1 = t1s(it1)
      t3 = t3s(it3)

      call set_goft_g(gn_t1, it1, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t2, it2, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t3, it3, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2, it1+it2, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t2t3, it2+it3, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2t3, it1+it2+it3, Ntsbi, gofts, ptn, t1s)
      
      ! assuming Ng = 1
      do f1 = 1, Nf  
    
      do ea = 1, Ne
      do eb = 1, Ne

                                 
          !if (oafac(ea,ea,f1) > minfac) then

              ! frequencies
              exparg = j1*(omef(eb,f1)+rwa)*t3 + j1*(omge(1,ea)+rwa)*t1 &
                       +j1*(omge(1,ea)-omge(1,eb))*t2


                       
              ! dephasing
              exparg = exparg + & 
                (-Kdge(1,ea)*t1 - Kdef(eb,f1)*t3)
            
              ! decay
              if (ea /= eb) then
                  exparg = exparg - Kdee(ea,eb)*t2
              else
                  ! dephasing rates for e1 == e2 contain -depopulation rates
                  exparg = exparg + Kdee(ea,eb)*t2              
              end if

              gg_ab_t2 = dot_product(zz2(:,ea,eb),gn_t2)
              gg_bb_t2t3 = dot_product(zz2(:,eb,eb),gn_t2t3)
              gg_fb_t2 = dot_product(a21(:,f1,eb),gn_t2)
              gg_fb_t2t3 = dot_product(a21(:,f1,eb),gn_t2t3)
              gg_aa_t1t2 = dot_product(zz2(:,ea,ea),gn_t1t2)
              gg_ab_t1 = dot_product(zz2(:,ea,eb),gn_t1)
              gg_ab_t3 = dot_product(zz2(:,ea,eb),gn_t3)
              gg_ab_t1t2t3 = dot_product(zz2(:,ea,eb),gn_t1t2t3)
              gg_af_t3 = dot_product(a21(:,f1,ea),gn_t3)
              gg_af_t1t2 = dot_product(a21(:,f1,ea),gn_t1t2)
              gg_af_t1t2t3 = dot_product(a21(:,f1,ea),gn_t1t2t3)
              gg_bf_t3 = dot_product(a21(:,f1,eb),gn_t3)
              gg_ff_t3 = dot_product(aa2(:,f1,f1),gn_t3)


!
! Quantathei Symbolic expression
!
! conjg(gg(a, b, t2)) - conjg(gg(b, b, t2 + t3)) - conjg(gg(f, b, t2))
! + conjg(gg(f, b, t2 + t3)) - gg(a, a, t1 + t2) - gg(a, b, t1)
! - gg(a, b, t3) + gg(a, b, t1 + t2 + t3) + gg(a, f, t3)
! + gg(a, f, t1 + t2) - gg(a, f, t1 + t2 + t3) + gg(b, f, t3)
! - gg(f, f, t3)

              !      
              ! line-shape functions
              exparg = exparg + &
                (conjg(gg_ab_t2) - conjg(gg_bb_t2t3) &
                 - conjg(gg_fb_t2) + conjg(gg_fb_t2t3) &
                 - gg_aa_t1t2 - gg_ab_t1 - gg_ab_t3 &
                 + gg_ab_t1t2t3 + gg_af_t3 + gg_af_t1t2 - gg_af_t1t2t3 &
                 + gg_bf_t3 - gg_ff_t3)
            
              ! exp
              prod = -oafac(ea,eb,f1)*exp(exparg)
              !prod = exp(exparg)
            
              r = r + prod
          
          !end if
    
      end do
      end do

      end do
    
      resp(it3,it1) = resp(it3,it1) + r
    
    end do
    end do

    !end do  ! f1

end subroutine nr3_r2fs_fic


!##############################################################################
!
!    ESA containing pathways, returned separated
!
!##############################################################################


subroutine nr3_r1fs_list_fic(orient_av, Ns, omge, omef, nnge, ddge, nnef, ddef, &
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
    use acetoaux
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
    complex(8), dimension(:,:), intent(inout) :: gofts
    integer, dimension(:,:) :: ptn
    real(8), dimension(:,:) :: SS1, SS2    
        
    ! original arguments           
    integer, intent(in) :: it2
    real(8), dimension(:), intent(in) :: t1s,t3s
    real(8), intent(in) :: rwa
    real(8), intent(in) :: rmin
    complex(8), dimension(:,:,:), intent(inout) :: resp
           
    ! local
    integer :: Ng, Ne, Nf
    integer :: Nt1, Nt3, Ntsbi
    integer :: it1, it3
    integer :: ea, eb, g1, f1
    real(8) :: t1,t2,t3
    real(8), dimension(:,:,:), allocatable :: oafac
    complex(8) :: r, prod, exparg, om
           
    real(8), dimension(:,:,:), allocatable :: zz2
    real(8), dimension(:,:,:), allocatable :: aa2
    real(8), dimension(:,:,:), allocatable :: a21
        
    ! lineshape function values
    complex(8), dimension(:), allocatable :: gn_t1, gn_t2, gn_t3
    complex(8), dimension(:), allocatable :: gn_t1t2, gn_t2t3, gn_t1t2t3

    complex(8) :: gg_aa_t1t2t3, gg_ab_t1, gg_ab_t1t2, gg_af_t1t2 
    complex(8) :: gg_af_t1t2t3, gg_af_t3, gg_ab_t3, gg_ab_t2t3
    complex(8) :: gg_bb_t2, gg_fb_t2, gg_fb_t3
    complex(8) :: gg_fb_t2t3, gg_ff_t3
    
    ! minimal value of the dipole factor
    real(8) :: minfac    
           
    
    Nt1 = size(t1s)
    Nt3 = size(t3s)
    Ntsbi = size(gofts,2)
    
    Ng = Ns(1)
    Ne = Ns(2)
    Nf = Ns(3)
    
    allocate(oafac(Ne,Ne,Nf))
    allocate(zz2(Ne,Ne,Ne))
    allocate(aa2(Ne,Nf,Nf))
    allocate(a21(Ne,Nf,Ne))
    allocate(gn_t1(Ne), gn_t2(Ne), gn_t3(Ne))
    allocate(gn_t1t2(Ne), gn_t2t3(Ne), gn_t1t2t3(Ne))

    call set_goft_mixing(SS1,zz2)
    call set_goft_mixing_22(SS2,aa2,Ne)
    call set_goft_mixing_21(SS2,SS1,a21)

    ! initial and final states (normally there is a sum over them)
    g1 = 1
    !f1 = 1 ! this points to ground state

    do f1 = 1, Nf
      call set_dipole_factor_f(g1, f1, "R1f0", orient_av, Ne, ddge, nnge, ddef, &
                             nnef, oafac(:,:,f1), &
                             rmin, minfac)    
    end do

    t2 = t1s(it2)
    
    do it1 = 1,Nt1
    do it3 = 1,Nt3
    
      r = 0.0d0
      t1 = t1s(it1)
      t3 = t3s(it3)

      call set_goft_g(gn_t1, it1, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t2, it2, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t3, it3, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2, it1+it2, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t2t3, it2+it3, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2t3, it1+it2+it3, Ntsbi, gofts, ptn, t1s)
      
      ! assuming Ng = 1

      do f1 = 1, Nf  
    
      do ea = 1, Ne
      do eb = 1, Ne

                                 
          !if (oafac(ea,ea,f1) > minfac) then

              ! frequencies
              exparg = j1*(omef(ea,f1)+rwa)*t3 - j1*(omge(1,ea)+rwa)*t1 &
                       -j1*(omge(1,ea)-omge(1,eb))*t2


                       
              ! dephasing
              exparg = exparg + & 
                (-Kdge(1,ea)*t1 - Kdef(ea,f1)*t3)
            
              ! decay
              if (ea /= eb) then
                  exparg = exparg - Kdee(ea,eb)*t2
              else
                  ! dephasing rates for e1 == e2 contain -depopulation rates
                  exparg = exparg + Kdee(ea,eb)*t2              
              end if

              gg_aa_t1t2t3 = dot_product(zz2(:,ea,ea),gn_t1t2t3)
              gg_ab_t1 = dot_product(zz2(:,ea,eb),gn_t1)
              gg_ab_t1t2 = dot_product(zz2(:,ea,eb),gn_t1t2)
              gg_af_t1t2 = dot_product(a21(:,f1,ea),gn_t1t2)
              gg_af_t1t2t3 = dot_product(a21(:,f1,ea),gn_t1t2t3)
              gg_ab_t3 = dot_product(zz2(:,ea,eb),gn_t3)
              gg_ab_t2t3 = dot_product(zz2(:,ea,eb),gn_t2t3)
              gg_af_t3 = dot_product(a21(:,f1,ea),gn_t3)
              gg_bb_t2 = dot_product(zz2(:,eb,eb),gn_t2)
              gg_fb_t2 = dot_product(a21(:,f1,eb),gn_t2)
              gg_fb_t3 = dot_product(a21(:,f1,eb),gn_t3)
              gg_fb_t2t3 = dot_product(a21(:,f1,eb),gn_t2t3)
              gg_ff_t3 = dot_product(aa2(:,f1,f1),gn_t3)

!
! Quantathei Symbolic expression
!
! -conjg(gg(a, a, t1 + t2 + t3)) - conjg(gg(a, b, t1)) 
! + conjg(gg(a, b, t1 + t2)) - conjg(gg(a, f, t1 + t2))
! + conjg(gg(a, f, t1 + t2 + t3)) - gg(a, b, t3) + gg(a, b, t2 + t3)
! + gg(a, f, t3) - gg(b, b, t2) + gg(f, b, t2) + gg(f, b, t3)
! - gg(f, b, t2 + t3) - gg(f, f, t3)  
!      
              ! line-shape functions
              exparg = exparg + &
                (-conjg(gg_aa_t1t2t3) -conjg(gg_ab_t1) &
                 + conjg(gg_ab_t1t2) - conjg(gg_af_t1t2) &
                 + conjg(gg_af_t1t2t3) - gg_ab_t3 + gg_ab_t2t3 &
                 + gg_af_t3 - gg_bb_t2 + gg_fb_t2 + gg_fb_t3 &
                 - gg_fb_t2t3 &
                 - gg_ff_t3)
            
              ! exp
              prod = -oafac(ea,eb,f1)*exp(exparg)
          
              r = r + prod
          
          !end if
    
      end do
      end do

      resp(f1, it3,it1) = resp(f1, it3,it1) + r
      r = 0.0d0

      end do
     
    end do
    end do

    !end do  ! f1

end subroutine nr3_r1fs_list_fic


subroutine nr3_r2fs_list_fic(orient_av, Ns, omge, omef, nnge, ddge, nnef, ddef, &
                      Kdge, Kdee, Kdef, &
                      gofts, ptn, SS1, SS2, &
                      it2, t1s, t3s, rwa, rmin, resp)
    !
    ! R2f* response of an three band multi-level system
    !
    ! Implementation with full interface
    !
    !
    !
    use acetodef 
    use acetolab
    use acetoaux
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
    complex(8), dimension(:,:), intent(inout) :: gofts
    integer, dimension(:,:) :: ptn
    real(8), dimension(:,:) :: SS1, SS2    
        
    ! original arguments           
    integer, intent(in) :: it2
    real(8), dimension(:), intent(in) :: t1s,t3s
    real(8), intent(in) :: rwa
    real(8), intent(in) :: rmin
    complex(8), dimension(:,:,:), intent(inout) :: resp
           
    ! local
    integer :: Ng, Ne, Nf
    integer :: Nt1, Nt3, Ntsbi
    integer :: it1, it3
    integer :: ea, eb, g1, f1
    real(8) :: t1,t2,t3
    real(8), dimension(:,:,:), allocatable :: oafac
    complex(8) :: r, prod, exparg, om
           
    real(8), dimension(:,:,:), allocatable :: zz2
    real(8), dimension(:,:,:), allocatable :: aa2
    real(8), dimension(:,:,:), allocatable :: a21
        
    ! lineshape function values
    complex(8), dimension(:), allocatable :: gn_t1, gn_t2, gn_t3
    complex(8), dimension(:), allocatable :: gn_t1t2, gn_t2t3, gn_t1t2t3

    complex(8) :: gg_ab_t2, gg_bb_t2t3, gg_fb_t2
    complex(8) :: gg_fb_t2t3, gg_aa_t1t2, gg_ab_t1
    complex(8) :: gg_ab_t3, gg_ab_t1t2t3, gg_af_t3, gg_af_t1t2 
    complex(8) :: gg_af_t1t2t3, gg_bf_t3, gg_ff_t3 
    
    ! minimal value of the dipole factor
    real(8) :: minfac      
    
    Nt1 = size(t1s)
    Nt3 = size(t3s)
    Ntsbi = size(gofts,2)
    
    Ng = Ns(1)
    Ne = Ns(2)
    Nf = Ns(3)
    
    allocate(oafac(Ne,Ne,Nf))
    allocate(zz2(Ne,Ne,Ne))
    allocate(aa2(Ne,Nf,Nf))
    allocate(a21(Ne,Nf,Ne))
    allocate(gn_t1(Ne), gn_t2(Ne), gn_t3(Ne))
    allocate(gn_t1t2(Ne), gn_t2t3(Ne), gn_t1t2t3(Ne))

    call set_goft_mixing(SS1,zz2)
    call set_goft_mixing_22(SS2,aa2,Ne)
    call set_goft_mixing_21(SS2,SS1,a21)

    ! initial and final states (normally there is a sum over them)
    g1 = 1
    !f1 = 1 ! this points to ground state

    do f1 = 1, Nf
      call set_dipole_factor_f(g1, f1, "R2f0", orient_av, Ne, ddge, nnge, ddef, &
                             nnef, oafac(:,:,f1), &
                             rmin, minfac)    
    end do

    t2 = t1s(it2)
    

    do it1 = 1,Nt1
    do it3 = 1,Nt3
    
      r = 0.0d0
      t1 = t1s(it1)
      t3 = t3s(it3)

      call set_goft_g(gn_t1, it1, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t2, it2, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t3, it3, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2, it1+it2, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t2t3, it2+it3, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2t3, it1+it2+it3, Ntsbi, gofts, ptn, t1s)
      
      ! assuming Ng = 1
      do f1 = 1, Nf  
    
      do ea = 1, Ne
      do eb = 1, Ne

                                 
          !if (oafac(ea,ea,f1) > minfac) then

              ! frequencies
              exparg = j1*(omef(eb,f1)+rwa)*t3 + j1*(omge(1,ea)+rwa)*t1 &
                       +j1*(omge(1,ea)-omge(1,eb))*t2


                       
              ! dephasing
              exparg = exparg + & 
                (-Kdge(1,ea)*t1 - Kdef(eb,f1)*t3)
            
              ! decay
              if (ea /= eb) then
                  exparg = exparg - Kdee(ea,eb)*t2
              else
                  ! dephasing rates for e1 == e2 contain -depopulation rates
                  exparg = exparg + Kdee(ea,eb)*t2              
              end if

              gg_ab_t2 = dot_product(zz2(:,ea,eb),gn_t2)
              gg_bb_t2t3 = dot_product(zz2(:,eb,eb),gn_t2t3)
              gg_fb_t2 = dot_product(a21(:,f1,eb),gn_t2)
              gg_fb_t2t3 = dot_product(a21(:,f1,eb),gn_t2t3)
              gg_aa_t1t2 = dot_product(zz2(:,ea,ea),gn_t1t2)
              gg_ab_t1 = dot_product(zz2(:,ea,eb),gn_t1)
              gg_ab_t3 = dot_product(zz2(:,ea,eb),gn_t3)
              gg_ab_t1t2t3 = dot_product(zz2(:,ea,eb),gn_t1t2t3)
              gg_af_t3 = dot_product(a21(:,f1,ea),gn_t3)
              gg_af_t1t2 = dot_product(a21(:,f1,ea),gn_t1t2)
              gg_af_t1t2t3 = dot_product(a21(:,f1,ea),gn_t1t2t3)
              gg_bf_t3 = dot_product(a21(:,f1,eb),gn_t3)
              gg_ff_t3 = dot_product(aa2(:,f1,f1),gn_t3)


!
! Quantathei Symbolic expression
!
! conjg(gg(a, b, t2)) - conjg(gg(b, b, t2 + t3)) - conjg(gg(f, b, t2))
! + conjg(gg(f, b, t2 + t3)) - gg(a, a, t1 + t2) - gg(a, b, t1)
! - gg(a, b, t3) + gg(a, b, t1 + t2 + t3) + gg(a, f, t3)
! + gg(a, f, t1 + t2) - gg(a, f, t1 + t2 + t3) + gg(b, f, t3)
! - gg(f, f, t3)

              !      
              ! line-shape functions
              exparg = exparg + &
                (conjg(gg_ab_t2) - conjg(gg_bb_t2t3) &
                 - conjg(gg_fb_t2) + conjg(gg_fb_t2t3) &
                 - gg_aa_t1t2 - gg_ab_t1 - gg_ab_t3 &
                 + gg_ab_t1t2t3 + gg_af_t3 + gg_af_t1t2 - gg_af_t1t2t3 &
                 + gg_bf_t3 - gg_ff_t3)
            
              ! exp
              prod = -oafac(ea,eb,f1)*exp(exparg)
              !prod = exp(exparg)
            
              r = r + prod
          
          !end if
    
      end do
      end do

      resp(f1,it3,it1) = resp(f1,it3,it1) + r
      r = 0.0d0
      
      end do
        
    end do
    end do

    !end do  ! f1

end subroutine nr3_r2fs_list_fic




!##############################################################################
!
!  Secular Markovian Energy Transfer Containing Pathways
!
!
!  The pathways below approximate the lineshape by Markov approximation
!
!
!##############################################################################


subroutine nr3_r2g_trans_fic(orient_av, Ns, omge, nnge, ddge, Kdge, Kdee, &
                      gofts, ptn, SS1, Ueet2, &
                      it2, t1s, t3s, rwa, rmin, resp)
    !
    ! R2g energy transfer response of an three band multi-level system
    !
    ! Implementation with full interface
    !
    !
    !
    use acetodef 
    use acetolab
    use acetoaux
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
    integer :: Nt1, Nt3, Ntsbi
    integer :: it1, it3
    integer :: e1, e2, g1, f1
    real(8) :: t1,t2,t3
    real(8), dimension(:,:), allocatable :: oafac
    complex(8) :: r, prod, exparg
           
    real(8), dimension(:,:,:), allocatable :: ss2
        
    ! lineshape function values
    complex(8), dimension(:), allocatable :: gn_t1, gn_t3

    complex(8) :: gg_11_t1, gg_22_t3

    ! minimal value of the dipole factor
    real(8) :: minfac    
           
        
    Nt1 = size(t1s)
    Nt3 = size(t3s)
    Ntsbi = size(gofts,2)
    
    Ng = Ns(1)
    Ne = Ns(2)
    
    allocate(oafac(Ne,Ne))
    allocate(ss2(Ne,Ne,Ne))
    allocate(gn_t1(Ne), gn_t3(Ne))


    call set_goft_mixing(SS1,ss2)

    ! initial and final states (normally there is a sum over them)
    g1 = 1
    f1 = 1 ! this points to ground state
    
    call set_dipole_factor_gt(g1, f1, "R2gt", orient_av, Ne, ddge, nnge, oafac, &
                           rmin, minfac)

    t2 = t1s(it2)


    do it1 = 1,Nt1
    do it3 = 1,Nt3
    
      r = 0.0d0
      t1 = t1s(it1)
      t3 = t3s(it3)

      call set_goft_g(gn_t1, it1, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t3, it3, Ntsbi, gofts, ptn, t1s)
      
      ! assuming Ng = 1
      do e1 = 1, Ne
      do e2 = 1, Ne
                                 
          !if (oafac(e1,e2) > minfac) then

              ! frequencies
              exparg = -j1*((omge(1,e1)+rwa)*t1 - (omge(1,e2)+rwa)*t3) 

              ! dephasing
              exparg = exparg + & 
                (-Kdge(1,e1)*t1 - Kdge(1,e2)*t3)
              
              gg_11_t1 = dot_product(ss2(:,e1,e1),gn_t1)
              gg_22_t3 = dot_product(ss2(:,e2,e2),gn_t3)
        
              ! line-shape functions
              exparg = exparg - conjg(gg_11_t1) - gg_22_t3 
            
              ! exp
              prod = oafac(e1,e2)*Ueet2(e2,e1)*exp(exparg)
          
              r = r + prod
          
          !end if
    
      end do
      end do
    
      resp(it3,it1) = resp(it3,it1) + r
    
    end do
    end do

end subroutine nr3_r2g_trans_fic

subroutine nr3_r1g_trans_fic(orient_av, Ns, omge, nnge, ddge, Kdge, Kdee, &
                      gofts, ptn, SS1, Ueet2, &
                      it2, t1s, t3s, rwa, rmin, resp)
    !
    ! R2g energy transfer response of an three band multi-level system
    !
    ! Implementation with full interface
    !
    !
    !
    use acetodef 
    use acetolab
    use acetoaux
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
    integer :: Nt1, Nt3, Ntsbi
    integer :: it1, it3
    integer :: e1, e2, g1, f1
    real(8) :: t1,t2,t3
    real(8), dimension(:,:), allocatable :: oafac
    complex(8) :: r, prod, exparg
           
    real(8), dimension(:,:,:), allocatable :: ss2
        
    ! lineshape function values
    complex(8), dimension(:), allocatable :: gn_t1, gn_t3

    complex(8) :: gg_11_t1, gg_22_t3

    ! minimal value of the dipole factor
    real(8) :: minfac    
           
        
    Nt1 = size(t1s)
    Nt3 = size(t3s)
    Ntsbi = size(gofts,2)
    
    Ng = Ns(1)
    Ne = Ns(2)
    
    allocate(oafac(Ne,Ne))
    allocate(ss2(Ne,Ne,Ne))
    allocate(gn_t1(Ne), gn_t3(Ne))


    call set_goft_mixing(SS1,ss2)

    ! initial and final states (normally there is a sum over them)
    g1 = 1
    f1 = 1 ! this points to ground state
    
    call set_dipole_factor_gt(g1, f1, "R1gt", orient_av, Ne, ddge, nnge, oafac, &
                           rmin, minfac)

    t2 = t1s(it2)


    do it1 = 1,Nt1
    do it3 = 1,Nt3
    
      r = 0.0d0
      t1 = t1s(it1)
      t3 = t3s(it3)

      call set_goft_g(gn_t1, it1, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t3, it3, Ntsbi, gofts, ptn, t1s)
      
      ! assuming Ng = 1
      do e1 = 1, Ne
      do e2 = 1, Ne
                                 
          !if (oafac(e1,e2) > minfac) then

              ! frequencies
              exparg = j1*((omge(1,e1)+rwa)*t1 + (omge(1,e2)+rwa)*t3) 

              ! dephasing
              exparg = exparg + & 
                (-Kdge(1,e1)*t1 - Kdge(1,e2)*t3)
              
              gg_11_t1 = dot_product(ss2(:,e1,e1),gn_t1)
              gg_22_t3 = dot_product(ss2(:,e2,e2),gn_t3)
        
              ! line-shape functions
              exparg = exparg - gg_11_t1 - gg_22_t3 
            
              ! exp
              prod = oafac(e1,e2)*Ueet2(e2,e1)*exp(exparg)
          
              r = r + prod
          
          !end if
    
      end do
      end do
    
      resp(it3,it1) = resp(it3,it1) + r
    
    end do
    end do

end subroutine nr3_r1g_trans_fic

subroutine nr3_r1fs_trans_fic(orient_av, Ns, omge, omef, nnge, ddge, nnef, &
                              ddef, Kdge, Kdee, Kdef, &
                              gofts, ptn, SS1, SS2, Ueet2, &
                              it2, t1s, t3s, rwa, rmin, resp)
    !
    ! R1f* energy transfer response of an three band multi-level system
    !
    ! Implementation with full interface
    !
    !
    !
    use acetodef 
    use acetolab
    use acetoaux
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
    complex(8), dimension(:,:), intent(inout) :: gofts
    integer, dimension(:,:) :: ptn
    real(8), dimension(:,:) :: SS1, SS2    
    ! population evolution matrix at time t2
    real(8), dimension(:,:) :: Ueet2
        
    ! original arguments           
    integer, intent(in) :: it2
    real(8), dimension(:), intent(in) :: t1s,t3s
    real(8), intent(in) :: rwa
    real(8), intent(in) :: rmin
    complex(8), dimension(:,:), intent(inout) :: resp
           
    ! local
    integer :: Ng, Ne, Nf
    integer :: Nt1, Nt3, Ntsbi
    integer :: it1, it3
    integer :: ea, eb, g1, f1
    real(8) :: t1,t2,t3
    real(8), dimension(:,:,:), allocatable :: oafac
    complex(8) :: r, prod, exparg, om
           
    real(8), dimension(:,:,:), allocatable :: zz2
    real(8), dimension(:,:,:), allocatable :: aa2
    real(8), dimension(:,:,:), allocatable :: a21
        
    ! lineshape function values
    complex(8), dimension(:), allocatable :: gn_t1, gn_t3

    complex(8) :: gg_aa_t1
    complex(8) :: gg_af_t1t2t3, gg_af_t3, gg_ab_t3, gg_ab_t2t3
    complex(8) :: gg_bb_t3, gg_bf_t3, gg_ff_t3
    
    ! minimal value of the dipole factor
    real(8) :: minfac    
           
    
    Nt1 = size(t1s)
    Nt3 = size(t3s)
    Ntsbi = size(gofts,2)
    
    Ng = Ns(1)
    Ne = Ns(2)
    Nf = Ns(3)
    
    allocate(oafac(Ne,Ne,Nf))
    allocate(zz2(Ne,Ne,Ne))
    allocate(aa2(Ne,Nf,Nf))
    allocate(a21(Ne,Nf,Ne))
    allocate(gn_t1(Ne), gn_t3(Ne))

    call set_goft_mixing(SS1,zz2)
    call set_goft_mixing_22(SS2,aa2,Ne)
    call set_goft_mixing_21(SS2,SS1,a21)

    ! initial and final states (normally there is a sum over them)
    g1 = 1

    do f1 = 1, Nf
      call set_dipole_factor_ft(g1, f1, "R1ft", orient_av, Ne, ddge, nnge, &
                                ddef, nnef, oafac(:,:,f1), rmin, minfac)    
    end do

    t2 = t1s(it2)
    
    do it1 = 1,Nt1
    do it3 = 1,Nt3
    
      r = 0.0d0
      t1 = t1s(it1)
      t3 = t3s(it3)

      call set_goft_g(gn_t1, it1, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t3, it3, Ntsbi, gofts, ptn, t1s)

      
      ! assuming Ng = 1

      do f1 = 1, Nf  
    
      do ea = 1, Ne
      do eb = 1, Ne

                                 
          !if (oafac(ea,ea,f1) > minfac) then

              ! frequencies
              exparg = j1*(omef(eb,f1)+rwa)*t3 - j1*(omge(1,ea)+rwa)*t1 
           
              ! dephasing
              exparg = exparg + & 
                (-Kdge(1,ea)*t1 - Kdef(ea,f1)*t3)
            
              gg_aa_t1 = dot_product(zz2(:,ea,ea),gn_t1)
              gg_bb_t3 = dot_product(zz2(:,eb,eb),gn_t3)
              gg_bf_t3 = dot_product(a21(:,f1,eb),gn_t3)
              gg_ff_t3 = dot_product(aa2(:,f1,f1),gn_t3)
                                  
              ! line-shape functions
              exparg = exparg + &
                (-conjg(gg_bb_t3)  + 2.0*real(gg_bf_t3) - gg_ff_t3) &
                 - conjg(gg_aa_t1)
            
              ! exp
              prod = -oafac(ea,eb,f1)*Ueet2(eb,ea)*exp(exparg)
          
              r = r + prod
          
          !end if
    
      end do
      end do

      end do
    
      resp(it3,it1) = resp(it3,it1) + r
    
    end do
    end do

    !end do  ! f1

end subroutine nr3_r1fs_trans_fic


subroutine nr3_r2fs_trans_fic(orient_av, Ns, omge, omef, nnge, ddge, nnef, &
                              ddef, Kdge, Kdee, Kdef, &
                              gofts, ptn, SS1, SS2, Ueet2, &
                              it2, t1s, t3s, rwa, rmin, resp)
    !
    ! R2f* energy transfer response of an three band multi-level system
    !
    ! Implementation with full interface
    !
    !
    !
    use acetodef 
    use acetolab
    use acetoaux
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
    complex(8), dimension(:,:), intent(inout) :: gofts
    integer, dimension(:,:) :: ptn
    real(8), dimension(:,:) :: SS1, SS2    
    ! population evolution matrix at time t2
    real(8), dimension(:,:) :: Ueet2
        
    ! original arguments           
    integer, intent(in) :: it2
    real(8), dimension(:), intent(in) :: t1s,t3s
    real(8), intent(in) :: rwa
    real(8), intent(in) :: rmin
    complex(8), dimension(:,:), intent(inout) :: resp
           
    ! local
    integer :: Ng, Ne, Nf
    integer :: Nt1, Nt3, Ntsbi
    integer :: it1, it3
    integer :: ea, eb, g1, f1
    real(8) :: t1,t2,t3
    real(8), dimension(:,:,:), allocatable :: oafac
    complex(8) :: r, prod, exparg, om
           
    real(8), dimension(:,:,:), allocatable :: zz2
    real(8), dimension(:,:,:), allocatable :: aa2
    real(8), dimension(:,:,:), allocatable :: a21
        
    ! lineshape function values
    complex(8), dimension(:), allocatable :: gn_t1, gn_t3

    complex(8) :: gg_aa_t1
    complex(8) :: gg_af_t1t2t3, gg_af_t3, gg_ab_t3, gg_ab_t2t3
    complex(8) :: gg_bb_t3, gg_bf_t3, gg_ff_t3
    
    ! minimal value of the dipole factor
    real(8) :: minfac    
           
    
    Nt1 = size(t1s)
    Nt3 = size(t3s)
    Ntsbi = size(gofts,2)
    
    Ng = Ns(1)
    Ne = Ns(2)
    Nf = Ns(3)
    
    allocate(oafac(Ne,Ne,Nf))
    allocate(zz2(Ne,Ne,Ne))
    allocate(aa2(Ne,Nf,Nf))
    allocate(a21(Ne,Nf,Ne))
    allocate(gn_t1(Ne), gn_t3(Ne))

    call set_goft_mixing(SS1,zz2)
    call set_goft_mixing_22(SS2,aa2,Ne)
    call set_goft_mixing_21(SS2,SS1,a21)

    ! initial and final states (normally there is a sum over them)
    g1 = 1

    do f1 = 1, Nf
      call set_dipole_factor_ft(g1, f1, "R2ft", orient_av, Ne, ddge, nnge, &
                                ddef, nnef, oafac(:,:,f1), rmin, minfac)    
    end do

    t2 = t1s(it2)
    
    do it1 = 1,Nt1
    do it3 = 1,Nt3
    
      r = 0.0d0
      t1 = t1s(it1)
      t3 = t3s(it3)

      call set_goft_g(gn_t1, it1, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t3, it3, Ntsbi, gofts, ptn, t1s)

      
      ! assuming Ng = 1

      do f1 = 1, Nf  
    
      do ea = 1, Ne
      do eb = 1, Ne

                                 
          !if (oafac(ea,ea,f1) > minfac) then

              ! frequencies
              exparg = j1*(omef(eb,f1)+rwa)*t3 + j1*(omge(1,ea)+rwa)*t1 
           
              ! dephasing
              exparg = exparg + & 
                (-Kdge(1,ea)*t1 - Kdef(ea,f1)*t3)
            
              gg_aa_t1 = dot_product(zz2(:,ea,ea),gn_t1)
              gg_bb_t3 = dot_product(zz2(:,eb,eb),gn_t3)
              gg_bf_t3 = dot_product(a21(:,f1,eb),gn_t3)
              gg_ff_t3 = dot_product(aa2(:,f1,f1),gn_t3)
                                  
              ! line-shape functions
              exparg = exparg + &
                (-conjg(gg_bb_t3)  + 2.0*real(gg_bf_t3) - gg_ff_t3) &
                 - gg_aa_t1
            
              ! exp
              prod = -oafac(ea,eb,f1)*Ueet2(eb,ea)*exp(exparg)
          
              r = r + prod
          
          !end if
    
      end do
      end do

      end do
    
      resp(it3,it1) = resp(it3,it1) + r
    
    end do
    end do

    !end do  ! f1

end subroutine nr3_r2fs_trans_fic

!##############################################################################
!
!  Non-Markovian Energy Transfer Containing Pathways of n-th Order
!
!
!  
!
!
!##############################################################################

subroutine nr3_r2g_trn_fic(orient_av, Ns, omge, nnge, ddge, Kdge, Kdee, &
                      gofts, ptn, SS1, Nmax, &
                      it2, t1s, t3s, rwa, rmin, resp)
    !
    ! R2g energy transfer response of first order of an three band
    ! multi-level system
    !
    ! Implementation with full interface
    !
    !
    !
    use acetodef 
    use acetolab
    use acetoaux
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
    !real(8), dimension(:,:) :: Ueet2
    integer :: Nmax  ! order of transfer to which we should integrate
        
    ! original arguments           
    integer, intent(in) :: it2
    real(8), dimension(:), intent(in) :: t1s,t3s
    real(8), intent(in) :: rwa
    real(8), intent(in) :: rmin
    complex(8), dimension(:,:), intent(inout) :: resp

    ! local
    integer :: Ng, Ne
    integer :: Nt1, Nt3, Ntsbi
    integer :: it1, it3
    integer :: e1, e2, g1, f1
    integer :: ino, tau, ec, ed
    real(8) :: t1,t2,t3,dt
    real(8), dimension(:,:), allocatable :: oafac
    complex(8) :: r, prod, exparg, pref, Knstt2
           
    real(8), dimension(:,:,:), allocatable :: ss2
        
    ! lineshape function values
    complex(8), dimension(:), allocatable :: gn_t1, gn_t3
    complex(8), dimension(:), allocatable :: gn_t2, gn_t1t2, gn_t2t3, gn_t1t2t3
    complex(8), dimension(:,:), allocatable :: gn_t2_td, gn_t2t3_td
    
    complex(8) :: gg_11_t1, gg_22_t3, gg_12_t2, gg_12_t2t3
    complex(8) :: gg_12_t1t2, gg_12_t1t2t3
    complex(8), dimension(:), allocatable :: gg_c1_t2_td, gg_d1_t2_td 
    complex(8), dimension(:), allocatable :: gg_c1_t2t3_td, gg_d1_t2t3_td
    
    complex(8), dimension(:,:,:), allocatable :: Gst
    complex(8), dimension(:,:), allocatable :: Knst, Knm1st, Kt2t3
    real(8), dimension(:), allocatable :: taus
 

    ! minimal value of the dipole factor
    real(8) :: minfac    
           
        
    ! nothing will be done if only zero's order is concidered
    if (Nmax <= 0) return        
        
    Nt1 = size(t1s)
    Nt3 = size(t3s)
    Ntsbi = size(gofts,2)
    
    Ng = Ns(1)
    Ne = Ns(2)
    
    allocate(oafac(Ne,Ne))
    allocate(ss2(Ne,Ne,Ne))
    allocate(gn_t1(Ne), gn_t2(Ne), gn_t3(Ne))
    allocate(gn_t1t2(Ne), gn_t2t3(Ne), gn_t1t2t3(Ne))
    allocate(gn_t2_td(Ne, it2),gn_t2t3_td(Ne, it2))
    allocate(Gst(Ne,Ne,it2),Knst(Ne,it2),Knm1st(Ne,it2),taus(it2))
    allocate(gg_c1_t2_td(it2), gg_c1_t2t3_td(it2))
    allocate(gg_d1_t2_td(it2), gg_d1_t2t3_td(it2))
    allocate(Kt2t3(Ne,Ne))
    
    call set_goft_mixing(SS1,ss2)

    ! initial and final states (normally there is a sum over them)
    g1 = 1
    f1 = 1 ! this points to ground state
    
    call set_dipole_factor_gt(g1, f1, "R2gt", orient_av, Ne, ddge, nnge, oafac, &
                           rmin, minfac)

    t2 = t1s(it2)
    dt = t1s(2)-t1s(1)
    
    print *, "T2 = ", t2, it2

    ! here we set lineshape functions over which we will integrate
    call set_goft_g_td(gn_t2_td, it2, it2, Ntsbi, gofts, ptn, t1s)

    taus = t1s(1:it2)    
    
    do it3 = 1,Nt3  
    
      call set_goft_g_td(gn_t2t3_td, it2, it2+it3, Ntsbi, gofts, ptn, t1s)
      
      t3 = t3s(it3)
      
      
      do e1 = 1, Ne
      do e2 = 1, Ne
      
              !
              ! integration: as many integrations as there are orders of
              ! transfer to be considered
              !

              Knst = 0
              Knst(e1,:) = 1.0d0  ! st (starting state is e1): this is Knst[a,tau] = 1[a,st]
              
              !
              ! Prepare lineshape function factors
              !
              Gst = 0.0d0

              
              if (Nmax == 1) then
                  
                  ! e1 is starting, e2 is final
                  gg_c1_t2_td(:) = matmul(ss2(:,e1,e2),gn_t2_td(:,:))
                  gg_c1_t2t3_td(:) = matmul(ss2(:,e1,e2),gn_t2t3_td(:,:))
                  gg_d1_t2_td(:) = matmul(ss2(:,e2,e2),gn_t2_td(:,:))
                  gg_d1_t2t3_td(:) = matmul(ss2(:,e2,e2),gn_t2t3_td(:,:))
                  
                  do tau = 1,it2
                      Gst(e2,e1,tau) = &
                      -(0.0d0,2.0d0)*imag(gg_c1_t2_td(tau)-gg_c1_t2t3_td(tau)) &
                      +(0.0d0,2.0d0)*imag(gg_d1_t2_td(tau)-gg_d1_t2t3_td(tau))
                  end do
                 

              else

                  do ec = 1, Ne
                  do ed = 1, Ne
                      gg_c1_t2_td(:) = matmul(ss2(:,ec,e2),gn_t2_td(:,:))
                      gg_c1_t2t3_td(:) = matmul(ss2(:,ec,e2),gn_t2t3_td(:,:))
                      gg_d1_t2_td(:) = matmul(ss2(:,ed,e2),gn_t2_td(:,:))
                      gg_d1_t2t3_td(:) = matmul(ss2(:,ed,e2),gn_t2t3_td(:,:))
                      
                      do tau = 1,it2
                          Gst(ec,ed,tau) = &
                          -(0.0d0,2.0d0)*imag(gg_c1_t2_td(tau)-gg_c1_t2t3_td(tau)) &
                          +(0.0d0,2.0d0)*imag(gg_d1_t2_td(tau)-gg_d1_t2t3_td(tau))
                      end do
                      
                  end do
                  end do
                  
              end if
              
              Knstt2 = 0.0d0
              
              !
              ! Integrate Nmax times
              !
              !
              ! We use diagonal rates with negative values 
              !
              do ino = 1, Nmax
              
                  Knm1st(:,:) = Knst(:,:)
                  
                  if (ino == Nmax) then  ! last integration
                      !do ec = 1, Ne
                        ec = e1
                        if (ec /= e2) then
                          !do tau = 1, it2 ! integration goes to t2
                            Knstt2 = Knstt2 &
                            + Kdee(e2,ec)*sum(exp(-(Kdee(e2,e2)-Kdee(ec,ec))*taus(1:it2)) &
                            *exp(Gst(e2,ec,1:it2))*Knm1st(ec,1:it2))*dt
                                           
                          !end do
                        end if
                      !end do
                  else           ! intermediate integrations
                      do ec = 1, Ne
                      do ed = 1, Ne
                        if (ec /= ed) then
                          do tau = 1, it2
                            Knst(ec,tau) = Knst(ec,tau) &
                            + Kdee(ec,ed)*sum(exp(-(Kdee(ec,ec)-Kdee(ed,ed))*taus(1:tau)) &
                            *exp(Gst(ec,ed,1:tau))*Knm1st(ed,1:tau))*dt
                          end do
                        end if
                      end do
                      end do
                 end if              
              
              end do

              Kt2t3(e2,e1) = Knstt2*exp(Kdee(e2,e2)*t2)
              
              !if (it2 == 1) then
              !    print *, e2, e1, Kt2t3(e2,e1)
              !    Kt2t3(e2,e1) = 0.0
              !end if
              !if (e2 == e1) then
              !    Kt2t3(e2,e1) = 1.0
              !else
              !    Kt2t3(e2,e1) = 0.0
              !end if
              
      end do
      end do



    do it1 = 1, Nt1
      
      r = 0.0d0
      t1 = t1s(it1)
      
      call set_goft_g(gn_t1, it1, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t2, it2, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t3, it3, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2, it1+it2, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t2t3, it2+it3, Ntsbi, gofts, ptn, t1s)
      call set_goft_g(gn_t1t2t3, it1+it2+it3, Ntsbi, gofts, ptn, t1s)  
          
      ! assuming Ng = 1
      do e1 = 1, Ne
      do e2 = 1, Ne
                 
                 
                           
          !if (oafac(e1,e2) > minfac) then

              ! frequencies
              exparg = -j1*((omge(1,e1)+rwa)*t1 - (omge(1,e2)+rwa)*t3) 

              ! dephasing
              exparg = exparg + & 
                (-Kdge(1,e1)*t1 - Kdge(1,e2)*t3)
              
              gg_11_t1 = dot_product(ss2(:,e1,e1),gn_t1)
              gg_22_t3 = dot_product(ss2(:,e2,e2),gn_t3)
              gg_12_t1t2 = dot_product(ss2(:,e1,e2),gn_t1t2)
              gg_12_t1t2t3 = dot_product(ss2(:,e1,e2),gn_t1t2t3)
              gg_12_t2 = dot_product(ss2(:,e1,e2),gn_t2)
              gg_12_t2t3 = dot_product(ss2(:,e1,e2),gn_t2t3)
                      
              ! Prefactor of the transfer response function
              pref = exp(conjg(- gg_11_t1 - gg_22_t3        &
                               - gg_12_t1t2 + gg_12_t1t2t3) &
                               + gg_12_t2 - gg_12_t2t3)   
                   
            

              Knstt2 = Kt2t3(e2,e1)
            
              ! final expression
              prod = oafac(e1,e2)*pref*exp(exparg)*Knstt2  !*Ueet2(e2,e1)
          
              r = r + prod
          
          !end if
    
      end do
      end do
    
      resp(it3,it1) = resp(it3,it1) + r
    
    end do
    end do

end subroutine nr3_r2g_trn_fic