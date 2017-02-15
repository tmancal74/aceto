!
! This file contains routines calculating response functions of third
! order  
!
!

subroutine nr3_r2g_fi(orient_av, Ns, omge, nnge, ddge, Kdge, Kdee, &
                      gofts, ptn, SS1, &
                      it2, t1s, t3s, rwa, resp)
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
    complex(8), dimension(:,:), intent(inout) :: resp
           
    ! local
    integer :: Ng, Ne
    integer :: Nt1, Nt3
    integer :: it1, it3
    integer :: e1, e2, g1, f1, k
    real(8) :: t1,t2,t3, tt1, tt2
    real(8), dimension(:,:), allocatable :: oafac
    complex(8) :: r, prod, exparg, aa, bb
           
    real(8), dimension(:,:,:), allocatable :: ss2
        
    ! lineshape function values
    complex(8), dimension(:), allocatable :: gn_t1, gn_t2, gn_t3
    complex(8), dimension(:), allocatable :: gn_t1t2, gn_t2t3, gn_t1t2t3

    complex(8) :: gg_21_t3, gg_21_t1, gg_21_t2, gg_22_t1t2
    complex(8) :: gg_11_t2t3, gg_21_t1t2t3
           
    ! reconstruction of the LAB object
    type(lab_settings) :: LAB
    
    real(8) :: minfac    
        
    LAB%orient_aver = orient_av
    
    minfac = 0.01
    
    Nt1 = size(t1s)
    Nt3 = size(t3s)
    
    Ng = Ns(1)
    Ne = Ns(2)
    
    ! initial and final states (normally there is a sum over them)
    g1 = 1
    f1 = 1 ! this points to ground state

    allocate(oafac(Ne,Ne))
    allocate(ss2(Ne,Ne,Ne))
    allocate(gn_t1(Ne), gn_t2(Ne), gn_t3(Ne))
    allocate(gn_t1t2(Ne), gn_t2t3(Ne), gn_t1t2t3(Ne))
    
    do e1 = 1, Ne
    do e2 = 1, ne
       ss2(:,e1,e2) = (SS1(:,e1)**2)*(SS1(:,e2)**2)
       oafac(e1,e2) = LAB%get_oafactor(nnge(:,g1,e1),nnge(:,g1,e2), &
                                       nnge(:,f1,e1),nnge(:,f1,e2))
       !print *, oafac(e1,e2)
    end do
    end do

    t2 = t1s(it2)


! The order of loops may be switched for better performace with larger
! systems
!      do e1 = 1, Ne
!      do e2 = 1, Ne

    do it1 = 1,Nt1
    do it3 = 1,Nt3
    
      r = 0.0d0
      t1 = t1s(it1)
      t3 = t3s(it3)

      do k = 1, Ne
          gn_t1(k) = gofts(ptn(k,k),it1)
          gn_t2(k) = gofts(ptn(k,k),it2)
          gn_t3(k) = gofts(ptn(k,k),it3)
      end do      
      if ((it1+it2) > Nt1) then
         do k = 1, Ne
             ! linear extrapolation
             aa = gofts(ptn(k,k),Nt1)
             bb = gofts(ptn(k,k),Nt1-1)
             tt2 = t1s(Nt1)
             tt1 = t1s(Nt1-1)
             gn_t1t2(k) = aa + (aa - bb)*((it1+it2-Nt1))
         end do
      else
         do k = 1, Ne
            gn_t1t2(k) = gofts(ptn(k,k),it1+it2)
         end do
      end if
      if ((it2+it3) > Nt1) then
         do k = 1, Ne
             ! linear extrapolation
             aa = gofts(ptn(k,k),Nt1)
             bb = gofts(ptn(k,k),Nt1-1)
             tt2 = t1s(Nt1)
             tt1 = t1s(Nt1-1)
             gn_t2t3(k) = aa + (aa - bb)*((it2+it3-Nt1)) 
         end do
      else
         do k = 1, Ne
            gn_t2t3(k) = gofts(ptn(k,k),it2+it3)
         end do
      end if
      if ((it1+it2+it3) > Nt1) then
         do k = 1, Ne
             ! linear extrapolation
             aa = gofts(ptn(k,k),Nt1)
             bb = gofts(ptn(k,k),Nt1-1)
             tt2 = t1s(Nt1)
             tt1 = t1s(Nt1-1)
             gn_t1t2t3(k) = aa + (aa - bb)*((it1+it2+it3-Nt1))           
         end do
      else
         do k = 1, Ne
            gn_t1t2t3(k) = gofts(ptn(k,k),it1+it2+it3)
         end do
      end if

!      ! assuming Ng = 1
      do e1 = 1, Ne
      do e2 = 1, Ne
          
          if (oafac(e1,e2) > minfac) then
             
          gg_21_t1 = dot_product(ss2(:,e2,e1),gn_t1)
          gg_21_t2 = dot_product(ss2(:,e2,e1),gn_t2)
          gg_21_t3 = dot_product(ss2(:,e2,e1),gn_t3)
          gg_22_t1t2 = dot_product(ss2(:,e2,e2),gn_t1t2)
          gg_11_t2t3 = dot_product(ss2(:,e1,e1),gn_t2t3)
          gg_21_t1t2t3 = dot_product(ss2(:,e2,e1),gn_t1t2t3)
          
          
          ! transition dipole moments
          prod = &
            oafac(e1,e2)*ddge(g1,e1)*ddge(g1,e2)*ddge(f1,e1)*ddge(f1,e2)
            
          ! frequencies
          exparg = -j1*((omge(1,e1)+rwa)*t1 - (omge(1,e2)+rwa)*t3) &
                   -j1*(omge(1,e2)-omge(1,e1))*t2
    
          ! dephasing
          exparg = exparg + & 
            (-Kdge(1,e1)*t1 - Kdge(1,e2)*t3)
            
          ! decay
          if (e1 /= e2) then
              exparg = exparg + &
                (- Kdee(e2,e1)*t2)
          end if
    
          ! line-shape functions
          exparg = exparg + &
            (-conjg(gg_21_t1) -conjg(gg_21_t3) + gg_21_t2 - conjg(gg_22_t1t2) &
             -gg_11_t2t3 + conjg(gg_21_t1t2t3))
            
          ! exp
          prod = prod*exp(exparg)
          
          r = r + prod
          
          end if
    
      end do
      end do
    
      resp(it1,it3) = resp(it1,it3) + r
    
    end do
    end do

!      end do
!      end do

end subroutine nr3_r2g_fi



subroutine nr3_r3g_fi(orient_av, Ns, omge, nnge, ddge, Kdge, Kdee, &
                      gofts, ptn, SS1, &
                      it2, t1s, t3s, rwa, resp)
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
    complex(8), dimension(:,:), intent(inout) :: resp
           
    ! local
    integer :: Ng, Ne
    integer :: Nt1, Nt3
    integer :: it1, it3
    integer :: e1, e2, g1, f1, k
    real(8) :: t1,t2,t3, tt1, tt2
    real(8) :: oafac
    complex(8) :: r, prod, exparg, aa, bb
           
    real(8), dimension(:,:,:), allocatable :: ss2
        
    ! lineshape function values
    complex(8), dimension(:), allocatable :: gn_t1, gn_t2, gn_t3
    complex(8), dimension(:), allocatable :: gn_t1t2, gn_t2t3, gn_t1t2t3

    complex(8) :: gg_21_t3, gg_21_t1, gg_21_t2, gg_22_t1t2
    complex(8) :: gg_11_t2t3, gg_21_t1t2t3
           
    ! reconstruction of the LAB object
    type(lab_settings) :: LAB
    
    LAB%orient_aver = orient_av
    
    Nt1 = size(t1s)
    Nt3 = size(t3s)
    
    Ng = Ns(1)
    Ne = Ns(2)
    
    ! initial and final states (normally there is a sum over them)
    g1 = 1
    f1 = 1 ! this points to ground state

    allocate(ss2(Ne,Ne,Ne))
    allocate(gn_t1(Ne), gn_t2(Ne), gn_t3(Ne))
    allocate(gn_t1t2(Ne), gn_t2t3(Ne), gn_t1t2t3(Ne))
    
    do e1 = 1, Ne
    do e2 = 1, ne
       ss2(:,e1,e2) = (SS1(:,e1)**2)*(SS1(:,e2)**2)
    end do
    end do

    t2 = t1s(it2)


! The order of loops may be switched for better performace with larger
! systems
!      do e1 = 1, Ne
!      do e2 = 1, Ne

    do it1 = 1,Nt1
    do it3 = 1,Nt3
    
      r = 0.0d0
      t1 = t1s(it1)
      t3 = t3s(it3)

      do k = 1, Ne
          gn_t1(k) = gofts(ptn(k,k),it1)
          gn_t2(k) = gofts(ptn(k,k),it2)
          gn_t3(k) = gofts(ptn(k,k),it3)
      end do      
      if ((it1+it2) > Nt1) then
         do k = 1, Ne
             ! linear extrapolation
             aa = gofts(ptn(k,k),Nt1)
             bb = gofts(ptn(k,k),Nt1-1)
             tt2 = t1s(Nt1)
             tt1 = t1s(Nt1-1)
             gn_t1t2(k) = aa + (aa - bb)*((it1+it2-Nt1))
         end do
      else
         do k = 1, Ne
            gn_t1t2(k) = gofts(ptn(k,k),it1+it2)
         end do
      end if
      if ((it2+it3) > Nt1) then
         do k = 1, Ne
             ! linear extrapolation
             aa = gofts(ptn(k,k),Nt1)
             bb = gofts(ptn(k,k),Nt1-1)
             tt2 = t1s(Nt1)
             tt1 = t1s(Nt1-1)
             gn_t2t3(k) = aa + (aa - bb)*((it2+it3-Nt1)) 
         end do
      else
         do k = 1, Ne
            gn_t2t3(k) = gofts(ptn(k,k),it2+it3)
         end do
      end if
      if ((it1+it2+it3) > Nt1) then
         do k = 1, Ne
             ! linear extrapolation
             aa = gofts(ptn(k,k),Nt1)
             bb = gofts(ptn(k,k),Nt1-1)
             tt2 = t1s(Nt1)
             tt1 = t1s(Nt1-1)
             gn_t1t2t3(k) = aa + (aa - bb)*((it1+it2+it3-Nt1))           
         end do
      else
         do k = 1, Ne
            gn_t1t2t3(k) = gofts(ptn(k,k),it1+it2+it3)
         end do
      end if

!      ! assuming Ng = 1
      do e1 = 1, Ne
      do e2 = 1, Ne

          gg_21_t1 = dot_product(ss2(:,e2,e1),gn_t1)
          gg_21_t2 = dot_product(ss2(:,e2,e1),gn_t2)
          gg_21_t3 = dot_product(ss2(:,e2,e1),gn_t3)
          gg_22_t1t2 = dot_product(ss2(:,e2,e2),gn_t1t2)
          gg_11_t2t3 = dot_product(ss2(:,e1,e1),gn_t2t3)
          gg_21_t1t2t3 = dot_product(ss2(:,e2,e1),gn_t1t2t3)
          

          oafac = LAB%get_oafactor(nnge(:,g1,e1),nnge(:,g1,e1), &
                                   nnge(:,f1,e2),nnge(:,f1,e2))
          
          ! transition dipole moments
          prod = &
            oafac*ddge(g1,e1)*ddge(g1,e2)*ddge(f1,e1)*ddge(f1,e2)
            
          ! frequencies
          exparg = -j1*((omge(1,e1)+rwa)*t1 - (omge(1,e2)+rwa)*t3) !&
                   !-j1*(omge(1,e2)-omge(1,e1))*t2
    
          ! dephasing
          exparg = exparg + & 
            (-Kdge(1,e1)*t1 - Kdge(1,e2)*t3)
            
          ! decay
          !if (e1 /= e2) then
          !    exparg = exparg + &
          !      (- Kdee(e2,e1)*t2)
          !end if
    
          ! line-shape functions
          exparg = exparg + &
            (-conjg(gg_21_t1) -conjg(gg_21_t3) + gg_21_t2 - conjg(gg_22_t1t2) &
             -gg_11_t2t3 + conjg(gg_21_t1t2t3))
            
          ! exp
          prod = prod*exp(exparg)
          
          r = r + prod
    
      end do
      end do
    
      resp(it1,it3) = resp(it1,it3) + r
    
    end do
    end do

!      end do
!      end do




end subroutine nr3_r3g_fi



