!
! This file contains routines calculating response functions of third
! order  
!
!

subroutine nr3_r2g_fi(orient_av, Ns, omge, nnge, ddge, Kdge, Kdee, &
                    t2, t1s, t3s, resp)
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
           
    ! original arguments           
    real(8), intent(in) :: t2
    real(8), dimension(:), intent(in) :: t1s,t3s
    complex(8), dimension(:,:), intent(inout) :: resp
           
    ! local
    integer :: Ng, Ne
    integer :: Nt1, Nt3
    integer :: it1, it3
    integer :: e1, e2, g1, f1
    real(8) :: t1,t3
    real(8) :: oafac
    complex(8) :: r, prod, exparg
    
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
    
    do it1 = 1,Nt1
    do it3 = 1,Nt3
    
      r = 0.0d0
      t1 = t1s(it1)
      t3 = t3s(it3)
     
      ! assuming Ng = 1
      do e1 = 1, Ne
        do e2 = 1, Ne
        
          oafac = LAB%get_oafactor(nnge(:,g1,e1),nnge(:,g1,e2), &
                                   nnge(:,f1,e1),nnge(:,f1,e2))
          
          ! transition dipole moments
          prod = &
            oafac*ddge(g1,e1)*ddge(g1,e2)*ddge(f1,e1)*ddge(f1,e2)
          
            
            
          ! frequencies
          exparg = -j1*(omge(e1,1)*t1 - omge(e2,1)*t3) &
                   +j1*(omge(e2,1)-omge(e1,1))*t2
    
          ! dephasing and decay
          exparg = exparg + & 
            (-Kdge(1,e1)*t1 - Kdge(1,e2)*t3) - Kdee(e2,e1)*t2      
    
          ! line-shape functions
          exparg = exparg + &
            0.0
            
          ! exp
          prod = prod*exp(exparg)
          
          r = r + prod
    
    
        end do
      end do
    
      resp(it1,it3) = resp(it1,it3) + r
    
    end do
    end do


end subroutine nr3_r2g_fi
