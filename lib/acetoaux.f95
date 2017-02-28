module acetoaux

contains

subroutine set_dipole_factor_g(g1, f1, pathw, orient_av, Ne, ddge, nnge, &
                             oafac, rtol, minfac)
    use acetolab
    implicit none

    integer, intent(in) :: g1, f1
    character(3) :: pathw
    real(8), dimension(:), intent(in) :: orient_av
    integer, intent(in) :: Ne
    real(8), dimension(:,:), intent(in) :: ddge
    real(8), dimension(:,:,:), intent(in) :: nnge
    
    real(8), dimension(:,:), intent(out) :: oafac
    real(8), intent(in)   :: rtol
    real(8), intent(out)  :: minfac
    
    integer :: e1, e2
    
    ! reconstruction of the LAB object
    type(lab_settings) :: LAB

    LAB%orient_aver = orient_av

    if (pathw == 'R1g') then
        do e1 = 1, Ne
        do e2 = 1, Ne
           oafac(e1,e2) = LAB%get_oafactor(nnge(:,g1,e1),nnge(:,g1,e2), &
                                           nnge(:,f1,e2),nnge(:,f1,e1))
           oafac(e1,e2) = oafac(e1,e2)*ddge(g1,e1)*ddge(g1,e2)* &
                                       ddge(f1,e2)*ddge(f1,e1)       
        end do
        end do       
    else if (pathw == 'R2g') then
        do e1 = 1, Ne
        do e2 = 1, ne
           oafac(e1,e2) = LAB%get_oafactor(nnge(:,g1,e1),nnge(:,g1,e2), &
                                           nnge(:,f1,e1),nnge(:,f1,e2))
           oafac(e1,e2) = oafac(e1,e2)*ddge(g1,e1)*ddge(g1,e2)* &
                                       ddge(f1,e1)*ddge(f1,e2)       
        end do
        end do
    else if (pathw == 'R3g') then
        do e1 = 1, Ne
        do e2 = 1, ne
           oafac(e1,e2) = LAB%get_oafactor(nnge(:,g1,e1),nnge(:,g1,e1), &
                                           nnge(:,f1,e2),nnge(:,f1,e2))
           oafac(e1,e2) = oafac(e1,e2)*ddge(g1,e1)*ddge(g1,e1)* &
                                       ddge(f1,e2)*ddge(f1,e2)       
        end do
        end do                        
    else if (pathw == 'R4g') then
        do e1 = 1, Ne
        do e2 = 1, ne
           oafac(e1,e2) = LAB%get_oafactor(nnge(:,g1,e1),nnge(:,g1,e1), &
                                           nnge(:,f1,e2),nnge(:,f1,e2))
           oafac(e1,e2) = oafac(e1,e2)*ddge(g1,e1)*ddge(g1,e1)* &
                                       ddge(f1,e2)*ddge(f1,e2)       
        end do
        end do                        
    else
        stop "Unsupported pathway"
    end if

    minfac = rtol*maxval(oafac)
    !print *, minfac
    
end subroutine set_dipole_factor_g

subroutine set_dipole_factor_gt(g1, f1, pathw, orient_av, Ne, ddge, nnge, &
                             oafac, rtol, minfac)
    use acetolab
    implicit none

    integer, intent(in) :: g1, f1
    character(4) :: pathw
    real(8), dimension(:), intent(in) :: orient_av
    integer, intent(in) :: Ne
    real(8), dimension(:,:), intent(in) :: ddge
    real(8), dimension(:,:,:), intent(in) :: nnge
    
    real(8), dimension(:,:), intent(out) :: oafac
    real(8), intent(in)   :: rtol
    real(8), intent(out)  :: minfac
    
    integer :: e1, e2
    
    ! reconstruction of the LAB object
    type(lab_settings) :: LAB

    LAB%orient_aver = orient_av

!    if (pathw == 'R1gt') then
!        do e1 = 1, Ne
!        do e2 = 1, ne
!           oafac(e1,e2) = LAB%get_oafactor(nnge(:,g1,e1),nnge(:,g1,e2), &
!                                           nnge(:,f1,e2),nnge(:,f1,e1))
!           oafac(e1,e2) = oafac(e1,e2)*ddge(g1,e1)*ddge(g1,e2)* &
!                                       ddge(f1,e2)*ddge(f1,e1)       
!        end do
!        end do       
!    else 
    if (pathw == 'R2gt') then
        do e1 = 1, Ne
        do e2 = 1, ne
           oafac(e1,e2) = LAB%get_oafactor(nnge(:,g1,e1),nnge(:,g1,e2), &
                                           nnge(:,f1,e1),nnge(:,f1,e2))
           oafac(e1,e2) = oafac(e1,e2)*ddge(g1,e1)*ddge(g1,e2)* &
                                       ddge(f1,e1)*ddge(f1,e2)       
        end do
        end do
    else if (pathw == 'R3gt') then
        do e1 = 1, Ne
        do e2 = 1, ne
           oafac(e1,e2) = LAB%get_oafactor(nnge(:,g1,e1),nnge(:,g1,e1), &
                                           nnge(:,f1,e2),nnge(:,f1,e2))
           oafac(e1,e2) = oafac(e1,e2)*ddge(g1,e1)*ddge(g1,e1)* &
                                       ddge(f1,e2)*ddge(f1,e2)       
        end do
        end do                        
!    else if (pathw == 'R4gt') then
!        do e1 = 1, Ne
!        do e2 = 1, ne
!           oafac(e1,e2) = LAB%get_oafactor(nnge(:,g1,e1),nnge(:,g1,e1), &
!                                           nnge(:,f1,e2),nnge(:,f1,e2))
!           oafac(e1,e2) = oafac(e1,e2)*ddge(g1,e1)*ddge(g1,e1)* &
!                                       ddge(f1,e2)*ddge(f1,e2)       
!        end do
!        end do                        
    else
        stop "Unsupported pathway"
    end if

    minfac = rtol*maxval(oafac)
    print *, minfac
    
end subroutine set_dipole_factor_gt

subroutine set_dipole_factor_f(g1, f1, pathw, orient_av, Ne, ddge, nnge, &
                               ddef, nnef, &
                               oafac, rtol, minfac)
    use acetolab
    implicit none

    integer, intent(in) :: g1, f1
    character(3) :: pathw
    real(8), dimension(:), intent(in) :: orient_av
    integer, intent(in) :: Ne
    real(8), dimension(:,:), intent(in) :: ddge
    real(8), dimension(:,:,:), intent(in) :: nnge
    real(8), dimension(:,:), intent(in) :: ddef
    real(8), dimension(:,:,:), intent(in) :: nnef    
    real(8), dimension(:,:), intent(out) :: oafac
    real(8), intent(in)   :: rtol
    real(8), intent(out)  :: minfac
    
    integer :: e1, e2, ff
    real(8) :: aux
    
    ! reconstruction of the LAB object
    type(lab_settings) :: LAB

    LAB%orient_aver = orient_av
    
    ff = f1
    !print *, "!!!!!!!!!! ff = ", ff
    if (pathw == 'R1f') then
        do e1 = 1, Ne
        do e2 = 1, Ne
           aux = LAB%get_oafactor(nnge(:,g1,e1),nnge(:,g1,e2), &
                                           nnef(:,e2,ff),nnef(:,e1,ff))
           oafac(e1,e2) = aux*ddge(g1,e1)*ddge(g1,e2)* &
                                       ddef(e2,ff)*ddef(e1,ff)  
        end do
        end do       
    else if (pathw == 'R2f') then
        do e1 = 1, Ne
        do e2 = 1, Ne
           aux = LAB%get_oafactor(nnge(:,g1,e1),nnge(:,g1,e2), &
                                           nnef(:,e1,ff),nnef(:,e2,ff))
           oafac(e1,e2) = aux*ddge(g1,e1)*ddge(g1,e2)* &
                                       ddef(e1,ff)*ddef(e2,ff)       
        end do
        end do
!    else if (pathw == 'R3f') then
!        do e1 = 1, Ne
!        do e2 = 1, ne
!           oafac(e1,e2) = LAB%get_oafactor(nnge(:,g1,e1),nnge(:,g1,e1), &
!                                           nnge(:,f1,e2),nnge(:,f1,e2))
!           oafac(e1,e2) = oafac(e1,e2)*ddge(g1,e1)*ddge(g1,e1)* &
!                                       ddge(f1,e2)*ddge(f1,e2)       
!        end do
!        end do                        
!    else if (pathw == 'R4f') then
!        do e1 = 1, Ne
!        do e2 = 1, ne
!           oafac(e1,e2) = LAB%get_oafactor(nnge(:,g1,e1),nnge(:,g1,e1), &
!                                           nnge(:,f1,e2),nnge(:,f1,e2))
!           oafac(e1,e2) = oafac(e1,e2)*ddge(g1,e1)*ddge(g1,e1)* &
!                                       ddge(f1,e2)*ddge(f1,e2)       
!        end do
!        end do                        
    else
        stop "Unsupported pathway"
    end if

    minfac = rtol*maxval(oafac)
    !print *, minfac, "(", rtol, maxval(oafac),")"
    
end subroutine set_dipole_factor_f

subroutine set_dipole_factor_ft(g1, f1, pathw, orient_av, Ne, ddge, nnge, &
                             oafac, rtol, minfac)
    use acetolab
    implicit none

    integer, intent(in) :: g1, f1
    character(4) :: pathw
    real(8), dimension(:), intent(in) :: orient_av
    integer, intent(in) :: Ne
    real(8), dimension(:,:), intent(in) :: ddge
    real(8), dimension(:,:,:), intent(in) :: nnge
    
    real(8), dimension(:,:), intent(out) :: oafac
    real(8), intent(in)   :: rtol
    real(8), intent(out)  :: minfac
    
    integer :: e1, e2
    
    ! reconstruction of the LAB object
    type(lab_settings) :: LAB

    LAB%orient_aver = orient_av

    if (pathw == 'R1ft') then
        do e1 = 1, Ne
        do e2 = 1, ne
           oafac(e1,e2) = LAB%get_oafactor(nnge(:,g1,e1),nnge(:,g1,e2), &
                                           nnge(:,f1,e2),nnge(:,f1,e1))
           oafac(e1,e2) = oafac(e1,e2)*ddge(g1,e1)*ddge(g1,e2)* &
                                       ddge(f1,e2)*ddge(f1,e1)       
        end do
        end do       
    else if (pathw == 'R2ft') then
        do e1 = 1, Ne
        do e2 = 1, ne
           oafac(e1,e2) = LAB%get_oafactor(nnge(:,g1,e1),nnge(:,g1,e2), &
                                           nnge(:,f1,e1),nnge(:,f1,e2))
           oafac(e1,e2) = oafac(e1,e2)*ddge(g1,e1)*ddge(g1,e2)* &
                                       ddge(f1,e1)*ddge(f1,e2)       
        end do
        end do
    else if (pathw == 'R3ft') then
        do e1 = 1, Ne
        do e2 = 1, ne
           oafac(e1,e2) = LAB%get_oafactor(nnge(:,g1,e1),nnge(:,g1,e1), &
                                           nnge(:,f1,e2),nnge(:,f1,e2))
           oafac(e1,e2) = oafac(e1,e2)*ddge(g1,e1)*ddge(g1,e1)* &
                                       ddge(f1,e2)*ddge(f1,e2)       
        end do
        end do                        
    else if (pathw == 'R4ft') then
        do e1 = 1, Ne
        do e2 = 1, ne
           oafac(e1,e2) = LAB%get_oafactor(nnge(:,g1,e1),nnge(:,g1,e1), &
                                           nnge(:,f1,e2),nnge(:,f1,e2))
           oafac(e1,e2) = oafac(e1,e2)*ddge(g1,e1)*ddge(g1,e1)* &
                                       ddge(f1,e2)*ddge(f1,e2)       
        end do
        end do                        
    else
        stop "Unsupported pathway"
    end if

    minfac = rtol*maxval(oafac)
    print *, minfac
    
end subroutine set_dipole_factor_ft


subroutine set_goft_g(gn, it, nmax, gofts, ptn, t1s)
    ! Assignment of line shapes in single exciton block
    !
    !
    implicit none
    complex(8), dimension(:), intent(out) :: gn
    integer :: nmax
    integer :: it
    complex(8), dimension(:,:), intent(in) :: gofts
    integer, dimension(:,:), intent(in) :: ptn
    real(8), dimension(:), intent(in) :: t1s
    ! local
    integer :: k, Ne
    real(8) :: tt1, tt2
    complex(8) :: aa, bb
           
    Ne = size(gn)
    if (it > nmax) then
        do k = 1, Ne
            ! linear extrapolation
            aa = gofts(ptn(k,k),nmax)
            bb = gofts(ptn(k,k),nmax-1)
            tt2 = t1s(nmax)
            tt1 = t1s(nmax-1)
            gn(k) = aa + (aa - bb)*((it-nmax))
        end do
    else
        do k = 1, Ne
            gn(k) = gofts(ptn(k,k),it)
        end do
    end if
                      
    
end subroutine set_goft_g

subroutine set_goft_f(gn, it, nmax, gofts, ptn, t1s)
    ! Assignment of line shapes in the two-exciton block
    !
    !
    !
    implicit none
    complex(8), dimension(:), intent(out) :: gn
    integer :: nmax
    integer :: it
    complex(8), dimension(:,:), intent(in) :: gofts
    integer, dimension(:,:), intent(in) :: ptn
    real(8), dimension(:), intent(in) :: t1s
    ! local
    integer :: k, l, Ne, a
    real(8) :: tt1, tt2
    complex(8) :: aa, bb, cc, dd
          
    ! pointer to sites has a size corresponding to number of sites
    Ne = size(ptn,1) 
    
    if (it > nmax) then
        do k = 1, Ne
        do l = k+1, ne
            ! linear extrapolation
            aa = gofts(ptn(k,k),nmax)
            bb = gofts(ptn(k,k),nmax-1)
            cc = gofts(ptn(l,l),nmax)
            dd = gofts(ptn(l,l),nmax-1)
            tt2 = t1s(nmax)
            tt1 = t1s(nmax-1)
            gn(k) = aa + (aa - bb)*((it-nmax))  &
                  + cc + (cc - dd)*((it-nmax))
        end do
        end do
    else
        a = 1
        do k = 1, Ne
        do l = k+1, ne
            gn(a) = gofts(ptn(k,k),it) + gofts(ptn(l,l),it)
            a = a + 1
        end do
        end do
    end if
                      
    
end subroutine set_goft_f

subroutine set_goft_mixing(SS1, ss2)
    real(8), dimension(:,:), intent(in) :: SS1
    real(8), dimension(:,:,:), intent(out) :: ss2 
    ! local
    integer :: Ne, e1, e2
    
    Ne = size(SS1,2)
    do e1 = 1, Ne
    do e2 = 1, ne
       ss2(:,e1,e2) = (SS1(:,e1)**2)*(SS1(:,e2)**2)
    end do
    end do

end subroutine set_goft_mixing

subroutine set_goft_mixing_22(SS2, A22, N1)
    ! Mixing matrix for calculation of two-exciton lineshape functions
    !
    !
    !
    !
    implicit none
    real(8), dimension(:,:), intent(in) :: SS2
    real(8), dimension(:,:,:), intent(out) :: A22
    integer, intent(in) :: N1
    ! local
    integer :: a, b, c, d
    integer :: e1, e2, f1, f2
    integer :: n, m, k, l
        
    A22 = 0.0d0
    
    ! exciton 1
    e1 = 1
    do a = 1, N1
    do b = a + 1, N1
      ! exciton 2
      e2 = 1
      do c = 1, N1
      do d = c + 1, N1
        ! local state 1
        f1 = 1
        do n = 1, N1
        do m = n + 1, N1
          ! local state 2
          f2 = 1
          do k = 1, N1
          do l = k + 1, N1
            
            if (n == k) then
              A22(n,e1,e2) = A22(n,e1,e2) +  &
                             SS2(e1,f1)*SS2(e1,f1)*SS2(e2,f2)*SS2(e2,f2)
            end if

            if (n == l) then
              A22(n,e1,e2) = A22(n,e1,e2) +  &
                             SS2(e1,f1)*SS2(e1,f1)*SS2(e2,f2)*SS2(e2,f2)
            end if
            
            if (m == k) then
              A22(m,e1,e2) = A22(m,e1,e2) +  &
                             SS2(e1,f1)*SS2(e1,f1)*SS2(e2,f2)*SS2(e2,f2)
            end if

            if (m == l) then
              A22(m,e1,e2) = A22(m,e1,e2) +  &
                             SS2(e1,f1)*SS2(e1,f1)*SS2(e2,f2)*SS2(e2,f2)
            end if
            
            f2 = f2 + 1            
          end do
          end do
          
          f1 = f1 + 1
        end do
        end do

        e2 = e2 + 1
      end do
      end do

      e1 = e1 + 1
    end do
    end do

end subroutine set_goft_mixing_22

subroutine set_goft_mixing_21(SS2, SS1, A21)
    ! FIXME: finish this
    implicit none
    real(8), dimension(:,:), intent(in) :: SS2
    real(8), dimension(:,:), intent(in) :: SS1
    real(8), dimension(:,:,:), intent(out) :: A21
    ! local
    integer :: N1, N2
    integer :: a, b, c, d
    integer :: e1, e2, f1
    integer :: n, m, k
    
    N1 = size(SS1,1)
    N2 = size(SS2,1)
    
    A21 = 0.0d0 

    e1 = 1
    do a = 1, N1
    do b = a+1, N1

      do e2 = 1, N1
      
         do k = 1, N1

           f1 = 1
           do n = 1, N1
           do m = n + 1, N1
             
             if (n == k) then
               A21(n,e1,e2) = A21(n,e1,e2) + &
                              SS2(e1,f1)*SS2(e1,f1)*SS1(e2,k)*SS1(e2,k)
             end if

             if (m == k) then
               A21(m,e1,e2) = A21(m,e1,e2) + &
                              SS2(e1,f1)*SS2(e1,f1)*SS1(e2,k)*SS1(e2,k)
             end if

           
             f1 = f1 + 1
           end do
           end do
    
         end do
         
      end do
      
      e1 = e1 + 1
    end do
    end do
    
end subroutine set_goft_mixing_21


end module acetoaux

