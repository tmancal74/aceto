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
        do e2 = 1, ne
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
    print *, minfac
    
end subroutine set_dipole_factor_g

subroutine set_dipole_factor_f(g1, f1, pathw, orient_av, Ne, ddge, nnge, &
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

    if (pathw == 'R1f') then
        do e1 = 1, Ne
        do e2 = 1, ne
           oafac(e1,e2) = LAB%get_oafactor(nnge(:,g1,e1),nnge(:,g1,e2), &
                                           nnge(:,f1,e2),nnge(:,f1,e1))
           oafac(e1,e2) = oafac(e1,e2)*ddge(g1,e1)*ddge(g1,e2)* &
                                       ddge(f1,e2)*ddge(f1,e1)       
        end do
        end do       
    else if (pathw == 'R2f') then
        do e1 = 1, Ne
        do e2 = 1, ne
           oafac(e1,e2) = LAB%get_oafactor(nnge(:,g1,e1),nnge(:,g1,e2), &
                                           nnge(:,f1,e1),nnge(:,f1,e2))
           oafac(e1,e2) = oafac(e1,e2)*ddge(g1,e1)*ddge(g1,e2)* &
                                       ddge(f1,e1)*ddge(f1,e2)       
        end do
        end do
    else if (pathw == 'R3f') then
        do e1 = 1, Ne
        do e2 = 1, ne
           oafac(e1,e2) = LAB%get_oafactor(nnge(:,g1,e1),nnge(:,g1,e1), &
                                           nnge(:,f1,e2),nnge(:,f1,e2))
           oafac(e1,e2) = oafac(e1,e2)*ddge(g1,e1)*ddge(g1,e1)* &
                                       ddge(f1,e2)*ddge(f1,e2)       
        end do
        end do                        
    else if (pathw == 'R4f') then
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
    
end subroutine set_dipole_factor_f

subroutine set_goft_g(gn, it, nmax, gofts, ptn, t1s)
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

end module acetoaux

