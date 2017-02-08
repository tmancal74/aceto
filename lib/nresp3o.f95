!
!
!
!
!
!
!

module nresp3o_mod

contains

subroutine R2g_2exnoR_dpc_none(i2, i3, i1, dts, Nen, om_eg, om_ee, &
                               dd_eg, gg_ee, rho0, RR)
  !
  ! At specific times given by indices i1, i2, and i3 calculates
  ! the value of the response function R1g
  !
  ! Naming of the variables follows the calling routine above
  !
  !use :: iso_c_binding
  implicit none
  
  integer, parameter :: c_double = kind(1.0d0)
  integer, parameter :: c_double_complex = kind((0.0d0,1.0d0))
  
  integer, intent(in) :: i2, i3, i1
  real(c_double), dimension(:), intent(in) :: dts
  integer, dimension(:), intent(in) :: Nen
  real(c_double), dimension(:,:), intent(in) :: om_eg, om_ee
  real(c_double), dimension(:,:), intent(in) :: dd_eg
  complex(c_double_complex), dimension(:,:,:), intent(in) :: gg_ee
  complex(c_double_complex), dimension(:,:), intent(in) :: rho0
  complex(c_double_complex), dimension(:,:,:) :: RR
  ! state indices
  integer :: ng1, ng2, ne1, ne2, Ntmax
  complex(c_double_complex) :: res
  real(c_double) :: t1, t2, t3, gam
  complex(c_double_complex), dimension(size(gg_ee,1),size(gg_ee,2)) :: &
          gg1,gg2,gg3,gg12,gg13,gg23,gg123

  ! values of time variables
  t1 = (i1-1)*dts(1)
  t2 = (i2-1)*dts(2)
  t3 = (i3-1)*dts(3)

  ! prepare line-shape functions for the given time
  Ntmax = size(gg_ee,3)
  ! long times will be interpolated
  if ((i1 + i2) > Ntmax) then
     gg12 = 0.0
  else
     gg12 = gg_ee(:,:,i1+i2)
  end if
  if ((i1 + i3) > Ntmax) then
     gg13 = 0.0
  else
     gg13 = gg_ee(i1+i3,:,:)
  end if
  if ((i2 + i3) > Ntmax) then
     gg23 = 0.0
  else
     gg23 = gg_ee(i2+i3,:,:)
  end if
  if ((i1 + i2 + i3) > Ntmax) then
     gg123 = 0.0
  else
     gg123 = gg_ee(i1+i2+i3,:,:)
  end if
  ! standard times
  gg1 = gg_ee(i1,:,:)
  gg2 = gg_ee(i2,:,:)
  gg3 = gg_ee(i3,:,:)
  
  ! ficticious dephasing for debugging
  gam = 1.0/100.0
  
  res = 0.0
  do ng1 = 1, Nen(1)
     do ne1 = 1, Nen(2)
        do ng2 = 1, Nen(1)
           do ne2 = 1, Nen(2)
           
!            Diagram R2g
!                                           
!            |ng2>   <ng2|
!       <----|-----------|
!            |ne1>   <ng2|U(ne1,ng1)*e^{-i*om(ne1,ng2)*t3}
!            |-----------|---->
!            |ne1>   <ne2|e^{-i*om(ne1,ne2)*t2}
!            |-----------|<----
!            |ne1>   <ng1| U(ne1,ng1)*e^{-i*om(ne1,ng1)*t1}
!       ---->|-----------|
!            |ng1>   <ng1|
!                                          

               ! add dipole section, frequency section, lineshape functions
               res = res + &
                   dd_eg(ne1,ng2)*dd_eg(ne2,ng2)* &
                   dd_eg(ne2,ng1)*dd_eg(ne1,ng1)*&
                   exp((0.0d0,-1.0d0)*(om_eg(ne1,ng1)*t1+om_ee(ne1,ne2)*t2 &
                                      +om_eg(ne1,ng2)*t3))* &
                   exp(-gam*t1 -gam*t3)* &
                   rho0(ng1,ng1)
              
           end do
        end do
     end do
  end do
  
  RR(i2,i3,i1) = RR(i2,i3,i1) + res

end subroutine R2g_2exnoR_dpc_none

end module nresp3o_mod





subroutine nresp3o_2exnoR_dpc_none(en, Nen, dd_eg, dd_fe, dephrates, relxrates, &
                                Nt, gg_ee, gg_ef, gg_ff, t1, t2, t3, rho0, RR, &
                                which_resp)
  use iso_c_binding

  use nresp3o_mod

  implicit none
  !
  ! Non-linear response of the 3rd order for an excitonic system with up to 2-exciton manifold
  ! Standard non-accellerated version
  !
  !
  ! Parameters
  ! ----------
  ! 
  ! en : double
  !    vector containg all energies
  !
  !
  ! number of ground state, one-exciton and double-exciton levels
  !
  !
  ! transition dipole moments (ground to one-exciton)
  !
  !
  ! transition dipole moments (one-exciton to two-exciton)
  ! 
  !
  ! vector containing dephasing rates
  !
  !
  ! vectors containing times on the t1, t2 and t3 axis
  !
  !
  ! initial density matrix of the ground state
  !
  !
  ! which responses to calculate. This is an integer array which contains indices of 
  ! responses that should be calculated. The mapping is the following
  !
  ! 1 = R1g
  ! 2 = R2g
  ! 3 = R3g
  ! 4 = R4g
  ! 5 = R1f*
  ! 6 = R2f
  ! and these not very common ones
  ! 7 = R3f*
  ! 8 = R4f
  !
  ! order of responses in the RR array corresponds to orders of indices in which_resp
  !
  !
  ! Respose function to be calculate
  !
  real(c_double), dimension(:), target, intent(in) :: en
  integer, dimension(3), intent(in) :: Nen
  real(c_double), dimension(:,:), intent(in) :: dd_eg
  real(c_double), dimension(:,:), intent(in) :: dd_fe
  real(c_double), dimension(:), intent(in) :: dephrates
  real(c_double), dimension(:), intent(in) :: relxrates
  integer, intent(in) :: Nt
  complex(c_double_complex), dimension(:,:,:), intent(in) :: gg_ee
  complex(c_double_complex), dimension(:,:,:), intent(in) :: gg_ef
  complex(c_double_complex), dimension(:,:,:), intent(in) :: gg_ff
  real(c_double), dimension(:), intent(in) :: t1
  real(c_double), dimension(:), intent(in) :: t2
  real(c_double), dimension(:), intent(in) :: t3
  complex(c_double_complex), dimension(:,:), intent(in) :: rho0
  integer, dimension(:), intent(in) :: which_resp
  complex(c_double_complex), dimension(:,:,:,:) :: RR


  !
  ! local variables
  !

  ! ground state energies
  real(c_double), dimension(:), pointer :: en_g
  ! one-exciton energies
  real(c_double), dimension(:), pointer :: en_e
  ! double-exciton energies
  real(c_double), dimension(:), pointer :: en_f
  ! lengths of the t1, t2 and t3 vectors
  integer :: Nt1, Nt2, Nt3

  real(c_double), dimension(Nen(2),Nen(1)) :: om_eg
  real(c_double), dimension(Nen(2),Nen(2)) :: om_ee
  real(c_double), dimension(3) :: dts

  ! time indices
  integer :: i1, i2, i3

  !
  integer :: k, j

  integer :: Nwr

  en_g => en(1:Nen(1))
  en_e => en(Nen(1)+1:Nen(1)+Nen(2))
  en_f => en(Nen(1)+Nen(2)+1:Nen(1)+Nen(2)+Nen(3))
  Nt1 = size(t1)
  Nt2 = size(t2)
  Nt3 = size(t3)

  Nwr = size(which_resp)
  
  do i1 = 1, Nen(2)
     do i2 = 1, Nen(1)
        om_eg(i1,i2) = en_e(i1)-en_g(i2)
     end do
  end do
  
  do i1 = 1, Nen(2)
     do i2 = 1, Nen(2)
        om_ee(i1,i2) = en_e(i1)-en_e(i2)
     end do
  end do
        

  RR = 0.0

  do i2 = 1, Nt2
     do i3 = 1, Nt3
        do i1 = 1, Nt1

           do k = 1, Nwr

              j = which_resp(k)
              
              !
              ! each call adds to the response funtion
              !
              select case(j)
                 
              case(1)
                 
                 !call R1g_up2ex_dpc_none(i2, i3, i1, Nen, &
                 !                        dd_eg, rho0, RR(k,:,:,:))
                 print *, "Not implemented"
              case(2)
                 
                 call R2g_2exnoR_dpc_none(i2, i3, i1, dts, Nen, om_eg, om_ee, &
                                         dd_eg, gg_ee, rho0, RR(k,:,:,:))
                 
              case(3)
                 
                 !call R1g_up2ex_dpc_none(i2, i3, i1, Nen, &
                 !                        dd_eg, rho0, RR(k,:,:,:))
                 print *, "Not implemented"
              case(4)
                 
                 !call R1g_up2ex_dpc_none(i2, i3, i1, Nen, &
                 !                        dd_eg, rho0, RR(k,:,:,:))
                 print *, "Not implemented"
              case(5)
                 
                 !call R1g_up2ex_dpc_none(i2, i3, i1, Nen, &
                 !                        dd_eg, rho0, RR(k,:,:,:))
                 print *, "Not implemented"
              case(6)
                 
                 !call R1g_up2ex_dpc_none(i2, i3, i1, Nen, &
                 !                        dd_eg, rho0, RR(k,:,:,:))
                 print *, "Not implemented"
              case(7)
                 
                 !call R1g_up2ex_dpc_none(i2, i3, i1, Nen, &
                 !                        dd_eg, rho0, RR(k,:,:,:))
                 print *, "Not implemented"
              case(8)
                 
                 !call R1g_up2ex_dpc_none(i2, i3, i1, Nen, &
                 !                        dd_eg, rho0, RR(k,:,:,:))
                 print *, "Not implemented"
              case default
                 
                 stop "Unknown response function type"
                 
              end select

           end do
               
        end do
     end do
  end do 


end subroutine nresp3o_2exnoR_dpc_none

