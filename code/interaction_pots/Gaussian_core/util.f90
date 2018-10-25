! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!=============================================================================!
!                                 U T I L                                     !
!=============================================================================!
!                                                                             !
! $Id: util.f90,v 1.1.1.1.2.1 2007/11/06 16:12:57 ccseac Exp $
!                                                                             !
!-----------------------------------------------------------------------------!
! Contains ancilliary utility routines used throughout the code. Many of      !
! these will be inlined manually when used for efficiency.                    !
!-----------------------------------------------------------------------------!
!                                                                             !
! $Log$
!
!
!=============================================================================!

module util

  use iso_c_binding
  implicit none

  ! Double precision
  integer,parameter :: dp = c_double

  ! pi
  real(kind=dp),parameter :: Pi = 3.141592653589793238462643383279502884197_dp
  real(kind=dp),parameter :: invPi = 1.0_dp/3.141592653589793238462643383279502884197_dp
  
  contains

  real(kind=dp) function util_determinant(matrix)
    !------------------------------------------------------------------------------!
    ! Computes the determinant of a 3x3 matrix.                                    !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2006                                                     !
    !------------------------------------------------------------------------------!     

    implicit none
    real(kind=dp),dimension(3,3) :: matrix

    real(kind=dp) :: det

    Det = Matrix(1,1)*(Matrix(2,2)*Matrix(3,3) - &
          Matrix(2,3)*Matrix(3,2))
    Det = Det - Matrix(1,2)*(Matrix(2,1)*Matrix(3,3) - &
          Matrix(2,3)*Matrix(3,1))
    Det = Det + Matrix(1,3)*(Matrix(2,1)*Matrix(3,2) - &
          Matrix(2,2)*Matrix(3,1))

    util_determinant = det

    !if ( det < tiny(0.0_dp) ) stop 'Error in util_determinant'

    return

  end function util_determinant

  subroutine util_recipmatrix(hmatrix,recip_matrix)
    !------------------------------------------------------------------------------!
    ! Calculates the matrix of reciprocal lattive vectors from the h_matrix        !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2006                                                     !
    !------------------------------------------------------------------------------! 
    implicit none
    real(kind=dp),dimension(3,3),intent(in)  :: hmatrix
    real(kind=dp),dimension(3,3),intent(out) :: recip_matrix
    real(kind=dp) :: vol

    ! invert hmatrix to get recip_matrix
    recip_matrix(1,1)=hmatrix(2,2)*hmatrix(3,3)-hmatrix(2,3)*hmatrix(3,2)
    recip_matrix(1,2)=hmatrix(2,3)*hmatrix(3,1)-hmatrix(2,1)*hmatrix(3,3)
    recip_matrix(1,3)=hmatrix(2,1)*hmatrix(3,2)-hmatrix(2,2)*hmatrix(3,1)
    
    recip_matrix(2,1)=hmatrix(1,3)*hmatrix(3,2)-hmatrix(1,2)*hmatrix(3,3)
    recip_matrix(2,2)=hmatrix(1,1)*hmatrix(3,3)-hmatrix(1,3)*hmatrix(3,1)
    recip_matrix(2,3)=hmatrix(1,2)*hmatrix(3,1)-hmatrix(1,1)*hmatrix(3,2)
    
    recip_matrix(3,1)=hmatrix(1,2)*hmatrix(2,3)-hmatrix(1,3)*hmatrix(2,2)
    recip_matrix(3,2)=hmatrix(1,3)*hmatrix(2,1)-hmatrix(1,1)*hmatrix(2,3)
    recip_matrix(3,3)=hmatrix(1,1)*hmatrix(2,2)-hmatrix(1,2)*hmatrix(2,1)
    
    ! Calculte cell volume
    vol =hmatrix(1,1)*recip_matrix(1,1) + &
         hmatrix(1,2)*recip_matrix(1,2) + &
         hmatrix(1,3)*recip_matrix(1,3)

    ! Scale reciprocal lattice by 2*pi/volume
    recip_matrix(:,:)=recip_matrix(:,:)*2.0_dp*Pi/vol

    return

  end subroutine util_recipmatrix


  subroutine util_images(v_ss,hmatrix,recip_matrix)
    !------------------------------------------------------------------------------!
    ! Computes the minimum image of a vector v_ss - Note that this is inlined for  !
    ! efficiency in most cases.                                                    !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2006                                                     !
    !------------------------------------------------------------------------------! 
    implicit none
    real(kind=dp),dimension(3),intent(inout) :: v_ss
    real(kind=dp),dimension(3,3),intent(in)  :: hmatrix,recip_matrix

    real(kind=dp) :: sx,sy,sz
    real(kind=dp) :: ssx,ssy,ssz

    integer :: idim

    ! compute fractional co-ordinates
    sx = recip_matrix(1,1)*v_ss(1) + &
         recip_matrix(2,1)*v_ss(2) + &
         recip_matrix(3,1)*v_ss(3)
    sy = recip_matrix(1,2)*v_ss(1) + &
         recip_matrix(2,2)*v_ss(2) + &
         recip_matrix(3,2)*v_ss(3)  
    sz = recip_matrix(1,3)*v_ss(1) + &
         recip_matrix(2,3)*v_ss(2) + &
         recip_matrix(3,3)*v_ss(3) 


    sx = sx*0.5_dp*invPi 
    sy = sy*0.5_dp*invPi
    sz = sz*0.5_dp*invPi 

    ! apply boundary conditions
    ssx = sx - floor(sx+0.5_dp,kind=dp)
    ssy = sy - floor(sy+0.5_dp,kind=dp)
    ssz = sz - floor(sz+0.5_dp,kind=dp)
    

    ! scale back up
    do idim=1,3
       v_ss(idim) = hmatrix(idim,1)*ssx + &
                    hmatrix(idim,2)*ssy + &
                    hmatrix(idim,3)*ssz
    end do   


    return

  end subroutine util_images

end module util


