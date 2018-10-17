! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!=============================================================================!
!                               E N E R G Y                                   !
!=============================================================================!
! Contains C-compatible Fortran 2008 routines which calculate the energy of a !
! system of particles interacting via the standard Lennard-Jones 12-6 model.  !
!                                                                             !
! Key variables                                                               !
!    pos(3,natoms)    : Array of 3d particle positions (absolute coordinates) !
!    hmatrix(3,3)     : 3 vectors which define the supercell perioic lattice  !
!                                                                             !
! These will correspond to numpy arrays when invoked via a Python interface.  !
!=============================================================================!
module energy

  use iso_c_binding
  implicit None                                 ! Impose strong typing

  integer,parameter :: dp = c_double            ! C compatible double precision

  Private                                       ! Everything is private
  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
                                                     !... unless exposed here.
  public :: create_lattice
  public :: compute_lattice_energy
  public :: compute_local_energy
  public :: compute_neighbour_list

  !---------------------------------------------------------------------------!
  !                        P u b l i c   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!
  
  !---------------------------------------------------------------------------!
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!

  ! Number of lattices
  integer,parameter :: num_lattices=2
  
  ! Constants defining the Lennard-Jones potential
  real(kind=dp),save :: lj_epsilon = 1.0_dp
  real(kind=dp),save :: lj_sigma   = 1.0_dp
  real(kind=dp),save :: lj_cut     = 3.5_dp
  real(kind=dp),save :: lj_lrc     = 0.0_dp
  real(kind=dp),save :: skin       = 0.5_dp 
  
  type lattice_t

     ! Set when number of atoms has been imposed
     logical :: initialised = .false.              

     ! Number of atoms in the lattices  (for consistency checking)
     integer :: natoms = -1
  
     ! Current energy due to the model Hamiltonian
     real(kind=dp) :: energy = 0.0_dp

     ! Lattice translation vectors used to compute interactions with images
     integer :: nivect
     real(kind=dp),allocatable,dimension(:,:) :: ivect

     ! Neighbour list. Maxneigh controls max no. neighbours per atom
     integer  :: maxneigh = 50
     integer,allocatable,dimension(:) :: nn
     integer,allocatable,dimension(:,:) ::jn,vn

  end type lattice_t

  ! Array of these lattice variables
  type(lattice_t),dimension(0:num_lattices-1) :: lattice
  
  ! Cache to avoid repeated computation
  real(kind=dp),save  :: sigma12,sigma6

  ! Which imaging method to use (image vectors vs minimum image)
  integer :: image_flag
  integer,parameter :: flag_iv = 1 , flag_mi = 2
  character(13),dimension(2) :: image_method = (/"Image vectors","Minimum image"/)
  

contains

  subroutine create_lattice(n,d,pos,dh2,dh1,hmatrix,ilat) bind(c)
    !==============================================================================!
    ! Initialises all arrays and shared variables used in computing the energy of  !
    ! the system (or of a single particle). Called on first invocation of either   !
    ! compute_model_energy or compute_local_real_energy.                           !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2017                                                     !
    !------------------------------------------------------------------------------!
    use util,       only : util_determinant,util_images,util_recipmatrix,pi
    implicit none

    ! input variables
    integer,intent(in) :: n,d,dh2,dh1,ilat
    real(kind=dp),dimension(d,n),intent(in) :: pos
    real(kind=dp),dimension(dh1,dh2),intent(in) :: hmatrix
    
    integer :: im,jm,km    ! loop counters
    integer :: ierr ! error flag
    real(kind=dp) :: s3,density
   
    if (ilat>=num_lattices) stop 'Error in create_lattice : ilat out of range'
    
    ! Set number of atoms for sanity checks on subsequent calls
    lattice(ilat)%natoms = n

    
    ! Cache to avoid repeated computation
    sigma6  = lj_sigma**6
    sigma12 = sigma6**2
    
    ! Invariant part of long range correction
    s3     = lj_sigma**3   
    lj_lrc = (8.0_dp/3.0_dp)*pi*lj_epsilon*s3*( (1.0_dp/3.0_dp)*(s3**3/lj_cut**9) - (s3/lj_cut**3) )


    !===============!
    ! Image vectors !
    !===============!

    ! Check which imaging method to use
    image_flag = flag_mi ! minimum image

    if ( 2.0*(lj_cut + skin) >= sqrt(dot_product(hmatrix(:,1),hmatrix(:,1)))) image_flag = flag_iv
    if ( 2.0*(lj_cut + skin) >= sqrt(dot_product(hmatrix(:,2),hmatrix(:,2)))) image_flag = flag_iv
    if ( 2.0*(lj_cut + skin) >= sqrt(dot_product(hmatrix(:,3),hmatrix(:,3)))) image_flag = flag_iv  

    write(0,'("Using image method : ",A13)')image_method(image_flag)
    if ( image_flag == flag_iv ) then
       write(0,'("Warning - system too small for minimum image convention.")')
       write(0,'("Using image method only suitable for solids.")')
    end if

    if ( image_flag == flag_iv ) then
    
       ! Find out how many lattice vectors we need to consider all images within cutoff
       im = floor((lj_cut+skin)/sqrt(dot_product(hmatrix(:,1),hmatrix(:,1))))+1
       jm = floor((lj_cut+skin)/sqrt(dot_product(hmatrix(:,2),hmatrix(:,2))))+1
       km = floor((lj_cut+skin)/sqrt(dot_product(hmatrix(:,3),hmatrix(:,3))))+1
       
       lattice(ilat)%nivect = (2*im+1)*(2*jm+1)*(2*km+1)
       !write(0,'("DEBUG - number of ivects ",I4)')lattice(ilat)%nivect
       
       if (allocated(lattice(ilat)%ivect)) then
           deallocate(lattice(ilat)%ivect, stat=ierr)
           if (ierr/=0) stop 'Error in create_lattice : Could not deallocate ivect'
       end if
       allocate(lattice(ilat)%ivect(1:3,1:lattice(ilat)%nivect),stat=ierr)
       if (ierr/=0) stop 'Error in create_lattice: Could not allocate ivect'
       
       ! compute current image translation vectors
       call compute_ivects(hmatrix,ilat)

    end if
    
    !============!
    ! Neighbours !
    !============!

    ! Compute maxneigh based on density and a safety factor of 20%
    density = lattice(ilat)%natoms/abs(util_determinant(hmatrix))
    lattice(ilat)%maxneigh = ceiling((density*4.0*pi*(lj_cut+skin)**3)/3.0_dp * 1.2_dp)

    if (allocated(lattice(ilat)%nn)) then
        deallocate(lattice(ilat)%nn,stat=ierr)
        if (ierr/=0) stop 'Error in create_lattice : Could not deallocate nn'
    end if
    allocate(lattice(ilat)%nn(1:lattice(ilat)%natoms),stat=ierr)
    if (ierr/=0) stop 'Error in create_lattice: Could not allocate nn array'

    if (allocated(lattice(ilat)%jn)) then
        deallocate(lattice(ilat)%jn,stat=ierr)
        if (ierr/=0) stop 'Error in create_lattice : Could not deallocate jn'
    end if
    allocate(lattice(ilat)%jn(1:lattice(ilat)%maxneigh,1:lattice(ilat)%natoms),stat=ierr)
    if (ierr/=0) stop 'Error in energy_init: Could not allocat jn array'
    
    if ( image_flag == flag_iv ) then
       if (allocated(lattice(ilat)%vn)) then
         deallocate(lattice(ilat)%vn,stat=ierr)
         if (ierr/=0) stop 'Error in create_lattice : Could not deallocate vn'
       end if
       allocate(lattice(ilat)%vn(1:lattice(ilat)%maxneigh,1:lattice(ilat)%natoms),stat=ierr)
       if (ierr/=0) stop 'Error in create_lattice : Could not allocate vn array'
    end if
       
    ! Set initialised flag
    lattice(ilat)%initialised = .true.

    ! Make sure first neighbour list is calculated 
    call compute_neighbour_list(lattice(ilat)%natoms,3,pos,3,3,hmatrix,ilat)

    ! Make sure energy if up-to-date
    lattice(ilat)%energy = compute_lattice_energy(lattice(ilat)%natoms,3,pos,3,3,hmatrix,ilat)
    
    return

  end subroutine create_lattice

  
  subroutine compute_ivects(hmatrix,ilat)
    !------------------------------------------------------------------------------!
    ! Computes the translation vectors needed to include all images such that each !
    ! atom sees all the images (including those of itself) within the cut off.     !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2017                                                     !
    !------------------------------------------------------------------------------!
    implicit none
    integer,intent(in) :: ilat
    real(kind=dp),dimension(3,3),intent(in) :: hmatrix
 
    !loop counters
    integer :: icell,jcell,kcell,im,jm,km,k
    real(kind=dp),dimension(3) :: sx,sy,sz

    im = floor((lj_cut+skin)/sqrt(dot_product(hmatrix(:,1),hmatrix(:,1))))+1
    jm = floor((lj_cut+skin)/sqrt(dot_product(hmatrix(:,2),hmatrix(:,2))))+1
    km = floor((lj_cut+skin)/sqrt(dot_product(hmatrix(:,3),hmatrix(:,3))))+1
    
    lattice(ilat)%nivect = (2*im+1)*(2*jm+1)*(2*km+1)

    ! we'd like the central cell to be entry 0
    ! we can flag it as non-self interacting
    lattice(ilat)%ivect(:,1) = 0.0_dp
    
    k = 2
    do icell = -im,im
       sx = real(icell,kind=dp)*hmatrix(:,1)
       do jcell = -jm,jm
          sy = real(jcell,kind=dp)*hmatrix(:,2) 
          do kcell = -km,km
             sz = real(kcell,kind=dp)*hmatrix(:,3)
             
             if ( abs(icell)+abs(jcell)+abs(kcell) == 0 ) cycle
             lattice(ilat)%ivect(:,k)  = sx + sy + sz
             k = k + 1

          end do
       end do
    end do

    return

  end subroutine compute_ivects

  real(kind=dp) function compute_local_energy(dmol,n,d,pos,dh2,dh1,hmatrix,ilat) bind(c)
    !------------------------------------------------------------------------------!
    ! Call either the minimum image or image vector routine to compute the energy  !
    ! of a single particle interacting with its neighbours.                        !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2017                                                     !
    !------------------------------------------------------------------------------!    
    implicit none
    integer,value,intent(in) :: dmol,d,n,dh1,dh2,ilat
    real(kind=dp),intent(in),dimension(d,n) :: pos
    real(kind=dp),intent(in),dimension(dh1,dh2) :: hmatrix

    integer :: imol
    real(kind=dp) :: loc_energy

    ! Sanity checks
    if (d/=3) stop 'Error in compute_local_real_energy: position array must be 3xn'
    if ((dh1/=3).or.(dh2/=3)) stop 'Error in compute_local_real_energy: cell array must be 3x3'

    if (ilat>=num_lattices) stop 'Error in compute_local_energy : ilat out of range'

    if (lattice(ilat)%initialised) then
       if (n/=lattice(ilat)%natoms) stop 'Error in compute_local_energy: Number of entries in position array has changed!'
    else
       call create_lattice(n,d,pos,dh2,dh1,hmatrix,ilat)
    end if

    ! We assume this routine will only ever be called from C/Python (never Fortran)
    ! and so we correct imol into 1-based indexing
    imol = dmol + 1
    if ( (imol<1).or.(imol>lattice(ilat)%natoms) ) stop 'Error in compute_local_real_energy : atom index out of range!'
    
    select case (image_flag)
    case (flag_iv)
       loc_energy = compute_local_energy_iv(imol,n,d,pos,dh2,dh1,hmatrix,ilat)
    case (flag_mi)
       loc_energy = compute_local_energy_mi(imol,n,d,pos,dh2,dh1,hmatrix,ilat)
    case default
       stop 'Error in compute_local_real_energy - unknown image flag'       
    end select

    compute_local_energy = loc_energy

  end function compute_local_energy

  real(kind=dp) function compute_lattice_energy(n,d,pos,dh2,dh1,hmatrix,ilat) bind(c)
    !------------------------------------------------------------------------------!
    ! Call either the minimum image or image vector routine to compute the energy  !
    ! of the entire system.                                                        !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2017                                                     !
    !------------------------------------------------------------------------------!    
    implicit none
    integer,value,intent(in) :: d,n,dh1,dh2,ilat
    real(kind=dp),intent(in),dimension(d,n) :: pos
    real(kind=dp),intent(in),dimension(dh1,dh2) :: hmatrix

    integer :: imol
    real(kind=dp) :: loc_energy

    ! Sanity checks
    if (d/=3) stop 'Error in compute_local_real_energy: position array must be 3xn'
    if ((dh1/=3).or.(dh2/=3)) stop 'Error in compute_local_real_energy: cell array must be 3x3'

    if (ilat>=num_lattices) stop 'Error in create_lattice : ilat out of range'

    if (.not.lattice(ilat)%initialised) then
        call create_lattice(n,d,pos,dh2,dh1,hmatrix,ilat)
    else
        if (n/=lattice(ilat)%natoms) stop 'Error in compute_lattice_energy : Number of entries in position array has changed!'
    end if

    select case (image_flag)
    case (flag_iv)
       loc_energy = compute_lattice_energy_iv(n,d,pos,dh2,dh1,hmatrix,ilat)
    case (flag_mi)
       loc_energy = compute_lattice_energy_mi(n,d,pos,dh2,dh1,hmatrix,ilat)
    case default
       stop 'Error in compute_lattice_energy - unknown image flag'       
    end select

    lattice(ilat)%energy   = loc_energy
    compute_lattice_energy = loc_energy

  end function compute_lattice_energy
  
    
  real(kind=dp) function compute_local_energy_mi(imol,n,d,pos,dh2,dh1,hmatrix,ilat) bind(c)
    !------------------------------------------------------------------------------!
    ! Calculates the real-space contribution to the energy due particle dmol.      !
    ! To be used when computing the changes in energy due to a trial move.         !
    ! This version uses the minimum image convention
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2017                                                     !
    !------------------------------------------------------------------------------!    
    implicit none
    integer,value,intent(in) :: imol,d,n,dh1,dh2,ilat
    real(kind=dp),intent(in),dimension(d,n) :: pos
    real(kind=dp),intent(in),dimension(dh1,dh2) :: hmatrix
    
    ! local variables
    real(kind=dp) :: La,Lb,Lc,invLa,invLb,invLc
    real(kind=dp)                :: Evdw,rcsq,invr2_ij
    real(kind=dp)                :: r2_ij,tmpE,ir6,ir12
    real(kind=dp),dimension(3)   :: tmpvect
    real(kind=dp),dimension(3)   :: ilj,jlj

    integer :: jmol,ln ! loop counters
    integer :: ji

    ! Assume orthorhomic boundary conditions
    La = sqrt(dot_product(hmatrix(:,1),hmatrix(:,1)))
    Lb = sqrt(dot_product(hmatrix(:,2),hmatrix(:,2)))
    Lc = sqrt(dot_product(hmatrix(:,3),hmatrix(:,3)))
    
    invLa = 1.0_dp/La
    invLb = 1.0_dp/Lb
    invLc = 1.0_dp/Lc
    
    Evdw  = 0.0_dp

    rcsq = lj_cut*lj_cut

    !------------------------------------------------!
    !         Lennard-Jones model                    !
    !------------------------------------------------!

    ilj(:) = pos(:,imol)

    do ln = 1,lattice(ilat)%nn(imol) ! loop over other molecule jmol
       
       jmol = lattice(ilat)%jn(ln,imol)   ! molecule
       
       jlj(:) = pos(:,jmol)

       ! compute separation vector
       tmpvect = jlj(:) - ilj(:)
       jlj(:) = pos(:,jmol)

       ! compute separation vector
       tmpvect = jlj(:) - ilj(:)

       tmpvect(1) = tmpvect(1) - La*anint(tmpvect(1)*invLa)
       tmpvect(2) = tmpvect(2) - Lb*anint(tmpvect(2)*invLb)
       tmpvect(3) = tmpvect(3) - Lc*anint(tmpvect(3)*invLc)
       
       r2_ij   = tmpvect(1)**2 + tmpvect(2)**2 + tmpvect(3)**2
             
       ! compute interactions if in range with the current image
       if ( r2_ij < rcsq ) then

          invr2_ij = 1.0_dp/r2_ij
          ir6  = invr2_ij*invr2_ij*invr2_ij
          ir12 = ir6*ir6 
          
          ! Pair interaction (COMPUTE tmpE)
          tmpE = sigma12*ir12 - sigma6*ir6
          
          Evdw  = Evdw + tmpE
             
          
       end if ! jmol within range of jmol inside image ji
       
    end do  ! end loop over neighbours of imol
    
    ! set return value
    compute_local_energy_mi = Evdw*4.0_dp*lj_epsilon
    
    return 

  end function compute_local_energy_mi

  
  real(kind=dp) function compute_local_energy_iv(imol,n,d,pos,dh2,dh1,hmatrix,ilat) bind(c)
    !------------------------------------------------------------------------------!
    ! Calculates the real-space contribution to the energy due particle dmol.      !
    ! To be used when computing the changes in energy due to a trial move.         !
    ! This version uses stored vectors to consider multiple images.                !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2017                                                     !
    !------------------------------------------------------------------------------    
    implicit none
    integer,value,intent(in) :: imol,d,n,dh1,dh2,ilat
    real(kind=dp),intent(in),dimension(d,n) :: pos
    real(kind=dp),intent(in),dimension(dh1,dh2) :: hmatrix
    
    ! local variables
    real(kind=dp)                :: Evdw,rcsq,invr2_ij
    real(kind=dp)                :: r2_ij,tmpE,ir6,ir12
    real(kind=dp),dimension(3)   :: tmpvect
    real(kind=dp),dimension(3)   :: ilj,jlj

    integer :: jmol,ln ! loop counters
    integer :: ji
    
    Evdw  = 0.0_dp

    rcsq = lj_cut*lj_cut

    !------------------------------------------------!
    !         Lennard-Jones model                    !
    !------------------------------------------------!

    ilj(:) = pos(:,imol)

    do ln = 1,lattice(ilat)%nn(imol) ! loop over other molecule jmol
       
       jmol = lattice(ilat)%jn(ln,imol)   ! molecule
       ji   = lattice(ilat)%vn(ln,imol)   ! image
       
       jlj(:) = pos(:,jmol) + lattice(ilat)%ivect(:,ji)

       ! compute separation vector
       tmpvect = jlj(:) - ilj(:)
       r2_ij   = tmpvect(1)**2 + tmpvect(2)**2 + tmpvect(3)**2
             
       ! compute interactions if in range with the current image
       if ( r2_ij < rcsq ) then

          invr2_ij = 1.0_dp/r2_ij
          ir6  = invr2_ij*invr2_ij*invr2_ij
          ir12 = ir6*ir6 
          
          ! Pair interaction (COMPUTE tmpE)
          tmpE = sigma12*ir12 - sigma6*ir6
          
          Evdw  = Evdw + tmpE
             
          
       end if ! jmol within range of jmol inside image ji
       
    end do  ! end loop over neighbours of imol
    
    ! set return value
    compute_local_energy_iv = Evdw*4.0_dp*lj_epsilon
    
    return 

  end function compute_local_energy_iv


  real(kind=dp) function compute_lattice_energy_mi(n,d,pos,dh2,dh1,hmatrix,ilat) bind(c)
    !------------------------------------------------------------------------------!
    ! Computes the energy of the entire simulation cell, using a pair potential.   !
    ! To be used when sampling the energy of the entire system or                  !
    ! when computing the energy change from a volume/cell trial move.              !
    ! This version uses the minimum image convention.                              !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2017                                                     !
    !------------------------------------------------------------------------------!
    use util, only : Util_determinant
    implicit none
    integer,value,intent(in) :: d,n,dh1,dh2,ilat
    real(kind=dp),intent(in),dimension(d,n) :: pos
    real(kind=dp),intent(in),dimension(dh1,dh2) :: hmatrix
    
    real(kind=dp) :: La,Lb,Lc,invLa,invLb,invLc
    real(kind=dp),dimension(3)   :: tmpvect
    real(kind=dp)                :: r2_ij,invr2_ij
    real(kind=dp)                :: Evdw
    real(kind=dp)                :: tmpE,rcsq,ir6,ir12
    real(kind=dp)                :: lattice_energy ! jmr

    real(kind=dp),dimension(3) :: ilj,jlj

    integer :: imol,jmol ! loop counters
    integer :: ji,ln

    ! Assume orthorhomic boundary conditions
    La = sqrt(dot_product(hmatrix(:,1),hmatrix(:,1)))
    Lb = sqrt(dot_product(hmatrix(:,2),hmatrix(:,2)))
    Lc = sqrt(dot_product(hmatrix(:,3),hmatrix(:,3)))
    
    invLa = 1.0_dp/La
    invLb = 1.0_dp/Lb
    invLc = 1.0_dp/Lc
    
    Evdw  = 0.0_dp

    rcsq = lj_cut*lj_cut

    !------------------------------------------------!
    !         Lennard-Jones model                    !
    !------------------------------------------------!
    do imol = 1,lattice(ilat)%natoms  ! loop over central molecule imol

       ilj(:) = pos(:,imol)

       do ln = 1,lattice(ilat)%nn(imol) ! loop over other molecule jmol

          jmol = lattice(ilat)%jn(ln,imol)   ! molecule

          jlj(:) = pos(:,jmol)

          ! compute separation vector
          tmpvect = jlj(:) - ilj(:)

          tmpvect(1) = tmpvect(1) - La*anint(tmpvect(1)*invLa)
          tmpvect(2) = tmpvect(2) - Lb*anint(tmpvect(2)*invLb)
          tmpvect(3) = tmpvect(3) - Lc*anint(tmpvect(3)*invLc)
          
          r2_ij   = tmpvect(1)**2 + tmpvect(2)**2 + tmpvect(3)**2
             
          ! compute interactions if in range with the current image
          if ( r2_ij < rcsq ) then

             invr2_ij = 1.0_dp/r2_ij
             ir6  = invr2_ij*invr2_ij*invr2_ij
             ir12 = ir6*ir6 
                                   
             ! Pair interaction (COMPUTE tmpE)
             tmpE = sigma12*ir12 - sigma6*ir6             
             Evdw  = Evdw + tmpE
             
                   
          end if ! jmol within range of jmol inside image ji
                       
       end do  ! end loop over neighbours of imol

    end do ! end loop over imol

    ! Correct for double counting (could avoid double counting but that
    ! is difficult if not using minimum image as here).
    lattice_energy = Evdw*4.0_dp*lj_epsilon*0.5_dp

    ! Add in long range correction
    lattice_energy = lattice_energy + lj_lrc*real(lattice(ilat)%natoms*lattice(ilat)%natoms,kind=dp)/util_determinant(hmatrix)
    
    compute_lattice_energy_mi = lattice_energy
    
    return 

  end function compute_lattice_energy_mi

  
  real(kind=dp) function compute_lattice_energy_iv(n,d,pos,dh2,dh1,hmatrix,ilat) bind(c)
    !------------------------------------------------------------------------------!
    ! Computes the energy of the entire simulation cell, using a pair potential.   !
    ! To be used when sampling the energy of the entire system or                  !
    ! when computing the energy change from a volume/cell trial move.              !
    ! This version uses stored vectors to consider multiple images.                !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2017                                                     !
    !------------------------------------------------------------------------------!
    use util, only : Util_determinant
    implicit none
    integer,value,intent(in) :: d,n,dh1,dh2,ilat
    real(kind=dp),intent(in),dimension(d,n) :: pos
    real(kind=dp),intent(in),dimension(dh1,dh2) :: hmatrix
    

    real(kind=dp),dimension(3)   :: tmpvect
    real(kind=dp)                :: r2_ij,invr2_ij
    real(kind=dp)                :: Evdw
    real(kind=dp)                :: tmpE,rcsq,ir6,ir12
    real(kind=dp)                :: lattice_energy ! jmr

    real(kind=dp),dimension(3) :: ilj,jlj

    integer :: imol,jmol ! loop counters
    integer :: ji,ln
    
    Evdw  = 0.0_dp

    rcsq = lj_cut*lj_cut

    !------------------------------------------------!
    !         Lennard-Jones model                    !
    !------------------------------------------------!
    do imol = 1,lattice(ilat)%natoms  ! loop over central molecule imol

       ilj(:) = pos(:,imol)

       do ln = 1,lattice(ilat)%nn(imol) ! loop over other molecule jmol

          jmol = lattice(ilat)%jn(ln,imol)   ! molecule
          ji   = lattice(ilat)%vn(ln,imol)   ! image

          jlj(:) = pos(:,jmol) + lattice(ilat)%ivect(:,ji)

          ! compute separation vector
          tmpvect = jlj(:) - ilj(:)
          r2_ij   = tmpvect(1)**2 + tmpvect(2)**2 + tmpvect(3)**2
             
          ! compute interactions if in range with the current image
          if ( r2_ij < rcsq ) then

             invr2_ij = 1.0_dp/r2_ij
             ir6  = invr2_ij*invr2_ij*invr2_ij
             ir12 = ir6*ir6 
                                   
             ! Pair interaction (COMPUTE tmpE)
             tmpE = sigma12*ir12 - sigma6*ir6             
             Evdw  = Evdw + tmpE
             
                   
          end if ! jmol within range of jmol inside image ji
                       
       end do  ! end loop over neighbours of imol

    end do ! end loop over imol

    ! Correct for double counting (could avoid double counting but that
    ! is difficult if not using minimum image as here).
    lattice_energy = Evdw*4.0_dp*lj_epsilon*0.5_dp

    ! Add in long range correction
    lattice_energy = lattice_energy + lj_lrc*real(lattice(ilat)%natoms*lattice(ilat)%natoms,kind=dp)/util_determinant(hmatrix)
    
    compute_lattice_energy_iv = lattice_energy
    
    return 

  end function compute_lattice_energy_iv


   subroutine  compute_neighbour_list(n,d,pos,dh2,dh1,hmatrix,ilat) bind(c)
    !------------------------------------------------------------------------------!
    ! Call either the minimum image or image vector routine to compute the energy  !
    ! of the entire system.                                                        !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2017                                                     !
    !------------------------------------------------------------------------------!    
    implicit none
    integer,value,intent(in) :: d,n,dh1,dh2,ilat
    real(kind=dp),intent(in),dimension(d,n) :: pos
    real(kind=dp),intent(in),dimension(dh1,dh2) :: hmatrix

    integer :: imol
    real(kind=dp) :: loc_energy

    ! Sanity checks
    if (d/=3) stop 'Error in compute_local_real_energy: position array must be 3xn'
    if ((dh1/=3).or.(dh2/=3)) stop 'Error in compute_local_real_energy: cell array must be 3x3'

    if (lattice(ilat)%initialised) then
       if (n/=lattice(ilat)%natoms) stop 'Error in compute_local_real_energy: Number of entries in position array has changed!'
    else
       call create_lattice(n,d,pos,dh2,dh1,hmatrix,ilat)
    end if

    select case (image_flag)
    case (flag_iv)
       call compute_neighbour_list_iv(n,d,pos,dh2,dh1,hmatrix,ilat)
    case (flag_mi)       
       call compute_neighbour_list_mi(n,d,pos,dh2,dh1,hmatrix,ilat)
    case default
       stop 'Error in compute_neighbour_list- unknown image flag'       
    end select

  end subroutine compute_neighbour_list
  
  subroutine compute_neighbour_list_iv(n,d,pos,dh2,dh1,hmatrix,ilat) bind(c)
    !------------------------------------------------------------------------------!
    ! Computes a Verlet neighbour list for use in calculating energy               !
    ! This version uses stored vectors to consider multiple images.                !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2017                                                     !
    !------------------------------------------------------------------------------!
    implicit none
    integer,value,intent(in) :: n,d,dh2,dh1,ilat
    real(kind=dp),dimension(d,n),intent(in) :: pos
    real(kind=dp),dimension(dh1,dh2),intent(in) :: hmatrix

    integer :: imol,jmol,k,ni
    real(kind=dp) :: r2_ij,rn
    real(kind=dp),dimension(3) :: v_ij,ilj,jlj,tmpvect

    ! Hard coded neighbour list cutoff with skin width 
    rn = lj_cut + skin
    
    call compute_ivects(hmatrix,ilat)

    do imol = 1,n

       ilj(:) = pos(:,imol)

       lattice(ilat)%nn(imol) = 0
       do jmol = 1,n

          jlj(:) = pos(:,jmol)

          v_ij(:)   = jlj(:) - ilj(:)

          do k = 1,lattice(ilat)%nivect
             if ( (k==1).and.(jmol==imol) )cycle

             tmpvect(:) = v_ij(:) + lattice(ilat)%ivect(:,k) ! apply image             
             r2_ij = dot_product(tmpvect,tmpvect)

             if (r2_ij<rn*rn) then
                lattice(ilat)%nn(imol) = lattice(ilat)%nn(imol) + 1
                ni = lattice(ilat)%nn(imol)        ! number of neighbours of imol
                lattice(ilat)%jn(ni,imol) = jmol   ! jmol is a neighbour of imol
                lattice(ilat)%vn(ni,imol) = k      ! in image k
             end if

          end do
          if (lattice(ilat)%nn(imol) > lattice(ilat)%maxneigh) stop 'Error in compute_neighbour_list: maxneigh exceeded'
       end do
       
       !write(0,'("Molecule ",I5," has ",I5," neighbours")')imol,nn(imol,ils)

    end do


  end subroutine compute_neighbour_list_iv


  subroutine compute_neighbour_list_mi(n,d,pos,dh2,dh1,hmatrix,ilat) bind(c)
    !------------------------------------------------------------------------------!
    ! Computes a Verlet neighbour list for use in calculating energy               !
    ! This version uses the minimum image convention.                              !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2017                                                     !
    !------------------------------------------------------------------------------!
    implicit none
    integer,value,intent(in) :: n,d,dh2,dh1,ilat
    real(kind=dp),dimension(d,n),intent(in) :: pos
    real(kind=dp),dimension(dh1,dh2),intent(in) :: hmatrix

    integer :: imol,jmol,k,ni
    real(kind=dp) :: r2_ij,rn
    real(kind=dp) :: La,Lb,Lc,invLa,invLb,invLc
    real(kind=dp),dimension(3) :: v_ij,ilj,jlj,tmpvect

    ! Hard coded neighbour list cutoff with skin width 
    rn = lj_cut + skin

    ! Assume orthorhomic boundary conditions
    La = sqrt(dot_product(hmatrix(:,1),hmatrix(:,1)))
    Lb = sqrt(dot_product(hmatrix(:,2),hmatrix(:,2)))
    Lc = sqrt(dot_product(hmatrix(:,3),hmatrix(:,3)))
    
    invLa = 1.0_dp/La
    invLb = 1.0_dp/Lb
    invLc = 1.0_dp/Lc
    
    do imol = 1,n

       ilj(:) = pos(:,imol)

       lattice(ilat)%nn(imol) = 0

       do jmol = 1,n

          if (imol==jmol) cycle

          jlj(:) = pos(:,jmol)

          v_ij(:)   = jlj(:) - ilj(:)
          
          tmpvect(1) = v_ij(1) - La*anint(v_ij(1)*invLa)
          tmpvect(2) = v_ij(2) - Lb*anint(v_ij(2)*invLb)
          tmpvect(3) = v_ij(3) - Lc*anint(v_ij(3)*invLc)
          
          r2_ij = dot_product(tmpvect,tmpvect)
          
          if (r2_ij<rn*rn) then
             lattice(ilat)%nn(imol) = lattice(ilat)%nn(imol) + 1
             ni = lattice(ilat)%nn(imol)        ! number of neighbours of imol
             lattice(ilat)%jn(ni,imol) = jmol   ! jmol is a neighbour of imol
          end if
             
          if (lattice(ilat)%nn(imol) > lattice(ilat)%maxneigh) stop 'Error in compute_neighbour_list: maxneigh exceeded'

       end do
       
    end do


  end subroutine compute_neighbour_list_mi
  

end module energy
