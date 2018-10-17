! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!=============================================================================!
!                               G A U S S                                     !
!=============================================================================!
! Contains C-compatible Fortran 2008 routines which calculate the energy of a !
! system of particles interacting via the standard Gaussian core model.       !
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
  integer,parameter  :: gauss_num_lattices=2
  
  ! Constants defining the Gaussian core potential
  real(kind=dp),save :: gauss_epsilon = 1.0_dp
  real(kind=dp),save :: gauss_sigma   = 1.0_dp
  real(kind=dp),save :: gauss_cut     = 3.7_dp
  real(kind=dp),save :: gauss_skin    = 0.5_dp

  type lattice_t
     
     ! Set when neighbours list and ivect arrays have been
     ! calculated for each lattice
     logical :: initialised = .false.              

     ! Number of atoms in the lattices  (for consistency checking)
     integer :: natoms = -1
  
     ! Current energy due to the model Hamiltonian
     real(kind=dp) :: energy = 0.0_dp

     ! Lattice translation vectors
     integer :: nivect
     real(kind=dp),allocatable,dimension(:,:) :: ivect
  
     ! Neighbour list. 
     integer :: maxneigh = 50 ! Maxneigh controls max no. neighbours per atom
     integer,allocatable,dimension(:)   :: nn      ! list of neighbour entries
     integer,allocatable,dimension(:,:) :: jn,vn   ! number and image of neighbout     

  end type lattice_t

  ! Array of these lattice variables
  type(lattice_t),dimension(0:gauss_num_lattices-1) :: lattice

contains

  subroutine create_lattice(n,d,pos,dh2,dh1,hmatrix,ilat) bind(c)
    !==============================================================================!
    ! Initialises a lattice object:                                                !
    !           pos(d,n) :: positions of n particles in d dimensions.              !
    !   hmatrix(dh2,dh1) :: dh2 x dh1 matrix containing the lattice vectors        !
    !               ilat :: index of lattice to initialise 
    !                                                                              !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2017                                                     !
    !------------------------------------------------------------------------------!
    use util,       only : util_determinant,util_images,util_recipmatrix,pi
    implicit none

    ! input variables
    integer,value,intent(in) :: n,d,dh2,dh1,ilat
    real(kind=dp),intent(in),dimension(d,n) :: pos
    real(kind=dp),intent(in),dimension(dh1,dh2) :: hmatrix

    ! local variables
    integer :: im,jm,km    ! loop counters
    integer :: ierr        ! error flag
    real(kind=dp) :: density

    !write(0,'("Entered create_lattice for lattice ",I5)')ilat

    if (ilat>=gauss_num_lattices) stop 'Error in create_lattice : ilat out of range'
    
    ! Set number of atoms
    lattice(ilat)%natoms = n

    !===============!
    ! Image vectors !
    !===============!
    
    ! Find out how many lattice vectors we need to consider all images within cutoff
    im = floor((gauss_cut+1.0_dp)/sqrt(dot_product(hmatrix(:,1),hmatrix(:,1))))+1
    jm = floor((gauss_cut+1.0_dp)/sqrt(dot_product(hmatrix(:,2),hmatrix(:,2))))+1
    km = floor((gauss_cut+1.0_dp)/sqrt(dot_product(hmatrix(:,3),hmatrix(:,3))))+1
    
    lattice(ilat)%nivect = (2*im+1)*(2*jm+1)*(2*km+1)

    if (allocated(lattice(ilat)%ivect)) then
       deallocate(lattice(ilat)%ivect, stat=ierr)
       if (ierr/=0) stop 'Error in create_lattice : Could not deallocate ivect'
    end if
    allocate(lattice(ilat)%ivect(1:d,1:lattice(ilat)%nivect),stat=ierr)
    if (ierr/=0) stop 'Error in create_lattice : Could not allocate ivect'
       
    ! compute current image translation vectors
    call compute_image_vectors(hmatrix,ilat)

    !============!
    ! Neighbours !
    !============!
    
    ! Compute maxneigh based on density and a safety factor of 20%
    density = lattice(ilat)%natoms/abs(util_determinant(hmatrix))
    lattice(ilat)%maxneigh = ceiling((density*4.0*pi*(gauss_cut+gauss_skin)**3)/3.0_dp * 1.2_dp)

    if (allocated(lattice(ilat)%nn)) then
       deallocate(lattice(ilat)%nn,stat=ierr)    
       if (ierr/=0) stop 'Error in create_lattice : Could not deallocate nn'
    end if
    allocate(lattice(ilat)%nn(1:lattice(ilat)%natoms),stat=ierr)
    if (ierr/=0) stop 'Error in create_lattice : Could not allocate nn array'
    
    if (allocated(lattice(ilat)%vn)) then
       deallocate(lattice(ilat)%vn,stat=ierr)    
       if (ierr/=0) stop 'Error in create_lattice : Could not deallocate vn'
    end if
    allocate(lattice(ilat)%vn(1:lattice(ilat)%maxneigh,1:lattice(ilat)%natoms),stat=ierr)
    if (ierr/=0) stop 'Error in create_lattice : Could not allocate vn array'
       
    if (allocated(lattice(ilat)%jn)) then
       deallocate(lattice(ilat)%jn,stat=ierr)    
       if (ierr/=0) stop 'Error in create_lattice : Could not deallocate jn'
    end if
    allocate(lattice(ilat)%jn(1:lattice(ilat)%maxneigh,1:lattice(ilat)%natoms),stat=ierr)
    if (ierr/=0) stop 'Error in create_lattice : Could not allocate jn array'


    ! Set initialised flag
    lattice(ilat)%initialised = .true.
    
    ! Make sure first neighbour list is calculated 
    call compute_neighbour_list(lattice(ilat)%natoms,d,pos,d,d,hmatrix,ilat)

    !write(*,'("Initialised lattice ",I5)')ilat
    
    ! Make sure energy if up-to-date
    lattice(ilat)%energy = compute_lattice_energy(n,d,pos,d,d,hmatrix,ilat)
    
    return

  end subroutine create_lattice

  subroutine compute_image_vectors(hmatrix,ilat) 
    !------------------------------------------------------------------------------!
    ! Computes the translation vectors needed to include all images such that each !
    ! atom sees all the images (including those of itself) within the cut-off      !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2017                                                     !
    !------------------------------------------------------------------------------!
    implicit none
    integer,intent(in) :: ilat
    real(kind=dp),dimension(3,3),intent(in) :: hmatrix
 
    !loop counters
    integer :: icell,jcell,kcell,im,jm,km,k
    real(kind=dp),dimension(3) :: sx,sy,sz

    if (ilat>=gauss_num_lattices) stop 'Error in create_lattice : ilat out of range'
    
    im = floor(gauss_cut/sqrt(dot_product(hmatrix(:,1),hmatrix(:,1))))+1
    jm = floor(gauss_cut/sqrt(dot_product(hmatrix(:,2),hmatrix(:,2))))+1
    km = floor(gauss_cut/sqrt(dot_product(hmatrix(:,3),hmatrix(:,3))))+1
    
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

  end subroutine compute_image_vectors


  real(kind=dp) function compute_local_energy(dmol,n,d,pos,dh2,dh1,hmatrix,ilat) bind(c)
    !------------------------------------------------------------------------------!
    ! Calculates the real-space contribution to the energy due particle dmol.      !
    ! To be used when computing the changes in energy due to a trial move.         !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2017                                                     !
    !------------------------------------------------------------------------------    
    implicit none
    integer,value,intent(in) :: dmol,d,n,dh1,dh2,ilat
    real(kind=dp),intent(in),dimension(d,n) :: pos
    real(kind=dp),intent(in),dimension(dh1,dh2) :: hmatrix
    
    ! local variables
    real(kind=dp)                :: Evdw,rcsq,invsigmasq
    real(kind=dp)                :: r2_ij,tmpE
    real(kind=dp),dimension(3)   :: tmpvect
    real(kind=dp),dimension(3)   :: ilj,jlj

    integer :: jmol,ln ! loop counters
    integer :: ji,imol

    ! Sanity checks
    if (d/=3) stop 'Error in compute_local_real_energy: position array must be 3xn'
    if ((dh1/=3).or.(dh2/=3)) stop 'Error in compute_local_real_energy: cell array must be 3x3'

    if (ilat>=gauss_num_lattices) stop 'Error in create_lattice : ilat out of range'
    
    if (lattice(ilat)%initialised) then
       if (n/=lattice(ilat)%natoms) stop 'Error in compute_local_real_energy: Number of entries in position array has changed!'
    else
       call create_lattice(n,d,pos,dh2,dh1,hmatrix,ilat)
    end if

    ! We assume this routine will only ever be called from C/Python (never Fortran)
    ! and so we correct imol into 1-based indexing
    imol = dmol + 1
    if ( (imol<1).or.(imol>lattice(ilat)%natoms) ) stop 'Error in compute_local_real_energy : atom index out of range!'
    
    Evdw  = 0.0_dp

    rcsq = gauss_cut*gauss_cut
    invsigmasq = 1.0_dp/(gauss_sigma*gauss_sigma)

    !------------------------------------------------!
    !         Gaussian core model                    !
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
          
          ! Pair interaction (COMPUTE tmpE)
          tmpE = exp(-r2_ij*invsigmasq)
          
          Evdw  = Evdw + tmpE
             
          
       end if ! jmol within range of jmol inside image ji
       
    end do  ! end loop over neighbours of imol
    
    ! set return value
    compute_local_energy = Evdw*gauss_epsilon
    
    return 

  end function compute_local_energy


  real(kind=dp) function compute_lattice_energy(n,d,pos,dh2,dh1,hmatrix,ilat) bind(c)
    !------------------------------------------------------------------------------!
    ! Computes the energy of the entire simulation cell, using a pair potential.   !
    ! To be used when sampling the energy of the entire system or                  !
    ! when computing the energy change from a volume/cell trial move.              !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2017                                                     !
    !------------------------------------------------------------------------------!
    implicit none
    integer,value,intent(in) :: d,n,dh1,dh2,ilat
    real(kind=dp),intent(in),dimension(d,n) :: pos
    real(kind=dp),intent(in),dimension(dh1,dh2) :: hmatrix
    
    real(kind=dp),dimension(3)   :: tmpvect
    real(kind=dp)                :: r2_ij
    real(kind=dp)                :: Evdw
    real(kind=dp)                :: tmpE,rcsq,invsigmasq

    real(kind=dp),dimension(3) :: ilj,jlj

    integer :: imol,jmol ! loop counters
    integer :: ji,ln

    ! Sanity checks
    if (d/=3) stop 'Error in compute_model_energy: position array must be 3xn'
    if ((dh1/=3).or.(dh2/=3)) stop 'Error in compute_model_energy: cell array must be 3x3'

    if (ilat>=gauss_num_lattices) stop 'Error in create_lattice : ilat out of range'
    
    if (.not.lattice(ilat)%initialised) then
       call create_lattice(n,d,pos,dh2,dh1,hmatrix,ilat)
    else
       if (n/=lattice(ilat)%natoms) stop 'Error in compute_model_energy: Number of entries in position array has changed!'
    end if
    
    Evdw  = 0.0_dp

    rcsq = gauss_cut*gauss_cut
    invsigmasq = 1.0_dp/(gauss_sigma*gauss_sigma)

    !------------------------------------------------!
    !         Gaussian core model                    !
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
                
             ! Pair interaction (COMPUTE tmpE)
             tmpE = exp(-r2_ij*invsigmasq)             
             Evdw  = Evdw + tmpE
             
                   
          end if ! jmol within range of jmol inside image ji
                       
       end do  ! end loop over neighbours of imol

    end do ! end loop over imol

    ! Correct for double counting (could avoid double counting but that
    ! is difficult if not using minimum image as here).
    lattice(ilat)%energy   = Evdw*gauss_epsilon*0.5_dp
    compute_lattice_energy = Evdw*gauss_epsilon*0.5_dp
    
    return 

  end function compute_lattice_energy

  subroutine compute_neighbour_list(n,d,pos,dh2,dh1,hmatrix,ilat) bind(c)
    !------------------------------------------------------------------------------!
    ! Computes a Verlet neighbour list for use in calculating energy               !
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

    ! Hard coded neighbour list cutoff, should probably be read from input
    rn = gauss_cut + gauss_skin

    ! Sanity checks
    if (d/=3) stop 'Error in compute_model_energy: position array must be 3xn'
    if ((dh1/=3).or.(dh2/=3)) stop 'Error in compute_model_energy: cell array must be 3x3'

    if (ilat>=gauss_num_lattices) stop 'Error in create_lattice : ilat out of range'
    
    if (.not.lattice(ilat)%initialised) then
       call create_lattice(n,d,pos,dh2,dh1,hmatrix,ilat)
    else
       if (n/=lattice(ilat)%natoms) stop 'Error in compute_model_energy: Number of entries in position array has changed!'
    end if

    ! Make sure image vectors are current
    call compute_image_vectors(hmatrix,ilat)

    do imol = 1,lattice(ilat)%natoms

       ilj(:) = pos(:,imol)

       lattice(ilat)%nn(imol) = 0
       
       do jmol = 1,lattice(ilat)%natoms

          jlj(:) = pos(:,jmol)

          v_ij(:) = jlj(:) - ilj(:)

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
          if (lattice(ilat)%nn(imol) > lattice(ilat)%maxneigh) &
               stop 'Error in compute_neighbour _list : maximum neighbours exceeded!'
       end do
       
       !write(0,'("Molecule ",I5," has ",I5," neighbours")')imol,nn(imol,ils)
       
    end do


  end subroutine compute_neighbour_list


end module energy
