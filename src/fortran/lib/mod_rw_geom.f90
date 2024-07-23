!!!#############################################################################
!!! Code written by Ned Thaddeus Taylor and Francis Huw Davies
!!! Code part of the ARTEMIS group (Hepplestone research group).
!!! Think Hepplestone, think HRG.
!!!#############################################################################
!!! Module made to read and write structure input files
!!! Currently supports:
!!!    -VASP
!!!    -Quantum Espresso
!!!    -CASTEP
!!!    -xyz (read only)
!!!#############################################################################
module rw_geom
  use constants, only: pi,real12
  use misc_raffle, only: to_upper,jump,Icount
  use misc_linalg, only: LUinv,modu
  implicit none


  private

  public :: igeom_input,igeom_output
  public :: bas_type, spec_type
  public :: clone_bas
  public :: convert_bas
  public :: geom_read,geom_write


  integer :: igeom_input=1,igeom_output=1

  type spec_type
     real(real12), allocatable ,dimension(:,:) :: atom
     real(real12) :: mass
     real(real12) :: charge
     character(len=3) :: name
     integer :: num
  end type spec_type
  type bas_type
     type(spec_type), allocatable, dimension(:) :: spec
     integer :: nspec
     integer :: natom
     real(real12) :: energy
     real(real12) :: lat(3,3) = 0._real12
     logical :: lcart = .false.
     logical, dimension(3) :: pbc = .true.
     character(len=1024) :: sysname
   contains
     procedure, pass(this) :: allocate_species
  end type bas_type

  

!!!updated  2024/06/25


contains

  subroutine allocate_species(this, num_species, species_symbols, species_count, atoms)
    implicit none
    class(bas_type), intent(inout) :: this
    integer, intent(in), optional :: num_species
    character(3), dimension(:), intent(in), optional :: species_symbols
    integer, dimension(:), intent(in), optional :: species_count
    real(real12), dimension(:,:), intent(in), optional :: atoms
    
    integer :: i, istart, iend

    if(present(num_species)) this%nspec = num_species

    if(allocated(this%spec)) deallocate(this%spec)
    allocate(this%spec(this%nspec))

    species_check: if(present(species_symbols))then
       if(size(species_symbols).ne.this%nspec) exit species_check
       this%spec(:)%name = species_symbols
    end if species_check

    natom_check: if(present(species_count))then
       if(size(species_count).ne.this%nspec) exit natom_check
       this%spec(:)%num = species_count
       istart = 1
       do i = 1, this%nspec
          iend = istart + this%spec(i)%num - 1
          allocate(this%spec(i)%atom(this%spec(i)%num,3))
          if(present(atoms))then
             this%spec(i)%atom = atoms(istart:iend,:3)
          end if
          istart = iend + 1
       end do
    end if natom_check
       
  end subroutine allocate_species

  
!!!#############################################################################
!!! sets up the name of output files and subroutines to read files
!!!#############################################################################
  subroutine geom_read(UNIT,basis,length)
    implicit none
    integer, intent(in) :: UNIT
    type(bas_type), intent(out) :: basis
    integer, optional, intent(in) :: length

    integer :: dim,i

    dim=3
    if(present(length)) dim=length

    select case(igeom_input)
    case(1)
       call VASP_geom_read(UNIT,basis,dim)
    case(2)
       call CASTEP_geom_read(UNIT,basis,dim)
    case(3)
       call QE_geom_read(UNIT,basis,dim)
    case(4)
       stop "ERROR: Not yet set up for CRYSTAL"
    case(5)
       call XYZ_geom_read(UNIT,basis,dim)
       write(0,'("WARNING: XYZ file format does not contain lattice data")')
    case(6)
       call extXYZ_geom_read(UNIT,basis,dim)
    end select
    if(dim.eq.4)then
       do i=1,basis%nspec
         basis%spec(i)%atom(:,4)=1._real12
       end do
    end if


  end subroutine geom_read
!!!#############################################################################


!!!#############################################################################
!!! sets up the name of output files and subroutines to read files
!!!#############################################################################
  subroutine geom_write(UNIT,basis)
    implicit none
    integer, intent(in) :: UNIT
    type(bas_type), intent(in) :: basis

!!! MAKE IT CHANGE HERE IF USER SPECIFIES LCART OR NOT
!!! AND GIVE IT THE CASTEP AND QE OPTION OF LABC !!!

    select case(igeom_output)
    case(1)
       call VASP_geom_write(UNIT,basis)
    case(2)
       call CASTEP_geom_write(UNIT,basis)
    case(3)
       call QE_geom_write(UNIT,basis)
    case(4)
       write(0,'("ERROR: ARTEMIS not yet set up for CRYSTAL")')
       stop
    case(5)   
       call XYZ_geom_write(UNIT,basis)
    case(6)
       call extXYZ_geom_write(UNIT,basis)
    end select


  end subroutine geom_write
!!!#############################################################################


!!!#############################################################################
!!! read the POSCAR or CONTCAR file
!!!#############################################################################
 subroutine VASP_geom_read(UNIT,basis,length)
   implicit none
   integer, intent(in) :: UNIT
   type(bas_type), intent(out) :: basis
   integer, intent(in), optional :: length

   integer :: pos,count,Reason
   real(real12) :: scal 
   character(len=100) :: lspec
   character(len=1024) :: buffer
   real(real12), dimension(3,3) :: reclat
   integer :: i,j,k,dim


!!!-----------------------------------------------------------------------------
!!! determines dimension of basis (include translation dimension for symmetry?)
!!!-----------------------------------------------------------------------------
   if(present(length))then
      dim = length
   else
      dim = 3
   end if

!!!-----------------------------------------------------------------------------
!!! read system name
!!!-----------------------------------------------------------------------------
   read(UNIT,'(A)',iostat=Reason) basis%sysname
   if(Reason.ne.0)then
      write(0,'(" The file is not in POSCAR format.")')
      write(0,'(" Exiting code ...")')
      call exit()
   end if
   read(UNIT,*) scal


!!!-----------------------------------------------------------------------------
!!! read lattice
!!!-----------------------------------------------------------------------------
   do i=1,3
      read(UNIT,*) (basis%lat(i,j),j=1,3)
   end do
   basis%lat=scal*basis%lat
   

!!!-----------------------------------------------------------------------------
!!! read species names and number of each atomic species
!!!-----------------------------------------------------------------------------
   read(UNIT,'(A)') lspec
   basis%nspec=Icount(lspec)
   allocate(basis%spec(basis%nspec))
   if(verify(lspec,' 0123456789').ne.0) then
      count=0;pos=1
      speccount: do
         i=verify(lspec(pos:), ' ')
         if (i.eq.0) exit speccount
         count=count+1
         pos=i+pos-1
         i=scan(lspec(pos:), ' ')
         if (i.eq.0) exit speccount
         basis%spec(count)%name=lspec(pos:pos+i-1)
         pos=i+pos-1
      end do speccount

      read(UNIT,*) (basis%spec(j)%num,j=1,basis%nspec)
   else !only numbers
      do count=1,basis%nspec
         write(basis%spec(count)%name,'(I0)') count
      end do
      read(lspec,*) (basis%spec(j)%num,j=1,basis%nspec)
   end if


!!!-----------------------------------------------------------------------------
!!! determines whether input basis is in direct or cartesian coordinates
!!!-----------------------------------------------------------------------------
   basis%lcart=.false.
   read(UNIT,'(A)') buffer
   if(verify(trim(buffer),'Direct').eq.0) basis%lcart=.false.
   if(verify(trim(buffer),'Cartesian').eq.0) then
      write(0,*) "NOT SURE IF CARTESIAN COORDINATES ARE SUPPORTED YET!"
      write(0,*) "PLEASE CHECK COORDINATES"
      basis%lcart=.true.
   end if


!!!-----------------------------------------------------------------------------
!!! read basis
!!!-----------------------------------------------------------------------------
   do i=1,basis%nspec
      allocate(basis%spec(i)%atom(basis%spec(i)%num,dim))
      basis%spec(i)%atom(:,:)=0._real12
      do j=1,basis%spec(i)%num
         read(UNIT,*) (basis%spec(i)%atom(j,k),k=1,3)
      end do
   end do


!!!-----------------------------------------------------------------------------
!!! convert basis if in cartesian coordinates
!!!-----------------------------------------------------------------------------
   if(basis%lcart)then
      reclat=transpose(LUinv(basis%lat))*2._real12*pi
      basis=convert_bas(basis,reclat)
   end if


!!!-----------------------------------------------------------------------------
!!! normalise basis to between 0 and 1 in direct coordinates
!!!-----------------------------------------------------------------------------
   do i=1,basis%nspec
      do j=1,basis%spec(i)%num
         do k=1,3
            basis%spec(i)%atom(j,k)=&
                 basis%spec(i)%atom(j,k)-floor(basis%spec(i)%atom(j,k))
         end do
      end do
   end do
   basis%natom=sum(basis%spec(:)%num)


 end subroutine VASP_geom_read
!!!#############################################################################


!!!#############################################################################
!!! writes out the structure in vasp poscar style format
!!!#############################################################################
 subroutine VASP_geom_write(UNIT,basis,lcart)
   implicit none
   integer, intent(in) :: UNIT
   type(bas_type), intent(in) :: basis
   logical, intent(in), optional :: lcart

   integer :: i,j
   character(100) :: fmt,string

   string="Direct"
   if(present(lcart))then
      if(lcart) string="Cartesian"
   end if

   write(UNIT,'(A)') trim(adjustl(basis%sysname))
   write(UNIT,'(F15.9)') 1._real12
   do i=1,3
      write(UNIT,'(3(F15.9))') basis%lat(i,:)
   end do
   write(fmt,'("(",I0,"(A,1X))")') basis%nspec
   write(UNIT,trim(adjustl(fmt))) (adjustl(basis%spec(j)%name),j=1,basis%nspec)
   write(fmt,'("(",I0,"(I0,5X))")') basis%nspec
   write(UNIT,trim(adjustl(fmt))) (basis%spec(j)%num,j=1,basis%nspec)
   write(UNIT,'(A)') trim(adjustl(string))
   do i=1,basis%nspec
      do j=1,basis%spec(i)%num
         write(UNIT,'(3(F15.9))') basis%spec(i)%atom(j,1:3)
      end do
   end do
   

 end subroutine VASP_geom_write
!!!#############################################################################


!!!#############################################################################
!!! read the QE geom file
!!!#############################################################################
 subroutine QE_geom_read(UNIT,basis,length)
   implicit none
   integer, intent(in) :: UNIT
   type(bas_type), intent(out) :: basis
   integer, intent(in), optional :: length

   integer Reason,i,j,k,dim,iline
   integer, dimension(1000) :: tmp_natom
   real(real12), dimension(3) :: tmpvec
   real(real12), dimension(3,3) :: reclat
   character(len=5) :: ctmp
   character(len=5), dimension(1000) :: tmp_spec
   character(len=1024) :: buffer,buffer2


!!!-----------------------------------------------------------------------------
!!! determines dimension of basis (include translation dimension for symmetry?)
!!!-----------------------------------------------------------------------------
   if(present(length))then
      dim = length
   else
      dim = 3
   end if

   basis%lcart=.false.
   basis%sysname="Converted_from_geom_file"


!!!-----------------------------------------------------------------------------
!!! read lattice
!!!-----------------------------------------------------------------------------
   rewind UNIT
   cellparam: do
      read(UNIT,'(A)',iostat=Reason) buffer
      if(Reason.ne.0)then
         write(0,'(" An issue with the QE input file format has been encountered.")')
         write(0,'(" Exiting code ...")')
         stop
      end if
      if(index(trim(buffer),"ibrav").ne.0)then
         write(0,'("ERROR: Internal error in QE_geom_read")')
         write(0,'(2X,"Subroutine not yet set up to read IBRAV lattices")')
         stop
      end if
      if(verify("CELL_PARAMETERS",buffer).eq.0) then
         exit cellparam
      end if
   end do cellparam
   do i=1,3
      read(UNIT,*) (basis%lat(i,j),j=1,3)
   end do


!!!-----------------------------------------------------------------------------
!!! determines whether input basis is in direct or cartesian coordinates
!!!-----------------------------------------------------------------------------
   iline=0
   rewind UNIT
   basfind: do
      read(UNIT,'(A)',iostat=Reason) buffer
      iline=iline+1
      if(verify("ATOMIC_POSITIONS",buffer).eq.0)then
         backspace(UNIT)
         read(UNIT,*) buffer,buffer2
         if(verify("crystal",buffer2).eq.0) basis%lcart=.false.
         if(verify("angstrom",buffer2).eq.0) basis%lcart=.true.
         exit basfind
      end if
   end do basfind


!!!-----------------------------------------------------------------------------
!!! read basis
!!!-----------------------------------------------------------------------------
   basis%natom=0
   basis%nspec=0
   tmp_natom=1
   basread: do
      read(UNIT,'(A)',iostat=Reason) buffer
      read(buffer,*) ctmp
      if(Reason.ne.0) exit
      if(trim(ctmp).eq.'') exit
      if(verify(buffer,' 0123456789').eq.0) exit
      basis%natom=basis%natom+1
      if(.not.any(tmp_spec(1:basis%nspec).eq.ctmp))then
         basis%nspec=basis%nspec+1
         tmp_spec(basis%nspec)=ctmp
      else
         where(tmp_spec(1:basis%nspec).eq.ctmp)
            tmp_natom(1:basis%nspec)=tmp_natom(1:basis%nspec)+1
         end where
      end if
   end do basread

   allocate(basis%spec(basis%nspec))
   basis%spec(1:basis%nspec)%name=tmp_spec(1:basis%nspec)
   do i=1,basis%nspec
      basis%spec(i)%num=0
      allocate(basis%spec(i)%atom(tmp_natom(i),dim))
   end do
   
   call jump(UNIT,iline)
   basread2: do i=1,basis%natom
      read(UNIT,*,iostat=Reason) ctmp,tmpvec(1:3)
      do j=1,basis%nspec
         if(basis%spec(j)%name.eq.ctmp)then
            basis%spec(j)%num=basis%spec(j)%num+1
            basis%spec(j)%atom(basis%spec(j)%num,1:3)=tmpvec(1:3)
            exit
         end if
      end do
   end do basread2


!!!-----------------------------------------------------------------------------
!!! convert basis if in cartesian coordinates
!!!-----------------------------------------------------------------------------
   if(basis%lcart)then
      reclat=transpose(LUinv(basis%lat))*2._real12*pi
      basis=convert_bas(basis,reclat)
   end if


!!!-----------------------------------------------------------------------------
!!! normalise basis to between 0 and 1 in direct coordinates
!!!-----------------------------------------------------------------------------
   do i=1,basis%nspec
      do j=1,basis%spec(i)%num
         do k=1,3
            basis%spec(i)%atom(j,k)=&
                 basis%spec(i)%atom(j,k)-floor(basis%spec(i)%atom(j,k))
         end do
      end do
   end do
   basis%natom=sum(basis%spec(:)%num)


 end subroutine QE_geom_read
!!!#############################################################################


!!!#############################################################################
!!! writes out the structure in QE geom style format
!!!#############################################################################
 subroutine QE_geom_write(UNIT,basis,lcart)
   implicit none
   integer, intent(in) :: UNIT
   type(bas_type), intent(in) :: basis
   logical, intent(in), optional :: lcart

   integer :: i,j
   character(10) :: string

   string="crystal"
   if(present(lcart))then
      if(lcart) string="angstrom"
   end if


   write(UNIT,'("CELL_PARAMETERS angstrom")')
   do i=1,3
      write(UNIT,'(3(F15.9))') basis%lat(i,:)
   end do
   write(UNIT,'("ATOMIC_SPECIES")')
   do i=1,basis%nspec
      write(UNIT,'(A)') trim(adjustl(basis%spec(i)%name))
   end do
   write(UNIT,'("ATOMIC_POSITIONS",1X,A)') trim(adjustl(string))
   do i=1,basis%nspec
      do j=1,basis%spec(i)%num
         write(UNIT,'(A5,1X,3(F15.9))') &
              basis%spec(i)%name,basis%spec(i)%atom(j,1:3)
      end do
   end do


 end subroutine QE_geom_write
!!!#############################################################################


!!!#############################################################################
!!! reads atoms from an CASTEP file
!!!#############################################################################
 subroutine CASTEP_geom_read(UNIT,basis,length)
   implicit none
   integer, intent(in) :: UNIT
   type(bas_type), intent(out) :: basis
   integer, intent(in), optional :: length

   integer :: Reason,itmp1
   integer :: i,j,k,dim,iline
   character(len=5) :: ctmp
   character(len=20) :: units
   character(len=200) :: buffer,store
   logical :: labc
   integer, dimension(1000) :: tmp_natom
   real(real12), dimension(3) :: abc,angle,dvtmp1
   real(real12), dimension(3,3) :: reclat
   character(len=5), dimension(1000) :: tmp_spec

   
!!!-----------------------------------------------------------------------------
!!! determines dimension of basis (include translation dimension for symmetry?)
!!!-----------------------------------------------------------------------------
   if(present(length))then
      dim = length
   else
      dim = 3
   end if


!!!-----------------------------------------------------------------------------
!!! reading loop of file
!!!-----------------------------------------------------------------------------
   tmp_spec=""
   tmp_natom=0
   iline=0
   basis%sysname="from CASTEP"
   rewind(UNIT)
   readloop: do
      iline=iline+1
      read(UNIT,'(A)',iostat=Reason) buffer
      if(Reason.ne.0) exit
      buffer=to_upper(buffer)
      if(scan(trim(adjustl(buffer)),'%').ne.1) cycle readloop
      if(index(trim(adjustl(buffer)),'%END').eq.1) cycle readloop
      read(buffer,*) store, buffer
      if(trim(buffer).eq.'') cycle readloop
      !!------------------------------------------------------------------------
      !! read lattice
      !!------------------------------------------------------------------------
      lattice_if: if(index(trim(buffer),"LATTICE").eq.1)then
         if(index(trim(buffer),"ABC").ne.0) labc=.true.
         if(index(trim(buffer),"CART").ne.0) labc=.false.
         store=""
         itmp1=0
         lattice_loop: do 
            itmp1=itmp1+1
            read(UNIT,'(A)',iostat=Reason) buffer
            if(Reason.ne.0) exit lattice_loop
            if(scan(trim(adjustl(buffer)),'%').eq.1) exit lattice_loop
            if(itmp1.eq.5)then
               write(0,'("ERROR: Too many lines in LATTICE block of structure file")')
               stop
            end if
            store=trim(store)//" "//trim(buffer)
         end do lattice_loop
         iline=iline+itmp1

         if(labc)then
            read(store,*) units,(abc(i),i=1,3), (angle(j),j=1,3)
            basis%lat=convert_abc_to_lat(abc,angle,.false.)
         else
            read(store,*) units,(basis%lat(i,:),i=1,3)
         end if
         cycle readloop
      end if lattice_if

      !!------------------------------------------------------------------------
      !! read basis
      !!------------------------------------------------------------------------
      basis_if: if(index(trim(buffer),"POSITIONS").eq.1) then
         if(index(trim(buffer),"ABS").ne.0) basis%lcart=.true.
         if(index(trim(buffer),"FRAC").ne.0) basis%lcart=.false.
         itmp1=0
         basis_loop1: do
            read(UNIT,'(A)',iostat=Reason) buffer
            if(Reason.ne.0) exit basis_loop1
            if(scan(trim(adjustl(buffer)),'%').eq.1) exit basis_loop1
            read(buffer,*) ctmp
            if(trim(ctmp).eq.'') exit
            if(verify(buffer,' 0123456789').eq.0) exit
            basis%natom=basis%natom+1
            if(.not.any(tmp_spec(1:basis%nspec).eq.ctmp))then
               basis%nspec=basis%nspec+1
               tmp_natom(basis%nspec)=1
               tmp_spec(basis%nspec)=ctmp
            else
               where(tmp_spec(1:basis%nspec).eq.ctmp)
                  tmp_natom(1:basis%nspec)=tmp_natom(1:basis%nspec)+1
               end where
            end if
         end do basis_loop1

         allocate(basis%spec(basis%nspec))
         basis%spec(1:basis%nspec)%name=tmp_spec(1:basis%nspec)
         do i=1,basis%nspec
            basis%spec(i)%num=0
            allocate(basis%spec(i)%atom(tmp_natom(i),dim))
         end do

         call jump(UNIT,iline)
         basis_loop2: do i=1,basis%natom
            read(UNIT,'(A)',iostat=Reason) buffer
            if(Reason.ne.0)then
               write(0,'("ERROR: Internal error in assigning the basis")')
               stop
            end if
            read(buffer,*) ctmp,dvtmp1(1:3)
            species_loop: do j=1,basis%nspec
               if(basis%spec(j)%name.eq.ctmp)then
                  basis%spec(j)%num=basis%spec(j)%num+1
                  basis%spec(j)%atom(basis%spec(j)%num,1:3)=dvtmp1(1:3)
                  exit species_loop
               end if
            end do species_loop
         end do basis_loop2

      end if basis_if
   
   end do readloop


!!!-----------------------------------------------------------------------------
!!! convert basis if in cartesian coordinates
!!!-----------------------------------------------------------------------------
   if(basis%lcart)then
      reclat=transpose(LUinv(basis%lat))*2._real12*pi
      basis=convert_bas(basis,reclat)
   end if


!!!-----------------------------------------------------------------------------
!!! normalise basis to between 0 and 1 in direct coordinates
!!!-----------------------------------------------------------------------------
   do i=1,basis%nspec
      do j=1,basis%spec(i)%num
         do k=1,3
            basis%spec(i)%atom(j,k)=&
                 basis%spec(i)%atom(j,k)-floor(basis%spec(i)%atom(j,k))
         end do
      end do
   end do
   basis%natom=sum(basis%spec(:)%num)


   return
 end subroutine CASTEP_geom_read
!!!#############################################################################


!!!#############################################################################
!!! writes lattice and basis in a CASTEP file format
!!!#############################################################################
 subroutine CASTEP_geom_write(UNIT,basis,labc,lcart)
   implicit none
   integer :: i,j,UNIT
   real(real12), dimension(3) :: abc,angle
   type(bas_type) :: basis
   character(4) :: string_lat,string_bas
   logical, intent(in), optional :: labc,lcart


   string_lat="CART"
   if(present(labc))then
      if(labc) string_lat="ABC"
   end if

   string_bas="FRAC"
   if(present(lcart))then
      if(lcart)then
         string_bas="ABS"
         write(0,'("ERROR: Internal error in CASTEP_geom_write")')
         write(0,'(2X,"Subroutine not yet set up to output cartesian &
              &coordinates")')
         stop
      end if
   end if

   write(UNIT,'("%block LATTICE_",A)') trim(string_lat)
   write(UNIT,'("ang")')
   if(present(labc))then
      if(labc)then
         do i=1,3
            abc(i)=modu(basis%lat(i,:))
         end do
         angle(1) = dot_product(basis%lat(2,:),basis%lat(3,:))/(abc(2)*abc(3))
         angle(2) = dot_product(basis%lat(1,:),basis%lat(3,:))/(abc(1)*abc(3))
         angle(3) = dot_product(basis%lat(1,:),basis%lat(2,:))/(abc(1)*abc(2))
         write(UNIT,'(3(F15.9))') abc
         write(UNIT,'(3(F15.9))') angle
         goto 10
      end if
   end if
   do i=1,3
      write(UNIT,'(3(F15.9))') basis%lat(i,:)
   end do

10 write(UNIT,'("%endblock LATTICE_",A)') trim(string_lat)

   write(UNIT,*)
   write(UNIT,'("%block POSITIONS_",A)') trim(string_bas)
   do i=1,basis%nspec
      do j=1,basis%spec(i)%num
         write(UNIT,'(A5,1X,3(F15.9))') basis%spec(i)%name,basis%spec(i)%atom(j,1:3)
      end do
   end do
   write(UNIT,'("%endblock POSITIONS_",A)') trim(string_bas)


 end subroutine CASTEP_geom_write
!!!#############################################################################


!!!#############################################################################
!!! reads atoms from an xyz file
!!!#############################################################################
 subroutine XYZ_geom_read(UNIT,basis,length)
   implicit none
   integer, intent(in) :: UNIT
   type(bas_type), intent(out) :: basis
   integer, intent(in), optional :: length

   integer :: Reason
   integer, allocatable, dimension(:) :: tmp_num
   real(real12), dimension(3) :: vec
   real(real12), allocatable, dimension(:,:,:) :: tmp_bas
   character(len=5) :: ctmp
   character(len=5), allocatable, dimension(:) :: tmp_spec
   integer :: i,j,dim

   dim=3
   if(present(length)) dim=length


   read(UNIT,*,iostat=Reason) basis%natom
   if(Reason.ne.0)then
      write(0,'(" The file is not in xyz format.")')
      write(0,'(" Exiting code ...")')
      call exit()
   end if
   read(UNIT,'(A)',iostat=Reason) basis%sysname


!!!-----------------------------------------------------------------------------
!!! read basis
!!!-----------------------------------------------------------------------------
   allocate(tmp_spec(basis%natom))
   allocate(tmp_num(basis%natom))
   allocate(tmp_bas(basis%natom,basis%natom,dim))
   tmp_num(:)=0
   tmp_spec=""
   tmp_bas=0
   basis%nspec=0
   do i=1,basis%natom
      read(UNIT,*,iostat=Reason) ctmp,vec(1:3)
      if(.not.any(tmp_spec(1:basis%nspec).eq.ctmp))then
         basis%nspec=basis%nspec+1
         tmp_spec(basis%nspec)=ctmp
         tmp_bas(basis%nspec,1,1:3)=vec(1:3)
         tmp_num(basis%nspec)=1
      else
         checkspec: do j=1,basis%nspec
            if(tmp_spec(j).eq.ctmp)then
               tmp_num(j)=tmp_num(j)+1
               tmp_bas(j,tmp_num(j),1:3)=vec(1:3)
               exit checkspec
            end if
         end do checkspec
      end if
   end do


!!!-----------------------------------------------------------------------------
!!! move basis from temporary basis to main basis.
!!! done to allow for correct allocation of number of and per species
!!!-----------------------------------------------------------------------------
   allocate(basis%spec(basis%nspec))
   basis%spec(1:basis%nspec)%name=tmp_spec(1:basis%nspec)
   do i=1,basis%nspec
      basis%spec(i)%name=tmp_spec(i)
      basis%spec(i)%num=tmp_num(i)
      allocate(basis%spec(i)%atom(tmp_num(i),dim))
      basis%spec(i)%atom(:,:)=0
      basis%spec(i)%atom(1:tmp_num(i),1:3)=tmp_bas(i,1:tmp_num(i),1:3)
   end do


 end subroutine XYZ_geom_read
!!!#############################################################################


!!!#############################################################################
!!! generates cartesian basis
!!!#############################################################################
 subroutine XYZ_geom_write(UNIT,basis)
   implicit none
   integer, intent(in) :: UNIT
   type(bas_type), intent(in) :: basis

   integer :: i,j


   write(UNIT,'("I0")') basis%natom
   write(UNIT,'("A")') basis%sysname
   do i=1,basis%nspec
      do j=1,basis%spec(i)%num
         write(UNIT,'(A5,1X,3(F15.9))') &
              basis%spec(i)%name,basis%spec(i)%atom(j,1:3)
      end do
   end do


 end subroutine XYZ_geom_write
!!!#############################################################################


!!!#############################################################################
!!! reads lattice and basis from an extended xyz file
!!!#############################################################################
 subroutine extXYZ_geom_read(UNIT,basis,length)
   implicit none
   integer, intent(in) :: UNIT
   type(bas_type), intent(out) :: basis
   integer, intent(in), optional :: length

   integer :: Reason
   integer :: index1, index2
   integer, allocatable, dimension(:) :: tmp_num
   real(real12), dimension(3) :: vec
   real(real12), allocatable, dimension(:,:,:) :: tmp_bas
   character(len=5) :: ctmp
   character(len=5), allocatable, dimension(:) :: tmp_spec
   character(len=1024) :: buffer
   integer :: i,j,dim

   dim=3
   basis%lcart=.true.
   if(present(length)) dim=length


!!!-----------------------------------------------------------------------------
!!! read system information
!!!-----------------------------------------------------------------------------
   read(UNIT,*,iostat=Reason) basis%natom
   if(Reason.ne.0)then
      write(0,'(" The file is not in xyz format.")')
      write(0,'(" Exiting code ...")')
      call exit()
   end if
   read(UNIT,'(A)',iostat=Reason) buffer
   if(Reason.ne.0)then
      write(0,'(" The file is not in xyz format.")')
      write(0,'(" Exiting code ...")')
      call exit()
   end if
   index1 = index(buffer,'Lattice="') + 9
   index2 = index(buffer(index1:),'"') + index1 - 2
   read(buffer(index1:index2),*) ( ( basis%lat(i,j), j = 1, 3), i = 1, 3)

   index1 = index(buffer,'free_energy=') + 12
   read(buffer(index1:),*) basis%energy


!!!-----------------------------------------------------------------------------
!!! read basis
!!!-----------------------------------------------------------------------------
   allocate(tmp_spec(basis%natom))
   allocate(tmp_num(basis%natom))
   allocate(tmp_bas(basis%natom,basis%natom,dim))
   tmp_num(:)=0
   tmp_spec=""
   tmp_bas=0
   basis%nspec=0
   do i=1,basis%natom
      read(UNIT,*,iostat=Reason) ctmp,vec(1:3)
      if(.not.any(tmp_spec(1:basis%nspec).eq.ctmp))then
         basis%nspec=basis%nspec+1
         tmp_spec(basis%nspec)=ctmp
         tmp_bas(basis%nspec,1,1:3)=vec(1:3)
         tmp_num(basis%nspec)=1
      else
         checkspec: do j=1,basis%nspec
            if(tmp_spec(j).eq.ctmp)then
               tmp_num(j)=tmp_num(j)+1
               tmp_bas(j,tmp_num(j),1:3)=vec(1:3)
               exit checkspec
            end if
         end do checkspec
      end if
   end do


!!!-----------------------------------------------------------------------------
!!! move basis from temporary basis to main basis.
!!! done to allow for correct allocation of number of and per species
!!!-----------------------------------------------------------------------------
   allocate(basis%spec(basis%nspec))
   basis%spec(1:basis%nspec)%name=tmp_spec(1:basis%nspec)
   basis%sysname = ""
   do i=1,basis%nspec
      basis%spec(i)%name=tmp_spec(i)
      basis%spec(i)%num=tmp_num(i)
      allocate(basis%spec(i)%atom(tmp_num(i),dim))
      basis%spec(i)%atom(:,:)=0
      basis%spec(i)%atom(1:tmp_num(i),1:3)=tmp_bas(i,1:tmp_num(i),1:3)
      write(buffer,'(I0,A)') basis%spec(i)%num,trim(basis%spec(i)%name)
      basis%sysname = basis%sysname//trim(buffer)
      if(i.lt.basis%nspec) basis%sysname = basis%sysname//"_"
   end do

 end subroutine extXYZ_geom_read
!!!#############################################################################


!!!#############################################################################
!!! generates cartesian basis
!!!#############################################################################
 subroutine extXYZ_geom_write(UNIT,basis)
   implicit none
   integer, intent(in) :: UNIT
   type(bas_type), intent(in) :: basis

   integer :: i,j


   write(UNIT,'("I0")') basis%natom
   write(UNIT,'(A,8(F0.8,1X),F0.8,A)', advance="no") &
        'Lattice="',((basis%lat(i,j),j=1,3),i=1,3),'"'
   write(UNIT,'(A,F0.8)', advance="no") ' free_energy=',basis%energy
   write(UNIT,'(A)', advance="no") ' pbc="T T T"'
   if(basis%lcart)then
      do i=1,basis%nspec
         do j=1,basis%spec(i)%num
            write(UNIT,'(A8,3(1X, F16.8))') &
               basis%spec(i)%name,basis%spec(i)%atom(j,1:3)
         end do
      end do
   else
      do i=1,basis%nspec
         do j=1,basis%spec(i)%num
            write(UNIT,'(A8,3(1X, F16.8))') basis%spec(i)%name, &
               matmul(basis%spec(i)%atom(j,1:3),basis%lat)
         end do
      end do
   end if

 end subroutine extXYZ_geom_write
!!!#############################################################################


!!!#############################################################################
!!! convert basis using latconv transformation matrix
!!!#############################################################################
 function convert_bas(basis,latconv) result(outbas)
   implicit none
   type(bas_type), intent(in) :: basis
   real(real12), dimension(3,3), intent(in) :: latconv
   
   integer :: is,ia,dim
   type(bas_type) :: outbas
   

   dim=size(basis%spec(1)%atom(1,:))
   allocate(outbas%spec(basis%nspec))
   outbas%natom=basis%natom
   outbas%nspec=basis%nspec
   outbas%sysname=basis%sysname
   outbas%lcart=.not.basis%lcart
   do is=1,basis%nspec
      allocate(outbas%spec(is)%atom(basis%spec(is)%num,dim))
      outbas%spec(is)=basis%spec(is)
      do ia=1,basis%spec(is)%num
         outbas%spec(is)%atom(ia,1:3)=&
              matmul(latconv,outbas%spec(is)%atom(ia,1:3))
      end do
   end do
   
 end function convert_bas
!!!#############################################################################


!!!#############################################################################
!!! converts lattice from abc and αβγ to lattice matrix
!!!#############################################################################
 function convert_abc_to_lat(abc,angle,radians) result(out_lat)
   use constants, only: pi
   implicit none
   real(real12), dimension(3) :: in_angle
   real(real12), dimension(3,3) :: out_lat

   real(real12), dimension(3), intent(in) :: abc,angle

   logical, optional, intent(in) :: radians


   if(present(radians))then
      if(.not.radians) in_angle=angle*pi/180._real12
   end if
!      in_angle=angle*pi/180._real12 ! this looks wrong, check it

   out_lat=0._real12

   out_lat(1,1)=abc(1)
   out_lat(2,:2)=(/abc(2)*cos(in_angle(3)),abc(2)*sin(in_angle(3))/)

   out_lat(3,1) = abc(3)*cos(in_angle(2))
   out_lat(3,2) = abc(3)*(cos(in_angle(1)) - cos(in_angle(2))*&
        cos(in_angle(3)))/sin(in_angle(3))
   out_lat(3,3) = sqrt(abc(3)**2._real12 - &
        out_lat(3,1)**2._real12 - &
        out_lat(3,2)**2._real12)


 end function convert_abc_to_lat
!!!#############################################################################


!!!#############################################################################
!!! converts lattice from matrix to abc and αβγ
!!!#############################################################################
 function convert_lat_to_abc(in_lat,radians) result(abc_angle)
   use constants, only: pi
   implicit none
   integer :: i
   real(real12), dimension(2,3) :: abc_angle

   real(real12), dimension(3,3), intent(in) :: in_lat

   logical, optional, intent(in) :: radians


   do i=1,3
      abc_angle(1,i)=modu(in_lat(i,:))
   end do
   do i=1,3
   end do
   abc_angle(2,1)=acos(dot_product(in_lat(2,:),in_lat(3,:))/&
        (abc_angle(1,2)*abc_angle(1,3)))
   abc_angle(2,3)=acos(dot_product(in_lat(1,:),in_lat(3,:))/&
        (abc_angle(1,1)*abc_angle(1,3)))
   abc_angle(2,3)=acos(dot_product(in_lat(1,:),in_lat(2,:))/&
        (abc_angle(1,1)*abc_angle(1,2)))

   if(present(radians))then
      if(.not.radians) abc_angle(2,:)=abc_angle(2,:)*180._real12/pi
   end if

 end function convert_lat_to_abc
!!!#############################################################################


!!!#############################################################################
!!! clones basis 1 onto basis 2
!!!#############################################################################
 subroutine clone_bas(inbas,outbas,trans_dim)
   implicit none
   integer :: i
   integer :: indim,outdim
   real(real12) :: val
   logical :: udef_trans_dim

   type(bas_type) :: inbas,outbas

   logical, optional, intent(in) :: trans_dim


!!!-----------------------------------------------------------------------------
!!! determines whether user wants output basis extra translational dimension
!!!-----------------------------------------------------------------------------
   indim = size(inbas%spec(1)%atom(1,:),dim=1)
   if(present(trans_dim))then
      udef_trans_dim = trans_dim
   elseif(indim.eq.4)then
      udef_trans_dim = .true.
   elseif(indim.eq.3)then
      udef_trans_dim = .false.
   end if


!!!-----------------------------------------------------------------------------
!!! sets up output basis atomic coordinates dimension
!!!-----------------------------------------------------------------------------
   if(udef_trans_dim)then
      outdim = 4
      val = 1._real12
   else
      outdim = 3
      val = 0._real12
   end if


!!!-----------------------------------------------------------------------------
!!! if already allocated, deallocates output basis
!!!-----------------------------------------------------------------------------
   if(allocated(outbas%spec))then
      do i=1,outbas%nspec
         if(allocated(outbas%spec(i)%atom)) deallocate(outbas%spec(i)%atom)
      end do
      deallocate(outbas%spec)
   end if


!!!-----------------------------------------------------------------------------
!!! allocates output basis and clones data from input basis to output basis
!!!-----------------------------------------------------------------------------
   allocate(outbas%spec(inbas%nspec))
   do i=1,inbas%nspec
      allocate(outbas%spec(i)%atom(&
           inbas%spec(i)%num,outdim))
      if(indim.eq.outdim)then
         outbas%spec(i)%atom(:,:indim) = inbas%spec(i)%atom(:,:indim)
      elseif(outdim.gt.indim)then
         outbas%spec(i)%atom(:,:indim) = inbas%spec(i)%atom(:,:indim)
         outbas%spec(i)%atom(:,outdim) = val
      else
         outbas%spec(i)%atom(:,:outdim) = inbas%spec(i)%atom(:,:outdim)
      end if
      outbas%spec(i)%mass = inbas%spec(i)%mass
      outbas%spec(i)%num = inbas%spec(i)%num
      outbas%spec(i)%name = inbas%spec(i)%name
   end do
!   outbas = inbas !using this will reallocate outbas to inbas
   outbas%nspec = inbas%nspec
   outbas%natom = inbas%natom
   outbas%lcart = inbas%lcart
   outbas%sysname = inbas%sysname
   outbas%energy = inbas%energy
   outbas%lat = inbas%lat


   return
 end subroutine clone_bas
!!!#############################################################################

end module rw_geom
