module vasp_file_handler
  use constants, only: real12
  use atomtype
  implicit none

  type unitcell 
     real(real12), dimension(3,3) :: cell 
  end type unitcell

contains

!!!#############################################################################
!!! repeat the atoms in the unit cell out to a -1,0,1 grid
!!!#############################################################################
  subroutine atomrepeater(structures,position,element,array,unit,atomnumber,length)
    use atomtype
    implicit none
    type (atom), dimension(:,:), allocatable :: array
    real(real12), dimension(3) :: position
    type(unitcell), dimension(:), allocatable :: unit

    integer :: j,atomnumber,length,x,y,z,structures,m
    character(3) :: element
    m=((atomnumber-1)*27)
    do x=-1,1
       do y=-1,1
          do z=-1,1
             m = m + 1
             do j=1, 3
                array(structures,m)%position(j)=position(j) + &
                     &(x*unit(structures)%cell(1,j)) + &
                     &(y*unit(structures)%cell(2,j)) + &
                     &(z*unit(structures)%cell(3,j))
             end do
             array(structures,m)%name=element
          end do
       end do
    end do
  end subroutine atomrepeater
!!!#############################################################################


!!!#############################################################################
!!! write INCAR file
!!!#############################################################################
  subroutine Incarwrite(filepath,nstep,bandno)
    implicit none
    character(1024) ::  name
    integer :: nstep, bandno
    character(*) :: filepath

    integer :: unit
    
    write(name,'(A,A6)') trim(filepath), "/INCAR"
    write(*,*) "The name of the new file will be", trim(adjustl(name))
    
    !write(*,*) trim(adjustl(name))
    open(newunit=unit,file=trim(adjustl(name)))!"pos/POSCAR_001")!trim(name))
    
    write(unit, *)"SYSTEM = RSS"
    write(unit, *)"### sys ###"
    write(unit, *)"GGA = PE "
    write(unit, *)"ISYM = 0"
    write(unit, *)"ENCUT = 500"
    write(unit, *)"ENAUG = 500"
    write(unit, *)"PREC = Accurate"
    write(unit, *)"ISTART = 1"
    write(unit, *)"ICHARG = 2"
    write(unit, *)"LWAVE = .FALSE."
    write(unit, *)"LAECHG =.FALSE."
    write(unit, *)"LMAXMIX = 6"
    write(unit, *)"LASPH = .TRUE."
    write(unit, *)"LMIXTAU = .TRUE."
    write(unit, *)"LREAL = .FALSE."
    write(unit, *)"NWRITE = 3"
    write(unit, *)"ALGO = N"
    write(unit, *)"LHAR = .FALSE."
    write(unit, *)"LVTOT = .FALSE."
    write(unit, *)
    
    write(unit, *)"### vdw ###"
    write(unit, *)"#IVDW = 11"
    write(unit, *)"#VDW_RADIUS = 80"
    write(unit, *)"#VDW_CNRADIUS = 50.0"
    write(unit, *)"#VDW_S8 = 0.722"
    write(unit, *)"#VDW_SR = 1.217"
    write(unit, *) 
    
    write(unit, *)"### elc ###"
    write(unit, *)"#AMIX = 0.6"
    write(name,'(I0.3)') nstep
    write(unit, *)"NELM = ", adjustL(trim(name))
    write(unit, *)"NELMIN = 5"
    write(unit, *)"NELMDL = -5"
    write(unit, '(1X,A,1X,I0.3)') "NBANDS =", bandno
    write(unit, *)"EDIFF = 10d-8"
    write(unit, *)"ISMEAR = 0"
    
    !!! THIS SHOULD BE CHANGED BASED ON INTUITION
    
    write(unit, *)"SIGMA = 0.2"
    write(unit, *) 
    
    write(unit, *)"### Mag ###"
    write(unit, *)"#ISPIN = 2"
    write(unit, *)"#MAGMOM = 1 -1 1 1 -1 -1"
    write(unit, *) 
    
    write(unit, *)"### ncl ###"
    write(unit, *)"#LNONCOLLINEAR = .TRUE."
    write(unit, *)"#LSORBIT = .TRUE."
    write(unit, *)"#GGA_COMPAT = .FALSE."
    write(unit, *) 
    
    write(unit, *)"### MPI ###"
    write(unit, *)"#NCORE = 24"
    write(unit, *)"#KPAR = 2"
    write(unit, *) 
    
    write(unit, *)"### Rlx ###"
    write(unit, *)"#LMAXPAW =-1"
    write(unit, *)"ADDGRID = .TRUE."
    write(unit, *)"#POTIM = 0.1"
    write(unit, *)"#NFREE = 15"
    write(unit, *)"#NSW = 150"
    write(unit, *)"#ISIF = 2"
    write(unit, *)"#IBRION = 1"
    write(unit, *)"#EDIFFG = -0.001"
    write(unit, *) 
    
    write(unit, *)"### dos ###"
    write(unit, *)"#EMIN = -3.0"
    write(unit, *)"#EMAX =  6.0"
    write(unit, *)"#NEDOS = 5000"
    write(unit, *)"#LORBIT = 11"
    
    close(unit)
  end subroutine Incarwrite
!!!#############################################################################


!!!#############################################################################
!!! write POTCAR file
!!!#############################################################################
  subroutine potwrite(filepath, elnames, eltot) 

    character(3), dimension(:), allocatable :: elnames
    character(1024) ::  name
    character(20) :: filepath
    integer :: eltot, i
    name=" "
    !write(*,*) filepath
    !write(*,*) eltot
    do i=1, eltot 
    
       !write(*,*) i
     if(i.eq.1) then
        !write(*,*) "I reached here"
        write(name, '(A3,1X,A6,A2,1X,A1,A,A)') "cat", "potcar", elnames(i),">",trim(adjustl(filepath)),"/POTCAR"
        !write(*,*) name
        call execute_command_line(name)
     else       
        write(name, '(A,1X,A,A,1X,A,A,1X,A,1X,A,A)') "cat ",trim(adjustl(filepath)),"/POTCAR", "potcar", elnames(i)," >>",&
             & trim(adjustl(filepath)),"/POTCAR1"   
    
        call execute_command_line(name)
        write(name,'(A,1X,A,A8,1X,A,A)') "mv ",trim(adjustl(filepath)),"/POTCAR1", trim(adjustl(filepath)),"/POTCAR"
        call execute_command_line(name)
     end if
    
    end do
    end subroutine potwrite
!!!#############################################################################


!!!#############################################################################
!!! write job and KPOINTS files
!!!#############################################################################
  subroutine Jobwrite(filepath, a,b,c)
    implicit none
    integer :: a,b,c
    character(1024) ::  name
    character(20) :: filepath
    integer :: unit

    name=" "
    !write(*,*) filepath
    write(name,*) "cp job_vasp_isca.in ", trim(filepath)
    call execute_command_line(name)
    write(name,'(A,A8)') trim(filepath), "/KPOINTS"
    open(newunit=unit,file=trim(adjustl(name)), status='new')!"pos/POSCAR_001")!trim(name))
    write(unit, '(A7)') "KPOINTS"
    write(unit, '(A1)') "0" 
    write(unit, '(A1)') "G" 
    write(unit, '(I0,X,I0,X,I0)') a, b, c
    write(unit, '(I0,X,I0,X,I0)') 0,0,0
    close(unit) 
  end subroutine Jobwrite
!!!#############################################################################


!!!#############################################################################
!!! write POSCAR file
!!!#############################################################################
  subroutine  poswrite(box,atomlist,len, structures, structno, prev_structures)
  
    integer :: i,j, len, reason, structures, structno, prev_structures
    type(atom), dimension(:,:) :: atomlist
    real(real12), dimension(3,3) :: box
    do i=1, len 
     if(i.eq.len) then 
        write(structures+10000,'(5X,A2)', IOStat=reason, advance="no")&
             &atomlist(structures,i)%name
     else 
        !if(i.eq.(len-1)) cycle
        if (atomlist(structures,i)%name.eq.atomlist(structures,i+1)%name) cycle
        write(structures+10000,'(5X,A2)', IOStat=reason, advance="no") &
             &atomlist(structures,i)%name      
     end if
    end do
    
    write(structures+10000, '(2X)')
    !write(structures+10000,*) "HERE"
    j=1
    do i=1, len 
     if(i.eq.len) then 
        write(structures+10000,'(1X,I3)', IOStat=reason, advance="no") j
     else
        if (atomlist(structures,i)%name.eq.atomlist(structures,i+1)%name) then 
           j=j+1
           cycle
        else 
           write(structures+10000,'(1X,I3)', IOStat=reason, advance="no") j
           j=1
    
        end if
     end if
    end do
    
    
    write(structures+10000,'(/,A)') "Cartesian"
    
    do i=1, len 
     write(structures+10000,'(1X,F0.16,1X,F0.16,1X,F0.16)') atomlist(structures,i)%position(:)
    end do
    
  end subroutine poswrite
!!!#############################################################################


!!!#############################################################################
!!! touch POSCAR subdirectories within pos/ directory
!!!#############################################################################
  subroutine touchposdir(structures,prev_structures)
    character(1024) :: name, tmp, command  
    logical :: dir_e
    integer :: structures, prev_structures 
    write(name,'(A11,I0.3,A7)')"pos/POSCAR_",structures+prev_structures,"/POSCAR"
    !!Calculates the new structure number, and writes it to tmp
    write(tmp,'(A11,I0.3)')"pos/POSCAR_",structures+prev_structures
    !!Checks if a directory to contain that file exsts already (It should never exist, could add warning)
    inquire(file=tmp, exist=dir_e)
    if(dir_e) then
    else
       !!Writes a command to create said directory
     write(command,'(A17,I0.3)')"mkdir pos/POSCAR_",structures+prev_structures
     call execute_command_line(trim(adjustl(command)))
    end if
  end subroutine touchposdir
!!!#############################################################################

!!!#############################################################################
!!! touch pos/ directory
!!!#############################################################################
  subroutine touchpos() 
    !> Creates a directory named "pos" if it does not already exist.
    !> This subroutine is used to ensure that the directory "pos" exists before performing any operations on it.
    !> If the directory already exists, no action is taken.
    !> If the directory does not exist, it is created using the "mkdir" command.
    !> This subroutine does not have any input or output parameters.
    character(1024) :: file, command 
    logical :: dir_e
  
    inquire(file="pos", exist=dir_e)
    if(dir_e) then
    else
       write(command,*) "mkdir pos"
       call execute_command_line(command)
    end if
  end subroutine touchpos
!!!#############################################################################


!!!#############################################################################
!!! Add POSCAR (directory or file????)
!!!#############################################################################
  subroutine addposcar(option_generate_files_only,option_filepath,stage, prev_structures_overwrite)
    type(unitcell), dimension(:), allocatable :: formula
    character(1024) :: tmp, name,option_filepath
    character(3), dimension(:), allocatable :: elnames
    real(real12) :: cellmultiplier
    real(real12), dimension(:), allocatable :: bondcutoff
    integer :: q,l,k,j,i,structno, ecount, eltot, leng,structures, stage
    integer :: prev_structures, option_generate_files_only, prev_structures_overwrite
    type (atom), dimension(:,:), allocatable :: tmplist,atomlist, alistrep
    integer, dimension(:), allocatable :: elno
    !! Wipes the randomly generated formula
    structures=1
    structno=1
    
    if (option_generate_files_only.eq.0) then
     prev_structures=structurecounter("pos")
     call touchpos()
     call touchposdir(structures,prev_structures)
    end if
    
    write(*,*) structures, prev_structures, "!!"
    allocate(formula(structno))
    
    write(6,*) "Please enter the filename you wish to add to the database"
    !read(*, *) name
    if (option_generate_files_only.eq.0) then
     write(name,'(A)') "POSCAR"
    else 
     write(*,*) "here"
     write(name,'(A,A)') trim(adjustl(option_filepath)),"/POSCAR"
    end if
    write(*,*) name
    open(50, file=name)
    write(*,*) "step 1"
    read(50, '(A)') tmp
    read(50, '(F16.0)') cellmultiplier
    
    
    do i=1, 3
     read(50,*) formula(1)%cell(i,1),formula(1)%cell(i,2),formula(1)%cell(i,3)
     formula(1)%cell(i,1)=formula(1)%cell(i,1)*cellmultiplier
     formula(1)%cell(i,2)=formula(1)%cell(i,2)*cellmultiplier
     formula(1)%cell(i,3)=formula(1)%cell(i,3)*cellmultiplier
     write(*,*) formula(1)%cell(i,1),formula(1)%cell(i,2),formula(1)%cell(i,3)
    
    end do
    
    
    
    
    if (option_generate_files_only.eq.0) then
    
     write(name,'(A11,I0.3,A7)')"pos/POSCAR_",structures+prev_structures,"/POSCAR"
     open(structures+10000, file=name,status="new")
     write(structures+10000,*) "Test"
     write(structures+10000,*) 1.0
     do i=1, 3
        write(structures+10000,*) formula(1)%cell(i,1), &
             &formula(1)%cell(i,2), formula(1)%cell(i,3)
     end do
    end if
    read(50,'(A)') tmp
    
    
    ecount=0
    eltot=0
    allocate(tmplist(structno,1000))
    allocate(elnames(1024))
    do i=1,len(tmp)-1
     if(i.eq.1) then 
        if((scan(tmp(i:i+1)," ").eq.0).or.&
             &((scan(tmp(i:i+1)," ").eq.2))) then
           eltot=eltot+1
           !write(*,*) tmp(i:i+1)
           if(scan(tmp(i:i+1)," ").eq.0) then
              elnames(eltot)=tmp(i:i+1)
           end if
           if((scan(tmp(i:i+1)," ").eq.2).and.(scan(tmp(i-1:i)," ").eq.1)) then
              elnames(eltot)=tmp(i:i)
           end if
        else
        end if
    
    
     else
        if((scan(tmp(i:i+1)," ").eq.0).or.&
             &((scan(tmp(i:i+1)," ").eq.2).and.(scan(tmp(i-1:i)," ")&
             &.eq.1))) then
           eltot=eltot+1
           !write(*,*) tmp(i:i+1)
           if(scan(tmp(i:i+1)," ").eq.0) then
              elnames(eltot)=tmp(i:i+1)
           end if
           if((scan(tmp(i:i+1)," ").eq.2).and.(scan(tmp(i-1:i)," ").eq.1)) then
              elnames(eltot)=tmp(i:i)
           end if
        else
        end if
     end if
    end do
    write(*,*) elnames(:), eltot
    allocate(bondcutoff(eltot))
    allocate(elno(eltot))
    !read(50,'(4X)', advance='no')
    k=0
    read(50,*) elno
    write(*,*) elno(:)
    read(50, *) tmp
    k=0
    do i=1, eltot
     do j=1, elno(i)
        k=k+1
    
        tmplist(structures,k)%name=elnames(i)
     end do
    end do
    write(*,*) k
    l=0
    do i=1, eltot
     do j=1, elno(i)
        l=l+1
        read(50, *)&
             &tmplist(structures,l)%position(1),tmplist(structures,l)%position(2),&
             &tmplist(structures,l)%position(3)
        write(*,*) "Warning! Direct POSCAR import ONLY"
        write(*,*) tmplist(structures,l)%position(1),tmplist(structures,l)%position(2),&
             &tmplist(structures,l)%position(3)
        wait(5)
        do q=1,3
           tmplist(structures,L)%position(q)=&
                &formula(structures)%cell(1,q)*tmplist(structures,L)%position(1)+&
                &formula(structures)%cell(2,q)*tmplist(structures,L)%position(2)+&
                &formula(structures)%cell(3,q)*tmplist(structures,L)%position(3)
           write(*,*) tmplist(structures,L)%position(q)
        end do
        !write(*,*) atomlist(structures,j)%position(:)
        !write(*,*) atomlist(1,j)%position(1),atomlist(1,j)%position(2),atomlist(1,j)%position(3\
    
     end do
    end do
    allocate(atomlist(1,k))
    do i=1, k
     atomlist(structures,i)%position=tmplist(structures,i)%position
     atomlist(structures,i)%name=tmplist(structures,i)%name
     !write(*,*) tmplist(structures,i)%position(:)
    
    end do
    deallocate(tmplist)
    if (option_generate_files_only.eq.0) then
     call poswrite(formula(structures)%cell,atomlist,k,1,1,prev_structures)
    end if
    allocate(alistrep(1,k*27))
    do i=1, k
     call atomrepeater(structures,atomlist(structures,i)%position,atomlist(structures,i)%name,&
          &alistrep,&
          &formula,i,k)
    end do
    !!! The problem is here. Firstly, you need to pass STAGE to this function so it can be correctly passed to the following functions 
    !!! Secondly, prev_structures shouldn't be calculated but also passed from the previous function, as loop number from calling loop
    !!!
    if (prev_structures_overwrite.eq.0) then 
     prev_structures=structurecounter("bon")
    else  
     prev_structures=prev_structures_overwrite-1
    end if
    
    !do i=1, k
    !   call generatebondfiles(stage,structures,prev_structures,atomlist,alistrep,eltot,elno,i,elnames)
    !end do
    bondcutoff(1)=1.0
    write(*,*) "MANUALLY EDIT"
    stop
    if (prev_structures_overwrite.eq.0) prev_structures=structurecounter("bad")
    !do i=1, k
    !   call generateanglefiles(stage,structures,prev_structures,atomlist,alistrep,eltot,elno,i,bondcutoff)
    !end do
    
    !do i=1, k
    !   call generate4files(stage,structures,prev_structures,atomlist,alistrep,eltot,elno,i,bondcutoff)
    !end do
    
    !do i=1, k 
    !   call generate4files(stage,structures,prev_structures,atomlist,alistrep,eltot,elno,i,bondcutoff)
    !end do
    
    write(tmp,'(A11,I0.3)')"pos/POSCAR_",structures+prev_structures
    call Incarwrite(trim(adjustl(tmp)),500, 20*k)
    call Jobwrite(tmp,1,1,1)
    call potwrite(tmp, elnames, eltot)
    
    
    
    
  end subroutine addposcar  
!!!#############################################################################


!!!#############################################################################
!!! read POSCAR file
!!!#############################################################################
  subroutine poscar_read(pathway,unit_cell,structure_elements,structure_stochiometry, structure_factor,&
    &coordinate_type,atomic_positions, header)
   use atomtype 
   implicit none 
   character(1024), intent(out) :: coordinate_type, header 
   character(1024) :: read_in, tmp
   character(1024), intent(in) :: pathway
   character(3), dimension(:), allocatable, intent(out) :: structure_elements
   character(3), dimension(:), allocatable :: temp_names
   integer, dimension(:), allocatable, intent(out) :: structure_stochiometry 
   integer :: loop, status, size,eltot, i, j,k,l,m
   integer, parameter :: line_buf_len= 1024*4
   character(LEN=line_buf_len) :: InS
   type(unitcell), dimension(1), intent(out) :: unit_cell 
   real(real12), dimension(:,:), intent(out), allocatable :: atomic_positions
   real(real12), dimension(:,:), allocatable :: temp_positions
   real(real12), intent(out) :: structure_factor 
   logical :: file_e  
   
   eltot=0
   allocate(structure_elements(1))
   inquire(file=trim(adjustl(pathway)), exist=file_e)
   write(*,*) file_e
   if(file_e) then 
    open(102,file=pathway)
   
    read(102,*) header
    read(102,*) structure_factor
    do loop=1, 3
       read(102,*) unit_cell(1)%cell(loop,1), unit_cell(1)%cell(loop,2), unit_cell(1)%cell(loop,3)
    end do
    read(102,'(A)') tmp
    write(*,*) tmp
    write(read_in,'(X,A)') trim(adjustl(tmp))
    write(*,*) read_in
    do i=1,len(read_in)-1
       if(i.eq.1) then
          if((scan(read_in(i:i+1)," ").eq.0).or.&
               &((scan(read_in(i:i+1)," ").eq.2))) then
             eltot=eltot+1
             if(scan(read_in(i:i+1)," ").eq.0) then
                allocate(temp_names(eltot))
   
                do j=1, eltot-1
                   temp_names(j)=structure_elements(j)
                end do
                temp_names(eltot)=read_in(i:i+1)
                deallocate(structure_elements)
                allocate(structure_elements(eltot))
                structure_elements=temp_names
                deallocate(temp_names)
   
             end if
             if((scan(read_in(i:i+1)," ").eq.2).and.(scan(read_in(i-1:i)," ").eq.1)) then
                allocate(temp_names(eltot))
   
                do j=1, eltot-1
                   temp_names(j)=structure_elements(j)
                end do
                temp_names(eltot)=read_in(i:i)
                deallocate(structure_elements)
                allocate(structure_elements(eltot))
                structure_elements=temp_names
                deallocate(temp_names)
             end if
          else
          end if
       else
          if((scan(read_in(i:i+1)," ").eq.0).or.&
               &((scan(read_in(i:i+1)," ").eq.2).and.(scan(read_in(i-1:i)," ")&
               &.eq.1))) then
             eltot=eltot+1
             !write(*,*) tmp(i:i+1)
             if(scan(read_in(i:i+1)," ").eq.0) then
                allocate(temp_names(eltot))
   
                do j=1, eltot-1
                   temp_names(j)=structure_elements(j)
                end do
                temp_names(eltot)=read_in(i:i+1)
                deallocate(structure_elements)
                allocate(structure_elements(eltot))
                structure_elements=temp_names
                deallocate(temp_names)
   
   
             end if
             if((scan(read_in(i:i+1)," ").eq.2).and.(scan(read_in(i-1:i)," ").eq.1)) then
                allocate(temp_names(eltot))
   
                do j=1, eltot-1
                   temp_names(j)=structure_elements(j)
                end do
                temp_names(eltot)=read_in(i:i)
                deallocate(structure_elements)
                allocate(structure_elements(eltot))
                structure_elements=temp_names
                deallocate(temp_names)
   
             end if
          else
          end if
       end if
    end do
    allocate(structure_stochiometry(eltot))
   
    read(102,*) structure_stochiometry 
    read(102,*) coordinate_type
    allocate(temp_positions(sum(structure_stochiometry),3))
    allocate(atomic_positions(sum(structure_stochiometry),3))
    if(coordinate_type.eq."Cartesian") then 
       do i=1, sum(structure_stochiometry)
          read(102,*) atomic_positions(i,1),atomic_positions(i,2),atomic_positions(i,3)
       end do
    else
       write(*,*) coordinate_type
       do i=1, sum(structure_stochiometry)
   
          read(102,*) temp_positions(i,1),temp_positions(i,2),temp_positions(i,3)
   
       end do
       do i=1, sum(structure_stochiometry)
          atomic_positions(i,:)=matmul(temp_positions(i,:),unit_cell(1)%cell) 
          !temp_positions(i,1)*unit_cell(1,1)+temp_positions(i,2)*unit_cell() 
          !write(*,*) atomic_positions(i,:)
          !write(*,*) temp_positions(i,:) 
          !write(*,*) unit_cell(i)%cell 
          !call sleep(5) 
       end do
    end if
   else
    write(*,*) "ERROR: NOT FOUND"
   end if
   
   end subroutine poscar_read
!!!#############################################################################


!!!#############################################################################
!!! count number of structures in a directory
!!!#############################################################################
  function structurecounter(dir_t)  result(n)
    logical :: file_exists
    character(1024) :: name, dir_name, dir_append
    character(3) :: dir_t
    integer :: n, gap_test
    n=1
    mainloop : do while(n.ne.0)
     if(dir_t.eq."don") dir_name="/DON_"     
     if(dir_t.eq."pos") dir_name="/POSCAR_"
     if(dir_t.eq."don") dir_append="/DON"
     if(dir_t.eq."pos") dir_append="/POSCAR"
     if(dir_t.eq."bon") dir_name="/BON_"
     if(dir_t.eq."bon") dir_append=""
     if(dir_t.eq."bad") dir_name="/BAD_"
     if(dir_t.eq."bad") dir_append=""
     write(name,'(A3,A,I0.3,A7)') dir_t,trim(adjustl(dir_name)), n,dir_append
     name=trim(adjustl(name))
    
     INQUIRE(FILE=name, EXIST=file_exists)
     if(file_exists) then 
        n=n+1
        cycle 
     else 
        do gap_test=1, 1000
           write(name,'(A3,A,I0.3,A7)') dir_t,trim(adjustl(dir_name)), n+gap_test,dir_append
           INQUIRE(FILE=name, EXIST=file_exists)
           if(file_exists) then
              n=n+gap_test
              cycle mainloop 
           end if
        end do
        n=n-1
    
    
        exit
     end if
    end do mainloop
  end function structurecounter
!!!#############################################################################
   
end module vasp_file_handler