module file_generator
  use constants, only: real12, pi
  use geom, only: get_bondlength, get_bondangle, get_dihedral_angle
  use vasp_file_handler, only: unitcell, &
       poscar_read
  use geom
  use atomtype, only: atom
  implicit none

contains


!!!#############################################################################
!!! generate file associated with bondlengths
!!!#############################################################################
  subroutine generatebondfiles( &
       stage,structures,prev_structures,&
       array,repeatedarray,eltot,stochio,atomnumber,structure_elements)
    implicit none
    character(1024) :: command,name,tmp, location_string, debug_string
    integer :: l,i,el,structures,prev_structures,stage, eltot, atomnumber, tmpint,m,j,k
    type (atom), dimension(:,:), allocatable :: array, repeatedarray
    integer, dimension(:), allocatable :: stochio
    real(real12), dimension(:), allocatable :: list, list_tmp
    character(3), dimension(:), allocatable :: structure_elements
    
    integer :: unit
    logical :: dir_e


    if (stage.eq.1) location_string="bon"
    if (stage.eq.2) location_string="bon2"
    if (stage.eq.3) location_string="bon3"
    
    inquire(file=adjustl(trim(location_string)),exist=dir_e)
    if(.not.dir_e)then
       write(name,*) "mkdir ", trim(adjustl(location_string))
       call execute_command_line(name)
    end if

    write(tmp,'(A,A,I0.3)') trim(adjustl(location_string)),"/BON_",structures+prev_structures
    inquire(file=tmp, exist=dir_e)
    if(.not.dir_e)then
       write(command,'(A,A,A,I0.3)') "mkdir ",trim(adjustl(location_string)),"/BON_",structures+prev_structures
       call execute_command_line(command)
    end if
    
    !! get index of atomnumber in structure relating to how many previous ...
    !! ... atoms of the same species have already been documented
    !! i.e. in the parent loop from which this is called, all atoms are ...
    !! ... looped over, and this handles the incrementing when the species ...
    !! ... changes
    i=1
    tmpint=atomnumber
    do while (i.ge.1)
       if((tmpint-stochio(i)).gt.0) then
          if(array(structures,atomnumber)%name.ne.&
               array(structures,tmpint-stochio(i))%name) then 
             tmpint=tmpint-stochio(i)
          end if
          i=i+1
          if(i.ge.size(stochio)) i=0
       else
          i=0
       end if
    end do

    !! open file
    write(name, '(A,A,I0.3,A1,A,A,I0.3)') trim(adjustl(location_string)),"/BON_",structures+prev_structures,"/",&
       &trim(adjustl(array(structures,atomnumber)%name)),"_",tmpint
    open(newunit=unit, file=name)

    !! loop over elements/species in the structure and print associated bondlengths
    m=0
    k=0
    do el=1, eltot
       m=0
       allocate(list(27*stochio(el)))
      
       !! loop over atoms in the structure and get bondlengths
       do i=1, stochio(el)
          do L=1, 27
             k=k+1
             if((L.eq.14).and.(get_bondlength(repeatedarray(structures,k)%position,&
                  &array(structures,atomnumber)%position)).eq.0) cycle
             if(get_bondlength(repeatedarray(structures,k)%position,&
                  &array(structures,atomnumber)%position).gt.5) cycle
      
             m=m+1
             list(m)=get_bondlength(repeatedarray(structures,k)%position,&
                  &array(structures,atomnumber)%position)
          end do
       end do
       !! check for duplicates
       do i=1, m, 1
          do j=i+1, m, 1
             if(i.eq.j) cycle
             if(abs(list(i)-list(j)).lt.0.1) then
                list(j)=0
             end if
          end do
       end do

       !! write species header
       write(unit,*) structure_elements(el)
    
       !! write bondlengths to file (if they are not too close to zero)
       do i=1, m
          if(list(i).lt.0.001) cycle
          write(unit,*) list(i), i
       end do
       deallocate(list)
    end do
    
    !! close file
    close(unit)

  end subroutine generatebondfiles
!!!#############################################################################
  

!!!#############################################################################
!!! generate file associated with bondangles
!!!#############################################################################
  subroutine generateanglefiles( &
       stage,structures,prev_structures,&
       array,repeatedarray,eltot,stochio,&
       atomnumber,bondlengthcutoff,structure_elements)
    character(1024) :: command,name,tmp,location_string, name2
    integer :: stage,prev_structures,structures,eltot,atomnumber,&
         bondint_1,bondint_2,bondint_3
    type (atom), dimension(:,:), allocatable :: array, repeatedarray
    integer, dimension(:), allocatable :: stochio
    real(real12), dimension(:,:), allocatable :: bondlengthcutoff
    character(3), dimension(:), allocatable :: structure_elements
   
    integer :: i,j,m,p,ip,l,k,n,x,tmpint
    integer :: unit1, unit2
    logical :: dir_e
   
    if (stage.eq.1) location_string="bad"
    if (stage.eq.2) location_string="bad2"
    if (stage.eq.3) location_string="bad3"
    inquire(file=location_string,exist=dir_e)
    if(.not.dir_e)then
     write(name,*) "mkdir ", trim(adjustl(location_string))
     call execute_command_line(name)
    end if

    write(tmp,'(A,A,I0.3)') trim(adjustl(location_string)),"/BAD_",structures+prev_structures
    inquire(file=tmp, exist=dir_e)
    if(.not.dir_e) then
     write(command,'(A,A,A,I0.3)')"mkdir ", trim(adjustl(location_string)),"/BAD_",structures+prev_structures
     call execute_command_line(command)
    end if


    !! get index of atomnumber in structure relating to how many previous ...
    !! ... atoms of the same species have already been documented
    !! i.e. in the parent loop from which this is called, all atoms are ...
    !! ... looped over, and this handles the incrementing when the species ...
    !! ... changes
    i=1
    tmpint=atomnumber
    do while (i.ge.1)
     write(*,*) tmpint, stochio(i)
     if((tmpint-stochio(i)).gt.0) then
        if(array(structures,atomnumber)%name.ne.(array(structures,tmpint-stochio(i))%name)) then
           tmpint=tmpint-stochio(i)
        end if
        i=i+1
        if(i.gt.size(stochio)) i=0
     else
        i=0
     end if
    end do
    write(*,*) tmpint, trim(adjustl(array(structures,atomnumber)%name))

    !! open file
    write(name, '(A,A,I0.3,A1,A,A,I0.3)') trim(adjustl(location_string)),"/BAD_",structures+prev_structures,"/",&
       &trim(adjustl(array(structures,atomnumber)%name)),"_",tmpint
    open(newunit=unit1, file=name)
   
    !! UNKNOWN
    do ip=1, size(structure_elements)
       if(structure_elements(ip).ne.array(structures,atomnumber)%name) then
          cycle
       else
          bondint_1=ip
       end if
    end do
   
     !! loop over elements/species in the structure and print associated bondangles
     x=0
     do i=1, eltot
      do j=1, stochio(i)
         do p=1, 27
            do ip=1, size(structure_elements)
               if(trim(adjustl(structure_elements(ip))).ne.&
                    &trim(adjustl(repeatedarray(structures,x+1)%name))) then
                  cycle
               else
                  bondint_2=ip
               end if
            end do
     
            m=0
            x=x+1
            if(x.eq.atomnumber) cycle
            if(get_bondlength(array(structures,atomnumber)%position,&
                 &repeatedarray(structures,x)%position).gt.bondlengthcutoff(bondint_1,bondint_2)) then 
               cycle
            end if
            if(get_bondlength(array(structures,atomnumber)%position,&
                 &repeatedarray(structures,x)%position).lt.0.01) cycle
     
            do k=1, eltot
               do L=1, stochio(k)
                  do n=1, 27
                     m=m+1
                     !write(*,*) i,j,p,k,l,n
                     if(m.eq.x) cycle
                     if(m.eq.atomnumber) cycle
     
                     do ip=1, size(structure_elements)
                        if(trim(adjustl(structure_elements(ip))).ne.&
                             &trim(adjustl(repeatedarray(structures,m)%name))) then
                           cycle
                        else
                           bondint_3=ip
                        end if
                     end do
     
     
     
                     if(get_bondlength(array(structures,atomnumber)%position,&
                          &repeatedarray(structures,m)%position).lt.0.01) cycle
                     if(get_bondlength(array(structures,atomnumber)%position,&
                          &repeatedarray(structures,m)%position).gt.bondlengthcutoff(bondint_1,bondint_3)) then 
                        cycle
                     end if
     
                     if(get_bondangle(repeatedarray(structures,m)%position,&
                          &array(structures,atomnumber)%position,&
                          &repeatedarray(structures,x)%position).ge.pi/2) then 
     
                        write(unit1,*) get_bondangle(repeatedarray(structures,m)%position,&
                             &array(structures,atomnumber)%position,&
                             &repeatedarray(structures,x)%position)
     
                        write(name2, '(A,A,I0.3,A1,A,A,A,A,A,A,I0.3)') &
                             &trim(adjustl(location_string)),"/BAD_",structures+prev_structures,"/",&
                             &trim(adjustl(structure_elements(bondint_2))),"_",&
                             &trim(adjustl(array(structures,atomnumber)%name)),"_",&
                             &trim(adjustl(structure_elements(bondint_3)))
                        open(unit2,file=name2, position="append")
                        write(unit2,*) pi - get_bondangle(repeatedarray(structures,m)%position,&
                             &array(structures,atomnumber)%position,&
                             &repeatedarray(structures,x)%position)
     
                     else
     
                        write(unit1,*) get_bondangle(repeatedarray(structures,m)%position,&
                             &array(structures,atomnumber)%position,&
                             &repeatedarray(structures,x)%position)
     
                        write(name2, '(A,A,I0.3,A1,A,A,A,A,A,A,I0.3)') &
                             &trim(adjustl(location_string)),"/BAD_",structures+prev_structures,"/",&
                             &trim(adjustl(structure_elements(bondint_2))),"_",&
                             &trim(adjustl(array(structures,atomnumber)%name)),"_",&
                             &trim(adjustl(structure_elements(bondint_3)))
                        open(newunit=unit2,file=name2, position="append")
                        write(unit2,*) get_bondangle(repeatedarray(structures,m)%position,&
                             &array(structures,atomnumber)%position,&
                             &repeatedarray(structures,x)%position)
                        close(unit2)
                     end if
    
                  end do
     
               end do
            end do
         end do
      end do
     end do
     write(*,*) "ANGLEFINISHED"
     close(unit1)

   end subroutine generateanglefiles
!!!#############################################################################

!!!#############################################################################
!!!
!!!#############################################################################
   subroutine generate4files(stage,structures,prev_structures,array,repeatedarray,&
        &eltot,stochio,atomnumber,bondlengthcutoff,structure_elements)
     character(1024) :: command,name,tmp,location_string
     integer :: stage,prev_structures,structures,eltot,atomnumber,tmpint,l,k,n,x,i,j,m,p,ii,jj,kk,LL, dim, tool, counter
     integer :: ip, bondint_1, bondint_2, bondint_3, bondint_4
     type (atom), dimension(:,:), allocatable :: array, repeatedarray
     integer, dimension(:), allocatable :: stochio
     real(real12) :: res
     real(real12), dimension(:,:), allocatable :: bondlengthcutoff
     real(real12), dimension(3) :: ab,ac,ad, t
     logical :: dir_e
     real(real12), dimension(:), allocatable :: list, list_tmp
     character(3), dimension(:), allocatable :: structure_elements
   
     if (stage.eq.1) location_string="4body"
     if (stage.eq.2) location_string="4body2"
     if (stage.eq.3) location_string="4body3"
     inquire(file=location_string,exist=dir_e)
     if(dir_e) then
     else
      write(name,*) "mkdir ", trim(adjustl(location_string))
      call execute_command_line(name)
     end if
     write(tmp,'(A,A,I0.3)') trim(adjustl(location_string)),"/4BOD_",structures+prev_structures
     inquire(file=tmp, exist=dir_e)
     if(dir_e) then
     else
      write(command,'(A,A,A,I0.3)')"mkdir ", trim(adjustl(location_string)),"/4BOD_",structures+prev_structures
      call execute_command_line(command)
     end if
     i=1
     tmpint=atomnumber
     do while (i.ge.1)
        if((tmpint-stochio(i)).gt.0) then
           if(array(structures,atomnumber)%name.ne.(array(structures,tmpint-stochio(i))%name)) then
              tmpint=tmpint-stochio(i)
           end if
       
           i=i+1
           if(i.gt.size(stochio)) i=0
        else
           i=0
        end if
     end do
     
     write(name, '(A,A,I0.3,A1,A,A,I0.3)') trim(adjustl(location_string)),"/4BOD_",structures+prev_structures,"/",&
        &trim(adjustl(array(structures,atomnumber)%name)),"_",tmpint
     open(101, file=name)
     counter=0
     allocate(list(1))
     x=0
   
     do ip=1, size(structure_elements)
        if(trim(adjustl(structure_elements(ip))).ne.&
             &trim(adjustl(repeatedarray(structures,atomnumber)%name))) then
           cycle
        else
           bondint_1=ip
        end if
     end do
   

     do i=1, eltot
      do j=1, stochio(i)
         do p=1, 27
     
            x=x+1
     
     
            do ip=1, size(structure_elements)
               if(trim(adjustl(structure_elements(ip))).ne.&
                    &trim(adjustl(repeatedarray(structures,x)%name))) then
                  cycle
               else
                  bondint_2=ip
               end if
            end do
    
            if(x.eq.atomnumber) cycle
            if(get_bondlength(array(structures,atomnumber)%position,&
                 &repeatedarray(structures,x)%position).gt.bondlengthcutoff(bondint_1,bondint_2)) cycle
     
            if(get_bondlength(array(structures,atomnumber)%position,&
                 &repeatedarray(structures,x)%position).lt.0.1) cycle
            if((P.eq.14).and.(get_bondlength(repeatedarray(structures,x)%position,&
                 &array(structures,atomnumber)%position)).eq.0) cycle
     
     
            m=0
            do k=1, eltot
               do L=1, stochio(k)
                  do n=1, 27
                     !write(*,*) i,j,p,k,l,n
                     !Change from atomnumber to x
                     m=m+1
     
     
                     do ip=1, size(structure_elements)
                        if(trim(adjustl(structure_elements(ip))).ne.&
                             &trim(adjustl(repeatedarray(structures,m)%name))) then
                           cycle
                        else
                           bondint_3=ip
                        end if
                     end do
     
     
                     if(m.eq.x) cycle
                     if(m.eq.atomnumber) cycle
                     if(get_bondlength(repeatedarray(structures,x)%position,&
                          &repeatedarray(structures,m)%position).lt.0.01) cycle
                     if(get_bondlength(repeatedarray(structures,x)%position,&
                          &repeatedarray(structures,m)%position).gt.bondlengthcutoff(bondint_2,bondint_3)) cycle
                     if((N.eq.14).and.(get_bondlength(repeatedarray(structures,m)%position,&
                          &repeatedarray(structures,x)%position)).eq.0) cycle
     
                     LL=0
                     do ii=1, eltot
                        do jj=1, stochio(ii)
                           do kk=1, 27
                              LL=LL+1
                              if(LL.eq.x) cycle
                              if(LL.eq.m) cycle
                              if(LL.eq.atomnumber) cycle
                              !change from atomnumber to x
     
                              do ip=1, size(structure_elements)
                                 if(trim(adjustl(structure_elements(ip))).ne.&
                                      &trim(adjustl(repeatedarray(structures,LL)%name))) then
                                    cycle
                                 else
                                    bondint_4=ip
                                 end if
                              end do
     
     
     
     
     
                              if(get_bondlength(repeatedarray(structures,x)%position,&
                                   &repeatedarray(structures,LL)%position).gt.bondlengthcutoff(bondint_2,bondint_4)) cycle
                              if((kk.eq.14).and.(get_bondlength(repeatedarray(structures,LL)%position,&
                                   &repeatedarray(structures,x)%position)).eq.0) cycle
                              ab=array(structures,atomnumber)%position-repeatedarray(structures,x)%position
                              ac=repeatedarray(structures,x)%position-repeatedarray(structures,m)%position
                              ad=repeatedarray(structures,x)%position-repeatedarray(structures,LL)%position
                              res=get_dihedral_angle(array(structures,atomnumber)%position,&
                                   &repeatedarray(structures,x)%position,&
                                   &repeatedarray(structures,m)%position,&
                                   &repeatedarray(structures,LL)%position)
     
     
                              if(isnan(res)) cycle
                              if(res.eq.1000) cycle
                              if(res.eq.0) cycle
                              counter=counter+1
                              !write(*,*) counter, size(list)
                              if(counter.ne.1) then
                                 allocate(list_tmp(counter-1))
                                 do tool=1, counter-1
                                    list_tmp(tool)=list(tool)
                                 end do
                                 deallocate(list)
                                 allocate(list(counter))
                                 do tool=1, counter-1
                                    list(tool)=list_tmp(tool)
                                 end do
                                 list(counter)=res
                                 deallocate(list_tmp)
                              else 
                                 list(1)=res
                              end if
                              !*get_bondlength(array(structures,atomnumber)%position,&
                              !                                 &repeatedarray(structures,LL)%position)
     
                           end do
                        end do
                     end do
    
     
                  end do
     
               end do
            end do
         end do
      end do
     end do

     !! write to file
     write_loop: do i = 1, counter 
        if(list(i).lt.0.001) cycle
        !! check for duplicates
        do j = 1, i-1, 1
          if(abs(list(i)-list(j)).lt.0.01)then
             list(i)=0
             cycle write_loop
          end if
        end do
      
        if(list(i).gt.pi/2) then 
           write(101,*) list(i)
        else
           write(101,*) list(i)
        end if
     end do write_loop
     
     close(101)
   end subroutine generate4files
!!!#############################################################################


!!!#############################################################################
!!! regenerate distribution files
!!!#############################################################################
   subroutine regenerate_distribution_files (prev_structures)
      implicit none
      character(1024) :: command,name,tmp, location_string
      integer :: l,i,structures,prev_structures,stage, eltot, atomnumber, tmpint,m,j,k, structure_loop
      type (atom), dimension(:,:), allocatable :: array, repeatedarray
      integer, dimension(:), allocatable :: stochio
      logical :: dir_e, file_e, empty
      
      
      
      character(1024) :: coordinate_type, header
      character(LEN=:), allocatable :: read_in
      character(1024) :: pathway, path_test
      character(3), dimension(:), allocatable :: structure_elements
      integer, dimension(:), allocatable :: structure_stochiometry
      integer :: loop, status, size, stat
      integer, parameter :: line_buf_len= 1024*4
      character(LEN=line_buf_len) :: InS
      type(unitcell), dimension(:), allocatable :: unit_cell
      real(real12), dimension(:,:), allocatable :: atomic_positions
      real(real12) :: structure_factor
      real(real12), dimension(:,:), allocatable :: bondcutoff
      logical ::  OK, set
      
      allocate(unit_cell(1))
      write(pathway,*) "rm -r bon bad 4body bon2 bad2 4body2 bon3 bad3 4body3"
      call execute_command_line(pathway)
      strucloop:do structure_loop=1, prev_structures
         !DEBUG
         !strucloop:do structure_loop=303, 304 
       write(pathway,'(A,I0.3,A)') "pos/POSCAR_",structure_loop,"/KPOINTS" 
       !write(*,*) pathway
       inquire(file=trim(adjustl(pathway)),exist=dir_e)
       if(dir_e.eqv..FALSE.) cycle strucloop
       write(pathway,'(A,I0.3,A)') "pos/POSCAR_",structure_loop,"/CONTCAR"
       inquire(file=trim(adjustl(pathway)),exist=dir_e)
       write(command,'(A,A,A,I0.3,A,I0.3,A,I0.3,A,I0.3,A)') &
            &"if [ -s ",trim(adjustl(pathway)), &
            &" ]; then rm -f  pos/POSCAR_",structure_loop,&
            &"/empty.txt; touch pos/POSCAR_",structure_loop&
            &,"/full.txt; else rm -f touch pos/POSCAR_",structure_loop&
            &,"/full.txt; touch pos/POSCAR_",structure_loop,&
            &"/empty.txt; fi"
       write(*,*) command
       call execute_command_line(command)
       write(path_test,'(A,I0.3,A)') "pos/POSCAR_",structure_loop,"/full.txt"
      
       inquire(file=trim(adjustl(path_test)),exist=dir_e)
      
       if(dir_e.eqv..FALSE.) then
          write(pathway,'(A,I0.3,A)') "pos/POSCAR_",structure_loop,"/POSCAR"
          inquire(file=trim(adjustl(pathway)),exist=dir_e)
          if(dir_e.eqv..FALSE.) cycle strucloop
       end if
      
      
       write(*,*) pathway, structure_loop 
      
      
       call poscar_read(pathway,unit_cell(1),structure_elements,structure_stochiometry, structure_factor,coordinate_type&
            &,atomic_positions, header)
       !write(*,*) pathway
       allocate(array(1,sum(structure_stochiometry)))
       allocate(repeatedarray(1,27*sum(structure_stochiometry)))
       k=0
       !write(*,*) size(array)
       do i=1, size(structure_stochiometry)
          do j=1, structure_stochiometry(i) 
             k=k+1
             array(1,k)%name=structure_elements(i)
             array(1,k)%position=atomic_positions(k,:)
          end do
      
      
       end do
      
      
       allocate(bondcutoff(size(structure_stochiometry),size(structure_stochiometry)))
       do i=1, size(structure_stochiometry)
          do j=1, size(structure_stochiometry)
             if(j.gt.i) cycle
             write(command,*) "echo $(sed -n 's|",trim(adjustl(structure_elements(i)))&
                  &," ",trim(adjustl(structure_elements(j)))," ||p'<chem.in)|",&
                  &"sed -n 's| .*||p'>tmp.out"
             call execute_command_line(command) 
             open(72,file="tmp.out")
             read(72,*,IOSTAT=stat) bondcutoff(i,j)
             IF(IS_IOSTAT_END(stat)) then 
                close(72)
                write(command,*) "echo $(sed -n 's|",trim(adjustl(structure_elements(j)))&
                     &," ",trim(adjustl(structure_elements(i)))," ||p'<chem.in)|",&
                     &"sed -n 's| .*||p'>tmp.out"
                call execute_command_line(command)
                open(72,file="tmp.out")
                read(72,*) bondcutoff(i,j)
                bondcutoff(i,j)=bondcutoff(i,j)*1.2
                bondcutoff(j,i)=bondcutoff(i,j)
                close(72)
             else
      
                close(72)
                bondcutoff(i,j)=bondcutoff(i,j)*1.2
                bondcutoff(j,i)=bondcutoff(i,j)
             end if
          end do
       end do
      
      
      
      
       ! write(*,*) "MANUALLY CHANGE"
       ! do i=1, sum(structure_stochiometry)
       !    call atomrepeater(1,array(1,i)%position,array(1,i)%name,&
       !         &repeatedarray,unit_cell,i,sum(structure_stochiometry)) 
       ! end do
      
      
      
      
      
       write(*,*) "YOU SHOULD ADD IN HERE A FUNCTION WHICH CERATES DYNAMIC BONDCUTOFF VALUES FOR DIFFERENT ELEMTNS."
       write(*,*) structure_stochiometry, sum(structure_stochiometry)
       do i=1, sum(structure_stochiometry) 
          call generatebondfiles(1,1,structure_loop-1,array,repeatedarray,size(structure_stochiometry),&             
               &structure_stochiometry,i,structure_elements)
          call generate4files(1,1,structure_loop-1,array,repeatedarray,size(structure_stochiometry),&
               &structure_stochiometry,i,bondcutoff,structure_elements)
      
          call generateanglefiles(1,1,structure_loop-1,array,repeatedarray,size(structure_stochiometry),&
               &structure_stochiometry,i,bondcutoff,structure_elements)
       end do
       deallocate(array)
       deallocate(repeatedarray)
       deallocate(bondcutoff)
      end do strucloop
      
      strucloop2:do structure_loop=1, prev_structures
         !  strucloop2:do structure_loop=303, 304
      
       write(pathway,'(A,I0.3,A)') "pos/POSCAR_",structure_loop,"/RELAX/CONTCAR" 
      
       write(command,'(A,A,A,A,A)') " if [ -s ", trim(adjustl(pathway))," ]; then :; else rm ",trim(adjustl(pathway))," ; fi"
       call execute_command_line(command)
       inquire(file=trim(adjustl(pathway)),exist=dir_e)
      
       if(dir_e.eqv..FALSE.) then 
          write(pathway,'(A,I0.3,A)') "pos/POSCAR_",structure_loop,"/RELAX/POSCAR"
      
          write(command,'(A,A,A,A,A)') " if [ -s ", trim(adjustl(pathway))," ]; then :; else rm ",trim(adjustl(pathway))," ; fi"
          call execute_command_line(command)
          inquire(file=trim(adjustl(pathway)),exist=dir_e)
      
          if(dir_e.eqv..FALSE.) cycle strucloop2
       else
          continue 
       end if
      
      
      
      
      
       call poscar_read(pathway,unit_cell(1),structure_elements,structure_stochiometry, structure_factor,coordinate_type&
            &,atomic_positions, header)
       write(*,*) pathway
       allocate(array(1,sum(structure_stochiometry)))
       allocate(repeatedarray(1,27*sum(structure_stochiometry)))
       k=0
       write(*,*) size(array)
       do i=1, size(structure_stochiometry)
          do j=1, structure_stochiometry(i) 
             k=k+1
             array(1,k)%name=structure_elements(i)
             array(1,k)%position=atomic_positions(k,:)
          end do
       end do
      
       allocate(bondcutoff(size(structure_stochiometry),size(structure_stochiometry)))
       do i=1, size(structure_stochiometry)
          do j=1, size(structure_stochiometry)
             if(j.gt.i) cycle
             write(command,*) "echo $(sed -n 's|",trim(adjustl(structure_elements(i)))&
                  &," ",trim(adjustl(structure_elements(j)))," ||p'<chem.in)|",&
                  &"sed -n 's| .*||p'>tmp.out"
             call execute_command_line(command)
             open(72,file="tmp.out")
             read(72,*,IOSTAT=stat) bondcutoff(i,j)
             IF(IS_IOSTAT_END(stat)) then
                close(72)
                write(command,*) "echo $(sed -n 's|",trim(adjustl(structure_elements(j)))&
                     &," ",trim(adjustl(structure_elements(i)))," ||p'<chem.in)|",&
                     &"sed -n 's| .*||p'>tmp.out"
                call execute_command_line(command)
                open(72,file="tmp.out")
                read(72,*) bondcutoff(i,j)
                bondcutoff(j,i)=bondcutoff(i,j)
                bondcutoff(i,j)=bondcutoff(i,j)*1.1
                close(72)
             else
      
                close(72)
                bondcutoff(j,i)=bondcutoff(i,j)
                bondcutoff(i,j)=bondcutoff(i,j)*1.1
             end if
          end do
       end do
      
       ! do i=1, sum(structure_stochiometry)
       !    call atomrepeater(1,array(1,i)%position,array(1,i)%name,&
       !         &repeatedarray,unit_cell,i,sum(structure_stochiometry)) 
       ! end do
      
       write(*,*) "YOU SHOULD ADD IN HERE A FUNCTION WHICH CERATES DYNAMIC BONDCUTOFF VALUES FOR DIFFERENT ELEMTNS."
       do i=1, sum(structure_stochiometry) 
          write(*,*) prev_structures,structure_loop,prev_structures+structure_loop
          call generatebondfiles(2,1,structure_loop-1,array,repeatedarray,size(structure_stochiometry),&
               &structure_stochiometry,i,structure_elements)
          write(*,*) "11"
          call generate4files(2,1,structure_loop-1,array,repeatedarray,size(structure_stochiometry),&
               &structure_stochiometry,i,bondcutoff,structure_elements)
          write(*,*) "22"
          call generateanglefiles(2,1,structure_loop-1,array,repeatedarray,size(structure_stochiometry),&
               &structure_stochiometry,i,bondcutoff,structure_elements)
          write(*,*) "33"
       end do
       deallocate(array)
       deallocate(repeatedarray)
       deallocate(bondcutoff)
      end do strucloop2
      strucloop3:do structure_loop=1, prev_structures
         ! strucloop3:do structure_loop=303, 304
      
       write(pathway,'(A,I0.3,A)') "pos/POSCAR_",structure_loop,"/RELAX/RELAX2/CONTCAR" 
       write(*,*) pathway
      
       write(command,'(A,A,A,A,A)') " if [ -s ", trim(adjustl(pathway))," ]; then :; else rm ",trim(adjustl(pathway))," ; fi"
       call execute_command_line (command)
      
       inquire(file=trim(adjustl(pathway)),exist=dir_e)
      
      
       if(dir_e.eqv..FALSE.) then 
          write(pathway,'(A,I0.3,A)') "pos/POSCAR_",structure_loop,"/RELAX/RELAX2/POSCAR"
      
          write(command,'(A,A,A,A,A)') " if [ -s ", trim(adjustl(pathway))," ]; then :; else rm ",trim(adjustl(pathway))," ; fi"
          call execute_command_line(command)
          inquire(file=trim(adjustl(pathway)),exist=dir_e)
      
          if(dir_e.eqv..FALSE.) cycle strucloop3
       else
          continue 
       end if
      
       call poscar_read(pathway,unit_cell(1),structure_elements,structure_stochiometry, structure_factor,coordinate_type&
            &,atomic_positions, header)
       write(*,*) pathway
       allocate(array(1,sum(structure_stochiometry)))
       allocate(repeatedarray(1,27*sum(structure_stochiometry)))
       k=0
       write(*,*) size(array)
       do i=1, size(structure_stochiometry)
          do j=1, structure_stochiometry(i) 
             k=k+1
             array(1,k)%name=structure_elements(i)
             array(1,k)%position=atomic_positions(k,:)
          end do
       end do
      
       allocate(bondcutoff(size(structure_stochiometry),size(structure_stochiometry)))
       do i=1, size(structure_stochiometry)
          do j=1, size(structure_stochiometry)
             if(j.gt.i) cycle
             write(command,*) "echo $(sed -n 's|",trim(adjustl(structure_elements(i)))&
                  &," ",trim(adjustl(structure_elements(j)))," ||p'<chem.in)|",&
                  &"sed -n 's| .*||p'>tmp.out"
             call execute_command_line(command)
             open(72,file="tmp.out")
             read(72,*,IOSTAT=stat) bondcutoff(i,j)
             IF(IS_IOSTAT_END(stat)) then
                close(72)
                write(command,*) "echo $(sed -n 's|",trim(adjustl(structure_elements(j)))&
                     &," ",trim(adjustl(structure_elements(i)))," ||p'<chem.in)|",&
                     &"sed -n 's| .*||p'>tmp.out"
                call execute_command_line(command)
                open(72,file="tmp.out")
                read(72,*) bondcutoff(i,j)
                bondcutoff(i,j)=bondcutoff(i,j)*1.1
      
                bondcutoff(j,i)=bondcutoff(i,j)
      
                close(72)
             else
      
                close(72)
                bondcutoff(i,j)=bondcutoff(i,j)*1.1
                bondcutoff(j,i)=bondcutoff(i,j)
             end if
          end do
       end do
       write(*,*) bondcutoff
       stop
      
       ! do i=1, sum(structure_stochiometry)
       !    call atomrepeater(1,array(1,i)%position,array(1,i)%name,&
       !         &repeatedarray,unit_cell,i,sum(structure_stochiometry)) 
       ! end do
      
       write(*,*) "YOU SHOULD ADD IN HERE A FUNCTION WHICH CERATES DYNAMIC BONDCUTOFF VALUES FOR DIFFERENT ELEMTNS."
       do i=1, sum(structure_stochiometry) 
          write(*,*) prev_structures,structure_loop,prev_structures+structure_loop
          call generatebondfiles(3,1,structure_loop-1,array,repeatedarray,size(structure_stochiometry),&
               &structure_stochiometry,i,structure_elements)
          write(*,*) "11"
          call generate4files(3,1,structure_loop-1,array,repeatedarray,size(structure_stochiometry),&
               &structure_stochiometry,i,bondcutoff,structure_elements)
          write(*,*) "22"
          call generateanglefiles(3,1,structure_loop-1,array,repeatedarray,size(structure_stochiometry),&
               &structure_stochiometry,i,bondcutoff,structure_elements)
          write(*,*) "33"
       end do
       deallocate(array)
       deallocate(repeatedarray)
       deallocate(bondcutoff)
      end do strucloop3
      
      end subroutine regenerate_distribution_files
!!!#############################################################################

end module file_generator