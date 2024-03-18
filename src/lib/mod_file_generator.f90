module file_generator
  use constants, only: real12, pi
  use atomtype
  use geom, only: get_bondlength, get_bondangle, get_dihedral_angle
  implicit none

contains


!!!#############################################################################
!!! generate file associated with bondlengths
!!!#############################################################################
  subroutine generatebondfiles( &
       stage,structures,prev_structures,&
       array,repeatedarray,eltot,stochio,atomnumber,structure_elements)
    use atomtype
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

end module file_generator