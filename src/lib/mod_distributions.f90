module distributions
  use constants, only: real12
  use geom
  use atomtype, only: atom
  implicit none

  type nearmatrix
     real(real12), dimension(3) :: position 
  end type nearmatrix


contains



subroutine angledistribution(atoma,resul)
  integer :: i,j,k,length
  character(3) :: atoma 
  character(1024) :: filename, buffer 
  logical :: dir_e 
  real(real12), dimension(:), allocatable :: resul
  
  if(allocated(resul)) deallocate(resul)
  write(filename,*) "angledist", trim(adjustl(atoma)),".in"
  filename=trim(adjustL(filename))
  inquire(file=trim(adjustl(filename)), exist=dir_e)
  if(dir_e) then 
     !write(*,*) "Found distribution for ", atoma
  else 
   write(*,*) "STOP. YOU VIOLATED THE LAW. YOU MUST HAVE DISTRIBUTIONS FOR", atoma
   stop 
  end if
  open(50,file=filename, status="old")
  j=0
  read(50,"(1024A)") buffer
  close(50)
  k=0
  length=0
  loop1: do i=1, 1024 
   if(i.lt.k) cycle
   if((scan(buffer(i:i)," ").eq.1)) cycle  
   loop2: do j=0, 1024-i
      if((scan(buffer(i:i+j)," ").gt.0)) then 
         k=i+j
         length=length+1
         exit loop2
      end if
  
   end do loop2
  end do loop1
  allocate(resul(length))
  length=0
  j=0
  k=0
  loop1p: do i=1, 1024
   if(i.lt.k) cycle
   if((scan(buffer(i:i)," ").eq.1)) cycle
   loop2p: do j=0, 1024-i
      if((scan(buffer(i:i+j)," ").gt.0)) then
         length=length+1
         read(buffer(i:i+j-1),*) resul(length)
  
         k=i+j
         exit loop2p
      end if
  
   end do loop2p
  end do loop1p
  
  end subroutine angledistribution

subroutine generate_bondangle_distribution (tmpvector,atomlist,alistrep,structures,bondcutoff, &
  &atom_number_new, new_position_needed&
  &,agausssamp,distribution,eltot, sigma1, sigma2,elnames, predicted_positions)
 character(3), dimension(:), allocatable :: elnames
 real(real12), dimension(3) :: tmpvector
 integer :: y, atom_number_new, structures, eltot,i,j,k,L, nlines 
 logical :: new_position_needed, file_exists
 type (atom), dimension(:,:), allocatable :: atomlist, alistrep, predicted_positions
 type (nearmatrix), dimension(:), allocatable :: nearneighbourmatrix
 real(real12) :: agausssamp, bondcutoff, bondangle_store, distribution, read_angle, sigma1, sigma2, step_significance
 real(real12) :: tmpvalue, total_significance
 character(1024) :: determined_anglefile
 real(real12), allocatable, dimension(:) :: initial_angle_peaks
 
 distribution=0
 !if(allocated(nearneighbourmatrix)) deallocate(nearneighbourmatrix) 
 allocate(nearneighbourmatrix(((atom_number_new-1)*27)+26))
 k=0
 
 do y=1, ((atom_number_new-1)*27)+26
  if(y.le.((atom_number_new-1)*27)) then 
     if(get_bondlength(tmpvector,alistrep(structures,y)%position).lt.bondcutoff) then
        k=k+1
        nearneighbourmatrix(k)%position=alistrep(structures,y)%position
     end if
  else
     if(get_bondlength(tmpvector,predicted_positions(1,y-(atom_number_new-1)*27)%position).lt.bondcutoff) then
        k=k+1
        nearneighbourmatrix(k)%position=predicted_positions(1,y-(atom_number_new-1)*27)%position
     end if
  end if
 end do
 
 
 !!! CHANGE alistrep etc to nearneighbourmatrix%position. SHould work.
 total_significance=0
 if(new_position_needed.neqv..TRUE.) then 
  do L=1, eltot
     !write(*,*) elnames(L), structures, atom_number_previous, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     !write(*,*) atomlist(structures,atom_number_previous+1)%name, "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
     if(trim(adjustl(elnames(L))).ne.trim(adjustl(atomlist(structures,atom_number_new)%name))) cycle
 
     do j=1, k 
        do y=1, k 
           if (j.eq.y) cycle 
           bondangle_store=get_bondangle(nearneighbourmatrix(j)%position,tmpvector,&
                &nearneighbourmatrix(k)%position)
           if (get_bondlength(tmpvector,nearneighbourmatrix(k)%position).lt.1.0) then 
              new_position_needed=.TRUE.
              cycle
           end if
           if (get_bondlength(tmpvector,nearneighbourmatrix(k)%position).lt.1.0) then
              new_position_needed=.TRUE.
              cycle
           end if
 
           if(new_position_needed.eqv..TRUE.) cycle 
 
           call angledistribution(elnames(L),initial_angle_peaks)
           do i=1, size(initial_angle_peaks)
              step_significance=1
              total_significance=total_significance+step_significance
              distribution=distribution+(1.0/2.0)*((erf((tmpvalue-initial_angle_peaks(i)+&
                   &(0.5*agausssamp))/(sigma2*sqrt(2.0)))-&
                   &erf((tmpvalue-initial_angle_peaks(i)-(0.5*agausssamp))/(sigma2*sqrt(2.0)))))
           end do
 
           write(determined_anglefile,'(A,A)') trim(adjustL(elnames(L))),"_evolved_angles"
           INQUIRE(FILE=determined_anglefile, EXIST=file_exists)
           if (file_exists.eqv..TRUE.) then
              nlines=0
              OPEN (11, file=trim(adjustl(determined_anglefile)))
              DO
                 READ (11,*, END=10) read_angle, step_significance
                 total_significance=total_significance+step_significance
                 distribution=distribution+step_significance*(1.0/2.0)*(&
                      &abs(erf((bondangle_store-read_angle+0.5*sigma1)/((sigma1)*2**dble(1.0/2.0)))-&
                      &erf((bondangle_store-read_angle-0.5*sigma1)/((sigma1)*2**dble(1.0/2.0)))))
                 nlines=nlines+1
              END DO
 10           CLOSE (11)
           end if
        end do
     end do
  end do
  distribution=abs(distribution/total_significance) 
 end if
 
 
 end subroutine generate_bondangle_distribution
 
 
 
 
 subroutine generate_bondlength_distribution (sigma1,sigma2,atom_number_proposed,distribution,atomlist,alistrep,&
  &tmpvector,elrad,eltot,structures,elnames,new_position_needed, predicted_positions)
 real(real12), dimension(:,:,:), allocatable :: elrad
 integer :: nlines, atom_number_proposed, atom_number_previous, atom_total_repeated,j,k,l,eltot,structures
 type (atom), dimension(:,:), allocatable :: atomlist, alistrep, predicted_positions
 real(real12), dimension(3) :: tmpvector
 real(real12) :: distribution, total_significance, step_significance, bondlength_store, sigma1,sigma2, read_bondlength 
 character(1024) :: determined_bondfile 
 logical :: file_exists, new_position_needed
 character(3), dimension(:), allocatable :: elnames
 
 atom_number_previous=atom_number_proposed-1
 atom_total_repeated=atom_number_previous*27
 total_significance=0
 new_position_needed=.FALSE.
 do while(new_position_needed.eqv..FALSE.)
  do j=1, atom_total_repeated
     do k=1, eltot 
        do L=1, eltot
           !! Two eltot loops in here to make sure the atom being added and the atom being checked are the correct species 
           if(atomlist(structures,atom_number_proposed)%name.ne.elnames(k)) cycle
           if(alistrep(structures,j)%name.ne.elnames(L)) cycle
           step_significance=1
           bondlength_store=get_bondlength(tmpvector,alistrep(structures,j)%position)
           if (bondlength_store.lt.0.7*elrad(1,L,k)) then 
              new_position_needed=.TRUE.
           end if
           total_significance=total_significance+step_significance
           distribution=distribution+(1.0/2.0)*(&
                &abs(erf((bondlength_store-elrad(1,k,L)+0.5*sigma2)/((sigma1)*2**dble(1.0/2.0)))-&
                &erf((bondlength_store-elrad(1,k,L)-0.5*sigma2)/((sigma1)*2**dble(1.0/2.0)))))
           file_exists=.FALSE.
 
           write(determined_bondfile,'(A,A,A,A)') &
                &trim(adjustL(elnames(k))),"_",trim(adjustL(elnames(L))),"_evolved_bondlength"
           INQUIRE(FILE=determined_bondfile, EXIST=file_exists)           
           if (file_exists .eqv..FALSE.) then  
              write(determined_bondfile,'(A,A,A,A)') &
                   &trim(adjustL(elnames(L))),"_",trim(adjustL(elnames(K))),"_evolved_bondlength"
              INQUIRE(FILE=determined_bondfile, EXIST=file_exists)
           end if
 
           if (file_exists .eqv..TRUE.) then 
              nlines=0
              OPEN (11, file=determined_bondfile)
              DO
                 READ (11,*, END=10) read_bondlength, step_significance
                 !           write(*,*) distribution, nlines, total_significance 
                 total_significance=total_significance+step_significance
                 distribution=distribution+step_significance*(1.0/2.0)*(&
                      &abs(erf((bondlength_store-read_bondlength+0.5*sigma2)/((sigma1)*2**dble(1.0/2.0)))-&
                      &erf((bondlength_store-read_bondlength-0.5*sigma2)/((sigma1)*2**dble(1.0/2.0)))))
                 nlines=nlines+1
              END DO
 10           CLOSE (11)
           end if
        end do
     end do
  end do
  !distribution=abs(distribution/total_significance)
  exit
 end do
 
 do while(new_position_needed.eqv..FALSE.)
  do j=1, 26
     do k=1, eltot
        !! Two eltot loops in here to make sure the atom being added and the atom being checked are the correct species
        if(atomlist(structures,atom_number_proposed)%name.ne.elnames(k)) cycle
        step_significance=1
        bondlength_store=get_bondlength(tmpvector,predicted_positions(1,j)%position)
        if (bondlength_store.lt.1.0) new_position_needed=.TRUE.
        total_significance=total_significance+step_significance
        distribution=distribution+(1.0/2.0)*(&
             &abs(erf((bondlength_store-elrad(1,k,k)+0.5*sigma2)/((sigma1)*2**dble(1.0/2.0)))-&
             &erf((bondlength_store-elrad(1,k,k)-0.5*sigma2)/((sigma1)*2**dble(1.0/2.0)))))
        file_exists=.FALSE.
 
        write(determined_bondfile,'(A,A,A,A)') &
             &trim(adjustL(elnames(k))),"_",trim(adjustL(elnames(k))),"_evolved_bondlength"
        INQUIRE(FILE=determined_bondfile, EXIST=file_exists)
        if (file_exists .eqv..TRUE.) then
           OPEN (11, file=determined_bondfile)
           DO
              READ (11,*, END=11) read_bondlength, step_significance
              !           write(*,*) distribution, nlines, total_significance
              total_significance=total_significance+step_significance
              distribution=distribution+step_significance*(1.0/2.0)*(&
                   &abs(erf((bondlength_store-read_bondlength+0.5*sigma2)/((sigma1)*2**dble(1.0/2.0)))-&
                   &erf((bondlength_store-read_bondlength-0.5*sigma2)/((sigma1)*2**dble(1.0/2.0)))))
              nlines=nlines+1
           END DO
 11        CLOSE (11)
        end if
     end do
  end do
  distribution=abs(distribution/total_significance)
 
  exit
 end do
 
 end subroutine generate_bondlength_distribution
 
 

end module distributions