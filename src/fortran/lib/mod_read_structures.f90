module read_structures
  use constants, only: real12
  use misc_raffle, only: grep
  use misc_linalg, only: modu
  use rw_geom, only: basis_type, geom_read, geom_write, igeom_input
  use rw_vasprun, only: get_energy_from_vasprun, get_structure_from_vasprun
  use evolver, only: gvector_container_type
#ifdef ENABLE_ATHENA
  use machine_learning, only: network_setup, &
       network_train, network_train_graph, &
       network_predict, network_predict_graph
  use athena, only: shuffle, random_setup, split, graph_type, edge_type
#endif
  implicit none


  private

  public :: get_evolved_gvectors_from_data
#ifdef ENABLE_ATHENA
  public :: get_graph_from_basis
#endif


contains

!!!#############################################################################
!!! read in the structures from the input directories and generate the gvectors
!!!#############################################################################
  function get_evolved_gvectors_from_data(input_dir, &
       file_format, gvector_container_template) &
       result(gvector_container)
    implicit none
    character(*), dimension(..), intent(in) :: input_dir
    type(gvector_container_type), intent(in), optional :: gvector_container_template
    type(gvector_container_type), allocatable :: gvector_container
    character(*), intent(in), optional :: file_format

    character(256) :: name
    integer :: i, j, ifile_format, num_structures
    real(real12) :: energy
    character(50) :: buffer
    character(256) :: format
    logical :: success

    integer :: xml_unit, unit, ierror
    integer :: num_files
    type(basis_type) :: basis
    character(256), dimension(:), allocatable :: structure_list
#ifdef ENABLE_ATHENA
    type(graph_type), dimension(:), allocatable :: graphs

    real(real12), dimension(:), allocatable :: labels, labels_train, labels_validate
    real(real12), dimension(:,:), allocatable :: dataset, data_train, data_validate
#endif


    if(present(gvector_container_template)) then
       gvector_container = gvector_container_template
    else
       gvector_container = gvector_container_type()
    end if

    if(present(file_format)) then
       select case(trim(adjustl(file_format)))
       case('vasprun.xml','xml','vasprun')
          ifile_format = 0
       case('POSCAR','OUTCAR')
          ifile_format = 1
          igeom_input = 1
       case('xyz','extxyz')
          ifile_format = 2
          igeom_input = 6
       case default
          write(*,*) "Unknown file format: ", file_format
          stop
       end select
    else
       ifile_format = 0
    end if

    !inquire(file='gvector_container.dat', exist=success)
    !if(success) then
    !   call gvector_container%read('gvector_container.dat')
    !   goto 100
    !end if
    !!! For each new run of the code, it should populate a new directory and ...
    !!! ... add to the existing ones.
    !!! And it should check that output_dir never equals any of the input_dirs (or database_dirs)
    select rank(input_dir)
    rank(0)
       structure_list = [ get_structure_list( input_dir, ifile_format ) ]
    rank(1)
       do i = 1, size(input_dir)
         if(i.eq.1)then
            structure_list = [ get_structure_list( input_dir(i), ifile_format ) ]
         else
            structure_list = [ structure_list, &
                               get_structure_list( input_dir(i), ifile_format )]
         end if
       end do
    end select


    num_structures = 0
#ifdef ENABLE_ATHENA
    allocate(graphs(0))
    allocate(labels(0))
#endif
    do i = 1, size(structure_list)

       write(*,*) "Reading structure: ", trim(adjustl(structure_list(i)))
       select case(ifile_format)
       case(0) ! vasprun.xml
          open(newunit=unit, file=trim(adjustl(structure_list(i)))//"/vasprun.xml")
          write(*,*) "Reading structures from xml"
          basis%energy = get_energy_from_vasprun(unit, success)
          if(.not.success) cycle
          rewind(unit)
          call get_structure_from_vasprun(unit, basis, success)
          if(.not.success) cycle
          close(unit)
       case(1)
          open(newunit=unit, file=trim(adjustl(structure_list(i)))//"/POSCAR")
          write(*,*) "Reading structures from POSCAR"
          call geom_read(unit, basis)
          close(unit)
          open(newunit=unit, file=trim(adjustl(structure_list(i)))//"/OUTCAR")
          call grep(unit, 'free  energy   TOTEN  =', lline=.false., success=success)
          if(.not.success) cycle
          backspace(unit)
          read(unit,*) buffer, buffer, buffer, buffer, energy
          close(unit)
          basis%energy = energy
       case(2)
          open(newunit=unit, file=trim(adjustl(structure_list(i))))
          write(*,*) "Reading structures from xyz"
          do
             read(unit,'(A)',iostat=ierror) buffer
             if(ierror .ne. 0) exit
             if(trim(buffer).eq."") cycle
             backspace(unit)
             call geom_read(unit, basis)
#ifdef ENABLE_ATHENA
             graphs = [ graphs, get_graph_from_basis(basis) ]
             labels = [ labels, basis%energy ]
#endif

             num_structures = num_structures + 1
             write(format,'("(""Found structure: "",I4,"" with energy: "",&
                  &F13.7,"" and formula:"",",I0,"(1X,A,I0))")') basis%nspec
             write(*,format) &
                  num_structures, &
                  basis%energy, &
                  ( trim(basis%spec(j)%name), basis%spec(j)%num, &
                  j=1, basis%nspec )
             call gvector_container%add(basis)
          end do
          cycle
       end select

       write(*,*) &
            "Found structure: ", trim(adjustl(structure_list(i))), &
            " with energy: ", basis%energy
       call gvector_container%add(basis)
       num_structures = num_structures + 1
      
       !!! STORE THE ENERGY IN AN ARRAY
       ! probably new structure format of crystal
       ! where crystal contains lattice, basis, and energy
       !!! DO SOMETHING ABOUT NESTED RELAXATIONS
    end do



100 call gvector_container%evolve(deallocate_systems_after_evolve=.false.)


#ifdef ENABLE_ATHENA
    call network_setup(num_inputs = 2, num_outputs = 1)
    call network_train_graph(graphs(:), labels(:), num_epochs = 200)
#endif

    igeom_input = 1

  end function get_evolved_gvectors_from_data
!!!#############################################################################


!!!#############################################################################
!!! get the list of structures from the input directories
!!!#############################################################################
  function get_structure_list(input_dir, ifile_format) result(structure_list)
    implicit none
    character(*), intent(in) :: input_dir
    integer, intent(in) :: ifile_format
    character(256), dimension(:), allocatable :: structure_list

    integer :: unit, ierror
    character(256) :: name
    logical :: file_exists, addit_file_exists, structures_found


    structures_found = .false.
    call execute_command_line( &
         "ls "//trim(adjustl(input_dir))//"/. >structures.txt", wait = .TRUE. )
    open(newunit=unit, file="structures.txt", status="old")
    do
       read(unit,'(A)',iostat=ierror) name
       if(is_iostat_end(ierror)) exit
       name = trim(adjustl(input_dir))//"/"//trim(adjustl(name))
       select case(ifile_format)
       case(0)
          inquire(file=trim(name)//"/vasprun.xml", exist=file_exists)
       case(1)
          inquire(file=trim(name)//"/POSCAR", exist=file_exists)
          inquire(file=trim(name)//"/OUTCAR", exist=addit_file_exists)
          file_exists = file_exists.and.addit_file_exists
       case(2)
          inquire(file=trim(name), exist=file_exists)
          file_exists = file_exists.and.(index(name, ".xyz") .gt. 0)
          if(.not.file_exists)then
             inquire(file=trim(name)//"/data.xyz", exist=file_exists)
             if(file_exists) name = trim(name)//"/data.xyz"
          end if
       end select
       if(.not.file_exists) cycle
       if(structures_found) then
          structure_list = [ structure_list, name ]
       else
          structure_list = [ name ]
          structures_found = .true.
       end if
    end do
    close(unit)

    if(.not.structures_found) then
       write(*,'("No structures found in directory: """,A, &
            &""" with format: ",I0)') &
            trim(input_dir), ifile_format
       stop 0
    end if

  end function get_structure_list
!!!#############################################################################


!!!#############################################################################
!!! 
!!!#############################################################################
#ifdef ENABLE_ATHENA
  function get_graph_from_basis(basis) result(graph)
    implicit none
    type(basis_type), intent(in) :: basis
    type(graph_type) :: graph

    integer :: is, ia, js, ja, i, j, k
    integer :: iatom, jatom
    integer :: amax, bmax, cmax
    type(edge_type) :: edge
    real(real12) :: rtmp1, cutoff_min, cutoff_max
    real(real12), dimension(3) :: diff, vtmp1

    
    graph%num_vertices = basis%natom
    graph%num_vertex_features = 2
    graph%num_edge_features = 1

    allocate(graph%vertex(graph%num_vertices))

    iatom = 0
    do is = 1, basis%nspec
       do ia = 1, basis%spec(is)%num
          iatom = iatom + 1
          allocate(graph%vertex(iatom)%feature(graph%num_vertex_features))
          graph%vertex(iatom)%feature = [ basis%spec(is)%charge / 100._real12, &
               basis%spec(is)%mass / 52._real12 ]
       end do
    end do

    cutoff_min = 0.5_real12
    cutoff_max = 6.0_real12
    amax = ceiling(cutoff_max/modu(basis%lat(1,:)))
    bmax = ceiling(cutoff_max/modu(basis%lat(2,:)))
    cmax = ceiling(cutoff_max/modu(basis%lat(3,:)))

    iatom = 0
    allocate(graph%edge(0))
    spec_loop1: do is=1,basis%nspec
       atom_loop1: do ia=1,basis%spec(is)%num
          iatom = iatom + 1
          jatom = 0
          spec_loop2: do js=is,basis%nspec
             atom_loop2: do ja=1,basis%spec(js)%num
                jatom = jatom + 1
                if(is.eq.js.and.ja.lt.ia) cycle atom_loop2
                diff = basis%spec(is)%atom(ia,:3) -  basis%spec(js)%atom(ja,:3)
                diff = diff - ceiling(diff - 0.5_real12)
                do i=-amax,amax+1,1
                   vtmp1(1) = diff(1) + real(i, real12)
                   do j=-bmax,bmax+1,1
                      vtmp1(2) = diff(2) + real(j, real12)
                      do k=-cmax,cmax+1,1
                         vtmp1(3) = diff(3) + real(k, real12)
                         rtmp1 = modu(matmul(vtmp1,basis%lat))
                         if( rtmp1 .gt. cutoff_min .and. &
                             rtmp1 .lt. cutoff_max )then
                            edge%index = [iatom,jatom]
                            edge%feature = [rtmp1]
                            graph%edge = [ graph%edge, edge ]
                         end if
                      end do
                   end do
                end do
             end do atom_loop2
          end do spec_loop2
       end do atom_loop1
    end do spec_loop1
    graph%num_edges = size(graph%edge)
    call graph%generate_adjacency()
    call graph%calculate_degree()


  end function get_graph_from_basis
#endif
!!!#############################################################################

end module read_structures