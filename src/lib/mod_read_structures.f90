module read_structures
  use constants, only: real12
  use misc_raffle, only: grep
  use misc_linalg, only: modu
  use rw_geom, only: bas_type, geom_read, geom_write, igeom_input
  use rw_vasprun, only: get_energy_from_vasprun, get_structure_from_vasprun
  use evolver, only: gvector_container_type, gvector_type
  use machine_learning, only: network_setup, &
       network_train, network_train_graph, &
       network_predict, network_predict_graph
  use athena, only: shuffle, random_setup, split, graph_type, edge_type
  implicit none


  private

  public :: get_evolved_gvectors_from_data
  public :: get_graph_from_basis


contains

!!!#############################################################################
!!! read in the structures from the input directories and generate the gvectors
!!!#############################################################################
  function get_evolved_gvectors_from_data(input_dir, &
       element_file, bond_file, element_list, &
       file_format, gvector_container_template) &
       result(gvector_container)
    implicit none
    character(*), dimension(..), intent(in) :: input_dir
    character(*), intent(in), optional :: element_file, bond_file
    type(gvector_container_type), intent(in), optional :: gvector_container_template
    type(gvector_container_type), allocatable :: gvector_container
    character(len=3), dimension(:), intent(in), optional :: element_list
    character(*), intent(in), optional :: file_format

    character(256) :: name
    integer :: i, j, ifile_format, num_structures
    real(real12) :: energy
    character(50) :: buffer
    character(256) :: format
    logical :: success

    integer :: xml_unit, unit, ierror
    integer :: num_files
    type(bas_type) :: basis
    type(gvector_type) :: gvector
    real(real12), dimension(3,3) :: lattice
    character(256), dimension(:), allocatable :: structure_list
    type(graph_type), dimension(:), allocatable :: graphs

    real(real12), dimension(:), allocatable :: labels, labels_train, labels_validate
    real(real12), dimension(:,:), allocatable :: dataset, data_train, data_validate


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
    allocate(graphs(0))
    allocate(labels(0))
    do i = 1, size(structure_list)

       write(*,*) "Reading structure: ", trim(adjustl(structure_list(i)))
       select case(ifile_format)
       case(0) ! vasprun.xml
          open(newunit=unit, file=trim(adjustl(structure_list(i)))//"/vasprun.xml")
          write(*,*) "Reading structures from xml"
          basis%energy = get_energy_from_vasprun(unit, success)
          if(.not.success) cycle
          rewind(unit)
          call get_structure_from_vasprun(unit, lattice, basis, success)
          if(.not.success) cycle
          close(unit)
       case(1)
          open(newunit=unit, file=trim(adjustl(structure_list(i)))//"/POSCAR")
          write(*,*) "Reading structures from POSCAR"
          call geom_read(unit, lattice, basis)
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
             call geom_read(unit, lattice, basis)
             call get_elements_masses_and_charges(basis)
             graphs = [ graphs, get_graph_from_basis(lattice, basis) ]
             labels = [ labels, basis%energy ]

             num_structures = num_structures + 1
             write(format,'("(""Found structure: "",I4,"" with energy: "",&
                  &F13.7,"" and formula:"",",I0,"(1X,A,I0))")') basis%nspec
             write(*,format) &
                  num_structures, &
                  basis%energy, &
                  ( trim(basis%spec(j)%name), basis%spec(j)%num, &
                  j=1, basis%nspec )
             call gvector_container%add(basis, lattice)
          end do
          cycle
       end select

       write(*,*) &
            "Found structure: ", trim(adjustl(structure_list(i))), &
            " with energy: ", basis%energy
       call gvector_container%add(basis, lattice)
       num_structures = num_structures + 1
      
       !!! STORE THE ENERGY IN AN ARRAY
       ! probably new structure format of crystal
       ! where crystal contains lattice, basis, and energy
       !!! DO SOMETHING ABOUT NESTED RELAXATIONS
    end do


100 if(present(element_file).and.present(element_list)) then
       call gvector_container%set_element_info(element_file, element_list)
    end if
    if(present(bond_file)) then
       call gvector_container%set_bond_info(bond_file)
    else
      call gvector_container%set_bond_info()
    end if



    call gvector_container%evolve(deallocate_systems_after_evolve=.false.)

    !!! do not deallocate structures
    !!! then load the athena library
    !!! set up the network
    !!! append 2, 3, and 4 body potentials
    !!! HOW DO WE HANDLE SPECIES?
    !!! A network for each species?
    !!! split dataset into train and test sets
    !!! train the network

   !  num_structures = size(gvector_container%system)
   !  write(*,*) "LOOKY", gvector_container%nbins
   !  allocate(dataset(sum(gvector_container%nbins), num_structures))
   !  do i = 1, num_structures
   !     dataset(1:gvector_container%nbins(1),i) = &
   !          sum(gvector_container%system(i)%df_2body,dim=2)
   !     dataset(gvector_container%nbins(1)+1:&
   !          sum(gvector_container%nbins(1:2)),i) = &
   !          sum(gvector_container%system(i)%df_3body,dim=2)
   !     dataset(sum(gvector_container%nbins(1:2))+1:&
   !          sum(gvector_container%nbins(1:3)),i) = &
   !          sum(gvector_container%system(i)%df_4body,dim=2)
   !  end do
   !  allocate(labels(num_structures))
   !  labels = [ gvector_container%system(:)%energy / gvector_container%system(:)%num_atoms ]

   !  call random_setup(1)
   !  call split( dataset, labels, &
   !              data_train, data_validate, &
   !              labels_train, labels_validate, &
   !              dim=2, left_size=0.8, right_size=0.2, shuffle=.true., seed=1)

   !  call network_setup(num_inputs = sum(gvector_container%nbins), &
   !       num_outputs = 1)
   !  call network_train(data_train, labels_train, num_epochs = 100)

   !  write(*,*) "predicting known"
   !  write(*,*) -1._real12 * network_predict(data_train(:,1:10)) * sqrt(dot_product(labels_train, labels_train))
   !  write(*,*) labels_train(1:10)
   !  write(*,*)

   !  write(*,*) "PREDICTING"
   !  write(*,*) "norm", sqrt(dot_product(labels_train, labels_train))
   !  write(*,*) -1._real12 * network_predict(data_validate) * sqrt(dot_product(labels_train, labels_train))
    
   !  write(*,*) labels_validate

    ! write(*,*) "LABELS"
    ! write(*,*) size(labels)
    ! write(*,*) labels
    ! call network_setup(num_inputs = 2, num_outputs = 1)
    ! call network_train_graph(graphs(:size(graphs)-10), labels(:size(graphs)-10), num_epochs = 100)

    
   !  write(*,*) "predicting known"
   !  write(*,*) network_predict_graph(graphs(:size(graphs)-10))
   !  write(*,*) labels(:size(graphs)-10)
   !  write(*,*)

   !  write(*,*) "PREDICTING"
   !  write(*,*) network_predict_graph(graphs(size(graphs)-10+1:))
   !  write(*,*) labels(size(graphs)-10+1:)
   !  write(*,*)


    call network_setup(num_inputs = 2, num_outputs = 1)
    call network_train_graph(graphs(:), labels(:), num_epochs = 200)

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
  function get_graph_from_basis(lattice, basis) result(graph)
    implicit none
    type(bas_type), intent(in) :: basis
    real(real12), dimension(3,3), intent(in) :: lattice
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
    amax = ceiling(cutoff_max/modu(lattice(1,:)))
    bmax = ceiling(cutoff_max/modu(lattice(2,:)))
    cmax = ceiling(cutoff_max/modu(lattice(3,:)))

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
                         rtmp1 = modu(matmul(vtmp1,lattice))
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
!!!#############################################################################


!!!#############################################################################
!!! get elements masses and charges
!!!#############################################################################
  subroutine get_elements_masses_and_charges(basis)
    implicit none
    type(bas_type), intent(inout) :: basis

    integer :: i
    real(real12) :: mass, charge

    do i = 1, basis%nspec
       select case(basis%spec(i)%name)
       case('H')
          mass = 1.00784_real12
          charge = 1.0_real12
       case('He')
          mass = 4.0026_real12
          charge = 2.0_real12
       case('Li')
          mass = 6.94_real12
          charge = 3.0_real12
       case('Be')
          mass = 9.0122_real12
          charge = 4.0_real12
       case('B')
          mass = 10.81_real12
          charge = 5.0_real12
       case('C')
          mass = 12.011_real12
          charge = 4.0_real12
       case('N')
          mass = 14.007_real12
          charge = 5.0_real12
       case('O')
          mass = 15.999_real12
          charge = 6.0_real12
       case('F')
          mass = 18.998_real12
          charge = 7.0_real12
       case('Na')
          mass = 22.989_real12
          charge = 1.0_real12
       case('Mg')
          mass = 24.305_real12
          charge = 2.0_real12
       case('Al')
          mass = 26.982_real12
          charge = 3.0_real12
       case('Si')
          mass = 28.085_real12
          charge = 4.0_real12
       case('P')
          mass = 30.974_real12
          charge = 5.0_real12
       case('S')  
          mass = 32.06_real12
          charge = 6.0_real12
       case('Cl')
          mass = 35.453_real12
          charge = 8.0_real12
       case('K')
          mass = 39.098_real12
          charge = 1.0_real12
       case('Ca')
          mass = 40.078_real12
          charge = 2.0_real12
       case('Sc')
          mass = 44.956_real12
          charge = 3.0_real12
       case('Ti')
          mass = 47.867_real12
          charge = 4.0_real12
       case('V')
          mass = 50.942_real12
          charge = 5.0_real12
       case('Cr')
          mass = 51.996_real12
          charge = 6.0_real12
       case('Mn')
          mass = 54.938_real12
          charge = 7.0_real12
       case('Fe')
          mass = 55.845_real12
          charge = 8.0_real12
       case('Co')
          mass = 58.933_real12
          charge = 9.0_real12
       case('Ni')
          mass = 58.693_real12
          charge = 10.0_real12
       case('Cu')
          mass = 63.546_real12
          charge = 11.0_real12
       case('Zn')
          mass = 65.38_real12
          charge = 12.0_real12
       case('Ga')
          mass = 69.723_real12
          charge = 13.0_real12
       case('Ge')
          mass = 72.63_real12
          charge = 14.0_real12
       case('As')
          mass = 74.922_real12
          charge = 15.0_real12
       case('Se')
          mass = 78.971_real12
          charge = 16.0_real12
       case('Br')
          mass = 79.904_real12
          charge = 17.0_real12
       case('Kr')
          mass = 83.798_real12
          charge = 18.0_real12
       case('Rb')
          mass = 85.468_real12
          charge = 19.0_real12
       case('Sr')
          mass = 87.62_real12
          charge = 20.0_real12
       case('Y')
          mass = 88.906_real12
          charge = 21.0_real12
       case('Zr')
          mass = 91.224_real12
          charge = 22.0_real12
       case('Nb')
          mass = 92.906_real12
          charge = 23.0_real12
       case('Mo')
          mass = 95.95_real12
          charge = 24.0_real12
       case('Tc')
          mass = 98.0_real12
          charge = 25.0_real12
       case('Ru')
          mass = 101.07_real12
          charge = 26.0_real12
       case('Rh')
          mass = 102.91_real12
          charge = 27.0_real12
       case('Pd')
          mass = 106.42_real12
          charge = 28.0_real12
       case('Ag')
          mass = 107.87_real12
          charge = 29.0_real12
       case('Cd')
          mass = 112.41_real12
          charge = 30.0_real12
       case('In')
          mass = 114.82_real12
          charge = 31.0_real12
       case('Sn')
          mass = 118.71_real12
          charge = 32.0_real12
       case('Sb')
          mass = 121.76_real12
          charge = 33.0_real12
       case('Te')
          mass = 127.6_real12
          charge = 34.0_real12
       case('I')
          mass = 126.9_real12
          charge = 35.0_real12
       case('Xe')
          mass = 131.29_real12
          charge = 36.0_real12
       case('Cs')
          mass = 132.91_real12
          charge = 37.0_real12
       case('Ba')
          mass = 137.33_real12
          charge = 38.0_real12
       case('La')
          mass = 138.91_real12
          charge = 39.0_real12
       case('Ce')
          mass = 140.12_real12
          charge = 40.0_real12
       case('Pr')
          mass = 140.91_real12
          charge = 41.0_real12
       case('Nd')
          mass = 144.24_real12
          charge = 42.0_real12
       case('Pm')
          mass = 145.0_real12
          charge = 43.0_real12
       case('Sm')
          mass = 150.36_real12
          charge = 44.0_real12
       case('Eu')
          mass = 152.0_real12
          charge = 45.0_real12
       case('Gd')
          mass = 157.25_real12
          charge = 46.0_real12
       case('Tb')
          mass = 158.93_real12
          charge = 47.0_real12
       case('Dy')
          mass = 162.5_real12
          charge = 48.0_real12
       case('Ho')
          mass = 164.93_real12
          charge = 49.0_real12
       case('Er')
          mass = 167.26_real12
          charge = 50.0_real12
       case('Tm')
          mass = 168.93_real12
          charge = 51.0_real12
       case('Yb')
          mass = 173.05_real12
          charge = 52.0_real12
       case('Lu')
          mass = 174.97_real12
          charge = 53.0_real12
       case('Hf')
          mass = 178.49_real12
          charge = 54.0_real12
       case('Ta')
          mass = 180.95_real12
          charge = 55.0_real12
       case('W')
          mass = 183.84_real12
          charge = 56.0_real12
       case('Re')
          mass = 186.21_real12
          charge = 57.0_real12
       case('Os')
          mass = 190.23_real12
          charge = 58.0_real12
       case('Ir')
          mass = 192.22_real12
          charge = 59.0_real12
       case('Pt')
          mass = 195.08_real12
          charge = 60.0_real12
       case('Au')
          mass = 196.97_real12
          charge = 61.0_real12
       case('Hg')
          mass = 200.59_real12
          charge = 62.0_real12
       case('Tl')
          mass = 204.38_real12
          charge = 63.0_real12
       case('Pb')
          mass = 207.2_real12
          charge = 64.0_real12
       case('Bi')
          mass = 208.98_real12
          charge = 65.0_real12
       case('Th')
          mass = 232.04_real12
          charge = 66.0_real12
       case('Pa')
          mass = 231.04_real12
          charge = 67.0_real12
       case('U')
          mass = 238.03_real12
          charge = 68.0_real12
       case('Np')
          mass = 237.0_real12
          charge = 69.0_real12
       case('Pu')
          mass = 244.0_real12
          charge = 70.0_real12
       case('Am')
          mass = 243.0_real12
          charge = 71.0_real12
       case('Cm')
          mass = 247.0_real12
          charge = 72.0_real12
       case('Bk')
          mass = 247.0_real12
          charge = 73.0_real12
       case('Cf')
          mass = 251.0_real12
          charge = 74.0_real12
       case('Es')
          mass = 252.0_real12
          charge = 75.0_real12
       case('Fm')
          mass = 257.0_real12
          charge = 76.0_real12
       case('Md')
          mass = 258.0_real12
          charge = 77.0_real12
       case('No')
          mass = 259.0_real12
          charge = 78.0_real12
       case('Lr')
          mass = 262.0_real12
          charge = 79.0_real12
       case('Rf')
          mass = 267.0_real12
          charge = 80.0_real12
       case('Db')
          mass = 270.0_real12
          charge = 81.0_real12
       case('Sg')
          mass = 271.0_real12
          charge = 82.0_real12
       case('Bh')
          mass = 270.0_real12
          charge = 83.0_real12
       case('Hs')
          mass = 277.0_real12
          charge = 84.0_real12
       case('Mt')
          mass = 276.0_real12
          charge = 85.0_real12
       case('Ds')
          mass = 281.0_real12
          charge = 86.0_real12
       case('Rg')
          mass = 280.0_real12
          charge = 87.0_real12
       case('Cn')
          mass = 285.0_real12
          charge = 88.0_real12
       case('Nh')
          mass = 284.0_real12
          charge = 89.0_real12
       case('Fl')
          mass = 289.0_real12
          charge = 90.0_real12
       case('Mc')
          mass = 288.0_real12
          charge = 91.0_real12
       case('Lv')
          mass = 293.0_real12
          charge = 92.0_real12
       case('Ts')
          mass = 294.0_real12
          charge = 93.0_real12
       case('Og')
          mass = 294.0_real12
          charge = 94.0_real12
       case default
          ! handle unknown element
          mass = 0.0_real12
          charge = 0.0_real12
       end select
       basis%spec(i)%mass = mass
       basis%spec(i)%charge = charge
    end do

  end subroutine get_elements_masses_and_charges
!!!#############################################################################

end module read_structures