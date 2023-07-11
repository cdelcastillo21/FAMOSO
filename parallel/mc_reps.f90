! Carlos del-Castillo-Negrete, April 2019
! MC_GRID module
! Module implementing functions corresponding to the mc_reps type.
module mc_reps
    use mpi
    use utilities
    use mc_grid 
    implicit none

    type particle_rep
            integer :: np                                   ! # of particles
            integer :: nd                                   ! # of dimensions 
            real(8), dimension(:,:), allocatable :: p_i     ! p_i(i,j) contains ith coordinate value of particle j 
    end type particle_rep

    type stat_rep
            type(grid) :: gd                                ! Grid over which f_yx is defined
            real(8), dimension(:,:), allocatable :: f_yx    ! Distribution function over given gd
    end type stat_rep

    abstract interface
        ! funx_x_i defines a general real valued function over nd reals 
        ! func_x_i(args(nd), args(nd-1), ..., args(1)) : R(nd) => R(1)
        function func_x_i (nd, args)
                real(8) :: func_x_i
                integer, intent(IN) :: nd
                real(8), dimension(nd), intent(IN) :: args
        end function func_x_i
    end interface

contains

    subroutine mc_destroy_srep(s_rep)
            type(stat_rep), intent(INOUT) :: s_rep

            ! Free memory associated with s_rep
            call mc_grid_destroy(s_rep%gd)
            deallocate(s_rep%f_yx)
    end subroutine mc_destroy_srep
    
    subroutine mc_destroy_prep(p_rep)
            type(particle_rep), intent(INOUT) :: p_rep

            ! Free memory assoicated with p_rep
            deallocate(p_rep%p_i)
    end subroutine mc_destroy_prep

    subroutine mc_prep_print( p_rep )
        
        type(particle_rep) :: p_rep

        integer :: i
        
        write (*,*) "-----------"
        write(*,*) "Printing ", p_rep%np, " particle positions:"
        do i=1,p_rep%np
            write(*,*) "p_i(", i,") = ", "(", p_rep%p_i(1,i), ", ", p_rep%p_i(2,i), ")"
        end do

    end subroutine mc_prep_print

    subroutine mc_srep_print( s_rep )
        type(stat_rep) :: s_rep

        write(*,*) "--------------------"
        write(*,*) "Printing srep:"
        call mc_grid_print( "s_rep defined over grid: ", s_rep%gd )
        call print_matrix("f_yx", s_rep%gd%szs(1), s_rep%gd%szs(2), s_rep%f_yx)
    end subroutine mc_srep_print

    ! Calculates integral of given distributuion function in s_rep over given grid
    function mc_srep_norm(s_rep)
        
        real(8) :: mc_srep_norm
        type(stat_rep), intent(IN) :: s_rep

        integer :: i, j, nx, ny, sx, sy
        real(8) :: norm, delta_v

        norm = 0.0
        delta_v = s_rep%gd%delta_x_i(1)*s_rep%gd%delta_x_i(2)
        nx = s_rep%gd%szs(2)
        ny = s_rep%gd%szs(1)
        sy = SIZE(s_rep%f_yx, DIM=1)
        sx = SIZE(s_rep%f_yx, DIM=2)

        call mc_grid_print("s_rep-grid in srep_norm", s_rep%gd)
        write(*,*) size(s_rep%f_yx)
        write(*,*) nx, ny

        do j=1,(nx-1)
            do i=1,(ny-1)
                ! calculate normalizing factor
                
                ! Quadrature to integrate f_yx over grid
                if (j<sx .and. i<sy) then
                  norm = norm+delta_v*(s_rep%f_yx(i,j)+&
                        &s_rep%f_yx(i+1,j)+s_rep%f_yx(i,j+1)+s_rep%f_yx(i+1,j+1))/4.0
                end if 
            end do
        end do

        mc_srep_norm = norm 
    end function mc_srep_norm 

    function mc_copy_srep( s_rep )
        type(stat_rep) :: mc_copy_srep, s_rep_copy
        type(stat_rep), intent(IN) :: s_rep

        integer :: i,j

        allocate(s_rep_copy%f_yx(s_rep%gd%szs(1),s_rep%gd%szs(2)))
        s_rep_copy%gd = mc_grid_copy(s_rep%gd)

        do j=1,s_rep%gd%szs(2)
            do i=1, s_rep%gd%szs(1)
                s_rep_copy%f_yx(i,j) = s_rep%f_yx(i,j)
            end do
        end do

        mc_copy_srep = s_rep_copy
    end function mc_copy_srep

    subroutine mc_e_thresh(f_yx, nx, ny, e_thresh)
        real(8), dimension(:,:), intent(INOUT) :: f_yx
        integer, intent(IN) :: nx,ny
        real(8), intent(IN) :: e_thresh

        integer :: i,j,max_idx, n_elements
        real(8) :: curr_norm, norm, max_el
        real(8), dimension(3) :: temp
        real(8), dimension(:,:), allocatable :: f_sort

        allocate(f_sort(3,nx*ny))

        ! Find norm of f_yx, store f_sort
        norm = 0.0
        do j=1,ny
            do i=1,nx
                norm = norm + f_yx(i,j)*f_yx(i,j)
                f_sort(1,(j-1)*nx + i) = i
                f_sort(2,(j-1)*nx + i) = j
                f_sort(3,(j-1)*nx + i) = f_yx(i,j)
            end do
        end do

        ! Set f_yx to zero
        f_yx(:,:) = 0.0

        ! Sort f_sort keeping track of idxs
        n_elements = nx*ny
        curr_norm = 0.0
        do i=1,nx*ny

            ! Set first in list as max element
            max_el = abs(f_sort(3,i))
            max_idx = i

            ! Find max element in remaining list
            do j=i+1,nx*ny
                if (abs(f_sort(3,j)) > max_el) then
                    max_idx = j
                    max_el = abs(f_sort(3,j))
                end if
            end do
            
            ! Swap to front of list
            temp(:) = f_sort(:,i)
            f_sort(:,i) = f_sort(:,max_idx)
            f_sort(:,max_idx) = temp(:)

            ! Accumulate current norm once next max element is added.
            ! Break and stop sorting if e_thresh exceeded. 
            curr_norm = curr_norm + max_el*max_el
            if (curr_norm/norm > e_thresh) then
                n_elements = i
                write(*,*) "n_elements truncated at ", n_elements
                exit
            end if 
        end do

        write(*,*) "n_elements set at ", n_elements
        ! Add sorted first n_elements of f_isort back to f_yx at right idxs
        do i=1,n_elements
            f_yx(INT(f_sort(1,i)),INT(f_sort(2,i))) = f_sort(3,i)
        end do
       
        deallocate(f_sort)
    end subroutine mc_e_thresh

    subroutine mc_reps_scale_all(s_rep_i, max_gd, N)
        type(stat_rep), dimension(:), intent(INOUT) :: s_rep_i
        type(grid), intent(INOUT) :: max_gd
        integer, intent(IN) :: N

        integer :: i,j,nd,max_sz
        integer, dimension(:,:), allocatable :: nz_idxs
        real(8), dimension(:,:), allocatable :: f_yx_temp

        max_gd = mc_grid_copy(s_rep_i(1)%gd)
        nd = max_gd%nd

        ! First find max grid that contains all s_reps.
        do i=1,N
            do j=1,nd
                if (max_gd%upper_bnds(j) < s_rep_i(i)%gd%upper_bnds(j)) then
                    max_gd%upper_bnds(j) = s_rep_i(i)%gd%upper_bnds(j)
                end if
                
                if (max_gd%lower_bnds(j) >  s_rep_i(i)%gd%lower_bnds(j)) then
                    max_gd%lower_bnds(j) = s_rep_i(i)%gd%lower_bnds(j)
                end if
            end do
        end do

        ! Calculate new # grid points to fit into new bounds while preserving old delta_x spacing 
        do j=1,nd
            max_gd%szs(j) = int(abs(max_gd%upper_bnds(j) - max_gd%lower_bnds(j))/max_gd%delta_x_i(j)) 
        end do 

        ! Allocate space for temporary dist_function
        allocate(f_yx_temp(max_gd%szs(1),max_gd%szs(2)))
        allocate(nz_idxs(2,nd))
        
        ! Next extend all distributuion functions to max grid size
        ! by setting value of f_yx outside original grid to 0. 
        do i=1,N

            ! Calculate indices in grid where f_yx for current s_rep will be non-zero
            ! in the new f_yx_temp defined over potentialy larger grid.
            do j=1,nd
                max_sz = SIZE(s_rep_i(i)%f_yx, DIM=j)
                nz_idxs(1,j) = int((s_rep_i(i)%gd%lower_bnds(j)-max_gd%lower_bnds(j))/max_gd%delta_x_i(j)) + 1
                nz_idxs(2,j) = max_gd%szs(j) - int((max_gd%upper_bnds(j)-s_rep_i(i)%gd%upper_bnds(j))/max_gd%delta_x_i(j))
                if ((1+nz_idxs(2,j)-nz_idxs(1,j))>max_sz) then
                  nz_idxs(2,j) = nz_idxs(2,j)-1
                end if 
            end do

            ! Copy over to temp, extrapolating to zero
            f_yx_temp(:,:) = 0.0
            f_yx_temp(nz_idxs(1,1):nz_idxs(2,1),nz_idxs(1,2):nz_idxs(2,2)) = s_rep_i(i)%f_yx(:,:)

            ! Reallocate space for new s_rep distribution function
            deallocate(s_rep_i(i)%f_yx)
            allocate(s_rep_i(i)%f_yx(max_gd%szs(1),max_gd%szs(2)))

            ! Copy over f_yx_temp
            s_rep_i(i)%f_yx(:,:) = f_yx_temp(:,:)

            ! Set s_rep_i's new grid
            s_rep_i(i)%gd = max_gd
        end do

        ! Free memory
        deallocate(f_yx_temp)
        deallocate(nz_idxs)

    end subroutine mc_reps_scale_all

    ! Constructs a monte-carlo particle representation (xp, yp) of a contium
    ! field given as a matrix (f) in a grid (xg, yg).
    function mc_f_to_xy(s_rep, np)
        ! Declare IN/OUT variables 			
        type(particle_rep) :: mc_f_to_xy, p_rep
        type(stat_rep) :: s_rep
        integer :: np
        
        ! Declare Auxiliary variables
        integer :: i,j,k,l
        integer :: nd, nx, ny, np_j
        real(8), dimension(:,:), allocatable :: np_ij
        real(8) :: delta_v, x_i, y_i, eta, np_in_cell, np_total, norm
        
        ! Initialize variables
        nx = s_rep%gd%szs(2)
        ny = s_rep%gd%szs(1)
        delta_v = s_rep%gd%delta_x_i(1)*s_rep%gd%delta_x_i(2)
        nd = s_rep%gd%nd

        ! Allocate space for particle arrays 
        allocate(np_ij(ny,nx))
        
        call random_seed()

        ! Cont # of particles per cell. Approximate
        ! the integral of f using simple quadrature:
        ! f_cm = (f(i,j)+f(i+1,j)+f(i,j+1)+f(i+1,j+1))
        np_total = 0
        np_ij(:,:) = 0.0
        do j=1,(nx-1)
            do i=1,(ny-1)
                ! calculate number of particles in cell
                np_in_cell = real(np)*delta_v*(s_rep%f_yx(i,j)+s_rep%f_yx(i+1,j)+s_rep%f_yx(i,j+1)+s_rep%f_yx(i+1,j+1))/4.0
                np_ij(i,j) = floor(np_in_cell)
                np_total = np_total + floor(np_in_cell)
            end do
        end do

        ! Allocate particle structure
        allocate(p_rep%p_i(2,INT(np_total)))
        p_rep%nd = nd
        p_rep%np = np_total

        ! Initialze particles with proper coordinates
        k=0
        x_i = s_rep%gd%lower_bnds(2)
        do j=1,(nx-1)
            y_i = s_rep%gd%lower_bnds(1)
            np_j = 0
            do i=1,(ny-1)
                if (np_ij(i,j) > 0) then
                    do l=1,INT(np_ij(i,j))
                        if ( (k+np_j+l) <= np_total ) then
                            call random_number(eta)
                            p_rep%p_i(1,k+np_j+l) = y_i + eta*s_rep%gd%delta_x_i(1)
                        end if 
                    end do
                end if
               
                y_i = y_i + s_rep%gd%delta_x_i(1)
                np_j = np_j + np_ij(i,j)
            end do

            ! Initialize x coordinates for current block of particles
            if (np_j > 0) then
                do l=1,np_j
                    if ( (k+l) <= np_total ) then 
                        call random_number(eta)
                        p_rep%p_i(2,k+l) = x_i + eta*s_rep%gd%delta_x_i(2)
                    end if 
                end do
            end if

            k = k+np_j
            x_i = x_i + s_rep%gd%delta_x_i(2)
        end do

        write(*,*) "mc_reps:: mc_f_to_xy  - np_total generated/np_desired  = ", np_total, "/", np
        deallocate(np_ij)
        
        mc_f_to_xy = p_rep
    end function mc_f_to_xy                

    ! Constructs an s_rep structure with a distribution function
    ! initialized to the given input function f_init over given gird gd. 
    function mc_build_srep (gd, f_init)

        ! Declare IN/OUT Variables
        type(stat_rep) :: mc_build_srep, s_rep
        type(grid), intent(IN) :: gd
        procedure(func_x_i), pointer, intent(IN) :: f_init

        ! Declare Auxiliary Variables
        integer :: i, j, nx, ny, nd
        real(8) :: delta_x, delta_y 
        real(8), dimension(:), allocatable :: x_i
       
        ! Initialize s_rep to contain given grid
        s_rep%gd = mc_grid_copy( gd )

        ! Initialize local variables
        nd = s_rep%gd%nd
        nx = s_rep%gd%szs(2)
        ny = s_rep%gd%szs(1)
        delta_x = (s_rep%gd%upper_bnds(2)-s_rep%gd%lower_bnds(2))/real(s_rep%gd%szs(2)-1) ! x coordinate grid spacing
        delta_y = (s_rep%gd%upper_bnds(1)-s_rep%gd%lower_bnds(1))/real(s_rep%gd%szs(1)-1) ! y coordinate grid spacing

        ! Allocate space for local array
        allocate(x_i(gd%nd))
        allocate(s_rep%f_yx(ny,nx))

        ! Allocate space for s_rep function
        ! at each point in the given grid and storing in s_rep
        x_i(2) = s_rep%gd%lower_bnds(2)
        do j=1,nx
                x_i(1) = s_rep%gd%lower_bnds(1)
                do i=1,ny
                    s_rep%f_yx(i,j) = f_init(nd,x_i)
                    x_i(1) = x_i(1) + delta_y
                end do
                x_i(2) = x_i(2) + delta_x
        end do

        ! Normalize the distribution function over the grid
        s_rep%f_yx = s_rep%f_yx/mc_srep_norm(s_rep)

        mc_build_srep = s_rep
    end function mc_build_srep 

    ! This function returns an s_rep structure that contains a matrix
    ! with the contiuum representation of a function on a given grid gd
    ! corresponding to the particle data in the given p_rep structure.
    function mc_xy_to_f (p_rep, gd)

        ! Declare IN/OUT variables
        type(stat_rep) :: mc_xy_to_f, s_rep
        type(particle_rep), intent(IN) :: p_rep
        type(grid), intent(IN) :: gd

        ! Declare Auxiliary Variables
        integer :: i, j, nx, ny, nd, np, flag
        integer, dimension(:), allocatable :: x_i
        real(8), dimension(:,:), allocatable :: np_ij

        ! MPI vars
        integer :: rank, mpi_size, ierror 
        call MPI_Comm_size(MPI_COMM_WORLD, mpi_size, ierror)
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

        ! Check input arguments 
        if (p_rep%nd /= gd%nd ) then
            write(*,*) "mc_xy_to_f: error in input arguments"
            return
        end if

        ! Initialize variables 
        np = p_rep%np/mpi_size
        nd = gd%nd
        nx = gd%szs(2)
        ny = gd%szs(1)

        write(*,*) "Task ", rank, " = NX - ", nx, ", NY - ", ny
       
        ! Allocate space for variables 
        allocate(x_i(nd))
        allocate(s_rep%f_yx(ny,nx))
        allocate(np_ij(ny, nx))
        
        ! Initialize arrays and structures
        s_rep%gd = mc_grid_copy(gd)
        np_ij(:,:) = 0
        
        ! Count number of particles in each cell grid
        do i=1,np
            ! Get index of what cell each particle coordinate lies in
            x_i(:) = floor((p_rep%p_i(:,i)-gd%lower_bnds(:)+gd%delta_x_i(:)/2.0)/gd%delta_x_i(:)) + 1
            
            ! Check that particle is in bounds
            flag = 1
            do j=1,nd
                if ( x_i(j) <= 0 .or. x_i(j) > gd%szs(j) ) then
                    flag = 0
                    exit           
                end if
            end do                        

            ! Only add to grid if in bounds 
            if ( flag == 1 ) then 
                np_ij(x_i(1),x_i(2)) = np_ij(x_i(1),x_i(2)) + 1.0
            end if
        end do

        ! Reduce np_ij accross processes and accumulate in master node
        if (rank == 0) then
            call MPI_Reduce(MPI_IN_PLACE, np_ij, ny*nx, MPI_REAL, & 
                        MPI_SUM, 0, MPI_COMM_WORLD, ierror)
        else 
            call MPI_Reduce(np_ij, np_ij, ny*nx, MPI_REAL, & 
                        MPI_SUM, 0, MPI_COMM_WORLD, ierror)
        end if

        ! Construct f_y at (y_x) grid from particle histogram count np_ij
        if (rank == 0) then
          s_rep%f_yx(:,:) = np_ij(:,:)
          s_rep%f_yx(1,:) = 2*s_rep%f_yx(1,:)
          s_rep%f_yx(ny,:) = 2*s_rep%f_yx(ny,:)
          s_rep%f_yx(:,1) = 2*s_rep%f_yx(:,1);
          s_rep%f_yx(:,nx) = 2*s_rep%f_yx(:,nx)

          ! Normalize f
          call mc_grid_print("s_rep-grid in xy_to_f", s_rep%gd)
          s_rep%f_yx = s_rep%f_yx / mc_srep_norm(s_rep)
        end if

        ! Free allocated memory
        deallocate(np_ij)
        deallocate(x_i)

        ! Set returned value
        mc_xy_to_f = s_rep
    end function mc_xy_to_f 

    ! Scans the p_rep structure to find the maximum size grid that contains all the particles.
    ! Scales this curr_gd to this new max grid but preserving the delta_x spacing of the curr_gd.
    subroutine mc_fit_grid_to_particles(curr_gd, p_rep)
        type(grid), intent(INOUT) :: curr_gd
        type(particle_rep), intent(IN) :: p_rep 

        integer :: i,j,np,nd
        real(8), dimension(:), allocatable :: new_upper, new_lower, new_szs

        ! MPI vars
        integer :: rank, ierror 
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

        ! Check arguments
        if (curr_gd%nd /= p_rep%nd) then
            write(*,*) "error"
            return
        end if
        
        ! Initialize local variables
        np = p_rep%np
        nd = curr_gd%nd

        ! Allocate space for variables
        allocate(new_lower(nd))
        allocate(new_upper(nd))

        do i=1,nd
            new_upper(i) = curr_gd%upper_bnds(i)
            new_lower(i) = curr_gd%lower_bnds(i)
        end do

        ! Find new upper and lower bounds that fit the given p_rep for each task
        do j=1,np
            do i=1,nd
                if (p_rep%p_i(i,j) > new_upper(i)) then
                    new_upper(i) = p_rep%p_i(i,j)
                end if

                if (p_rep%p_i(i,j) < new_lower(i)) then
                    new_lower(i) = p_rep%p_i(i,j)
                end if
            end do
        end do

        ! Round up to multiples of grid spacing to perserve it
        do i=1,nd
          ! Round upper bnd up 
          new_upper(i) = ROUND_UP(new_upper(i),curr_gd%delta_x_i(i))
                                   
          ! Round lower bnd down
          new_lower(i) = ROUND_DOWN(new_lower(i),curr_gd%delta_x_i(i))
        end do

        call MPI_ALLreduce(new_lower, curr_gd%lower_bnds, nd, MPI_REAL, &
                                MPI_MIN, MPI_COMM_WORLD, ierror);
        call MPI_ALLreduce(new_upper, curr_gd%upper_bnds, nd, MPI_REAL, &
                                MPI_MAX, MPI_COMM_WORLD, ierror);

        ! Compute new grid spacing size
        do i=1,nd
          curr_gd%szs(i) = int(abs(curr_gd%upper_bnds(i) - &
                                  curr_gd%lower_bnds(i))/curr_gd%delta_x_i(i)) 
        end do
        
        ! Free allocated memory
        deallocate(new_upper)
        deallocate(new_lower)
    end subroutine mc_fit_grid_to_particles

end module mc_reps

