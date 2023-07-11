! Carlos del-Castillo-Negrete, April, 2019
! MONTECARLO module - Standard Monte Carlo Algorithm
module montecarlo
    use mpi
    use utilities
    use mc_grid
    use mc_reps
    implicit none

    type mc_sim_data
            integer :: N
            type(stat_rep), dimension(:), allocatable :: s_rep_i 
            real(8), dimension(:), allocatable :: t_i
            real(8), dimension(:), allocatable :: cpu_t_i
    end type mc_sim_data 

    type mc_problem
            type(grid) :: gd
            procedure(func_x_i), pointer, nopass :: f_init 
            procedure(stoch_diff_dx), pointer, nopass :: dX
            real(8) :: t_i, t_f, delta_t
            type(stat_rep) :: s_rep_init
            type(mc_sim_data) :: mc_data 
    end type mc_problem

    ! Define abstract format of functions used in the montecarlo module 
    abstract interface
        ! stoch_diff_dx defines a function that returns the dx incremental
        ! value of a given particle with the nd coordinates definend in x_i over
        ! the time step delta_t. d_eta is an array of nd uniformily distributed
        ! random numbers in the intervale [0,1] that this function uses as the
        ! random pertubations of the stochastic process modeled. 
        function stoch_diff_dx (nd, delta_t, x_i, d_eta)
                integer, intent(IN) :: nd
                real(8), dimension(nd) :: stoch_diff_dx
                real(8), intent(IN) :: delta_t
                real(8), dimension(nd), intent(IN) :: x_i, d_eta
        end function stoch_diff_dx
    end interface

    ! Declare Parameters of module 
    integer, parameter :: max_dim = 2

contains

    function mc_prob_construct(t_i, t_f, delta_t, gd, f_init, dX)
        type(mc_problem) :: mc_prob_construct, mc_prob
        type(grid), intent(IN) :: gd 
        real(8), intent(IN) :: t_i, t_f, delta_t
        procedure(func_x_i), pointer, intent(IN) :: f_init 
        procedure(stoch_diff_dx), pointer, intent(IN) :: dX
       
        mc_prob%gd = mc_grid_copy( gd )
        mc_prob%t_i = t_i
        mc_prob%t_f = t_f
        mc_prob%delta_t = delta_t
        mc_prob%f_init => f_init
        mc_prob%dX => dX

        mc_prob_construct = mc_prob
    end function mc_prob_construct

    ! mc_step iterates the montecarlo particle model for nt time
    ! time steps of size delta_t according to the stochastic 
    ! differential equation f_dv.
    subroutine mc_step (p_rep, nt, delta_t, f_dv, max_gd)
        ! Declare IN/OUT variables
        type(particle_rep), intent(INOUT) :: p_rep
        integer, intent(IN) :: nt                               ! # time 
        real(8), intent(IN) :: delta_t                          ! time step size
        procedure(stoch_diff_dx), pointer, intent(IN) :: f_dv    ! external subroutine that evaluates stochastic diff eq.
        type(grid), optional, intent(INOUT) :: max_gd
        
        ! Declare variables for microscopic integration
        integer :: i, j, nd, np 
        real(8), dimension(:), allocatable :: dv        ! 
        real(8), dimension(:), allocatable :: d_eta     ! Random variables

        ! Initialize local variables
        nd = p_rep%nd
        np = p_rep%np

        ! Allocate space for local variables
        allocate (dv(nd))              
        allocate (d_eta(nd))
       
        ! Initiate random number generator
        call init_random_seed()
        
        ! Perform microscopic integration
        do i = 1,nt
            ! For each particle calculate new position at next time step
            do j=1,np
                ! Generate random numbers for f_dv function to use
                call random_number(d_eta)

                ! Calculate increments dv to each coordinate x_i by calling 
                ! the passed external function
                dv = f_dv(nd, delta_t, p_rep%p_i(:,j), d_eta)

                ! Add increment to current position
                p_rep%p_i(:,j) = p_rep%p_i(:,j) + dv(:)
            end do
        end do                
      
        ! if max_gd argument is specified, then compute new max gd
        ! Note - this will synchronize all processes via Allreduce call
        if (present(max_gd)) then
            call mc_fit_grid_to_particles(max_gd,p_rep)                        
        end if

        ! Free memory 
        deallocate (dv)              
        deallocate (d_eta)
    end subroutine mc_step 
  
    function mc_solve(mc_prob, np)
        ! Input/Ouput variables
        type(stat_rep) :: mc_solve
        type(mc_problem), intent(IN) :: mc_prob
        integer, intent(IN) :: np

        ! Variables of problem
        integer :: nt, np_i
        type(grid) :: max_gd
        type(particle_rep) :: p_rep
        type(stat_rep) :: s_rep_init, s_rep_final 

        ! MPI vars
        integer :: rank, mpi_size, ierror 

        ! Initialize MPI
        call MPI_Init(ierror)
        call MPI_Comm_size(MPI_COMM_WORLD, mpi_size, ierror)
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

        ! Divide particles amongst all processes
        np_i = np / mpi_size

        ! Max grid that contains all particles after iterating montecarlo module
        max_gd = mc_grid_copy( mc_prob%gd )

        call mc_grid_print( "curr_grid", mc_prob%gd )

        ! Calculate number of time steps to iterate model for 
        nt = nint((mc_prob%t_f-mc_prob%t_i)/mc_prob%delta_t)

        ! Build s_rep and initialize to init distribution function
        s_rep_init = mc_build_srep(mc_prob%gd, mc_prob%f_init)

        ! Construct particle representation from given initial s_rep 
        p_rep = mc_f_to_xy (s_rep_init, np_i)

        ! Iterate microscopic integration to solve system for given number of timesteps
        call mc_step(p_rep, nt, mc_prob%delta_t, mc_prob%dX, max_gd)

        ! Construct final distribution function
        s_rep_final = mc_xy_to_f(p_rep, max_gd)

        ! Free used structures
        call mc_destroy_prep(p_rep)
        call mc_grid_destroy(max_gd)
        call mc_destroy_srep(s_rep_init)

        call MPI_Finalize(ierror)

        mc_solve = s_rep_final
    end function mc_solve

    subroutine  mc_solve_and_record(mc_prob, np, N, out_data_file)
        ! Declare IN/OUT variables
        type(mc_problem), intent(INOUT) :: mc_prob
        integer, intent(IN) :: np, N
        character(*), optional, intent(INOUT) :: out_data_file

        ! Variables method 
        integer :: i,nt, np_n
        real :: t_1, t_2
        type(mc_sim_data) :: mc_data
        type(particle_rep) :: p_rep

        ! MPI vars
        integer :: rank, mpi_size, ierror 

        ! Initialize MPI
        call MPI_Init(ierror)
        call MPI_Comm_size(MPI_COMM_WORLD, mpi_size, ierror)
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

        ! Calculate number of particles to processed by this node
        np_n = np / mpi_size

        ! Allocate space for local structures and arrays
        allocate(mc_data%s_rep_i(N+1))
        allocate(mc_data%t_i(N+1))
        allocate(mc_data%cpu_t_i(N+1))

        ! calculate descrete number of time steps
        nt = nint((mc_prob%t_f-mc_prob%t_i)/(mc_prob%delta_t*N))
        mc_data%N = N+1

        ! Build initial s_rep distributuion. If initial s_rep provided
        ! start from that. If not resort to f_init fuction
        if ( allocated(mc_prob%s_rep_init%f_yx) ) then
            mc_data%s_rep_i(1) = mc_copy_srep(mc_prob%s_rep_init)
        else
            mc_data%s_rep_i(1) = mc_build_srep(mc_prob%gd, mc_prob%f_init)
        end if  
        mc_data%t_i(1) = mc_prob%t_i
        mc_data%cpu_t_i(1) = 0.0

        ! Construct particle representation from given initial s_rep 
        p_rep = mc_f_to_xy(mc_data%s_rep_i(1), np_n)

        do i=1,N
            ! Iterate microscopic integration to solve system for given number of timesteps
            t_1 = MPI_Wtime()
            call mc_step(p_rep, nt, mc_prob%delta_t, mc_prob%dX, mc_prob%gd)
            t_2 = MPI_Wtime()

            ! Construct and store reconstructed distributuion function 
            ! Note - Only master process gets s_rep snapshots
	    if(rank == 0) then
	      call mc_grid_print("max grid in solve_and_record", mc_prob%gd)
	    end if 
            mc_data%s_rep_i(i+1) = mc_xy_to_f(p_rep, mc_prob%gd)
            mc_data%t_i(i+1) = mc_data%t_i(i)+nt*mc_prob%delta_t 
            mc_data%cpu_t_i(i+1) = t_2-t_1
        end do

        ! Master - Process -> Store MC data by aggregating to past data if necessary.
        ! If not in master process, destroy mc_prob data structure as it is not needded.
        if (rank == 0) then
          if (allocated(mc_prob%mc_data%s_rep_i)) then
              call mc_merge_data(mc_prob%mc_data, mc_data)
          else
              mc_prob%mc_data = mc_copy_data( mc_data )
              call mc_destroy_sim_data(mc_data)
          end if

          if (present(out_data_file)) then
            call mc_data_output(out_data_file, mc_prob)
          end if
        else
          call mc_destroy_sim_data(mc_data)
        end if

        ! Destroy particle rep (only used in this method)
        call mc_destroy_prep(p_rep)  

        call MPI_Finalize(ierror)
    end subroutine mc_solve_and_record

    ! Returns an allocated copy of the data stored in the mc_data structure.
    function mc_copy_data(mc_data)
        type(mc_sim_data) :: mc_copy_data, mc_copy
        type(mc_sim_data), intent(IN) :: mc_data

        integer :: i
        
        mc_copy%N = mc_data%N 
        
        allocate(mc_copy%s_rep_i(mc_data%N))
        allocate(mc_copy%t_i(mc_data%N))
        allocate(mc_copy%cpu_t_i(mc_data%N))

        do i=1,mc_data%N
            mc_copy%t_i(i) = mc_data%t_i(i)
            mc_copy%cpu_t_i(i) = mc_data%cpu_t_i(i)
            mc_copy%s_rep_i(i) = mc_copy_srep(mc_data%s_rep_i(i))
        end do

        mc_copy_data = mc_copy
    end function mc_copy_data

    ! Merges data from mc_data_1 and mc_data_2 and returns merged data.
    ! NOTE: this routine assumes mc_data_1 contains s_reps at times t_i < mc_data_2. 
    subroutine mc_merge_data(mc_data_1, mc_data_2)
        type(mc_sim_data), intent(INOUT) :: mc_data_1, mc_data_2 
        type(mc_sim_data) :: merged_data

        merged_data%N = mc_data_1%N+mc_data_2%N
        
        allocate(merged_data%s_rep_i(merged_data%N))
        allocate(merged_data%t_i(merged_data%N))
        allocate(merged_data%cpu_t_i(merged_data%N))
       
        merged_data%t_i(1:mc_data_1%N) = mc_data_1%t_i(1:mc_data_1%N)  
        merged_data%t_i((mc_data_1%N+1):merged_data%N) = mc_data_2%t_i(1:mc_data_2%N)  
        merged_data%cpu_t_i(1:mc_data_1%N) = mc_data_1%cpu_t_i(1:mc_data_1%N)  
        merged_data%cpu_t_i((mc_data_1%N+1):merged_data%N) = mc_data_2%cpu_t_i(1:mc_data_2%N)  
        merged_data%s_rep_i(1:mc_data_1%N) = mc_data_1%s_rep_i(1:mc_data_1%N)  
        merged_data%s_rep_i((mc_data_1%N+1):merged_data%N) = mc_data_2%s_rep_i(1:mc_data_2%N)  

        ! deallocate old structures
        call mc_destroy_sim_data(mc_data_1)
        call mc_destroy_sim_data(mc_data_2)

        mc_data_1 = mc_copy_data(merged_data)

        call mc_destroy_sim_data(merged_data)
    end subroutine mc_merge_data

    ! Prints the contents of an mc_data structure
    subroutine mc_data_print(mc_data)
        type(mc_sim_data) :: mc_data

        integer :: i

        write(*,*) "-------------------------"
        write(*,*) "Printing mc_data - ", mc_data%n, " data points"
        do i=1,mc_data%n
            write(*,*) "data point at t = ", mc_data%t_i(i)
            write(*,*) "Computation time = ", mc_data%cpu_t_i(i)
            call mc_srep_print( mc_data%s_rep_i(i) )
        end do
    end subroutine mc_data_print

    ! Outputs data stroed in mc_data structure to file with the name stored in desc
    ! Can load data into matlab/octave with the script load_mc_data.m
    subroutine mc_data_output(desc, mc_prob)
        character(*), intent(IN) :: desc
        type(mc_problem), intent(INOUT) :: mc_prob 

        integer :: i,j,nx,ny,N,k
        type(grid) :: max_gd
        real(8), dimension(:,:,:), allocatable :: A 
      
        character(8) :: ext = ".mc_data"
        character(100) :: file_name

        ! Open file with given name
        file_name=trim(desc)//ext
        open(unit=1,file=trim(adjustl(desc//ext)),blank='ZERO')

        ! Make sure all the s_reps in the mc_data structure are of same size
        call mc_reps_scale_all(mc_prob%mc_data%s_rep_i, max_gd, mc_prob%mc_data%N)

        N = mc_prob%mc_data%N
        ny = max_gd%szs(1) 
        nx = max_gd%szs(2)

        write(*,*) "Sizes - ", N, ny, nx

        ! Write data to file
        write(1,*) mc_prob%delta_t
        write(1,*) max_gd%szs
        write(1,*) max_gd%lower_bnds
        write(1,*) max_gd%upper_bnds
        write(1,*) N
        write(1,*) mc_prob%mc_data%cpu_t_i
        write(1,*) mc_prob%mc_data%t_i

        do k=1,N
            do j=1,nx
                do i=1,ny
                    write(1,*)  mc_prob%mc_data%s_rep_i(k)%f_yx(i,j)
                end do
            end do
        end do

        close(1)
    end subroutine mc_data_output

    ! Destroy allocated memory associated with an mc_data structure
    subroutine mc_destroy_sim_data(mc_data)
            type(mc_sim_data) :: mc_data
            integer :: i

            ! Free allocated memory associated with p_rep
            do i=1,mc_data%N
                    call mc_destroy_srep(mc_data%s_rep_i(i))
            end do
            deallocate(mc_data%s_rep_i)
            deallocate(mc_data%t_i)
            deallocate(mc_data%cpu_t_i)
    end subroutine mc_destroy_sim_data
end module montecarlo

