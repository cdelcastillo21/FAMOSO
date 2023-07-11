! Carlos del-Castillo-Negrete 
! test_mc - Program tests the standard Monte Carlo solver
! routine from the FAMOSO library. Data of simulation outputted
! to file in ../data/ director.
program test_mc
    use mpi
    use montecarlo
    use test_functions

    implicit none

    ! Parameters of system
    integer, parameter  :: nd = 2
    real(8), parameter :: delta_t = 0.005, t_init = 0.0, t_final = 32.0 
    real(8), dimension(nd), parameter :: delta_x_i = (/0.1, 0.1/)
    real(8), dimension(nd), parameter :: lower_bnds = (/ 0.0, -3.0/)
    real(8), dimension(nd), parameter :: upper_bnds = (/ 4.0, 3.0/)

    ! Variables of problem
    type(mc_problem) :: mc_prob
    type(grid) :: init_gd
    type(stat_rep) :: s_rep_final 

    ! Procedure pointers
    procedure(func_x_i), pointer :: f_init=>f_init_2d_double_hump
    procedure(stoch_diff_dx), pointer :: f_dv=>dV_vpar_vperp

    ! Input Arguments - Num particles and Num snapshots
    integer :: num_args, ix, np, n_dump
    character(len=20), dimension(:), allocatable :: args
    character(:), allocatable :: run_name

    num_args = command_argument_count()
    allocate(args(num_args))
    do ix = 1, num_args
      call get_command_argument(ix,args(ix))
      write(*,*) args(ix)
    end do

    read(args(1),"(I)") np
    read(args(2),"(I)") n_dump

    ! Construct initial grid
    init_gd = mc_grid_build(nd, delta_x_i, lower_bnds, upper_bnds)

    ! Construct montecarlo problem structure within defined parameters, functions and initial grid
    mc_prob = mc_prob_construct(t_init, t_final, delta_t, init_gd, f_init, f_dv)
    
    write(*,*) LEN(args(3)) 
    run_name = TRIM(args(3))
    write(*,*) LEN(run_name) 
    call mc_solve_and_record(mc_prob, np, n_dump, run_name)
end program test_mc


