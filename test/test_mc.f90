! Carlos del-Castillo-Negrete 
! April, 2015
! test_mc - Program tests the standard Monte Carlo solver
! routine from the FAMOSO library. Data of simulation outputted
! to file in ../data/ director.
program test_mc
    use montecarlo
    use test_functions

    implicit none

    ! Parameters of system
    integer, parameter  :: np=1000, nd=2, n_dump = 16
    real(8), parameter :: delta_t = 0.005, t_init = 0.0, t_final = 32.0 
    real(8), dimension(nd), parameter :: delta_x_i = (/0.1, 0.1/)
    real(8), dimension(nd), parameter :: lower_bnds = (/ 0.0, -3.0/)
    real(8), dimension(nd), parameter :: upper_bnds = (/ 4.0, 3.0/)

    ! Variables of problem
    type(mc_problem) :: mc_prob
    type(grid) :: init_gd

    ! Procedure pointers
    procedure(func_x_i), pointer :: f_init=>f_init_2d_double_hump
    procedure(stoch_diff_dx), pointer :: f_dv=>dV_vpar_vperp

    ! Construct initial grid
    init_gd = mc_grid_build(nd, delta_x_i, lower_bnds, upper_bnds)

    ! Construct montecarlo problem structure within defined parameters, functions and initial grid
    mc_prob = mc_prob_construct(t_init, t_final, delta_t, init_gd, f_init, f_dv)
    
    ! solve the system storing n_dump snapshots of distribution function along the way.  
    call mc_solve_and_record(mc_prob, np, n_dump)
    
    call mc_data_output("../data/mc_test", mc_prob)
end program test_mc


