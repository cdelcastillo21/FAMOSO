! Carlos del-Castillo-Negrete
! April, 2015
! test_pjmc - Program to test PMC algorithm implemented in FAMOSO library.
! Simulation data outputed to a given file in the ../data/ director. 
program test_pjmc

    use projective_montecarlo
    use mc_grid
    use test_functions

    implicit none

    ! Parameters of system
    integer, parameter :: r1 = 5, r2 = 5
    integer, parameter  :: np=100000, nd=2, num_proj_iter = 1
    real(8), parameter :: delta_t = 0.005, t_init = 0.0, t_final = 8
    real(8), parameter :: pj_ratio = 0.75, pj_aff = 0.2, e_thresh = 0.9
    real(8), dimension(nd), parameter :: delta_x_i = (/0.1, 0.1/)
    real(8), dimension(nd), parameter :: lower_bnds = (/ 0.0, -3.0/)
    real(8), dimension(nd), parameter :: upper_bnds = (/ 4.0, 3.0/)

    ! Variables of problem
    type(pjmc_problem) :: pjmc_prob
    type(grid) :: init_gd

    ! Procedure pointers
    procedure(func_x_i), pointer :: f_init=>f_init_2d_double_hump
    procedure(stoch_diff_dx), pointer :: f_dv=>dV_vpar_vperp

    ! Construct initial grid
    init_gd = mc_grid_build(nd, delta_x_i, lower_bnds, upper_bnds)

    ! Construct montecarlo problem structure within defined parameters, functions and initial grid
    pjmc_prob = pjmc_prob_construct(t_init, t_final, delta_t, num_proj_iter, init_gd, f_init, f_dv, pj_ratio, pj_aff, e_thresh)
    
    ! solve the system and get the final distribution function
    call pjmc_solve(pjmc_prob, np, r1, r2, init_gd)
   
    ! output data 
    call pjmc_data_output("../data/pjmc_test", pjmc_prob)
     
end program test_pjmc
