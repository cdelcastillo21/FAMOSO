! Carlos del-Castillo-Negrete, April, 2015
! FAMOSO - FAst MOnte Carlo SOLver - PROJECTIVE MONTECARLO module
! Module implementing functions corresponding to projective monte carlo solver routines.
module projective_montecarlo

    use montecarlo
    use glram

    type :: pjmc_problem
        type(grid) :: gd
        procedure(func_x_i), pointer, nopass :: f_init 
        procedure(stoch_diff_dx), pointer, nopass :: dX
        real(8) :: t_i, t_f, delta_t
        type(stat_rep) :: s_rep_init
        type(mc_sim_data) :: mc_data 
        integer :: num_iter
        real(8) :: pj_ratio 
        real(8) :: pj_aff
        real(8) :: e_thresh
        type(glram_rep) :: g_rep
    end type pjmc_problem

    integer :: DEBUG = 1 
contains
    
    function pjmc_prob_construct(t_i, t_f, delta_t, num_iter, gd, f_init, dX, pj_ratio, pj_aff, e_thresh)
        type(pjmc_problem) :: pjmc_prob_construct, pjmc_prob
        type(grid), intent(IN) :: gd 
        integer, intent(IN) :: num_iter
        real(8), intent(IN) :: t_i, t_f, delta_t, pj_ratio, pj_aff, e_thresh
        procedure(func_x_i), pointer, intent(IN) :: f_init 
        procedure(stoch_diff_dx), pointer, intent(IN) :: dX
       
        pjmc_prob%gd = mc_grid_copy( gd )
        pjmc_prob%t_i = t_i
        pjmc_prob%t_f = t_f
        pjmc_prob%delta_t = delta_t
        pjmc_prob%pj_ratio = pj_ratio
        pjmc_prob%pj_aff = pj_aff
        pjmc_prob%e_thresh = e_thresh
        pjmc_prob%num_iter = num_iter
        pjmc_prob%f_init => f_init
        pjmc_prob%dX => dX

        pjmc_prob_construct = pjmc_prob
    end function pjmc_prob_construct

    subroutine pjmc_solve(pjmc_prob, np, r1, r2, max_gd)
        type(pjmc_problem), intent(INOUT) :: pjmc_prob
        integer, intent(IN) :: np, r1, r2
        type(grid), optional, intent(IN) :: max_gd

        integer :: i, n_data, timer
        real(8) :: time_interval, t_i_mic, t_f_mic, t_proj, cpu_proj_time

        ! Calculate number of time steps to iterate regular montecarlo
        ! integration for. 
        time_interval = (pjmc_prob%t_f-pjmc_prob%t_i)/real(pjmc_prob%num_iter)
        n_data = int(1.0/pjmc_prob%pj_aff)
        t_i_mic = pjmc_prob%t_i
        t_f_mic = pjmc_prob%t_i + time_interval*(1.0-pjmc_prob%pj_ratio) 
        t_proj = pjmc_prob%t_i + time_interval
        cpu_proj_time = 0.0

        ! If diagnostic variable set output parameters of projective solver:
        if ( DEBUG == 1 ) then
            write(*,*) " ------ Starting Projective Monte Carlo Solver ----- "
            write(*,*) "PARAMETERS:"
            write(*,*) "-----------"
            write(*,*) "t_i = ", pjmc_prob%t_i, ", t_f = ", pjmc_prob%t_f
            write(*,*) "# of projective iterations = ", pjmc_prob%num_iter
            write(*,*) "projection step time interval = ", time_interval
            write(*,*) "projection ratio = ", pjmc_prob%pj_ratio, "projection aff = ", pjmc_prob%pj_aff
            write(*,*) "t_i_mic = ", t_i_mic, ", t_f_mic = ", t_f_mic, " t_proj = ", t_proj
            if (present(max_gd)) then
                call mc_grid_print("USING MAX GRID:", max_gd)
            end if
        end if

        do i=1,pjmc_prob%num_iter
            ! Integrate microscopically. Record snapshots for projection
            pjmc_prob%t_i = t_i_mic
            pjmc_prob%t_f = t_f_mic
       
            call mc_grid_print("curr_grid:", pjmc_prob%gd) 
            write(*,*) "microscopic integration from ", t_i_mic, " to ", t_f_mic
            write(*,*) "recording ", n_data, " snapshots"
            call start_timer(timer)
            if (present(max_gd)) then
                call pjmc_microscopic_integrate(pjmc_prob, np, n_data, cpu_proj_time, max_gd)
            else
                call pjmc_microscopic_integrate(pjmc_prob, np, n_data, cpu_proj_time)
            end if
            write(*,*) "..........done with mic_integration - ", get_time(timer), " seconds"
  
            ! Construct projection and set as initialization s_rep for next iteration
            if ( allocated(pjmc_prob%s_rep_init%f_yx) ) then
                call mc_destroy_srep(pjmc_prob%s_rep_init)
            end if

            ! Check if projection step necessary
            if (t_proj > t_f_mic) then
                write(*,*) "**** Projecting from ", t_f_mic, " to ", t_proj
                write(*,*) "**** Using rank ", r1 , " by ", r2, " approximation "
                call start_timer(timer)
                pjmc_prob%s_rep_init = pjmc_project_glram(pjmc_prob%mc_data,t_proj,n_data,r1,r2,pjmc_prob%e_thresh)
                cpu_proj_time = get_time(timer)
                write(*,*) "..... DONE in ", cpu_proj_time, " seconds"
            else
                pjmc_prob%s_rep_init = mc_copy_srep(pjmc_prob%mc_data%s_rep_i(pjmc_prob%mc_data%N))
            end if

            ! Set new time interval for next iteration
            t_i_mic = t_i_mic + time_interval
            t_f_mic = t_f_mic + time_interval
            t_proj = t_proj + time_interval 
        end do 
       
    end subroutine pjmc_solve
    
    subroutine  pjmc_microscopic_integrate(pjmc_prob, np, N, cpu_t_0, max_gd)
        ! Declare IN/OUT variables
        type(pjmc_problem), intent(INOUT) :: pjmc_prob
        integer, intent(IN) :: np, N
        real(8), intent(IN) :: cpu_t_0
        type(grid), optional, intent(IN) :: max_gd

        ! Variables method 
        integer :: i,nt,timer
        type(mc_sim_data) :: mc_data
        type(particle_rep) :: p_rep
        type(stat_rep) :: s_rep_recon

        ! Allocate space for local structures and arrays
        allocate(mc_data%s_rep_i(N+1))
        allocate(mc_data%t_i(N+1))
        allocate(mc_data%cpu_t_i(N+1))

        ! calculate descrete number of time steps
        nt = nint((pjmc_prob%t_f-pjmc_prob%t_i)/(pjmc_prob%delta_t*N))
        mc_data%N = N+1

        ! Build initial s_rep distributuion. If initial s_rep provided
        ! start from that. If not resort to f_init fuction
        if ( allocated(pjmc_prob%s_rep_init%f_yx) ) then
            write(*,*) "    ***Restarting Microscopic Integration***"
            write(*,*) "    initializing microscopic integration from previous dist_function"
            mc_data%s_rep_i(1) = mc_copy_srep(pjmc_prob%s_rep_init)
        else
            mc_data%s_rep_i(1) = mc_build_srep(pjmc_prob%gd, pjmc_prob%f_init)
            write(*,*) "    initializing microscopic integration from init dist_function"
        end if  
        mc_data%t_i(1) = pjmc_prob%t_i
        mc_data%cpu_t_i(1) = cpu_t_0

        ! Construct particle representation from given initial s_rep 
        p_rep = mc_f_to_xy(mc_data%s_rep_i(1), np)
        
        do i=1,N
            write(*,*) "-------- starting microscopic stepper iteration ", i, "------"
            ! Iterate microscopic integration to solve system for given number of timesteps
            call start_timer(timer)
            call mc_step(p_rep, nt, pjmc_prob%delta_t, pjmc_prob%dX)

            ! Update grid of problem to contain all the particles  
            if (present(max_gd)) then
                call mc_fit_grid_to_particles(pjmc_prob%gd, p_rep, max_gd)
            else
                call mc_fit_grid_to_particles(pjmc_prob%gd, p_rep)
            end if

            ! construct and store reconstructed distributuion function 
            mc_data%s_rep_i(i+1) = mc_xy_to_f(p_rep, pjmc_prob%gd)
            mc_data%t_i(i+1) = mc_data%t_i(i)+nt*pjmc_prob%delta_t 
            mc_data%cpu_t_i(i+1) = get_time(timer)
            
            write(*,*) "...... DONE ..... ", mc_data%cpu_t_i(i+1), " seconds."
        end do

        if (allocated(pjmc_prob%mc_data%s_rep_i)) then
            call mc_merge_data(pjmc_prob%mc_data, mc_data)
        else
            pjmc_prob%mc_data = mc_copy_data( mc_data )
            call mc_destroy_sim_data(mc_data)
        end if

        ! Destroy particle rep (only used in this method)
        call mc_destroy_prep(p_rep)  

    end subroutine pjmc_microscopic_integrate

    function pjmc_project_glram( mc_data, t_project, N, l1, l2, e_thresh)
        type(stat_rep) :: pjmc_project_glram, s_rep_proj
        type(mc_sim_data), intent(INOUT) :: mc_data
        real(8), intent(IN) :: t_project
        integer(4), intent(IN) :: N, l1, l2
        real(8), intent(IN) :: e_thresh

        integer :: i, j, N_tot
        real(8), dimension(:), allocatable :: t_i
        real(8), dimension(:,:), allocatable :: M_proj, temp
        type(glram_rep) :: g_rep

        integer, parameter :: n_for_projection = 2

        ! Allocate space for local variables
        allocate(M_proj(l1,l2))
        allocate(t_i(N))

        ! Compute glram compressed rep of n_data most recent data points
        ! stored in the mc_data structure. Use r1*r2 coefficient matrices
        g_rep = pjmc_compute_glram(mc_data, l1, l2, N)

        ! Initialize variablesV
        N_tot = mc_data%N
        t_i = mc_data%t_i((N_tot-N+1):N_tot)
        s_rep_proj%gd = mc_grid_copy(mc_data%s_rep_i(mc_data%N)%gd)

        ! Now repeat but filtering M_i using e_thresh
        do i=1,g_rep%N
            call mc_e_thresh(g_rep%M_i(:,:,i), l1, l2, e_thresh)
        end do

        ! Project each coefficient according to linear model
        do j=1,l2
            do i=1,l1
                M_proj(i,j) = linear_project_simple( N, t_i, g_rep%M_i(i,j,:), t_project)
            end do
        end do

        ! Reconstruct projected matrix
        s_rep_proj%f_yx = glram_reconstruct(M_proj, g_rep)

        ! Find minimum entry and add to function -> don't want negative #s
        do j=1,g_rep%c
            do i=1,g_rep%r
                if (s_rep_proj%f_yx(i,j) < 0.0) then
                    s_rep_proj%f_yx(i,j) = 0.0
                end if
            end do
        end do

        ! Normalize
        s_rep_proj%f_yx(:,:) = s_rep_proj%f_yx(:,:)/mc_srep_norm(s_rep_proj)
        
        deallocate(M_proj)
        deallocate(t_i)
        call glram_destroy(g_rep)

        pjmc_project_glram = s_rep_proj 
    end function pjmc_project_glram

    function pjmc_compute_glram(mc_data, r1, r2, N)
        type(glram_rep) :: pjmc_compute_glram, g_rep
        type(mc_sim_data), intent(INOUT) :: mc_data
        integer, intent(IN) :: r1, r2 
        integer, intent(IN) :: N 

        integer :: i, j, start_idx, r, c, N_tot
        type(grid) :: max_gd
        real(8), dimension(:,:,:), allocatable :: A_i

        ! Make sure all s_rep_i's are defined over same grid
        ! Also get back max_gd over which all s_rep_i's are defined
        call mc_reps_scale_all(mc_data%s_rep_i, max_gd, mc_data%N)

        ! Initialize variables
        r = max_gd%szs(1)
        c = max_gd%szs(2)
        N_tot = mc_data%N
        start_idx = N_tot - N + 1

        ! Allocate space for temp matrices to use to compute glram compression
        allocate(A_i(r,c,N))

        ! Copy all s_reps distribution functions into matrices
        j = 1 
        do i=start_idx,N_tot
                A_i(:,:,j) = mc_data%s_rep_i(i)%f_yx(:,:)
                j = j + 1
        end do  

        ! Call external function from glram module to run algorithm
        g_rep = glram_construct( A_i, r, c, N, r1, r2)

        ! Free allocated memory for arrays and structures
        deallocate(A_i)
        call mc_grid_destroy(max_gd)

        ! Set return value 
        pjmc_compute_glram = g_rep
    end function pjmc_compute_glram

    ! Linear projection
    ! m = (N*sum(x_i*y_i) - sum(x_i)*sum(y_i)) / (N*sum(x_i^2) - (sum(x_i))^2)
    ! b = (sum(y_i)*sum(x_i^2) - sum(x_i)*sum(x_i*y_i))/(n*sum(x_i^2)-(sum(x_i))^2)
    ! y_proj = m*x_i + b        
    function linear_project_simple(n, x_i, y_i, x_proj)
        real(8) :: linear_project_simple, y_proj
        integer, intent(IN) :: n
        real(8), intent(IN) :: x_proj
        real(8), dimension(n), intent(IN) :: x_i, y_i

        real(8) :: slope 

        ! Calculate slope between last two points
        slope = (y_i(n) - y_i(n-1))/(x_i(n)-x_i(n-1))

        ! Calculate projected value
        y_proj = y_i(n) + (x_proj-x_i(n))*slope

        linear_project_simple = y_proj
    end function linear_project_simple

    ! Linear projection
    ! m = (N*sum(x_i*y_i) - sum(x_i)*sum(y_i)) / (N*sum(x_i^2) - (sum(x_i))^2)
    ! b = (sum(y_i)*sum(x_i^2) - sum(x_i)*sum(x_i*y_i))/(n*sum(x_i^2)-(sum(x_i))^2)
    ! y_proj = m*x_i + b        
    function linear_project(n, x_i, y_i, x_proj)
        real(8) :: linear_project, y_proj
        integer, intent(IN) :: n
        real(8), intent(IN) :: x_proj
        real(8), dimension(n), intent(IN) :: x_i, y_i

        integer :: i
        real(8) :: m, b
        real(8) :: sum_xi, sum_yi, sum_xiyi, sum_xi_2, sum_yi_2

        sum_xi = 0.0
        sum_yi = 0.0
        sum_xiyi = 0.0
        sum_xi_2 = 0.0
        sum_yi_2 = 0.0

        ! calculate values for linear regression formula
        do i=1,n
            sum_xi = sum_xi + x_i(i)
            sum_yi = sum_yi + y_i(i)
            sum_xiyi = sum_xiyi + x_i(i)*y_i(i)
            sum_xi_2 = sum_xi_2 + x_i(i)*x_i(i)
            sum_yi_2 = sum_yi_2 + y_i(i)*y_i(i)
        end do 

        ! calculate linear coefficients
        m = (n*sum_xiyi - sum_xi*sum_yi)/(n*sum_xi_2 - sum_xi*sum_xi)
        b = (sum_yi*sum_xi_2 - sum_xi*sum_xiyi)/(n*sum_xi_2 - sum_xi*sum_xi)

        ! Calculate projected value
        y_proj = m*x_proj + b

        linear_project = y_proj
    end function linear_project

    subroutine pjmc_data_output(desc, pjmc_prob)
        character(*), intent(IN) :: desc
        type(pjmc_problem), intent(INOUT) :: pjmc_prob 

        integer :: i,j,nx,ny,N,k
        type(grid) :: max_gd
        real(8), dimension(:,:,:), allocatable :: A 
      
        character(10) :: ext = ".pjmc_data"
        character(100) :: file_name

        ! Only output data if it exists
        if (allocated(pjmc_prob%mc_data%s_rep_i)) then

            ! Open file with given name and extension
            file_name=desc//ext
            open(unit=1,file=trim(adjustl(desc//ext)))

            ! Make sure all the s_reps in the mc_data structure are of same size
            call mc_reps_scale_all(pjmc_prob%mc_data%s_rep_i, max_gd, pjmc_prob%mc_data%N)

            N = pjmc_prob%mc_data%N
            ny = max_gd%szs(1) 
            nx = max_gd%szs(2)

            allocate(A(ny,nx,N))

            do k=1,N
                do j=1,nx
                    do i=1,ny
                        A(i,j,k) = pjmc_prob%mc_data%s_rep_i(k)%f_yx(i,j)
                    end do
                end do
            end do
            
            ! Write data to file
            write(1,*) pjmc_prob%delta_t, pjmc_prob%num_iter, pjmc_prob%pj_ratio, pjmc_prob%pj_aff, &
                !&pjmc_prob%g_rep%N pjmc_prob%g_rep%l1, pjmc_prob%g_rep%l2, mc_data%t_i, &
                &N, max_gd%szs, max_gd%lower_bnds, max_gd%upper_bnds, pjmc_prob%mc_data%t_i, pjmc_prob%mc_data%cpu_t_i, A 

            deallocate(A)

            close(1)
        end if

    end subroutine pjmc_data_output
    
end module projective_montecarlo
