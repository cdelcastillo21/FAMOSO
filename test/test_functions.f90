! Carlos del-Castillo-Negrete
! April, 2015
! test_functions - Modulue containing an example iniital distribution function
! and an example stochastic differential equation function that can be used by the 
! Monte Carlo solver routines in FAMOSO. These particular functions correspond to a
! plasma physics particle simulation using Fokker-Planck's equations. 
module test_functions

    real(8), parameter :: v_perp_min = 0.001

contains

    ! Initial distribution function - Initial state of system
    function f_init_2d_double_hump(nd, x_i)

        real(8) :: f_init_2d_double_hump, f_init
        integer, intent(IN) :: nd
        real(8), dimension(nd), intent(IN) :: x_i

        ! Parameters of function
        real(8), parameter :: cons1 = 2.0/sqrt(3.14159265359) 
        
        if (nd == 2) then
            f_init = cons1*x_i(1)*exp(-(x_i(2)-2.0)**2-(x_i(1)-2.0)**2) + cons1*x_i(1)*exp(-(x_i(2)+2.0)**2-(x_i(1)-2.0)**2)
        else
            f_init = 0.0
        end if

        f_init_2d_double_hump = f_init

    end function f_init_2d_double_hump

    ! Set of stochastic differential equations giving dV increment to v_par, v_perp
    function dV_vpar_vperp(nd, delta_t, x_i, dw)
            
        integer, intent(IN) :: nd
        real(8), dimension(nd) :: dV_vpar_vperp, dv
        real(8), intent(IN) :: delta_t
        real(8), dimension(nd), intent(IN) :: x_i, dw

        ! Declare Parameters
        real(8), parameter :: cons1 = 2.0/sqrt(3.14159265359) 
        real(8), parameter :: alpha = 1.0

        ! Declare variables for microscopic integration
        real(8) :: vpar, vperp
        real(8) :: v, H, nu_d, erfv, A
        real(8), dimension(nd) :: d_eta 

        ! Initialize variables
        d_eta = dw
        vpar = x_i(2)
        vperp = x_i(1)

        ! Check if particle in illegal area
        if ( vperp < v_perp_min ) then
            vperp = v_perp_min
        end if

        ! Calculate erfv, v, H, nu_d, B, A
        v = sqrt(vpar**2.0 + vperp**2.0)
        erfv = erf(v)
        H = (erfv - cons1*v*exp(-v**2.0)) / (v**3.0)
        nu_d = (erfv / v**3.0) - (0.5*H / v**2.0)
        A = -(0.5/v + alpha*v)*H + cons1*(1.0/v)*exp(-v**2.0)

!        ! Fit Given random numbers to proper distributuion 
        d_eta = sqrt(delta_t)*sign(1.0_8,1.0-2.0*d_eta)

        ! Calculate v_par and v_perp at next time step
        dv(2) = (-nu_d*vpar+A*vpar/v)*delta_t + (-vperp*sqrt(nu_d))*d_eta(1) + (vpar*sqrt(H)/v)*d_eta(2)
        dv(1) = (0.5*nu_d*(vpar**2-vperp**2)/vperp+A*vperp/v)*delta_t+(vpar*sqrt(nu_d))*d_eta(1) + (sqrt(H)*vperp/v)*d_eta(2)

        ! Do not let particle into illegal area 
        if ( (x_i(1) + dv(1)) < v_perp_min ) then
            dv(1) = -(x_i(1) - v_perp_min)
        end if

        dV_vpar_vperp = dv
    end function dV_vpar_vperp 

end module test_functions
