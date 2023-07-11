! Carlos del-Castillo-Negrete, March 23, 2015
! FAMOSO - FAst MOnte Carlo SOLver - MC_GRID module
! Module implementing functions corresponding to the mc_grid type.
module mc_grid
 
        use utilities

        type grid
                integer :: nd                                        ! # Dimensions in grid
                integer, dimension(:), allocatable :: szs            ! # grid points for each dimension
                real(8), dimension(:), allocatable :: delta_x_i      ! width of each grid cell
                real(8), dimension(:), allocatable :: upper_bnds     ! upper bounds of grid 
                real(8), dimension(:), allocatable :: lower_bnds     ! lowerb bounds of grid
        end type grid

contains

        ! Returns an allocated grid structure with the given paraemeters. 
        function mc_grid_build(nd, delta_x_i, lower_bnds, upper_bnds)
                type(grid) :: mc_grid_build, gd
                integer, intent(IN) :: nd
                real(8), dimension(nd), intent(IN) :: delta_x_i, lower_bnds, upper_bnds

                integer :: i 

                ! Allocate space for grid arrays
                allocate(gd%szs(nd))
                allocate(gd%upper_bnds(nd))
                allocate(gd%lower_bnds(nd))
                allocate(gd%delta_x_i(nd))

                ! Copy over grid variable
                gd%nd = nd
                
                do i=1,nd
                    gd%delta_x_i(i) = delta_x_i(i) 
                    gd%lower_bnds(i) = ROUND_DOWN(lower_bnds(i),delta_x_i(i))
                    gd%upper_bnds(i) = ROUND_UP(upper_bnds(i),delta_x_i(i))
                    gd%szs(i) =  int((gd%upper_bnds(i)-gd%lower_bnds(i))/gd%delta_x_i(i))
                end do

                ! Return grid structure
                mc_grid_build = gd
        end function mc_grid_build

        ! Prints the message contained in desc and the grid information contained in gd to stdout.
        subroutine mc_grid_print( desc, gd )
            character(*) :: desc
            type(grid), intent(IN) :: gd

            integer :: i 

            write(*,*) desc 
            do i=1,gd%nd
                write(*,*) "    nd = ", gd%nd
                write(*,*) "    delta_x(", i, ") = ", gd%delta_x_i(i)
                write(*,*) "    lower_bound(", i, ") = ", gd%lower_bnds(i)
                write(*,*) "    upper_bound(", i, ") = ", gd%upper_bnds(i)
                write(*,*) "    szs(", i, ") = ", gd%szs(i)
            end do

        end subroutine mc_grid_print

        ! Returns an allocated grid structure that is an exact copy of the given grid in gd.
        function mc_grid_copy( gd )
                type(grid) :: new_gd, mc_grid_copy
                type(grid), intent(IN) :: gd

                new_gd = mc_grid_build(gd%nd, gd%delta_x_i, gd%lower_bnds, gd%upper_bnds)
           
                mc_grid_copy = new_gd 
        end function mc_grid_copy


        ! Rescales the current grid to fit the new_upper and new_lower bounds in each dimension
        ! but without modifying the original grid scpacing delta_x_i.
        subroutine mc_grid_rescale(curr_gd, new_upper, new_lower)
                type(grid), intent(INOUT) :: curr_gd
                real(8), dimension(:), intent(INOUT) :: new_upper, new_lower 

                integer :: i,nd
                integer, dimension(:), allocatable :: new_szs

                ! Check arguments
 
                ! Initialize local variables
                nd = curr_gd%nd

                ! Allocate space for local variables
                allocate(new_szs(nd))

                ! Scale up grid sizes to preserve delta_x_i spacing in each dimension
                ! Round upper/lower bounds to closest multiple of delta_x_i
                do i=1,nd
                    ! Round upper bnd up 
                    new_upper(i) = ROUND_UP(new_upper(i),curr_gd%delta_x_i(i))
                                             
                    ! Round lower bnd down
                    new_lower(i) = ROUND_DOWN(new_lower(i),curr_gd%delta_x_i(i))

                    ! Compute new grid spacing 
                    new_szs(i) = int(abs(new_upper(i) - new_lower(i))/curr_gd%delta_x_i(i)) 
                end do

                ! update curr_gd
                curr_gd%upper_bnds(:) = new_upper(:)
                curr_gd%lower_bnds(:) = new_lower(:)
                curr_gd%szs(:) = new_szs(:)

                ! Free used memory
                deallocate(new_szs)
        end subroutine mc_grid_rescale

        ! Frees allocated memory associated with a gd structure. 
        subroutine mc_grid_destroy(gd)
                type(grid) :: gd

                deallocate(gd%szs)
                deallocate(gd%upper_bnds)
                deallocate(gd%lower_bnds)
        end subroutine mc_grid_destroy

end module mc_grid
