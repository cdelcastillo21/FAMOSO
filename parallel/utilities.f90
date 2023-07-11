! Carlos del-Castillo-Negrete, April 2019
! UTILITIES module
! Conatains utility subroutines and functions to 
! facilitate the manipulation and input/output of matrices
module utilities
        use ifport 
     
contains

        subroutine start_timer(t)
            integer, intent(OUT) :: t
            
            call system_clock(t)
        end subroutine start_timer

        function get_time(t)
            real(8) :: get_time
            integer, intent(IN) :: t
            integer :: curr_t, clck_rate

            call system_clock(curr_t,clck_rate)

            get_time = real(curr_t-t)/real(clck_rate)
        end function get_time

        ! Rounds x down to nearest multiple of y
        function ROUND_DOWN(x, y)

            real(8) :: ROUND_DOWN
            real(8), intent(IN) :: x,y
            
            if ( mod(x,y) /= 0.0 ) then
                ROUND_DOWN = x - (y - abs(mod(x,y)))
            else
                ROUND_DOWN = x
            end if
        end function ROUND_DOWN

        ! Rounds x up to nearest multiple of y
        function ROUND_UP(x, y)
            real(8) :: ROUND_UP
            real(8), intent(IN) :: x,y 

            if ( mod(x,y) /= 0.0 ) then
                ROUND_UP = x + (y - abs(mod(x,y)))
            else
                ROUND_UP = x
            end if

        end function ROUND_UP

        ! Auxiliary routine to print a matrix A to stdout                                                     
        subroutine PRINT_MATRIX( desc, m, n, A)
                
                character(*) :: desc
                integer :: m, n
                integer :: i, j
                real(8), dimension(:,:) :: A
                
                write(*,*)
                write(*,*) desc

                do i=1,m
                        write(*,9998) ( A(i,j), j=1,n )
                end do

                9998 FORMAT(100f10.5)
                return
        end subroutine print_matrix

        ! Auxiliary routine to print a set of N matrices to stdout                                                              
        subroutine PRINT_MATRICES_N( desc, r, c, N, A)
                
                character(*) :: desc
                character(10) :: sub_desc
                integer :: r, c, N
                integer :: k
                real(8), dimension(:,:,:) :: A
               
                write(*,*) "-------------------------" 
                write(*,'(a4, a2, i2, a3, i2, a3, i2, a2)') desc, "[ ", r, " X ", c, " X ", N,  " ]"

                do k=1,N
                        write(sub_desc,9997) "(:,:,", k, ")"
                        call print_matrix(sub_desc, r, c, A(:,:,k))
                end do

                write(*,*) "-------------------------" 

                9997 FORMAT(a5, i2, a1)
                return
        end subroutine print_matrices_n

        ! Auxiliary routine to print a vector v to stdout
        subroutine PRINT_VECTOR( desc, l, v)

                character(*) :: desc
                integer :: l, i
                real(8), dimension(:) :: v

                write(*,*)
                write(*,*) desc

                do i=1,l
                        write(*,'(f10.5)') v(i)
                end do

                return
        end subroutine 

        ! Outputs a given matrix A to a .data file that 
        ! should be readable by matlab/octave
        subroutine OUTPUT_MATRIX( desc, m, n, A)
                
                character(*) :: desc
                integer :: i, j, m, n
                real(8), dimension(:,:) :: A

                ! Open file with appropriate name -> overwrite existing
                open(unit=1,file=desc)

                ! Print formatted matrix
                do i=1,m
                        write(1,9998) ( A(i,j), j=1,n )
                end do

                ! Close file
                close(1)
                
                9998 FORMAT(100f10.5)
                return
        end subroutine OUTPUT_MATRIX

        ! Initializes a random number generator. Parallel safe
        
        subroutine INIT_RANDOM_SEED()
            use iso_fortran_env, only: int64
            implicit none
            integer, allocatable :: seed(:)
            integer :: i, n, un, istat, dt(8), pid
            integer(int64) :: t
          
            call random_seed(size = n)
            allocate(seed(n))
            ! First try if the OS provides a random number generator
            open(newunit=un, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old", iostat=istat)
            if (istat == 0) then
               read(un) seed
               close(un)
            else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(t)
               if (t == 0) then
                  call date_and_time(values=dt)
                  t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24_int64 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
               end if
               pid = getpid()
               t = ieor(t, int(pid, kind(t)))
               do i = 1, n
                  seed(i) = lcg(t)
               end do
            end if
            call random_seed(put=seed)
          contains
            ! This simple PRNG might not be good enough for real work, but is
            ! sufficient for seeding a better PRNG.
            function lcg(s)
              integer :: lcg
              integer(int64) :: s
              if (s == 0) then
                 s = 104729
              else
                 s = mod(s, 4294967296_int64)
              end if
              s = mod(s * 279470273_int64, 4294967291_int64)
              lcg = int(mod(s, int(huge(0), int64)), kind(0))
            end function lcg
        end subroutine INIT_RANDOM_SEED       

end module utilities

