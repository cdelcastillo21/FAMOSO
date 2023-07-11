! Carlos del-Castillo-Negrete, March 23, 2015
! FAMOSO - FAst MOnte Carlo SOLver - GLRAM module
! Module implementing Generalized Low Rank Approximation of Matrices
! algorithm (GLRAM). Defines a glram representation data structure
! that consists of the compressed matrices M_i, the transformation 
! matrices LT and RT, and the dimensions/size of the data set. 
module glram
        
        use utilities

        ! A full data set of N rXc matrices A_i can be 
        ! compressed into this data structure using GLRAM decomposition
        type glram_rep
                integer :: r, c, N, l1, l2
                real(8), dimension(:,:), allocatable :: LT, RT
                real(8), dimension(:,:,:), allocatable :: M_i
        end type glram_rep

        ! Constants of glram algorithm
        integer, parameter :: max_iter = 20                  
        real(8), parameter :: eps = 0.000001                
contains
        
        ! Destroys the given g_rep by freeing associated memory
        ! g_rep must have been already allocated (i.e. returned from glram_construt())
        subroutine glram_destroy(g_rep)
                type(glram_rep), intent(INOUT) :: g_rep

                deallocate(g_rep%LT)
                deallocate(g_rep%RT)
                deallocate(g_rep%M_i)
        end subroutine glram_destroy

        ! Prints to stdout the components of a g_rep structure.
        subroutine glram_print(g_rep)
            type(glram_rep), intent(IN) :: g_rep
            integer :: i

            write(*,*) "-----------------"
            write(*,*) "N = ", g_rep%N, ", r = ", g_rep%r, ", c = ", g_rep%c
            write(*,*) "l1 = ", g_rep%l1, ", l2 = ", g_rep%l2 
            call print_matrix("LT: ", g_rep%r, g_rep%l1, g_rep%LT)
            call print_matrix("RT: ", g_rep%c, g_rep%l2, g_rep%RT)
            do i=1,g_rep%N
                call print_matrix("M_i", g_rep%l1, g_rep%l2, g_rep%M_i(:,:,1))
            end do 

        end subroutine glram_print

        ! Returns glram_rep structure containing the compressed representation
        ! of the N, rXc matrices stored in A_i. Compressed matrices are of dimension l1xl2
        function glram_construct( A_i, r, c, N, l1, l2)
                
                implicit none

                type(glram_rep) :: glram_construct
                type(glram_rep) :: g_rep
                integer, intent(IN) :: r, c, N, l1, l2
                real(8), dimension(r,c,N), intent(IN) :: A_i

                ! Declaration of subroutine variables
                integer :: i, j, k, num_eigen, lwork, lwork_max, info, il, iu
                integer, dimension(:), allocatable :: iwork, ifail
                real(8) :: diff_lt, diff_rt 
                real(8), dimension(:), allocatable :: eigen_vals, work, temp_v
                real(8), dimension(:,:), allocatable :: M_r, M_l, TR, TL, LT_prev, RT_prev, temp_init

                ! Declare parameters 
                real(8), parameter :: abstol = -1.0             ! -1.0 corresponds to default value
                real(8), parameter :: alpha_1 = 1.0, beta_0 = 0.0, beta_1 = 1.0

                ! Declare external functions to compute svd and matrix products
                external DGEMM
                external DGESVX
               
                ! intrinsic functions
                intrinsic MIN, MAX, INT

                ! Initialieze g_rep structure 
                g_rep%r = r
                g_rep%c = c
                g_rep%N = N
                g_rep%l1 = l1
                g_rep%l2 = l2
                
                ! Allocate space for matrices and compressed rep
                allocate(g_rep%LT(r,l1))
                allocate(g_rep%RT(c,l2))
                allocate(g_rep%M_i(l1,l2,N))

                ! Set LAPACK lwork_max variable 
                lwork_max = 8*MAX(r,c)

                ! Allocate data structures
                allocate(temp_v(MAX(r,c)))
                allocate(temp_init(r,c))
                allocate(M_r(c,c))
                allocate(M_l(r,r))
                allocate(TR(c,l1))
                allocate(TL(r,l2))
                allocate(LT_prev(r,l1))
                allocate(RT_prev(c,l2))
                allocate(eigen_vals(MAX(l1,l2)))
                allocate(iwork(5*MAX(r,c)))
                allocate(ifail(MAX(r,c)))
                allocate(work(lwork_max))


                ! Set upper and lower boundary of eigenvalues we want -> largest l1 eigenvectors of A(1) 
                iu = r 
                il = r-l1+1 
                num_eigen =iu-il+1

                ! First query for optimal workspace -> lwork = -1 for query
                lwork = -1
                temp_init(:,:) = A_i(:,:,1)
                CALL DSYEVX('V','I','U', r, temp_init, r, 0.0, 0.0, il, iu, abstol, num_eigen, eigen_vals(1:l1), &
                           &g_rep%LT, r, work, lwork, iwork, ifail, info) 
                lwork = MIN(lwork_max, INT(work(1)))
 
                ! Perform actual computation 
                CALL DSYEVX('V','I','U', r, temp_init, r, 0.0, 0.0, il, iu, abstol, num_eigen, eigen_vals(1:l1), &
                           &g_rep%LT, r, work, lwork, iwork, ifail, info) 

                ! Iterative algorithm to find LT,RT
                do i=1,max_iter
                        ! Reset temp arrays
                        TR(:,:) = 0
                        TL(:,:) = 0
                        M_r(:,:) = 0
                        M_l(:,:) = 0

                        ! Copy previous values of transformations
                        LT_prev(:,:) = g_rep%LT(:,:);
                        RT_prev(:,:) = g_rep%RT(:,:);

                        ! Form matrix M_r
                        do j=1,N
                                ! Compute TL = A_i*LT
                                CALL DGEMM('T','N', c, l1, r, alpha_1, A_i(:,:,j), r, g_rep%LT, r, beta_0, TR, c)
                                
                                ! Add on to M_R -> M_R = M_R + (TR)TR*
                                CALL DGEMM('N','T', c, c, l1, alpha_1, TR, c, TR, c, beta_1, M_r, c)         
                        end do
                        
                        ! Compute RT by finding eigenvectors of M_r
                        ! Set upper and lower boundary of eigenvalues we want -> largest l2 eigenvectors
                        iu = c
                        il = c-l2+1
                        num_eigen = iu-il+1

                        ! First query for optimal workspace -> lwork = -1 for query
                        lwork = -1
                        CALL DSYEVX('V','I','U', c, M_r, c, 0.0, 0.0, il, iu, abstol, num_eigen, eigen_vals(1:l2), &
                                   &g_rep%RT, c, work, lwork, iwork, ifail, info) 
                        lwork = MIN(lwork_max, INT(work(1)))

                        ! Perform actual computation
                        CALL DSYEVX('V','I','U', c, M_r, c, 0.0, 0.0, il, iu, abstol, num_eigen, eigen_vals(1:l2), &
                                   &g_rep%RT, c, work, lwork, iwork, ifail, info) 
   

                        ! Form matrix M_l using new R_i
                        do j=1,N
                                ! Compute TL = AR CHECK THIS
                                CALL DGEMM('N','N', r, l2, c, alpha_1, A_i(:,:,j), r, g_rep%RT, c, beta_0, TL, r)
                                
                                ! Add on to M_L -> M_L = M_L + (TL)TL*
                                CALL DGEMM('N','T', r, r, l2, alpha_1, TL, r, TL, r, beta_1, M_l, r)         
                        end do
                        
                        ! Compute LT by finding eigenvectors of M_l
                        ! Set upper and lower boundary of eigenvalues we want -> largest l1 eigenvectors of M_l
                        iu = r 
                        il = r-l1+1 
                        num_eigen =iu-il+1

                        ! First query for optimal workspace -> lwork = -1 for query
                        lwork = -1
                        CALL DSYEVX('V','I','U', r, M_l, r, 0.0, 0.0, il, iu, abstol, num_eigen, eigen_vals(1:l1), &
                                   &g_rep%LT, r, work, lwork, iwork, ifail, info) 
                        lwork = MIN(lwork_max, INT(work(1)))

                        ! Perform actual computation
                        CALL DSYEVX('V','I','U', r, M_l, r, 0.0, 0.0, il, iu, abstol, num_eigen, eigen_vals(1:l1), &
                                   &g_rep%LT, r, work, lwork, iwork, ifail, info) 

                        ! Check for convergence by finding net change between RT, LT and see if its < eps
                        ! Only check after the first iteration of the algorithm. 
                        if ( i > 1 ) then
                            diff_lt = 0.0
                            do j=1,l1
                                do k=1,r
                                    diff_lt = diff_lt + (g_rep%LT(k,j)-LT_prev(k,j))**2
                                end do
                            end do
                            diff_rt = 0.0
                            do j=1,l2
                                do k=1,c
                                diff_rt = diff_rt + (g_rep%RT(k,j)-RT_prev(k,j))**2
                                end do
                            end do

                            if ( diff_lt < eps .and. diff_rt > eps ) then
                                write(*,*) "glram:: glram_construct - LT converged (", diff_lt, "), rt has not on iteration ", i
                            else if ( diff_rt < eps .and. diff_lt > eps ) then
                                write(*,*) "glram:: glram_construct - RT converged (", diff_rt, "), rt has not on iteration ", i
                            else if ( diff_lt < eps .and. diff_rt < eps ) then 
                                write(*,*) "glram:: glram_construct - RT (", diff_rt, ") and LT (", diff_lt, &
                                    &") converged on iteraetion ", i
                                EXIT
                            end if
                        end if 

                end do

                ! Rearrange LT and RT -> Column vectors must be in opposite
                do i=1,l2/2
                        ! Copy vector and swap
                        temp_v(1:c) = g_rep%RT(:,i)
                        g_rep%RT(:,i) = g_rep%RT(:,l2-i+1)
                        g_rep%RT(:,l2-i+1) = temp_v(1:c)        
                end do
                do i=1,l1/2
                        temp_v(1:r) = g_rep%LT(:,i)
                        g_rep%LT(:,i) = g_rep%LT(:,l1-i+1)
                        g_rep%LT(:,l1-i+1) = temp_v(1:r)
                end do

                ! Compute M_i -> M_i = LT*A_iR
                do i=1,N
                        ! Compute A_iR -> store temporarily in TL
                        CALL DGEMM('N','N', r, l2, c, alpha_1, A_i(:,:,i), r, g_rep%RT, c, beta_0, TL, r) 

                        ! Now M_i = LT*TL
                        CALL DGEMM('T', 'N', l1, l2, r, alpha_1, g_rep%LT, r, TL, r, beta_0, g_rep%M_i(:,:,i), l1) 
                end do

                deallocate(M_r)
                deallocate(M_l)
                deallocate(TR)
                deallocate(TL)
                deallocate(eigen_vals)
                deallocate(iwork)
                deallocate(ifail)
                deallocate(work)
                deallocate(temp_v)
                
                glram_construct = g_rep
        end function glram_construct

        ! Reconstructs the given coefficient matrix M into an original rXc sized matrix using
        ! the transformations stored in the g_rep structure. 
        function glram_reconstruct(M, g_rep)
            real(8), dimension(:,:), allocatable :: glram_reconstruct
            real(8), dimension(:,:), intent(IN) :: M
            type(glram_rep) :: g_rep

            real(8), dimension(:,:), allocatable :: temp, A_recon
            real(8), parameter :: alpha_1 = 1.0, beta_0 = 0.0

            allocate(temp(g_rep%l1,g_rep%c))
            allocate(A_recon(g_rep%r,g_rep%c))

            ! Reconstruct projection to original data matrix sizeSet constants
            ! A_proj = LT(M_proj)RT^T

            ! Compute temp = (M_proj)RT^T
            CALL DGEMM('N', 'T', g_rep%l1, g_rep%c, g_rep%l2, alpha_1, M, g_rep%l1, g_rep%RT, g_rep%c, beta_0, temp, g_rep%l1)

            ! Compute A_proj = LT(temp)
            CALL DGEMM('N', 'N', g_rep%r, g_rep%c, g_rep%l1, alpha_1, g_rep%LT, g_rep%r, temp, g_rep%l1, beta_0, A_recon, g_rep%r)

            deallocate(temp)

            glram_reconstruct = A_recon
        end function glram_reconstruct

        ! Returns reconstructed data in A_i according to stored glram representaiton in g_rep
        ! Note - subroutine assumes memory for A_i already exists and is of proper size
        subroutine glram_reconstruct_all(A_i, g_rep)
                real(8), intent(INOUT), dimension(:,:,:) :: A_i
                type(glram_rep), intent(IN) :: g_rep

                ! Declaration of subroutine variables
                integer :: i, N, r, c, l1, l2
                real(8), dimension(:,:), allocatable :: temp
 
                ! Declare parameters 
                real(8), parameter :: alpha_1 = 1.0, beta_0 = 0.0

                ! Set constants
                N = g_rep%N
                r = g_rep%r
                c = g_rep%c
                l1 = g_rep%l1
                l2 = g_rep%l2

                ! Allocate space for temp variables
                allocate(temp(l1,c))

                ! Reconstruct original data points according to A_i = LT(M_i)(RT*)                
                do i=1,N
                
                        ! Compute temp = (M_i)RT*
                        CALL DGEMM('N', 'T', l1, c, l2, alpha_1, g_rep%M_i(:,:,i), l1, g_rep%RT, c, beta_0, temp, l1)

                        ! Compute A_i = LT(temp)
                        CALL DGEMM('N', 'N', r, c, l1, alpha_1, g_rep%LT, r, temp, l1, beta_0, A_i(:,:,i), r)
                end do
               
                ! Free used memory 
                deallocate(temp)

                RETURN
        end subroutine glram_reconstruct_all

end module glram
