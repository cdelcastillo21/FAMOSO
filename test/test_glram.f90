! Carlos del-Castillo-Negrete 
! April, 2015
! test_glram - Program to test the glram algorithm implemented in the FAMOSO library.
! Function computes tensor product decomposition of a random set of matrices and prints
! out the results of the decomposition. 
program test_glram
        
        use utilities
        use glram
        
        implicit none

        ! Declare local variables
        integer, parameter :: r = 4, c = 4, N = 3
        integer, parameter :: l1 = 2, l2 = 2
        real(8), dimension(r,c,N) :: A_i, A_i_reconstructed
        type(glram_rep) :: g_rep
        character(100) :: full_path
               
        write(*,*) "--------------------TEST GLRAM PROCEDURE---------------"
        write(*,*) "r = ", r, "c = ",  c, "N = ",  N, "l1 = ", l1, "l2 = ", l2

        ! Initiate random number generator
        call init_random_seed()

        ! Initialize A_i's to random numbers
        call random_number(A_i)

        ! Call glram function to construct compressed representation of A_i
        g_rep = glram_construct( A_i, r, c, N, l1, l2)

        ! Reconstruct A_i
        CALL glram_reconstruct_all(A_i_reconstructed, g_rep)

        ! Output results to stdout
        CALL PRINT_MATRIX("RT", c, l2, g_rep%RT)
        CALL PRINT_MATRIX("LT", r, l1, g_rep%LT)
        CALL PRINT_MATRICES_N("M_i", l1, l2, N, g_rep%M_i)
        CALL PRINT_MATRICES_N("A_i", r, c, N, A_i)
        CALL PRINT_MATRICES_N("A_R", r, c, N, A_i)

        ! Destroy allocated memory
        CALL glram_destroy(g_rep)
end program

