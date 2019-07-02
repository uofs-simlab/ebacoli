! Copyright (c) 2013, Paul Muir, Jack Pew
! Paul Muir, Mathematics and Computing Science, Saint Mary's University.
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
! * Redistributions of source code must retain the above copyright
!   notice, this list of conditions and the following disclaimer.
! * Redistributions in binary form must reproduce the above copyright
!   notice, this list of conditions and the following disclaimer in the
!   documentation and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

program test_toy2_par_plus_ode

    ! A minimal example driver for the EBACOLI95 wrapper of EBACOLI.
    ! See the comments within ebacoli95.f95 for details regarding use.

    use ebacoli95_mod, only: ebacoli95_init, ebacoli95, ebacoli95_vals, &
                             ebacoli95_sol, ebacoli95_sol_teardown

    implicit none
    integer, parameter  :: dp = kind(0d0)

    ! Declare data structure that will hold solution
    type(ebacoli95_sol) :: sol

    ! Subsystem sizes
    integer,  parameter :: nu = 2, nv = 2, nw = 0
    integer,  parameter, dimension(3) :: npde_sub = (/nu,nv,nw/)
    integer,  parameter :: npde = sum(npde_sub)

    ! Set boundary locations
    double precision, parameter :: xL = 1._dp, xR = 2._dp

    ! Maximum number of subintervals to allow
    integer,  parameter :: nint_max = 500

    ! Set (problem dependent) output time to 1
    real(dp), parameter :: tout = 1.D0

    ! Assume output at 11 points over (problem specific) spatial domain,
    integer,  parameter :: nout = 11

    integer :: i, k, ij

    external f, bndxa, bndxb, uinit, derivf, difbxa, difbxb, truu

    ! BEGIN executable statements

    ! Initialization (set spatial domain = [xL,xR]); a default uniform
    ! spatial mesh having 10 subintervals will be constructed.
    call ebacoli95_init(sol, npde_sub, (/xL,xR/), nint_max = nint_max)

    call write_subsystem_sizes(sol%npde_sub)

    ! Compute solution at tout
    call ebacoli95(sol, tout, f, bndxa, bndxb, uinit, derivf, difbxa, difbxb)

    ! Output solution over uniformly spaced points if solution successful
    print '("idid=",i5)',sol%idid
    if (sol%idid > 0) then
       print '("idid > 0 => Successful computation")'
       call write_solution_uniformly_spaced(sol,nout)
    end if

    ! Display actual error
    call compare_L2_error(sol,truu)

    call ebacoli95_sol_teardown(sol)

end program
