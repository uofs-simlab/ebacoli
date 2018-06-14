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
! This file contains the BACOLI source code (including the new code
! that implements the interpolation based spatial error estimates.)

program curve_example_driver

    use ebacoli95_mod, only: ebacoli95_init, ebacoli95, ebacoli95_vals
    use ebacoli95_mod, only: ebacoli95_sol, ebacoli95_sol_teardown
    use ebacoli95_mod, only: ebacoli95_output_Bsplines

    implicit none
    integer, parameter     :: dp = kind(0d0)
    type(ebacoli95_sol)     :: sol

    integer,  parameter    :: npde = 4
    integer,  parameter    :: nu   = 4
    integer,  parameter    :: nderiv = 2
    real(dp), parameter    :: xa = 0, xb = 1
    integer,  parameter    :: nin = 10

    integer                :: nout, ntout
    real(dp), allocatable  :: xout(:), uout(:,:,:)
    real(dp)               :: tout, tstart, tstop, atol(npde), rtol(npde), xin(nin)

    integer                :: i, j, ier

    external f, bndxa, bndxb, uinit, derivf, difbxa, difbxb


    ! The following parameters are used in the rcdsys.f example.
    real(dp) :: a1,a2,d1,d2,r,c,n,pe1,pe2
    common /rcdstuff/ a1,a2,d1,d2,r,c,n,pe1,pe2
    a1 = 30._dp
    a2 = 30._dp
    d1 = 1.5_dp
    d2 = 1.2_dp
    r = 1000._dp
    c = .96_dp
    n = 1._dp
    pe1 = 1e4_dp
    pe2 = 1e4_dp


    !-------------------------------------------------------------------
    ! Write out value of npde to allow user to confirm that its value
    ! is appropriate for the problem to be solved.

    write(6,*) 'The total number of equations is ', npde
    write(6,*) 'The size of the parabolic subsystem is ', nu
    write(6,*)

    ! Set output time
    tstop = 72.d0
    ntout = 288

    ! Set tolerance
    atol = 1d-4
    rtol = atol
    !    atol(4) = 1d-1

    ! Initialize a grid to pass to ebacoli
    xin = (/(xa+i*(xb-xa)/(nin-1), i=0,nin-1)/)

    !-------------------------------------------------------------------
    ! Initialization: Allocate storage and set problem parameters.
    call ebacoli95_init(sol, npde, nu, xin, atol=atol, rtol=rtol)
    tstart = sol%t0
    tout = tstart + 1.d-12


    !-------------------------------------------------------------------
    ! Integrate solution from t=0 to t=tout.
    print '(/"THE INPUT IS")'
    print 900, sol%kcol, sol%nint, sol%npde, tstop
    print 901, sol%atol(1), sol%rtol(1), "LOI"

    ! Compute the solution at tout
    call ebacoli95(sol, tout, f, bndxa, bndxb, uinit, derivf, difbxa, difbxb)

    call ebacoli95_output_Bsplines(sol, 0, tstart)

    do j = 1, ntout

       tout = tstart + (j)*(tstop-tstart)/(ntout)

       ! Compute the solution at tout
       call ebacoli95(sol, tout, f, bndxa, bndxb, uinit, derivf, difbxa, difbxb)
       if (sol%idid <= 0) goto 800
    !-------------------------------------------------------------------
    ! Output results.
       if (sol%idid > 0) then
          call ebacoli95_output_Bsplines(sol, j, tout)
       end if

    end do

    !-------------------------------------------------------------------
    ! The end.
    call ebacoli95_sol_teardown(sol) ; stop
600 print '("Error: Improperly formatted input")' ; stop
700 print '("Error: Could not allocate storage")' ; stop
800 print '("Error: Was not able to integrate to tstop")' ; stop

    !-------------------------------------------------------------------
    ! Formats
900 format("kcol = ",i2,", nint0 = ",i4,", npde = ",i3,", tout = ",es7.1)
901 format("atol = ",es7.1,", rtol = ",es7.1,",",17x,a3)
902 format("Number of subintervals in the current mesh:",i8)
end program
