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

program colbeck_driver

    use ebacoli95_mod, only: ebacoli95_init, ebacoli95, ebacoli95_vals
    use ebacoli95_mod, only: ebacoli95_sol, ebacoli95_sol_teardown
    use ebacoli95_mod, only: ebacoli95_output_Bsplines

    implicit none
    integer, parameter     :: dp = kind(0d0)
    type(ebacoli95_sol)     :: sol

    integer,  parameter, dimension(3) :: npde_sub = (/2,0,0/)
    integer,  parameter    :: npde = 2
    integer,  parameter    :: nu   = 2
    real(dp), parameter    :: xa = 0d0, xb = 1.d0

    ! Specify the two sections of the input mesh:
    ! - 1 from xa to xsplit, nint1 number of intervals
    ! - 2 from xsplit to xb, nint2 number of intervals
    real(dp), parameter    :: xsplit = 0.5d0
    integer,  parameter    :: nint1 = 10, nint2 = 10
    integer,  parameter    :: nin = nint1+nint2+1
    real(dp), parameter    :: h1 = (xsplit-xa)/nint1, h2 = (xb-xsplit)/nint2
    real(dp)               :: xin(nin)

    integer                :: nout, ntout, kcol
    real(dp), allocatable  :: xout(:), uout(:,:,:)
    real(dp)               :: tout, tstart, tstop, atol(npde), rtol(npde)

    integer                :: i, j, ier

    external f, bndxa, bndxb, uinit, derivf, difbxa, difbxb

    ! parameters
    double precision rhoair,rholiq,rhoice,cair,cliq,cice,Lfus,omega,Tfrz,thetar,ksnow,eps
    integer ns
    common /colbeck/ rhoair,rholiq,rhoice,cair,cliq,cice,Lfus,omega,Tfrz,thetar,ksnow,eps,ns
    ! parameters not needed in function evaluation
    doubleprecision porosity, permeability

    double precision inton,smstep

    rhoair = 1.2754d0
    rholiq = 1000.d0
    rhoice = 917.d0

    cair = 1005.d0
    cliq = 4181.d0
    cice = 2114.d0

    Lfus = rhoice*cice/(rholiq*0.0058d0)

    Tfrz = 273.15d0
    omega = 50.d0

    porosity = 0.673
    thetar = 0.07d0*porosity  ! Swi*phi

    permeability = 2.967d-8
    ksnow = rholiq*9.81*permeability/1.781e-3

    ns = 3

    eps = 0.1d0

    ! write(*,*) smstep(0.d0,1.d0)
    ! write(*,*) smstep(3.d0*3600.d0-0.d0,1.d0)
    ! write(*,*) 0.1d0*inton(0.d0,3000.d0,3.d0*3600.d0,1.d0)
    !-------------------------------------------------------------------
    ! Write out value of npde to allow user to confirm that its value
    ! is appropriate for the problem to be solved.

    write(6,*) 'The total number of equations is ', npde
    write(6,*) 'The size of the parabolic subsystem is ', nu
    write(6,*)

    ! Set output time
    tstop = 1000.d0
    ntout = 100

    ! Set tolerance
    atol = 1d-5
    rtol = atol

    ! Set kcol (order of expansion is p = kcol+1)
    ! kcol=4 is default
    kcol = 4

    ! Initialize a grid to pass to ebacoli
    xin(1) = xa
    do i = 1, nint1
       xin(i+1) = xa+i*h1
    end do
    do i = 1, nint2
       xin(nint1+i+1) = xsplit + i*h2
    end do

    !-------------------------------------------------------------------
    ! Initialization: Allocate storage and set problem parameters.
    call ebacoli95_init(sol, npde_sub, xin, atol=atol, kcol=kcol, rtol=rtol, nint_max=1000)
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

    ! Write out some of the solver data
    print '("final time      = ",e7.2)', sol%time
    print '("final step size = ",e7.2)', sol%prev_time_step_size
    print '("steps taken     = ",i10)', sol%num_accepted_time_steps
    print '("num remeshings  = ",i10)', sol%num_remeshings
    print '("ini remeshings  = ",i10)', sol%num_ini_remeshings
    print '("cold restarts   = ",i10)', sol%num_cold_restarts

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
end program ! colbeck_driver
