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
! This file contains the EBACOLI source code (including the new code
! that implements the interpolation based spatial error estimates.)

!-----------------------------------------------------------------------
! This is an example driver for running EBACOLI through the F95 wrapper.
! See the comments within ebacoli95.f95 for details regarding use.

! This driver writes output for each PDE using EBACOLI's spatial mesh
! points at ntout uniform points along the temporal domain. tstart is 0,
! tstop is the end time, and the spatial domain is [xa,xb]. A number of
! parameters are hard coded and should be edited as needed: npde (the
! number of PDEs), xa, xb. Several other parameters take default values
! through the F95 wrapper.
!
! The output format used here allows the user to take advantage of
! EBACOLI's spatial adaptivity in their data visualisation, to some
! extent, but the trade-off is that data points are slightly more
! difficult to work with than those over a uniform grid.
!
! If returning after uniform time steps is not desirable,
! it is possible to modify this code to have ebacoli return after
! a specific number of time steps (as performed by DASSL).
! See the documentation inside ebacoli95.f95 for details.
!
! After a successful computation, this program will have written files
! named Points001, ..., PointsNPDE containing data points as ordered
! triples, one per line. The coordinates are ordered (X, T, U).
!
! This driver should be linked with with ebacoli.f, ebacoli-aux.f,
! ebacoli95.f95 and a problem definition file, eg, burg1.f.
!
!--------------------------------------------------------------------
! Example Python plotting code:
!
!   import matplotlib as mpl
!   mpl.use('AGG')  # for systems not running a GUI
!   from mpl_toolkits.mplot3d import Axes3D
!   from matplotlib import cm
!   import matplotlib.pyplot as plt
!   import numpy as np
!
!   styling = {
!     'cmap': cm.coolwarm,
!     'linewidth': 0,
!     'antialiased': True
!   }
!
!   x, t, u = np.loadtxt('Points1', unpack=True)
!   fig = plt.figure()
!   ax = fig.add_subplot(111, projection='3d')
!   ax.plot_trisurf(x, t, u, **styling)
!
!   ax.set_xlabel('$x$')
!   ax.set_ylabel('$t$')
!   ax.set_zlabel('$u(t,x)$')
!
!   plt.savefig('trimesh.png')
!
! axes3d.plot_trisurf() was added to matplotlib in version 1.2.0.
! It drapes a surface over a triangularization of the data.
!-----------------------------------------------------------------------

program trimesh_example_driver

    use ebacoli95_mod, only: ebacoli95_init, ebacoli95, ebacoli95_vals
    use ebacoli95_mod, only: ebacoli95_sol, ebacoli95_sol_teardown

    implicit none
    integer, parameter     :: dp = kind(0d0)
    type(ebacoli95_sol)     :: sol

    ! Choose npde to be consistent with the problem specific source
    ! code; see burg1.f, burg2.f, CahnAllen.f, rcdsys.f, sinmads.f,
    ! and steady.f.

    integer,  parameter, dimension(3)    :: npde_sub = (/2,2,2/)
    integer,  parameter    :: npde = sum(npde_sub)
    integer,  parameter    :: nu   = npde_sub(1), nv = npde_sub(2), nw = npde_sub(3)
    real(dp), parameter    :: xa = 0, xb = 70
    real(dp), allocatable  :: xinit(:)
    real(dp), allocatable  :: uout(:,:)
    real(dp)               :: tout, tstart, tstop, atol(npde), rtol(npde)

    integer                :: i, j, k, ier, ntout, nint_init
    character(len=32)      :: fname, npde_str

    external f, bndxa, bndxb, uinit, derivf, difbxa, difbxb

    ! System parameters
    double precision    epsilon, beta, gamma,  sigmai, sigmae
    common /fhn/       epsilon, beta, gamma, sigmai, sigmae

    epsilon = 0.1D0
    beta    = 1.D0
    gamma   = 0.5D0
    sigmai = 1.D0
    sigmae = 1.D0

    !-------------------------------------------------------------------
    ! Write out value of npde to allow user to confirm that its value
    ! is appropriate for the problem to be solved.
    call write_subsystem_sizes(npde_sub)

    ! Set solution parameters
    tstop = 30.d0
    ntout = 100
    atol(1) = 1d-2
    rtol(1) = atol(1)
    atol(4) = 1d-2
    rtol(4) = atol(4)

    ! set initial mesh
    nint_init = 40
    allocate(xinit(nint_init), stat=ier)
    do i = 1, 30
       xinit(i) = (i-1)*0.5
    end do
    do i = 31, 40
       xinit(i) = 15 + (i-30)*5.5
    end do
    ! write(*,*) xinit
    ! call exit(0)

    !-------------------------------------------------------------------
    ! Initialization: Allocate storage and set problem parameters.
   ! call ebacoli95_init(sol, npde_sub, (/xa,xb/), atol=atol, rtol=rtol, nint_max=10000)
    call ebacoli95_init(sol, npde_sub, xinit, atol=atol, rtol=rtol, nint_max=10000)

    allocate(uout(npde,sol%nint_max+1), stat=ier)
    if (ier /= 0 .or. sol%idid == -1000) goto 700
    tstart = sol%t0

    !-------------------------------------------------------------------
    ! Open files for output, 1 file per field
    do k = 1, npde
        write(npde_str,*) k
        fname = 'Points' // adjustl(trim(npde_str))
        open(unit=10+k,file=fname)
    end do

    !-------------------------------------------------------------------
    ! Integrate solution from t=0 to t=tout.
    print '(/"THE SOLVE PARAMETERS ARE")'
    print 900, sol%kcol, sol%nint, sol%npde, tstop
    print 901, sol%atol(1), sol%rtol(1), "LOI"

    do j = 2, ntout
        tout = tstart + (j-1)*(tstop-tstart)/(ntout-1)

        ! ebacoli call signature:
!        subroutine ebacoli95(sol, tout, f, bndxa, bndxb, uinit, &
!                derivf, difbxa, difbxb, tstop, nsteps)
         call ebacoli95(sol, tout, f, bndxa, bndxb, uinit, derivf, difbxa, difbxb)
        if (sol%idid <= 0) goto 800

        if (j == 2) then
            do i = 1, sol%nint+1
                call uinit(sol%x(i), uout(1,i), npde)
            end do
            do k = 1, npde
                do i = 1, sol%nint+1
                    write(10+k,*) sol%x(i), tstart, uout(k,i)
                end do
            end do
        end if

        call ebacoli95_vals(sol, sol%x(1:sol%nint+1), uout)

        do k = 1, npde
            do i = 1, sol%nint+1
                write(10+k,*) sol%x(i), sol%t0, uout(k,i)
            end do
        end do
    end do

    print '("IDID         = ",i2)', sol%idid
    print '("nsteps       = ",i10)', sol%num_accepted_time_steps
    print '("nint (tstop) = ",i10)', sol%nint

    !-------------------------------------------------------------------
    ! Cleanup and terminate
    call ebacoli95_sol_teardown(sol)
    stop

    !-------------------------------------------------------------------
    ! Errors
700 print '("Error: Could not allocate storage")' ; error stop
800 print '("Error: Was not able to integrate to tstop")' ; error stop

    !-------------------------------------------------------------------
    ! Formats!
900 format("kcol = ",i2,", nint0 = ",i4,", npde = ",i3,", tout = ",es7.1)
901 format("atol = ",es7.1,", rtol = ",es7.1,",",17x,a3)

end program
