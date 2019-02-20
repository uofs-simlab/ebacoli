program simple_example_driver

    ! A minimal example driver for the eBACOLI95 wrapper of eBACOLI.
    use ebacoli95_mod, only: ebacoli95_init, ebacoli95, ebacoli95_vals, &
                             ebacoli95_sol, ebacoli95_sol_teardown
    implicit none
    integer, parameter :: dp = kind(0d0)

    ! Declare data structure that will hold solution
    type(ebacoli95_sol) :: sol

    ! Choose npde to be consistent with the problem specific source
    ! Assume output at 11 points over (problem specific) spatial domain,
    ! and select a maximum number of subintervals.
    integer,  parameter :: npde = 4, nu = 2, nout = 11, nint_max = 500
    real(dp), parameter :: xa = 1, xb = 2

    ! Set (problem dependent) output time to 1
    real(dp), parameter :: tout = 1

    ! Declare output points and output solution values arrays
    real(dp) :: xout(nout), uout(npde*nout)

    integer :: i, k, ij

    ! Fortran77 subroutines defining the problem
    external f, bndxl, bndxr, uinit, derivf, difbxl, difbxr

    ! Initialization (set spatial domain = [xa,xb]); a default uniform
    ! spatial mesh having 10 subintervals will be constructed.
    call ebacoli95_init(sol, npde, nu, (/xa,xb/), nint_max = nint_max)

    ! Compute solution at tout
    call ebacoli95(sol, tout, f, bndxl, bndxr, uinit, derivf, difbxl, difbxr)

    ! Output idid to check for a successful computation
    print '("idid=",i5)',sol%idid
    if (sol%idid > 0) print '("idid > 0 => Successful computation")'

    ! Output solution at tout for nout values of x uniformly
    ! distributed over spatial domain
    if (sol%idid > 0) then
        xout = (/(xa+i*(xb-xa)/(nout-1), i=0,nout-1)/)
        call ebacoli95_vals(sol, xout, uout)

        print '("At t=",f4.2)', sol%t0
        write(*,'(/a)') 'the solution is'
        write(*,'(a13,a27)') 'XOUT', 'UOUT'
        do i = 1, nout
           ij = (i-1)*npde
           write(*,*) xout(i), (uout(ij+k), k = 1, npde)
        end do
    end if

    ! Clean up ebacoli
    call ebacoli95_sol_teardown(sol)

end program
