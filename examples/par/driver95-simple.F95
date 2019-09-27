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

#ifndef MAC_NU
#define MAC_NU 0
#endif
#ifndef MAC_NV
#define MAC_NV 0
#endif
#ifndef MAC_NW
#define MAC_NW 0
#endif

program simple_example_driver

    ! A minimal example driver for the EBACOLI95 wrapper of EBACOLI.
    ! See the comments within ebacoli95.f95 for details regarding use.

    use ebacoli95_mod, only: ebacoli95_init, ebacoli95, ebacoli95_vals, &
                             ebacoli95_sol, ebacoli95_sol_teardown

    implicit none
    integer, parameter  :: dp = kind(0d0)

    ! Declare data structure that will hold solution
    type(ebacoli95_sol)  :: sol

    ! Choose npde to be consistent with the problem specific source
    ! code; see burg1.f, burg2.f, CahnAllen.f, rcdsys.f, sinmads.f,
    ! and steady.f.
    ! Assume output at 11 points over (problem specific) spatial domain,
    ! and maximum number of subintervals = 50.
    integer,  parameter, dimension(3) :: npde_sub = (/MAC_NU,MAC_NV,MAC_NW/)
    integer,  parameter :: npde = sum(npde_sub), nout = 11, nint_max = 500

    ! Set (problem dependent) output time to 1
    real(dp), parameter :: tout = 1

    ! Declare output points and output solution values arrays
    real(dp) :: xout(nout), uout(npde*nout)

    integer :: i, k, ij

    external f, bndxa, bndxb, uinit, derivf, difbxa, difbxb

    ! eps is used in the burg1.f, burg2.f and Cahn_Allen.f examples.
    real(dp) :: eps
    common /burger/ eps
    ! The following parameters are used in the rcdsys.f example.
    real(dp) :: a1,a2,d1,d2,r,c,n,pe1,pe2
    common /rcdstuff/ a1,a2,d1,d2,r,c,n,pe1,pe2

!-------------------------------------------------------------------------
    ! Set problem dependent parameters
    eps = 1e-3_dp
    a1 = 30._dp
    a2 = 30._dp
    d1 = 1.5_dp
    d2 = 1.2_dp
    r = 1000._dp
    c = .96_dp
    n = 1._dp
    pe1 = 1e4_dp
    pe2 = 1e4_dp

    ! Write out value of npde to allow user to confirm that its value
    ! is appropriate for the problem to be solved.

    write(6,*) 'The number of PDEs is assumed to be ', npde
    write(6,*)

    ! Initialization (set spatial domain = [0,1]); a default uniform
    ! spatial mesh having 10 subintervals will be constructed.
    call ebacoli95_init(sol, npde_sub, (/0._dp,1._dp/), nint_max = nint_max)

    ! Compute solution at tout
    call ebacoli95(sol, tout, f, bndxa, bndxb, uinit, derivf, difbxa, difbxb)

    ! Output idid to check for a successful computation
    print '("idid=",i5)',sol%idid
    if (sol%idid > 0) print '("idid > 0 => Successful computation")'

    ! Output solution at tout for nout values of x uniformly
    ! distributed over spatial domain
    if (sol%idid > 0) then
        xout = (/(i*0.1_dp, i=0,nout-1)/)
        call ebacoli95_vals(sol, xout, uout)

        print '("At t=",f4.2)', sol%t0
        write(*,'(/a)') 'the solution is'
        write(*,'(a13,a27)') 'XOUT', 'UOUT'
        do i = 1, nout
           ij = (i-1)*npde
           write(*,*) xout(i), (uout(ij+k), k = 1, npde)
        end do
    end if

    call ebacoli95_sol_teardown(sol)

end program
