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

subroutine read_concrete_parameters(inFile)
! Subroutine for reading concrete rewetting parameters
  character(len=*), intent(in) :: inFile

  integer, parameter :: Npar = 40
  character value_str(Npar)
  character(120) line

  integer i, j, stat, pos, posp1
  double precision par(Npar)

  ! parameters
  double precision cg01 ,cg02 ,cg03 ,cg04 ,cg05 ,cg06 ,cg07 ,cg08 ,cg09 ,cg10 &
         ,cg11 ,cg12 ,cg13 ,cg14 ,cg15 ,cg16 ,cg17 ,cg18 ,cg19 ,cg20 &
         ,cg21 ,cg22 ,cg23 ,cg24 ,cg25 ,cg26 ,cg27 ,cg28 ,cg29 ,cg30 &
         ,cg31 ,cg32 ,cg33 ,cg34 ,cg35 ,cg36 ,cg37 ,cg38 ,cg39 ,cg40

  common /concrete/ cg01 ,cg02 ,cg03 ,cg04 ,cg05 ,cg06 ,cg07 ,cg08 ,cg09 ,cg10 &
         ,cg11 ,cg12 ,cg13 ,cg14 ,cg15 ,cg16 ,cg17 ,cg18 ,cg19 ,cg20 &
         ,cg21 ,cg22 ,cg23 ,cg24 ,cg25 ,cg26 ,cg27 ,cg28 ,cg29 ,cg30 &
         ,cg31 ,cg32 ,cg33 ,cg34 ,cg35 ,cg36 ,cg37 ,cg38 ,cg39 ,cg40

  open(unit=10, file=inFile, status='old', iostat=stat)
  if (stat .ne. 0)then
    write (*,*) inFile, ' can not be opened'
    go to 99
  end if

  read(10,*) ! discard header line
  do i = 1, Npar
     read(10,'(a)') line
     call StripSpaces(line)
     ! convert the fourth column into double precision number
     pos = 0
     do j = 1,3
        pos = pos+index(line(pos+1:),',')
     end do
     posp1 = pos+index(line(pos+1:),',')
     read(line(pos+1:posp1-1), *) par(i)
  end do

  cg01 = par(1)
  cg02 = par(2)
  cg03 = par(3)
  cg04 = par(4)
  cg05 = par(5)
  cg06 = par(6)
  cg07 = par(7)
  cg08 = par(8)
  cg09 = par(9)
  cg10 = par(10)
  cg11 = par(11)
  cg12 = par(12)
  cg13 = par(13)
  cg14 = par(14)
  cg15 = par(15)
  cg16 = par(16)
  cg17 = par(17)
  cg18 = par(18)
  cg19 = par(19)
  cg20 = par(20)
  cg21 = par(21)
  cg22 = par(22)
  cg23 = par(23)
  cg24 = par(24)
  cg25 = par(25)
  cg26 = par(26)
  cg27 = par(27)
  cg28 = par(28)
  cg29 = par(29)
  cg30 = par(30)
  cg31 = par(31)
  cg32 = par(32)
  cg33 = par(33)
  cg34 = par(34)
  cg35 = par(35)
  cg36 = par(36)
  cg37 = par(37)
  cg38 = par(38)
  cg39 = par(39)
  cg40 = par(40)

! close files
  99 continue
  close(10)

end subroutine

subroutine StripSpaces(string)
  character(len=*) :: string
  integer :: stringLen
  integer :: last, actual

  stringLen = len (string)
  last = 1
  actual = 1

  do while (actual < stringLen)
     if (string(last:last) == ' ') then
        actual = actual + 1
        string(last:last) = string(actual:actual)
        string(actual:actual) = ' '
     else
        last = last + 1
        if (actual < last) &
             actual = last
     endif
  end do

end subroutine StripSpaces

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
    integer,  parameter :: npde = 5, nu = 4, nout = 11, nint_max = 50000, nin=1000

    ! Set (problem dependent) output time to 1
    real(dp), parameter :: tout = 0.0000001D0

    ! Declare output points and output solution values arrays
    real(dp) :: xin(nin), xout(nout), uout(npde*nout)
    real(dp) :: xa, xb

    integer :: i, k, ij

    external f, bndxa, bndxb, uinit, derivf, difbxa, difbxb

    ! parameters
    double precision cg01 ,cg02 ,cg03 ,cg04 ,cg05 ,cg06 ,cg07 ,cg08 ,cg09 ,cg10 &
         ,cg11 ,cg12 ,cg13 ,cg14 ,cg15 ,cg16 ,cg17 ,cg18 ,cg19 ,cg20 &
         ,cg21 ,cg22 ,cg23 ,cg24 ,cg25 ,cg26 ,cg27 ,cg28 ,cg29 ,cg30 &
         ,cg31 ,cg32 ,cg33 ,cg34 ,cg35 ,cg36 ,cg37 ,cg38 ,cg39 ,cg40

    common /concrete/ cg01 ,cg02 ,cg03 ,cg04 ,cg05 ,cg06 ,cg07 ,cg08 ,cg09 ,cg10 &
         ,cg11 ,cg12 ,cg13 ,cg14 ,cg15 ,cg16 ,cg17 ,cg18 ,cg19 ,cg20 &
         ,cg21 ,cg22 ,cg23 ,cg24 ,cg25 ,cg26 ,cg27 ,cg28 ,cg29 ,cg30 &
         ,cg31 ,cg32 ,cg33 ,cg34 ,cg35 ,cg36 ,cg37 ,cg38 ,cg39 ,cg40

    call read_concrete_parameters('concrete_parameters.csv')

    xa          = 0.D0
    ! xb          = 0.01D0
    xb          = cg34

    ! Write out value of npde to allow user to confirm that its value
    ! is appropriate for the problem to be solved.
    write(6,*) 'The total number of equations is ', npde
    write(6,*) 'The size of the parabolic subsystem is ', nu
    write(6,*)

    ! Initialize a grid to pass to ebacoli
    xin = (/(xa+i*(xb-xa)/(nin-1), i=0,nin-1)/)

    ! Initialization (set spatial domain = [0,1]); a default uniform
    ! spatial mesh having 10 subintervals will be constructed.
    call ebacoli95_init(sol, npde, nu, xin, atol=(/1d-2/), &
         rtol=(/1d-2/), nint_max = nint_max)

    ! Compute solution at tout
    call ebacoli95(sol, tout, f, bndxa, bndxb, uinit, derivf, difbxa, difbxb)

    ! Output idid to check for a successful computation
    print '("idid=",i5)',sol%idid
    if (sol%idid > 0) print '("idid > 0 => Successful computation")'

    ! Output solution at tout for nout values of x uniformly
    ! distributed over spatial domain
    if (sol%idid > -10) then
        xout = (/(xa+i*(xb-xa)/(nout-1), i=0,nout-1)/)
        call ebacoli95_vals(sol, xout, uout)

        print '("At t = ",f7.2)', sol%t0
        write(*,'(/a)') 'the solution is'
        write(*,'(a13,a27)') 'XOUT', 'UOUT'
        do i = 1, nout
           ij = (i-1)*npde
           write(*,*) xout(i), (uout(ij+k), k = 1, npde)
        end do
    end if

    call ebacoli95_sol_teardown(sol)

end program
