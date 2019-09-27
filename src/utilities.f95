! Utility functions

! ----------------------------------------------------------------------------
subroutine write_subsystem_sizes(npde_sub)

  ! Write the total number of PDEs and the size of each subsystem to stdout

  implicit none

  ! size of each type of subsystem
  integer npde_sub(3)

  write(6,*) 'Total number of PDEs: ', sum(npde_sub)
  write(6,*) '  Size of parabolic:', npde_sub(1)
  write(6,*) '  Size of ODE      :', npde_sub(2)
  write(6,*) '  Size of elliptic :', npde_sub(3)
  write(6,*)

end subroutine write_subsystem_sizes


! ----------------------------------------------------------------------------
subroutine write_solution_uniformly_spaced(sol,nout)

  ! Write solution (to stdout) on uniformly space grid

  use ebacoli95_mod, only: ebacoli95_sol, ebacoli95_vals

  implicit none

  type(ebacoli95_sol) :: sol
  integer nout

  ! Declare output points and output solution values arrays
  double precision :: xout(nout), uout(sol%npde*nout)

  integer i, ij, k

  ! Output solution at uniformly spaced coordinates
  xout = (/(sol%x(1)+i*(sol%x(sol%nint+1)-sol%x(1))/(nout-1), i=0,nout-1)/)
  call ebacoli95_vals(sol, xout, uout)

  print '("At t=",f4.2)', sol%t0
  write(*,'(/a)') 'the solution is'
  write(*,'(a13,a27)') 'XOUT', 'UOUT'
  do i = 1, nout
     ij = (i-1)*sol%npde
     write(*,*) xout(i), (uout(ij+k), k = 1, sol%npde)
  end do

end subroutine write_solution_uniformly_spaced


! ----------------------------------------------------------------------------
subroutine write_discrete_L2_error(sol,truu,npoints)

  ! Computes the discrete L2 error of the solution at npoints equally spaced
  ! points. Writes the result to screen

  use ebacoli95_mod, only: ebacoli95_sol, ebacoli95_vals

  implicit none

  type(ebacoli95_sol) :: sol
  integer npoints

  ! Declare output points and output solution values arrays
  double precision :: xout(npoints), uout(sol%npde*npoints), uactual(sol%npde*npoints)

  double precision error

  integer i, ij, k

  ! Compute the discrete L2 error relative to known solution truu(t,x)

  ! Solution at uniformly spaced coordinates
  xout = (/(sol%x(1)+i*(sol%x(sol%nint+1)-sol%x(1))/(npoints-1), i=0,npoints-1)/)
  call ebacoli95_vals(sol, xout, uout)

  error = 0
  do i = 1, npoints
     ij = 1+(i-1)*sol%npde
     call truu(sol%t0,xout(i),uactual(ij),sol%npde)
  end do

  error = norm2(uactual-uout)

  write(*,*) "Discrete L2 error at output points:", error

end subroutine write_discrete_L2_error
