! Utility functions for testing


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
subroutine compare_L2_error(sol,truu)

  ! Compute the L2 error relative to known solution truu(t,x)

  write(*,*) 'L2 @ U. lolz'

end subroutine compare_L2_error
