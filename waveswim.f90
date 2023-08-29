! waveswim.f90
! Integration of the equations for an elongated, gyrotactic swimmer in the velocity field produced by a linear wave.
! 2D case, equations as on
! Ventrella et al. "Microswimmer trapping in surface waves with shear", eqs 2.4
! variables Y(1)->x, Y(2)->z, Y(3)->phi
! Input file: input.dat
! Output file: girotax.dat
! Author: Francesco Ventrella (2023)

module global
  implicit none
  save
  !Integration parameters
  real*8  :: ti , tf , h
  !phisical parameters
  real*8  :: a , lamda, alfa, nu, Psi, xo, zo, phio, zh, hz, vg, beta, sigma
  real*8, parameter :: pi = 3.14159265358979d0, g = 9.81d0
  real*8 :: k, w
  integer :: M
end module

!MAIN PROGRAM
program girotax
  use global
  implicit none

  real*8  :: Y(3), t
  integer :: N, Neq=3, i, j, iout

  !read input parameters
  open(3, file = 'input.dat')
  read(3,*) ti !starting time
  read(3,*) tf !ending time
  read(3,*) h !integration time step
  read(3,*) beta !characteristic depth (negative value)
  read(3,*) a !wave's amplitude
  read(3,*) lamda !elongation parameter
  read(3,*) alfa !steepness
  read(3,*) nu !adimesional swimming speed
  read(3,*) xo !x initial condition
  read(3,*) zo !z initial condition
  read(3,*) phio !angular initial condition
  read(3,*) Psi !nondimensional gyrotactic time
  read(3,*) hz !step for initial condition of trajectories on vertical axis
  read(3,*) M !number of trajectories
  read(3,*) vg !settling velocity
  read(3,*) sigma !nondimensional shear parameter
  read(3,*) iout !prints swimmer positions every iou steps
  close(3)

  !wave's parameters
  k = alfa/a !wave number
  w = dsqrt(g*k) !frequency

  N = NINT((tf-ti)/h) !number of integration steps

  !output file
  !4 columns
  ! time, x, y,swimming angle
  open(1, file = 'girotax.dat')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do j=0,M,1

    !set initial conditions
    Y(1) = xo
    Y(2) = zo-hz*j
    Y(3) = phio
    t=ti

    do i=1,N

      call rk4(t,Y,Neq)
      t=ti+i*h

      !save X,Z,Phi
      if (mod(i,iout) .eq. 0) then
        write(1,'(4g20.8)') t,Y(1),Y(2),Y(3)
      endif

      ! logic control on surface z=0 and foundal z=beta
      if ( Y(2) > 0.d0 ) then
        exit
      endif
      if ( Y(2) < beta ) then
        exit
      endif

    end do

    write(1,*) NEW_LINE(' ') !for organizing data plotting

  end do

  close(1)

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine rhs(t,Y,R)
  !compute right hand side for rk4
  use global
  implicit none

  real*8 t, x, z, phi
  real*8 Y(3),R(3)

  x   = Y(1)
  z   = Y(2)
  phi = Y(3)

  !ref. to eq. (2.4)
  R(1) = alfa*dexp(z)*dcos(x-t) + nu*dsin(phi) + sigma*(dabs(beta)+z)
  R(2) = alfa*dexp(z)*dsin(x-t) + nu*dcos(phi) - vg
  R(3) = lamda*alfa*dexp(z)*dcos(x-t+2*phi) - 1.d0/(2.d0*Psi)*dsin(phi) + sigma*0.5d0*(1.d0+lamda*dcos(2.d0*phi))

end subroutine

subroutine rk4(x,Y,Neq)
  !Runge-Kutta algorithm fourth order
  use global
  implicit none

  real*8 x
  integer Neq,i
  real*8::Y(Neq)
  real*8::Yl(Neq),k1(Neq),K2(Neq),k3(Neq),k4(Neq)

  call rhs(x,Y,k1)

  do i=1,Neq
    Yl(i)=Y(i)+0.5d0*h*K1(i)
  end do

  call rhs(x+0.5d0*h,Yl,k2)

  do i=1,Neq
    Yl(i)=Y(i)+0.5d0*h*K2(i)
  end do

  call rhs(x+0.5d0*h,Yl,k3)

  do i=1,Neq
    Yl(i)=Y(i)+h*K3(i)
  end do

  call rhs(x+h,Yl,k4)

  do i=1,Neq
    Y(i)=Y(i)+h/6.d0*(k1(i)+2.d0*k2(i)+2.d0*k3(i)+k4(i))
  end do


end subroutine
