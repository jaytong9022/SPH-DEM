MODULE Constants
  implicit none
  private
  public :: is2D
  public :: cubic_spline, cubic_spline_grad
  public :: pi, h, dt, m, g, vis_eff, alpha, dx, h_tank, l_tank, h_water, l_water
  
  logical :: is2D = .true.
  real(8), parameter :: pi = 3.141592653589793d0
  real(8), parameter :: h = 0.012d0         ! Kernel distance
  real(8), parameter :: dt = 0.001d0        ! Time step
  real(8), parameter :: g = -9.81d0         ! Gravity acceleration
  real(8), parameter :: vis_eff = 1.0e-6    ! Laminar viscosity
  real(8), parameter :: alpha = 0.2d0       ! Courant coefficient
  real(8), parameter :: dx = 0.01d0         ! Particle initial spacing
  real(8), parameter :: m = 1000.d0*dx*dx   ! Particle mass
  real(8), parameter :: h_tank = 3.0d0
  real(8), parameter :: l_tank = 5.0d0
  real(8), parameter :: h_water = 1.0d0
  real(8), parameter :: l_water = 0.5d0
  
contains

function cubic_spline(r_) result(w)
  real(8), intent(in) :: r_
  real(8) :: w, q, sigma

  q = r_ / h

  if (is2D) then
    sigma = 10.0d0 / (7.0d0 * pi * h**2)
  else
    sigma = 15.0d0 / (pi * h**3)
  end if

  if (q >= 0.0d0 .and. q <= 1.0d0) then
    w = sigma * (1.0d0 - 1.5d0*q**2 + 0.75d0*q**3)
  else if (q > 1.0d0 .and. q <= 2.0d0) then
    w = sigma * 0.25d0 * (2.0d0 - q)**3
  else
    w = 0.0d0
  end if
end function cubic_spline

function cubic_spline_grad(r_vec) result(grad_w)
  real(8), intent(in) :: r_vec(3)
  real(8) :: grad_w(3)

  real(8) :: r, q, sigma, dw_dr
  integer :: dim, d

  r = sqrt(sum(r_vec**2))
  q = r / h
  dim = merge(2, 3, is2D)

  if (is2D) then
    sigma = 10.0d0 / (7.0d0 * pi * h**2)
  else
    sigma = 15.0d0 / (pi * h**3)
  end if

  if (q >= 0.0d0 .and. q <= 1.0d0) then
    dw_dr = sigma * (-3.0d0*q + 2.25d0*q**2) / h
  else if (q > 1.0d0 .and. q <= 2.0d0) then
    dw_dr = -0.75d0 * sigma * (2.0d0 - q)**2 / h
  else
    dw_dr = 0.0d0
  end if

  if (r < 1.0d-12) then
    grad_w = 0.0d0
  else
    do d = 1, dim
      grad_w(d) = dw_dr * r_vec(d) / r
    end do
    do d = dim+1, 3
      grad_w(d) = 0.0d0
    end do
  end if
end function cubic_spline_grad

end module Constants
