MODULE SPH_Formulation
  use Constants
  use Particle_System
  use omp_lib              ! OpenMP 支持
  implicit none
  
contains
    
  subroutine compute_density(ps_)
  type(ParticleSystem), intent(inout) :: ps_

  integer :: i, j
  integer :: dim
  real(8) :: rho_i, rij

  dim = merge(2, 3, is2D)

  !$omp parallel do private(i, j, rho_i, rij) schedule(static)
  do i = 1, ps_%n_particles
    rho_i = m * cubic_spline(0.0d0)
    do j = 1, ps_%particles(i)%neighbor_list%ncount
      rij = ps_%particles(i)%neighbor_list%pair(j)%r
      rho_i = rho_i + m * cubic_spline(rij)
    end do
    ps_%particles(i)%rho_0 = rho_i
  end do
  !$omp end parallel do
  end subroutine compute_density

!----------------------------------------------------------------------------------------
subroutine compute_viscosity_acceleration(ps_)
  type(ParticleSystem), intent(inout) :: ps_

  integer :: i, j, d
  integer :: dim, idx2
  real(8) :: r2, dot, rho_avg
  real(8) :: eta2, v_diff(3)
  real(8) :: vrec_temp(3), wg(3)

  eta2 = (0.01d0 * h)**2
  dim = merge(2, 3, is2D)

  !$omp parallel do private(i,j,d,idx2,r2,dot,rho_avg,v_diff) schedule(static)
  do i = 1, ps_%n_particles
    ps_%particles(i)%a_visc = 0.0d0

    do j = 1, ps_%particles(i)%neighbor_list%ncount
      idx2 = ps_%particles(i)%neighbor_list%pair(j)%j
      r2   = ps_%particles(i)%neighbor_list%pair(j)%r2
      vrec_temp = ps_%particles(i)%neighbor_list%pair(j)%r_vec
      wg = cubic_spline_grad(vrec_temp)
      
      ! v_diff = v_j - v_i
      v_diff(1:dim) = ps_%particles(idx2)%v(1:dim) - ps_%particles(i)%v(1:dim)
      
      ! dot(r_ij, ∇W)
      dot = sum(vrec_temp(1:dim) * wg(1:dim))

      rho_avg = 0.5d0 * (ps_%particles(i)%rho_0 + ps_%particles(idx2)%rho_0)
      if (rho_avg > 0.0d0 .and. r2 + eta2 > 0.0d0) then
        do d = 1, dim
          ps_%particles(i)%a_visc(d) = ps_%particles(i)%a_visc(d) + &
               m * (4.0d0 * vis_eff / rho_avg) * v_diff(d) * dot / (r2 + eta2)
        end do
      end if
    end do
  end do
  !$omp end parallel do
  end subroutine compute_viscosity_acceleration


end module SPH_Formulation