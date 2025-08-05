program TestParticleSystem
  use Particle_System
  use Constants
  use SPH_Formulation
  implicit none

  type(ParticleSystem) :: ps
  type(Box) :: sim_box
  type(Particle) :: p
  integer :: i, j, idx, np, key, idx2
  real(8) :: hcell, x_temp, y_temp, z_temp
  real(8) :: dist, dx_, dy_
  integer :: neighbor_count
  real(8) :: x0(3), x1(3)
  real(8) :: polyline(3, 4)

  polyline(:, 1) = (/ 0.0d0, h_tank, 0.0d0 /)
  polyline(:, 2) = (/ 0.0d0, 0.0d0, 0.0d0 /)
  polyline(:, 3) = (/ l_tank, 0.0d0, 0.0d0 /)
  polyline(:, 4) = (/ l_tank, h_tank, 0.0d0 /)

  ! |                                                                       |
  ! | ------------|                                                         |
  ! |             |       DAMBREAK MODEL TANK                               |
  ! |             |                                                         |  
  ! |             |                                                         |
  ! |             |                                                         |
  ! |-----------------------------------------------------------------------|
  
  hcell = 2.0*h
  neighbor_count = 0
  sim_box%x_min = (/-5*h, -5*h, 0.0d0/)
  sim_box%x_max = (/l_tank+2*h, h_tank+2*h, 0.0d0/)
  x0 = (/ dx, dx, 0.0d0/)
  x1 = (/ l_water, h_water, 0.0d0 /)
  
  call init_particle_system(ps, hcell, sim_box)
  call generate_polyline_particles(ps, dx, polyline, 4, 2)   ! 2 for wall boundary particles
  call generate_orthogonal_distribution(ps, dx, x0, x1, 1)   ! 1 for fluid particles
  call build_cell_list(ps)  
  call build_neighbor_list(ps)
  call compute_density(ps)
  
  open(unit=20, file='density.csv', status='replace')
  write(20, '(A)') 'id,x,y,z,rho,type'
  do i = 1, ps%n_particles
    write(20, '(I4,",",F12.6,",",F12.6,",",F12.6,","F12.6,",",I4)') i, ps%particles(i)%x(1), ps%particles(i)%x(2),ps%particles(i)%x(3), ps%particles(i)%rho_0, ps%particles(i)%p_type
  end do
  close(20)

print *, "Particle densities written to 'density.csv'"
    
end program TestParticleSystem
