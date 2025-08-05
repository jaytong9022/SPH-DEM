module Particle_System
  use omp_lib
  use Constants
  implicit none

  integer, parameter :: MAX_PARTICLES_PER_CELL = 100
  integer, parameter :: MAX_NEIGHBORS = 30
  
  type :: Box
    real(8) :: x_min(3)
    real(8) :: x_max(3)
  end type Box

  type NeighborPair
    integer :: j
    real(8) :: r_vec(3)
    real(8) :: r
    real(8) :: r2
  end type NeighborPair

  type Nlist
    integer :: ncount
    type(NeighborPair) :: pair(MAX_NEIGHBORS)
  end type Nlist
  
  ! type = 0 for dummy, 1 for fluid, 2 for wall boundary, 3 for DEM paricles
  type :: Particle
    real(8) :: x(3)
    real(8) :: v(3)
    real(8) :: rho_0
    real(8) :: rho_star
    real(8) :: a_visc(3)  ! 粘性加速度
    real(8) :: a_p(3)     ! 压力加速度
    integer :: p_type  
    type(Nlist) :: neighbor_list
  end type Particle

  type :: Cell
    integer :: count = 0
    integer :: particle_indices(MAX_PARTICLES_PER_CELL)
  end type Cell
  
  type :: ParticleSystem
    real(8) :: hcell
    type(Box) :: domain
    integer :: n_particles = 0
    integer :: capacity = 0
    integer :: key_lim(3)
    integer :: key_max
    integer, allocatable :: key_diff(:)

    type(Particle), allocatable :: particles(:)
    type(Cell), allocatable :: cell_list(:)
    type(Cell) :: removal_cell
  end type ParticleSystem
  
contains

  subroutine init_particle_system(ps_, hcell_, domain_)
  type(ParticleSystem), intent(out) :: ps_
  real(8), intent(in) :: hcell_
  type(Box), intent(in) :: domain_

  integer :: i, j, k, count
  real(8) :: x_min(3), x_max(3)

  ps_%hcell = hcell_
  ps_%domain = domain_
  x_min = domain_%x_min
  x_max = domain_%x_max
  is2D = (x_max(3) - x_min(3) == 0.0d0)

  ! 使用 ceiling，避免 ceil 报错
  do i = 1, 3
    ps_%key_lim(i) = ceiling((x_max(i) - x_min(i)) / hcell_)
  end do

  if (is2D) ps_%key_lim(3) = 1

  ps_%key_max = ps_%key_lim(1) * ps_%key_lim(2) * ps_%key_lim(3)

  ! 邻域偏移 key_diff
  if (is2D) then
    allocate(ps_%key_diff(9))
    count = 0
    do j = -1, 1
      do i = -1, 1
        count = count + 1
        ps_%key_diff(count) = i + ps_%key_lim(1) * j
      end do
    end do
  else
    allocate(ps_%key_diff(27))
    count = 0
    do k = -1, 1
      do j = -1, 1
        do i = -1, 1
          count = count + 1
          ps_%key_diff(count) = i + ps_%key_lim(1) * (j + ps_%key_lim(2) * k)
        end do
      end do
    end do
  end if

  allocate(ps_%cell_list(ps_%key_max))
  end subroutine init_particle_system
  
!--------------------------------------------------------------------------------------------! 
  subroutine build_cell_list(ps_)
  type(ParticleSystem), intent(inout) :: ps_
  integer :: i, key, tid, nthreads
  integer :: j, k, idx2
  integer, allocatable :: local_counts(:,:)
  integer, allocatable :: local_lists(:,:,:)
  logical, parameter :: debug_mode = .false.

  nthreads = omp_get_max_threads()
  
  ! 初始化主cell
  do i = 1, ps_%key_max
    ps_%cell_list(i)%count = 0
  end do

  ! 分配线程私有cell buffer
  allocate(local_counts(ps_%key_max, nthreads))
  allocate(local_lists(MAX_PARTICLES_PER_CELL, ps_%key_max, nthreads))
  local_counts = 0

  !$omp parallel private(i, key, tid)
    tid = omp_get_thread_num() + 1  ! Fortran索引从1开始

  !$omp do schedule(static)
    do i = 1, size(ps_%particles)
      ! 防止非法 i（理论上不会发生，但防御式编程）
      if (i < 1 .or. i > ps_%n_particles) then
        if (debug_mode) print *, "WARNING: Invalid i =", i
        cycle
      end if

      key = find_key(ps_, ps_%particles(i)%x)
      if (key > 0 .and. key <= ps_%key_max) then
        if (local_counts(key, tid) < MAX_PARTICLES_PER_CELL) then
          local_counts(key, tid) = local_counts(key, tid) + 1
          local_lists(local_counts(key, tid), key, tid) = i
        else if (debug_mode) then
          print *, "WARNING: local cell overflow at cell", key, " thread", tid
        end if
      end if
    end do
  !$omp end do
  !$omp end parallel

  ! 汇总线程私有cell数据到全局
  do j = 1, ps_%key_max
    ps_%cell_list(j)%count = 0
    do tid = 1, nthreads
      do k = 1, local_counts(j, tid)
        idx2 = local_lists(k, j, tid)

        ! 防止非法索引
        if (idx2 < 1 .or. idx2 > ps_%n_particles) then
          if (debug_mode) then
            print *, "ERROR: Invalid particle idx2 =", idx2, "in cell", j, "from thread", tid
          end if
          cycle
        end if

        if (ps_%cell_list(j)%count < MAX_PARTICLES_PER_CELL) then
          ps_%cell_list(j)%count = ps_%cell_list(j)%count + 1
          ps_%cell_list(j)%particle_indices(ps_%cell_list(j)%count) = idx2
        else
          print *, "WARNING: Cell ", j, " exceeds max capacity during merge."
          exit
        end if
      end do
    end do
  end do

  deallocate(local_counts, local_lists)
  end subroutine build_cell_list
  
!--------------------------------------------------------------------------------------------! 
  subroutine add_particle(ps_, p_)
    type(ParticleSystem), intent(inout) :: ps_
    type(Particle), intent(in) :: p_
    integer :: key, i, j, neighbor_key, pid
    logical :: is_duplicate
    real(8), parameter :: eps2 = 1.0d-20

    key = find_key(ps_, p_%x)
    if (key <= 0 .or. key > ps_%key_max) return

    is_duplicate = .false.
    do j = 1, size(ps_%key_diff)
        neighbor_key = key + ps_%key_diff(j)
        if (neighbor_key < 1 .or. neighbor_key > ps_%key_max) cycle
        do i = 1, ps_%cell_list(neighbor_key)%count
            pid = ps_%cell_list(neighbor_key)%particle_indices(i)
            if (sum((ps_%particles(pid)%x - p_%x)**2) < eps2) then
                is_duplicate = .true.
                exit
            end if
        end do
        if (is_duplicate) exit
    end do

    if (is_duplicate) return

    if (.not. allocated(ps_%particles)) then
        ps_%capacity = 1024
        allocate(ps_%particles(ps_%capacity))
    else if (ps_%n_particles >= ps_%capacity) then
        call resize_particles(ps_, ps_%capacity * 2)
    end if

    ps_%n_particles = ps_%n_particles + 1
    ps_%particles(ps_%n_particles) = p_

    if (ps_%cell_list(key)%count >= MAX_PARTICLES_PER_CELL) then
      print *, "Error: cell ", key, " overflow! Too many particles."
    stop
    end if

  ps_%cell_list(key)%count = ps_%cell_list(key)%count + 1
  ps_%cell_list(key)%particle_indices(ps_%cell_list(key)%count) = ps_%n_particles
  end subroutine add_particle

  subroutine resize_particles(ps_, new_capacity)
      type(ParticleSystem), intent(inout) :: ps_
      integer, intent(in) :: new_capacity
      type(Particle), allocatable :: tmp(:)

      ! 为粒子数组分配新空间，并将已有粒子移动过去
      allocate(tmp(new_capacity))
      if (ps_%n_particles > 0) tmp(1:ps_%n_particles) = ps_%particles(1:ps_%n_particles)
      call move_alloc(tmp, ps_%particles)
      ps_%capacity = new_capacity
  end subroutine resize_particles


  subroutine generate_orthogonal_distribution(ps_, dx_, x0, x1, type_id)
  type(ParticleSystem), intent(inout) :: ps_
  real(8), intent(in) :: dx_
  real(8), intent(in) :: x0(3), x1(3)
  integer, intent(in) :: type_id

  real(8) :: dy_, xpos, ypos, zpos
  integer :: i, j
  type(Particle) :: p_

  dy_ = dx_ 
  j = 0
  ypos = x0(2)
  zpos = x0(3)
  do while (ypos <= x1(2))
    xpos = x0(1)
    do while (xpos <= x1(1))
      p_%x = (/ xpos, ypos, zpos /)
      p_%p_type = type_id
      call add_particle(ps_, p_)
      xpos = xpos + dx_
    end do
      ypos = ypos + dy_
      j = j + 1
  end do
  end subroutine generate_orthogonal_distribution
  
  subroutine generate_polyline_particles(ps_, dx_, vertices, n_vert, type_id)
  type(ParticleSystem), intent(inout) :: ps_
  real(8), intent(in) :: dx_
  real(8), intent(in) :: vertices(3, n_vert)  ! 每列一个点
  integer, intent(in) :: n_vert
  integer, intent(in) :: type_id

  integer :: i, j, n_seg_points
  real(8) :: len, dir(3), p0(3), p1(3), pos(3)
  type(Particle) :: p_

  ! 如果顶点数小于2，直接返回
  if (n_vert < 2) return

  ! 遍历每对相邻的顶点
  do i = 1, n_vert - 1
    p0 = vertices(:, i)
    p1 = vertices(:, i+1)

    ! 计算两点之间的方向向量和距离
    dir = p1 - p0
    len = sqrt(sum(dir**2))

    ! 如果长度为0，则跳过该对点
    if (len < 1.0d-12) cycle  ! 增加一个阈值，避免计算无效点对

    ! 归一化方向向量
    dir = dir / len

    ! 根据长度计算需要的分段数
    n_seg_points = max(1, int(len / dx_ + 0.5))  ! 四舍五入，避免过小分段数

    ! 生成粒子
    do j = 0, n_seg_points
        pos = p0 + dx_ * j * dir
        p_%x = pos
        p_%p_type = type_id
        call add_particle(ps_, p_)
    end do
  end do
  end subroutine generate_polyline_particles
  
 !--------------------------------------------------------------------------------------------- 
  subroutine build_neighbor_list(ps_)
  type(ParticleSystem), intent(inout) :: ps_

  integer :: i, j, k, idx, idx2, nk, dim, key, n, pid, nid
  real(8) :: r_vec(3), r, r2, r_temp

  dim = merge(2, 3, is2D)

  !$omp parallel do private(i,j,k,idx2,nk,key,r_vec,r2,grad_w) schedule(static)
  do idx = 1, ps_%n_particles
    ps_%particles(idx)%neighbor_list%ncount = 0

    key = find_key(ps_, ps_%particles(idx)%x)
    if (key < 1 .or. key > ps_%key_max) cycle

    do k = 1, size(ps_%key_diff)
      j = key + ps_%key_diff(k)
      if (j < 1 .or. j > ps_%key_max) cycle

      do i = 1, ps_%cell_list(j)%count
        idx2 = ps_%cell_list(j)%particle_indices(i)
        if (idx2 == idx) cycle

        ! 计算 r_vec, r², r
        r_vec(1:dim) = ps_%particles(idx2)%x(1:dim) - ps_%particles(idx)%x(1:dim)
        r2 = sum(r_vec(1:dim)**2)
        r = sqrt(r2)
        if(r <= 2.0d0*h) then 
          nk = ps_%particles(idx)%neighbor_list%ncount + 1
          if (nk <= MAX_NEIGHBORS) then
            ps_%particles(idx)%neighbor_list%ncount = nk
            ps_%particles(idx)%neighbor_list%pair(nk)%j = idx2
            ps_%particles(idx)%neighbor_list%pair(nk)%r_vec = r_vec
            ps_%particles(idx)%neighbor_list%pair(nk)%r = r
            ps_%particles(idx)%neighbor_list%pair(nk)%r2 = r2
          end if
        endif
      end do
    end do
  end do
  !$omp end parallel do
  
  open(unit=99, file='neigbourlist.csv', status='replace', action='write')

    write(99, '(A)') 'particle_id, neighbor_id, distance'

    do i = 1, ps_%n_particles
        n = ps_%particles(i)%neighbor_list%ncount
        do j = 1, n
            pid = i
            nid = ps_%particles(i)%neighbor_list%pair(j)%j
            r_temp = ps_%particles(i)%neighbor_list%pair(j)%r
            write(99, '(I10,",", I10,",", F15.8)') pid, nid, r_temp
        end do
    end do

    close(99)
  end subroutine build_neighbor_list

!--------------------------------------------------------------------------------------------!
  function find_key(ps_, x) result(key)
    type(ParticleSystem), intent(in) :: ps_
    real(8), intent(in) :: x(3)
    integer :: key, i, j, k

    ! 使用相对于 domain%x_min 的坐标定位 cell，边界粒子归入最后一个cell
    i = min(ps_%key_lim(1), 1 + int((x(1) - ps_%domain%x_min(1)) / ps_%hcell))
    j = min(ps_%key_lim(2), 1 + int((x(2) - ps_%domain%x_min(2)) / ps_%hcell))
    k = min(ps_%key_lim(3), 1 + int((x(3) - ps_%domain%x_min(3)) / ps_%hcell))

    if (i < 1 .or. j < 1 .or. k < 1 .or. &
        i > ps_%key_lim(1) .or. j > ps_%key_lim(2) .or. k > ps_%key_lim(3)) then
      key = -1
    else
      key = i + ps_%key_lim(1)*(j-1) + ps_%key_lim(1)*ps_%key_lim(2)*(k-1)
    end if
  end function find_key
  
end module Particle_System
