! BiCGSTAB求解器
subroutine BICG(n, nnz, val, row, col, b, x)

integer, parameter :: max_iter = 1000  ! 最大迭代次数
DOUBLE PRECISION, parameter :: tol = 1.0d-6    ! 收敛容差
integer, intent(in) :: n, nnz
integer, intent(in) :: row(nnz), col(nnz)
DOUBLE PRECISION, intent(in) :: val(nnz), b(n)
DOUBLE PRECISION, intent(inout) :: x(n)
DOUBLE PRECISION, allocatable :: r(:), r_tilde(:), p(:), v(:), s(:), t(:)
DOUBLE PRECISION :: alpha, beta, omega, rho_old, rho_new, resid, b_norm, temp
integer :: iter
    
allocate(r(n), r_tilde(n), p(n), v(n), s(n), t(n))
    
! 计算初始残差 r = b - A*x
call coo_matrix_vector_multiply(n, nnz, row, col, val, x, v)
r = b - v
    
! 初始化r_tilde为r
r_tilde = r
    
! 计算初始残差范数和b的范数
b_norm = norm2(b)
resid = norm2(r) / b_norm
    
! 如果初始残差已经很小，直接返回
if (resid <= tol) then
    iter = 0
    deallocate(r, r_tilde, p, v, s, t)
    return
end if
    
! 初始化参数
rho_old = 1.0
alpha = 1.0
omega = 1.0
p = 0.0
    
! BiCGSTAB主循环
do iter = 1, max_iter
    rho_new = dot_product(r_tilde, r)
        
    ! 检查重启条件
    if (rho_new == 0.0 .or. omega == 0.0) then
        exit
    end if
        
    beta = (rho_new / rho_old) * (alpha / omega)
    p = r + beta * (p - omega * v)
        
    ! 计算 v = A*p
    call coo_matrix_vector_multiply(n, nnz, row, col, val, p, v)
        
    alpha = rho_new / dot_product(r_tilde, v)
        
    ! 计算 s = r - alpha*v
    s = r - alpha * v
        
    ! 检查收敛
    resid = norm2(s) / b_norm
    if (resid <= tol) then
        x = x + alpha * p
        exit
    end if
        
    ! 计算 t = A*s
    call coo_matrix_vector_multiply(n, nnz, row, col, val, s, t)
        
    ! 计算 omega
    omega = dot_product(t, s) / dot_product(t, t)
        
    ! 更新解向量
    x = x + alpha * p + omega * s
        
    ! 更新残差
    r = s - omega * t
        
    ! 检查收敛
    resid = norm2(r) / b_norm
    if (resid <= tol) then
        exit
    end if
        
    rho_old = rho_new
        
        
    ! 使用已声明的变量
    temp = dot_product(r_tilde, r)  ! <-- 现在temp已正确声明
        
    ! 检查是否达到最大迭代次数
    if (iter >= max_iter) then
        exit
    end if
end do
deallocate(r, r_tilde, p, v, s, t)
contains
    subroutine coo_matrix_vector_multiply(n, nnz, row, col, val, x, y)   
        integer, intent(in) :: n, nnz
        integer, intent(in) :: row(nnz), col(nnz)
        DOUBLE PRECISION, intent(in) :: val(nnz), x(n)
        DOUBLE PRECISION, intent(out) :: y(n)
        integer :: i, k
        
        y = 0.0
        
        do k = 1, nnz
            y(row(k)) = y(row(k)) + val(k) * x(col(k))
        end do
    end subroutine
end subroutine



