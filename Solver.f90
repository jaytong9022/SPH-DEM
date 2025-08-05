! BiCGSTAB�����
subroutine BICG(n, nnz, val, row, col, b, x)

integer, parameter :: max_iter = 1000  ! ����������
DOUBLE PRECISION, parameter :: tol = 1.0d-6    ! �����ݲ�
integer, intent(in) :: n, nnz
integer, intent(in) :: row(nnz), col(nnz)
DOUBLE PRECISION, intent(in) :: val(nnz), b(n)
DOUBLE PRECISION, intent(inout) :: x(n)
DOUBLE PRECISION, allocatable :: r(:), r_tilde(:), p(:), v(:), s(:), t(:)
DOUBLE PRECISION :: alpha, beta, omega, rho_old, rho_new, resid, b_norm, temp
integer :: iter
    
allocate(r(n), r_tilde(n), p(n), v(n), s(n), t(n))
    
! �����ʼ�в� r = b - A*x
call coo_matrix_vector_multiply(n, nnz, row, col, val, x, v)
r = b - v
    
! ��ʼ��r_tildeΪr
r_tilde = r
    
! �����ʼ�в����b�ķ���
b_norm = norm2(b)
resid = norm2(r) / b_norm
    
! �����ʼ�в��Ѿ���С��ֱ�ӷ���
if (resid <= tol) then
    iter = 0
    deallocate(r, r_tilde, p, v, s, t)
    return
end if
    
! ��ʼ������
rho_old = 1.0
alpha = 1.0
omega = 1.0
p = 0.0
    
! BiCGSTAB��ѭ��
do iter = 1, max_iter
    rho_new = dot_product(r_tilde, r)
        
    ! �����������
    if (rho_new == 0.0 .or. omega == 0.0) then
        exit
    end if
        
    beta = (rho_new / rho_old) * (alpha / omega)
    p = r + beta * (p - omega * v)
        
    ! ���� v = A*p
    call coo_matrix_vector_multiply(n, nnz, row, col, val, p, v)
        
    alpha = rho_new / dot_product(r_tilde, v)
        
    ! ���� s = r - alpha*v
    s = r - alpha * v
        
    ! �������
    resid = norm2(s) / b_norm
    if (resid <= tol) then
        x = x + alpha * p
        exit
    end if
        
    ! ���� t = A*s
    call coo_matrix_vector_multiply(n, nnz, row, col, val, s, t)
        
    ! ���� omega
    omega = dot_product(t, s) / dot_product(t, t)
        
    ! ���½�����
    x = x + alpha * p + omega * s
        
    ! ���²в�
    r = s - omega * t
        
    ! �������
    resid = norm2(r) / b_norm
    if (resid <= tol) then
        exit
    end if
        
    rho_old = rho_new
        
        
    ! ʹ���������ı���
    temp = dot_product(r_tilde, r)  ! <-- ����temp����ȷ����
        
    ! ����Ƿ�ﵽ����������
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



