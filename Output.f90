!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  THIS SUBROUTINE DOES NOT NEED TO BE CHANGED   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	VISUALIZATION
SUBROUTINE output(time_step)
USE mdata
!!IMPLICIT NONE
!!INTEGER::I
!!PARTICLE SNAPSHOTS
!!OPEN(2,FILE='xy')
!!REWIND(2)
!!    DO I=1,N
!!        WRITE(2,14)X1(I),Y1(I),IFLAG(I),P(I)
!!        14	FORMAT(2f10.5,I5,f10.5)
!!    ENDDO
!!CLOSE(2)
!!RETURN
implicit none
integer, intent(in) :: time_step
integer :: i, iunit, stat
character(len=256) :: filename, dirname, time_str
logical :: dir_exists

! 定义输出目录和文件名
dirname = 'output'

! 检查并创建 output 目录（如果不存在）
inquire(file=trim(dirname)//'/.', exist=dir_exists)  ! 检查目录是否存在
if (.not. dir_exists) then
    call system('mkdir "' // trim(dirname) // '" > nul 2>&1')  ! 屏蔽所有输出
end if

! 将时间步格式化为4位数字（如 5 → "0005"）
write(time_str, '(I6.6)') time_step

! 生成文件名（显式处理路径分隔符）
filename = trim(dirname) // '\sph_output_' // trim(time_str) // '.vtu'


! 打开文件（带错误处理）
open(newunit=iunit, file=filename, status='replace', action='write', iostat=stat)
if (stat /= 0) then
    print *, '错误: 无法创建文件 '//trim(filename)
    stop
end if

! 写入VTK头部
write(iunit, '(a)') '<?xml version="1.0" encoding="UTF-8"?>'  ! <<<
write(iunit, '(a)') '<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian">'
write(iunit, '(a)') '<UnstructuredGrid>'
write(iunit, '(a,i0,a,i0,a)') '  <Piece NumberOfPoints="', N, '" NumberOfCells="', N, '">'

write(iunit, '(A)') '      <PointData>'
write(iunit, '(A)') '        <DataArray type="Int32" Name="Type" format="ascii">'
write(iunit, '(10I2.1)') IFLAG(1: N)  ! 假设ptype已初始化为整数
write(iunit, '(A)') '        </DataArray>'
write(iunit, '(A)') '      </PointData>'

! 写入粒子位置（2D转3D，z=0）
write(iunit, '(A)') '      <Points>'
write(iunit, '(A)') '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
do i = 1, N
    write(iunit, '(3F10.6)') X1(i), Y1(i), 0.0
end do
write(iunit, '(A)') '        </DataArray>'
write(iunit, '(A)') '      </Points>'

! 速度（矢量）
write(iunit, '(A)') '      <PointData Scalars="Pressure" Vectors="Velocity">'
write(iunit, '(A)') '        <DataArray type="Float32" Name="Velocity" NumberOfComponents="3" format="ascii">'
do i = 1, N
    write(iunit, '(3F10.6)') U1(i), V1(i), 0.0
end do
write(iunit, '(A)') '        </DataArray>'
  
! 压力（标量）
write(iunit, '(A)') '        <DataArray type="Float32" Name="Pressure" format="ascii">'
write(iunit, '(10F10.2)')  P(1:N)
write(iunit, '(A)') '        </DataArray>'
write(iunit, '(A)') '      </PointData>'

! 写入粒子单元（每个粒子作为一个顶点）
write(iunit, '(A)') '      <Cells>'
write(iunit, '(A)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
write(iunit, '(10I8)') (i-1, i=1,N)  ! VTK索引从0开始
write(iunit, '(A)') '        </DataArray>'
write(iunit, '(A)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
write(iunit, '(10I8)') (i, i=1,N)    ! 每个单元的偏移量
write(iunit, '(A)') '        </DataArray>'
write(iunit, '(A)') '        <DataArray type="UInt8" Name="types" format="ascii">'
write(iunit, '(10I4)') (1, i=1,N)    ! 1=VTK_VERTEX
write(iunit, '(A)') '        </DataArray>'
write(iunit, '(A)') '      </Cells>'

! 结束文件
write(iunit, '(A)') '    </Piece>'
write(iunit, '(A)') '  </UnstructuredGrid>'
write(iunit, '(A)') '</VTKFile>'

close(iunit)
END SUBROUTINE output