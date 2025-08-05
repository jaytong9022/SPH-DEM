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

! �������Ŀ¼���ļ���
dirname = 'output'

! ��鲢���� output Ŀ¼����������ڣ�
inquire(file=trim(dirname)//'/.', exist=dir_exists)  ! ���Ŀ¼�Ƿ����
if (.not. dir_exists) then
    call system('mkdir "' // trim(dirname) // '" > nul 2>&1')  ! �����������
end if

! ��ʱ�䲽��ʽ��Ϊ4λ���֣��� 5 �� "0005"��
write(time_str, '(I6.6)') time_step

! �����ļ�������ʽ����·���ָ�����
filename = trim(dirname) // '\sph_output_' // trim(time_str) // '.vtu'


! ���ļ�����������
open(newunit=iunit, file=filename, status='replace', action='write', iostat=stat)
if (stat /= 0) then
    print *, '����: �޷������ļ� '//trim(filename)
    stop
end if

! д��VTKͷ��
write(iunit, '(a)') '<?xml version="1.0" encoding="UTF-8"?>'  ! <<<
write(iunit, '(a)') '<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian">'
write(iunit, '(a)') '<UnstructuredGrid>'
write(iunit, '(a,i0,a,i0,a)') '  <Piece NumberOfPoints="', N, '" NumberOfCells="', N, '">'

write(iunit, '(A)') '      <PointData>'
write(iunit, '(A)') '        <DataArray type="Int32" Name="Type" format="ascii">'
write(iunit, '(10I2.1)') IFLAG(1: N)  ! ����ptype�ѳ�ʼ��Ϊ����
write(iunit, '(A)') '        </DataArray>'
write(iunit, '(A)') '      </PointData>'

! д������λ�ã�2Dת3D��z=0��
write(iunit, '(A)') '      <Points>'
write(iunit, '(A)') '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
do i = 1, N
    write(iunit, '(3F10.6)') X1(i), Y1(i), 0.0
end do
write(iunit, '(A)') '        </DataArray>'
write(iunit, '(A)') '      </Points>'

! �ٶȣ�ʸ����
write(iunit, '(A)') '      <PointData Scalars="Pressure" Vectors="Velocity">'
write(iunit, '(A)') '        <DataArray type="Float32" Name="Velocity" NumberOfComponents="3" format="ascii">'
do i = 1, N
    write(iunit, '(3F10.6)') U1(i), V1(i), 0.0
end do
write(iunit, '(A)') '        </DataArray>'
  
! ѹ����������
write(iunit, '(A)') '        <DataArray type="Float32" Name="Pressure" format="ascii">'
write(iunit, '(10F10.2)')  P(1:N)
write(iunit, '(A)') '        </DataArray>'
write(iunit, '(A)') '      </PointData>'

! д�����ӵ�Ԫ��ÿ��������Ϊһ�����㣩
write(iunit, '(A)') '      <Cells>'
write(iunit, '(A)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
write(iunit, '(10I8)') (i-1, i=1,N)  ! VTK������0��ʼ
write(iunit, '(A)') '        </DataArray>'
write(iunit, '(A)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
write(iunit, '(10I8)') (i, i=1,N)    ! ÿ����Ԫ��ƫ����
write(iunit, '(A)') '        </DataArray>'
write(iunit, '(A)') '        <DataArray type="UInt8" Name="types" format="ascii">'
write(iunit, '(10I4)') (1, i=1,N)    ! 1=VTK_VERTEX
write(iunit, '(A)') '        </DataArray>'
write(iunit, '(A)') '      </Cells>'

! �����ļ�
write(iunit, '(A)') '    </Piece>'
write(iunit, '(A)') '  </UnstructuredGrid>'
write(iunit, '(A)') '</VTKFile>'

close(iunit)
END SUBROUTINE output