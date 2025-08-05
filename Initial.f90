!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    THIS IS THE INITIALIZATION SUBROUTINE - NEED TO BE CHANGED ACCORDING TO DIFFERENT APPLICATIONS     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE init
USE mdata

!!Added By TongJ
Xrt = 0.8
Xlt = -0.1
Xa1 = 0.0
Xa2 = 0.7
Yup = 0.5
Ydn = 0.0
Xp = 0.2
Yp = 0.3
Sp = 0.01
!! INITIALIZATION OF WALL PARTICLES
!! LEFT WALL......
N=1
X1(1)=0.0
Y1(1)=Yup
U1(1)=0.0
V1(1)=0.0
IFLAG(1)=1
DO WHILE(ABS(Y1(N)-Sp)>=1.0E-4)
N=N+1
Y1(N)=Y1(N-1)-Sp
X1(N)=X1(N-1)
U1(N)=0.0
V1(N)=0.0
IFLAG(N)=1
END DO
!!  HORIZONTAL WALL......
N=N+1
X1(N)=Xlt
Y1(N)=0.0
U1(N)=0.0
V1(N)=0.0
IFLAG(N)=1
DO WHILE(ABS(X1(N)-(Xrt+Sp))>=1.0E-4)
N=N+1
Y1(N)=Y1(N-1)
X1(N)=X1(N-1)+Sp
U1(N)=0.0
V1(N)=0.0
    IF((X1(N)>0.0).AND.(0.7-X1(N)>1.0E-4)) THEN
    IFLAG(N)=3
    ELSE
    IFLAG(N)=1
    ENDIF
END DO
!!	RIGHT WALL
X1(N)=Xa2
Y1(N)=Sp !!Careful
U1(N)=0.0
V1(N)=0.0
IFLAG(N)=1
DO WHILE(ABS(Y1(N)-(Yup))>=1.0E-4)
N=N+1
Y1(N)=Y1(N-1)+Sp
X1(N)=X1(N-1)
U1(N)=0.0
V1(N)=0.0
IFLAG(N)=1
END DO
!! WRITING INITIAL FILE
!! OPEN(1,FILE='O.DAT')
!! WRITE(1,21) (I,X1(I),Y1(I),IFLAG(I),I=1,N)
!! 21 FORMAT(I5,2F10.3,I5)
!! CLOSE(1)

!! INITIALIZATION OF INNER PARTICLES
!! FROM LEFT BOUNDARY TO RIGHT BOUNDARY

DO I=INT((Yup-Ydn)/Sp+(Xa1-Xlt)/Sp+2),INT((Yup-Ydn)/Sp+(Xa1-Xlt+Xp)/Sp+1)
	N=N+1
	X1(N)=X1(I)
	Y1(N)=Y1(I)+Sp
    U1(N)=0.0
	V1(N)=0.0
	IFLAG(N)=2
    DO WHILE(Y1(N)<(Yp-Sp))
    N=N+1
	X1(N)=X1(I)
	Y1(N)=Y1(N-1)+Sp
	U1(N)=0.0
	V1(N)=0.0
	IFLAG(N)=2
    END DO
END DO


!!OPEN(1,FILE='O.DAT')
!!WRITE(1,21) (I,X1(I),Y1(I),IFLAG(I),I=1,N)
!!21  FORMAT(I5,2F10.4,I5)
!!CLOSE(1)

!!  INITIALIZATION OF DUMMY PARTICLE NUMBERS AND COORDINATES
!! FIRST LINE OF DUMMY PARTICELS
!! LEFT WALL......
    DO I=1,INT((Yup-Ydn)/Sp)
	N=N+1
	X1(N)=X1(I)-Sp
	Y1(N)=Y1(I)
    IFLAG(N)=0
	END DO
!!  HORIZONTAL WALL......
    DO I=INT((Yup-Ydn)/Sp)+1,INT((Yup-Ydn)/Sp+(Xrt-Xlt)/Sp)+1
	N=N+1
	X1(N)=X1(I)
	Y1(N)=Y1(I)-Sp
    IFLAG(N)=0
	END DO
!!  RIGHT WALL......
    DO I=INT((Yup-Ydn)/Sp+(Xrt-Xlt)/Sp)+2,INT((Yup-Ydn)*2/Sp+(Xrt-Xlt)/Sp)+1
	N=N+1
	X1(N)=X1(I)+Sp
	Y1(N)=Y1(I)
    IFLAG(N)=0
	END DO
!!  SECOND LINE OF DUMMY PARTICELS
!!  LEFT WALL......
    DO I=1,INT((Yup-Ydn)/Sp)
	N=N+1
	X1(N)=X1(I)-2*Sp
	Y1(N)=Y1(I)
    IFLAG(N)=0
	END DO
!!  HORIZONTAL WALL
    DO I=INT((Yup-Ydn)/Sp)+1,INT((Yup-Ydn)/Sp+(Xrt-Xlt)/Sp)+1
	N=N+1
	X1(N)=X1(I)
	Y1(N)=Y1(I)-2*Sp
    IFLAG(N)=0
	END DO
!!  RIGHT WALL......
    DO I=INT((Yup-Ydn)/Sp+(Xrt-Xlt)/Sp)+2,INT((Yup-Ydn)*2/Sp+(Xrt-Xlt)/Sp)+1
	N=N+1
	X1(N)=X1(I)+2*Sp
	Y1(N)=Y1(I)
    IFLAG(N)=0
	END DO

    !!  THIRD LINE OF DUMMY PARTICELS
!!  LEFT WALL......
    DO I=1,INT((Yup-Ydn)/Sp)
	N=N+1
	X1(N)=X1(I)-2*Sp
	Y1(N)=Y1(I)
    IFLAG(N)=0
	END DO
!!  HORIZONTAL WALL
    DO I=INT((Yup-Ydn)/Sp)+1,INT((Yup-Ydn)/Sp+(Xrt-Xlt)/Sp)+1
	N=N+1
	X1(N)=X1(I)
	Y1(N)=Y1(I)-2*Sp
    IFLAG(N)=0
	END DO
!!  RIGHT WALL......
    DO I=INT((Yup-Ydn)/Sp+(Xrt-Xlt)/Sp)+2,INT((Yup-Ydn)*2/Sp+(Xrt-Xlt)/Sp)+1
	N=N+1
	X1(N)=X1(I)+2*Sp
	Y1(N)=Y1(I)
    IFLAG(N)=0
	END DO

    !!  FORTH LINE OF DUMMY PARTICELS
!!  LEFT WALL......
    DO I=1,INT((Yup-Ydn)/Sp)
	N=N+1
	X1(N)=X1(I)-2*Sp
	Y1(N)=Y1(I)
    IFLAG(N)=0
	END DO
!!  HORIZONTAL WALL
    DO I=INT((Yup-Ydn)/Sp)+1,INT((Yup-Ydn)/Sp+(Xrt-Xlt)/Sp)+1
	N=N+1
	X1(N)=X1(I)
	Y1(N)=Y1(I)-2*Sp
    IFLAG(N)=0
	END DO
!!  RIGHT WALL......
    DO I=INT((Yup-Ydn)/Sp+(Xrt-Xlt)/Sp)+2,INT((Yup-Ydn)*2/Sp+(Xrt-Xlt)/Sp)+1
	N=N+1
	X1(N)=X1(I)+2*Sp
	Y1(N)=Y1(I)
    IFLAG(N)=0
	END DO

 !!OPEN(2,FILE='O2.DAT')
 !!WRITE(2,22) (I,XDUMMY(I),YDUMMY(I),I=1,N)
 !!22 FORMAT(I5,2F10.3)
 !!CLOSE(2)

!!  MASS OF INNER PARTICLE
	DO I=1,N
	PM(I)=1000*Sp*Sp
	END DO
!!  MASS OF DUMMY PARTICLE
	PMD=1000*Sp*Sp
!!  DENSITY OF INNER PARTICLE
	DO I=1,N
	ROU(I)=1000.0
	END DO
!!  DENSITY OF DUMMY PARTICLE
	ROUDUMMY=1000.0
!!  FOR LAMINAR FLOW EDDY VISCOSITY IS ZERO
	DO I=1,N
	VISEFF(I)=0.0
	END DO
!!  END OF INITIALIZATION
RETURN
END SUBROUTINE init