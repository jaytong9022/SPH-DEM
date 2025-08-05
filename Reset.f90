!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     THIS SUBROUTIME DOES NOT NEED TO BE CHANGED    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE reset
USE mdata
!! RESET PARTICLE VELOCITY AND POSITION TO PREVIOUS STEP
!! RESET VELOCITIES
U1=U2
V1=V2
!! RESET COORDINATES
X1=X2
Y1=Y2
RETURN
END SUBROUTINE reset