CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PROGRAMA  P R O J B A N D A
C
C    Agrupa datos x,z a intervalos dx y calcula el valor medio y la 
C    desviacion estandar en ese intervalo. Puede promediar con una funcion peso. 
C
C    INPUT COLUMNS: X,Z[,WEIGHT]   (X must increase monotonously)
C    OUPUT COLUMNS: XSAM,ZMEAN,SD
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DIMENSION X(100000),Y(100000),Z(100000),P(100000),ZP(100000)
      CHARACTER*40 FILEIN,FILEOUT
      A=-0.00036
      B=1.
      PMIN=0.1
      PRINT *,'ENTER INPUT FILE: '
      READ(5,'(A)') FILEIN
      PRINT *,'ENTER OUTPUT FILE: '
      READ(5,'(A)') FILEOUT
      PRINT *
      PRINT *,'ENTER SAMPLING INTERVAL: '
      READ(5,*) DX
      PRINT *
      PRINT *,'ENTER WINDOW WIDTH: '
      READ(5,*) DXW
      OPEN(1,FILE=FILEIN,STATUS='OLD')
      OPEN(2,FILE=FILEOUT,STATUS='UNKNOWN')
      I=1
    1 READ(1,*,END=9) X(I),Z(I),Y(I)
C     write(6,*) C1,C2,X(I),Z(I),Y(I)
      I=I+1
      GOTO 1
    9 CONTINUE
      I=I-1
      XSAM=X(1)
      GOTO 4
    3 XSAM=XSAM+DX
      IF(XSAM .GT. X(I))   GOTO 99
    4 N=0
      SP=0.
      SZ=0.
      SZZ=0.
      SZP=0.
      XMIN=XSAM-DXW/2.
      XMAX=XSAM+DXW/2.
      DO 2 J=1,I
        IF(X(J) .GE. XMIN .AND. X(J) .LE. XMAX)   THEN
            N=N+1
            XD=ABS(XSAM-X(J))
            YD=ABS(Y(J))
            D=SQRT(XD*XD+YD*YD)
            P(J)=A*D*D+B
            IF(P(J) .LT. PMIN)  P(J)=PMIN
            ZP(J)=Z(J)*P(J)
            SZP=SZP+ZP(J)
            SZ=SZ+Z(J)
            SZZ=SZZ+Z(J)*Z(J)
            SP=SP+P(J)
        ENDIF
    2 CONTINUE
      ZMEAN=SZP/SP
      SD=SQRT((SZZ-SZ*SZ/N)/N)
      WRITE(2,*) XSAM,ZMEAN,SD
      GOTO 3
   99 CONTINUE
      CLOSE(1)
      CLOSE(2)
      STOP
      END
