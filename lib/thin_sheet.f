CC **********************************************************************
CC  Set of subroutines to calculate:
CC   -The velocity field of a thin sheet (subroutine velocity_field)
CC      from the boundary conditions (subroutine READ_BC)
CC   -Thickening of the sheet (subroutine THICKEN)
CC **********************************************************************

      SUBROUTINE velocity_field (TEMPSMa, Dt, AX, BY, n, m, nn,
     +			vissup, visinf, visTer, vis, nincogn,
     +			tallmax, alfa, nitermax, average_Pressure, u, v,
     +			BCfile, nbanda)

CC  TEMPSMa[Ma]	If != 0 then uses the previous velocity field.
CC  Dt [s]   	Time increment, if the amount of deformation is too 
CC			large then Dt will be reduced to fit 'quotadef'
CC  AX [m]	domain length in x
CC  BY [m]	domain length in y
CC  n, m	Number of nodes (in x,y directions) minus 1: column index 
CC			runs from 0 (west) to n (east)
CC  nn		Number of nodes (n+1)*(m+1)
CC  nincogn	2*nn
CC  nbanda	4(n+1)+7 
CC  nitermax	max. number of iterations for convergence between 
CC			viscosity and velocity (if exceeded => divergence).
CC			If 0 => viscosity doesn't depend on strain rate.
CC  tallmax	Criterium of convergence. If SUM(Dvel/vel) < tallmax =>
CC			convergence. (p.e., =0.01).
CC  alfa	0 => explicit ;  alfa=1 => implicit. 0.5 recommended
CC  visTer [Pa]	Array of thermal term of viscosity: 
CC			1/2 * stregth[Pa*m] / layer.thickness[m]
CC  vis [Pa*s]	Array of effective viscosity of the layer:
CC			visTer / ref.strain.rate
CC			Returns the calculated viscosity (a function of 
CC			the new strain rate).
CC  vissup [Pa*s]  Imposed upper limit for calculated viscosity 'vis'. (~10^25 Pa*s)
CC  visinf [Pa*s]  Imposed lower limit for calculated viscosity 'vis'. (~10^22 Pa*s)
CC  average_Pressure [Pa=N/m2] Integral of (P(z)*dz)/thickness_layer between the surface (z=topo) and Zcomp.
CC                      Where P(z) is the lithostatic pressure (integral of rho(z)*g*dz), 
CC			Zcomp is the depth of compensation, and thickness_layer=Zcomp+elevation. 
CC			The independent term of the thin_sheet equilibrium equation is the lateral gradient of -average_Pressure.
CC			Example: in a thickenned crust and in a ridge (thinner lithosphere) this average_Pressure is lower.
!!			The velocity field goes from the lowest to the highest average_Pressure.
CC  u, v [m/s]	Arrays of velocity in x, y directions in the previous time step.
CC			Sorted from west to east and south to north.
CC			Returns the new calculated velocity.
CC  BCfile [char*84]  Boundary conditions file name 

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*84 TITOL_BC, BCfile
      PARAMETER (pi=3.1415926535897932D0, FACTEMP=3.1536D7,  
     +		quotadef=8.D-2, DP_sup=10000.D6, DP_inf=-10000.D6)
      DIMENSION a(nincogn,nbanda), u(nn), v(nn), b(nincogn),
     +		average_Pressure(nn), visTer(nn), vis(nn), visold(nn),
     +		vel(nincogn), velold(nincogn)

CC  Parametres per a les Condicions de Contorn
      PARAMETER (IBC=2,IFAULT=1)   
C	vr=10/3.1536D10
        
      tall=0.D0
      Dx=AX/n
      Dy=BY/m
      kvel=1
      DO kixy=1,nn
          visold(kixy)=vis(kixy)
          vel(kvel)=0.D0
          vel(kvel+1)=0.D0
          velold(kvel)=u(kixy)
          velold(kvel+1)=v(kixy)
          kvel=kvel+2
      END DO
CC *************** CONSTRUCCIO DE LA MATRIU *********************
      Li=2*(n+1)+3
      Ls=Li
      Ld=Li+1
      kp=Ld+2*(n+1) 
      kn=Ld-2*(n+1)
       
      niter=0
 111   CONTINUE 
       niter=niter+1
       a=0.D0
       b=0.D0

CC  ----------  Condicions de Contorn --------------------
CC  POL DE ROTACIO: 1:NUVEL-1A.  2:Argus. 3:RM2.  4:P071.  5:Mckenzie
		IF(IBC.EQ.1) THEN
			dlonpol=-20.6D0
			dlatpol=21.0D0
			omegada=0.12D0
		ENDIF
		IF(IBC.EQ.2) THEN
			dlonpol=-20.3D0
			dlatpol=18.8D0
			omegada=0.104D0
		ENDIF
		IF(IBC.EQ.3) THEN
			dlonpol=-21.19D0
			dlatpol=25.23D0
			omegada=0.1D0
		ENDIF
		IF(IBC.EQ.4) THEN
			dlonpol=-23.5D0
			dlatpol=29.2D0
			omegada=0.14D0
		ENDIF
		IF(IBC.EQ.5) THEN
			dlonpol=-28.2D0
			dlatpol=22.7D0
			omegada=0.27D0
		ENDIF
C       CALL fixCCSud (a,b,m,n,nincogn,nbanda,Dx,Dy,vr,
C     +			dlonpol,dlatpol,omegada)
C       CALL fixCCWest (a,b,m,n,nincogn,nbanda,Dx,Dy,vr,
C     +			dlonpol,dlatpol,omegada,vis,nn)
C       CALL fixCCEst (a,b,m,n,nincogn,nbanda,Dx,Dy,vr,
C     +			dlonpol,dlatpol,omegada,IFAULT)
C       CALL fixCCNort (a,b,m,n,nincogn,nbanda,Dx,Dy,vr)
       
      CALL READ_BC (TEMPSMa,niter,n,m,nn,nincogn,nbanda,Dx,Dy,vis,
     +			a,b,TITOL_BC, BCfile)

C  ------------------------------------------------------
C   PUNTS INTERIORS
       DP_max=-1.D100
       DP_min=1.D100
       Dvis_max=-1.D100
       Dvis_min=1.D100
      Point_iy: DO iy=1,m-1
         Point_ix: DO ix=1,n-1
               k=ix+1+iy*(n+1)
                GLv=vis(k)
                G1=vis(k+1)
                G2=vis(k-1)
                G3=vis(k+n+1)
                G4=vis(k-n-1)
C  ---------------------------------------------------------------------
               GLvx=(G1-G2)/(2.D0*Dx)
               GLvy=(G3-G4)/(2.D0*Dy)
	       Dvis_max=MAX(Dvis_max,GLvx,GLvy)
	       Dvis_min=MIN(Dvis_min,GLvx,GLvy)
	       DP_x=(average_Pressure(k+1)-average_Pressure(k-1))/
     +	       		(2.D0*Dx)
	       DP_y=(average_Pressure(k+n+1)-average_Pressure(k-n-1))/
     +	       		(2.D0*Dy)
CC	LIMIT DEL GRADIENT DE LA average_Pressure
C		DP_x=MIN(DP_x,DP_sup)
C	       	DP_x=MAX(DP_x,DP_inf)
C		DP_y=MIN(DP_y,DP_sup)
C	       	DP_y=MAX(DP_y,DP_inf)
               DP_max=MAX(DP_max,DP_x,DP_y)
               DP_min=MIN(DP_min,DP_x,DP_y)
	       Tx=-DP_x
	       Ty=-DP_y
               l=2*ix+1+2*iy*(n+1)
               DIAGONAL1=-(8.D0*GLv)/(Dx*Dx)-(2.D0*GLv)/(Dy*Dy)
               a(l,Ld+2)=((4.D0*GLv)/(Dx*Dx)+(2.D0*GLvx)/Dx)/DIAGONAL1
               a(l,kp)=(GLv/(Dy*Dy)+GLvy/(2.D0*Dy))/DIAGONAL1
               a(l,Ld)=1.D0
               a(l,kn)=(GLv/(Dy*Dy)-GLvy/(2.D0*Dy))/DIAGONAL1
               a(l,Ld-2)=((4.D0*GLv)/(Dx*Dx)-(2.D0*GLvx)/Dx)/DIAGONAL1
               a(l,kp+3)=((3.D0*GLv)/(4.D0*Dx*Dy))/DIAGONAL1
               a(l,kn+3)=-a(l,kp+3)
               a(l,kp-1)=-a(l,kp+3)
               a(l,kn-1)=a(l,kp+3)
               a(l,kp+1)=(GLvx/Dy)/DIAGONAL1
               a(l,kn+1)=-a(l,kp+1)
               a(l,Ld+3)=(GLvy/(2.D0*Dx))/DIAGONAL1
               a(l,Ld-1)=-a(l,Ld+3)
               b(l)=Tx/DIAGONAL1
               l=2*ix+2+2*iy*(n+1)
               DIAGONAL2=-2.D0*GLv*(1.D0/(Dx*Dx)+4.D0/(Dy*Dy))
               a(l,Ld+2)=(GLv/(Dx*Dx)+GLvx/(2.D0*Dx))/DIAGONAL2
               a(l,kp)=((4.D0*GLv)/(Dy*Dy)+(2.D0*GLvy)/Dy)/DIAGONAL2
               a(l,Ld)=1.D0
               a(l,kn)=((4.D0*GLv)/(Dy*Dy)-(2.D0*GLvy)/Dy)/DIAGONAL2
               a(l,Ld-2)=(GLv/(Dx*Dx)-GLvx/(2.D0*Dx))/DIAGONAL2
               a(l,kp+1)=((3.D0*GLv)/(4.D0*Dx*Dy))/DIAGONAL2
               a(l,kn+1)=-a(l,kp+1)
               a(l,kp-3)=-a(l,kp+1)
               a(l,kn-3)=a(l,kp+1)
               a(l,Ld+1)=(GLvy/Dx)/DIAGONAL2
               a(l,Ld-3)=-a(l,Ld+1)
               a(l,kp-1)=(GLvx/(2.D0*Dy))/DIAGONAL2
               a(l,kn-1)=-a(l,kp-1)
               b(l)=Ty/DIAGONAL2
         END DO Point_ix 
      END DO Point_iy 

      CALL sistbanda(a,b,nincogn,nbanda,Li,vel)
      
      IF(TEMPSMa.EQ.0.D0.AND.niter.EQ.1) GOTO 19 
      IF(nitermax.LE.1) GOTO 19
      vel=alfa*vel+(1.D0-alfa)*velold
C         DO 20 kvel=1,nincogn 
C  20             vel(kvel)=alfa*vel(kvel)+(1.D0-alfa)*velold(kvel)

 19    CONTINUE
       epmig=0.D0
       cepmig=0.D0
       epmax=0.D0
       NVISSUP=0
       NVISINF=0
       IF(nitermax.EQ.0) THEN
!	  PRINT*,'DELTA PRESSIO:minim =',DP_min/1.D6,'MPa/m(10**6 N/m3)'
!	  PRINT*,'              maxim =',DP_max/1.D6,'MPa/m(10**6 N/m3)' 
!	    PRINT*,' Dvis_min =',Dvis_min, ', Dvis_max =',Dvis_max 
            PRINT*,'   No Iteration. The viscosity is the same'
              DO iy=1,m-1
                  DO ix=1,n-1
                      kxy=ix+1+iy*(n+1)
                      kd=2*ix+1+2*iy*(n+1)
                      kdp=kd+2*(n+1)
                      kdn=kd-2*(n+1)
                      ux=(vel(kd+2)-vel(kd-2))/(2.D0*Dx)
                      vx=(vel(kd+3)-vel(kd-1))/(2.D0*Dx)
                      vy=(vel(kdp+1)-vel(kdn+1))/(2.D0*Dy)
                      uy=(vel(kdp)-vel(kdn))/(2.D0*Dy)
                      epuntzz=-(ux+vy)
	    	      Esecond=((ux*ux+vy*vy+epuntzz*epuntzz)/2.D0)+
     +                            ((1.D0/4.D0)*(uy+vx)*(uy+vx))	
     	    	      epeffec=DSQRT(Esecond)
		      epmax=MAX(ABS(epeffec),epmax)
	    	      epmig=epmig+ABS(epeffec)
		      cepmig=cepmig+1
                  END DO  
              END DO  
	    epmig=epmig/cepmig
            vismax=MAXVAL(vis)
            vismin=MINVAL(vis)
            GOTO 999
       ENDIF
      IF(niter.EQ.1) WRITE(6,60) visinf,vissup
 60   FORMAT(/4X,'Iteration',31X,'__ nodes limited for viscosity __'/
     +  4X'number',5X,'mean Dv',6X,'mean Dvis',7X,
     +  1P,G9.2,' Pa.s',4X,1P,G9.2,' Pa.s') 

C ****************  TROBO LA NOVA VISCOSITAT *********************
      Vis_iy: DO iy=1,m-1
         Vis_ix: DO ix=1,n-1
	    kxy=ix+1+iy*(n+1)
            kd=2*ix+1+2*iy*(n+1)
            kdp=kd+2*(n+1)
            kdn=kd-2*(n+1)
            ux=(vel(kd+2)-vel(kd-2))/(2.D0*Dx)
            vx=(vel(kd+3)-vel(kd-1))/(2.D0*Dx)
            vy=(vel(kdp+1)-vel(kdn+1))/(2.D0*Dy)
            uy=(vel(kdp)-vel(kdn))/(2.D0*Dy)
            epuntzz=-(ux+vy)
C            Epunt=DSQRT(2.D0*(ux*ux+vy*vy+ux*vy+
C     +                            (1.D0/4.D0)*(uy+vx)*(uy+vx)))

!!!  Controlar quin strain rate utilitza per calcular la viscositat !!!
	    Esecond=((ux*ux+vy*vy+((ux+vy)*(ux+vy)))/2.D0)+
     +                            ((1.D0/4.D0)*(uy+vx)*(uy+vx))	
     	    epeffec=DSQRT(Esecond)
	    epmig=epmig+ABS(epeffec)
	    cepmig=cepmig+1
	    epmax=MAX(ABS(epeffec),epmax)
	    visnova=ABS(visTer(kxy)/epeffec)
            vis(kxy)=alfa*visnova+(1.D0-alfa)*visold(kxy)
            IF(vis(kxy).GT.vissup) THEN
            	vis(kxy)=vissup
           	NVISSUP=NVISSUP+1
            ENDIF
            IF(vis(kxy).LT.visinf) THEN
            	vis(kxy)=visinf
           	NVISINF=NVISINF+1
            ENDIF
         END DO Vis_ix
      END DO Vis_iy
      
      epmig=epmig/cepmig
      Dumig=0.D0
!      Dumax=0.D0
      velmax=0.D0
      Dvismig=0.D0
      Dvismax=0.D0
      NPU=0
      DO iy=1,m-1
         DO ix=1,n-1
              kd=ix+1+iy*(n+1)
              k=2*ix+1+2*iy*(n+1)
              velmod=DSQRT(vel(k)*vel(k)+vel(k+1)*vel(k+1))
              veloldm=DSQRT(velold(k)*velold(k)+velold(k+1)*velold(k+1))
              Du=(ABS(velmod-veloldm))/velmod  
              Dumig=Dumig+Du  
!              Dumax=MAX(Du,Dumax)
              Dvis=(ABS(visold(kd)-vis(kd)))/vis(kd)
              Dvismig=Dvismig+Dvis
              Dvismax=MAX(Dvis,Dvismax) 
              velmax=MAX(velmod,velmax) 
	      NPU=NPU+1
         END DO
      END DO
      Dvismig=Dvismig/NPU 
      Dumig=Dumig/NPU
      vismax=MAXVAL(vis)
      vismin=MINVAL(vis)

      IF(velmax.EQ.0.D0.OR.vismax.EQ.0.D0) THEN
              PRINT*,'             EL CAMP DE VELOCITATS ES NUL'
              PRINT*,' Terme maxim de la velocitat (velmax) =',velmax,
     +           '     Terme maxim de la viscositat (vismax) =',vismax
              GOTO 999 
      ENDIF
       tall=Dumig
       tallvis=Dvismig
       WRITE(6,65) niter,tall,tallvis,NVISINF,NVISSUP 
 65    FORMAT(4X,I3,5X,F10.5,4X,F10.5,12X,I4,14X,I4)  
       if(niter.ge.nitermax) THEN    
             PRINT*,'   Maximum number of iterations!'
             GOTO 999
       ENDIF
       IF(tall.GT.tallmax) then
           k=1
	   visold=vis
           velold=vel
           GOTO 111
       ENDIF   
       WRITE(6,"(4X,'velocity converged')")     
999    CONTINUE 

CC    FIXO LA VISCOSITAT A LES VORES 
           DO ix=1,n-1
                 kxys=ix+1
                 kxyn=ix+1+m*(n+1)
                 vis(kxys)=vis(kxys+n+1)
                 vis(kxyn)=vis(kxyn-n-1)
           END DO
           DO iy=0,m
                 kxyw=1+iy*(n+1)
                 kxye=n+1+iy*(n+1)
                 vis(kxyw)=vis(kxyw+1)
                 vis(kxye)=vis(kxye-1)
           END DO

C ----------------------------------------------------------------------
      strainmx=epmax*Dt
      WRITE(6,61) vismin,vismax,epmig,epmax,epmig*Dt,strainmx
 61   FORMAT(
     +  4X,'Viscosity:            minim:',1P,G12.4,' Pa.s,   maxim:',
     +    1P,G12.4,' Pa.s'/
     +  4X,'Effect. strain rate:  medium:',1P,G12.4,' s-1,   maximum:',
     +    1P,G12.4,' s-1'/
     +  4X,'Effect. strain:       medium:',1P,G12.4,'        maximum:',
     +    1P,G12.4)
     	 
          IF(strainmx.GT.quotadef) THEN
               Dt=(quotadef/epmax)-5.D4
               Dtany=Dt/FACTEMP
               PRINT*,'TOO MUCH DEFORMATION -> DECREASE '
               PRINT*,'THE TIME INTERVAL,  Dt =',Dtany,' anys'
               strainmx=epmax*Dt
               PRINT*,'         NEW MAXIMUM STRAIN :',strainmx
          ENDIF

CC	Divideixo el vector vel (velocitat) en dos vectors:
CC	u (velocitat en x) i v (velocitat en y).
	kvel=1
      DO iy=0,m
	 DO ix=0,n
	     kxy=ix+1+iy*(n+1)
	     u(kxy)=vel(kvel)
	     v(kxy)=vel(kvel+1)
	     kvel=kvel+2
         END DO
      END DO

      
      RETURN
      END SUBROUTINE velocity_field
       
CC ********************************************************************
CC ********************************************************************

      SUBROUTINE READ_BC (TEMPSMa,niter,n,m,nn,nincogn,nbanda,Dx,Dy,
     +				vis,a,b,TITOL_BC,BCfile)   

CC   Read the Boundary Conditions from the file: BC.in
CC   First line: Title
CC	ix,iy,ITBC,t1,t2 
CC	ITBC: Boundary Condition type.
CC	t1,t2: Boundary Condition 1 and 2.

CC  ITBC=1    velocity fixed 				->  vel_x=t1  and  vel_y=t2 
CC  ITBC=12   stress xx and xy fixed (Est, West free)	-> tau_xx=t1  and  tau_xy=t2 
CC  ITBC=13   stress xx and yy fixed   		   	-> tau_xx=t1  and  tau_yy=t2 
CC  ITBC=23   stress xy and yy fixed (North, South free)-> tau_xy=t1  and  tau_yy=t2 
CC  ITBC=4    free slip (vel normal=0, tau_xy=0)	-> don't use t1 and t2

CC  Stress free -> normal stress and xy zero.

    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION a(nincogn,nbanda),b(nincogn),vis(nn)
      CHARACTER*84 TITOL_BC,BCfile
      PARAMETER (NBCmax=600,FACVEL=3.1536D10)

      OPEN(4,file=BCfile)
      NBC=2*(n+1)+2*(m-1)
      Li=2*(n+1)+3
      Ls=Li
      Ld=Li+1
      kp=Ld+2*(n+1)
      kn=Ld-2*(n+1)

       READ(4,1001) TITOL_BC
 1001  FORMAT(A80)
!       IF(TEMPSMa.EQ.0.AND.niter.EQ.1) WRITE (6,3) TITOL_BC
 3     FORMAT (/2X,'READ BOUNDARY CONDITIONS: ',/4X,A80/
     +    4X,'(ix, iy, BC type:   u,v / tau_1,tau_2 )')
 33     FORMAT (1X,3I4,1P,2G14.6)
 
       NBC_read=0
       DO 7 inBC=1,NBCmax
           READ(4,*,IOSTAT=IOS) ix,iy,ITBC,t1,t2
           IF (IOS < 0) THEN
               GO TO 120
           ELSE IF (IOS > 0) THEN
               WRITE (*, 90) NDATA
 90            FORMAT (/' ERROR: Encountered bad data after ',
     +                  I6,' successful READs.')
               STOP
           ENDIF
!       	IF(TEMPSMa.EQ.0.AND.niter.EQ.1) WRITE (6,33) ix,iy,ITBC,t1,t2
           NBC_read=NBC_read+1
	   leq1=2*ix+1+2*iy*(n+1)
	   leq2=2*ix+2+2*iy*(n+1)
	   kxy=ix+1+iy*(n+1)
	   kp=Ld+2*(n+1)
           kn=Ld-2*(n+1)
CCC ====================================================
CCCC	ITBC=1  fixo (u,v)
	    IF(ITBC.EQ.1) THEN
C	    	PRINT*,'A: ix,iy:',ix,iy,'  ITBC:',ITBC
	   	vel_x=t1
		vel_y=t2
             	a(leq1,Ld)=1.D0
             	b(leq1)=vel_x
             	a(leq2,Ld)=1.D0
             	b(leq2)=vel_y
		GOTO 99
	    ENDIF	   
CCC ====================================================
CCCC	ITBC=12  fixo:  tau_xx (eq.1) i  tau_xy (eq.2)
CCCC	ITBC=13  fixo:  tau_xx (eq.1) i  tau_yy (eq.2)
CCCC	ITBC=23  fixo:  tau_xy (eq.1) i  tau_yy (eq.2)
	    IF(ITBC.GT.11.AND.ITBC.LT.20) THEN
	   	tau_xx=t1
		GOTO 411
	    ENDIF
 412	    CONTINUE
	    IF(ITBC.EQ.12) THEN
		tau_xy=t2
		ieq=1
		GOTO 522
	    ENDIF
	    IF(ITBC.EQ.13) THEN
		tau_yy=t2
		GOTO 813
	    ENDIF
	    IF(ITBC.EQ.23) THEN
		tau_xy=t1
		ieq=0
		GOTO 522
	    ENDIF
 525	    CONTINUE
	    IF(ITBC.EQ.23) THEN
		tau_yy=t2
		GOTO 813
	    ENDIF
CCC ====================================================
CCCC	ITBC=4  free slip (v norm=0, dv(tang)/dx(tang) =0)
	    IF(ITBC.EQ.4) THEN
		IF(iy.EQ.0) THEN
             	    a(leq1,kp)=1.D0
             	    a(leq1,Ld)=-1.D0
             	    b(leq1)=0.D0
             	    a(leq2,Ld)=1.D0
		    b(leq2)=0.D0
		    GOTO 99
	    	ENDIF		
		IF(iy.EQ.m) THEN
             	    a(leq1,Ld)=1.D0
             	    a(leq1,kn)=-1.D0
             	    b(leq1)=0.D0
             	    a(leq2,Ld)=1.D0
		    b(leq2)=0.D0
		    GOTO 99
	    	ENDIF		
		IF(ix.EQ.0) THEN
             	    a(leq1,Ld)=1.D0
             	    b(leq1)=0.D0
             	    a(leq2,Ld+2)=1.D0
             	    a(leq2,Ld)=-1.D0
		    b(leq2)=0.D0
		    GOTO 99
	    	ENDIF		
		IF(ix.EQ.n) THEN
             	    a(leq1,Ld)=1.D0
             	    b(leq1)=0.D0
             	    a(leq2,Ld)=1.D0
             	    a(leq2,Ld-2)=-1.D0
		    b(leq2)=0.D0
		    GOTO 99
	    	ENDIF		
	    ENDIF		

CCCC---------------------------------------------------------
CCCC---------------------------------------------------------
CC  Condition   tau_xx:  Equation 1
 411		CONTINUE
		IF(ix.NE.0.AND.ix.NE.n) THEN	
		    a(leq1,Ld+2)=1.D0
		    a(leq1,Ld-2)=-1.D0
             	    b(leq1)=(Dx*tau_xx)/vis(kxy)
		ENDIF
		IF(ix.EQ.0) THEN	
		    a(leq1,Ld+2)=1.D0
		    a(leq1,Ld)=-1.D0
             	    b(leq1)=(Dx*tau_xx)/(2.D0*vis(kxy))
		ENDIF
		IF(ix.EQ.n) THEN	
		    a(leq1,Ld)=1.D0
		    a(leq1,Ld-2)=-1.D0
             	    b(leq1)=(Dx*tau_xx)/(2.D0*vis(kxy))
		ENDIF
		GOTO 412
CCCC---------------------------------------------------------
CCCC---------------------------------------------------------
CC  Condition	tau_xy:  Equation 1 (ieq=0),	  Equation 2 (ieq=1)
 522		CONTINUE
		leq=leq1+ieq
		IF(iy.NE.0.AND.iy.NE.m) THEN
		    a(leq,kp-ieq)=1.D0/(2.D0*Dy)
		    a(leq,kn-ieq)=-1.D0/(2.D0*Dy)
		ENDIF
		IF(iy.EQ.0) THEN
		    a(leq,kp-ieq)=1.D0/Dy
		    a(leq,Ld-ieq)=-1.D0/Dy
		ENDIF
		IF(iy.EQ.m) THEN
		    a(leq,Ld-ieq)=1.D0/Dy
		    a(leq,kn-ieq)=-1.D0/Dy
		ENDIF
		IF(ix.NE.0.AND.ix.NE.n) THEN
		    a(leq,Ld+3-ieq)=1.D0/(2.D0*Dx)
		    a(leq,Ld-1-ieq)=-1.D0/(2.D0*Dx)
		ENDIF		    
		IF(ix.EQ.0) THEN
		    a(leq,Ld+3-ieq)=1.D0/Dx
		    a(leq,Ld+1-ieq)=-1.D0/Dx
		ENDIF
		IF(ix.EQ.n) THEN
		    a(leq,Ld+1-ieq)=1.D0/Dx
		    a(leq,Ld-1-ieq)=-1.D0/Dx
		ENDIF
             	b(leq)=tau_xy/vis(kxy)
	    	IF(ieq.EQ.1) GOTO 99
	    	IF(ITBC.EQ.23) GOTO 525
CCCC---------------------------------------------------------
CCCC---------------------------------------------------------
CC  Condition   tau_yy:  Equation 2
 813		CONTINUE
		IF(iy.NE.0.AND.iy.NE.m) THEN	
		    a(leq2,kp)=1.D0
		    a(leq2,kn)=-1.D0
             	    b(leq2)=(Dy*tau_yy)/vis(kxy)
		ENDIF
		IF(iy.EQ.0) THEN	
		    a(leq2,kp)=1.D0
		    a(leq2,Ld)=-1.D0
             	    b(leq2)=(Dy*tau_yy)/(2.D0*vis(kxy))
		ENDIF
		IF(iy.EQ.m) THEN	
		    a(leq2,Ld)=1.D0
		    a(leq2,kn)=-1.D0
             	    b(leq2)=(Dy*tau_yy)/(2.D0*vis(kxy))
		ENDIF
	        GOTO 99
CCCC---------------------------------------------------------
		
 99	CONTINUE
 7      CONTINUE

 100  CONTINUE
      IF(NBC_read.EQ.0) PRINT*,' FILE NO READ. number of data:',NBC_read
      WRITE(6,"(/' Dataset contained more than ',I5,
     +           ' data; Check the dimensions.')") NBCmax
 120  CONTINUE
!      WRITE(6,"(/5X,'Reading of data completed:',I7,' data points'/ )")
!     +            NBC_read
 	
      IF(NBC_read.NE.NBC) THEN
           WRITE(6,"(' El numero de punts llegits:',I7,' no coincideix',
     &         'amb 2(n+1)+2(m-1)=',I7//'PROGRAMA ATURAT')")NBC_read,NBC
           STOP
      ENDIF
      CLOSE(4)

      RETURN
      END SUBROUTINE READ_BC

CC ********************************************************************
CC ********************************************************************

      SUBROUTINE vertical_strain_rate (AX, BY, n, m, nn, u, v, epuntzz) 

CC    Calculation of the vertical strain rate and control that 
CC		it should be lower than ep_limit.

CC  epuntzz [1/s]  Array returning the vertical strain rate.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ep_limit=5.D-15)
      DIMENSION u(nn), v(nn), epuntzz(nn)

      Dx=AX/n
      Dy=BY/m
      epuntzz_min=1.D30
      epuntzz_max=-1.D30
     
       DO 42 iy=1,m-1
         DO 42 ix=1,n-1
	    	kxy=ix+1+iy*(n+1)
		ux=(u(kxy+1)-u(kxy-1))/(2.D0*Dx)
		vy=(v(kxy+n+1)-v(kxy-n-1))/(2.D0*Dy)
		epuntzz(kxy)=-(ux+vy)
          	epuntzz_max=MAX(epuntzz(kxy),epuntzz_max)
              	epuntzz_min=MIN(epuntzz(kxy),epuntzz_min)
 42   CONTINUE

CC    CONTINIUM VERTICAL STRAIN RATE TO THE BOUNDARIES
             DO 35 ix=1,n-1
                 kxys=ix+1
                 kxyn=ix+1+m*(n+1)
                 epuntzz(kxys)=epuntzz(kxys+n+1)
                 epuntzz(kxyn)=epuntzz(kxyn-n-1)
 35          CONTINUE
             DO 37 iy=0,m
                 kxyw=1+iy*(n+1)
                 kxye=n+1+iy*(n+1)
                 epuntzz(kxyw)=epuntzz(kxyw+1)
                 epuntzz(kxye)=epuntzz(kxye-1)
 37          CONTINUE
CC	LIMIT TO THE VERTICAL STRAIN RATE
          DO 125 kxy=1,nn
C  		limit the calculated vert.str.rate 'epuntzz'.
    		IF(epuntzz(kxy).GT.ep_limit) epuntzz(kxy)=ep_limit
            	IF(epuntzz(kxy).LT.(-1.D0*ep_limit)) 
     +			epuntzz(kxy)=-1.D0*ep_limit
          	epuntzz_max=MAX(epuntzz(kxy),epuntzz_max)
              	epuntzz_min=MIN(epuntzz(kxy),epuntzz_min)
 125	CONTINUE   
 
      WRITE(6,61) epuntzz_min,epuntzz_max,ep_limit
 61   FORMAT(4X,'Vertical strain rate: minim:',1P,G12.4,
     +  ' s-1,    maxim:',1P,G12.4,' s-1,  limit:',1P,G9.2,' s-1')

      RETURN
      END SUBROUTINE vertical_strain_rate

CC **********************************************************************
CC **********************************************************************

       SUBROUTINE thicken (Dt, n, m, nn, AX, BY, u, v, epuntzz, 
     +			       thickness, IBC_thicken)

!  New thickness of a layer with velocity field (u,v) after a time 
!	interval Dt(s)
!  (u,v) (m/s) : array of horizontal velocities.
!	d(thickness)/dt=
!		thickness*epuntzz-(u*d(thickness)/dx+v*d(thickness)/dy)

!      Amb els logaritmes dona problemes quan algun gruix es zero:
!	d[ln(thickness)]/dt=
!		epuntzz-(u*d[ln(thickness)]/dx+v*d[ln(thickness)]/dy)

! IBC_thicken=0 -> No temporal thickening variations on the boundaries,
!			thickness=thickness_old
! IBC_thicken!=0 -> No lateral variations of the thickening on the boundaries,
!			d(thickness)/dx=d(thickness)/dy=0


       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION u(nn),v(nn), epuntzz(nn), 
     +		thickness(nn),thickness_old(nn)

      Dx=AX/n
      Dy=BY/m
      thickness_old(:)=thickness(:)

      Column: DO iy=1,m-1
         Row: DO ix=1,n-1
            kxy=ix+1+iy*(n+1)
            TERME1=thickness_old(kxy)*epuntzz(kxy)
            sx=(thickness_old(kxy+1)-thickness_old(kxy-1))/(2.D0*Dx)
            sy=(thickness_old(kxy+n+1)-thickness_old(kxy-n-1))/(2.D0*Dy)
!!           TERME2=0.D0 		! LAGRANGIAN sist. (grid moves with the material)
            TERME2=u(kxy)*sx+v(kxy)*sy	! EULERIAN sist. (grid fixed)
            thickness(kxy)=thickness_old(kxy)+Dt*(TERME1-TERME2)	    
	    IF(thickness(kxy).LT.0.0) thickness(kxy)=0.D0
         END DO Row
      END DO Column
!  Boundary Conditions     
      IF(IBC_thicken==0) THEN  ! ix=0,ix=n,iy=0,iy=m -> thickness=thickness_old
	  DO iy=1,m-1
             kxy0=1+iy*(n+1)
             kxyn=n+1+iy*(n+1)
             thickness(kxy0)=thickness_old(kxy0)
             thickness(kxyn)=thickness_old(kxyn)
          END DO
          DO ix=0,n
             kxy0=ix+1
             kxym=ix+1+m*(n+1)
             thickness(kxy0)=thickness_old(kxy0)
             thickness(kxym)=thickness_old(kxym)
          END DO
       ELSE   ! ix=0,ix=n,iy=0,iy=m -> d(thickness)/dx=d(thickness)/dy=0
          DO iy=1,m-1
             kxy0=1+iy*(n+1)
             kxyn=n+1+iy*(n+1)
             thickness(kxy0)=thickness(kxy0+1)
             thickness(kxyn)=thickness(kxyn-1)
          END DO
          DO ix=0,n
             kxy0=ix+1
             kxym=ix+1+m*(n+1)
             thickness(kxy0)=thickness(kxy0+n+1)
             thickness(kxym)=thickness(kxym-n-1)
          END DO
      ENDIF 
CC --------------------------------------------------------------------
      thickmin=MINVAL(thickness)
      thickmax=MAXVAL(thickness)
      WRITE(6,65) thickmin,thickmax
 65   FORMAT(4X,'Thickness:',12X,'minim:',F10.2,' m',9X,
     +        'maxim: ',F10.2,' m')

      RETURN
      END SUBROUTINE thicken
