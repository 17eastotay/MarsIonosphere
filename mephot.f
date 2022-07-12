      PROGRAM MEPHOTMAIN
C     ===============  EPHOT  ========================================
C     040207
C     This program is now modified to work for Mars
C     Ina Robertson
C
C     Program for the calculation of the photoelectron production energy
C     spectrum, at specified altitudes.
C     CEPHOT: for a cometary atmosphere.
CC
C     150509
C     Modified a switching of Radial case and Parabolic case.
C     Changed photon flux data from Viking data to HESSR data.
C     (HESSR: http://www.galactitech.net/hessrdata/Mars/Spectra/)
C     Modified bin number of wave length for HESSR.
C     Phfluxbin0p5a and Phfluxbin1nm are connected in 6 nm.
C     (I used Phfluxbin0p5a below 6 nm, and Phfluxbin1nm above it.)
C     The unit of photon flux is /cm2/s/1nm.
C     (Note that same unit is used in Phfluxbin0p5a. You need to multiply
C     solar flux by 0.05 below 6 nm.)
C
C     Shotaro Sakai
C
C     150522
C     Modefied auger electron part.
C     I increased the number of peak of oxygen and nitrogen auger electron
C     to three.
C     COEF(I) is newly introduced for energies of auger electron.
C
C     Shotaro Sakai
C
C
C     21/08/09
C     Added other mid-major neutral species, as well as Nitrogen chemistry
C     O2, NO, N, Ar, He, H
C     only account for simple primary ionization (NO+, N+, Ar+, He+, H+)
C     DxP=Ix*nx*e^(-tau_ave)
C
C     Antonio Renzaglia, Taylor Easton
C
C     IMPLICIT NONE
CGL  -----------   parameter declaration. -----------------------
      PARAMETER(MXE=200,MXL=132,MXZ=500,MXN=11,MXI=8,MXA=5,MXC=310)
      INTEGER SW,choice,NNN(MXN),ll 	! switch for energy bin structure.

      REAL Z(MXZ),ZKM(MXZ),N2(MXZ),CO2(MXZ),CO(MXZ),O(MXZ),HE(MXZ),
     &     n_NO(MXZ),n_N(MXZ),n_Ar(MXZ),n_H(MXZ),n_O2(MXZ),n_C(MXZ),
     $     TN(MXZ),TDEN(MXZ),
     &     EXT(MXZ),COLUMN(MXN),PROD(MXN),XN(MXN,MXZ),
     &     PHE1(MXL),PHE2(MXL),ZFLUX(MXL),FLUX(MXL),RAD_cs(MXL),
     &     TPOT(MXN,MXI),SIGABS(MXN,MXL),SIGION(MXN,MXL),
     &     PROB(MXN,MXI,MXL),ENER(MXE),DEL(MXE),SPECT(MXE),
     &
     &     ELSPEC(MXZ,MXE),PRION(MXN,MXI),
     &     EPSIL1(MXN,MXI,MXL),EPSIL2(MXN,MXI,MXL),PEHEAT(MXN,MXI),
     &
     &     AGSIG(MXA,MXL),AGLW(MXA,MXZ),TAUFAC(MXZ),CHI,
     &     N2BR(MXN,MXL),CO2BR(MXN,MXL),COBR(MXN,MXL),OBR(MXN,MXL),
     &     O2BR(MXN,MXL),NOBR(MXN,MXL),wavelen(MXC),AGCS(MXN,MXL),
     &     PHE0(MXL),PHE01(MXL),PHE02(MXL),ZFLUX0(MXL)
      CHARACTER*40 FLNM12,FLNM4,FLNM18,FLNM14,FLNM20
      CHARACTER*1 YN

C     COMET VARIABLES
      REAL DHEL, QZ, FIONZ, TEMP, NEUTDENS, COEF(5)
      INTEGER COUNT,PARA
      REAL S(500),ramagl(500),sza(500),rad(500),rkm1(500)
      REAL sza1,DELS,Rm,hsubs,ramdir,cosx,arclth,radius,chi2
      REAL D_N2(MXZ),D_CO2(MXZ),D_CO(MXZ),D_O(MXZ),D_HE(MXZ)
      REAL D_NO(MXZ),D_N(MXZ),D_O2(MXZ),D_Ar(MXZ),D_H(MXZ),D_C(MXZ)
      real n2_b,n2_o, TAUAVE(MXZ), IO2,INO,IN,IAr,IH,IHe
      real zt(mxz),tnm(mxz)
      DHEL = 1.524

C     tauave is average optical thickness at a specific altitude
C     Ionization rate [1/s] Ix of species x

C----------------------------------------------------------------------------
C     NOTE: This programm is dimensioned for seven species but is set up for the
C     following two species only N2, CH4.  CORRECTION:  THIS WAS FOR TITAN, but
C     for Mars the number of species has been changed to 5:  N2, CO2, CO, O, HE
C     ADDED 5 more neutrals: NO, N, O2, Ar, H
C     ADDED another neutral (C) and photoabsorbtion/ionization cross sections (O2, N, C, H, NO, Ar)
C--------------------------------------------------------------------------

C     read the input parameters ----------------------

      PRINT*, 'PROGRAM = MEPHOT'
      PRINT*, 'CALCULATES THE PRIMARY PRODUCTION RATES OF'
      PRINT*, 'THE MAJOR IONS'


      OPEN(5,FILE='MEPBRNEW.DAT', STATUS='OLD') !photoionization cross sections
      OPEN(10,FILE='AIRPHTNEW.DAT', STATUS='OLD') !airglow cross sections (For N2)
      OPEN(12,FILE='neutraldens_reordered.txt',STATUS='OLD') !! 12: neutral density
      OPEN(14,FILE='PhotonFlux/gj699_flare_matchmaven_resortedAU.txt',
     &     STATUS='UNKNOWN') !! 14: solar flux
      open(200,file='CrossSection/wavelength.txt',status='old')
      open(201,file='CrossSection/absorption_apr22.txt',status='old')
      open(202,file='CrossSection/ionization_apr22.txt',status='old')
      open(203,file='CrossSection/n2branching_fism.txt',status='old')
      open(204,file='CrossSection/co2branching_fism.txt',status='old')
      open(205,file='CrossSection/cobranching_fism.txt',status='old')
      open(206,file='CrossSection/obranching_fism.txt',status='old')
      open(207,file='CrossSection/airglow_fism.txt',status='old')
      open(208,file='tpot.txt',status='old')
      open(209,file='mtgcm_neutraltemp.txt',status='old')
      open(210,file='CrossSection/o2branching.txt',status='old')
      open(211,file='CrossSection/nobranching.txt',status='old')

c      OPEN(4, FILE='pespect_m.dat',STATUS='UNKNOWN') !! 4:pespect.dat
c      OPEN(16, FILE='PARAM.CHK',STATUS='UNKNOWN')
      OPEN(8, FILE='EBIN.CHK', STATUS='UNKNOWN')
c      OPEN(9, FILE='MOREDAT',STATUS='UNKNOWN')
      OPEN(18,FILE='primprod.dat',STATUS='UNKNOWN')	!! 18: pho.ion prod
c      OPEN(20,FILE='elecprod.dat',STATUS='UNKNOWN') !! 20: elect production
c      OPEN(50,FILE='TEST.DAT',STATUS='UNKNOWN')
      open(100,file='peprod.txt')
      open(24,FILE='checktauave.dat',STATUS='UNKNOWN')    !!check tauave value by altitude
      open(26,FILE='checkwavelength.dat',STATUS='UNKNOWN')  !!check what wavelengths are used
      open(28,FILE='checkflux.dat',STATUS='UNKNOWN')      !!check flux vs wavelength
      open(29,FILE='checkabsorb.dat',STATUS='UNKNOWN')    !!check absorption CS vs wavelength
      open(30,FILE='checkioniz.dat',STATUS='UNKNOWN')     !!check ionization CS vs wavelength
      open(31,FILE='calcionizrate.dat',STATUS='UNKNOWN')  !!outputs for calculating ionization rates
C     PRINT*, 'DONE CREATING FILES'
      open(110,file='neutral_temperature.txt',status='unknown')
C     *** END MODIFICATION

      DATA NNN/5,7,6,6,1,5,1,1,1,4,1/
      DATA DEL/20*0.5,20*0.5,20*1.0,20*3.0,20*5.0,20*10.0,20*10.0,
     &     20*20.0,20*50.0,20*400.0/
      DATA COEF/0.91,0.93,0.93,0.93,0.9/
C     STARTIM=TIMER()
C     CPUSPT1=STARTIM

      print *,'Step 1: input from files MEPBRNEW.DAT '
C      READ(5,2005) SW,IMAX,JMAX,LMAX,LMIN,NX1,NAX
      READ(5,*) SW,IMAX,JMAX,LMAX,LMIN,NX1,NAX
      READ(5,201) Z(1),DELZ,RT,GT
      READ(5,617) NBINS
      PRINT *,'NBINS = ',NBINS
      PRINT*, IMAX
      PRINT*, JMAX
      PRINT*, LMAX
c      WRITE(4,200) CHI,SW,IMAX,JMAX,LMAX,LMIN
c      WRITE(4,209) Z(1),DELZ,RT,GT
c      WRITE(16,200) CHI,SW,IMAX,JMAX,LMAX,LMIN
c      WRITE(16,209) Z(1),DELZ,RT,GT
c      WRITE(18,200) CHI,SW,IMAX,JMAX,LMAX,LMIN
c      WRITE(18,209) Z(1),DELZ,RT,GT
c      WRITE(18,9999) FLNM14
c      WRITE(18,9999) FlNM12
      WRITE(18,210)


 210  FORMAT(' Z(km)    N+,    N2+,    CO2+,    CO+,     O+',
     &     ',      C+,     He+,     NO+,   O2+,    Ar+',
     &     ',      H+,    N++,    CO2++,    CO++,    O++')

C     CHI: is zenith angle in radians;
C     SW:  is switch: for total ionization only set greater than two;
C     For equal 1ev boxes set SW=0.  For unequal boxes, 0.5ev up to 10ev
C     and 1ev boxes for energies greater than 10ev - set SW=1 <== ***
C     IMAX: =number of species considered;
C     LMAX: =wavelength increments;
C     Z(1): is lowest altitude;
C     DELZ: is the altitude increment
C     JMAX:= no of altitude grids to be used in the calculations;
C     NX1 & NAX: 1st grid no. and the no. of grids skipped (NAX-1) between
C     calculated pts.
C     RT : = Mars radius; GT = gravity at the surface
C     !!!be sure that densities are provided down to grazing height!!!
C     XN(1) is N2, XN(2) is CO2, XN(3) is CO, XN(4) is O, XN(5) is He.
C     XN(6) is O2, XN(7) is N, XN(8) is C, XN(9) is H, XN(10) is NO, XN(11) is Ar

c      READ(14,232)(PHE1(L),PHE2(L),ZFLUX(L),L=1,LMAX)
C  ------ input of solar photon flux
      READ(14,*)(PHE0(L),ZFLUX0(L),L=1,LMAX)
      do L=1,LMAX
         PHE01(L) = PHE0(L) + 0.5
         PHE02(L) = PHE0(L) - 0.5
         PHE01(L) = PHE01(L)*10.0
         PHE02(L) = PHE02(L)*10.0
         ZFLUX0(L) = ZFLUX0(L)*1.0e-4/1.6e-19/(12397.7/phe0(l)/10.0)
      enddo
      do l=1,lmax
         ll = lmax - (l-1)
         PHE1(l) = PHE01(ll)
         PHE2(l) = PHE02(ll)
         ZFLUX(l) = ZFLUX0(ll)
      enddo
      write(28,*) lmax
      do l=1,lmax
         write(28,*)phe1(l),ZFLUX(l)
      enddo
c      stop
c      WRITE(16,232)(PHE1(L),PHE2(L),ZFLUX(L),L=1,LMAX)
C     PHE1,PHE2 are the boundaries of a wavelength bin. ZFLUX is the solar
C     flux strength.

C  ------ Input of neutral temperature
      read(209,*)(zt(i),tnm(i),i=1,500)

C  ------ Parabora check, the length of field line for Parabolic case
      OPEN(101,FILE='parameter.txt', STATUS='OLD')
      read(101,*)PARA,chi,DELS
      S(1) = 1.0E5
      Rm = 3400.0E5
      chi = chi*3.141592/180.0
      if (para .eq. 1) then
         hsubs = 100.0E5 + Rm*cos(chi)
      else
         hsubs = 1.0E5
      endif
      DO J=2,JMAX
         S(J) = S(J-1) + DELS
      ENDDO

C  ------ initialising the branching ratios
      DO L=1,LMAX
         DO I=1,IMAX
            DO JJ=1,8
               PROB(I,JJ,L)=0.0
               EPSIL1(I,JJ,L)=0.0
               EPSIL2(I,JJ,L)=0.0
            ENDDO
c            SIGION(I,L)=-1.0
         ENDDO
      ENDDO

      DO II=1,100
         DO JJ=1,100
            ELSPEC(II,JJ)=0.0
         ENDDO
      ENDDO

C-----------------------------------------------------
C     input for absorbtion and ionization cross sections, and branch ratio.
      read(201,*)
      read(202,*)
      read(210,*)
      read(211,*)
	do l=1,lmax
      read(201,*)RAD_cs(l),(SIGABS(i,l),i=1,IMAX)
      read(202,*)RAD_cs(l),(SIGION(i,l),i=1,IMAX)
	enddo
      read(203,*)((N2BR(n,l),n=1,5),l=1,lmax)
      read(204,*)((CO2BR(n,l),n=1,7),l=1,lmax)
      read(205,*)((COBR(n,l),n=1,6),l=1,lmax)
      read(206,*)((OBR(n,l),n=1,6),l=1,lmax)
      read(207,*)(agcs(1,l),l=1,lmax)
      read(210,*)((O2BR(n,l),n=1,5),l=1,lmax)
      read(211,*)((NOBR(n,l),n=1,4),l=1,lmax)

         write(29,*)"Wavelength,N2,CO2,CO,O,He,O2,N,C,H,NO,Ar"
         write(30,*)"Wavelength,N2,CO2,CO,O,He,O2,N,C,H,NO,Ar"
      do L=1,LMAX
          WRITE(29,*)PHE1(L),(SIGABS(I,L),I=1,IMAX)
          WRITE(30,*)PHE1(L),(SIGION(I,L),I=1,IMAX)
      enddo

c      do l=1,lmax
c         write(*,*)l,sigion(2,l)
c      enddo


      DO LL=1,LMAX

C     -------            solar flux ZFLUX scaling factor 0.011 for Titan
C     NOTE:  FOR COMETS THE SCALING FACTOR IS (1AU/DHEL)**2
C     THE DISTANCE FROM MARS TO THE SUN (DHEL) IS 1.524 AU.

c         ZFLUX(LL)=ZFLUX(LL)*1.E9/DHEL/DHEL
         PHE1(LL)=12397.7/PHE1(LL)
         PHE2(LL)=12397.7/PHE2(LL)
         DO II=1,IMAX
            SIGABS(II,LL)=SIGABS(II,LL)*1.E-18
         ENDDO
      ENDDO
c      do i=1,imax
c         do l=1,lmax
c            write(1300,*)i,l,sigabs(i,l)
c         enddo
c      enddo
c      do l=1,lmax
c         write(*,*)l,sigion(1,l),sigion(2,l),sigion(3,l)
c      enddo

      DO L=1,LMAX
         DO I=1,IMAX
            SIGION(I,L)=SIGION(I,L)*1.0E-18
            IF(SIGION(I,L).LT.0.0)SIGION(I,L)=SIGABS(I,L)
         ENDDO
      ENDDO
c      do i=1,imax
c         do l=1,lmax
c            write(1100,*)i,l,sigion(i,l)
c         enddo
c      enddo
c      do l=1,lmax
c         write(*,*)l,sigion(1,l),sigion(2,l),sigion(3,l)
c      enddo
c      stop
C     ----- ionization potential for L1 final state, for the K1 ion (eV units)
      DO K1=1,IMAX
         K2=NNN(K1)
c         READ(5,92) (TPOT(K1,L1),L1=1,K2)
         READ(208,92) (TPOT(K1,L1),L1=1,K2)
c         WRITE(16,92) (TPOT(K1,L1),L1=1,K2)
      ENDDO

CG -------   Preset arrays
      DO II=1,MXA
         DO L=1,MXL
            AGSIG(II,L)=0.0
 785     ENDDO
 783  ENDDO

CG --- read from AIRPHT, the airglow cross sections.
c      READ(10,205)NEXA,LBMIN,LBMAX
c      DO L=LBMIN,LBMAX
c         READ(10,202)(AGSIG(I1,L),I1=1,NEXA)
c 787  ENDDO
      do l=1,lmax
         AGSIG(1,L) = agcs(1,l)
      enddo
c      do l=1,lmax
c         write(1200,*)l,agsig(1,l)
c      enddo
c      stop
      DO L=1,LMAX
         DO I=1,IMAX
            K1=NNN(I)
            DO K=1,K1
C     TPOT(I,K) is the ionization potential in eV
C     for the L1 final state for species K.
               EPSIL1(I,K,L)=PHE1(L)-TPOT(I,K)
               EPSIL2(I,K,L)=PHE2(L)-TPOT(I,K)
            ENDDO
         ENDDO
 93   ENDDO

c      WRITE(16,216) (SIGION(1,L),L=1, LMAX)
c      WRITE(16,216) (SIGION(2,L),L=1, LMAX)
c      WRITE(16,216) (SIGION(3,L),L=1, LMAX)
c      WRITE(16,216) (SIGION(4,L),L=1, LMAX)
c      WRITE(16,216) (SIGION(5,L),L=1, LMAX)

      DO I=1,IMAX
         NIFS=NNN(I)
C---  reading in the branching ratio for species I,
C     --    		ion final state N and wavelength interval L.
C     READ(5,206)((PROB(I,N,L),N=1,NIFS),L=1,LMAX)
C     *** MODIFIED BY CRB
C     READ(5,206)((PROB(I,N,L),N=1,NIFS),L=1,LMAX)
c         IF(NIFS.EQ.5)READ(5,206)((PROB(I,N,L),N=1,NIFS),L=1,LMAX)
c         IF(NIFS.EQ.7)READ(5,9206)((PROB(I,N,L),N=1,NIFS),L=1,LMAX)
c         IF(NIFS.EQ.6)READ(5,8206)((PROB(I,N,L),N=1,NIFS),L=1,LMAX)
c         IF(NIFS.EQ.1)READ(5,7206)((PROB(I,N,L),N=1,NIFS),L=1,LMAX)

C	Order: 1=N2+, 2=CO2+, 3=CO+, 4=O+, 5=He+, 6=O2+, 7=N+, 8=C+, 9=H+
         DO N=1,NIFS
            DO L=1,LMAX
               IF ( I .EQ. 1 ) THEN
                  PROB(I,N,L) = N2BR(N,L)
               ELSEIF ( I .EQ. 2 ) THEN
                  PROB(I,N,L) = CO2BR(N,L)
               ELSEIF ( I .EQ. 3 ) THEN
                  PROB(I,N,L) = COBR(N,L)
               ELSEIF ( I .EQ. 4 ) THEN
                  PROB(I,N,L) = OBR(N,L)
               ELSEIF ( I .EQ. 5 ) THEN
                  PROB(I,N,L) = 1.0
	       ELSEIF ( I .EQ. 6 ) THEN
		  PROB(I,N,L) = O2BR(N,L)
               ELSEIF ( I .EQ. 7 ) THEN
                  PROB(I,N,L) = 1.0
               ELSEIF ( I .EQ. 8 ) THEN
                  PROB(I,N,L) = 1.0
               ELSEIF ( I .EQ. 9 ) THEN
                  PROB(I,N,L) = 1.0
	       ELSEIF ( I .EQ. 10 ) THEN
                  PROB(I,N,L) = NOBR(N,L)
	       ELSEIF ( I .EQ. 11 ) THEN
                  PROB(I,N,L) = 1.0
               ENDIF
            ENDDO
         ENDDO
C     *** END MODIFICATION

 917  ENDDO

c      WRITE(16,206)((PROB(1,N,L),N=1,5),L=1,LMAX)
c      WRITE(16,9206)((PROB(2,N,L),N=1,7),L=1,LMAX)
c      WRITE(16,8206)((PROB(3,N,L),N=1,6),L=1,LMAX)
c      WRITE(16,8206)((PROB(4,N,L),N=1,6),L=1,LMAX)
c      WRITE(16,7206)((PROB(5,N,L),N=1,1),L=1,LMAX)
!      do i=1,imax
!         do l=1,lmax
!            write(1000,*)i,l,SIGABS(I,L)
!            write(1100,*)i,l,SIGION(I,L)
!         enddo
!      enddo
!      do l=1,lmax
!         write(1200,*)l,AGSIG(1,L)
!      enddo
!      do i=1,imax
!         NIFS=NNN(I)
!         do n=1,nifs
!            do l=1,lmax
!               write(1300,*)i,n,l,prob(i,n,l)
!            enddo
!         enddo
!      enddo
      IF(SW .LE.0) THEN
C---  This section sets up the energy bins for 195*1 eV bins
         ENER(1)=0.5
         DO M=2,NBINS
            ENER(M)=ENER(M-1)+0.5*(DEL(M-1)+DEL(M))
 1492    ENDDO
C     GO TO 738
      ELSE

CG-------  This section sets up the unequal energy bins with the following
CG-------  structure.
C
C     i=     20     90    100    120     130    150    160      170      185
C     20x.5   70x1.  10x2.  20x5.  10x10. 20x10. 10x50.  10x100.  15x200.
C     +------+------+------+------+-------+------+------+--------+--------+
C     E=0eV   10     80.   100.   200.    300.   500.   1000.    2000.    5000.
C

C     731     CONTINUE

         ENER(1)=0.25
         DO M=2,NBINS
            ENER(M)=ENER(M-1)+0.5*(DEL(M-1)+DEL(M))
         ENDDO

         WRITE(8,221)(ENER(M),M=1,NBINS)
         WRITE(8,221)(DEL(M),M=1,NBINS)
 221     FORMAT(1F8.2)

      ENDIF

      print *,' Step 2: input from file neutral atmosphere'

      READ(12,*)

      DO J=1,JMAX
C-  ....................reading in the neutral species.

C *** MODIFIED BY JCV 8/00 TO READ FROM MODIFIED NEUTRAL FILES
C
c         READ(12,615) RKM, N2(J), CO2(J), CO(J), O(J), HE(J)
	READ(12,*) RKM,N2(J),CO2(J),CO(J),O(J),HE(J),
     $     n_O2(J),n_N(J),n_C(J),n_H(J),n_NO(J),n_Ar(J)

         RKM1(J) = RKM*1.0e5 + rm*cos(chi)

C--   ....................Assign neutral temperature.
C     FOR MARS THE NEUTRAL TEMPERATURE IS SET TO 200 K.
c         TN(J)=200.
c      WRITE(50,*)ZKM(J),N2(J), CO2(J), CO(J), O(J), HE(J)
c      WRITE(16,615)ZKM(J),N2(J),CO2(J),CO(J),O(J),HE(J)

CG- ....................set up the spatial grid.
         IF (para .eq. 0) then
            ZKM(J)=Z(1)*(1.E-5)+DELZ*(1.E-5)*FLOAT(J-1)
            Z(J)=ZKM(J)*1.E5
            IF (RKM.NE.ZKM(J))PRINT *,'!!!!CHECK ZKM(J), J=',J
            write(*,*)Z(j)
            sza(J) = chi
         ENDIF
C.....................................preset airglow (AGLW)
         DO IA1=1,MXA
            AGLW(IA1,J)=0.0
         ENDDO
      ENDDO

CG- ....................set up the spatial grid for parabolic case
      DO J=1,JMAX
         ARCLTH=S(J)
         IF (PARA.EQ.1) THEN
            CALL PARA2(HSUBS,ARCLTH,ramdir,radius)
            ramagl(J) = ramdir
            rad(J) = radius
            ZKM(J) = (rad(j) - rm*cos(chi))*1.0e-5
            Z(J) = ZKM(J)*1.0E5
            chi2 = asin(Rm*sin(chi)
     &           /sqrt(Rm*Rm + 2.0*(rad(j)-Rm*cos(chi))*Rm*cos(chi)
     &           +(rad(j)-Rm*cos(chi))*(rad(j)-Rm*cos(chi))))
            cosx = cos(chi2)*cos(ramagl(J)*3.141592/180.0)
            sza(J) = acos(cosx)
            call extrap(500,rkm1,N2,rad(J),D_N2(J),1,1)
            call extrap(500,rkm1,CO2,rad(J),D_CO2(J),1,1)
            call extrap(500,rkm1,CO,rad(J),D_CO(J),1,1)
            call extrap(500,rkm1,O,rad(J),D_O(J),1,1)
            call extrap(500,rkm1,HE,rad(J),D_HE(J),1,1)
            call extrap(500,rkm1,n_O2,rad(J),D_O2(J),1,1)
            call extrap(500,rkm1,n_N,rad(J),D_N(J),1,1)
            call extrap(500,rkm1,n_C,rad(J),D_C(J),1,1)
	    call extrap(500,rkm1,n_H,rad(J),D_H(J),1,1)
            call extrap(500,rkm1,n_NO,rad(J),D_NO(J),1,1)
            call extrap(500,rkm1,n_Ar,rad(J),D_Ar(J),1,1)
            call extrap(500,rkm1,tnm,rad(J),TN(J),1,1)
         ENDIF
         write(110,*)rad(j),tn(j)
      ENDDO

      if (para .eq. 1) then
         DO J=1,JMAX
            N2(J) = D_N2(J)
            CO2(J) = D_CO2(J)
            CO(J) = D_CO(J)
            O(J) = D_O(J)
            HE(J) = D_HE(J)
	    n_O2(J) = D_O2(J)
            n_N(J) = D_N(J)
            n_C(J) = D_C(J)
	    n_H(J) = D_H(J)
            n_NO(J) = D_NO(J)
            n_Ar(J) = D_Ar(J)
         ENDDO
      endif

      LGT=1

CG=============!===========================!=====================
CG   ----------!     Begin Altitude loop.  !---------------------
C--------------!---------------------------!---------------------

      DO 7 J=NX1,JMAX,NAX 	!!!! alt. loop !!!!!!!!!!!!!!!
         SZA1 = SZA(J)
         DO I=1,IMAX
            PROD(I)=0.
            DO K=1,MXI
               PEHEAT(I,K)=0.
               PRION(I,K)=0.
            ENDDO
         ENDDO

         XN(1,J)=N2(J)
         XN(2,J)=CO2(J)
         XN(3,J)=CO(J)
         XN(4,J)=O(J)
         XN(5,J)=HE(J)
	 XN(6,J)=n_O2(J)
	 XN(7,J)=n_N(J)
         XN(8,J)=n_C(J)
	 XN(9,J)=n_H(J)
	 XN(10,J)=n_NO(J)
	 XN(11,J)=n_Ar(J)

         chi2 = asin(Rm*sin(chi)
     &       /sqrt(Rm*Rm + 2.0*(rad(j)-Rm*cos(chi))*Rm*cos(chi)
     &       +(rad(j)-Rm*cos(chi))*(rad(j)-Rm*cos(chi))))
        write(*,*)j,Z(J)*1.0e-5,XN(1,J),XN(2,J),XN(3,J),XN(4,J),XN(5,J),
     $            XN(6,J),XN(7,J),XN(8,J),XN(9,J),XN(10,J),XN(11,J),
     $            sza1*180.0/3.141592,chi2*180.0/3.141593
C--------------Call sub. to calculate column densities

         CALL RCOLUM(SZA1,IMAX,Z,JMAX,J,TN,N2,CO2,CO,
     $        O,HE,n_NO,n_N,n_O2,n_Ar,n_H,n_C,COLUMN,
     $        RT,GT,XN,LGT,PARA,CHI)

 82      DO M= 1,NBINS
            SPECT(M)= 0.0
 13      ENDDO

C------------------!----------------------------!-------------------
CG ----------------!   Starts wavelength loop   !-------------------
C  ----------------!----------------------------!-------------------
         DO 6 L=LMIN,LMAX
            TAU=0.
C-----computation of optical depth for absorption.
            DO 4 I=1,IMAX
               TAU=TAU+SIGABS(I,L)*COLUMN(I)
		write(26,*) J, L, TAU
 4          CONTINUE
            FLUX(L)=ZFLUX(L)*EXP(-TAU)
C-----!! flux decay due to absorption.

C-----find average optical thickness at each altitude
C-----want tau where wavelength is 50nm
	    IF (L == 50) THEN
		TAUAVE(J) = TAU
	    ENDIF

CCC---------    airglow volume emission rate for L calc here
CCC------------ at moment only N+/N2 is calculated
c            DO I1=1,NEXA
            DO I1=1,1
               AGLW(I1,J)=AGLW(I1,J)+N2(J)*AGSIG(I1,L)*FLUX(L)*1.0E-18
 794        ENDDO

            PTOT=0.

C----------------!-------------------------!-----------------
C----------------! Loop over species       !-----------------
C----------------!-------------------------!-----------------
            DO 304 I=1,IMAX
               NIFS=NNN(I)
               DO 302 K=1,NIFS
                  E1= EPSIL1(I,K,L)
                  E2= EPSIL2(I,K,L)
                  IF(E2.LT.0.) GO TO 302
                  IF(E1.LT.0.) GO TO 378

CG          -----------------------------------------------
                  DSPECT= XN(I,J)*SIGION(I,L)*FLUX(L)*PROB(I,K,L)
CG          -----------------------------------------------
                  GO TO 379
 378              E1=0.
                  DSPECT=(XN(I,J)*SIGION(I,L)*FLUX(L)*PROB(I,K,L))*
     &                 (E2/(E2-E1))
 379              CONTINUE 	!! go to pass !!
                  Y=E2-E1
CG          -----------------------------------------------------
CG          PRION(I,K):   production of K ion state for species I
CG          PEHEAT(I,K):  energy given to the photoelectron resulting
CG                        from the production of the K ion state for
CG                        I species.
CG          -----------------------------------------------
                  PROD(I)=PROD(I)+DSPECT
                  PRION(I,K)=PRION(I,K)+DSPECT
                  PEHEAT(I,K)=PRION(I,K)*(E1+E2)*0.5
CG             hv-IP(k) :should actually be hv - IP(k) - POT(k)
                  PTOT=PTOT+PEHEAT(I,K)

CG          -----------------------------------------------
CG 	    retrieve the energy box number for E1 and E2
CG          --------------------------------------------
                  IF(SW.GT.2) GO TO 302
                  IF(SW .GT. 0) GO TO 735
 301              CALL BOXNUM(E1,E2,M1,M2,R1,R2)
                  GO TO 736
 735              CALL BOXNN(E1,E2,M1,M2,R1,R2)
 736              IF(M1.GT.NBINS)GO TO 302
                  IF(M2.GT.NBINS)GO TO 302
C.................................     fill the boxes between M1 and M2
		  IF(M2-M1 .lt. 0) GO TO 51
		  IF(M2-M1 .eq. 0) GO TO 52
		  IF(M2-M1 .gt. 0) GO TO 53

c 51               WRITE(16, 106) L, I, K
 51					print *,'error'
c                  WRITE(8,98)L,I,K,E1,E2
                  GO TO 302
 52               SPECT(M1)= SPECT(M1)+ DSPECT
                  GO TO 302
 53               DO 60 M=M1,M2
		     IF(M-M1 .lt. 0) GO TO 43
		     IF(M-M1 .eq. 0) GO TO 54
		     IF(M-M1 .gt. 0) GO TO 55
CG                                   43: write error msg.
 54                  SPECT(M)= SPECT(M)+ DSPECT*(R1- E1)/Y
                     GO TO 60
 55 		     IF(M-M2 .lt. 0) GO TO 56
		     IF(M-M2 .eq. 0) GO TO 57
		     IF(M-M2 .gt. 0) GO TO 44
 56                  SPECT(M)= SPECT(M)+ DSPECT*DEL(M)/Y
                     GO TO 60
 57                  SPECT(M)= SPECT(M)+ DSPECT*(E2- R2)/Y
 60               CONTINUE
 302           CONTINUE 	!!! exit pass.
 304        CONTINUE            !!! end loop over I species.
 6       CONTINUE 		!!! end loop over L energies.

C Auger Electron Section
C   BE of K-shell electron for nitrogen Auger from N2: 402 eV?
C   BE of K-shell electron for oxygen Auger from CO2, CO and O: 532 eV
CC I=1 N2; I=2 CO2; I=3 CO; I=4 O; I=5 He; I=6 O2; I=7 N; I=8 C; I=9 H; I=10 NO; I=11 Ar
          DO I=1,IMAX
             NXX=NNN(I)
             EAUG1=COEF(I)*TPOT(I,NXX)
C             EAUG1=0.9*TPOT(I,NXX)
             EAUG2=EAUG1+0.1
             CALL BOXNN(EAUG1,EAUG2,M1,M2,R1,R2)
             IF (I .EQ. 5) THEN
                SPECT(M1)=SPECT(M1)+PRION(I,NXX)
             ELSE
                SPECT(M1-1)=SPECT(M1-1)+PRION(I,NXX)/3.0
                SPECT(M1)=SPECT(M1)+PRION(I,NXX)/3.0
                SPECT(M1+1)=SPECT(M1+1)+PRION(I,NXX)/3.0
             ENDIF
          ENDDO

C NOW SKIP I =1 (N2) AND DO THE SECOND AUGER ELECTRON FOR
C CO2 AND CO.  O AND HE ARE ALSO SKIPPED
CC BE of K-shell electron for carbon Auger from CO2 and CO: 284 eV

          DO I = 2, 3
             NXX=NNN(I) - 1
             EAUG1=0.9*TPOT(I,NXX)
             EAUG2=EAUG1+0.1
             CALL BOXNN(EAUG1,EAUG2,M1,M2,R1,R2)
             SPECT(M1)=SPECT(M1)+PRION(I,NXX)
          ENDDO

C----------------------------------------------
C----     Chemistry section.   ----------------
C----------------------------------------------

C ARR - Add simple ionization for NO, N, O2, Ar, H (cut to just NO and Ar - 2/2/22)
C ARR - have added cross sections for all above species, no simple ionization
C These values are for solar minimum conditions - no longer need these
C	INO = 3.9611E-07
C	IAr = 3.1517E-07

C  Ion grouping
c  revised for mars. 11/3/11

	DN2P=PROD(1)-AGLW(1,J)-PRION(1,5)	  	 		! N2+ = sum(N2+:X,A,B,S) - N+ -Ni++
	DCO2P=PRION(2,1)			         		! CO2+
	DCOP=PROD(3)-PRION(3,5)-PRION(3,6)+PRION(2,2)	 		! CO+
	DOP=PROD(4)-PRION(4,6)+PRION(2,3)+PRION(6,3)+PRION(10,2)	! O+
	DHEP=PROD(5)					 		! He+
        DO2P=PRION(6,5)                                 		! O2+
	DNP=AGLW(1,J)+PROD(7)+PRION(10,3)                     		! N+ = airglow from N2 + photoionization of N
	DCP=PROD(8)+PRION(2,4)		                 		! C+
        DHP=PROD(9)                                      		! H+
	DNOP=PRION(10,1)                                 		! NO+
        DArP=PROD(11)                                    		! Ar+

C Added Ni++ and CHi++, which is the final products due to the Auger electron emmision
	DNPP=PRION(1,5)				    !N++
	DCO2PP=PRION(2,5)+PRION(2,6)+PRION(2,7)     !CO2++
	DCOPP=PRION(3,5)+PRION(3,6)		    !CO++
	DOPP=PRION(4,6)				    !O++

      WRITE(18,1003) ZKM(J),DNP,DN2P,DCO2P,DCOP,DOP,DCP,DHEP,
     &      DNOP,DO2P,DArP,DHP,DNPP,DCO2PP,DCOPP,DOPP

      WRITE(24,*) ZKM(J), TAUAVE(J)

	IF (J == 452) THEN
		write(31,*) "Altitude(km):   ", ZKM(J)
		write(31,*) "Species A, ", "A density(1/cm2), ", "A+ PR(1/cm2/s)"
		write(31,*) "N2  ", XN(1,J), (PROD(1)-AGLW(1,J)-PRION(1,5))
		write(31,*) "CO2 ", XN(2,J), (PRION(2,1))
		write(31,*) "CO  ", XN(3,J), (PROD(3)-PRION(3,5)-PRION(3,6))
		write(31,*) "O   ", XN(4,J), (PROD(4)-PRION(4,6))
		write(31,*) "He  ", XN(5,J), (PROD(5))
		write(31,*) "O2  ", XN(6,J), (PRION(6,5))
		write(31,*) "N   ", XN(7,J), (PROD(7))
		write(31,*) "C   ", XN(8,J), (PROD(8))
		write(31,*) "H   ", XN(9,J), (PROD(9))
		write(31,*) "NO  ", XN(10,J), (PRION(10,1))
		write(31,*) "Ar  ", XN(11,J), (PROD(11))
	ENDIF

      DO M=1,NBINS
         IM=NBINS+1-M
         ELSPEC(J,IM)=SPECT(M)
      ENDDO
c      WRITE(9,100) J,Z(J)
c      WRITE(9,108)(PROD(I),I=1,IMAX)

c      ENUM=0.0
c      ENERGY=0.0
c      DO M=1,NBINS
c         ENERGY=ENERGY+ENER(M)*SPECT(M)
c         ENUM=ENUM+SPECT(M)
c      ENDDO
C * * * this if is for when zenith angle ge to 90 deg.  * * *
c      IF(ENUM.LT.1.E-35) GO TO 6123
c      EBAR=ENERGY/ENUM
c      GO TO 6124
c6123  EBAR=1.E35
c6124  CONTINUE

CG======================================================
CG  ======   output section               ==============
CG------------------------------------------------------
C      WRITE(9,1463)
C      WRITE(9,1464)(I,(PRION(I,K),K=1,NNN(I)),I=1,IMAX)
C      WRITE(9,1464)(I,(PEHEAT(I,K),K=1,NNN(I)),I=1,IMAX)

c         WRITE(50,1464) Z(J)*1E-5,PROD(3),(PRION(3,K),K=1,NNN(3))

C	*** MODIFIED BY CRB 7/99
C	WRITE(9,1464) 1,(PRION(1,K),K=1,NNN(1))
C	WRITE(9,9464) 2,(PRION(2,K),K=2,NNN(2))
C	WRITE(9,1464) 1,(PEHEAT(1,K),K=1,NNN(1))
C	WRITE(9,9464) 2,(PEHEAT(2,K),K=1,NNN(2))
C	*** END MODIFICATION

C      WRITE(9,1494) ENUM,ENERGY,EBAR,PTOT,FACCC

C      CPUSPT2=TIMER()
C      CPUSPT=CPUSPT2-CPUSPT1
C      CPUSPT1=CPUSPT2
C      IF (MOD(J,20).EQ.1) PRINT *,' J CPU(S):',J,CPUSPT
CG  ----------
    7 CONTINUE


C      DO K = 1, NNN(2)
C         WRITE(50,*) 'CO2', K, PROD(2),PRION(2,K)
C      ENDDO
C      DO K = 1, NNN(3)
C         WRITE(50,*) 'CO', K, PROD(3),PRION(3,K)
C      ENDDO

CG  ========== End J loop, the spatial grid ==============


C	WRITE(18,321)(TAUFAC(J),J=1,JMAX)
c      WRITE(4,401)JMAX,NBINS
C	WRITE(6,1601)
C 1601   FORMAT(' Write large data set to 4: PESPECT',/5X,$)
      DO L=1,NBINS
c        WRITE(4,*) (ELSPEC(J,L),J=1,JMAX)
        WRITE(100,*) (ELSPEC(J,L),J=1,JMAX)
	   IF (MOD(L,20).EQ.0) THEN
		PERCENT=100.0*L/FLOAT(NBINS)
C		WRITE(6,1602)PERCENT
 1602		FORMAT('+',F4.0,'% ',$)
	   ENDIF
      ENDDO

      PRINT *,'OUTPUT: 4.PESPECT    18.PHOION    9.MOREDAT'
      PRINT *,'        6.PARAM.CHK  8.EBIN.CHK'

      STOP

CG===================================================
CG    write error message
CG    -------------------
c 42    WRITE(8,98)I,K,L
 42			print *,'error'
c      GO TO 302
c 43    WRITE(8,99)I,K,L,E1,E2
 43			print *,'error'
c      GO TO 302
c 44    WRITE(8,45)I,K,L,E1,E2
 44			print *,'error'
c      GO TO 302
45    FORMAT('ERROR 3',2X,I3,2X,I3,2X,I3,2X,E10.3,2X,E10.3)
92    FORMAT(8F7.2)
98    FORMAT(' ERROR 1',2X,I3,2X,I3,2X,I3,2E13.4)
99    FORMAT(' ERROR 2',2X,I3,2X,I3,2X,I3,2X,E10.3,2X,E10.3)

CG ======= format list
  100 FORMAT(/,3X,'LEVEL',I3,9X,'HEIGHT',1PE13.2,' CM.')
  101 FORMAT(10X,'CHI= ',F4.1,10X,'COLUMN DENSITIES',//8X,'CO2',1PE11.3,
     &7X,'CO',1PE11.3,7X,'O',1PE11.3 ,7X,'N2',1PE11.3,7X,'HE',1PE11.3//)
  102 FORMAT(('0', I3, 1PE12.4, 3(I6, E12.4)))
  103 FORMAT('0', 10X, 'COMPUTATION COMPLETE - LOWER BOUNDARY HAS NOT BE
     &EN REACHED' )
  104 FORMAT('0', 10X, 'NORMAL COMPLETION OF COMPUTATION' )
  105 FORMAT(35X,'PHOTOELECTRON ENERGY SPECTRUM'//('0',8(I4, 1PE11.3)))
  106 FORMAT('0', 15X, 'L=', I3, 5X, 'I=', I2, 5X, 'K=', I2, 5X,
     &'DATA ERROR')
  107 FORMAT(8E10.4)
 108  FORMAT(1X,'ELECTRON PRODUCTION RATES FOR:'/'  N2',1PE15.7,'  CH4'
     &, 1PE15.7)
  200 FORMAT(F6.3,7I4)
 2005 FORMAT(7I4)
  201 FORMAT(4E10.3)
  202 FORMAT(10F8.3)
  203 FORMAT(2(3I3, 2F9.2, F8.3))
  204 FORMAT(8F7.2)
 205  FORMAT(1X,10I5)
C  206 FORMAT(1X,<NIFS>F7.3)
C	*** MODIFIED BY CRB
 206  FORMAT(1X,5F7.4)
 7206 FORMAT(1X,1F7.4)
 8206 FORMAT(1X,6F7.4)
 9206 FORMAT(1X,7F7.4)
C	*** END MODIFICATION
  207 FORMAT(1X,F6.2,3E10.3)
  209 FORMAT(1P4E10.3)
C  210 FORMAT(3X,'J Z(km)',5X,'N2+',5X,'CH4+',5X,'CH3+',5X,'CH2+',5X,'CH+',
C	& 5X,'H+',5X,'H2+',5X,'N+',5X,'C2H+',5X,'C2H2+',5X,'C2H3+',
C	& 5X,'C2H4+',5X,'C2H6+')
  212 FORMAT(1X,5F8.4)
  216 FORMAT(8E9.2)
C 232 FORMAT(2F8.2,F10.4)
c  232 FORMAT(2F9.2,E15.3)                              !USED FOR NEWSFLX.DAT
  232 FORMAT(F9.2,E15.3)                              !USED FOR NEWSFLX.DAT
  321 FORMAT(/' PHOTOABSORPTION SCALING FACTOR J=1,JMAX',/,(1p8e10.3))
 401   FORMAT(2I6)
 402   FORMAT(1X,1P10E11.3)
 9402  FORMAT(1X,7E11.3)
 614  FORMAT(1X,F4.0,10E9.3)
  615 FORMAT(1X,F7.0,2X,1P5E11.3)
 617  FORMAT(I4)
 678  FORMAT(2X,1PE12.3)
  874 FORMAT(//30X,14HSPECIES NUMBER,2I6,//)
  875 FORMAT(/,2X,'ALTITUDE *','VOLUME EMISSION RATES IONIZATION NONE')
  879 FORMAT(5X,F6.1,1P8E11.3)
 1003 FORMAT(1X,F7.3,1P15E9.2)
 1004 FORMAT(1X,1P7E9.2)
 1463 FORMAT(5X,'excited ion prod. rates & energy given to the e-')
 1464 FORMAT(1X,F6.0,7E9.2)
 9464 FORMAT(10X,I4,5X,8E13.6)
 1494 FORMAT(5X,'ENUM=',1PE10.2,'  ENERGY=',1PE10.2,'  EBAR=',
     &   1PE10.2,' TOTAL HEAT INPUT=',1PE10.2,' FACCC=',1PE10.2)
 1601 FORMAT(' WRITE LARGE DATA SET TO 4: PESPECT',/5X,$)

        END

C==============================================================
      SUBROUTINE RCOLUM(CHI,IMAX,Z,JMAX,J,TN,N2,CO2,CO,O,HE,NO,
     $     N,O2,Ar,H,C,COLUMN,RT,GT,XN,LGT,PARA,CHI0)

C  This subroutine is set for all five neutrals. The part for SZA greater
C  than 90 deg. is not used. In case it is needed, reexam that part.

      DIMENSION Z(JMAX),TN(JMAX),N2(JMAX),CO2(JMAX),CO(JMAX),O(JMAX),
     $     HE(JMAX),NO(JMAX),N(JMAX),O2(JMAX),Ar(JMAX),H(JMAX),C(JMAX),
     $     COLUMN(IMAX),XN(IMAX,JMAX),VERTCL(11),GHN(7)
      REAL CHI,Z,TN,XN,COLUMN,RT,GT,N2,CO2,CO,O,HE,NO,N,O2,Ar,H,C,
     $     CHI0,HEI
      INTEGER IMAX,JMAX,J,LGT,PARA

      DATA PI/3.141592654/
      COLUMN(1)=1.E35
      COLUMN(2)=1.E35
      COLUMN(3)=1.E35
      COLUMN(4)=1.E35
      COLUMN(5)=1.E35
      COLUMN(6)=1.E35
      COLUMN(7)=1.E35
      COLUMN(8)=1.E35
      COLUMN(9)=1.E35
      COLUMN(10)=1.E35
      COLUMN(11)=1.E35

      if (para .eq. 0) then
         HEI = RT+Z(J)
      elseif (para .eq. 1) then
         HEI = sqrt(RT*RT + 2.0*Z(J)*RT*cos(CHI0) + Z(J)*Z(J))
      endif
c      write(*,*)RT+Z(J),HEI
C  ***  MODIFIED BY CRB
C      IF(CHI.GT.1.964) GO TO 41
C      IF CHI .GT. 116 DEG, THEN DONT DO ANY CALCULATIONS
      IF (CHI.LT.2.025) THEN
C  ***  END MODIFICATIONS
         RHO=(N2(J)*28.+CO2(J)*44.+CO(J)*28.+O(J)*16.+HE(J)*4.+
     $        NO(J)*30.+N(J)*14.+O2(J)*32.+Ar(J)*40.+H(J)+C(J)*6.)
     $        *1.662E-24
c         GR=GT*(RT/(RT+Z(J)))**2
         GR=GT*(RT/HEI)**2
         HN=1.38E-16*TN(J)*(N2(J)+CO2(J)+CO(J)+O(J)+HE(J)+NO(J)+
     $         N(J)+O2(J)+Ar(J)+H(J)+C(J))/(RHO*GR)
C                HN is the scale height for the mean mass.
c         HG=(RT+Z(J))/HN
         HG=HEI/HN
         HF=0.5*HG*(COS(CHI)**2)
         TT=SQRT(HF)
         HN2=(1.38E-16*TN(J))/(28.*1.662E-24*GR)

C *** MODIFIED BY CRB
         HCO2=HN2*28.0/44.0
         HCO =HN2
         HO = HN2*28.0/16.0
         HHE= HN2*28.0/4.0
         HNO= HN2*28.0/30.0
         HN = HN2*28.0/14.0
         HO2= HN2*28.0/32.0
         HAr= HN2*28.0/40.0
         HH = HN2*28.0/1.0
         HC = HN2*28.0/6.0
C
C      IF(CHI-1.571) 20,20,30
C  20  SECZI=SQRT(0.5*PI*HG)*SPERFC(TT)

         VERTCL(1)=N2(J)*HN2
         VERTCL(2)=CO2(J)*HCO2
         VERTCL(3)=CO(J)*HCO
         VERTCL(4)=O(J)*HO
         VERTCL(5)=HE(J)*HHE
         VERTCL(6)=NO(J)*HNO
         VERTCL(7)=N(J)*HN
         VERTCL(8)=O2(J)*HO2
         VERTCL(9)=Ar(J)*HAr
         VERTCL(10)=H(J)*HH
	 VERTCL(11)=C(J)*HC

         IF (CHI.LE.1.571) THEN ! FOR CHI LESS THAN 90 DEG
 20         SECZI=SQRT(0.5*PI*HG)*SPERFC(TT)
            COLUMN(1)=VERTCL(1)*SECZI
            COLUMN(2)=VERTCL(2)*SECZI
            COLUMN(3)=VERTCL(3)*SECZI
            COLUMN(4)=VERTCL(4)*SECZI
            COLUMN(5)=VERTCL(5)*SECZI
            COLUMN(6)=VERTCL(6)*SECZI
            COLUMN(7)=VERTCL(7)*SECZI
            COLUMN(8)=VERTCL(8)*SECZI
            COLUMN(9)=VERTCL(9)*SECZI
            COLUMN(10)=VERTCL(10)*SECZI
	    COLUMN(11)=VERTCL(11)*SECZI
         ELSE                   ! FOR CHI GREATER THAN 90 DEG
c            GHRG=(RT+Z(J))*SIN(CHI) !R_0, RADIUS AT 90 DEG
            GHRG=HEI*SIN(CHI) !R_0, RADIUS AT 90 DEG
            GHZ=GHRG-RT
            IF (GHZ.GT.0.) THEN
               HF0=0.5*HG*(COS(PI-CHI)**2)
               TT0=SQRT(HF0)
               SECZI0=SQRT(0.5*PI*HG)*SPERFC(TT0)
               SECZI=SQRT(0.5*PI*GHRG/HN)*SPERFC(0.)
c               VAR1=EXP((RT+Z(J)-GHRG)/HN2)
c               VAR2=EXP((RT+Z(J)-GHRG)/HCO2)
c               VAR3=EXP((RT+Z(J)-GHRG)/HCO )
c               VAR4=EXP((RT+Z(J)-GHRG)/HO)
c               VAR5=EXP((RT+Z(J)-GHRG)/HHE)
               VAR1=EXP((HEI-GHRG)/HN2)
               VAR2=EXP((HEI-GHRG)/HCO2)
               VAR3=EXP((HEI-GHRG)/HCO)
               VAR4=EXP((HEI-GHRG)/HO)
               VAR5=EXP((HEI-GHRG)/HHE)
               VAR6=EXP((HEI-GHRG)/HNO)
               VAR7=EXP((HEI-GHRG)/HN)
               VAR8=EXP((HEI-GHRG)/HO2)
               VAR9=EXP((HEI-GHRG)/HAr)
               VAR10=EXP((HEI-GHRG)/HH)
	       VAR11=EXP((HEI-GHRG)/HC)
               COLUMN(1)=VERTCL(1)*(2.*SECZI*VAR1-SECZI0)
               COLUMN(2)=VERTCL(2)*(2.*SECZI*VAR2-SECZI0)
               COLUMN(3)=VERTCL(3)*(2.*SECZI*VAR3-SECZI0)
               COLUMN(4)=VERTCL(4)*(2.*SECZI*VAR4-SECZI0)
               COLUMN(5)=VERTCL(5)*(2.*SECZI*VAR5-SECZI0)
               COLUMN(6)=VERTCL(6)*(2.*SECZI*VAR6-SECZI0)
               COLUMN(7)=VERTCL(7)*(2.*SECZI*VAR7-SECZI0)
               COLUMN(8)=VERTCL(8)*(2.*SECZI*VAR8-SECZI0)
               COLUMN(9)=VERTCL(9)*(2.*SECZI*VAR9-SECZI0)
               COLUMN(10)=VERTCL(10)*(2.*SECZI*VAR10-SECZI0)
	       COLUMN(11)=VERTCL(11)*(2.*SECZI*VAR11-SECZI0)
            ENDIF
         ENDIF
      ENDIF

C      DO 21 I=1,IMAX
C   21 COLUMN(I)=VERTCL(I)*SECZI
C      GO TO 40
C
C  30  GHRG=(RT+Z(J))*SIN(CHI)
C      GHZ=GHRG-RT
C      IF(GHZ) 41,41,50
C  50  GHG=GT*(RT/(RT+GHZ))**2.
C      IF(GHZ.GT.121.E5) GO TO 60
C      GHTN=300.-1.4*GHZ*1.E-5
C      GHHG=1.38E-16*GHTN/7.3E-23/GHG
C      GHNT=N2(1)*EXP((120.E5-GHZ)/GHHG)
C      GHN(1)=.9992*GHNT
C      GHN(2)=.0004*GHNT
C      GHXG=GHRG/GHHG
C      DO 22 I=1,IMAX
C  22  COLUMN(I)=SQRT(0.5*PI*GHXG)*GHHG*(2.*GHN(I)-XN(I,J)*SPERFC(TT))
C      GO TO 40
C
C  ***  THIS IS FOR GHZ .GT. 120 KM.  ***
C  60  IF(CHI.GT.1.919) GO TO 42
C ** FOR CHI=100 DEG., PT. P IS 24 LEVELS HIGHER THAN PT. G **
C      K=LGT+24
C      GO TO 43
C  42  K=LGT+94
C  43  COLUMN(1)=SQRT(0.5*PI*GHRG/HN2)*HN2*(2.*N2(LGT)-N2(K)*
C     K          SPERFC(TT))
C      COLUMN(2)=SQRT(0.5*PI*GHRG/HCH4)*HCH4*(2.*CH4(LGT)-
C     &          CH4(K)*SPERFC(TT))
C      LGT=LGT+1
C      GO TO 40
C  41  DO 23 I=1,IMAX
C  23  COLUMN(I)=1.E35
C   40 CONTINUE

C  *** END MODIFICATION
      RETURN
      END

C==========================================
      SUBROUTINE BOXNUM(E1, E2,M1,M2,R1,R2)
C        this subroutine finds the box numbers corresponding to the
C    energies E1 and E2, and CALLS them M1 and M2
C  this subroutine sets up one ev boxes between 1 and 100ev
   31 M1= E1+ 1.0
      M2= E2+ 1.0
      R1=M1
      R2=(M2-1)
      RETURN
      END

C=========================================================================
      SUBROUTINE BOXNN(E1,E2,M1,M2,R1,R2)
C        This subroutine finds the energy bin numbers in which
C   energies E1 and E2 reside, and calls them M1 and M2. R1 is the
C   upper boundary energy of bin M1. R2 is the lower boundary energy
C   of bin M2.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	IF (E1.LT.10.0) THEN
	  M1=1+E1/0.5
	  R1=M1*0.5
	ELSE IF (E1.LT.20.0) THEN
	  M1=21+(E1-10.)/0.5
	  R1=10.+(M1-20.)*0.5
	ELSE IF (E1.LT.40.0) THEN
	  M1=41+(E1-20.)/1.
	  R1=20.+(M1-40)*1.
	ELSE IF (E1.LT.100.0) THEN
	  M1=61+(E1-40.)/3.
	  R1=40.+(M1-60)*3.
	ELSE IF (E1.LT.200.0) THEN
	  M1=81+(E1-100.)/5.
	  R1=100.+(M1-80)*5.
	ELSE IF (E1.LT.400.0) THEN
      M1=101+(E1-200.)/10.
	  R1=200.+(M1-100)*10.
	ELSE IF (E1.LT.600.0) THEN
	  M1=121+(E1-400.)/10.
	  R1=400.+(M1-120)*10.
	ELSE IF (E1.LT.1000.0) THEN
	  M1=141+(E1-600.)/20.
	  R1=600.+(M1-140)*20.
	ELSE IF (E1.LT.2000.0) THEN
	  M1=161+(E1-1000.)/50.
	  R1=1000.+(M1-160)*50.
	ELSE
	  M1=181+(E1-2000.)/400.
	  R1=2000.+(M1-180)*400.
	ENDIF

	IF (E2.LT.10.0) THEN
	  M2=1+E2/0.5
	  R2=(M2-1)*0.5
	ELSE IF (E2.LT.20.0) THEN
	  M2=21+(E2-10.)/0.5
	  R2=10.+(M2-21.)*0.5
	ELSE IF (E2.LT.40.0) THEN
	  M2=41+(E2-20.)/1.
	  R2=20.+(M2-41)*1.
	ELSE IF (E2.LT.100.0) THEN
	  M2=61+(E2-40.)/3.
	  R2=40.+(M2-61)*3.
	ELSE IF (E2.LT.200.0) THEN
	  M2=81+(E2-100.)/5.
	  R2=100.+(M2-81)*5.
	ELSE IF (E2.LT.400.0) THEN
	  M2=101+(E2-200.)/10.
	  R2=200.+(M2-101)*10.
	ELSE IF (E2.LT.600.0) THEN
	  M2=121+(E2-400.)/10.
	  R2=400.+(M2-121)*10.
	ELSE IF (E2.LT.1000.0) THEN
	  M2=141+(E2-600.)/20.
	  R2=600.+(M2-141)*20.
	ELSE IF (E2.LT.2000.0) THEN
	  M2=161+(E2-1000.)/50.
	  R2=1000.+(M2-161)*50.
	ELSE
	  M2=181+(E2-2000.)/400.
	  R2=2000.+(M2-181)*400.
	ENDIF

   43 RETURN
	END

CG==========================
      FUNCTION SPERFC(DUMMY)
      IF (DUMMY-8. .lt. 0) GO TO 10
      IF (DUMMY-8. .eq. 0) GO TO 10
      IF (DUMMY-8. .gt. 0) GO TO 20
   10 SPERFC=(1.0606963+0.55643831*DUMMY)/(1.0619896+1.72 45609*DUMMY
     &   +DUMMY*DUMMY)
      GO TO 40
   20 SPERFC=0.56498823/(0.06651874+DUMMY)
   40 CONTINUE
      RETURN
      END

       SUBROUTINE NEUTD(FLAG,QQ,DD,RR,DENS)
       REAL UN, TAU1, TAU2, TAU3, LAMBDA1, LAMBDA2
       REAL LAMBDA3, DENSTOT, DENS, PI
       INTEGER FLAG

C  	FLAG = 1 -> H2O
C	FLAG = 2 -> CO2
C	FLAG = 3 -> CO

       PI = 2.*ASIN(1.0)
       UN = 1.0E+5
       TAU1 = 1.0E+6*DD*DD
       TAU2 = TAU1
       TAU3 = TAU1
       LAMBDA1 = TAU1*UN
       LAMBDA2 = TAU2*UN
       LAMBDA3 = TAU3*UN
       DENSTOT = QQ/(RR*RR*UN*4*PI)

       IF(FLAG.EQ.1) THEN
         DENS = 0.85*DENSTOT*EXP(-RR/LAMBDA1)
       ELSEIF(FLAG.EQ.2) THEN
         DENS = 0.08*DENSTOT*EXP(-RR/LAMBDA2)
       ELSE
         DENS = 0.07*DENSTOT*EXP(-RR/LAMBDA3)
       ENDIF

       RETURN
       END

C============================
	subroutine para2(x0,s,alpha,r)
c	common /angle/theta,rp,alpha		! for checking
        REAL x1,x0,s,rp,r,theta,alpha,rflank

	x1=x0*0.5
C theta, rp are angle and radius with respect to the focal point of parabola.
	call parasr(x1,s,rp)
c	call parasr(x0,s,r)
	theta=2.0*acos(sqrt(x1/rp))

	r=rp**2+x1**2+rp*x0*cos(theta)
	r=sqrt(r)

C alpha is the sza of (s,r) point measured from the center of the planet.
	alpha=asin(rp/r*sin(theta))
	alpha=alpha*180.0/3.1415926
	rflank=sqrt(2.0)*x0
	if (r.gt.rflank) alpha=180.0-alpha
c	write(50,*) alpha, r
	return
	end

CCCC   Subroutine to calculate the radius of each spatial grid
CCCC   calls at line ???.
      SUBROUTINE PARASR(X0,S,R)
C         X0:   subsolar altitude of the magnetic field line.
C         S:    arc length between the subsolar point and grid.
C         R:    radius of grid point.
      REAL X0,S,R,S1,R1,S2,R2
      INTEGER I
      IF (S .LE. (X0*0.3)) GOTO 404
C     first guess.
        R=X0+S
        S1=((R-X0)*R)**0.5+X0*ALOG(((R-X0)**0.5+R**0.5)/X0**0.5)
        R1=(S-X0*ALOG(((R-X0)**0.5+R**0.5)/X0**0.5))**2/R+X0
       DO 401 I=1,50
        S2=((R1-X0)*R1)**0.5+X0*ALOG(((R1-X0)**0.5+R1**0.5)/X0**0.5)
        IF (ABS((S2-S)/S) .LT. 1.E-4) GOTO 402
      IF (ABS((S2-S1)/S) .LE. 1.E-4) GOTO 402
        R2=(R1-R)*(S-S1)/(S2-S1)+R
        IF (R2 .LE. X0) R2=X0
        R=R1
        S1=S2
        R1=R2
 401   CONTINUE
 402   CONTINUE
       R=R1
       GOTO 405
 404   R=0.2492*S**2/X0+X0
 405   CONTINUE
      RETURN
      END
C==============================
	SUBROUTINE EXTRAP(N,X,Y,XP,YP,IEXP,IRESET)
C ......   Program to interpolate data at xp. X(N),Y(N) are N dimensional
C ......   data set. Given a x point at xp, the program returns yp for y
C ......   value.
C ......   IEXP=0 : linear interpolation.
C ......   IEXP=1 : exponential-exponential linear interpolation.
C ......   IRESET=0  : keep old J to start search.
C ......   IRESET=1  : reset J to 1.

	REAL X(N),Y(N)
	INTEGER JEXTRAP
	IF (IRESET.EQ.1) J=1
c	J = 1
	CALL FINDJ(N,X,XP,J,IFLAG)
	IF (IEXP.EQ.0) THEN
	  CALL LINEAR(N,X,Y,XP,J,IFLAG,YP)
	ELSE IF (IEXP.EQ.1) THEN
	  CALL EXPO(N,X,Y,XP,J,IFLAG,YP)
	ENDIF
	RETURN
	END

	SUBROUTINE LINEAR(N,X,Y,XP,J,IFLAG,YP)
	REAL X(N),Y(N)
	IF (IFLAG .EQ. 1) THEN
	  YP=Y(J)
	ELSE
 		IF (IFLAG .EQ. -1) THEN
		  JL=1
		  JR=2
		ELSE IF (IFLAG .EQ. -2) THEN
		  JL=N-1
		  JR=N
		ELSE IF (IFLAG .EQ. 0) THEN
		  JL=J
		  JR=J+1
		ELSE
		    PRINT *,'IFLAG not returned right'
		ENDIF
		  YP=(Y(JR)-Y(JL))*(XP-X(JL))/(X(JR)-X(JL))+Y(JL)
	ENDIF
	RETURN
	END

	SUBROUTINE EXPO(N,X,Y,XP,J,IFLAG,YP)
	REAL X(N),Y(N),XX1,XX2,YR,YL,XR,XL,YP,XPL
	INTEGER JL,JR
	IF (IFLAG .EQ. 1) THEN
	  YP=Y(J)
	ELSE
 		IF (IFLAG .EQ. -1) THEN
		  JL=1
		  JR=2
		ELSE IF (IFLAG .EQ. -2) THEN
		  JL=N-1
		  JR=N
		ELSE IF (IFLAG .EQ. 0) THEN
		  JL=J
		  JR=J+1
		ELSE
		    PRINT *,'IFLAG not returned right'
		ENDIF

		if ( (y(jr).gt.1.e-30) .and. (y(jl).gt.1.e-30)
     &		.and. (x(jr).gt.1.e-30) .and. (x(jl).gt.1.e-30) ) then
		  YR=ALOG(Y(JR))
	   	  YL=ALOG(Y(JL))
		  XR=ALOG(X(JR))
	 	  XL=ALOG(X(JL))
		  XPL=ALOG(XP)
		  YP=(YR-YL)*(XPL-XL)/(XR-XL)+YL
	 	  YP=EXP(YP)
		else
		  YP=0.0
		endif
	ENDIF

	RETURN
	END

	SUBROUTINE FINDJ(NPT,XX,X,J,IFLAG)
	PARAMETER (NMX=500)
	REAL XX(NPT)
	IF (X .LT. XX(1))	GOTO 1
      IF (X .GT. XX(NPT)) GOTO 2
C	J=1
 13	IF (X .EQ. XX(J))	GOTO 14
	IF ((X .GT. XX(J)) .AND. (X .LT. XX(J+1))) GOTO 15
	J=J+1
	IF (J .LE. NPT) GOTO 13
	PRINT *,' ## INFINITE LOOP IN FINDJ, X=',X
	STOP
 1	IFLAG=-1
	RETURN
 2	IFLAG=-2
	RETURN
 14	IFLAG=1
	RETURN
 15	IFLAG=0
	RETURN
	END
