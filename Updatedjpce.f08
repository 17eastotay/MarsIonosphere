program JPCE
!*******************************************************************************
!* Edited by Antonio Renzaglia and Joseph W. Wimmergren 10.02.19
!         This code is taken from Stephen's version on Jupiter and adapted
!     to work on Mars.
!*******************************************************************************
!* Edited by Stephen J. Houston 2.25.18
!*******************************************************************************
!	  This code is taken from the original tcpe.f code designed for TITAN
!     I am modifying it to work in Jupiter for the auroral region.
!     This was taken from the tpce96.f code worked by M. Richard and C. Wyllie
!	  Modification to Jupiter Sept. 2011
!     I have removed all the comments that are irrelevant to the Jupiter work
!     For more details on the code or on earlier versions please refer
!     to the Titan versions.

!     Diagnostcs have been updated to print out the values for the 38
!     neutrals and 74 ions. 10.7.07 M Richard

!     MODifIED CODE TO read COMPLETE DENSITY FILE
!       NOTE if NIONS IS CHANGED GO TO GAUSS JORDAN

!  THIS IS NED'S IMPROVED CODE - ORIGINAL VERSION.
!  22Nov2005 - modified code to work on kumhd and with telect.  IPR
!  4July2005 - wrap up comments on this code - used in 1996 model that produced
!  the 1998 Planet & Space Sci paper - Keller, Anicich, & Cravens.
!  25 May 97 - add 5 new neutral photoionization rates from new ephotib2.for
!  for a total of 23 photoionization rates. Electron ionization rates are still
!  read in (only 18 of them from old 1992 two-stream code BUT the five new
!  photoionization rates are complemented by assuming that e-imp rates are the
!  same as photoionization rates.
!  17 Oct 96 - set up diagnostic routines to look at prod/loss.
!  10 Oct 96 - Make overhaul in program structure to simplify in preparation
!  for adding more ions.  NGrid=229=>725-3005x10 km.

!     Neutral Code:
!     1=CO2,2=CO, 3=O, 4=N2, 5=H2, 6=H, 7=Ar, 8=He, 9=H2O, 10=NO, 11=O2, 12=N, 13=C


!     Ion Code:
!     1=CO2+, 2=O2+, 3=O+, 4=CO+, 5=HCO+, 6=H+, 7=N2+, 8=Ar+, 9=H2+
!     10=NO+, 11=H3+, 12=OH+, 13=ArH+, 14=N2H+, 15=HCO2+, 16=HOC+, 17=C+, 18=N+, 19=He+, 20=NH+


!     Input Files Needed:
!	1. Neutral atmosphere;
!      	2. Photoionization Rates from JEphot;
!     	3. Magsphere electron impact ionization rates from two stream code;
!   	4. Temperature
!     Units in this code are all cgs.
!*******************************************************************************

implicit real*8(a-h,o-z)

!*******************************************************************************
integer NIons,NNeut,NGrid, NReact
character*1 Blank
parameter(NIons=21,NNeut=13,NGrid=500,NReact=122,Blank=' ')
integer ILen,I,J,K,L,M,N,MaxIter,IDiag,IT,NDiagPoint,o
integer NDiagIons,ID,JD,NProdTerms,NLossTerms,JJ,KK
integer LOrder(100),KOrder(100),LTemp,KTemp,JIonMax
integer IDiagPoints(10),JDiagIons(100),DIMaxIndx(12)
integer IYr,IMon,IDay,IHr,IMin,I100th,COUNT
integer JKJJ, JKJ,DUM1,DUM2,DUM3,DUM4
real Tolerance,ToldX,TolF,DNTmp,DInit,DeltaZ,ZLower,ISec
real TE,TF,PTerm,LTerm,ErrF,ErrdX,ProdTemp,LossTemp,ProdPer,LossPer
real,dimension(NIons) :: Prod,Loss,TrueLoss,PNet,Beta,DelX,ERCoeff,ERecomb,TD
real,dimension(NGrid) :: Alt,Rad,DE,TEGr,alt_te,te_orig
real,dimension(NGrid) :: D_CO2,D_CO,D_O,D_N2,D_H2,D_H,D_Ar,D_He,D_H2O,D_NO,D_O2,D_N
real,dimension(NIons,NIons) :: Alpha
real,dimension(NGrid,NIons) :: DI,ProdPh
real,dimension(NGrid,NNeut) :: DN, DNREF
real PhIon(NGrid,23),EImpPr(NGrid,18),ProdOrder(100),LossOrder(100)
real R(NNeut,NIons,NIons)
real DIMax(12)
real MassDensity(36)
real FEXP, FEXPA, AMU, GT, chi, elect, Tvib
real BOLTZ, MASSC2H2, MASSC2H4, MASSC2H6, MASSC4H2, SH(NGrid,4)
integer JJ1, pe
real Ionflux
parameter(BOLTZ=1.38E-23,MASSC2H2=26.03728,MASSC2H4=28.05316,MASSC2H6=30.06904)
parameter(MASSC4H2= 50.05868,AMU=1.66E-27,GT=2.479e3,Ionflux=1)
logical DiagFlag,SkipToDiagFlag
character*50 OutFile,AtmosFil,PhIonFil,EImpPrFil,TElecFil
character*70 ParamFile,DiagFile,ReactFil
!     Comment out the vax vms way of obtaining the date & time - msfortran
!     uses a bunch of integers to get these data.
character*9 DateStr
character*8 TimeStr
character*13 NName(NNeut)
character*13 IName(NIons)
character*3 sza
character(len=100) filename
!****************************** Data Declaration *******************************
!********************************* Initialize **********************************
R=0.0;TD=0.0;DN=0.0;TeGr=0.0;ERCoeff=0.0;PhIon=0.0;Rad=0.0;EImpPr=0.0
Tolerance=0;Pterm=0.0;LTerm=0.0;ErrF=0.0;ErrdX=0.0;ProdTemp=0.0;LossTemp=0.0
Prod=0.0;Loss=0.0;TrueLoss=0.0;PNet=0.0;Beta=0.0;DelX=0.0;ERecomb=0.0
Alpha=0.0;DI=0.0;ProdPH=0.0;ProdOrder=0.0;LossOrder=0.0
!******************************** Main Program *********************************
! Get current date and time for this run to put into output data file.
Call Date_and_Time(DateStr,TimeStr)
read(DateStr(1:4),'(i4)') IYr
read(DateStr(5:6),'(i2)') IMon
read(DateStr(7:8),'(i2)') IDay
read(TimeStr(1:2),'(i2)') IHr
read(TimeStr(3:4),'(i2)') IMin
read(TimeStr(5:8),'(F4.2)') ISec
! Open a file of parameters containing input/output file names,
!    Newton-Raphson tolerances, etc.
!	Print*, 'Please enter file name of the parameter file.'
!	read(5,51) ParamFile
ILen = (Index(ParamFile,Blank))-1
!      Open(Unit=25,File=ParamFile(1:ILen),Status='OLD')
! OPEN(3,FILE='input/cephot/szangles.par',STATUS='UNKNOWN')
! read(3,*)CHI
CHI=0
open(unit=124, file='dnref.txt', status='unknown')
open(unit=125, file='dnrefplox.txt',status='unknown')

open(unit=123, FILE='neutralstest.txt', status='unknown')

write(*,*) 'SZA = ', CHI
! write(filename,"('common/configs/jpce_parameters',I0,'.cfg')") int(CHI)
! open(unit=25,file=filename,status='unknown')
		chi =0.0 !Uncomment and comment previous two lines if you want a fixed SZA
      if(chi .eq. 0)then
      	SZA = '00/'
      	Open(Unit=25,FILE='./configs/jpce_parameters00test.cfg',&
            Status='old')
      elseif(chi .eq. 60)then
      	SZA = '60/'
      	Open(Unit=25,FILE='./configs/jpce_parameters60.cfg',&
            Status='OLD')
      elseif(chi .eq. 80)then
      	SZA = '80/'
      	Open(Unit=25,FILE='./configs/jpce_parameters80.cfg',&
            Status='OLD')
      elseif(chi .eq. 90)then
      	SZA = '90/'
      	Open(Unit=25,FILE='./configs/jpce_parameters90.cfg',&
            Status='OLD')
      elseif(chi .eq. 100)then
      	SZA = 'ME/'
      	Open(Unit=25,FILE='./configs/jpce_parametersMonoE.cfg',&
           Status='OLD')
      elseif(chi .eq. 200)then
      	SZA = 'SE/'
      	Open(Unit=25,FILE='./configs/jpce_parametersSecEl.cfg',&
           Status='OLD')
      end if

! read names from param file and open input and output files.
! First read five comment lines at top of parameter file, then
! read the flag indicating whether or not diagnostic output is desired.
do I=1,11 !
	read(25,*)
end do
!* Atmosphere file. SJH - Well-mixed atmosphere of Jupiter
write(AtmosFil,"('./input/neutraldens_expanded_202108.txt')") !change here JW
open(unit=30,file=AtmosFil,status='old')
!* Reaction file
write(ReactFil,"('./input/marschemistry.dat')") !change here JW
open(unit=60,file=ReactFil,status='old')
!* Ion production from heavy ion precipitation (labeled as photoion here)
write(PhIonFil,"('./input/primprod_flare.dat')")  !int(CHI)
open(unit=35,file=PhIonFil,status='unknown')
!* Ion production from secondary electron impact
! write(EImpPrFil,"('./common/output/',I0,'/secprod.dat')") !int(CHI)
! open(unit=40,file=EImpPrFil,status='unknown')
!* Electron temp
 write(TElecFil,"('./input/electemp.dat')")
 open(unit=45,file=TElecFil,status='old') !Electron temp.


!* Output log file
 write(filename,"('./output/jpce.dat')") !int(CHI)
open(unit=50,file=filename)
!* Output atmosphere file of all densities used in the code
open(unit=57,file='jpce_density.dat')
!* Electron recombination cross-sections
open(unit=56,file='ERECOMB.DAT')
!* Lifetime file
write(filename,"('./output/Lifetime.dat')") !int(CHI)
open(unit=58,file=filename)
!* JPCE production file
write(filename,"('./output/jpce_prodin.dat')") !int(CHI)
open(unit=59,file=filename)
! write(filename,"('./common/output/',I0,'/jpce80.out')") int(CHI)
! open(unit=70,file=filename,status='old') !jpce 80 output file
!* JPCE ion densities
open(unit=71,file='jpceIonDensin.dat')
!* check electon temperature
open(unit=72,file='mpceETemp.dat')

! write to output file information on which input files were used.
write(50,6000)'! Neutral Atmosphere used = ',trim(AtmosFil),&
						 	'DATE = ',IMon,'-',IDay,'-',IYr
write(50,6001)'! Photoionization Rates Used = ',trim(PhIonFil),&
							'TIME = ',IHr,':',IMin,':',ISec
write(50,*)'! Electron Impact Ionization Rates Used = ',EImpPrFil
write(50,*)'! Temperatures used = ',trim(TElecFil)
6000 format(1x,A,A,3x,A,I2,A,I2,A,I4)
6001 format(1x,A,A,13x,A,I2,A,I2,A,F4.1)

! read in the electron TEMPERATURE profile
! Discard first 4 lines of comments and labels.  Rad is dummy.
 do i=1,Ngrid
 	read(45,*) alt_te(i),te_orig(i)
 end do

! read in the NEUTRAL atmosphere at each grid point.
read(30,*)
do  I=1,NGRID !Assuming a well mixed atmosphere
	read(30,*) Rad(I),D_N2(I),D_CO2(I),D_CO(I),D_O(I),D_He(I),&
	 D_NO(I),D_N(I),D_O2(I),D_Ar(I),D_H(I)		!JW changed Rad(i) to Alt(i)
end do
!Rad=Rad*1e5 !ARR commented out. Keep rad in km, not cm
! Use the H2 density to find the H2vib density by multiplying by exp{-21960/T(z)}  !JW aren't using H2vib?
! do I=1,NGrid
! 	Tvib= 1470				!commented this line, the line before, and the line after because we aren't using hvib?
! 	DN(I,9) = DN(I,1)*EXP(-(21960)/(2*TeGr(i)))  !(TeGr(I)*1))!Use *2 for photoelectrons
! end do
! Calculate the neutral densities for the other hydrocarbons
! For altitudes below the homopause (350 km for Maurellis atmos) use
! mixing ratios taken from: Perry et al. 1999 and Gladstone et al. 1996
! for higher altitudes use static atmosphere laws. Calculate first the scale height
! for each element and use the density at 350km as n_0
! Then n = n_0 * exp(-z/H)
do I=1,NGrid !!!   need H2 and H2O, scale to CO2 density
	! if(I .lt. 17) then
		D_H2(I) = D_He(I)        !H2 from Haider et al., 2011
		D_H2O(I) = D_NO(I)        !H2O from Haider et al., 2011
! 	elseif(I .ge. 17)then
! ! scale heights
! 		SH(I,1) = (BOLTZ* TeGr(I))/(MASSC2H2*AMU * GT/100)/1000
! 		SH(I,2) = (BOLTZ* TeGr(I))/(MASSC2H4*AMU * GT/100)/1000
! 		SH(I,3) = (BOLTZ* TeGr(I))/(MASSC2H6*AMU * GT/100)/1000
! 		SH(I,4) = (BOLTZ* TeGr(I))/(MASSC4H2*AMU * GT/100)/1000
! 	! Densities
! 		DN(I,5) = DN(16,5) * EXP(-(RAD(I)*1E-5-350)/SH(I,1))
! 		DN(I,6) = DN(16,6) * EXP(-(RAD(I)*1E-5-350)/SH(I,2))
! 		DN(I,7) = DN(16,7) * EXP(-(RAD(I)*1E-5-350)/SH(I,3))
! 		DN(I,8) = DN(16,8) * EXP(-(RAD(I)*1E-5-350)/SH(I,4))
! 	end if
end do

! read In Neutral Species Number Assigned and Name
read(60,*)
do DUM1=1, NNeut
	read(60,*) DUM2, NName(DUM1)
	write (123,*) DUM2, nName(DUM1) !changed by jw
end do

! read in Ion Species Name, ERCoeff and TD
read(60,*)
read(60,*)
do DUM1=1, NIons
	read(60,*) DUM2, IName(DUM1), ERCoeff(DUM1), TD(DUM1)
end do

! read in REACTION LIST from file in the form Neutral X + Ion Y --> Ion Z
read(60,*)
read(60,*)
read(60,*)
do DUM1=1, NReact
	read(60,*) DUM2, DUM3, DUM4, R(DUM2,DUM3,DUM4)
end do

!   Print Off Reaction Rates
open(37,File='./output/REACTION_LIST.DAT')
write(37,*) 'NEUTRAL LIST'
Do DUM1=1,NNEUT
	write(37,*) DUM1, trim(NName(DUM1))
end do
write(37,*)
write(37,*) 'ION LIST ERCoeff'
do DUM1=1,NIONS
	write(37,'(I3,A,A13,2ES10.2)') DUM1,' ',trim(IName(DUM1)),ERCoeff(DUM1),TD(DUM1)
end do
write(37,*)
write(37,*) 'REACTION LIST'
write(37,*) 'Neutral X + Ion Y --> Ion Z'
Do DUM1=1,NNEUT
	Do DUM2=1,NIONS
  	Do DUM3=1,NIONS
			if (R(DUM1,DUM2,DUM3) .NE. 0.0)then
				write(37,'(I5,I5,I5,ES10.3,A,A7,A,A13,A,A13)') DUM1,DUM2,DUM3,&
				R(DUM1,DUM2,DUM3),' ',trim(NName(DUM1)),' + ',trim(IName(DUM2)),' --> ',trim(IName(DUM3))
			end if
		end do
	end do
end do
!* Uncomment the next statements if you want to turn off any photoelec contrib.
!* Check to see if photoelectrons are being taken into consideration
!	  open(1, file='input/elect/jelect.par',STATUS='old')
!	  read(1,*) pe !pe = 0 -> no photoe; pe = 1 -> include photoe
!	  rewind(1)
pe=1 !pe=2 sets the primary production to ion production from ion precipitation

if(pe .eq. 0)then !set all primary production to zero if there are no photoelec considered
	open(2,file='./Desktop/Fortran/JPCE/primprod.dat',status='unknown') !input????? But listed as an output earlier in the file
	do i=1,NGrid
		read(2,*) k,(PhIon(i,j),j=1,15)
	  if(PhIon(i,j) .eq. 0.0)then
	  	PhIon(I,J) = 1e-20
	  endif
  enddo
elseif(pe .eq. 1)then
	! read in PHOTIONIZATION production rates at each grid point.
	! Discard first 3 lines of comments and labels.  Rad is dummy.
	!do I=1,1  !get rid of rows of ion names/other stuff
		read(35,*)
	!end do
	write(*,*) 'photoelectron input'
	do I=1,NGrid
		read(35,*)Alt(i),(PhIon(I,J),J=1,15)
	  if(PhIon(I,J) .eq. 0.0)then
	    PhIon(i,j) = 1e-20
	  endif
	end do
elseif(pe .eq. 2)then
! read the prim prod from ion precip
	read(35,*) !Skip the header line
	write(*,*) 'ion precip input'
	!skip 5 total
	read (35,*)
	read (35,*)
	read (35,*)
	read (35,*)
	do I=NGrid,1,-1 !Needs to be flipped - SJH
		read(35,*)Rad(i),(PhIon(I,J),J=1,8)
		! write(*,*)Rad(i),(PhIon(I,J),J=1,8)
    if(PhIon(I,J) .eq. 0.0)then
    	PhIon(i,j) = 1e-20
		endif
  end do
else
	write(*,*) 'wrong photoelectron criteria'
endif
! read in ELECTRON IMPACT-ion production rate at each grid point.
! Discard first lines of comments and labels.  Rad is dummy.
! read(40,*)
! do I=1,NGrid
! 	read(40,*)Rad(i),(EImpPr(I,J),J=1,8)
!   if(EImpPr(I,J) .eq. 0.0)then
!   	EImpPr(i,j) = 1e-20
!   endif
! write(*,*) NGRID, RAD
! end do
!PhIon=PhIon*1e6 !SJH - multiplying by 10^6 to make it a more realistic input !commented out by jw
! EImpPr=EImpPr*1e6

!!!convert neutral density altitudes to match primprod
!!!also want to convert electron temp altitudese to match primprod
!call extrap(500,rkm1,N2,rad(J),D_N2(J),1,1)

DO J=1,NGRID
	call extrap(500,RAD,D_CO2,Alt(J),DN(J,1),1,1)
	call extrap(500,RAD,D_CO,Alt(J),DN(J,2),1,1)
	call extrap(500,RAD,D_O,Alt(J),DN(J,3),1,1)
	call extrap(500,RAD,D_N2,Alt(J),DN(J,4),1,1)
	call extrap(500,RAD,D_H2,Alt(J),DN(J,5),1,1)
	call extrap(500,RAD,D_H,Alt(J),DN(J,6),1,1)
	call extrap(500,RAD,D_Ar,Alt(J),DN(J,7),1,1)
	call extrap(500,RAD,D_He,Alt(J),DN(J,8),1,1)
	call extrap(500,RAD,D_H2O,Alt(J),DN(J,9),1,1)
	call extrap(500,RAD,D_NO,Alt(J),DN(J,10),1,1)
	call extrap(500,RAD,D_O2,Alt(J),DN(J,11),1,1)
	call extrap(500,RAD,D_N,Alt(J),DN(J,12),1,1)

	call extrap(500,RAD,Te_orig,Alt(J),TeGr(J),1,1)
ENDDO


!testing the neutral densities

write(57,*) '# Alt [km]    CO2        CO         O        N2      H2      ',& !JW Changed these
		'H      Ar      He      H2O      NO      O2      N  [cm^-3]'

do  I=NGrid,1,-1 !Assuming a well mixed atmosphere
	write(123,*) Rad(I), DN(I,1),DN(I,2),DN(I,3),DN(I,4),DN(I,5),DN(I,6),DN(I,7),&
		DN(I,8),DN(I,9),DN(I,10),DN(I,11),DN(I,12)
end do

do I = 1,NGRID
  write(57,7115) Alt(i), (DN(I,J),J=1,12)
	do o=1,12
		DNREF(I,o)=DN(I,o)
	end do
7115 format(1x,10(1PE10.3))
end do

! read the parameters needed for the Newton-Raphson routine.
! MaxTrial = the maximum number of iterations allowed even if
!              the procedure does not converge.
! TolF and ToldX are the tolerances allowed in the sum of the functions
!    or sum of the changes in the roots in order to verify that convergence
!    has occurred.
! do i=1, 10
! 	read (25,*)    !JW put here to get rid of first 11 lines to get to tolerance
! end do

read(25,*) Tolerance
read(25,*) MaxIter
TolF= Tolerance
ToldX= Tolerance

! read initial guesses for the ion densities DI.
read(25,*) DInit
!* commented to be read in from photoe guesses
do J=1,NIons
	DI(1,J) = DInit
end do
write(50,46)Tolerance,DInit
46 format(1x,'! Tolerance = ',1PE8.1,' Initial Densities=',1PE8.1)

! Set up the altitude grid. SJH - No longer use this altitude grid JW: uncommented this
! ZLower = 100.0
! DeltaZ = 3.0
! do I=1,NGrid
! 	Alt(I) = ZLower + DeltaZ*(I-1)
!end do
!Alt=Rad/1e5 !commented out by jw
! Set up massive if statement to cancel out real calcs to jump to
! diagnostics for testing them.
SkipToDiagFlag = .FALSE.
if(.NOT.(SkipToDIagFlag)) then
! Start the MAIN LOOP computing the densities at each gridpoint.
	! do k = 1,7
  ! 	read(70,*)
  ! end do
  do I=1,NGrid !Loop through every altitude
! Set up the guesses for the ion densities at this grid point.  If this is
! the lowest grid point, use the values input in the parameter file.  If
! this is any other altitude point - use the results from previous point
! as initial guesses for this grid point.
!  		read(70,*) rad(i), elect, (DI(I,J),J = 1,NIons)
		if (I.NE.1) then
!	  		read(70,*) rad(i), elect, (DI(I,J),J = 1,NIons)
!	  		write(71,*)rad(i), elect, (DI(I,J),J = 1,NIons)
! nom - mod to read previous output files
			do 70 J=1,NIons
					DI(I,J) = DI((I-1),J)
					70 continue
		end if
! Extract the electron temperature at this grid point.
		TE = TEGr(I)
! Compute the ion production terms which come from photon and electron
!    ionization.  Normally this is done by simply adding the rates from
!    the input files of this data.
! Since all the ions do NOT have sources from photo/magsphere-ionization
! there is not an easy way to set up a vector for these rates
! and simply read in the data directly to the rate vector.
! Therefore this set of correlation relations is needed. 1,3,4,7,
		ProdPh(I,1) = PhIon(I,3)  !+ 1*EImpPr(I,2)    !CO2+ prod
		ProdPh(I,2) = PhIon(I,9)  !+ 1*EImpPr(I,3)    !O2+ prod
		ProdPh(I,3) = PhIon(I,5)  !+ 1*EImpPr(I,1)    !O+ prod
		ProdPh(I,4) = PhIon(I,4)  !+ 1*EImpPr(I,1)    !CO+ prod
		ProdPh(I,6) = PhIon(I,11) !+ 1*EImpPr(I,1)    !H+ prod
		ProdPh(I,7) = PhIon(I,2)  !+ 1*EImpPr(I,1)    !N2+ prod
		ProdPh(I,8) = PhIon(I,10) !+ 1*EImpPr(I,1)    !Ar+ prod
		ProdPh(I,10) = PhIon(I,8) !+ 1*EImpPr(I,1)    !NO+ prod
		ProdPh(I,17) = PhIon(I,6) !+ 1*EImpPr(I,1)    !C+ prod
		ProdPh(I,18) = PhIon(I,1) !+ 1*EImpPr(I,1)    !N+ prod
		ProdPh(I,19) = PhIon(I,7) !+ 1*EImpPr(I,1)    !He+ prod


! Multiply the ProdPh by the total # of ions
		do k=1,21
			ProdPh(i,k) = Ionflux*ProdPh(i,k)
		enddo
		write(59,*) i,ProdPh(I,1),ProdPh(I,3),ProdPh(I,4),&
			ProdPh(I,7)
! Compute the rate constants which depend on the electron temp.
! This if for recomb react rates of the form k= ERCoeff(300/Te)**TD
! TF = Temp factor for elec recomb rxs
		write(72,*) Alt(i), TE
		TF = 3.E2/TE
! TD and ERCoeff are obtained from the chem rates input file
		do J=1,NIons
			ERecomb(J) = ERCoeff(J)*TF**TD(J)
			if (J == 21 .and. TE > 800) then
				ERecomb(J) = 4.73E-05*TF**(-0.74)
			endif
		end do
		if(i.eq.39)then
			do JJ1=1,NIons
      	write(56,*)'ERECOMB(',JJ1,') =',ERecomb(JJ1)
			end do
		end if
! Start Newton-Raphson algorithm iterations to find the
! photochemical solutions for DI at the given gridpoint I.
		do 900 IT=1,MaxIter
! Zero the vector and matrix used in calls to Gauss-Jordan code as a part of
!     the Newton-Raphson algorithm.  Just for safety.
			Beta=0.0;Alpha=0.0
			DE(I) = 0.0 ! Index I is the altitude (NGrids)
			do J=1,NIons !find electron density at this altitude & iter
 				DE(I) = DE(I) + DI(I,J)
			end do
			do J=1,NIons
	    	Prod(J) = ProdPh(I,J)    !start production with photo and e-imp prd
	    	Loss(J) = ERecomb(J)*DE(I)  !start loss with e-recomb
	    	do L=1,NNeut !sum all prod/loss terms for the jth ion at this alt
	    		do K=1,NIons
						PTerm = DI(I,K)*DN(I,L)*r(L,K,J)  !r(lkj)=> n(L)+i(K)-->i(J)
						LTerm = DN(I,L)*r(L,J,K)          !r(ljk)=> n(L)+i(J)-->i(K)
						Prod(J) = Prod(J) + PTerm
						Loss(J) = Loss(J) + LTerm
					end do ! K=1,NIons
				end do ! L=1,NNeut
	    	TrueLoss(J) = Loss(J)*DI(I,J)  !Loss is useful for finding lifetime
				! write(*,*) DI(I,J),Loss(J),r(2,J,2)
	    	PNet(J) = Prod(J)-TrueLoss(J) !    but we need true loss in pnet
	    	Beta(J) = -PNet(J)
			end do ! J=1,NIons
! Check to see if we have hit the roots by checking to see if the sum of all
! the functions (-Beta's) just computed is less than the given convergence
! tolerance.  If the sum is less, then we have the roots and can go to the
! next altitude point.
			ErrF = 0.0
			do J=1,NIons
				if(IT .gt. 400)then
					write(*,*) 'nion=', j, 'beta= ', Beta(j)
				end if
				ErrF = ErrF + ABS(Beta(J))
				! write(*,*) j, Errf, beta(J)
			end do
			write(*,*) 'I=',I,'  Alt=',Rad(I),' Iteration=',IT,' ErrF=',ErrF
			if (ErrF .LE. TolF) then
      	goto 999
			end if
! Compute the alpha matrix for the NR algorithm.  This is a matrix of
! derivatives of each function with respect to each root.
! Alpha(j,k)=d(PNet)j/(dx)k - There are standard expressions for diagonal
!  and off-diagonal elements.
			do J=1,NIons
				do K=1,NIons
					Alpha(J,K) = -ERecomb(J)*DI(I,J)
					do L=1,NNeut
						Alpha(J,K) = Alpha(J,K) + DN(I,L)*r(L,K,J)
					end do !L=1,NNeut
				end do !K=1,NIons through all elements of the matrix for this ion
				Alpha(J,J) = -(Loss(J)) - (DI(I,J)*ERecomb(J))  !reset diagonal ele
			end do !J=1,NIons through all the ions

! Now we have the values for the alpha matrix and the beta vector.  These are
! input to the Gauss-Jordan matrix inversion routine, which returns the inverse
! of alpha in the same memory location as the input matrix, and returns the
! solution (i.e. dX) of the matrix equation (Alpha)*(dX)=Beta in the same
! location as the vector beta was input.
			Call GaussJ(Alpha,NIons,NIons,Beta,1,1)
open(17, file="conv.txt", status="unknown", position="append", action="write") !JW
! With this new information from GaussJ we use the NR algorithm to set up
! new guesses for the roots.  The algorithm states: X(n+1)=X(n) + dX
			do J=1,NIons
				DelX(J) = Beta(J)
				DI(I,J) = DI(I,J) + DelX(J)
				if (J.eq.2)then
					if (I.eq.3)then
					write(17,*) DI(I,J)
				end if
			end if

! nom Modification: There are species that end up with SMALL negative numbers
! I will force them to be positive numbers to make sure we have no calculation
! issues
				if(DI(I,J) .LT. 0.0)then
					DI(I,J)= 0.0
				end if
			end do
! Now check to see if we hit the roots by summing the values of dX.  This
! is a second possible criterion for convergence.  As we approach the roots,
! the values of dX for each variable approach zero.  If we have hit the roots
! we go to the next altitude point, otherwise we go back and set up and run
! another NR iteration with the new guesses for ion densities.  If we have
! hit the roots we compute the electron density corresponding to this set
! of ion densities, so it will be consistent with this set of ion densities,
! not with the previous iterations set of densities used to calculate DE above.
			ErrdX = 0.0
			do J=1,NIons
				ErrdX = ErrdX + ABS(DelX(J))
			end do

			if(ErrdX.LE.ToldX) then
				DE(I) = 0.0
				do J=1,NIons
					DE(I) = DE(I) + DI(I,J)
				end do




				goto 999 !break out of the iteration loop and go to next alt
			end if  !whether or not have found roots in dX
		900 Continue    ! with NR) iterations until converged for this alt
		999 if (I .eq. 1) then !XRDB change this to our neutrals?
      write(58,*) '    Alt            CO2               CO                O                N2                &
			              H2               H              Ar             He                  H2O'
    end if
		write(58,*) Alt(I),1.0/Loss(1),1.0/Loss(2),1.0/Loss(3), 1.0/Loss(4), &
		1.0/Loss(5), 1.0/Loss(6), 1.0/Loss(7), 1.0/Loss(8), 1.0/Loss(9)
  !  write(58,1023) Alt(I),1.0/Loss(1),1.0/Loss(2),1.0/Loss(3), 1.0/Loss(4), &
		!1.0/Loss(5), 1.0/Loss(6), 1.0/Loss(7), 1.0/Loss(8), 1.0/Loss(9)
		!1023 format(F8.2,4(2x,ES10.4))




		write (125,*) DNREF(i,1),DNREF(i,2),DNREF(i,3),DNREF(i,4),DNREF(i,5),DNREF(i,6),DNREF(i,7),DNREF(i,8), DNREF(i,9)


	if (DN(I,8) .NE. DNREF(I,8)) then !FLAG
		write(124,*) "Problem at altitude", I, DN(I,8), DNREF(I,8)
	else
		write(124,*) I, "NO ERROR"

	end if


	end do !End I=1,NGrid Loop through every altitude
! Finished with the computation of the densities as a function of altitude.
! Now print the results to the output file in the specified format.
! First write a header line.
  write(50,*)'!'
  write(50,1024) 'Z(km)','e-',(trim(INAME(DUM1)),DUM1=1,NIons)
	1024 format(a5,a8,2x,27(a10,2x))

  do 524 I=1,NGrid!was NGrid
    write(50,1025) Alt(I),DE(I),(DI(I,J),J=1,NIons)
    write(101,*) DI(I,1), ALT(I)
    write(102,*) DI(I,2), ALT(I)
    write(103,*) DI(I,3), ALT(I)
    write(104,*) DI(I,4), ALT(I)
    write(105,*) DI(I,5), ALT(I)
    write(106,*) DI(I,6), ALT(I)
    write(107,*) DI(I,7), ALT(I)
    write(108,*) DI(I,8), ALT(I)
    write(109,*) DI(I,9), ALT(I)
    write(110,*) DI(I,10), ALT(I)
    write(111,*) DI(I,11), ALT(I)
    write(112,*) DI(I,12), ALT(I)
    write(113,*) DI(I,13), ALT(I)
    write(114,*) DI(I,14), ALT(I)
    write(115,*) DI(I,15), ALT(I)
    write(116,*) DI(I,16), ALT(I)
    write(117,*) DI(I,17), ALT(I)
    write(118,*) DI(I,18), ALT(I)
    write(119,*) DI(I,19), ALT(I)
    write(120,*) DI(I,20), ALT(I)
    write(121,*) DI(I,21), ALT(I)
















		1025 format(F5.0,'	',28(1PE10.3,2x))
	524 Continue
 end if  !skip calculations to check diagnostic section

! Diagnostic section - get various diagnostics for desired altitudes and ions.
	read(25,*) DiagFlag,NDiagPoint,NDiagIons
	if(DiagFlag) then
		read(25,*) (IDiagPoints(I),I=1,NDiagPoint) !get indices of alt points
	  if (NDiagIons .EQ. 1000) then   !want to get all the ions w/ diagnostics
	  	do J=1,NIons
				JDiagIons(J) = J     !set the ions desired to each ion species
			end do
	  	read(25,*)      !read the next line but ignore it
	  else
	  	read(25,*) (JDiagIons(J),J=1,NDiagIons)   !get indices of ion species
	  end if     !desire to get diagnostics for every ion
		write(DiagFile,"('./jpce/jpce',I0,'.diag')") int(CHI)
		open(unit=55,file=DiagFile)
		do I=1,NDiagPoint !Loop through specific altitudes
	  	ID = IDiagPoints(I)
	  	TE = TeGr(ID)
	    write(55,7100) IDay,IMon,IYr,IHr,IMin,ISec,Alt(ID)
7100 format(1x,'Diagnostics for Mars PCE on ',I2,'-',I2,'-',I4,' at ',I2,&
		 ':',I2,':',F4.1,/,5x,'at an altitude of ',F7.1,&
		 ' km where the neutral densities are:',/,3x,'CO2',8x,'CO',7x,'O',8x,'N2',8x,'H2',8x,'H',8x,'Ar',8x,'He',8x,'H20',7x,&
			'NO',8x,'O2',8x,'N') !JOE WIMMERGREN need to these to the 12 Neutral's for Mars
	  	write(55,7110) (DN(ID,J),J=1,12) !AR changed from 1 to 12, and included the neutrals as well as 8x spaces between
7110 format(1x,4(1PE10.3))

	  	write(55,*)'At this altitude the electron and ion densities are:'
	  	write(55,7120)
7120 format(3X,'e-',8x,'1-CO2+',5X,'2-O2+',6X,'3-O+',4X,'4-CO+',4X,'5-HCO+',4X,&
		 '6-H+',5X,'7-N2+',6X,'8-Ar+',5X,'9-H2+',5X,'10-NO+',5x,'11-H3+')
	  	write(55,7022) DE(ID),(DI(ID,J),J=1,11) !JW: Don't know why these are split up in this way
7022 format(1X,12(1PE10.3))
	  	write(55,7024)
7024 format(3X,'12-OH+',4X,'13-ArH+',3X,'14-N2H+',4X,'15-HCO2+',3X,&
		 '16-HOC+',3X,'17-C+',4X,'18-N+',4X,'19-He+',3X,'20-NH+',3X,'21-H3O+')
	  	write(55,7022) (DI(ID,J),J=12,21)
	  	! write(55,7026) !commented this and the next line out as we only have 16 ions.
! 7026 format(3X,'24-CnHm+',2X,'25-c-C3H6+',1x,'26-HC3H+',3X,'27-C2H+')
! 	  	write(55,7022) (DI(ID,J),J=24,27) !jw commented this out and only had 16 ions
	  	write(55,*)'The electron temperature = ',Te,' ==> TFact = ',(300.0/Te),&
			' ==> Sqrt(TFact) = ',(SQRT((300.0/Te)))

! Sort out the top twelve species at this altitude.
	  do K=1,12
	    DIMaxIndx(K) = 0
      DIMax(K) =  0.0
		end do
	  do J=1,NIons    !find the largest density first
	    if (DI(ID,J).GE.DIMax(1)) then
	      DIMax(1) = DI(ID,J)
	      DIMaxIndx(1) = J
	    end if
		end do
	  do K=2,12
	    do J=1,NIons
	      if((DI(ID,J).GE.DIMax(K)).AND.(DI(ID,J).LT.DIMax(K-1)))then
					DIMax(K) = DI(ID,J)
					DIMaxIndx(K) = J
				end if
			end do
		end do
	  write(55,7070)
7070 format(' The twelve most densely populated species at this altitude are:')
	  do K=1,12
			write(55,7073) DIMaxIndx(K),trim(IName(DIMaxIndx(K))),DIMax(K)
		end do
7073 format(1x,'Ion numbered ',I2,'(',A,')',' with density of ',1PE10.3)
	  if (NDiagIons.EQ.1000) then
	    JIonMax = NIons
	  else
	    JIonMax = NDiagIons
	  end if
	  do J=1,JIonMax
	    JD = JDiagIons(J)    !pick out the next desired ion species
	    write(55,7101) JD,trim(IName(JD))
7101 format(1x,/,1x,'The production terms for ion species numbered ',I3,'(',A,&
		 ')','  are:')
		 	write(55,7102)ProdPh(ID,JD)
7102 format(' Production from photon and electron impact = ',1PE10.3,' cm-3/s')
	    Prod(JD) = ProdPh(ID,JD)   !start prod for this ion at photoprod
			LOrder=0;KOrder=0;ProdOrder=0.0 !zero vector storing prod terms for sort
	    ! do JJ=1,100     !zero out vector storing prod terms for sort
	    !   LOrder(JJ) = 0
	    !   KOrder(JJ) = 0
			! 	ProdOrder(JJ) = 0.0
			! end do
	    NProdTerms = 0    !count number of non-zero production terms
	    do L=1,NNeut !loop through all possible ion-neutral rxs - get prod
	      do K=1,NIons
					PTerm = DI(ID,K)*DN(ID,L)*r(L,K,JD)
					if (PTerm.NE.0.0) then
		 				NProdTerms = NProdTerms + 1
					 	ProdOrder(NProdTerms) = PTerm
					 	LOrder(NProdTerms) = L   !track the indices that go with
					 	KOrder(NProdTerms) = K   !   this prod term
					 	Prod(JD) = Prod(JD) + PTerm   !sum all prod to get total
					end if !have a non-zero production with this reaction
				end do !End K=1,NIons
			end do !End L=1,NNeut
!loop through all possible ion-neutral rxs - get prod
	    do JJ=1,NProdTerms    !sort prod terms from largest to small
	      do KK=(JJ+1),NProdTerms !sorting of production terms by magnitude
					if(ProdOrder(JJ).LE.ProdOrder(KK)) then
		  			ProdTemp = ProdOrder(JJ)
		  			ProdOrder(JJ) = ProdOrder(KK)
		  			ProdOrder(KK) =ProdTemp
		  			LTemp = LOrder(JJ)
		  			LOrder(JJ) = LOrder(KK)
		  			LOrder(KK) = LTemp
		  			KTemp = KOrder(JJ)
		  			KOrder(JJ) = KOrder(KK)
		  			KOrder(KK) = KTemp
				 	end if    !Need to switch this pair
			 	end do !End KK=(JJ+1),NProdTerms
		 	end do !End JJ=1,NProdTerms
			do JJ=1,NProdTerms
				ProdPer = (ProdOrder(JJ)/Prod(JD))*100.0
				write(55,7160) JD,trim(IName(JD)),LOrder(JJ),trim(NName(LOrder(JJ))),KOrder(JJ),&
					trim(IName(KOrder(JJ))),ProdOrder(JJ),ProdPer
			end do !End JJ=1,NProdTerms
7160 format(1x,'Production of ion ',I2,'(',A,')',' from neutral ',I2,'(',A,')',&
		 ' reacting with ion ',I2,'(',A,') = ',1PE11.3,' cm-3/s','=',0PF8.3,'%')
		 !The 0P MUST be w/ F or it is multiplied by ten
		 	write(55,7162)JD,Prod(JD)
7162 format(1x,'The TOTAL PRODUCTION for ion ',I3,' = ',1PE10.3,' cm-3/s')
	    write(55,7106) JD,trim(IName(JD))
7106 format(1x,/,1x,'The modified loss terms for ion species numbered ',I3,&
		 '(',A,')','  are:')
	    ERecomb(JD) = ERCoeff(JD)*((300.0/Te)**TD(JD))
	    write(55,7108)(ERecomb(JD)*DE(ID))
7108 format(1x,'ModLoss from electron recombination is:',1PE10.3,' /s')
	    Loss(JD) = ERecomb(JD)*DE(ID)   !start loss with e-recomb
			LOrder=0;KOrder=0;ProdOrder=0.0 !zero vector storing prod terms for sort
	    ! do JJ=1,100     !zero out vector storing loss terms for sort
	    !   LOrder(JJ) = 0
	    !   KOrder(JJ) = 0
			! 	LossOrder(JJ) = 0.0
			! end do
	    NLossTerms = 0    !count number of non-zero loss terms
	    do L=1,NNeut !loop through all possible ion-neutral rxs - get losses
	      do K=1,NIons
					LTerm = DN(ID,L)*r(L,JD,K)
					if (LTerm.NE.0.0) then
		 				NLossTerms = NLossTerms + 1
		 				LossOrder(NLossTerms) = LTerm
		 				LOrder(NLossTerms) = L   !track the indices that go with
		 				KOrder(NLossTerms) = K   !   this loss term
		 				Loss(JD) = Loss(JD) + LTerm   !sum all losses to get total
					end if    !have a non-zero loss with this reaction
				end do !End K=1,NIons
			end do !End L=1,NNeut
			do JJ=1,NLossTerms    !sort loss terms from largest to small
				do KK=(JJ+1),NLossTerms !sorting of loss terms by magnitude
					if(LossOrder(JJ).LE.LossOrder(KK)) then
						LossTemp = LossOrder(JJ)
						LossOrder(JJ) = LossOrder(KK)
						LossOrder(KK) =LossTemp
						LTemp = LOrder(JJ)
						LOrder(JJ) = LOrder(KK)
						LOrder(KK) = LTemp
						KTemp = KOrder(JJ)
						KOrder(JJ) = KOrder(KK)
						KOrder(KK) = KTemp
					end if    !Need to switch this pair
				end do !End KK=(JJ+1),NProdTerms
			end do !End JJ=1,NProdTerms
			do JJ=1,NLossTerms
				LossPer = (LossOrder(JJ)/Loss(JD))*100.0
				write(55,7175) JD,trim(IName(JD)),LOrder(JJ),trim(NName(LOrder(JJ))),KOrder(JJ),&
			 		trim(IName(KOrder(JJ))),LossOrder(JJ),LossPer
			end do
7175  format(1x,'Loss of ion ',I3,'(',A,')',' from neutral ',I3,'(',A,') ',&
		  'reacting to produce ion ',I3,'(',A,')',' = ',1PE10.3,'/s','=',0PF8.3,'%')
	    write(55,7172)JD,Loss(JD)
7172  format(1x,'The TOTAL modified loss for ion ',I3,' = ',1PE10.3,' /s')
	   	write(55,7180)JD,(1.0/Loss(JD)),(Loss(JD)*DI(ID,JD))
7180  format(1x,'The lifetime for ion ',I3,' = ',1PE10.3,'s  and',&
      ' the actual loss rate = ',1PE10.3,' cm-3/s',/,1x,/)
	  end do !End J=1,JIonMax loop through all desired ion species

! Compute density in each mass bin at this altitude and output by bins.
!	  MassDensity(1) = DI(ID,11)             ! 1 amu
!	  MassDensity(2) = DI(ID,10)             ! 2 amu
!	  MassDensity(3) = DI(ID,9)              ! 3 amu
!	  MassDensity(4) = DI(ID,8)              !12 amu
!	  MassDensity(5) = DI(ID,7)              !13 amu
!	  MassDensity(6) = DI(ID,6)+DI(ID,2)     !14 amu
!	  MassDensity(7) = DI(ID,5)+DI(ID,12)    !15 amu
!	  MassDensity(8) = DI(ID,4)              !16 amu
!	  MassDensity(9) = DI(ID,3)              !17 amu
!	  MassDensity(10) = DI(ID,27)            !18 amu
!	  MassDensity(11) = DI(ID,26)            !19 amu
!	  MassDensity(12) = DI(ID,22)            !25 amu
!	  MassDensity(13) = DI(ID,14)+DI(ID,21)  !26 amu
!	  MassDensity(14) = DI(ID,15)+DI(ID,20)  !27 amu
!	  MassDensity(15) = DI(ID,1)+DI(ID,16)+DI(ID,19)   !28amu
!	  MassDensity(16) = DI(ID,13)+DI(ID,18)  !29 amu
!	  MassDensity(17) = DI(ID,17)            !30 amu
!	MassDensity(18) = DI(ID,28)            !37 amu
!	MassDensity(19) = DI(ID,29)+DI(ID,50)  !38 amu
!	MassDensity(20) = DI(ID,23)+DI(ID,25)  !39 amu
!	MassDensity(21) = DI(ID,30)            !40 amu
!	MassDensity(22) = DI(ID,24)            !41 amu
!	MassDensity(23) = DI(ID,31)            !42 amu
!	MassDensity(24) = DI(ID,32)            !43 amu
!	MassDensity(25) = DI(ID,46)            !50 amu
!	MassDensity(26) = DI(ID,40)+DI(ID,47)  !51 amu
!	MassDensity(27) = DI(ID,41)            !52 amu
!	MassDensity(28) = DI(ID,33)            !53 amu
!	MassDensity(29) = DI(ID,34)            !55 amu
!	MassDensity(30) = DI(ID,35)            !63 amu
!	MassDensity(31) = DI(ID,36)            !65 amu
!	MassDensity(32) = DI(ID,48)            !67 amu
!	MassDensity(33) = DI(ID,49)            !69 amu
!	MassDensity(34) = DI(ID,51)            !77 amu
!	MassDensity(35) = DI(ID,37)+DI(ID,45)  !79 amu
!	MassDensity(36) = DI(ID,38)            !91 amu

!	write(55,7520) Alt(ID)
!7520    format(1x,'Density per amu bin for altitude of: ',F5.0)
!	write(55,7521)
!7521    format(1x,2x,'1amu',6x,'2amu',6x,'3amu',6x,'12amu',5x,'13amu',
!     &5x,'14amu',5x,'15amu',5x,'16amu',5x,'17amu',5x,'18amu',5x,
!     &'19amu',5x ,'25amu')
!	write(55,7522) (MassDensity(J),J=1,12)
!7522    format(1x,12(1PE10.3))
!	write(55,7525)
!7525    format(1x,2x,'26amu',5x,'27amu',5x,'28amu',5x,'29amu',5x,
!     &'30amu',5x,'37amu',5x,'38amu',5x,'39amu',5x,'40amu',5x,'41amu',
!     &5x,'42amu',5x,'43amu')
!	write(55,7522) (MassDensity(J),J=13,24)
!	write(55,7529)
!7529    format(1x,2x,'50amu',5x,'51amu',5x,'52amu',5x,'53amu',5x,
!     &'55amu',5x,'63amu',5x,'65amu',5x,'67amu',5x,'69amu',5x,'77amu',
!     &5x,'79amu',5x,'91amu')
!	write(55,7522) (MassDensity(J),J=25,36)

! write a couple of lines to separate this alt from next one.
		write(55,7531)
7531 format(1x,/,1x,/,1x)
	end do !End I=1,NDiagPoint loop through all desired altitude points
end if !diagnostics desired


! do i=1, NGRID
! 	write (125,*) DNREF(I,1),DNREF(I,2),DNREF(I,3),DNREF(I,4),DNREF(I,5),DNREF(I,6),DNREF(I,7),DNREF(I,8), DNREF(I,9)
! end do

close(unit=123)
close(unit=124)
close(unit=125)
close(unit=11)
Close(Unit=25)
Close(Unit=30)
Close(Unit=35)
Close(Unit=40)
! Close(Unit=45)
Close(Unit=50)
Close(Unit=55)

51 format(A)
end program


!...5...10...15...20...25...30...35...40...45...50...55...60...65...70...75...80
SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)

!* From Press, et. al. - pp 28-29.
!* Linear equation solution by Gauss-Jordan elimination to solve equation
!* of the form: A . x = b.  A is an input matrix of N x N elements, stored
!* in an array of physical dimentions NP by NP.  B is an input matrix of N
!* x M containing the M right-hand side vectors, stored in an array of
!* physical dimensions NP by MP.  On output, A is replaced by its matrix
!* inverse, and B is replaced by the corresponding set of solution vectors (X).
parameter (NMAX = 27) !have to change?
DIMENSION A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
!* The integer arrays IPIV,INDXR, and INDXC are used for bookeeping on the
!* pivoting. NMAX should be as large as the largest anticipated value of N.
do J=1,N
	IPIV(J) = 0
end do
do I=1,N
	BIG = 0.0
	do J=1,N
		if(IPIV(J).NE.1) then
		  do K=1,N
				if (IPIV(K).EQ.0) then
					if (ABS(A(J,K)).GE.BIG) then
						BIG=ABS(A(J,K))
						IROW=J
						ICOL=K
					end if
				else if (IPIV(K).GT.1) then
					STOP 'Singular matrix at 12'
				end if
			end do
		end if
	end do
	IPIV(ICOL) = IPIV(ICOL) + 1
!* We now have the pivot element, so we interchange rows, if needed, to put
!* the pivot element on the diagonal.  The columns are not physically
!* interchanged, only relabled: INDX(I), the column of the Ith pivot
!* element, is the Ith column that is reduced, while INDXR(I) is the row in
!* which that pivot element was originially located.  If INDXR(I) is not
!* equal to INDXC(I) there is an implied column interchange.  With this
!* form of bookeeping, the solution B's will end up in the correct order,
!* and the inverse matrix will be scrambled by columns.
	if (IROW.NE.ICOL) then
		do L=1,N
		  DUM=A(IROW,L)
		  A(IROW,L)=A(ICOL,L)
		  A(ICOL,L)=DUM
		end do
		do L=1,M
		  DUM=B(IROW,L)
		  B(IROW,L)=B(ICOL,L)
		  B(ICOL,L)=DUM
		end do
	end if
	INDXR(I)=IROW
	INDXC(I)=ICOL
	if (A(ICOL,ICOL).EQ.0.) STOP 'Singular matrix at 15'
	PIVINV=1./A(ICOL,ICOL)
	A(ICOL,ICOL)=1.
	do L=1,N
		A(ICOL,L)=A(ICOL,L)*PIVINV
	end do
	do L=1,M
		B(ICOL,L)=B(ICOL,L)*PIVINV
	end do
	do LL=1,N
		if (LL.NE.ICOL) then
			DUM=A(LL,ICOL)
			A(LL,ICOL)=0.
			do L=1,N
				A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
			end do
		  do L=1,M
				B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
			end do
		end if
	end do
end do
do L=N,1,-1
	if (INDXR(L).NE.INDXC(L)) then
		do K=1,N
			DUM=A(K,INDXR(L))
		  A(K,INDXR(L))=A(K,INDXC(L))
		  A(K,INDXC(L))=DUM
		end do
	end if
end do
return
end subroutine

!==============================
	SUBROUTINE EXTRAP(N,X,Y,XP,YP,IEXP,IRESET)
! ......   Program to interpolate data at xp. X(N),Y(N) are N dimensional
! ......   data set. Given a x point at xp, the program returns yp for y
! ......   value.
! ......   IEXP=0 : linear interpolation.
! ......   IEXP=1 : exponential-exponential linear interpolation.
! ......   IRESET=0  : keep old J to start search.
! ......   IRESET=1  : reset J to 1.

	REAL X(N),Y(N)
	INTEGER JEXTRAP
	IF (IRESET.EQ.1) J=1
!	J = 1
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

		if ( (y(jr).gt.1.e-30) .and. (y(jl).gt.1.e-30)&
     		.and. (x(jr).gt.1.e-30) .and. (x(jl).gt.1.e-30) ) then
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
!	J=1
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
