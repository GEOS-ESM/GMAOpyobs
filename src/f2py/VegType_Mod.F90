Module VegType_Mod

!BOP
!
! !MODULE: VegType_Mod --- Determines vegetation types at (lat,lon)
!
! !DESCRIPTION: 
!
!  Contain routines to determine the vegetation type given a list of
!  longitudes and latitudes.
!
!  Code derived from a stand alone program provided by Saulo Freitas.
!
!EOP

 Implicit NONE
 PRIVATE

 !PUBLIC METHODS:

  Public VegType_Initialize    ! Get it going
  Public VegType_GetDetailed   ! IGBP land cover types
  Public VegType_GetSimple     ! Simplified, aggregated types
  Public VegType_Finalize    

 !PUBLIC DATA TYPES:
  Public VegType
  type VegType 
     character(len=512) :: DirName ! directory where database resides
  end type VegType
   

  CONTAINS

  subroutine VegType_Initialize(this, DirName)
    implicit NONE
    type(VegType),    intent(out) :: this
    character(len=*), intent(in)  :: DirName
    this%DirName = DirName
  end subroutine VegType_Initialize

  subroutine VegType_GetDetailed(this, lons, lats, nobs, iVeg)
    implicit NONE
    type(VegType), intent(in)  :: this
    integer,       intent(in)  :: nobs
    real,          intent(in)  :: lons(nobs)
    real,          intent(in)  :: lats(nobs)
    integer,       intent(out) :: iVeg(nobs)
    call define_veg(nobs,lats,lons,iVeg,this%DirName)
  end subroutine VegType_GetDetailed

  subroutine VegType_GetSimple(this, lons, lats, nobs, iVeg)
    implicit NONE
    type(VegType), intent(in)  :: this
    integer,       intent(in)  :: nobs
    real,          intent(in)  :: lons(nobs)
    real,          intent(in)  :: lats(nobs)
    integer,       intent(out) :: iVeg(nobs)
    integer :: iVeg_(nobs)
    call define_veg(nobs,lats,lons,iVeg_,this%DirName)
    call agreg_veg(nobs,lats,lons,iVeg_,iVeg)
  end subroutine VegType_GetSimple

  subroutine VegType_Finalize(this)
    implicit NONE
    type(VegType), intent(inout) :: this
    this%DirName = '/dev/null'
  end subroutine VegType_Finalize

!---------------------------------------------------------------

!                   ---------------------
!                   Code by Saulo Freitas
!                   ---------------------
!

!-------------------------------------------------------------------
!
subroutine  agreg_veg(nfocos,qlat,qlon,qveg,qveg_agreg)
implicit none
integer j, ifoc,nfocos
real,    dimension(nfocos) :: qlat,qlon
integer, dimension(nfocos) :: qveg, qveg_agreg
integer :: catb(1:17)
data catb/ &
           2, 1, 2, 2                 & !floresta tropical 2 and 4 / extra trop fores 1,3,5
         , 2, 3, 3, 3, 3              & !cerrado/woody savanna :6 a 9
         , 4, 4, 4, 4, 4, -15, 4, -17 / !pastagem/lavouras: 10 ...

integer :: agreg(nfocos) ! local workspace

do ifoc=1,nfocos
   j = qveg(ifoc)
   if ( j > 0 .and. j < 18 ) then 
      qveg_agreg(ifoc) = catb( j )
   else
      qveg_agreg(ifoc) = 0
   end if
enddo

! classify the IGBP 'Evergreen Broadleaf' as extra-tropical 
! if it is outside of the [30S, 30N] zone
agreg = qveg_agreg
where ( (agreg==1) .and. ((qlat < -30.) .or. (qlat > 30.)) )
   qveg_agreg = 2 ! extra-tropical forests
end where

return
end subroutine agreg_veg

!----------------------------------------------------------------------------

subroutine  define_veg(nfocos,qlat,qlon,qveg,OFN)
integer                    :: nfocos
real,    dimension(nfocos) :: qlat,qlon
real * 4                   :: deltax,deltay,deltallo
real * 4                   :: deltaxp,deltayp
integer, dimension(nfocos) :: qveg
integer, allocatable :: DATO(:)
CHARACTER*512 TITLE
CHARACTER (len=*)::  OFN
character(len=1), allocatable:: cdato(:)
integer :: LB,iblksizo,no,isbego,iwbego,np,mof,iodim

!--------Diretorio onde se encontram os arquivos de vegetacao e o prefixo'
!OFN='./data_veg/IGBP'

LB=lastchar(OFN)
TITLE=OFN(1:LB)//'HEADER'
lb=lastchar(title)
call vegopenf(29,title(1:lb),'formatted','old','read',0)
read(29,2)iblksizo,no,isbego,iwbego
2    format(4i5)
close(29)


deltax=1000. ! 1 km resolucao
deltay=1000. ! 1 km resolucao

deltallo=float(iblksizo)/float(no-1)
np=min(10,max(1,int(deltax/(deltallo*111000.))))
deltaxp=deltax/float(np)
deltayp=deltay/float(np)
mof=4
iodim=no*no*mof
allocate (dato(iodim))
allocate (cdato(no*no))

 call le_veg(nfocos,NO,MOF,NP,DELTALLO,DELTAXP,DELTAYP,IBLKSIZO  &
   ,ISBEGO,IWBEGO,DATO,qveg,OFN,cdato,qlat,qlon)

deallocate (dato,cdato)

RETURN
END subroutine define_veg

!------------------------------------------------------------------------

subroutine le_veg(nfocos,NO,MOF,NP,DELTALLO,DELTAXP,DELTAYP,IBLKSIZO  &
                 ,ISBEGO,IWBEGO,DATO,qveg,OFN,cdato,qlat,qlon)

integer NO,NP,MOF,nfocos,IBLKSIZO,ISBEGO,IWBEGO,NONO
integer nofr,iof,ifoc
integer DATO(NO,NO,MOF)
character*1 cdato(no,no)
integer, dimension(nfocos) :: qveg
real,    dimension(nfocos) :: qlat,qlon
real * 4                   :: deltax,deltay,deltallo
real * 4                   :: deltaxp,deltayp,glatp,glonp
real * 4                   :: RIO,RJO
integer                    :: IO,JO

integer iso(100),iwo(100)
character(len=*) ofn
character*512 title3
character*3 title1
character*4 title2
integer datp,ISOC,IWOC,IOFR,JOFR,ISOCPT,ISOCPO,IWOCPH,IWOCPTIWOCPO,LB
integer i,j,IWOCPT,IWOCPO
logical l1,l2

NONO=NO*NO
NOFR=0
DO IOF=1,MOF
   ISO(IOF)=0
   IWO(IOF)=0
ENDDO

do ifoc=1,nfocos
 glatp=qlat(ifoc)
 glonp=qlon(ifoc)
            ISOC=(INT((GLATP-FLOAT(ISBEGO))/FLOAT(IBLKSIZO)+200.)  &
               -200)*IBLKSIZO+ISBEGO
            IWOC=(INT((GLONP-FLOAT(IWBEGO))/FLOAT(IBLKSIZO)+400.)  &
               -400)*IBLKSIZO+IWBEGO

            DO IOFR=1,NOFR
               JOFR=IOFR
               IF(ISO(IOFR).EQ.ISOC.AND.IWO(IOFR).EQ.IWOC)GO TO 10
            ENDDO
            ISOCPT=ABS(ISOC)/10
            ISOCPO=ABS(ISOC)-ISOCPT*10
            IWOCPH=ABS(IWOC)/100
            IWOCPT=(ABS(IWOC)-IWOCPH*100)/10
            IWOCPO=ABS(IWOC)-IWOCPH*100-IWOCPT*10
            IF(ISOC.GE.0)THEN
               WRITE(TITLE1,'(2I1,A1)')ISOCPT,ISOCPO,'N'
            ELSE
               WRITE(TITLE1,'(2I1,A1)')ISOCPT,ISOCPO,'S'
            ENDIF
            IF(IWOC.GE.0)THEN
               WRITE(TITLE2,'(3I1,A1)')IWOCPH,IWOCPT,IWOCPO,'E'
            ELSE
               WRITE(TITLE2,'(3I1,A1)')IWOCPH,IWOCPT,IWOCPO,'W'
            ENDIF
            LB=lastchar(OFN)
            TITLE3=OFN(1:LB)//TITLE1//TITLE2
            LB=lastchar(TITLE3)
            INQUIRE(FILE=TITLE3(1:LB),EXIST=L1,OPENED=L2)
            IF(.NOT.L1)THEN
!                     PRINT*, ' FILE',TITLE3(1:LB),' DOES NOT EXIST '
!     +                     'WATER IS ASSUMED'

               DATP = 0
               GO TO 20
            ENDIF
            IF(NOFR.GE.MOF)THEN
               DO IOF=1,MOF
                  ISO(IOF)=0
                  IWO(IOF)=0
               ENDDO
               NOFR=0
            ENDIF
            NOFR=NOFR+1
            JOFR=NOFR

!print *,' Calling veopen ********************'
!print*, title3(1:lb)//char(0),'r'//char(0)
!print*, title3(1:lb)//char(0)
!print*, char(0)
!print*, 'r'//char(0)


!!!            print *, 'file = ', trim(title3(1:lb)//char(0))
            call vegopen(title3(1:lb)//char(0),'r'//char(0))
            call vegchar(4,no*no,cdato)
            call vegclose() 

            do j=1,no
               do i=1,no
                  dato(i,j,nofr) = ichar(cdato(i,j))
               enddo
            enddo

            ISO(NOFR)=ISOC
            IWO(NOFR)=IWOC
10             CONTINUE
            RIO=(GLONP-FLOAT(IWOC))/DELTALLO+1.5
            RJO=(GLATP-FLOAT(ISOC))/DELTALLO+1.5
!
            IO=INT(RIO)
            JO=INT(RJO)

            DATP=DATO(IO,JO,JOFR)

!	    print*,title3(1:lb),GLATP,GLONP,DATP
!           if(ifoc.lt.17)then
!	    print*,ifoc,IO,JO,DATP,no,IWOC,ISOC
!            endif

20               CONTINUE

            qveg(ifoc)=datp

!srf- corrige veg=0
  	    if(qveg(ifoc)==0)  qveg(ifoc)= 17	   

enddo

return
end subroutine le_veg
!----------------------------------------------------------------------------------------
subroutine define_bas_olson(bas_olson,nveg_olson,nfocos,bas_by_fire,qveg_olson)
implicit none
integer nveg_olson,nfocos,i
real, dimension(nveg_olson) :: bas_olson
real,    dimension(4,0:96) :: olson_data
real,    dimension(nfocos) :: bas_by_fire
integer, dimension(nfocos) :: qveg_olson
data olson_data/  &
!------------------------------------------------------------------------------------------------------
!Medium carbon density (Kg C/m2) 	
!Revised medium carbon density (Kg C/m2) 	
!Minimum carbon density (Kg C/m2) 	
!Maximum carbon density (Kg C/m2) 	
!Ecosystem codes
!Ecosystem complexes 
!------------------------------------------------------------------------------------------------------
! MED   REV     MIN    MAX     Olson Global Ecosystem Legend (glcc/usgs versao 2)
!                             !X = adaptado, nao existente no dado original 
                              !    (http://cdiac.ornl.gov)
 0.00,  0.00,	0.0,   0.0, & !X 0 INTERRUPTED AREAS (GLOBAL GOO
 1.00,  0.80,	0.6,   2.0, & !X 1 URBAN			
 1.00,  0.80,	0.6,   2.0, & !X 2 LOW SPARSE GRASSLAND 	
16.00, 13.00,  12.0,  20.0, & !X 3 CONIFEROUS FOREST		
10.00,  7.00,	6.0,  14.0, & !X 4 DECIDUOUS CONIFER FOREST	
10.00,  9.00,	8.0,  14.0, & !X 5 DECIDUOUS BROADLEAF FOREST	
10.00,  9.00,	8.0,  14.0, & !X 6 EVERGREEN BROADLEAF FORESTS  
 3.00,  3.00,	2.0,   5.0, & !X 7 TALL GRASSES AND SHRUBS	
 0.40,  0.30,	0.2,   1.0, & !X  8 BARE DESERT  		
 0.50,  0.50,	0.0,   1.2, & !X 9 UPLAND TUNDRA		
 2.00,  2.00,	1.0,   3.0, & !X 10 IRRIGATED GRASSLAND  	
 0.40,  0.30,	0.2,   1.0, & !X 11 SEMI DESERT  		
 0.00,  0.00,	0.0,   0.0, & !  12 GLACIER ICE  		
 0.00,  0.00,	0.0,   0.0, & !  13 WOODED WET SWAMP		
 0.00,  0.00,	0.0,   0.0, & !  14 INLAND WATER 		
 0.00,  0.00,	0.0,   0.0, & !  15 SEA WATER			
1.30,  0.90,	0.5,   3.0, & !X 16 SHRUB EVERGREEN		
1.30,  0.90,	0.5,   3.0, & !X 17 SHRUB DECIDUOUS		
10.00,  7.00,	6.0,  14.0, & !X 18 MIXED FOREST AND FIELD	
!15.00, 12.00,	4.0,  25.0, & !X 19 EVERGREEN FOREST AND FIELDS 
15.00, 14.62,	4.0,  25.0, & !X 19 EVERGREEN FOREST AND FIELDS 
 8.00,  6.00,	4.0,  11.0, & ! 20 COOL RAIN FOREST			      
 8.00,  6.00,	4.0,  11.0, & ! 21 CONIFER BOREAL FOREST		      
16.00, 13.00,  12.0,  20.0, & ! 22 COOL CONIFER FOREST  		      
10.00,  7.00,	6.0,  14.0, & ! 23 COOL MIXED FOREST			         
10.00,  7.00,	6.0,  14.0, & ! 24 MIXED FOREST 			         
10.00,  9.00,	8.0,  14.0, & ! 25 COOL BROADLEAF FOREST		      
10.00,  9.00,	8.0,  14.0, & ! 26 DECIDUOUS BROADLEAF FOREST		      
16.00, 13.00,  12.0,  20.0, & ! 27 CONIFER FOREST			      
 5.00,  5.00,	1.0,  15.0, & ! 28 MONTANE TROPICAL FORESTS		      
15.00, 12.00,	4.0,  25.0, & ! 29 SEASONAL TROPICAL FOREST		       
 1.00,  0.70,	0.4,   2.0, & ! 30 COOL CROPS AND TOWNS 		      
 1.00,  0.80,	0.6,   2.0, & ! 31 CROPS AND TOWN			      
 7.00,  6.00,	5.0,   9.0, & ! 32 DRY TROPICAL WOODS			      
!15.00, 12.00,	4.0,  25.0, & ! 33 TROPICAL RAINFOREST  		       
15.00, 14.62,	4.0,  25.0, & ! 33 TROPICAL RAINFOREST  		       
10.00,  7.00,	6.0,  14.0, & ! 34 TROPICAL DEGRADED FOREST		      
 3.00,  3.00,	2.0,   4.0, & ! 35 CORN AND BEANS CROPLAND		      
 3.00,  3.00,	2.0,   4.0, & ! 36 RICE PADDY AND FIELD 		      
 2.00,  2.00,	1.0,   3.0, & ! 37 HOT IRRIGATED CROPLAND		       
 2.00,  2.00,	1.0,   3.0, & ! 38 COOL IRRIGATED CROPLAND		       
 2.00,  2.00,	1.0,   3.0, & ! 39 COLD IRRIGATED CROPLAND		       
 1.00,  0.80,	0.6,   2.0, & ! 40 COOL GRASSES AND SHRUBS		      
 1.30,  0.90,	0.5,   3.0, & ! 41 HOT AND MILD GRASSES AND SHRUBS	      
 1.00,  1.00,	0.5,   4.0, & ! 42 COLD GRASSLAND			      
! 3.00,  3.00,	2.0,   5.0, & ! 43 SAVANNA (WOODS)			      
 3.00,  3.10,	2.0,   5.0, & ! 43 SAVANNA (WOODS)			      
 2.00,  2.00,	1.0,   6.0, & ! 44 MIRE, BOG, FEN			      
 3.00,  2.00,	1.0,   6.0, & ! 45 MARSH WETLAND			      
 4.00,  3.00,	2.0,   8.0, & ! 46 MEDITERRANEAN SCRUB  		      
 4.00,  3.00,	2.0,   8.0, & ! 47 DRY WOODY SCRUB			      
 5.00,  4.00,	2.0,  10.0, & ! 48 DRY EVERGREEN WOODS  		      
 0.40,  0.30,	0.2,   1.0, & ! 49 VOLCANIC ROCK			       
 0.05,  0.05,	0.0,   0.2, & ! 50 SAND DESERT  			      
 0.40,  0.30,	0.2,   1.0, & ! 51 SEMI DESERT SHRUBS			       
 0.60,  0.60,	0.3,   1.0, & ! 52 SEMI DESERT SAGE			      
 0.50,  0.50,	0.0,   1.2, & ! 53 BARREN TUNDRA			      
 0.50,  0.50,	0.0,   1.2, & ! 54 COOL SOUTHERN HEMISPHERE MIXED FORESTS     
 4.00,  3.00,	2.0,   5.0, & ! 55 COOL FIELDS AND WOODS		       
 5.00,  4.00,	4.0,   8.0, & ! 56 FOREST AND FIELD			       
 5.00,  4.00,	4.0,   8.0, & ! 57 COOL FOREST AND FIELD		       
! 4.00,  3.00,	2.0,   5.0, & ! 58 FIELDS AND WOODY SAVANNA		       
 4.00,  1.9,	2.0,   5.0, & ! 58 FIELDS AND WOODY SAVANNA !srf : palacios-orueta
! 4.00,  3.00,	2.0,   6.0, & ! 59 SUCCULENT AND THORN SCRUB		      
 4.00,  3.70,	2.0,   6.0, & ! 59 SUCCULENT AND THORN SCRUB		      
11.00,  8.00,	6.0,  14.0, & ! 60 SMALL LEAF MIXED WOODS		      
11.00,  8.00,	6.0,  14.0, & ! 61 DECIDUOUS AND MIXED BOREAL FOREST	      
 5.00,  5.00,	2.0,   8.0, & ! 62 NARROW CONIFERS			      
 2.00,  2.00,	1.0,   5.0, & ! 63 WOODED TUNDRA			      
 1.50,  1.00,	1.0,   2.0, & ! 64 HEATH SCRUB  			      
 3.00,  3.00,	0.0,  10.0, & ! 65 COASTAL WETLAND - NW 		       
 3.00,  3.00,	0.0,  10.0, & ! 66 COASTAL WETLAND - NE 		       
 3.00,  3.00,	0.0,  10.0, & ! 67 COASTAL WETLAND - SE 		       
 3.00,  3.00,	0.0,  10.0, & ! 68 COASTAL WETLAND - SW 		       
 0.50,  0.50,	0.0,   1.2, & ! 69 POLAR AND ALPINE DESERT		      
 0.50,  0.50,	0.0,   1.2, & ! 70 GLACIER ROCK 			      
 0.40,  0.30,	0.2,   1.0, & ! 71 SALT PLAYAS  			       
 3.00,  2.00,	1.0,   6.0, & ! 72
 0.00,  0.00,	0.0,   0.0, & !   73 WATER AND ISLAND FRINGE		 
 0.00,  0.00,	0.0,   0.0, & !   74 LAND, WATER, AND SHORE		 
 0.00,  0.00,	0.0,   0.0, & !   75 LAND AND WATER, RIVERS		 
 0.00,  0.00,	0.0,   0.0, & !   76 CROP AND WATER MIXTURES		 
10.00,  7.00,   6.0,  14.0, & !X  77 SOUTHERN HEMISPHERE CONIFERS	  
10.00,  7.00,   6.0,  14.0, & !X  78 SOUTHERN HEMISPHERE MIXED FOREST   
10.00,  7.00,   6.0,  14.0, & !X  79 WET SCLEROPHYLIC FOREST		  
 0.00,  0.00,	0.0,   0.0, & !   80 COASTLINE FRINGE		 
 0.00,  0.00,	0.0,   0.0, & !   81 BEACHES AND DUNES  	 
 0.00,  0.00,	0.0,   0.0, & !   82 SPARSE DUNES AND RIDGES		 
 0.00,  0.00,	0.0,   0.0, & !   83 BARE COASTAL DUNES 	 
 0.00,  0.00,	0.0,   0.0, & !   84 RESIDUAL DUNES AND BEACHES  
 0.00,  0.00,	0.0,   0.0, & !   85 COMPOUND COASTLINES		 
 0.00,  0.00,	0.0,   0.0, & !   86 ROCKY CLIFFS AND SLOPES		 
 0.00,  0.00,	0.0,   0.0, & !   87 SANDY GRASSLAND AND SHRUBS  
10.00,  7.00,   6.0,  14.0, & !X  88 BAMBOO				 
10.00,  7.00,   6.0,  14.0, & !X  89 MOIST EUCALYPTUS  		 
!15.00, 12.00,   4.0,  25.0,& !X  90 RAIN GREEN TROPICAL FOREST	 
15.00, 14.62,   4.0,  25.0, & !X  90 RAIN GREEN TROPICAL FOREST	 
!3.00,  3.00,	2.0,   5.0, & !X  91 WOODY SAVANNA			
 3.00,  5.00,	2.0,   5.0, & !X  91 WOODY SAVANNA			
 3.00,  3.00,	2.0,   5.0, & !X  92 BROADLEAF CROPS			
 1.00,  0.80,	0.6,   2.0, & !X  93 GRASS CROPS			
 1.00,  0.80,	0.6,   2.0, & !X  94 CROPS, GRASS, SHRUBS		
 3.00,  3.00,	2.0,   5.0, & !X  95 EVERGREEN TREE CROP
 3.00,  3.00,	2.0,   5.0  / !X  96 DECIDUOUS TREE CROP
!-----------------------------------------------------------------------------------------------

bas_by_fire(:)=0.

do i=1,nveg_olson-1
 bas_olson(i)   = 2.*olson_data(2,i)  ! total biomass ~ 2. * Carbon density 
enddo

!- define densidade de biomassa acima do solo no array bas_by_fire
do i=1,nfocos

 if(qveg_olson(i) > 96 ) then
  !print*, 'erro nos dados de vegetacao de OLSON - mudando para oceano'
  qveg_olson(i)=15
 endif
 bas_by_fire(i) = bas_olson(qveg_olson(i))
 !if(bas_by_fire(i) < 0.001)  print*,'OLSON=',i,qveg_olson(i),bas_by_fire(i)
enddo

return
end subroutine define_bas_olson
!------------------------------------------------------------------------------------
subroutine vegopenf(iunit, filenm, formt, stat, act, iclob)

! replaces old jclopen and jclget
! files are overwritten unless iclob (ICLOBBER) set to 1

implicit none

integer iunit, iclob
character(len=*) filenm, formt, stat, act
logical exans
integer, external :: lastchar

inquire(FILE=filenm,EXIST=exans)
if(exans.and.iclob.eq.0.and.  &
     ((act == 'WRITE') .or. (act == 'write'))) then
   print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   print*,'!!!   trying to open file name :'
   print*,'!!!       ',filenm
   print*,'!!!   but it already exists. run is ended.'
   print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   call exit(1)
endif

! print*,'filenm,formt,stat2=',filenm,'  ',formt,'  ',stat
open(iunit,STATUS=stat,FILE=filenm,FORM=formt)
!print*,'F_open - ',filenm(1:lastchar(filenm))

return
end subroutine vegopenf

!-------------------------------------------------------------------
integer function lastchar(str)
character(len=*) str
integer ln,n

! returns last non-blank character position from a string

ln=len(str)
do n=ln,1,-1
   if(str(n:n).ne.' ') then
      lastchar=n
      return
   endif
enddo
lastchar=0

return
end function lastchar

!-------------------------------------------------------------------

end Module VegType_Mod
