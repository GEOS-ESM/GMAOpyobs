! **********************************************
!
!  Simple f77 wrapper for the Python interface.
!
! **********************************************

 subroutine GetSimpleVeg (lons,lats,DirName,nobs,iVeg)

   use VegType_Mod

   implicit none

   character(len=*), intent(in)  :: DirName
   integer,          intent(in)  :: nobs
   real,             intent(in)  :: lons(nobs)
   real,             intent(in)  :: lats(nobs)
   integer,          intent(out) :: iVeg(nobs)

!                       ----
 
 ! Local variables
 ! ---------------
   type(VegType) :: veg

   call VegType_Initialize(veg, DirName)
   call VegType_GetSimple(veg, lons, lats, nobs, iVeg)
   call VegType_Finalize(veg)
 
 end subroutine GetSimpleVeg

!..........................................................................

 subroutine GetDetailedVeg (lons,lats,DirName,nobs,iVeg)

   use VegType_Mod

   implicit none

   character(len=*), intent(in)  :: DirName
   integer,          intent(in)  :: nobs
   real,             intent(in)  :: lons(nobs)
   real,             intent(in)  :: lats(nobs)
   integer,          intent(out) :: iVeg(nobs)

!                       ----
 
 ! Local variables
 ! ---------------
   type(VegType) :: veg

   call VegType_Initialize(veg, DirName)
   call VegType_GetDetailed(veg, lons, lats, nobs, iVeg)
   call VegType_Finalize(veg)
 
 end subroutine GetDetailedVeg

