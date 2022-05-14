!
! Simple f2py interface
!

subroutine sgp4Track(lon, lat, n, name, nymd, nhms, dtsec, rc )

  use SGP4_Mod
  implicit NONE

  character(len=*), intent(in) :: name      ! Satellite name or TLE file name
  integer, intent(in)          :: n         ! size of output
  integer, intent(in)          :: nymd(2)   ! Beginning/ending date: YYYYMMDD
  integer, intent(in)          :: nhms(2)   ! Beginning/ending time: HHMMSS
  INTEGER, intent(in)          :: dtsec     ! Time step [secs]

! !OUTPUT PARAMETERS:

  real*8, intent(out)        :: lon(n)   ! Ground track longitudes [degrees]  
  real*8, intent(out)        :: lat(n)   ! Ground track latitudes  [degrees]  
  integer, intent(out)         :: rc       ! Error code = 0 all is well

  integer :: m
  real*8, pointer :: tlon(:), tlat(:)

  CALL SGP4_Track(tlon, tlat, name, nymd, nhms, dtsec, rc)

  m = size(tlon)
  if ( m > n ) then
     rc = 1
     m = n
  end if

  lon = -9999.0
  lat = -9999.0
  lon(1:m) = tlon(1:m)
  lat(1:m) = tlat(1:m)

  deallocate(tlon,tlat)

end subroutine sgp4Track


