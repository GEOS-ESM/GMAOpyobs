!  $Id: glint_mod.F90,v 1.2 2017/08/24 15:08:03 adasilva Exp $

  MODULE glint_mod 
!BOP

! !MODULE: glint_mod  --- Glint angle calculation module 


      IMPLICIT NONE

     PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:

      PUBLIC CALC_sun_zen_az   ! Generate sun azimuth and zenith angles 
      PUBLIC calc_glint        ! Generates glint angle
 
#ifdef __PROTEX__

  !DESCRIPTION:  
  
   This module provides necessary functions for calculating
   Glint angle for satellite and sun. 
 
#endif

!  !REVISION HISTORY: 

!     10Dec2009   Albayrak  Initial implementation.
!                 See individual routines for further detail.

!EOP


! !REVISION HISTORY
!
!  10Dec2009   Albayrak  Initial implementation 
!              (see indivudual subroutines for further information).


 

!     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Internal constants (not public!)
!     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
       INTEGER, PARAMETER :: dp=SELECTED_REAL_KIND(15,307)


!     2-Other Coef
       REAL(dp), PARAMETER :: pi=3.141592653589793
!       REAL(dp), PARAMETER :: earth_radius=6371.0  !sphere
       REAL(dp), PARAMETER :: earth_radius=6378.137D0  !km -- sphere
       REAL(dp), PARAMETER :: f = 1.D0 / 298.257       ! earth flattening factor
       REAL(dp), PARAMETER :: myeps = 0.0000001        !for zero checks
       REAL(dp), PARAMETER :: mu = 398600.4415         !km^3/s^2  Gravitational constant of Earth
       REAL(dp), PARAMETER :: cdeg2km = earth_radius * (pi / 180.0) ! constant for deg2km 
       REAL(dp), PARAMETER :: ckm2deg = (1.0/earth_radius) * (180.0/pi)
       REAL(dp), PARAMETER :: cr2d = 180/pi ! radian2degree
       REAL(dp), PARAMETER :: cd2r = pi/180 ! deg2rad
       REAL(dp), PARAMETER :: cd2s = 24.0*60.0*60.0 ! day2second

!     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Definition of overloaded functions 
!     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

        interface calc_glint
              module procedure calc_glint1  ! sat pos is given in lat lon
              module procedure calc_glint2  ! sat pos is given in ECEF
        end interface            
 
!     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Definition of structs ...
!     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


      type julian_cal
          real(dp) :: day
          real(dp) :: ephemeris_day
          real(dp) :: century 
          real(dp) :: ephemeris_century 
          real(dp) :: ephemeris_millenium
      end type julian_cal


      type earth_heliocen_pos
          real(dp) :: longitude
          real(dp) :: latitude
          real(dp) :: radius
      end type earth_heliocen_pos


     type sun_geocentric_pos
          real(dp) :: longitude
          real(dp) :: latitude    
     end type sun_geocentric_pos

     type topocentric_sun_pos
          real(dp) :: rigth_ascension_parallax
          real(dp) :: rigth_ascension
          real(dp) :: declination
     end type topocentric_sun_pos

     TYPE vector
          REAL :: x, y, z
     END TYPE vector

           
!     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     ERROR codes 
!     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!      0 - OK
!      3 - ALLOCATE (memory) problem
!      5 - satelite is not defined
!     90 - Swath width is not given correctly
!     99 - Satellite name is not correct




CONTAINS




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP
 
! !IROUTINE:  calc_glint --- Calculates glint angle
!
! !IIROUTINE: calc_glint  --- Main function to calculate Glint angle

!
! !INTERFACE:
!

SUBROUTINE calc_glint1(GLINT_ANGLE, nymd, nhms, long, lat, alt, tlon, tlat, talt)
       IMPLICIT NONE
! !ARGUMENTS:
!
       real(dp), intent(out) :: GLINT_ANGLE ! Glint angle (deg)
       integer,  intent(in)  :: nymd        ! Date 
       integer,  intent(in)  :: nhms        ! Time 
       real(dp), intent(in)  :: long        ! Longitude of the pixel (deg), 
       real(dp), intent(in)  :: lat         ! Latitude of the pixel (deg), 
       real(dp), intent(in)  :: alt         ! Altitude of the pixel (meters) above mean sea level
       real(dp), intent(in)  :: tlon        ! Longitude of sat (deg), 
       real(dp), intent(in)  :: tlat        ! Latitude of sat (deg), 
       real(dp), intent(in)  :: talt        ! Altitude of sat (km) above mean sea level


       real(dp), dimension(1:3)  :: pos  ! Satellite position vec (3x1)
       real(dp) :: x, y, z

#ifdef ___PROTEX___
!
   !DESCRIPTION:

   This function calls sun location, and sensor location routines, 
   and obtaines related variables. After that calculates glint angle
#endif
!
!EOP



!      Other variables
       REAL(dp)      :: fraction2day
       integer, DIMENSION(1:6)  :: timevec       ! (/Year, Month, Day, Hour, Minute, Seconds/)
       REAL(dp)  :: sunzenith1, sunazimuth1
       REAL(dp)  :: senzenith1, senazimuth1
       REAL(dp)  :: MTHET, MPHI, MTHET0, MPHI0, MDPHI,  MSCATT




!translate -- tlat, tlon to xyz and then assign to pos vector

       print*, "tlat, tlan, talt ", tlat, tlon, talt
       CALL LLA2ECEF(tlat*cd2r, tlon*cd2r, talt*1000, x, y, z)
       pos(1) = x/1000.0
       pos(2) = y/1000.0
       pos(3) = z/1000.0 
       print*, " **** __ ", pos 


       CALL get_timevec(nymd, nhms, fraction2day,  timevec)
       CALL CALC_sun_zen_az(sunzenith1, sunazimuth1, nymd, nhms, long, lat, alt)
       CALL calc_sensor_zen_az(senzenith1, senazimuth1, long,lat, alt, pos)

!      I assume here
!      theta  sensor zenith
!      MPHI   sensor azimuth
!      theta0 solar zenith
!      MPHI0  solar azimuth

       MTHET  = senzenith1
       MPHI   = senazimuth1
       MTHET0 = sunzenith1
       MPHI0  = sunazimuth1

!      Relative Azimuth
       MDPHI = abs(MPHI0 - MPHI -180.0)
       if (MDPHI>360.0) then
           MDPHI=mod(MDPHI,360.0)
       end if
       if (MDPHI>180.0) then
           MDPHI=360.0-MDPHI
       end if

       if ( (MTHET0>0.0).AND.(MTHET>0.0).AND.(MDPHI>0.0) ) then
          ! Scattering
          MSCATT = -cos(MTHET0*cd2r)*cos(MTHET*cd2r) +           &
                    sin(MTHET0*cd2r)*sin(MTHET*cd2r) * cos(MDPHI*cd2r)
          MSCATT = acos(MSCATT)*cr2d
          ! Glint:
          GLINT_ANGLE=(cos(MTHET0*cd2r))*(cos(MTHET*cd2r)) +     &
                      ((sin(MTHET0*cd2r))*(sin(MTHET*cd2r)) * ( cos(MDPHI*cd2r)))
          GLINT_ANGLE = (acos(GLINT_ANGLE))*cr2d
       end if

END SUBROUTINE calc_glint1




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP
 
! !IROUTINE:  calc_glint --- Calculates glint angle
!
! !IIROUTINE: calc_glint  --- Main function to calculate Glint angle

!
! !INTERFACE:
!

SUBROUTINE calc_glint2(GLINT_ANGLE, nymd, nhms, long, lat, alt, pos)
       IMPLICIT NONE
! !ARGUMENTS:
!
       real(dp), intent(out) :: GLINT_ANGLE ! Glint angle (deg)
       integer,  intent(in)  :: nymd        ! Date 
       integer,  intent(in)  :: nhms        ! Time 
       real(dp), intent(in)  :: long        ! Longitude of the pixel (deg), 
       real(dp), intent(in)  :: lat         ! Latitude of the pixel (deg), 
       real(dp), intent(in)  :: alt         ! Altitude of the pixel (meters) above mean sea level
       real(dp), dimension(1:3), intent(in)  :: pos  ! Satellite position vec (3x1)

#ifdef ___PROTEX___
!
   !DESCRIPTION:

   This function calls sun location, and sensor location routines, 
   and obtaines related variables. After that calculates glint angle
#endif
!
!EOP



!      Other variables
       REAL(dp)      :: fraction2day
       integer, DIMENSION(1:6)  :: timevec       ! (/Year, Month, Day, Hour, Minute, Seconds/)
       REAL(dp)  :: sunzenith1, sunazimuth1
       REAL(dp)  :: senzenith1, senazimuth1
       REAL(dp)  :: MTHET, MPHI, MTHET0, MPHI0, MDPHI,  MSCATT

       CALL get_timevec(nymd, nhms, fraction2day,  timevec)
       CALL CALC_sun_zen_az(sunzenith1, sunazimuth1, nymd, nhms, long, lat, alt)
       CALL calc_sensor_zen_az(senzenith1, senazimuth1, long,lat, alt, pos)

!      I assume here
!      theta  sensor zenith
!      MPHI   sensor azimuth
!      theta0 solar zenith
!      MPHI0  solar azimuth

       MTHET  = senzenith1
       MPHI   = senazimuth1
       MTHET0 = sunzenith1
       MPHI0  = sunazimuth1

!      Relative Azimuth
       MDPHI = abs(MPHI0 - MPHI -180.0)
       if (MDPHI>360.0) then
           MDPHI=mod(MDPHI,360.0)
       end if
       if (MDPHI>180.0) then
           MDPHI=360.0-MDPHI
       end if

       if ( (MTHET0>0.0).AND.(MTHET>0.0).AND.(MDPHI>0.0) ) then
          ! Scattering
          MSCATT = -cos(MTHET0*cd2r)*cos(MTHET*cd2r) +           &
                    sin(MTHET0*cd2r)*sin(MTHET*cd2r) * cos(MDPHI*cd2r)
          MSCATT = acos(MSCATT)*cr2d
          ! Glint:
          GLINT_ANGLE=(cos(MTHET0*cd2r))*(cos(MTHET*cd2r)) +     &
                      ((sin(MTHET0*cd2r))*(sin(MTHET*cd2r)) * ( cos(MDPHI*cd2r)))
          GLINT_ANGLE = (acos(GLINT_ANGLE))*cr2d
       end if

END SUBROUTINE calc_glint2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP

! !IIROUTINE: calc_sensor_zen_az --- Calculates sensor (satellite) zenith, and azimuth

! !INTERFACE:
!
 SUBROUTINE calc_sensor_zen_az(senz, sena, long, lat, alt, pos)
 
 
        IMPLICIT NONE

!
! !ARGUMENTS:
! 
       real(dp), intent(out) :: senz, sena        ! sensor zenith and azimuth (deg) 
       real(dp), intent(in)  :: long, lat, alt    ! degrees
       real(dp), dimension(:), intent(in)  :: pos ! satellite position vec (3x1) ECEF

#ifdef ___PROTEX___

    !DESCRIPTION:

    Translated from IDL code to fortran by Arif Albayrak

   

   {\bf pro lonlat2viewangle}
   Calculates solar and sensor view angles from geodetic longitude 
   and latitude, spacecraft position vector, and solar unit vector.
 
   Inputs: \\
   lon	- longitude of pixel in degrees \\
   lat	- longitude of pixel in degrees \\
   pos(3)	- vector from earth center to spacecraft in ECR, km 
 
   Outputs: \\
   senz	- sensor zenith angle of pixel in degrees \\
   sena	- sensor azimuth angle of pixel in degrees 
 
   Notes: \\
   Currently only handles one lon/lat position per call.
 
   Written By: \\
   Bryan Franz, SAIC GSC, NASA/SIMBIOS Project, April 1998. \\
   (with much help from geonav.pro by Fred Patt)
 

#endif 


!EOP


!      Other variables
       REAL(dp), DIMENSION(:,:) :: rmat(1:3,1:3)
       REAL(dp), DIMENSION(1:3) :: up, ea, no, gvec 
       REAL(dp), DIMENSION(1,1:3) :: scvec
       REAL(dp), DIMENSION(1:3,1):: senl
       REAL(dp)                 :: rlon, rlat, upxy, phi, r


       rmat = 0.0 

       rlon = long * pi / 180.D0
       rlat = lat * pi / 180.D0

!      first, we must define the axes (in ecr space) of a
!      pixel-local coordinate system, with the z-axis along
!      the geodetic pixel normal, x-axis pointing east, and
!      y-axis pointing north.

! !   !up = d1_array(cos(rlat) .* cos(rlon),cos(rlat) .* sin(rlon),sin(rlat));
       up(1) = cos(rlat) * cos(rlon)
       up(2) = cos(rlat) * sin(rlon)
       up(3) = sin(rlat)

       !up = [ cos(rlat) .* cos(rlon),cos(rlat) .* sin(rlon),sin(rlat)];
       upxy = sqrt(up(0 +1) * up(0 +1) + up(1 +1) * up(1 +1))
!      ea = d1_array(-up(1 +1) ./ upxy,up(0 +1) ./ upxy,0.0);
       ea(1) = -up(1 +1) / upxy
       ea(2) = up(0 +1)  / upxy
       ea(3) = 0.0
       !no = cross(up,ea)
       no(1) = up(2) * ea(3) - up(3) * ea(2)
       no(2) = up(3) * ea(1) - up(1) * ea(3)
       no(3) = up(1) * ea(2) - up(2) * ea(1)

  !    now we need the pixel-to-spacecraft vector in ecr,
  !    compute geocentric pixel location vector in km and subtract
  !    from spacecraft position.

       !phi = atann(tan(rlat) .* (1 - f).**2);  ! geocentric latitude
       phi = atan(tan(rlat) * (1 - f)**2)    ! geocentric latitude
       r = earth_radius * (1 - f) / sqrt(1 - (2 - f) * f * (cos(phi)**2))  ! dist to earth surface

!      gvec = r .* [cos(phi) .* cos(rlon),cos(phi) .* sin(rlon),sin(phi)];
       gvec(1) = r * (cos(phi) * cos(rlon))
       gvec(2) = r * (cos(phi) * sin(rlon))
       gvec(3) = r * sin(phi)
       scvec(1,1:3) = pos(1:3) - gvec(1:3)

!      now we can transform the pixel-to-spacecraft and sun 
!      vectors into the local frame.
       rmat(0 +1,:) = ea
       rmat(1 +1,:) = no
       rmat(2 +1,:) = up
       senl = matmul(rmat , TRANSPOSE(scvec))    
!      compute the sensor zenith and azimuth

       senz = atan2(sqrt(senl(1,1 ) * senl(1,1) + senl(2,1) * senl(2,1)),senl(3,1)) * cr2d

       if ((senz > 0.05d0)) then
          sena = atan2(senl(1,1), senl(2,1)) * cr2d
       else
          sena = 0.0d0;
       end if

       if ((sena < 0.d0)) then
          sena = sena + 360.0d0;
          !    print,pos
          !    print,scvec
          !    print,gvec
          !    print,lon,lat
          !    print,solz,sola
          !    print,senz,sena
       end if

  END SUBROUTINE calc_sensor_zen_az


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP

! !IIROUTINE: CALC_sun_zen_az --- Calculates sun zenith, and azimuth

! !INTERFACE:
!
SUBROUTINE CALC_sun_zen_az(sun_zen, sun_az, nymd, nhms, long, lat, alt)

       IMPLICIT NONE
!
! !ARGUMENTS:
!
       real(dp), intent(out) :: sun_az, sun_zen  ! degrees 
       integer,  intent(in)  :: nymd             ! date (Ex:20081231)
       integer,  intent(in)  :: nhms             ! time (Ex:015057)
       real(dp), intent(in)  :: long, lat        ! pixel location --> (deg), (deg)
       real(dp), intent(in)  :: alt              ! meters above mean sea level



#ifdef ___PROTEX___

    !DESCRIPTION:

    This function is translated from C++ to Fortran90 by Arif Albayrak

    Original code documentation: \\
    This function computes the sun position (zenith and azimuth angle at the observer 
    location) as a function of the observer local time and position. \\
    It is an implementation of the algorithm presented by Reda et Andreas in: 
    Reda, I., Andreas, A. (2003) Solar position algorithm for solar 
    radiation application. National Renewable Energy Laboratory (NREL) 
    Technical report NREL/TP-560-34302. 
    This document is avalaible at www.osti.gov/bridge \\
    This algorithm is based on numerical approximation of the exact equations.
    The authors of the original paper state that this algorithm should be
    precise at +/- 0.0003 degrees. I have compared it to NOAA solar table
   (http://www.srrb.noaa.gov/highlights/sunrise/azel.html) and to USNO solar
   table (http://aa.usno.navy.mil/data/docs/AltAz.html) and found very good
   correspondance (up to the precision of those tables), except for large
   zenith angle, where the refraction by the atmosphere is significant 
   (difference of about 1 degree). Note that in this code the correction 
   for refraction in the atmosphere as been implemented for a temperature 
   of 10C (283 kelvins) and a pressure of 1010 mbar. See the subfunction 
   $sun_topocentric_zenith_angle_calculation$ for a possible modification 
   to explicitely model the effect of temperature and pressure as describe
   in Reda and Andreas (2003). 

   History \\
   09/03/2004  Original creation by Vincent Roy (vincent.roy@drdc-rddc.gc.ca) \\
   10/03/2004  Fixed a bug in $julian_calculation$ subfunction (was
               incorrect for year 1582 only), Vincent Roy\\
   18/03/2004  Correction to the header (help display) only. No changes to   
               the code. (changed the elevation field in location structure
               information to altitude), Vincent Roy \\
   13/04/2004  Following a suggestion from Jody Klymak (jklymak@ucsd.edu),
               allowed the 'time' input to be passed as a Matlab time string. \\ 
   22/08/2005  Following a bug report from Bruce Bowler
               (bbowler@bigelow.org), modified the $julian_calculation$ function. Bug
               was 'MATLAB has allowed structure assignment  to a non-empty non-structure 
               to overwrite the previous value.  This behavior will continue in this release, 
               but will be an error in a  future version of MATLAB.  For advice on how to 
               write code that  will both avoid this warning and work in future versions of 
               MATLAB,  see R14SP2 Release Notes'. Script should now be
               compliant with futher release of Matlab...

#endif



!EOP

!      Other variables
       REAL(dp)                 :: JulD, fraction2day  !, ntime_day
       integer, DIMENSION(1:6)  :: timevec       ! (/Year, Month, Day, Hour, Minute, Seconds/)
       REAL(dp)                 :: nutation_long, nutation_obliq
       REAL(dp)                 :: true_obliquity, aberration_correction
       REAL(dp)                 :: apparent_sun_longitude, JD, JC
       REAL(dp)                 :: mean_stime, apparent_stime_at_greenwich
       REAL(dp)                 :: argument_numerator, argument_denominator
       REAL(dp)                 :: sun_rigth_ascension
       REAL(dp)                 :: argument, sun_geocentric_declination
       REAL(dp)                 :: observer_local_hour, topocentric_local_hour

       type(julian_cal)         :: julian
       type(earth_heliocen_pos) :: earth_heliocentric_position
       type(sun_geocentric_pos) :: sun_geocentric_position
       type(topocentric_sun_pos):: topocentric_sun_position
!      ......................


       CALL get_timevec(nymd, nhms, fraction2day,  timevec)
!      calculated at UTC ... not local time
       CALL JDay_improved ( timevec(1),timevec(2), timevec(3), timevec(4), timevec(5), timevec(6), julian ) 
       CALL calc_earth_heliocentric_pos(julian, earth_heliocentric_position)
       CALL calc_sun_geocentric_position(sun_geocentric_position, earth_heliocentric_position)
!      4. Calculate the nutation in longitude and obliquity (in degrees).
       CALL calc_nutation(nutation_long, nutation_obliq, julian) ! degree
!      5. Calculate the true obliquity of the ecliptic (in degrees)
       CALL calc_true_obliquity(true_obliquity, julian, nutation_long, nutation_obliq) 
!      6. Calculate the aberration correction (in degrees)
!      This function compute the aberration_correction, as a function of the
!      earth-sun distance. 
       aberration_correction = -20.4898/(3600.0 * earth_heliocentric_position%radius)
!      7. (in degrees)
       apparent_sun_longitude = sun_geocentric_position%longitude + nutation_long + aberration_correction
!      8. compute the apparent sideral time at Greenwich (in degrees)
!      Mean sideral time, in degrees
       JD = julian%day
       JC = julian%century
       mean_stime = 280.46061837 + (360.98564736629*(JD-2451545.0)) + (0.000387933*JC**2) - (JC**3.0/38710000.0)
!      Limit the range to [0-360]
       CALL set_to_range(mean_stime, 0.0, 360.0)
       apparent_stime_at_greenwich = mean_stime + (nutation_long * cos(true_obliquity * pi/180.0))

!      9. Calculate the sun rigth ascension (in degrees)
       argument_numerator = (sin(apparent_sun_longitude * pi/180.0) * cos(true_obliquity * pi/180.0)) -  &
       (tan(sun_geocentric_position%latitude * pi/180.0) * sin(true_obliquity * pi/180.0))
       argument_denominator = cos(apparent_sun_longitude * pi/180.0)
! 
       sun_rigth_ascension = atan2(argument_numerator, argument_denominator) * 180.0/pi
!      Limit the range to [0,360];
       CALL set_to_range(sun_rigth_ascension, 0.0, 360.0)

!      10. Calculate the geocentric sun declination (in degrees). Positive or
!      negative if the sun is north or south of the celestial equator. 
       argument = (sin(sun_geocentric_position%latitude * pi/180.0) *     &
                  cos(true_obliquity * pi/180.0)) +                       &
                  (cos(sun_geocentric_position%latitude * pi/180.0) *     &
                  sin(true_obliquity * pi/180.0) * sin(apparent_sun_longitude * pi/180.0))
       sun_geocentric_declination = asin(argument) * 180.0/pi

!      11. Calculate the observer local hour angle (in degrees, westward from south)
       observer_local_hour = apparent_stime_at_greenwich + long - sun_rigth_ascension
!      Set the range to [0-360]
       CALL set_to_range(observer_local_hour, 0.0, 360.0)




!       12. Calculate the topocentric sun position (rigth ascension, declination and
!       rigth ascension parallax in degrees) with respect to the observer local position 
!       at the Earth surface
       CALL calc_topocentric_sun_position(topocentric_sun_position,         &
                                          earth_heliocentric_position,      &
                   lat, long, alt, observer_local_hour,                     &
                   sun_rigth_ascension, sun_geocentric_declination)

!      13. Calculate the topocentric local hour angle (in degrees)
       topocentric_local_hour = observer_local_hour - topocentric_sun_position%rigth_ascension_parallax



       CALL calc_suntopocen_zenith_angle(sun_az, sun_zen, lat, long,      &
            alt, topocentric_sun_position, topocentric_local_hour)

print*, "in side function == ", sun_az, sun_zen


END SUBROUTINE CALC_sun_zen_az





SUBROUTINE calc_suntopocen_zenith_angle(sun_az, sun_zen, lat, long,      &
            alt, topocentric_sun_position, topocentric_local_hour)

    implicit none



       REAL(dp), intent(out) :: sun_az, sun_zen
       real(dp), intent(in)  :: long, lat, alt, topocentric_local_hour
       Type(topocentric_sun_pos), intent(inout) :: topocentric_sun_position

!      Other variables
       REAL(dp)  :: argument, true_elevation, refraction_corr
       REAL(dp)  :: apparent_elevation
       REAL(dp)  :: nominator, denominator

!      This function compute the sun zenith angle, taking into account the
!      atmospheric refraction. A default temperature of 283K and a
!      default pressure of 1010 mbar are used. 

!      Topocentric elevation, without atmospheric refraction
       argument = (sin(lat * pi/180.0) *                                      &
                  sin(topocentric_sun_position%declination * pi/180.0)) +     &
                  (cos(lat * pi/180.0) * cos(topocentric_sun_position%declination * pi/180.0)  & 
                  * cos(topocentric_local_hour * pi/180.0))
       true_elevation = asin(argument) * 180.0/pi

!      Atmospheric refraction correction (in degrees)
       argument = true_elevation + (10.3/(true_elevation + 5.11));
       refraction_corr = 1.02 / (60.0 * tan(argument * pi/180.0));





!      For exact pressure and temperature correction, use this, 
!      with P the pressure in mbar amd T the temperature in Kelvins:
!      refraction_corr = (P/1010) * (283/T) * 1.02 / (60 * tan(argument * pi/180));

!      Apparent elevation
       if(true_elevation > -5.0) then
          apparent_elevation = true_elevation + refraction_corr
       else
          apparent_elevation = true_elevation
       end if

       sun_zen = 90.0 - apparent_elevation

!      Topocentric azimuth angle. The +180 conversion is to pass from astronomer
!      notation (westward from south) to navigation notation (eastward from
!      north);
       nominator = sin(topocentric_local_hour * pi/180.0)
       denominator = (cos(topocentric_local_hour * pi/180.0)      &
                     * sin(lat * pi/180.0)) -       & 
            (tan(topocentric_sun_position%declination * pi/180.0) & 
            * cos(lat * pi/180.0))
       sun_az = (atan2(nominator, denominator) * 180.0/pi) + 180.0
!      Set the range to [0-360]
       CALL set_to_range(sun_az, 0.0, 360.0)

END SUBROUTINE calc_suntopocen_zenith_angle




SUBROUTINE calc_topocentric_sun_position(this, earth_heliocentric_position, &
                   lat, long, alt, observer_local_hour,                     &
                   sun_rigth_ascension, sun_geocentric_declination)

    implicit none


       Type(topocentric_sun_pos), intent(inout) :: this
       !Type(sun_geocentric_pos), intent(in) :: sun_geocentric_position
       Type(earth_heliocen_pos), intent(in) :: earth_heliocentric_position
       real(dp), intent(in) :: long, lat, alt
       REAL(dp), intent(in) :: sun_rigth_ascension, sun_geocentric_declination
       REAL(dp), intent(in) :: observer_local_hour

!      Other variables
       REAL(dp) :: eq_horizontal_parallax, u, x, y
       REAL(dp) :: nominator, denominator, sun_rigth_ascension_parallax


!      Equatorial horizontal parallax of the sun in degrees
       eq_horizontal_parallax = 8.794 / (3600.0 * earth_heliocentric_position%radius)
!      Term u, used in the following calculations (in radians)
       u = atan(0.9966479 * tan(lat * pi/180.0))
!      Term x, used in the following calculations
       x = cos(u) + ((alt/6378140.0) * cos(lat * pi/180.0))

!      Term y, used in the following calculations
       y = (0.99664719 * sin(u)) + ((alt/6378140.0) * sin(lat * pi/180.0))

!      Parallax in the sun rigth ascension (in radians)
       nominator = -x * sin(eq_horizontal_parallax * pi/180.0) *       &
                   sin(observer_local_hour * pi/180.0)
       denominator = cos(sun_geocentric_declination * pi/180.0) -      &
                     (x * sin(eq_horizontal_parallax * pi/180.0) *     &
                     cos(observer_local_hour * pi/180.0))
       sun_rigth_ascension_parallax = atan2(nominator, denominator)
!      Conversion to degrees. 
       this%rigth_ascension_parallax = sun_rigth_ascension_parallax * 180.0/pi
!      Topocentric sun rigth ascension (in degrees)
       this%rigth_ascension = sun_rigth_ascension + (sun_rigth_ascension_parallax * 180.0/pi)
!      Topocentric sun declination (in degrees)
       nominator = (sin(sun_geocentric_declination * pi/180.0) -      &
                   (y*sin(eq_horizontal_parallax * pi/180.0))) *      &
                   cos(sun_rigth_ascension_parallax);
       denominator = cos(sun_geocentric_declination * pi/180.0) -     &
                   (x*sin(eq_horizontal_parallax * pi/180.0)) *       &
                   cos(observer_local_hour * pi/180.0)
       this%declination = atan2(nominator, denominator) * 180.0/pi
END SUBROUTINE calc_topocentric_sun_position

SUBROUTINE calc_true_obliquity(true_obliquity, julian, nutation_long, nutation_obliq)

    implicit none

       REAL(dp), intent(out)          :: true_obliquity
       REAL(dp), intent(in)           :: nutation_long, nutation_obliq
       type(julian_cal), intent(in)   :: julian

!      Other variables
       REAL(dp), DIMENSION(1:11)        :: p
       REAL(dp)                         :: U, mean_obliquity


       DATA p/2.45, 5.79, 27.87, 7.12, -39.05, -249.67, -51.38, 1999.25, &
              -1.55, -4680.93, 84381.448/
!      mean_obliquity = polyval(p, julian.ephemeris_millenium/10);
       U = julian%ephemeris_millenium/10.0;
       mean_obliquity = p(1)*U**10.0 + p(2)*U**9.0 + p(3)*U**8.0      &
                      + p(4)*U**7.0 + p(5)*U**6.0 + p(6)*U**5.0 +     &
                        p(7)*U**4.0 + p(8)*U**3.0 + p(9)*U**2.0 + p(10)*U + p(11)
       true_obliquity = (mean_obliquity/3600) + nutation_obliq

END SUBROUTINE calc_true_obliquity




SUBROUTINE calc_nutation(nutation_long, nutation_obliq, julian)
    implicit none


       REAL(dp), intent(out)           :: nutation_long, nutation_obliq
       type(julian_cal), intent(in)    :: julian

!      Other variables
       INTEGER, DIMENSION(1:63, 1:5)   :: Y_terms
       REAL(dp), DIMENSION(1:63, 1:4)  :: nutation_terms
       REAL(dp), DIMENSION(1:63, 1:1)  :: tab_arg
       REAL(dp), DIMENSION(1:63)       :: delta_longitude, delta_obliquity
       REAL(dp), DIMENSION(1:5,1)      :: Xi
       REAL(dp), PARAMETER             :: temp1 = 1.0/189474.0 
       REAL(dp), PARAMETER             :: temp2 = -(1.0/300000.0)
       REAL(dp), PARAMETER             :: temp3 = -(1.0/56250.0)
       REAL(dp), PARAMETER             :: temp4 =  (1.0/327270.0)
       REAL(dp), PARAMETER             :: temp5 =  (1.0/450000.0)
       REAL(dp)                        :: JCE, X0, X1, X2, X3, X4
       REAL(dp), DIMENSION(1:4)        :: p, p1, p2, p3, p4
       INTEGER                         :: i,j



!      This function compute the nutation in longtitude and in obliquity, in
!      degrees. 

!      All Xi are in degrees. 
       JCE = julian%ephemeris_century
!      1. Mean elongation of the moon from the sun 
       !temp = 1.0/189474.0
       DATA p/temp1, -0.0019142, 445267.11148, 297.85036/
!      X0 = polyval(p, JCE);
       X0 = p(1) * JCE**3.0 + p(2) * JCE**2.0 + p(3) * JCE + p(4) ! This is faster than polyval...

!      2. Mean anomaly of the sun (earth)
       DATA p1/temp2, -0.0001603, 35999.05034, 357.52772/
!      X1 = polyval(p, JCE);
       X1 = p1(1) * JCE**3.0 + p1(2) * JCE**2.0 + p1(3) * JCE + p1(4) 

!      3. Mean anomaly of the moon
       DATA p2/temp3, 0.0086972, 477198.867398, 134.96298/
!      X2 = polyval(p, JCE);
       X2 = p2(1) * JCE**3.0 + p2(2) * JCE**2.0 + p2(3) * JCE + p2(4) 

!      4. Moon argument of latitude
       DATA p3/temp4, -0.0036825, 483202.017538, 93.27191/
!      X3 = polyval(p, JCE);
       X3 = p3(1) * JCE**3.0 + p3(2) * JCE**2.0 + p3(3) * JCE + p3(4) 

!      5. Longitude of the ascending node of the moon's mean orbit on the
!      ecliptic, measured from the mean equinox of the date
       DATA p4/temp5, 0.0020708, -1934.136261, 125.04452/
!      X4 = polyval(p, JCE);
       X4 = p4(1) * JCE**3.0 + p4(2) * JCE**2.0 + p4(3) * JCE + p4(4) 

!      Y tabulated terms from the original code
       DATA (( Y_terms(i,j), j=1,5), i=1,63) &      
             / 0, 0, 0, 0, 1, -2, 0, 0, 2, 2, 0, 0, 0, 2, 2,   & 
               0, 0, 0, 0, 2,  0, 1, 0, 0, 0, 0, 0, 1, 0, 0,   &
               -2, 1, 0, 2, 2, 0, 0, 0, 2, 1, 0, 0, 1, 2, 2,   &
               -2, -1, 0, 2, 2,-2, 0, 1, 0, 0,-2, 0, 0, 2, 1,  &
               0, 0, -1, 2, 2, 2, 0, 0, 0, 0,  0, 0, 1, 0, 1,  &
               2, 0, -1, 2, 2,  0, 0, -1, 0, 1, 0, 0, 1, 2, 1, & 
              -2, 0, 2, 0, 0,  0, 0, -2, 2, 1, 2, 0, 0, 2, 2,  &
               0, 0, 2, 2, 2,   0, 0, 2, 0, 0, -2, 0, 1, 2, 2, & 
               0, 0, 0, 2, 0,  -2, 0, 0, 2, 0, 0, 0, -1, 2, 1, & 
               0, 2, 0, 0, 0,  2, 0, -1, 0, 1, -2, 2, 0, 2, 2, & 
               0, 1, 0, 0, 1, -2, 0, 1, 0, 1,  0, -1, 0, 0, 1, & 
               0, 0, 2, -2, 0, 2, 0, -1, 2, 1, 2, 0, 1, 2, 2,  &
               0, 1, 0, 2, 2, -2, 1, 1, 0, 0, 0, -1, 0, 2, 2,  &
               2, 0, 0, 2, 1,  2, 0, 1, 0, 0, -2, 0, 2, 2, 2,  &
               -2, 0, 1, 2, 1, 2, 0, -2, 0, 1, 2, 0, 0, 0, 1,  &
               0, -1, 1, 0, 0,-2, -1, 0, 2, 1,-2, 0, 0, 0, 1,  &
               0, 0, 2, 2, 1, -2, 0, 2, 0, 1, -2, 1, 0, 2, 1,  &
               0, 0, 1, -2, 0,-1, 0, 1, 0, 0, -2, 1, 0, 0, 0,  &
               1, 0, 0, 0, 0,  0, 0, 1, 2, 0, 0, 0, -2, 2, 2,  & 
               -1, -1, 1, 0, 0, 0, 1, 1, 0, 0,0, -1, 1, 2, 2,  &
               2, -1, -1, 2, 2, 0, 0, 3, 2, 2,2, -1, 0, 2, 2/


       DATA (( nutation_terms(i,j), j=1,4), i=1,63) &      
      /-171996.0, -174.2, 92025.0, 8.9, -13187.0, -1.6, 5736.0, -3.1, -2274.0, -0.2, 977.0, -0.5, & 
        2062.0, 0.2, -895.0, 0.5, 1426.0, -3.4, 54.0, -0.1, 712.0, 0.1, -7.0, 0.0,                &
        -517.0, 1.2, 224.0, -0.6, -386.0, -0.4, 200.0, 0.0, -301.0, 0.0, 129.0, -0.1,             &
        217.0, -0.5, -95.0, 0.3, -158.0, 0.0, 0.0, 0.0, 129.0, 0.1, -70.0, 0.0,                   &
        123.0, 0.0, -53.0, 0.0, 63.0, 0.0, 0.0, 0.0, 63.0, 0.1, -33.0, 0.0,                       &
        -59.0, 0.0, 26.0, 0.0, -58.0, -0.1, 32.0, 0.0, -51.0, 0.0, 27.0, 0.0,                     &
        48.0, 0.0, 0.0, 0.0, 46.0, 0.0, -24.0, 0.0, -38.0, 0.0, 16.0, 0.0,                        &
        -31.0, 0.0, 13.0, 0.0, 29.0, 0.0, 0.0, 0.0, 29.0, 0.0, -12.0, 0.0,                        &
        26.0, 0.0, 0.0, 0.0, -22.0, 0.0, 0.0, 0.0, 21.0, 0.0, -10.0, 0.0,                         &
        17.0, -0.1, 0.0, 0.0, 16.0, 0.0, -8.0, 0.0, -16.0, 0.1, 7.0, 0.0,                         &
        -15.0, 0.0, 9.0, 0.0, -13.0, 0.0, 7.0, 0.0, -12.0, 0.0, 6.0, 0.0,                         &
        11.0, 0.0, 0.0, 0.0, -10.0, 0.0, 5.0, 0.0,  -8.0, 0.0, 3.0, 0.0,                          &
        7.0, 0.0, -3.0, 0.0, -7.0, 0.0, 0.0, 0.0, -7.0, 0.0, 3.0, 0.0,                            &
        -7.0, 0.0, 3.0, 0.0, 6.0, 0.0, 0.0, 0.0 , 6.0, 0.0, -3.0, 0.0 ,                           &
        6.0, 0.0, -3.0, 0.0, -6.0, 0.0, 3.0, 0.0, -6.0, 0.0, 3.0, 0.0,                            &
        5.0, 0.0, 0.0, 0.0 , -5.0, 0.0, 3.0, 0.0, -5.0, 0.0, 3.0, 0.0,                            &
        -5.0, 0.0, 3.0, 0.0,  4.0, 0.0, 0.0, 0.0,  4.0, 0.0, 0.0, 0.0,                            &
        4.0, 0.0, 0.0, 0.0,  -4.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0,                            &
        -4.0, 0.0, 0.0, 0.0,  3.0, 0.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0,                            &
        -3.0, 0.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0,                            &
        -3.0, 0.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0/



!      Using the tabulated values, compute the delta_longitude and
!      delta_obliquity. 
       Xi(1,1) =   X0
       Xi(2,1) =   X1
       Xi(3,1) =   X2
       Xi(4,1) =   X3
       Xi(5,1) =   X4


       tab_arg = matmul( real(Y_terms,dp), Xi) * pi/180.0        !tabulated_argument
       delta_longitude(1:63) = ((nutation_terms(1:63,1) + (nutation_terms(1:63,2) * JCE))) 
       !delta_longitude(1:63) = ((nutation_terms(1:63,1) + (nutation_terms(1:63,2) * JCE))) * sin(tab_arg(1:63,1)
       DO i = 1,63
         delta_longitude(i) = delta_longitude(i) * sin(tab_arg(i,1))
       END DO
     

       delta_obliquity(1:63) = ((nutation_terms(1:63,3) + (nutation_terms(1:63,4) * JCE))) 
       DO i = 1,63
         delta_obliquity(i) = delta_obliquity(i) * cos(tab_arg(i,1))
       END DO
!      Nutation in longitude
       nutation_long = sum(delta_longitude) / 36000000.0

!      Nutation in obliquity
       nutation_obliq = sum(delta_obliquity) / 36000000.0

END SUBROUTINE calc_nutation





SUBROUTINE calc_sun_geocentric_position(this, earth_heliocentric_position)
    implicit none


       Type(sun_geocentric_pos), intent(inout) :: this
       Type(earth_heliocen_pos), intent(in) :: earth_heliocentric_position


       this%longitude = earth_heliocentric_position%longitude + 180.0
!      Limit the range to [0,360];
       CALL set_to_range(this%longitude, 0.0, 360.0)

       this%latitude = -earth_heliocentric_position%latitude
!      Limit the range to [0,360]
       CALL set_to_range(this%latitude, 0.0, 360.0)

END SUBROUTINE calc_sun_geocentric_position

 SUBROUTINE calc_earth_heliocentric_pos(julian, this)
! This function compute the earth position relative to the sun, using
! tabulated values. 

! Tabulated values for the longitude calculation
! L terms  from the original code.
       Type(julian_cal), intent(in) :: julian
       Type(earth_heliocen_pos), intent(inout) :: this

!      Other variables 


      REAL(dp), DIMENSION(:), ALLOCATABLE  :: A0, B0, C0, A1, B1, C1, A2, B2, C2
      REAL(dp), DIMENSION(:), ALLOCATABLE  :: A3, B3, C3, A4, B4, C4



       REAL(dp), DIMENSION(1:64, 1:3)  :: L0_terms
       REAL(dp), DIMENSION(1:34, 1:3)  :: L1_terms
       REAL(dp), DIMENSION(1:20, 1:3)  :: L2_terms
       REAL(dp), DIMENSION(1:7,  1:3)  :: L3_terms
       REAL(dp), DIMENSION(1:3,  1:3)  :: L4_terms
       REAL(dp), DIMENSION(1:1,  1:3)  :: L5_terms   
       REAL(dp), DIMENSION(1:5,  1:3)  :: B0_terms
       REAL(dp), DIMENSION(1:2,  1:3)  :: B1_terms


       REAL(dp), DIMENSION(1:40,  1:3) :: R0_terms
       REAL(dp), DIMENSION(1:10,  1:3) :: R1_terms
       REAL(dp), DIMENSION(1:6,  1:3)  :: R2_terms
       REAL(dp), DIMENSION(1:2,  1:3)  :: R3_terms
       REAL(dp), DIMENSION(1:1,  1:3)  :: R4_terms




!        REAL(dp), DIMENSION(1:SIZE(L0_terms,1)) :: A0, B0, C0
!        REAL(dp), DIMENSION(1:SIZE(L1_terms,1)) :: A1, B1, C1
!        REAL(dp), DIMENSION(1:SIZE(L2_terms,1)) :: A2, B2, C2
!        REAL(dp), DIMENSION(1:SIZE(L3_terms,1)) :: A3, B3, C3
!        REAL(dp), DIMENSION(1:SIZE(L4_terms,1)) :: A4, B4, C4
        REAL(dp)                                :: A5, B5, C5
       REAL(dp)                                :: L0, L1, L2, L3, L4, L5, JME
       INTEGER                    :: i,j
 

       DATA (( L0_terms(i,j), j=1,3), i=1,64) &      
                     /175347046.0, 0.0, 0.0, 3341656.0, 4.6692568, 6283.07585,  &  
                     34894.0, 4.6261, 12566.1517, 3497.0, 2.7441, 5753.3849,    &  
                     3418.0, 2.8289, 3.5231, 3136.0, 3.6277, 77713.7715,        &  
                     2676.0, 4.4181, 7860.4194, 2343.0, 6.1352, 3930.2097,      &  
                     1324.0, 0.7425, 11506.7698, 1273.0, 2.0371, 529.691,       &
                     1199.0, 1.1096, 1577.3435, 990.0, 5.233, 5884.927,         &
                     902.0, 2.045, 26.298, 857.0, 3.508, 398.149,               &
                     780.0, 1.179, 5223.694, 753.0, 2.533, 5507.553,            &
                     505.0, 4.583, 18849.228, 492.0, 4.205, 775.523,            &
                     357.0, 2.92, 0.067, 317.0, 5.849, 11790.629,               &
                     284.0, 1.899, 796.298, 271.0, 0.315, 10977.079,            &
                     243.0, 0.345, 5486.778, 206.0, 4.806, 2544.314,            &
                     205.0, 1.869, 5573.143, 202.0, 2.4458, 6069.777,           &
                     156.0, 0.833, 213.299, 132.0, 3.411, 2942.463,             &
                     126.0, 1.083, 20.775, 115.0, 0.645, 0.98,                  &
                     103.0, 0.636, 4694.003, 102.0, 0.976, 15720.839,           &
                     102.0, 4.267, 7.114, 99.0, 6.21, 2146.17,                  &
                     98.0, 0.68, 155.42,  86.0, 5.98, 161000.69,                &
                     85.0, 1.3, 6275.96, 85.0, 3.67, 71430.7,                   &
                     80.0, 1.81, 17260.15, 79.0, 3.04, 12036.46,                &
                     71.0, 1.76, 5088.63 , 74.0, 3.5, 3154.69,                  &
                     74.0, 4.68, 801.82, 70.0, 0.83, 9437.76,                   &
                     62.0, 3.98, 8827.39, 61.0, 1.82, 7084.9,                   &
                     57.0, 2.78, 6286.6,  56.0, 4.39, 14143.5,                  &
                     56.0, 3.47, 6279.55, 52.0, 0.19, 12139.55,                 &
                     52.0, 1.33, 1748.02, 51.0, 0.28, 5856.48,                  &
                     49.0, 0.49, 1194.45, 41.0, 5.37, 8429.24,                  &
                     41.0, 2.4, 19651.05, 39.0, 6.17, 10447.39,                 &
                     37.0, 6.04, 10213.29, 37.0, 2.57, 1059.38,                 &
                     36.0, 1.71, 2352.87, 36.0, 1.78, 6812.77,                  &
                     33.0, 0.59, 17789.85, 30.0, 0.44, 83996.85,                &
                     30.0, 2.74, 1349.87, 25.0, 3.16, 4690.48/


       DATA (( L1_terms(i,j), j=1,3), i=1,34) &    
                     /628331966747.0, 0.0, 0.0, 206059.0, 2.678235, 6283.07585, &  
                      4303.0, 2.6351, 12566.1517, 425.0, 1.59, 3.523,           &  
                      119.0, 5.796, 26.298, 109.0, 2.966, 1577.344,             &  
                      93.0, 2.59, 18849.23, 72.0, 1.14, 529.69,                 &
                      68.0, 1.87, 398.15,  67.0, 4.41, 5507.55,                 &
                      59.0, 2.89, 5223.69, 56.0, 2.17, 155.42,                  &
                      45.0, 0.4, 796.3,  36.0, 0.47, 775.52,                    &
                      29.0, 2.65, 7.11, 21.0, 5.34, 0.98,                       &
                      19.0, 1.85, 5486.78,  19.0, 4.97, 213.3,                  &
                      17.0, 2.99, 6275.96, 16.0, 0.03, 2544.31,                 &
                      16.0, 1.43, 2146.17,  15.0, 1.21, 10977.08,               &
                      12.0, 2.83, 1748.02,                                      &
                      12.0, 3.26, 5088.63, 12.0, 5.27, 1194.45,                 &
                      12.0, 2.08, 4694.0,  11.0, 0.77, 553.57,                  &
                      10.0, 1.3, 3286.6,  10.0, 4.24, 1349.87,                  &
                      9.0, 2.7, 242.73,  9.0, 5.64, 951.72,                     &
                      8.0, 5.3, 2352.87,  6.0, 2.65, 9437.76,                   &
                      6.0, 4.67, 4690.48/ 


       DATA (( L2_terms(i,j), j=1,3), i=1,20) &    
                     /52919.0, 0.0, 0.0, 8720.0, 1.0721, 6283.0758,             &
                      309.0, 0.867, 12566.152, 27.0, 0.05, 3.52,                &
                      16.0, 5.19, 26.3, 16.0, 3.68, 155.42,                     &
                      10.0, 0.76, 18849.23,  9.0, 2.06, 77713.77,               &
                      7.0, 0.83, 775.52,  5.0, 4.66, 1577.34,                   &
                      4.0, 1.03, 7.11, 4.0, 3.44, 5573.14,                      &
                      3.0, 5.14, 796.3,  3.0, 6.05, 5507.55,                    &
                      3.0, 1.19, 242.73, 3.0, 6.12, 529.69,                     &
                      3.0, 0.31, 398.15, 3.0, 2.28, 553.57,                     &
                      2.0, 4.38, 5223.69, 2.0, 3.75, 0.98/ 


       DATA (( L3_terms(i,j), j=1,3), i=1,7) &    
                     /289.0, 5.844, 6283.076, 35.0, 0.0, 0.0,                   &  
                      17.0, 5.49, 12566.15, 3.0, 5.2, 155.42,                   &
                      1.0, 4.72, 3.52, 1.0, 5.3, 18849.23,                      &
                      1.0, 5.97, 242.73/ 

      DATA (( L4_terms(i,j), j=1,3), i=1,3) &    
                     /114.0, 3.142, 0.0, 8.0, 4.13, 6283.08,                    & 
                      1.0, 3.84, 12566.15/


      DATA (( L5_terms(i,j), j=1,3), i=1,1) &    
                     /1.0, 3.14, 0.0/


!        print*, size(L0_terms,1)
! 
!        print*, L2_terms(1,1), L2_terms(1,2), L2_terms(1,3)
!        print*, L2_terms(2,1), L2_terms(2,2), L2_terms(2,3)
!        print*, L0_terms(3,1), L0_terms(3,2), L0_terms(3,3)
!        print*, L0_terms(64,1), L0_terms(64,2), L0_terms(64,3)


       Allocate ( A0(1:SIZE(L0_terms,1)), B0(1:SIZE(L0_terms,1)), C0(1:SIZE(L0_terms,1)) )
       Allocate ( A1(1:SIZE(L1_terms,1)), B1(1:SIZE(L1_terms,1)), C1(1:SIZE(L1_terms,1)) )
       Allocate ( A2(1:SIZE(L2_terms,1)), B2(1:SIZE(L2_terms,1)), C2(1:SIZE(L2_terms,1)) )
       Allocate ( A3(1:SIZE(L3_terms,1)), B3(1:SIZE(L3_terms,1)), C3(1:SIZE(L3_terms,1)) )
       Allocate ( A4(1:SIZE(L4_terms,1)), B4(1:SIZE(L4_terms,1)), C4(1:SIZE(L4_terms,1)) )



       A0(1:SIZE(L0_terms,1)) = L0_terms(1:SIZE(L0_terms,1),1)
       B0(1:SIZE(L0_terms,1)) = L0_terms(1:SIZE(L0_terms,1),2)
       C0(1:SIZE(L0_terms,1)) = L0_terms(1:SIZE(L0_terms,1),3)
       A1(1:SIZE(L1_terms,1)) = L1_terms(1:SIZE(L1_terms,1),1);
       B1(1:SIZE(L1_terms,1)) = L1_terms(1:SIZE(L1_terms,1),2);
       C1(1:SIZE(L1_terms,1)) = L1_terms(1:SIZE(L1_terms,1),3);
       A2(1:SIZE(L2_terms,1)) = L2_terms(1:SIZE(L2_terms,1),1);
       B2(1:SIZE(L2_terms,1)) = L2_terms(1:SIZE(L2_terms,1),2);
       C2(1:SIZE(L2_terms,1)) = L2_terms(1:SIZE(L2_terms,1),3);
       A3(1:SIZE(L3_terms,1)) = L3_terms(1:SIZE(L3_terms,1),1);
       B3(1:SIZE(L3_terms,1)) = L3_terms(1:SIZE(L3_terms,1),2);
       C3(1:SIZE(L3_terms,1)) = L3_terms(1:SIZE(L3_terms,1),3);
       A4(1:SIZE(L4_terms,1)) = L4_terms(1:SIZE(L4_terms,1),1);
       B4(1:SIZE(L4_terms,1)) = L4_terms(1:SIZE(L4_terms,1),2);
       C4(1:SIZE(L4_terms,1)) = L4_terms(1:SIZE(L4_terms,1),3);
       A5 = L5_terms(1,1);
       B5 = L5_terms(1,2);
       C5 = L5_terms(1,3);

       JME = julian%ephemeris_millenium

!      Compute the Earth Heliochentric longitude from the tabulated values. 


       !L0 = C0 * JME
       L0 = sum(A0 * cos(B0 + (C0 * JME)))
       L1 = sum(A1 * cos(B1 + (C1 * JME)))
       L2 = sum(A2 * cos(B2 + (C2 * JME)))
       L3 = sum(A3 * cos(B3 + (C3 * JME)))
       L4 = sum(A4 * cos(B4 + (C4 * JME)))
       L5 = A5 * cos(B5 + (C5 * JME))

       !print*, L0, L1, L2, L3, L4, L5

       this%longitude = (L0 + (L1 * JME) + (L2 * JME**2.0) + (L3 * JME**3.0) + (L4 * JME**4.0) + (L5 * JME**5.0)) / 1e8
!      Convert the longitude to degrees. 
       this%longitude = this%longitude * 180/pi
!      Limit the range to [0,360[;
       CALL set_to_range(this%longitude, 0.0, 360.0)
       deallocate(A0, B0, C0, A1, B1, C1, A2, B2, C2, A3, B3, C3, A4, B4, C4)

!      *************************************************************


       DATA (( B0_terms(i,j), j=1,3), i=1,5) &      
                     /280.0, 3.199, 84334.662, 102.0, 5.422, 5507.553,     &
                      80.0, 3.88, 5223.69, 44.0, 3.7, 2352.87,             & 
                      32.0, 4.0, 1577.34/

       DATA (( B1_terms(i,j), j=1,3), i=1,2) &      
                     /9.0, 3.9, 5507.55, 6.0, 1.73, 5223.69/


       Allocate ( A0(1:SIZE(B0_terms,1)), B0(1:SIZE(B0_terms,1)), C0(1:SIZE(B0_terms,1)) )
       Allocate ( A1(1:SIZE(B1_terms,1)), B1(1:SIZE(B1_terms,1)), C1(1:SIZE(B1_terms,1)) )


       A0(1:SIZE(B0_terms,1)) = B0_terms(1:SIZE(B0_terms,1),1)
       B0(1:SIZE(B0_terms,1)) = B0_terms(1:SIZE(B0_terms,1),2)
       C0(1:SIZE(B0_terms,1)) = B0_terms(1:SIZE(B0_terms,1),3)
       A1(1:SIZE(B1_terms,1)) = B1_terms(1:SIZE(B1_terms,1),1)
       B1(1:SIZE(B1_terms,1)) = B1_terms(1:SIZE(B1_terms,1),2)
       C1(1:SIZE(B1_terms,1)) = B1_terms(1:SIZE(B1_terms,1),3)

       L0 = sum(A0 * cos(B0 + (C0 * JME)))
       L1 = sum(A1 * cos(B1 + (C1 * JME)))


       this%latitude = (L0 + (L1 * JME)) / 1e8; 
!      Convert the latitude to degrees. 
       this%latitude = this%latitude * 180.0/pi;
!      Limit the range to [0,360];
       CALL set_to_range(this%latitude, 0.0, 360.0);

       deallocate (A0, B0, C0, A1, B1, C1)
!      *************************************************************


      DATA (( R0_terms(i,j), j=1,3), i=1,40) &  
                         /100013989.0, 0.0, 0.0, 1670700.0, 3.0984635, 6283.07585,  &  
                          13956.0, 3.05525, 12566.1517, 3084.0, 5.1985, 77713.7715, &  
                          1628.0, 1.1739, 5753.3849, 1576.0, 2.8469, 7860.4194,     &  
                          925.0, 5.453, 11506.77, 542.0, 4.564, 3930.21,            &  
                          472.0, 3.661, 5884.927,  346.0, 0.964, 5507.553,          & 
                          329.0, 5.9, 5223.694, 307.0, 0.299, 5573.143,             & 
                          243.0, 4.273, 11790.629, 212.0, 5.847, 1577.344,          &  
                          186.0, 5.022, 10977.079, 175.0, 3.012, 18849.228,         &  
                          110.0, 5.055, 5486.778, 98.0, 0.89, 6069.78,              &  
                          86.0, 5.69, 15720.84, 86.0, 1.27, 161000.69,              & 
                          85.0, 0.27, 17260.15, 63.0, 0.92, 529.69,                 &  
                          57.0, 2.01, 83996.85, 56.0, 5.24, 71430.7,                & 
                          49.0, 3.25, 2544.31, 47.0, 2.58, 775.52,                  &
                          45.0, 5.54, 9437.76, 43.0, 6.01, 6275.96,                 &  
                          39.0, 5.36, 4694.0,  38.0, 2.39, 8827.39,                 &
                          37.0, 0.83, 19651.05, 37.0, 4.9, 12139.55,                &
                          36.0, 1.67, 12036.46, 35.0, 1.84, 2942.46,                &
                          33.0, 0.24, 7084.9, 32.0, 0.18, 5088.63,                  &
                          32.0, 1.78, 398.15, 28.0, 1.21, 6286.6,                   &
                          28.0, 1.9, 6279.55, 26.0, 4.59, 10447.39/

      DATA (( R1_terms(i,j), j=1,3), i=1,10) &  
                         /103019.0, 1.10749, 6283.07585, 1721.0, 1.0644, 12566.1517,&  
                          702.0, 3.142, 0.0, 32.0, 1.02, 18849.23,                  & 
                          31.0, 2.84, 5507.55, 25.0, 1.32, 5223.69,                 &
                          18.0, 1.42, 1577.34, 10.0, 5.91, 10977.08,                &
                          9.0, 1.42, 6275.96,  9.0, 0.27, 5486.78/


      DATA (( R2_terms(i,j), j=1,3), i=1,6) &  
                         /4359.0, 5.7846, 6283.0758, 124.0, 5.579, 12566.152,       & 
                          12.0, 3.14, 0.0, 9.0, 3.63, 77713.77,                     &
                          6.0, 1.87, 5573.14, 3.0, 5.47, 18849.0/  


      DATA (( R3_terms(i,j), j=1,3), i=1,2) &  
                         /145.0, 4.273, 6283.076, 7.0, 3.92, 12566.15/

      DATA (( R4_terms(i,j), j=1,3), i=1,1) &  
                         /4.0, 2.56, 6283.08/

       Allocate ( A0(1:SIZE(R0_terms,1)), B0(1:SIZE(R0_terms,1)), C0(1:SIZE(R0_terms,1)) )
       Allocate ( A1(1:SIZE(R1_terms,1)), B1(1:SIZE(R1_terms,1)), C1(1:SIZE(R1_terms,1)) )
       Allocate ( A2(1:SIZE(R2_terms,1)), B2(1:SIZE(R2_terms,1)), C2(1:SIZE(R2_terms,1)) )
       Allocate ( A3(1:SIZE(R3_terms,1)), B3(1:SIZE(R3_terms,1)), C3(1:SIZE(R3_terms,1)) )
       Allocate ( A4(1:SIZE(R4_terms,1)), B4(1:SIZE(R4_terms,1)), C4(1:SIZE(R4_terms,1)) )


       A0(1:SIZE(R0_terms,1)) = R0_terms(1:SIZE(R0_terms,1),1)
       B0(1:SIZE(R0_terms,1)) = R0_terms(1:SIZE(R0_terms,1),2)
       C0(1:SIZE(R0_terms,1)) = R0_terms(1:SIZE(R0_terms,1),3)

       A1(1:SIZE(R1_terms,1)) = R1_terms(1:SIZE(R1_terms,1),1);
       B1(1:SIZE(R1_terms,1)) = R1_terms(1:SIZE(R1_terms,1),2);
       C1(1:SIZE(R1_terms,1)) = R1_terms(1:SIZE(R1_terms,1),3);
       A2(1:SIZE(R2_terms,1)) = R2_terms(1:SIZE(R2_terms,1),1);
       B2(1:SIZE(R2_terms,1)) = R2_terms(1:SIZE(R2_terms,1),2);
       C2(1:SIZE(R2_terms,1)) = R2_terms(1:SIZE(R2_terms,1),3);
       A3(1:SIZE(R3_terms,1)) = R3_terms(1:SIZE(R3_terms,1),1);
       B3(1:SIZE(R3_terms,1)) = R3_terms(1:SIZE(R3_terms,1),2);
       C3(1:SIZE(R3_terms,1)) = R3_terms(1:SIZE(R3_terms,1),3);
       A4(1:SIZE(R4_terms,1)) = R4_terms(1:SIZE(R4_terms,1),1);
       B4(1:SIZE(R4_terms,1)) = R4_terms(1:SIZE(R4_terms,1),2);
       C4(1:SIZE(R4_terms,1)) = R4_terms(1:SIZE(R4_terms,1),3);
!        A5 = L5_terms(1,1);
!        B5 = L5_terms(1,2);
!        C5 = L5_terms(1,3);
       L0 = sum(A0 * cos(B0 + (C0 * JME)))
       L1 = sum(A1 * cos(B1 + (C1 * JME)))
       L2 = sum(A2 * cos(B2 + (C2 * JME)))
       L3 = sum(A3 * cos(B3 + (C3 * JME)))
       L4 = A4(1) * cos(B4(1) + (C4(1) * JME))
 
!        print*, "--<> ", A4(1), B4(1), C4(1)
!        print*, "--<> ", L0, L1, L2, L3, L4

!       this%longitude = (L0 + (L1 * JME) + (L2 * JME**2.0) + (L3 * JME**3.0) + (L4 * JME**4.0) + (L5 * JME**5.0)) / 1e8
!      Convert the longitude to degrees. 
!       this%longitude = this%longitude * 180/pi
!      Limit the range to [0,360[;
!       CALL set_to_range(this%longitude, 0.0, 360.0)
!       print*, this%longitude
       deallocate(A0, B0, C0, A1, B1, C1, A2, B2, C2, A3, B3, C3, A4, B4, C4)
       this%radius = (L0 + (L1 * JME) + (L2 * JME**2.0) + (L3 * JME**3.0) + (L4 * JME**4.0)) / 1e8
!       print*, ";;", this%radius
 ENDSUBROUTINE calc_earth_heliocentric_pos



 SUBROUTINE set_to_range (var, min_interval, max_interval)
      implicit NONE

      REAL(dp), INTENT(INOUT) :: var
      REAL, INTENT(IN) :: min_interval, max_interval




! if(var>0)
!     var = var - max_interval * floor(var/max_interval);
! else
!     var = var - max_interval * ceil(var/max_interval);
! end
! 
! if(var<min_interval)
!     var = var + max_interval;
! end

       var = var - max_interval * floor(var/max_interval);
       if(var<min_interval) then
          var = var + max_interval
       end if

      return
END SUBROUTINE set_to_range



SUBROUTINE get_timevec(nymd, nhms, fraction2day, tt)
!SUBROUTINE get_time(iSat, nymd, nhms, fraction2day, ntime_day, tt)
        IMPLICIT  NONE
        INTEGER, INTENT(IN)     :: nymd, nhms
        REAL(dp), INTENT(OUT)   :: fraction2day  !, ntime_day
        INTEGER,  DIMENSION(1:6), INTENT(OUT)   :: tt
!       Other
        INTEGER :: Year, Month, Day, Hour, Minute, Seconds 

        Year   = int(nymd / 10000 )
        Month  = mod ( nymd,  10000 ) / 100
        Day    = mod ( nymd,    100 )
        Hour = int(nhms / 10000 )
        Minute = mod ( nhms,  10000 ) / 100
        Seconds = mod ( nhms,    100 )
        fraction2day = (Hour*60*60+Minute*60+Seconds)/(24.0*60.0*60.0)  !fraction
                                         !calculated from seconds  
        !ntime_day=(-ODS_Julian(RefDate(iSat))+ODS_Julian(nymd))+fraction2day
        tt = (/Year, Month, Day, Hour, Minute, Seconds/)
        

END SUBROUTINE get_timevec


SUBROUTINE JDay_improved ( Year,Mon,Day,Hr,Min, Sec, this )

        IMPLICIT NONE
        INTEGER, INTENT(IN)             :: Year, Mon, Day, Hr, Min, Sec
        Type(julian_cal), intent(inout) :: this
        real(dp)                        :: delta_t


        this%day = 367.0D0 * Year                                       &
             - INT( (7* (Year+INT ( (Mon+9)/12) ) ) * 0.25D0 )   &
             + INT( 275*Mon / 9 )                                &
             + Day + 1721013.5D0                                 &
             + ( (Sec/60.0D0 + Min ) / 60.0D0 + Hr ) / 24.0D0    
!            - 0.5D0*DSIGN(1.0D0, 100.0D0*Year + Mon - 190002.5D0) + 0.5D0
        delta_t = 0.0 ! 33.184
        this%ephemeris_day = this%day + (delta_t/86400.0)
        this%century = (this%day - 2451545.0) / 36525.0 
        this%ephemeris_century = (this%ephemeris_day - 2451545.0) / 36525.0
        this%ephemeris_millenium = this%ephemeris_century / 10.0 
END SUBROUTINE JDay_improved


integer function ODS_Julian ( CalDate )
      implicit NONE
      integer  CalDate  ! Calendar date in the format YYYYMMDD
                        !   where YYYY is the year, MM is the
                        !   month and DD is the day.  A negative
                        !   number implies that the year is B.C.
!     Other variables
!     ---------------
      integer     Year
      integer     Month
      integer     Day
      integer     iGreg  ! Gregorian Calendar adopted Oct 12, 1582
      parameter ( iGreg = 15 + 31 * ( 10 + 12 * 1582 ) )
      integer     JulDay
      integer     jy, jm, ja

      Year   =       CalDate / 10000
      Month  = mod ( CalDate,  10000 ) / 100
      Day    = mod ( CalDate,    100 )
!     Change year 0 to year 1
      if ( Year  .eq. 0 ) Year = 1
!     Account for the nonexisting year 0
      if ( Year  .lt. 0 ) Year = Year + 1
      if ( Month .gt. 2 ) then
         jy = Year
         jm = Month + 1
      else
         jy = Year  - 1
         jm = Month + 13
      endif
      JulDay = int ( 365.25  * jy )       &
             + int ( 30.6001 * jm )       &
             + Day + 1720995
!     Test whether to change to Gregorian Celendar
      if ( Day + 31 * ( Month + 12 * Year ) .ge. iGreg) then
        ja     = int ( 0.01 * jy )
        Julday = JulDay + 2 - ja + int ( 0.25 * ja )
      endif
      ODS_Julian = JulDay
      return
end function ODS_Julian


SUBROUTINE LLA2ECEF(lat, lon, alt, x, y, z)
!       % ECEF3LLA - convert earth-centered earth-fixed (ECEF)
!       %            cartesian coordinates to latitude, longitude,
!       %            and altitude
!       %
!       % USAGE:
!       % [lat,lon,alt] = ecef2lla(x,y,z)
!       %
!       % lat = geodetic latitude (radians)
!       % lon = longitude (radians)
!       % alt = height above WGS84 ellipsoid (m)
!       % x = ECEF X-coordinate (m)
!       % y = ECEF Y-coordinate (m)
!       % z = ECEF Z-coordinate (m)
!       %
!       % Notes: (1) This function assumes the WGS84 model.
!       %        (2) Latitude is customary geodetic (not geocentric).
!       %        (3) Inputs may be scalars, vectors, or matrices of the same
!       %            size and shape. Outputs will have that same size and shape.
!       %        (4) Tested but no warranty; use at your own risk.
!       %        (5) Michael Kleder, April 2006
       
       REAL(dp), INTENT(in) :: lat,lon,alt
       REAL(dp), INTENT(out)  :: x,y,z
!      % WGS84 ellipsoid constants:
       REAL(dp), PARAMETER  :: a = 6378137 ! those should be defined outside
       REAL(dp), PARAMETER  :: e = 8.1819190842622e-2
       REAL(dp), PARAMETER  :: pi = 3.14159265358979323846 
       REAL(dp)             :: b, ep, p, th, N
       LOGICAL k



!      % WGS84 ellipsoid constants:


       ! intermediate calculation
       ! (prime vertical radius of curvature)
       N = a / sqrt(1 - e**2 * sin(lat)**2)
       ! results:
       x = (N+alt) * cos(lat) * cos(lon)
       y = (N+alt) * cos(lat) * sin(lon)
       z = ((1-e**2) * N + alt) * sin(lat)


END SUBROUTINE LLA2ECEF


END MODULE glint_mod ! ........END........ 

! ------------------------------------------------------
! ------------------------------------------------------
! ------------------------------------------------------
! 
! 
! SUBROUTINE zeroTo360(x, unit, y)
!        IMPLICIT NONE
! 
!        INTEGER,                INTENT(IN)  :: unit
!        REAL(dp), INTENT(INOUT)  :: x
!        REAL(dp), INTENT(OUT) :: y
!        REAL(dp)   :: deg
!        INTEGER    :: j
! 
! !     ............................................................................
! 
!        if (unit.EQ.1) then
!            deg = 2.0*pi
!        else 
!            deg = 360.0
!        endif
!        if (x >= deg) then 
!           x = x - aint(x/deg)*deg
!        elseif (x < 0) then
!           x = x - (aint(x/deg) - 1)*deg
!        endif
!        y = x
! END SUBROUTINE zeroTo360
! 
! 
! SUBROUTINE Orbits_Track0(iSat, p)
!      IMPLICIT NONE
!  
!        !character(len=*), intent(in) :: Sat_name ! Satellite name
!        REAL(dp), DIMENSION(:), INTENT(out) :: p(num_satparam)
!        integer, intent(in)          :: iSat
!        integer                      :: rc ! error code
!        type (sat_param) :: SatParam   
! !      Other parameters
!        real(dp)   :: co
! 
! 
! !     ............................................................................
! 
! 
! !        INTEGER, PARAMETER :: in=10
! ! 
! ! 
! ! !         Character Typerun, typeinput
! ! !         Integer Code, NumSats, whichconst
! ! !         REAL*8 startmfe, stopmfe, deltamin
! ! ! 
! ! ! * ----------------------------  Locals  -------------------------------
! ! !         REAL*8 J2, mu, RadiusEarthKm,VKmPerSec, xke, tumin
! ! !         REAL*8 BC,EPDay, sec, xpdotp, j3, j4, j3oj2 
! ! !         REAL*8 startsec, stopsec, startdayofyr, stopdayofyr, jdstart, 
! ! !      &         jdstop
! ! !         INTEGER startyear, stopyear, startmon, stopmon, startday, 
! ! !      &          stopday, starthr, stophr, startmin, stopmin 
! ! !         INTEGER Yr,Mon,Day,Hr,Minute,  ICrdno,nexp,bexp, error
! ! !         CHARACTER Show
! ! 
! !          INTEGER  :: ICrdno
! ! !         REAL(dp) :: 
! !         character(len=130)  :: LongStr1, LongStr2 
! ! 
! ! 
! ! ! *****************************
! !      ! OPEN (UNIT=in, FILE='dataaqua_20091026', STATUS='OLD', ACTION="READWRITE", POSITION="REWIND", IOSTAT=ios)
! ! 
! ! 
! !         OPEN(UNIT=in, FILE = 'dataaqua_20091026', STATUS='OLD', ACCESS = 'SEQUENTIAL' )
! ! 
! !        LongStr1 = ' '
! !        !READ(10,'(a130)',END=999) LongStr1
! ! 50     READ(in,'(a130)') LongStr1
! ! ! 999    rc = 999
! ! !        return
! ! 
! ! 
! !        IF(LongStr1(1:1) .eq. '#') GOTO 50 ! Commented line of text, skip
! !        print*, LongStr1
! ! ! 
! ! !         READ(LongStr1,500) ICRDNO,SatNum,SatName,EpochYr,EpDay,
! ! !      &                       NDot,NDDot,nexp,BStar,bexp,EPHTYP,ELNO
! ! !   500   FORMAT( I1,1X,I5,1X,A10,I2,D12.0,1X,D10.0,1X,
! ! !      &          F6.5,I2,1X,F6.5,I2,1X,I1,1X,I4 )
! ! 
! ! 
! ! 
! !       CLOSE(UNIT=in,STATUS='KEEP')
! 
! 
! 
! 
! ! * ----------- READ THE SECOND LINE OF ELEMENT SET AND TIME ------------
! !         LongStr2 = ' '
! !    51   READ(10,'(a130)',END=999) LongStr2
! !         IF(LongStr2(1:1) .eq. '#') GOTO 51 ! Commented line of text, skip
! ! 
! !         IF (Typerun.eq.'V') THEN
! !           READ(LongStr2,502) ICRDNO,Inclo,nodeo,Ecco,Argpo,Mo,No,REVI,
! !      &              startmfe, stopmfe, DeltaMin
! !          else
! !           READ(LongStr2,501) ICRDNO,Inclo,nodeo,Ecco,Argpo,Mo,No,REVI
! !          endif
! !   501   FORMAT( I1,7X,D8.0,1X,D8.0,1X,F7.7,1X,D8.0,1X,D8.0,1X,D11.0,I5)
! !   502   FORMAT( I1,7X,D8.0,1X,D8.0,1X,F7.7,1X,D8.0,1X,D8.0,1X,D11.0,I5,
! !      &          1X,F12.6,F12.6,F12.6 )
! 
! ! *********************************
! 
! 
! 
! 
!        RefDate  = 20091009 ! this is date for following parameters
!        RefTime  = 021057   ! this is time 
! 
!        !CALL satname2int(Sat_name, iSat, rc)
!        if (iSat==1) then !aqua
!            SatParam%eccentricity    =  0.0000672   ! eccentricity
!            SatParam%anomaly         =  304.6193    ! mean anamoly
!            SatParam%inclination     =  98.1767     ! sattelite inclination angle
!            SatParam%ascension       =  69.3045     ! right ascension of node
!            SatParam%omega           =  55.5082     ! degree argument of perigee
!            SatParam%nrevolution     =  14.57121647  ! num of rev per data (nbar) mean motion
!            SatParam%altitude        =  702         ! km
!            SatParam%semimajor       =  -9999999.0  ! km later will be calculated
!            SatParam%period          =  -9999999.0  ! period in seconds 
!        elseif (iSat==2) then !calipso
!            SatParam%eccentricity    =  0.0001054   ! eccentricity
!            SatParam%anomaly         =  276.1217    ! mean anamoly
!            SatParam%inclination     =  98.1903     ! sattelite inclination angle
!            SatParam%ascension       =  237.3156    ! right ascension of node
!            SatParam%omega           =  84.0117     ! degree argument of perigee
!            SatParam%nrevolution     =  14.571294  ! num of rev per data (nbar) mean motion
!            SatParam%altitude        =  702         ! km
!            SatParam%semimajor       =  -9999999.0  ! km later will be calculated
!            SatParam%period          =  -9999999.0  ! period in seconds 
!        elseif (iSat==3) then !cloudsat
!            SatParam%eccentricity    =  0.0000864   ! eccentricity
!            SatParam%anomaly         =  270.1434    ! mean anamoly
!            SatParam%inclination     =  98.1903     ! sattelite inclination angle
!            SatParam%ascension       =  236.8604    ! right ascension of node
!            SatParam%omega           =  89.9856     ! degree argument of perigee
!            SatParam%nrevolution     =  14.5712970  ! num of rev per data (nbar) mean motion
!            SatParam%altitude        =  702         ! km
!            SatParam%semimajor       =  -9999999.0  ! km later will be calculated
!            SatParam%period          =  -9999999.0  ! period in seconds 
!        elseif (iSat==4) then !aura
!            SatParam%eccentricity    =  0.0001545   ! eccentricity
!            SatParam%anomaly         =  268.2933    ! mean anamoly
!            SatParam%inclination     =  98.1883     ! sattelite inclination angle
!            SatParam%ascension       =  237.1841    ! right ascension of node
!            SatParam%omega           =  91.8445     ! degree argument of perigee
!            SatParam%nrevolution     =  14.5712533  ! num of rev per data (nbar) mean motion
!            SatParam%altitude        =  702         ! km
!            SatParam%semimajor       =  -9999999.0  ! km later will be calculated
!            SatParam%period          =  -9999999.0  ! period in seconds 
!        elseif (iSat==5) then !terra
!            SatParam%eccentricity    =  0.0001307   ! eccentricity
!            SatParam%anomaly         =  293.0040    ! mean anamoly
!            SatParam%inclination     =  98.2183     ! sattelite inclination angle
!            SatParam%ascension       =  359.8416    ! right ascension of node
!            SatParam%omega           =  8.8523      ! degree argument of perigee
!            SatParam%nrevolution     =  14.57102212 ! num of rev per data (nbar) mean motion
!            SatParam%altitude        =  702         ! km
!            SatParam%semimajor       =  -9999999.0  ! km later will be calculated
!            SatParam%period          =  -9999999.0  ! period in seconds 
!        else
!               rc = 5 !satellite is not defined
!               return
!        endif 
! 
!        p(1) = SatParam%eccentricity
!        p(2) = SatParam%anomaly
!        p(3) = SatParam%inclination
!        p(4) = SatParam%ascension
!        P(5) = SatParam%omega
!        P(6) = SatParam%nrevolution
!        P(7) = SatParam%altitude
!       !P(8) = SatParam%semimajor
! 
! 
! ! Transformations
!        co   = 8681660.4;  !n~8681660.4 * a^(-3/2) Space mission analysis and design
!        P(8) = (SatParam%nrevolution/co)**(-2.0/3.0);   !7083.4456; % sattellite semi-major axis -- km
! 
!        co   = 331.24915;  ! coefficient for 
!        P(9) = (P(8)/331.24915)**(3.0/2.0)*60; ! Sat Period in seconds
!        P(9) = 16*24*60*60
!        ! print*, p
! 
!        ! print*, "       "
!        ! print*, "----------------------------------------------------------"
!        ! print*, "satellite name: ", Sat_name, " and carresponding Sat number: ", iSat
!        ! print*, "Learning is performed on:", RefDate, " at:", 000000, " hour/min/sec"
!        ! print*, "This is the ", Day2year, " th day of the year"   
!  
!        ! print*, "----------------------------------------------------------"
!        ! print*, "      "
!  
! END SUBROUTINE Orbits_Track0
! 
! SUBROUTINE Orbits_Track1(name, p)
! !SUBROUTINE Orbits_Track1(name, Satcoef_vec1, Normcoef_vec1)
!      ! See interface
!      IMPLICIT NONE
!  
!        character(len=*), intent(in) :: name ! Satellite name
!        !real(dp), intent(out), dimension(1:Num_coef) :: Satcoef_vec1
!        !real(dp), intent(out), dimension(1:Num_normcoef) :: Normcoef_vec1
!        REAL(dp), DIMENSION(:), INTENT(OUT)       :: p(num_satparam)
!        !integer                               :: start_ind, end_ind
!        !integer                               :: Normcoef_ind
!        integer                               :: iSat
!        integer                               :: rc
!  !     ...................................
! 
!         CALL satname2int(name, iSat, rc)
!         CALL Orbits_Track0(iSat, p)
!  
!  END SUBROUTINE Orbits_Track1
!  
! SUBROUTINE Orbits_Track2(iSat, p)
! !SUBROUTINE Orbits_Track2(iSat, Satcoef_vec, Normcoef_vec)
!      ! See interface
!      IMPLICIT NONE
!  
!        integer, intent(in) :: iSat ! this is the sat number ex: Aqua=1
!        REAL(dp), DIMENSION(:), INTENT(OUT)       :: p(num_satparam)
!        !real(dp), intent(out), dimension(1:Num_coef) :: Satcoef_vec
!        !real(dp), intent(out), dimension(1:Num_normcoef) :: Normcoef_vec
!        !integer                               :: start_ind, end_ind
!        !integer                               :: Normcoef_ind
!  
!  !     ...................................
!  
!        CALL Orbits_Track0(iSat, p)
! 
! END SUBROUTINE Orbits_Track2
! 
! 
! 
! 
! ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! 
! !-------------------------------------------------------------------------
! 
! 
! !                         -----------------
! !                         INTERNAL ROUTINES
! !                         -----------------
! !-------------------------------------------------------------------------
! 
! Subroutine satname2int(name, satnum, rc)
!        IMPLICIT NONE
!        character(len=*), intent(in) :: name ! Satellite name
!        integer, intent(out) :: satnum       ! Satellite number
!        integer, intent(out) :: rc
! 
!        rc = 0
!        if ((name.EQ."AQUA").or.(name.EQ."aqua").or.(name.EQ."Aqua")) then         
!            satnum = 1
!        elseif ((name.EQ." calipso").or.(name.EQ."CALIPSO").or.(name.EQ."Calipso")) then         
!            satnum = 2
!        elseif ((name.EQ." cloudsat").or.(name.EQ."CLOUDSAT").or.(name.EQ."CloudSat")) then         
!            satnum = 3
!        elseif ((name.EQ." aura").or.(name.EQ."AURA").or.(name.EQ."Aura")) then         
!            satnum = 4
!        elseif ((name.EQ." terra").or.(name.EQ."TERRA").or.(name.EQ."Terra")) then         
!            satnum = 5
!        else 
!           rc =99
!           return
!        endif
! end Subroutine satname2int
! 
! 
! ! SUBROUTINE  ECI2ECEF(iSat,fraction2day,  ECI_est, ECEF_est)
! !         IMPLICIT NONE
! !         INTEGER, INTENT(IN)  :: iSat
! !         REAL(dp), DIMENSION(:,:), INTENT(IN)  :: ECI_est
! !         REAL(dp), DIMENSION(:,:), INTENT(OUT) :: ECEF_est
! !         REAL(dp), DIMENSION(3,3) :: RotM        ! rotation matrix from In2ECEF         
! !         INTEGER   :: t1 
! !         REAL(dp)      :: Goft, a1, a2
! !         REAL(dp)      :: fraction2day 
! ! 
! !         t1 =  Day2year(iSat)
! !         Goft = g0 + g1p*t1 + g2p(iSat) * fraction2day
! !         a1 = cos(Goft);
! !         a2 = sin(Goft);
! !         RotM = 0.0
! !         RotM(1,1) = a1
! !         RotM(1,2) = a2
! !         RotM(2,1) = -a2
! !         RotM(2,2) = a1
! !         RotM(3,3) = 1
! !         ECEF_est = MATMUL(RotM, ECI_est)
! ! END SUBROUTINE ECI2ECEF  
! 
! 
! 
! SUBROUTINE ECEF2LLA(x, y, z, lat, lon, alt)
! !       % ECEF3LLA - convert earth-centered earth-fixed (ECEF)
! !       %            cartesian coordinates to latitude, longitude,
! !       %            and altitude
! !       %
! !       % USAGE:
! !       % [lat,lon,alt] = ecef2lla(x,y,z)
! !       %
! !       % lat = geodetic latitude (radians)
! !       % lon = longitude (radians)
! !       % alt = height above WGS84 ellipsoid (m)
! !       % x = ECEF X-coordinate (m)
! !       % y = ECEF Y-coordinate (m)
! !       % z = ECEF Z-coordinate (m)
! !       %
! !       % Notes: (1) This function assumes the WGS84 model.
! !       %        (2) Latitude is customary geodetic (not geocentric).
! !       %        (3) Inputs may be scalars, vectors, or matrices of the same
! !       %            size and shape. Outputs will have that same size and shape.
! !       %        (4) Tested but no warranty; use at your own risk.
! !       %        (5) Michael Kleder, April 2006
!        REAL(dp), INTENT(IN)  :: x,y,z
!        REAL(dp), INTENT(OUT) :: lat,lon,alt
! !      % WGS84 ellipsoid constants:
!        REAL(dp), PARAMETER  :: a = 6378137 ! those should be defined outside
!        REAL(dp), PARAMETER  :: e = 8.1819190842622e-2
!        REAL(dp), PARAMETER  :: pi = 3.14159265358979323846 
!        REAL(dp)             :: b, ep, p, th, N
!        LOGICAL k
! 
!        b   = sqrt(a**2 * (1-e**2))
!        ep  = sqrt((a**2-b**2)/b**2)
!        p   = sqrt(x**2+y**2)
!        th  = atan2(a*z,b*p)
!        lon = atan2(y,x)
!        lat = atan2((z+ep**2*b*sin(th)**3),(p-e**2*a*cos(th)**3))
!        N   = a/sqrt(1-e**2*sin(lat)**2)
!        alt = p/cos(lat)-N
! !%     !return lon in range [0,2*pi)
!        !lon = mod(lon,2*pi) ! does not the same as matlab
!        lon = lon-floor(lon/(2*pi))*(2*pi)
! !      % correct for numerical instability in altitude near exact poles:
! !      % (after this correction, error is about 2 millimeters, which is about
! !      % the same as the numerical precision of the overall function)
!        !k= (abs(x)<1.AND.abs(y)<1 )
!        !print *, k
!        !!alt(k) = abs(z(k))-b
! END SUBROUTINE ECEF2LLA
! 
! 
! 
! 
! 
! integer function ODS_Julian ( CalDate )
!       implicit NONE
!       integer  CalDate  ! Calendar date in the format YYYYMMDD
!                         !   where YYYY is the year, MM is the
!                         !   month and DD is the day.  A negative
!                         !   number implies that the year is B.C.
! !     Other variables
! !     ---------------
!       integer     Year
!       integer     Month
!       integer     Day
!       integer     iGreg  ! Gregorian Calendar adopted Oct 12, 1582
!       parameter ( iGreg = 15 + 31 * ( 10 + 12 * 1582 ) )
!       integer     JulDay
!       integer     jy, jm, ja
! 
!       Year   =       CalDate / 10000
!       Month  = mod ( CalDate,  10000 ) / 100
!       Day    = mod ( CalDate,    100 )
! !     Change year 0 to year 1
!       if ( Year  .eq. 0 ) Year = 1
! !     Account for the nonexisting year 0
!       if ( Year  .lt. 0 ) Year = Year + 1
!       if ( Month .gt. 2 ) then
!          jy = Year
!          jm = Month + 1
!       else
!          jy = Year  - 1
!          jm = Month + 13
!       endif
!       JulDay = int ( 365.25  * jy )       &
!              + int ( 30.6001 * jm )       &
!              + Day + 1720995
! !     Test whether to change to Gregorian Celendar
!       if ( Day + 31 * ( Month + 12 * Year ) .ge. iGreg) then
!         ja     = int ( 0.01 * jy )
!         Julday = JulDay + 2 - ja + int ( 0.25 * ja )
!       endif
!       ODS_Julian = JulDay
!       return
! end function ODS_Julian
! 
! REAL(dp) FUNCTION get_fraction(t1)
!        IMPLICIT NONE
!        REAL(dp), INTENT(IN) :: t1
!        REAL(dp) :: fraction2day, temp3, temp4 
! 
!        if (t1.EQ.0.0) then
!            fraction2day = 0.0
!        else if (t1.lt.1.0 .AND. t1.gt.0) then
!            fraction2day = t1
!        else
!            temp3 = floor(t1)
!            temp4 = t1
!            if (temp3 == temp4) then 
!              fraction2day = temp3 ! either is ok 
!            else
!              fraction2day = MOD(temp4,temp3)
!            end if
!        end if
!        get_fraction = fraction2day
! END FUNCTION get_fraction
! 
! 
! REAL(dp) FUNCTION timeday2sidereal(tm)
!        IMPLICIT NONE
! !      Ref:  Fundementals of astrodynamics and Applications Vallado p192  
! !      output in degrees
!        INTEGER, DIMENSION(1:6), INTENT(IN) :: tm
!        INTEGER :: yr, mo, day, hr, mi, sec 
!        INTEGER, PARAMETER :: hr2sec = 60*60
!        REAL(dp) :: temp1, JD, TUT1, TETAGMST
! !timem = [1992 08  20 12 14 00]
!        yr   = tm(1);
!        mo   = tm(2);
!        day  = tm(3);
!        hr   = tm(4);
!        mi   = tm(5);
!        sec  = tm(6);
! !      calc jul date 
!        temp1 = ((real(sec,dp)/60.0 + real(mi,dp))/60.0) + real(hr,dp);
!        JD = 367 * real(yr,dp) - floor(  7*  (real(yr,dp) + floor( (real(mo,dp)+9.0)/12.0))  / 4.0)    +   &
!             floor(275*real(mo,dp)/9.0)  + real(day,dp) + 1721013.5 + temp1/24.0   !2.455113590937500e+06
! !      Julian centuries for particular epoch (J200)
!        TUT1 = (JD - 2451545.0) / 36525.0;
! !      find GMST
!        TETAGMST = 67310.54841 + (876600.0 * real(hr2sec,dp) + 8640184.812866)*TUT1 +  &
!                   0.093104 * (TUT1**2) - 6.2*(10**(-6)) * (TUT1**3);  ! need extended precision  
!       TETAGMST=mod(TETAGMST,86400.0);
!       TETAGMST=TETAGMST/240; !seconds to degrees ! tothours = TETAGMST/15
!       timeday2sidereal = TETAGMST
! END FUNCTION timeday2sidereal
! 
! 
! 
! 
! 
! SUBROUTINE jd2sse      ( jd, Direction, sse )
!        IMPLICIT NONE
!        REAL(dp), INTENT(INOUT) ::  sse, jd
!        CHARACTER*4 Direction
! 
!         ! --------------------  Implementation   ----------------------
!         IF ( Direction.eq.'FROM' ) THEN
!             jd = 2451544.5D0 + sse/86400.0D0
!           ELSE
!             sse = (jd - 2451544.5D0) * 86400.0D0
!           ENDIF
!       RETURN
! END SUBROUTINE jd2sse
! 
! 
! subroutine built_vecsize(sim_start, sim_end, time_day, say)
!        implicit none
!        real(dp), intent(in) :: sim_start, sim_end, time_day
!        integer          :: say
!        real(dp)             :: temp1
! 
!        say = 1
!        temp1 = sim_start
! 
!        !print*, "inside simend-----", sim_end, temp1, time_day 
!        !WRITE(*,"(1F10.2)")  sim_end
! 
!        do while (sim_end.GE.temp1)
!          temp1 = temp1 + time_day
!          if (temp1.LE.sim_end) then
!           say = say + 1 
!            end if
!             enddo
! end subroutine built_vecsize
! 
! subroutine built_vec(sim_start, sim_end, time_day, t1)
!        implicit none
!        real(dp), intent(in) :: sim_start, sim_end, time_day
!        integer          :: say
!        real(dp)             :: temp1
!        real(dp),  dimension(:) :: t1 
! 
!        say = 1         
!        temp1 = sim_start ! garante that it lesim start !time
!        t1(1) = sim_start
!        do while (sim_end.GE.temp1)
!          temp1 = temp1 + time_day
!           say = say + 1 
!           if (say>size(t1)) then
!              exit 
!           end if
!           t1(say) = temp1
!        enddo
! end subroutine built_vec
! 
! 
! 
! 
!       
! END MODULE orbits_mod1 ! ........END........ 
! 
! 
! 
! 
! ! ****************************************************
! ! ****************************************************
! SPA.h
! 
! 
! ! --------------------------------------------------
! MODULE SPAH_1         ! global declarations
! ! --------------------------------------------------
! INCLUDE 'C2F.FD' 
! INCLUDE 'C2F.FI' 
!   INTEGER,PARAMETER :: SPA_ZA=0
!   INTEGER,PARAMETER :: SPA_ZA_INC=0+1
!   INTEGER,PARAMETER :: SPA_ZA_RTS=0+1+1
!   INTEGER,PARAMETER :: SPA_ALL=0+1+1+1
!   TYPE :: SPA_DATA_T 
!     INTEGER :: year 
!     INTEGER :: month 
!     INTEGER :: day 
!     INTEGER :: hour 
!     INTEGER :: minute 
!     INTEGER :: second 
!     REAL(8) :: delta_t 
!     REAL(8) :: timezone 
!     REAL(8) :: longitude 
!     REAL(8) :: latitude 
!     REAL(8) :: elevation 
!     REAL(8) :: pressure 
!     REAL(8) :: temperature 
!     REAL(8) :: slope 
!     REAL(8) :: azm_rotation 
!     REAL(8) :: atmos_refract 
!     INTEGER :: function 
!     REAL(8) :: jd 
!     REAL(8) :: jc 
!     REAL(8) :: jde 
!     REAL(8) :: jce 
!     REAL(8) :: jme 
!     REAL(8) :: l 
!     REAL(8) :: b 
!     REAL(8) :: r 
!     REAL(8) :: theta 
!     REAL(8) :: beta 
!     REAL(8) :: x0 
!     REAL(8) :: x1 
!     REAL(8) :: x2 
!     REAL(8) :: x3 
!     REAL(8) :: x4 
!     REAL(8) :: del_psi 
!     REAL(8) :: del_epsilon 
!     REAL(8) :: epsilon0 
!     REAL(8) :: epsilon 
!     REAL(8) :: del_tau 
!     REAL(8) :: lamda 
!     REAL(8) :: nu0 
!     REAL(8) :: nu 
!     REAL(8) :: alpha 
!     REAL(8) :: delta 
!     REAL(8) :: h 
!     REAL(8) :: xi 
!     REAL(8) :: del_alpha 
!     REAL(8) :: delta_prime 
!     REAL(8) :: alpha_prime 
!     REAL(8) :: h_prime 
!     REAL(8) :: e0 
!     REAL(8) :: del_e 
!     REAL(8) :: e 
!     REAL(8) :: eot 
!     REAL(8) :: srha 
!     REAL(8) :: ssha 
!     REAL(8) :: sta 
!     REAL(8) :: zenith 
!     REAL(8) :: azimuth180 
!     REAL(8) :: azimuth 
!     REAL(8) :: incidence 
!     REAL(8) :: suntransit 
!     REAL(8) :: sunrise 
!     REAL(8) :: sunset 
!   END TYPE
!   TYPE (SPA_DATA_T) :: spa_data
! 
! END MODULE
! INCLUDE 'C2F_LIB.F90'  
! !  0 Errors detected
! 
! 
! 
! 
! 
! 
! ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! SPA.c
! 
! 
! ! --------------------------------------------------
! MODULE SPA_1         ! global declarations
! ! --------------------------------------------------
! INCLUDE 'C2F.FD' 
! INCLUDE 'C2F.FI' 
!   INTEGER,PARAMETER :: SPA_ZA=0
!   INTEGER,PARAMETER :: SPA_ZA_INC=0+1
!   INTEGER,PARAMETER :: SPA_ZA_RTS=0+1+1
!   INTEGER,PARAMETER :: SPA_ALL=0+1+1+1
!   INTEGER,PARAMETER :: TERM_A=0
!   INTEGER,PARAMETER :: TERM_B=0+1
!   INTEGER,PARAMETER :: TERM_C=0+1+1
!   INTEGER,PARAMETER :: TERM_COUNT=0+1+1+1
!   INTEGER,PARAMETER :: TERM_X0=0
!   INTEGER,PARAMETER :: TERM_X1=0+1
!   INTEGER,PARAMETER :: TERM_X2=0+1+1
!   INTEGER,PARAMETER :: TERM_X3=0+1+1+1
!   INTEGER,PARAMETER :: TERM_X4=0+1+1+1+1
!   INTEGER,PARAMETER :: TERM_X_COUNT=0+1+1+1+1+1
!   INTEGER,PARAMETER :: TERM_PSI_A=0
!   INTEGER,PARAMETER :: TERM_PSI_B=0+1
!   INTEGER,PARAMETER :: TERM_EPS_C=0+1+1
!   INTEGER,PARAMETER :: TERM_EPS_D=0+1+1+1
!   INTEGER,PARAMETER :: TERM_PE_COUNT=0+1+1+1+1
!   INTEGER,PARAMETER :: JD_MINUS=0
!   INTEGER,PARAMETER :: JD_ZERO=0+1
!   INTEGER,PARAMETER :: JD_PLUS=0+1+1
!   INTEGER,PARAMETER :: JD_COUNT=0+1+1+1
!   INTEGER,PARAMETER :: SUN_TRANSIT=0
!   INTEGER,PARAMETER :: SUN_RISE=0+1
!   INTEGER,PARAMETER :: SUN_SET=0+1+1
!   INTEGER,PARAMETER :: SUN_COUNT=0+1+1+1
!   TYPE :: SPA_DATA_T 
!     INTEGER :: year 
!     INTEGER :: month 
!     INTEGER :: day 
!     INTEGER :: hour 
!     INTEGER :: minute 
!     INTEGER :: second 
!     REAL(8) :: delta_t 
!     REAL(8) :: timezone 
!     REAL(8) :: longitude 
!     REAL(8) :: latitude 
!     REAL(8) :: elevation 
!     REAL(8) :: pressure 
!     REAL(8) :: temperature 
!     REAL(8) :: slope 
!     REAL(8) :: azm_rotation 
!     REAL(8) :: atmos_refract 
!     INTEGER :: function 
!     REAL(8) :: jd 
!     REAL(8) :: jc 
!     REAL(8) :: jde 
!     REAL(8) :: jce 
!     REAL(8) :: jme 
!     REAL(8) :: l 
!     REAL(8) :: b 
!     REAL(8) :: r 
!     REAL(8) :: theta 
!     REAL(8) :: beta 
!     REAL(8) :: x0 
!     REAL(8) :: x1 
!     REAL(8) :: x2 
!     REAL(8) :: x3 
!     REAL(8) :: x4 
!     REAL(8) :: del_psi 
!     REAL(8) :: del_epsilon 
!     REAL(8) :: epsilon0 
!     REAL(8) :: epsilon 
!     REAL(8) :: del_tau 
!     REAL(8) :: lamda 
!     REAL(8) :: nu0 
!     REAL(8) :: nu 
!     REAL(8) :: alpha 
!     REAL(8) :: delta 
!     REAL(8) :: h 
!     REAL(8) :: xi 
!     REAL(8) :: del_alpha 
!     REAL(8) :: delta_prime 
!     REAL(8) :: alpha_prime 
!     REAL(8) :: h_prime 
!     REAL(8) :: e0 
!     REAL(8) :: del_e 
!     REAL(8) :: e 
!     REAL(8) :: eot 
!     REAL(8) :: srha 
!     REAL(8) :: ssha 
!     REAL(8) :: sta 
!     REAL(8) :: zenith 
!     REAL(8) :: azimuth180 
!     REAL(8) :: azimuth 
!     REAL(8) :: incidence 
!     REAL(8) :: suntransit 
!     REAL(8) :: sunrise 
!     REAL(8) :: sunset 
!   END TYPE
!   TYPE (SPA_DATA_T) :: spa_data
!   INTEGER,PARAMETER :: l_subcount(6 )= (/64,34,20,7,3,1/)
!   INTEGER,PARAMETER :: b_subcount(2 )= (/5,2/)
!   INTEGER,PARAMETER :: r_subcount(5 )= (/40,10,6,2,1/)
!   REAL(8),PARAMETER :: L_TERMS(TERM_C)= (/175347046.0,0,0/)
!   REAL(8),PARAMETER :: L_TERMS(TERM_C)= (/206059.0,2.678235,6283.07585/)
!   REAL(8),PARAMETER :: L_TERMS(TERM_C)= (/8720.0,1.0721,6283.0758/)
!   REAL(8),PARAMETER :: L_TERMS(TERM_C)= (/35,0,0/)
!   REAL(8),PARAMETER :: L_TERMS(TERM_C)= (/8,4.13,6283.08/)
!   REAL(8),PARAMETER :: L_TERMS(TERM_C)= (/1,3.14,0,  &
!   REAL(8),PARAMETER :: B_TERMS(TERM_C)= (/280.0,3.199,84334.662/)
!   REAL(8),PARAMETER :: B_TERMS(TERM_C)= (/6,1.73,5223.69/)
!   REAL(8),PARAMETER :: R_TERMS(TERM_C)= (/100013989.0,0,0/)
!   REAL(8),PARAMETER :: R_TERMS(TERM_C)= (/1721.0,1.0644,12566.1517/)
!   REAL(8),PARAMETER :: R_TERMS(TERM_C)= (/124.0,5.579,12566.152/)
!   REAL(8),PARAMETER :: R_TERMS(TERM_C)= (/7,3.92,12566.15/)
!   REAL(8),PARAMETER :: R_TERMS(TERM_C)= (/4,2.56,6283.08,  &
!   REAL(8),PARAMETER :: L_TERMS(TERM_C)= (/8,4.13,6283.08/)
!  
! ! --------------------------------------------------
! FUNCTION rad2deg(radians)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: radians
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = (180.0/ 3.1415926535897932384626433832795028841971 )*radians
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION deg2rad(degrees)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: degrees
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = (3.1415926535897932384626433832795028841971 /180.0)*degrees
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION limit_degrees(degrees_)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: degrees_,degrees
! ! - - - local declarations - - -
!   REAL(8) :: limited
! 
! ! - - - begin - - -
!   degrees = degrees_
!   degrees = degrees / 360.0
!   limited = 360.0*(degrees-floor(degrees))
!   IF (limited < 0) limited = limited + 360.0
!   output_8 = limited
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION limit_degrees180pm(degrees_)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: degrees_,degrees
! ! - - - local declarations - - -
!   REAL(8) :: limited
! 
! ! - - - begin - - -
!   degrees = degrees_
!   degrees = degrees / 360.0
!   limited = 360.0*(degrees-floor(degrees))
!   IF (limited < -180.0) THEN
!     limited = limited + 360.0
!   ELSE IF (limited > 180.0) THEN
!     limited = limited - 360.0
!   END IF
!   output_8 = limited
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION limit_degrees180(degrees_)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: degrees_,degrees
! ! - - - local declarations - - -
!   REAL(8) :: limited
! 
! ! - - - begin - - -
!   degrees = degrees_
!   degrees = degrees / 180.0
!   limited = 180.0*(degrees-floor(degrees))
!   IF (limited < 0) limited = limited + 180.0
!   output_8 = limited
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION limit_zero2one(value)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: value
! ! - - - local declarations - - -
!   REAL(8) :: limited
! 
! ! - - - begin - - -
!   limited = value - floor(value)
!   IF (limited < 0) limited = limited + 1.0
!   output_8 = limited
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION limit_minutes(minutes)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: minutes
! ! - - - local declarations - - -
!   REAL(8) :: limited
! 
! ! - - - begin - - -
!   limited = minutes
!   IF (limited < -20.0) THEN
!     limited = limited + 1440.0
!   ELSE IF (limited > 20.0) THEN
!     limited = limited - 1440.0
!   END IF
!   output_8 = limited
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION dayfrac_to_local_hr(dayfrac,timezone)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: dayfrac,timezone
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = 24.0*limit_zero2one(dayfrac + timezone/24.0)
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION third_order_polynomial(a,b,c,d,x)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: a,b,c,d,x
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = ((a*x + b)*x + c)*x + d
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION validate_inputs(spa)  RESULT (output_4)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   INTEGER :: output_4                                       
!   TYPE (SPA_DATA_T) :: spa
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   IF ((spa%year < -2000) .OR. (spa%year > 6000)) THEN
!     output_4 = 1
!     RETURN
!   END IF
!     IF ((spa%month < 1 ) .OR. (spa%month > 12 )) THEN
!       output_4 = 2
!       RETURN
!     END IF
!     IF ((spa%day < 1 ) .OR. (spa%day > 31 )) THEN
!       output_4 = 3
!       RETURN
!     END IF
!     IF ((spa%hour < 0 ) .OR. (spa%hour > 24 )) THEN
!       output_4 = 4
!       RETURN
!     END IF
!     IF ((spa%minute < 0 ) .OR. (spa%minute > 59 )) THEN
!       output_4 = 5
!       RETURN
!     END IF
!     IF ((spa%second < 0 ) .OR. (spa%second > 59 )) THEN
!       output_4 = 6
!       RETURN
!     END IF
!     IF ((spa%pressure < 0 ) .OR. (spa%pressure > 5000)) THEN
!       output_4 = 12
!       RETURN
!     END IF
!     IF ((spa%temperature <= -273) .OR. (spa%temperature > 6000)) THEN
!       output_4 = 13
!       RETURN
!     END IF
!     IF ((spa%hour == 24 ) .AND. (spa%minute > 0 )) THEN
!       output_4 = 5
!       RETURN
!     END IF
!     IF ((spa%hour == 24 ) .AND. (spa%second > 0 )) THEN
!       output_4 = 6
!       RETURN
!     END IF
!     IF (ABS(spa%delta_t) > 8000 ) THEN
!       output_4 = 7
!       RETURN
!     END IF
!     IF (ABS(spa%timezone) > 18 ) THEN
!       output_4 = 8
!       RETURN
!     END IF
!     IF (ABS(spa%longitude) > 180 ) THEN
!       output_4 = 9
!       RETURN
!     END IF
!     IF (ABS(spa%latitude) > 90 ) THEN
!       output_4 = 10
!       RETURN
!     END IF
!     IF (ABS(spa%atmos_refract) > 5 ) THEN
!       output_4 = 16
!       RETURN
!     END IF
!     IF ( spa%elevation < -6500000) THEN
!       output_4 = 11
!       RETURN
!     END IF
!   IF ((spa%function == SPA_ZA_INC) .OR. (spa%function == SPA_ALL)) THEN
!     IF (ABS(spa%slope) > 360) THEN
!       output_4 = 14
!       RETURN
!     END IF
!       IF (ABS(spa%azm_rotation) > 360) THEN
!         output_4 = 15
!         RETURN
!       END IF
!   END IF
!   output_4 = 0
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION julian_day(year_,month_,day,hour,minute,second,tz)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   INTEGER :: year_,year,month_,month,day,hour,minute,second
!   REAL(8) :: tz
! ! - - - local declarations - - -
!   REAL(8) :: day_decimal,julian_day,a
! 
! ! - - - begin - - -
!   year = year_
!   month = month_
!   day_decimal = day + (hour - tz + (minute + second/60.0)/60.0)/24.0
!   IF (month < 3) THEN
!     month = month + 12
!     year = year-1
!   END IF
!   julian_day = floor(365.25*(year+4716.0)) + floor(30.6001*(month+1)) + day_decimal - 1524.5
!   IF (julian_day > 2299160.0) THEN
!     a = floor(year/100)
!     julian_day = julian_day +( (2 - a + floor(a/4)))
!   END IF
!   output_8 = julian_day
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION julian_century(jd)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: jd
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = (jd-2451545.0)/36525.0
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION julian_ephemeris_day(jd,delta_t)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: jd,delta_t
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = jd+delta_t/86400.0
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION julian_ephemeris_century(jde)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: jde
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = (jde - 2451545.0)/36525.0
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION julian_ephemeris_millennium(jce)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: jce
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = (jce/10.0)
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION earth_periodic_term_summation(terms,count,jme)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   INTEGER :: count
!   REAL(8) :: jme,terms(TERM_COUNT,)
! ! - - - local declarations - - -
!   INTEGER :: i
!   REAL(8) :: sum
! 
! ! - - - begin - - -
!   sum = 0
!   DO i = 0,count-1
!     sum = sum + terms(TERM_A+1,i+1)*COS(terms(TERM_B+1,i+1)+terms(TERM_C+1,i+1)*jme)
!   END DO
!   output_8 = sum
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION earth_values(term_sum,count,jme)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   INTEGER :: count
!   REAL(8) :: jme,term_sum(*)
! ! - - - local declarations - - -
!   INTEGER :: i
!   REAL(8) :: sum
! 
! ! - - - begin - - -
!   sum = 0
!   DO i = 0,count-1
!     sum = sum + term_sum(i+1)*jme**i
!   END DO
!   sum = sum / 1.0e8
!   output_8 = sum
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION earth_heliocentric_longitude(jme)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: jme
! ! - - - local declarations - - -
!   INTEGER :: i
!   REAL(8) :: sum(6 )
! 
! ! - - - begin - - -
!   DO i = 0,6-1
!     sum(i+1) = earth_periodic_term_summation(L_TERMS(i+1), l_subcount(i+1), jme)
!   END DO
!   output_8 = limit_degrees(rad2deg(earth_values(sum, 6 , jme)))
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION earth_heliocentric_latitude(jme)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: jme
! ! - - - local declarations - - -
!   INTEGER :: i
!   REAL(8) :: sum(2 )
! 
! ! - - - begin - - -
!   DO i = 0,2-1
!     sum(i+1) = earth_periodic_term_summation(B_TERMS(i+1), b_subcount(i+1), jme)
!   END DO
!   output_8 = rad2deg(earth_values(sum, 2 , jme))
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION earth_radius_vector(jme)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: jme
! ! - - - local declarations - - -
!   INTEGER :: i
!   REAL(8) :: sum(5 )
! 
! ! - - - begin - - -
!   DO i = 0,5-1
!     sum(i+1) = earth_periodic_term_summation(R_TERMS(i+1), r_subcount(i+1), jme)
!   END DO
!   output_8 = earth_values(sum, 5 , jme)
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION geocentric_longitude(l)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: l
! ! - - - local declarations - - -
!   REAL(8) :: theta
! 
! ! - - - begin - - -
!   theta = l + 180.0
!   IF (theta >= 360.0) theta = theta - 360.0
!   output_8 = theta
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION geocentric_latitude(b)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: b
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = -b
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION mean_elongation_moon_sun(jce)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: jce
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = third_order_polynomial(1.0/189474.0, -0.0019142, 445267.11148, 297.85036, jce)
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION mean_anomaly_sun(jce)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: jce
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = third_order_polynomial(-1.0/300000.0, -0.0001603, 35999.05034, 357.52772, jce)
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION mean_anomaly_moon(jce)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: jce
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = third_order_polynomial(1.0/56250.0, 0.0086972, 477198.867398, 134.96298, jce)
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION argument_latitude_moon(jce)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: jce
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = third_order_polynomial(1.0/327270.0, -0.0036825, 483202.017538, 93.27191, jce)
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION ascending_longitude_moon(jce)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: jce
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = third_order_polynomial(1.0/450000.0, 0.0020708, -1934.136261, 125.04452, jce)
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION xy_term_summation(i,x)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   INTEGER :: i
!   REAL(8) :: x(*)
! ! - - - local declarations - - -
!   INTEGER :: j
!   REAL(8) :: sum
! 
! ! - - - begin - - -
!   sum = 0
!   DO j = 0,TERM_X_COUNT-1
!     sum = sum + x(j+1)*Y_TERMS(j+1,i+1)
!   END DO
!   output_8 = sum
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! SUBROUTINE nutation_longitude_and_obliquity(jce,x,del_psi,del_epsilon)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: jce,x(*),del_psi,del_epsilon
! ! - - - local declarations - - -
!   INTEGER :: i
!   REAL(8) :: xy_term_sum,sum_psi,sum_epsilon
! 
! ! - - - begin - - -
!   sum_psi = 0
!   sum_epsilon = 0
!   DO i = 0,63-1
!     xy_term_sum = deg2rad(xy_term_summation(i, x))
!     sum_psi = sum_psi +( (PE_TERMS(TERM_PSI_A+1,i+1) + jce*PE_TERMS(TERM_PSI_B+1,i+1))*SIN(xy_term_sum))
!     sum_epsilon = sum_epsilon +( (PE_TERMS(TERM_EPS_C+1,i+1) + jce*PE_TERMS(TERM_EPS_D+1,i+1))*COS(xy_term_sum))
!   END DO
!   del_psi = sum_psi / 36000000.0
!   del_epsilon = sum_epsilon / 36000000.0
! 
! END SUBROUTINE
!  
! ! --------------------------------------------------
! FUNCTION ecliptic_mean_obliquity(jme)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: jme
! ! - - - local declarations - - -
!   REAL(8) :: u
! 
! ! - - - begin - - -
!   u = jme/10.0
!   output_8 =
!   RETURN
! (*( -39.05 + u*( 7.12 + u*( 27.87 + u*( 5.79 + u*2.45)))))))))
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION ecliptic_true_obliquity(delta_epsilon,epsilon0)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: delta_epsilon,epsilon0
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = delta_epsilon + epsilon0/3600.0
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION aberration_correction(r)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: r
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = -20.4898 / (3600.0*r)
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION apparent_sun_longitude(theta,delta_psi,delta_tau)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: theta,delta_psi,delta_tau
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = theta + delta_psi + delta_tau
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION greenwich_mean_sidereal_time(jd,jc)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: jd,jc
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 =
!   RETURN
! jc*(0.000387933 - jc/38710000.0))
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION greenwich_sidereal_time(nu0,delta_psi,epsilon)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: nu0,delta_psi,epsilon
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = nu0 + delta_psi*COS(deg2rad(epsilon))
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION geocentric_sun_right_ascension(lamda,epsilon,beta)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: lamda,epsilon,beta
! ! - - - local declarations - - -
!   REAL(8) :: lamda_rad,epsilon_rad
! 
! ! - - - begin - - -
!   lamda_rad = deg2rad(lamda)
!   epsilon_rad = deg2rad(epsilon)
!   output_8 =
!   RETURN
! CALL tan(deg2rad(beta))
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION geocentric_sun_declination(beta,epsilon,lamda)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: beta,epsilon,lamda
! ! - - - local declarations - - -
!   REAL(8) :: beta_rad,epsilon_rad
! 
! ! - - - begin - - -
!   beta_rad = deg2rad(beta)
!   epsilon_rad = deg2rad(epsilon)
!   output_8 =
!   RETURN
! CALL COS(beta_rad)
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION observer_hour_angle(nu,longitude,alpha_deg)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: nu,longitude,alpha_deg
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = limit_degrees(nu + longitude - alpha_deg)
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION sun_equatorial_horizontal_parallax(r)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: r
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = 8.794 / (3600.0 * r)
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! SUBROUTINE sun_right_ascension_parallax_and_topocentric_dec(latitude,elevation,xi,h,delta,delta_alpha,delta_prime)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: latitude,elevation,xi,h,delta,delta_alpha,delta_prime
! ! - - - local declarations - - -
!   REAL(8) :: delta_alpha_rad,lat_rad,xi_rad,h_rad,delta_rad,u,y,x
! 
! ! - - - begin - - -
!   lat_rad = deg2rad(latitude)
!   xi_rad = deg2rad(xi)
!   h_rad = deg2rad(h)
!   delta_rad = deg2rad(delta)
!   u = atan(0.99664719 * tan(lat_rad))
!   y = 0.99664719 * SIN(u) + elevation*SIN(lat_rad)/6378140.0
!   x = COS(u) + elevation*COS(lat_rad)/6378140.0
!   delta_alpha_rad = atan2( - x*SIN(xi_rad) *SIN(h_rad),  &
!       COS(delta_rad) - x*SIN(xi_rad) *COS(h_rad))
!   delta_prime = rad2deg(atan2((SIN(delta_rad) - y*SIN(xi_rad))*COS(delta_alpha_rad),  &
!       COS(delta_rad) - x*SIN(xi_rad) *COS(h_rad)))
!   delta_alpha = rad2deg(delta_alpha_rad)
! 
! END SUBROUTINE
!  
! ! --------------------------------------------------
! FUNCTION topocentric_sun_right_ascension(alpha_deg,delta_alpha)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: alpha_deg,delta_alpha
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = alpha_deg + delta_alpha
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION topocentric_local_hour_angle(h,delta_alpha)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: h,delta_alpha
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = h - delta_alpha
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION topocentric_elevation_angle(latitude,delta_prime,h_prime)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: latitude,delta_prime,h_prime
! ! - - - local declarations - - -
!   REAL(8) :: lat_rad,delta_prime_rad
! 
! ! - - - begin - - -
!   lat_rad = deg2rad(latitude)
!   delta_prime_rad = deg2rad(delta_prime)
!   output_8 =
!   RETURN
! CALL COS(lat_rad)
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION atmospheric_refraction_correction(pressure,temperature,atmos_refract,e0)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: pressure,temperature,atmos_refract,e0
! ! - - - local declarations - - -
!   REAL(8) :: del_e
! 
! ! - - - begin - - -
!   del_e = 0
!   IF (e0 >= -1*(0.26667  + atmos_refract))  &
!   del_e = (pressure / 1010.0) * (283.0 / (273.0 + temperature)) *  &
!       1.02 / (60.0 * tan(deg2rad(e0 + 10.3/(e0 + 5.11))))
!   output_8 = del_e
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION topocentric_elevation_angle_corrected(e0,delta_e)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: e0,delta_e
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = e0 + delta_e
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION topocentric_zenith_angle(e)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: e
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = 90.0 - e
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION topocentric_azimuth_angle_neg180_180(h_prime,latitude,delta_prime)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: h_prime,latitude,delta_prime
! ! - - - local declarations - - -
!   REAL(8) :: h_prime_rad,lat_rad
! 
! ! - - - begin - - -
!   h_prime_rad = deg2rad(h_prime)
!   lat_rad = deg2rad(latitude)
!   output_8 =
!   RETURN
! CALL COS(h_prime_rad)
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION topocentric_azimuth_angle_zero_360(azimuth180)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: azimuth180
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = azimuth180 + 180.0
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION surface_incidence_angle(zenith,azimuth180,azm_rotation,slope)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: zenith,azimuth180,azm_rotation,slope
! ! - - - local declarations - - -
!   REAL(8) :: zenith_rad,slope_rad
! 
! ! - - - begin - - -
!   zenith_rad = deg2rad(zenith)
!   slope_rad = deg2rad(slope)
!   output_8 =
!   RETURN
! CALL SIN(slope_rad)
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION sun_mean_longitude(jme)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: jme
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 =
!   RETURN
! (*(1/49931.0 + jme*(-1/15300.0 + jme*(-1/2000000.0))))))
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION eot(m,alpha,del_psi,epsilon)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: m,alpha,del_psi,epsilon
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = limit_minutes(4.0*(m - 0.0057183 - alpha + del_psi*COS(deg2rad(epsilon))))
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION approx_sun_transit_time(alpha_zero,longitude,nu)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: alpha_zero,longitude,nu
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 = (alpha_zero - longitude - nu) / 360.0
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION sun_hour_angle_at_rise_set(latitude,delta_zero,h0_prime)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: latitude,delta_zero,h0_prime
! ! - - - local declarations - - -
!   REAL(8) :: h0,latitude_rad,delta_zero_rad,argument
! 
! ! - - - begin - - -
!   h0 = -99999
!   latitude_rad = deg2rad(latitude)
!   delta_zero_rad = deg2rad(delta_zero)
!   argument = (SIN(deg2rad(h0_prime)) - SIN(latitude_rad)*SIN(delta_zero_rad)) /
!   (COS(latitude_rad)*COS(delta_zero_rad))
!   IF (ABS(argument) <= 1) h0 = limit_degrees180(rad2deg(acos(argument)))
!   output_8 = h0
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! SUBROUTINE approx_sun_rise_and_set(m_rts,h0)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: h0,m_rts(*)
! ! - - - local declarations - - -
!   REAL(8) :: h0_dfrac
! 
! ! - - - begin - - -
!   h0_dfrac = h0/360.0
!   m_rts(SUN_RISE+1) = limit_zero2one(m_rts(SUN_TRANSIT+1) - h0_dfrac)
!   m_rts(SUN_SET+1) = limit_zero2one(m_rts(SUN_TRANSIT+1) + h0_dfrac)
!   m_rts(SUN_TRANSIT+1) = limit_zero2one(m_rts(SUN_TRANSIT+1))
! 
! END SUBROUTINE
!  
! ! --------------------------------------------------
! FUNCTION rts_alpha_delta_prime(ad,n)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: n,ad(*)
! ! - - - local declarations - - -
!   REAL(8) :: a,b
! 
! ! - - - begin - - -
!   a = ad(JD_ZERO] - ad(JD_MINUS+1)
!   b = ad(JD_PLUS] - ad(JD_ZERO+1)
!   IF (ABS(a) >= 2.0) a = limit_zero2one(a)
!     IF (ABS(b) >= 2.0) b = limit_zero2one(b)
!   output_8 = ad(JD_ZERO+1) + n * (a + b + (b-a)*n)/2.0
!   RETURN
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION rts_sun_altitude(latitude,delta_prime,h_prime)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   REAL(8) :: latitude,delta_prime,h_prime
! ! - - - local declarations - - -
!   REAL(8) :: latitude_rad,delta_prime_rad
! 
! ! - - - begin - - -
!   latitude_rad = deg2rad(latitude)
!   delta_prime_rad = deg2rad(delta_prime)
!   output_8 =
!   RETURN
! CALL COS(latitude_rad)
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! FUNCTION sun_rise_and_set(m_rts,h_rts,delta_prime,latitude,h_prime,h0_prime,sun)  RESULT (output_8)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   REAL(8) :: output_8                                       
!   INTEGER :: sun
!   REAL(8) :: latitude,h0_prime,m_rts,h_rts,delta_prime(*),h_prime(*)
! ! - - - local declarations - - -
! 
! ! - - - begin - - -
!   output_8 =
!   RETURN
! (360.0*COS(deg2rad(delta_prime(sun+1)))*COS(deg2rad(latitude))*SIN(deg2rad(h_prime(sun+1))))
! 
! END FUNCTION
!  
! ! --------------------------------------------------
! SUBROUTINE calculate_geocentric_sun_right_ascension_and_declination(spa)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   TYPE (SPA_DATA_T) :: spa
! ! - - - local declarations - - -
!   REAL(8) :: x(TERM_X_COUNT)
! 
! ! - - - begin - - -
!   spa%jc = julian_century(spa%jd)
!   spa%jde = julian_ephemeris_day(spa%jd, spa%delta_t)
!   spa%jce = julian_ephemeris_century(spa%jde)
!   spa%jme = julian_ephemeris_millennium(spa%jce)
!   spa%l = earth_heliocentric_longitude(spa%jme)
!   spa%b = earth_heliocentric_latitude(spa%jme)
!   spa%r = earth_radius_vector(spa%jme)
!   spa%theta = geocentric_longitude(spa%l)
!   spa%beta = geocentric_latitude(spa%b)
!   spa%x0 = mean_elongation_moon_sun(spa%jce)
!   x(TERM_X0+1) = spa%x0
!   spa%x1 = mean_anomaly_sun(spa%jce)
!   x(TERM_X1+1) = spa%x1
!   spa%x2 = mean_anomaly_moon(spa%jce)
!   x(TERM_X2+1) = spa%x2
!   spa%x3 = argument_latitude_moon(spa%jce)
!   x(TERM_X3+1) = spa%x3
!   spa%x4 = ascending_longitude_moon(spa%jce)
!   x(TERM_X4+1) = spa%x4
!   CALL nutation_longitude_and_obliquity(spa%jce, x, (spa%del_psi))
!   spa%epsilon0 = ecliptic_mean_obliquity(spa%jme)
!   spa%epsilon = ecliptic_true_obliquity(spa%del_epsilon, spa%epsilon0)
!   spa%del_tau = aberration_correction(spa%r)
!   spa%lamda = apparent_sun_longitude(spa%theta, spa%del_psi, spa%del_tau)
!   spa%nu0 = greenwich_mean_sidereal_time(spa%jd, spa%jc)
!   spa%nu = greenwich_sidereal_time(spa%nu0, spa%del_psi, spa%epsilon)
!   spa%alpha = geocentric_sun_right_ascension(spa%lamda, spa%epsilon, spa%beta)
!   spa%delta = geocentric_sun_declination(spa%beta, spa%epsilon, spa%lamda)
! 
! END SUBROUTINE
!  
! ! --------------------------------------------------
! SUBROUTINE calculate_eot_and_sun_rise_transit_set(spa)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   TYPE (SPA_DATA_T) :: spa
! ! - - - local declarations - - -
!   INTEGER :: i,a%srha
!   REAL(8) :: nu,m,h0,n,h0_prime,alpha(JD_COUNT),delta(JD_COUNT),m_rts(SUN_COUNT) 
!   REAL(8) :: nu_rts(SUN_COUNT),h_rts(SUN_COUNT),alpha_prime(SUN_COUNT) 
!   REAL(8) :: delta_prime(SUN_COUNT),h_prime(SUN_COUNT)
!   TYPE (SPA_DATA_T) :: sun_rts
! 
! ! - - - begin - - -
!   h0_prime = -1*(0.26667  + spa%atmos_refract)
!   sun_rts = spa
!   m = sun_mean_longitude(spa%jme)
!   spa%eot = eot(m, spa%alpha, spa%del_psi, spa%epsilon)
!   sun_rts%second = 0
!   sun_rts%minute = sun_rts%second
!   sun_rts%hour = sun_rts%minute
!   sun_rts%timezone = 0.0
!   sun_rts%jd = julian_day(sun_rts%year, sun_rts%month, sun_rts%day,  &
!       sun_rts%hour, sun_rts%minute, sun_rts%second, sun_rts%timezone)
!   CALL calculate_geocentric_sun_right_ascension_and_declination(sun_rts)
!   nu = sun_rts%nu
!   sun_rts%delta_t = 0
!   sun_rts%jd = sun_rts%jd-1
!   DO i = 0,JD_COUNT-1
!     CALL calculate_geocentric_sun_right_ascension_and_declination(sun_rts)
!     alpha(i+1) = sun_rts%alpha
!     delta(i+1) = sun_rts%delta
!     sun_rts%jd = sun_rts%jd+1
!   END DO
!   m_rts(SUN_TRANSIT+1) = approx_sun_transit_time(alpha(JD_ZERO+1), spa%longitude, nu)
!   h0 = sun_hour_angle_at_rise_set(spa%latitude, delta(JD_ZERO+1), h0_prime)
!   IF (h0 >= 0) THEN
!     CALL approx_sun_rise_and_set(m_rts, h0)
!     DO i = 0,SUN_COUNT-1
!       nu_rts(i+1) = nu + 360.985647*m_rts(i+1)
!       n = m_rts(i+1) + spa%delta_t/86400.0
!       alpha_prime(i+1) = rts_alpha_delta_prime(alpha, n)
!       delta_prime(i+1) = rts_alpha_delta_prime(delta, n)
!       h_prime(i+1) = limit_degrees180pm(nu_rts(i+1) + spa%longitude - alpha_prime(i+1))
!       h_rts(i+1) = rts_sun_altitude(spa%latitude, delta_prime(i+1), h_prime(i+1))
!     END DO
!     spa%srha = h_prime(SUN_RISE+1)
!     spa%ssha = h_prime(SUN_SET+1)
!     spa%sta = h_rts(SUN_TRANSIT+1)
!     spa%suntransit = dayfrac_to_local_hr(m_rts(SUN_TRANSIT+1) - h_prime(SUN_TRANSIT+1) / 360.0,  &
!         spa%timezone)
!     spa%sunrise = dayfrac_to_local_hr(sun_rise_and_set(m_rts, h_rts, delta_prime,  &
!         spa%latitude, h_prime, h0_prime, SUN_RISE), spa%timezone)
!     spa%sunset = dayfrac_to_local_hr(sun_rise_and_set(m_rts, h_rts, delta_prime,  &
!         spa%latitude, h_prime, h0_prime, SUN_SET), spa%timezone)
!   ELSE
!     a%srha = spa%ssha = spa%sta = spa%suntransit = spa%sunrise = spa%sunset = -99999
!   END IF
! 
! END SUBROUTINE
!  
! ! --------------------------------------------------
! FUNCTION spa_calculate(spa)  RESULT (output_4)
! ! --------------------------------------------------
! USE SPA_1
! IMPLICIT NONE
! ! - - - arg types - - -
!   INTEGER :: output_4                                       
!   TYPE (SPA_DATA_T) :: spa
! ! - - - local declarations - - -
!   INTEGER :: result
! 
! ! - - - begin - - -
!   result = validate_inputs(spa)
!   IF (result == 0) THEN
!     spa%jd = julian_day(spa%year, spa%month, spa%day,  &
!         spa%hour, spa%minute, spa%second, spa%timezone)
!     CALL calculate_geocentric_sun_right_ascension_and_declination(spa)
!     spa%h = observer_hour_angle(spa%nu, spa%longitude, spa%alpha)
!     spa%xi = sun_equatorial_horizontal_parallax(spa%r)
!     CALL sun_right_ascension_parallax_and_topocentric_dec(spa%latitude, spa%elevation, spa%xi, spa%h, spa%delta, (spa%del_alpha))
!     spa%alpha_prime = topocentric_sun_right_ascension(spa%alpha, spa%del_alpha)
!     spa%h_prime = topocentric_local_hour_angle(spa%h, spa%del_alpha)
!     spa%e0 = topocentric_elevation_angle(spa%latitude, spa%delta_prime, spa%h_prime)
!     spa%del_e = atmospheric_refraction_correction(spa%pressure, spa%temperature,  &
!         spa%atmos_refract, spa%e0)
!     spa%e = topocentric_elevation_angle_corrected(spa%e0, spa%del_e)
!     spa%zenith = topocentric_zenith_angle(spa%e)
!     spa%azimuth180 = topocentric_azimuth_angle_neg180_180(spa%h_prime, spa%latitude,  &
!         spa%delta_prime)
!     spa%azimuth = topocentric_azimuth_angle_zero_360(spa%azimuth180)
!     IF ((spa%function == SPA_ZA_INC) .OR. (spa%function == SPA_ALL))  &
!     spa%incidence = surface_incidence_angle(spa%zenith, spa%azimuth180,  &
!         spa%azm_rotation, spa%slope)
!     IF ((spa%function == SPA_ZA_RTS) .OR. (spa%function == SPA_ALL)) CALL calculate_eot_and_sun_rise_transit_set(spa)
!   END IF
!   output_4 = result
!   RETURN
! 
! END FUNCTION
! INCLUDE 'C2F_LIB.F90'  
! !  0 Errors detected
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
