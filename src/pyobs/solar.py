"""
DEMO TO COMPUTE SOLAR ZENITH ANGLE
version 6 April 2017
by Antti Lipponen

Slightly modified by Arlindo da Silva, NASA/GSFC (October 2020)

Copyright (c) 2017 Antti Lipponen

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import numpy as np  # numerics & matrix algebra
from datetime import datetime  # dates and times


#################################################################################
# FUNCTION TO COMPUTE SOLAR AZIMUTH AND ZENITH ANGLE
# translated to Python from http://www.psa.es/sdg/sunpos.htm
#################################################################################
def solar_angles(day, dLongitude, dLatitude):
    """
        inputs: day: datetime object
                dLongitude: longitudes (scalar or Numpy array) in deg
                dLatitude: latitudes (scalar or Numpy array) in deg
        output: solar angles: (zenith,azimuth) in deg
    """
    dHours, dMinutes, dSeconds = day.hour, day.minute, day.second
    iYear, iMonth, iDay = day.year, day.month, day.day

    dEarthMeanRadius = 6371.01
    dAstronomicalUnit = 149597890

    ###################################################################
    # Calculate difference in days between the current Julian Day
    # and JD 2451545.0, which is noon 1 January 2000 Universal Time
    ###################################################################
    # Calculate time of the day in UT decimal hours
    dDecimalHours = dHours + (dMinutes + dSeconds / 60.) / 60.
    # Calculate current Julian Day
    liAux1 = int((iMonth - 14.) / 12.)
    liAux2 = int((1461. * (iYear + 4800. + liAux1)) / 4.) + int((367. * (iMonth - 2. - 12. * liAux1)) / 12.) - int((3. * int((iYear + 4900. + liAux1) / 100.)) / 4.) + iDay - 32075.
    dJulianDate = liAux2 - 0.5 + dDecimalHours / 24.
    # Calculate difference between current Julian Day and JD 2451545.0
    dElapsedJulianDays = dJulianDate - 2451545.0

    ###################################################################
    # Calculate ecliptic coordinates (ecliptic longitude and obliquity of the
    # ecliptic in radians but without limiting the angle to be less than 2*Pi
    # (i.e., the result may be greater than 2*Pi)
    ###################################################################
    dOmega = 2.1429 - 0.0010394594 * dElapsedJulianDays
    dMeanLongitude = 4.8950630 + 0.017202791698 * dElapsedJulianDays  # Radians
    dMeanAnomaly = 6.2400600 + 0.0172019699 * dElapsedJulianDays
    dEclipticLongitude = dMeanLongitude + 0.03341607 * np.sin(dMeanAnomaly) + 0.00034894 * np.sin(2. * dMeanAnomaly) - 0.0001134 - 0.0000203 * np.sin(dOmega)
    dEclipticObliquity = 0.4090928 - 6.2140e-9 * dElapsedJulianDays + 0.0000396 * np.cos(dOmega)

    ###################################################################
    # Calculate celestial coordinates ( right ascension and declination ) in radians
    # but without limiting the angle to be less than 2*Pi (i.e., the result may be
    # greater than 2*Pi)
    ###################################################################
    dSin_EclipticLongitude = np.sin(dEclipticLongitude)
    dY = np.cos(dEclipticObliquity) * dSin_EclipticLongitude
    dX = np.cos(dEclipticLongitude)
    dRightAscension = np.arctan2(dY, dX)
    if dRightAscension < 0.0:
        dRightAscension = dRightAscension + 2.0 * np.pi
    dDeclination = np.arcsin(np.sin(dEclipticObliquity) * dSin_EclipticLongitude)

    ###################################################################
    # Calculate local coordinates ( azimuth and zenith angle ) in degrees
    ###################################################################
    dGreenwichMeanSiderealTime = 6.6974243242 + 0.0657098283 * dElapsedJulianDays + dDecimalHours
    dLocalMeanSiderealTime = (dGreenwichMeanSiderealTime * 15. + dLongitude) * (np.pi / 180.)
    dHourAngle = dLocalMeanSiderealTime - dRightAscension
    dLatitudeInRadians = dLatitude * (np.pi / 180.)
    dCos_Latitude = np.cos(dLatitudeInRadians)
    dSin_Latitude = np.sin(dLatitudeInRadians)
    dCos_HourAngle = np.cos(dHourAngle)
    dZenithAngle = (np.arccos(dCos_Latitude * dCos_HourAngle * np.cos(dDeclination) + np.sin(dDeclination) * dSin_Latitude))
    dY = -np.sin(dHourAngle)
    dX = np.tan(dDeclination) * dCos_Latitude - dSin_Latitude * dCos_HourAngle
    dAzimuth = np.arctan2(dY, dX)
    dAzimuth[dAzimuth < 0.0] = dAzimuth[dAzimuth < 0.0] + 2.0 * np.pi
    dAzimuth = dAzimuth / (np.pi / 180.)
    # Parallax Correction
    dParallax = (dEarthMeanRadius / dAstronomicalUnit) * np.sin(dZenithAngle)
    dZenithAngle = (dZenithAngle + dParallax) / (np.pi / 180.)

    return dAzimuth - 180.0, dZenithAngle

def get_cosz (tyme, dLongitude, dLatitude):
    """
        inputs: tyme: datetime object
                dLongitude: longitudes (scalar or Numpy array) in deg
                dLatitude: latitudes (scalar or Numpy array) in deg
        output: cosize of solar zenith angle
    """
    dHours, dMinutes, dSeconds = tyme.hour, tyme.minute, tyme.second
    iYear, iMonth, iDay = tyme.year, tyme.month, tyme.day

    dEarthMeanRadius = 6371.01
    dAstronomicalUnit = 149597890

    ###################################################################
    # Calculate difference in days between the current Julian Day
    # and JD 2451545.0, which is noon 1 January 2000 Universal Time
    ###################################################################
    # Calculate time of the day in UT decimal hours
    dDecimalHours = dHours + (dMinutes + dSeconds / 60.) / 60.
    # Calculate current Julian Day
    liAux1 = int((iMonth - 14.) / 12.)
    liAux2 = int((1461. * (iYear + 4800. + liAux1)) / 4.) + int((367. * (iMonth - 2. - 12. * liAux1)) / 12.) - int((3. * int((iYear + 4900. + liAux1) / 100.)) / 4.) + iDay - 32075.
    dJulianDate = liAux2 - 0.5 + dDecimalHours / 24.
    # Calculate difference between current Julian Day and JD 2451545.0
    dElapsedJulianDays = dJulianDate - 2451545.0

    ###################################################################
    # Calculate ecliptic coordinates (ecliptic longitude and obliquity of the
    # ecliptic in radians but without limiting the angle to be less than 2*Pi
    # (i.e., the result may be greater than 2*Pi)
    ###################################################################
    dOmega = 2.1429 - 0.0010394594 * dElapsedJulianDays
    dMeanLongitude = 4.8950630 + 0.017202791698 * dElapsedJulianDays  # Radians
    dMeanAnomaly = 6.2400600 + 0.0172019699 * dElapsedJulianDays
    dEclipticLongitude = dMeanLongitude + 0.03341607 * np.sin(dMeanAnomaly) + 0.00034894 * np.sin(2. * dMeanAnomaly) - 0.0001134 - 0.0000203 * np.sin(dOmega)
    dEclipticObliquity = 0.4090928 - 6.2140e-9 * dElapsedJulianDays + 0.0000396 * np.cos(dOmega)

    ###################################################################
    # Calculate celestial coordinates ( right ascension and declination ) in radians
    # but without limiting the angle to be less than 2*Pi (i.e., the result may be
    # greater than 2*Pi)
    ###################################################################
    dSin_EclipticLongitude = np.sin(dEclipticLongitude)
    dY = np.cos(dEclipticObliquity) * dSin_EclipticLongitude
    dX = np.cos(dEclipticLongitude)
    dRightAscension = np.arctan2(dY, dX)
    if dRightAscension < 0.0:
        dRightAscension = dRightAscension + 2.0 * np.pi
    dDeclination = np.arcsin(np.sin(dEclipticObliquity) * dSin_EclipticLongitude)

    ###################################################################
    # Calculate local coordinates ( azimuth and zenith angle ) in degrees
    ###################################################################
    dGreenwichMeanSiderealTime = 6.6974243242 + 0.0657098283 * dElapsedJulianDays + dDecimalHours
    dLocalMeanSiderealTime = (dGreenwichMeanSiderealTime * 15. + dLongitude) * (np.pi / 180.)
    dHourAngle = dLocalMeanSiderealTime - dRightAscension
    dLatitudeInRadians = dLatitude * (np.pi / 180.)
    dCos_Latitude = np.cos(dLatitudeInRadians)
    dSin_Latitude = np.sin(dLatitudeInRadians)
    dCos_HourAngle = np.cos(dHourAngle)
    dZenithAngle = (np.arccos(dCos_Latitude * dCos_HourAngle * np.cos(dDeclination) + np.sin(dDeclination) * dSin_Latitude))
  
    # Parallax Correction
    dParallax = (dEarthMeanRadius / dAstronomicalUnit) * np.sin(dZenithAngle)
    dZenithAngle = (dZenithAngle + dParallax)
    

    return np.cos(dZenithAngle)

def test_angles():
    
#################################################################################
# COMPUTE ZENITH ANGLES AND AZIMUTHS
#################################################################################

    # time
    t = datetime(2017, 4, 6, 9, 30, 0)  # 6th April 2017 09:30:00 UTC

    # grid dimensions
    Nlats = 90
    Nlons = 90

    # coordinates in grid
    lat, lon = np.linspace(-90.0, 90.0, Nlats), \
               np.linspace(-180.0, 180.0, Nlons)  

    # compute solar zenith angle and azimuth (be careful with azimuth: I haven't checked this at all)
    saz, sza = solar_angles(t, lon.ravel(), lat.ravel())

    print "Solar Zenith Angle: ", np.cos(np.pi*sza/180.)

#.................................................................................
if __name__ == "__main__":

    from grads import GrADS

    rad = '/home/adasilva/opendap/m2/opendap/tavg1_2d_rad_Nx'
    ga = GrADS(Window=False,Echo=False)

    fh = ga.open(rad)
    solar = ga.exp('swtnt') # (ny,nx)

    lon = solar.grid.lon
    lat = solar.grid.lat
    tyme = solar.grid.tyme

    Lon = ga.exp('lon')
    Lat = ga.exp('lat')

    cosz = get_cosz(tyme[0],Lon,Lat)

    I = solar<=0  # night
    J = solar>0   # day

    cosz_ngt = cosz[I]
    cosz_day = cosz[J]

    print "cosz night: ", cosz_ngt.min(), cosz_ngt.max()
    print "cosz day: ",   cosz_day.min(), cosz_day.max()
    


    

    
