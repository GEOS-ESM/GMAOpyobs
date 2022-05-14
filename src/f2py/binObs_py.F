      subroutine screenObs2D(obs,sObs,im,jm,missing)
      implicit NONE
      integer, intent(in)  :: im, jm
      real, intent(in)     :: obs(im,jm)
      real, intent(in)     :: missing
      real, intent(out)    :: sObs(im,jm)
!
! Screen obs. Assigns missing  values to pixel if nearby obs is MISSING.
!
!
      integer              :: i, j, ii, jj
      integer              :: iis, iie
      integer              :: jjs, jje

      sObs = obs

      do i = 1, im
         do j = 1, jm
            if ( abs(obs(i,j)-missing) .gt. 0.01*abs(missing) ) then
               ! Look at all the pixels around it
               if ( i .eq. 1 ) then
                  iis = i
               else
                  iis = i - 1
               end if
               if ( i .eq. im ) then
                  iie = im
               else
                  iie = i + 1
               end if
               if ( j .eq. 1 ) then
                  jjs = j
               else
                  jjs = j - 1
               end if
               if ( j .eq. jm ) then
                  jje = j
               else
                  jje = j + 1
               end if

               do ii = iis, iie
                  do jj = jjs, jje
                     if ( abs(obs(ii,jj)-missing) .le.
     &                    0.01*abs(missing) ) then
                       sObs(i,j) = missing
                    end if
                  end do
               end do
            end if
         end do
      end do

      end

!-----------
      subroutine binObs2D(lon,lat,obs,nobs,gObs,im,jm,missing)
!vb      use MAPL_BaseMod, only: MAPL_GetHorzijIndex
      implicit NONE
      integer, intent(in) :: im, jm, nobs
      real, intent(in)  :: lon(nobs)
      real, intent(in)  :: lat(nobs)
      real, intent(in)  :: obs(nobs)
      real, intent(in)  :: missing
      real, intent(out) :: gObs(im,jm)
!
! Bins obs. It assumes a global 2D GEOS-5 A-Grid:
!
!  Longitudes in [-180,10)
!  Latitudes  in [-90,90]
!
!

      integer :: i, j, n
      real :: dLon, dLat, xLon
      real :: nsample(im,jm)
      integer :: ics(nobs),jcs(nobs)
      real, parameter :: radToDeg = 57.2957795

      gObs = 0.0
      nSample = 0.0
      if (jm == 6*im) then
         print*, "warning, function under dev"
!vb         call MAPL_GetHorzIJIndex(im,jm,nobs,ics,jcs,lon=lon
!vb     &   /radToDeg,lat=lat/radToDeg)

!vb         do n= 1,nobs
!vb            if (abs(obs(n)-missing) > 0.01*abs(missing)) then
!vb            gObs(ics(n),jcs(n)) = gObs(ics(n),jcs(n)) + obs(n)
!vb            nSample(ics(n),jcs(n)) = nSample(ics(n),jcs(n)) + 1.0
!vb            end if
!vb         enddo

      else

         dLon = 360. / im
         dLat = 180. / ( jm - 1.)

         do n = 1, nobs

            if ( (abs(lon(n))>180.) .OR. (abs(lat(n))>90.) ) cycle

            xLon = lon(n)
            if ( xLon >= 180. ) xLon = xLon - 360.
    
            i = 1 + nint((xlon + 180.  ) / dLon)
            j = 1 + nint((lat(n) + 90. ) / dLat)

            if ( i>im ) i = i - im
            if ( i<1  ) i = i + im
            
            if (abs(obs(n)-missing) > 0.01*abs(missing)) then
            gObs(i,j) = gObs(i,j) + obs(n)
            nSample(i,j) = nSample(i,j) + 1.0
            end if
         end do

      end if

      where ( nsample > 0 ) 
              gObs = gObs / nSample
      elsewhere
              gObs = missing
      end where

      end

!-----------

      subroutine binObs3D(lon,lat,obs,nobs,gObs,im,jm,km,missing)
!vb      use MAPL_BaseMod, only: MAPL_GetHorzijIndex
      implicit NONE
      integer, intent(in) :: im, jm, km,nobs
      real, intent(in)  :: lon(nobs)
      real, intent(in)  :: lat(nobs)
      real, intent(in)  :: obs(nobs,km)
      real, intent(in)  :: missing
      real, intent(out) :: gObs(im,jm,km)
!
! Bins obs. It assumes a global 3D GEOS-5 A-Grid:
!
!  Longitudes in [-180,10]
!  Latitudes  in [-90,90]
!
!

      integer :: i, j, k, n
      real :: dLon, dLat, xLon
      real :: nSample(im,jm)
      integer :: ics(nobs),jcs(nobs)
      real, parameter :: radToDeg = 57.2957795

      dLon = 360. / im
      dLat = 180. / ( jm - 1.)

      gObs = 0.0

      do k = 1, km
         nSample = 0.0

         if (jm == 6*im)then
               print*, "warning, function under dev"
!vb               call MAPL_GetHorzijIndex(im,jm,nobs,ics,jcs,
!vb     &         lon=lon/radToDeg,lat=lat/radToDeg)

!vb               if (abs(obs(n,k)-missing) > 0.01*abs(missing)) then
!vb               gObs(ics(n),jcs(n),k) = gObs(ics(n),jcs(n),k) + obs(n,k)
!vb               nSample(ics(n),jcs(n)) = nSample(ics(n),jcs(n)) + 1.0
!vb               end if

         else
         
            do n = 1,nobs
                    
               if ( (abs(lon(n))>180.) .OR. (abs(lat(n))>90.) ) cycle

               xLon = lon(n)
               if ( xLon >= 180. ) xLon = xLon - 360.
               
               i = 1 + nint((xlon + 180.  ) / dLon)
               j = 1 + nint((lat(n) + 90. ) / dLat)

               if ( i>im ) i = i - im
               if ( i<1  ) i = i + im
          
   !           if (abs(obs(n,k)-missing) > 0.1) then
               if (abs(obs(n,k)-missing) > 0.01*abs(missing)) then
                  gObs(i,j,k) = gObs(i,j,k) + obs(n,k)
                  nSample(i,j) = nSample(i,j) + 1.0
               end if

            end do

         end if

         where ( nsample > 0 ) 
            gObs(:,:,k) = gObs(:,:,k) / nSample
         elsewhere
            gObs(:,:,k) = missing
         end where
      
      end do

      end
!-------------
      subroutine binRms2D(lon,lat,obs,nobs,gObs,im,jm,missing)
!vb      use MAPL_BaseMod, only: MAPL_GetHorzijIndex
      implicit NONE
      integer, intent(in) :: im, jm, nobs
      real, intent(in)  :: lon(nobs)
      real, intent(in)  :: lat(nobs)
      real, intent(in)  :: obs(nobs)
      real, intent(in)  :: missing
      real, intent(out) :: gObs(im,jm)
!
! Bins RMS.  It assumes a global GEOS-5 A-Grid:
!
!  Longitudes in [-180,10)
!  Latitudes  in [-90,90]
!
!

      integer :: i, j, n
      real :: dLon, dLat, xLon
      real :: nsample(im,jm)
      integer :: ics(nobs),jcs(nobs)
      real, parameter :: radToDeg = 57.2957795

      dLon = 360. / im
      dLat = 180. / ( jm - 1.)

      gObs = 0.0
      nSample = 0.0
      if (jm == 6*im) then
         print*, "warning, function under dev"
!vb         call MAPL_GetHorzIJIndex(im,jm,nobs,ics,jcs,
!vb     &   lon=lon/radToDeg,lat=lat/radToDeg)

!vb         do n= 1,nobs
!vb            if (abs(obs(n)-missing) > 0.01*abs(missing)) then
!vb            gObs(ics(n),jcs(n)) = gObs(ics(n),jcs(n)) + obs(n)
!vb            nSample(ics(n),jcs(n)) = nSample(ics(n),jcs(n)) + 1.0
!vb            end if
!vb         enddo

      else

         do n = 1, nobs

            if ( (abs(lon(n))>180.) .OR. (abs(lat(n))>90.) ) cycle

            xLon = lon(n)
            if ( xLon >= 180. ) xLon = xLon - 360.
    
            i = 1 + nint((xlon + 180.  ) / dLon)
            j = 1 + nint((lat(n) + 90. ) / dLat)

            if ( i>im ) i = i - im
            if ( i<1  ) i = i + im
            
   !        if (abs(obs(n)-missing) > 0.1) then
            if (abs(obs(n)-missing) > 0.01*abs(missing)) then
            gObs(i,j) = gObs(i,j) + obs(n)*obs(n)
            nSample(i,j) = nSample(i,j) + 1.0
            end if
         end do

      end if

      where ( nsample > 0 ) 
              gObs = sqrt(abs(gObs / nSample))
      elsewhere
              gObs = missing
      end where

      end

!-----------

      subroutine binRms3D(lon,lat,obs,nobs,gObs,im,jm,km,missing)
!vb      use MAPL_BaseMod, only: MAPL_GetHorzijIndex
      implicit NONE
      integer, intent(in) :: im, jm, km,nobs
      real, intent(in)  :: lon(nobs)
      real, intent(in)  :: lat(nobs)
      real, intent(in)  :: obs(nobs,km)
      real, intent(in)  :: missing
      real, intent(out) :: gObs(im,jm,km)
!
! Bins RMS. It assumes a global 3D GEOS-5 A-Grid:
!
!  Longitudes in [-180,10]
!  Latitudes  in [-90,90]
!
!

      integer :: i, j, k, n
      real :: dLon, dLat, xLon
      real :: nSample(im,jm)
      integer :: ics(nobs),jcs(nobs)
      real, parameter :: radToDeg = 57.2957795

      dLon = 360. / im
      dLat = 180. / ( jm - 1.)

      gObs = 0.0

      do k = 1, km
         nSample = 0.0
         if (jm == 6*im)then
               print*, "warning, function under dev"
!vb               call MAPL_GetHorzijIndex(im,jm,nobs,ics,jcs,
!vb     &         lon=lon/radToDeg,lat=lat/radToDeg)

!vb               if (abs(obs(n,k)-missing) > 0.01*abs(missing)) then
!vb               gObs(ics(n),jcs(n),k) = gObs(ics(n),jcs(n),k) + obs(n,k)
!vb               nSample(ics(n),jcs(n)) = nSample(ics(n),jcs(n)) + 1.0
!vb               end if
         else 
            do n = 1,nobs
                    
               if ( (abs(lon(n))>180.) .OR. (abs(lat(n))>90.) ) cycle

               xLon = lon(n)
               if ( xLon >= 180. ) xLon = xLon - 360.
               
               i = 1 + nint((xlon + 180.  ) / dLon)
               j = 1 + nint((lat(n) + 90. ) / dLat)

               if ( i>im ) i = i - im
               if ( i<1  ) i = i + im
          
   !           if (abs(obs(n,k)-missing) > 0.1) then
               if (abs(obs(n,k)-missing) > 0.01*abs(missing)) then
                  gObs(i,j,k) = gObs(i,j,k) + obs(n,k)*obs(n,k)
                  nSample(i,j) = nSample(i,j) + 1.0
               end if

            end do

         end if

         where ( nsample > 0 ) 
            gObs(:,:,k) = sqrt(abs(gObs(:,:,k) / nSample))
         elsewhere
            gObs(:,:,k) = missing
         end where
      
      end do

      end
!--------------
      subroutine binObs3DP(lon,lat,obs,pobs,pe,nobs,nlev,
     & gObs,im,jm,km,missing)
!vb      use MAPL_BaseMod, only: MAPL_GetHorzijIndex
      implicit NONE
      integer, intent(in) :: im, jm, km,nobs
      real, intent(in)  :: lon(nobs)
      real, intent(in)  :: lat(nobs)
      real, intent(in)  :: obs(nobs,nlev)
      real, intent(in)  :: pobs(nobs,nlev)
      real, intent(in)  :: pe(im,jm,km)
      real, intent(in)  :: missing
      real, intent(out) :: gObs(im,jm,km-1)
!



!
! Bins obs. It assumes a global 3D GEOS-5 A-Grid:
!
!  Longitudes in [-180,10]
!  Latitudes  in [-90,90]
!
! 
      integer :: i, j, k, n, P, nlev, m
      real :: dLon, dLat, xLon
      real :: nSample(im,jm,km)
      
      dLon = 360. / im
      dLat = 180. / ( jm - 1.)

      
      gObs = 0.0
      m = 0
      do k = 1, km-1
         nSample = 0.0
         do n = 1,nobs
                    
            if ( (abs(lon(n))>180.) .OR. (abs(lat(n))>90.) ) cycle

            xLon = lon(n)
            if ( xLon >= 180. ) xLon = xLon - 360.
            
            i = 1 + nint((xlon + 180.  ) / dLon)
            j = 1 + nint((lat(n) + 90. ) / dLat)

            if ( i>im ) i = i - im
            if ( i<1  ) i = i + im
            
            do P = 1 , nlev

               if (k == km-1) then    ! if CALIPSO pressure > GEOS-5 P for edge 73 -> 
                                      ! put CALIPSO values in level 72 
                  if (pobs(n,P) > pe(i,j,k)) then
                   
                    if (abs(obs(n,P)-missing) > 0.01*abs(missing)) then
                                            
                      gObs(i,j,k) = gObs(i,j,k) + obs(n,P)
                      nSample(i,j,k) = nSample(i,j,k) + 1.0
                    end if
                  end if
               
               else if ( (pobs(n,P) > pe(i,j,k)).AND.
     &            (pobs(n,P) < pe(i,j,k+1) ) )then
                                      
                   if (abs(obs(n,P)-missing) > 0.01*abs(missing)) then

                      gObs(i,j,k) = gObs(i,j,k) + obs(n,P)
                      nSample(i,j,k) = nSample(i,j,k) + 1.0
                   end if
               end if
               
            end do
         end do
         
         where ( nsample(:,:,k) > 0 )            
            gObs(:,:,k) = gObs(:,:,k) / nSample(:,:,k)
         elsewhere
            gObs(:,:,k) = missing
         end where
      
      end do

      end

!--------------
      subroutine binObs3DH(lon,lat,obs,hobs,he,nobs,nlev,
     & im,jm,km,missing, gObs)
      implicit NONE
      integer, intent(in) :: im, jm, km,nobs
      real, intent(in)  :: lon(nobs)
      real, intent(in)  :: lat(nobs)
      real, intent(in)  :: obs(nobs,nlev)
      real, intent(in)  :: hobs(nlev)
      real, intent(in)  :: he(km,nobs)
      real, intent(in)  :: missing
      real, intent(out) :: gObs(im,jm,km)
!
!
! Bins obs. It assumes a global 3D GEOS-5 A-Grid:
!
!  Longitudes in [-180,10]
!  Latitudes  in [-90,90]
!
! 
      integer :: i, j, k, n, h, nlev, m,t
      real :: dLon, dLat, xLon
      real :: nSample(im,jm,km)
      
      dLon = 360. / im
      dLat = 180. / ( jm - 1.)

      
      gObs = 0.0
      m = 0
      do k = 1, km
         nSample = 0.0
         do n = 1,nobs
                    
            if ( (abs(lon(n))>180.) .OR. (abs(lat(n))>90.) ) cycle

            xLon = lon(n)
            if ( xLon >= 180. ) xLon = xLon - 360.
            
            i = 1 + nint((xLon + 180.  ) / dLon)
            j = 1 + nint((lat(n) + 90. ) / dLat)
            
            if ( i>im ) i = i - im
            if ( i<1  ) i = i + im
            

            do h = 1 , nlev                 
              if (k == km) then 
                 if  (hobs(h) < he(km,n)) then
                    
                    if (abs(obs(n,h)-missing) > 0.01*abs(missing)) then
                      gObs(i,j,k) = gObs(i,j,k) + obs(n,h)
                      
                      nSample(i,j,k) = nSample(i,j,k) + 1.0 
                      
                    end if
                 end if  
              else if ( (hobs(h) < he(k,n)).AND.
     &            (hobs(h) > he(k+1,n) )) then
                                 
                   if (abs(obs(n,h)-missing) > 0.01*abs(missing)) then

                      gObs(i,j,k) = gObs(i,j,k) + obs(n,h)
                      nSample(i,j,k) = nSample(i,j,k) + 1.0   
                              
                   end if
               end if
              
               
            end do
         end do
         
         where ( nsample(:,:,k) > 0 )            
            gObs(:,:,k) = gObs(:,:,k) / nSample(:,:,k)
         elsewhere
            gObs(:,:,k) = missing
         end where
      
      end do
      
      end
!--------------

      subroutine binAreas(lon,lat,area,nobs,gArea,im,jm)
!vb      use MAPL_BaseMod, only: MAPL_GetHorzijIndex
      implicit NONE
      integer, intent(in) :: im, jm, nobs
      real, intent(in)  :: lon(nobs)
      real, intent(in)  :: lat(nobs)
      real, intent(in)  :: area(nobs)
      real, intent(out) :: gArea(im,jm)
!
! Bins areas of no-fire pixels. It assumes a global GEOS-5 A-Grid:
!
!  Longitudes in [-180,180)
!  Latitudes  in [-90,90]
!
!
      integer :: i, j, n
      real :: dLon, dLat, xLon
      integer ics(nobs),jcs(nobs)
      real, parameter :: radToDeg = 57.2957795

      gArea = 0.0

      if (jm == 6*im) then
         print*, "warning, function under dev"
!vb         call MAPL_GetHorzIJIndex(im,jm,nobs,ics,jcs,
!vb     &   lon=lon/radToDeg,lat=lat/radToDeg)

!vb         do n= 1,nobs
!vb            gArea(ics(n),jcs(n)) = gArea(ics(n),jcs(n)) + area(n)
!vb         enddo

      else

         dLon = 360. / im
         dLat = 180. / ( jm - 1.)

         do n = 1, nobs

            xLon = lon(n)
            if ( xLon >= 180. ) xLon = xLon - 360.
    
            i = 1 + nint((xlon + 180.  ) / dLon)
            j = 1 + nint((lat(n) + 90. ) / dLat)

            if ( i>im ) i = i - im
            if ( i<1  ) i = i + im

            gArea(i,j) = gArea(i,j) + area(n)

         end do
      end if

      end

!--------------
      subroutine binAreasNR(lon,lat,area,nobs,gArea,im,jm)
      implicit NONE
      integer, intent(in) :: im, jm, nobs
      real, intent(in)  :: lon(nobs)
      real, intent(in)  :: lat(nobs)
      real, intent(in)  :: area(nobs)
      real, intent(out) :: gArea(im,jm)
!
! Bins areas of no-fire pixels. It assumes a uniform DEPE grid.
!
!  Longitudes in (-180,180)
!  Latitudes  in (-90,90)
!
!
      integer :: i, j, n
      real :: dLon, dLat, xLon

      dLon = 360. / im
      dLat = 180. / jm

      gArea = 0.0
      do n = 1, nobs

         xLon = lon(n)
         if ( xLon >= 180. ) xLon = xLon - 360.
 
         i = nint((xLon   + 180.0) / dLon + 0.5)
         j = nint((lat(n) +  90.0) / dLat + 0.5)

         if ( i>im ) then
             i = i - im
!            print *, '[w] adjustment i>im', xLon, lat(n)
         end if
         if ( i<1  ) then
             i = i + im
!            print *, '[w] adjustment i<1', xLon, lat(n)
         end if

         if ( j>jm ) then
             j = jm
!            print *, '[w] adjustment j>jm', xLon, lat(n)
         end if
         if ( j<1  ) then
             j = 1
!            print *, '[w] adjustment j<1', xLon, lat(n)
         end if

         gArea(i,j) = gArea(i,j) + area(n)

      end do

      end

!..................................................................

      subroutine findBox(lon,lat,nobs,im,jm,gridbox)
      implicit NONE
      integer, intent(in)  :: im, jm, nobs
      real, intent(in)     :: lon(nobs)
      real, intent(in)     :: lat(nobs)
      integer, intent(out) :: gridbox(nobs)
!
!     Returns gridbox number for each lon,lat on input.
!
      print *, 'findBox: not implemented yet.'
      gridbox = 0

      end

!..................................................................

      subroutine decimateSwath(a, im, jm, b, ih, jh, idel, missing, rc )
      implicit NONE
      integer, intent(in) :: im     
      integer, intent(in) :: jm     
      integer, intent(in) :: ih     
      integer, intent(in) :: jh     
      integer, intent(in) :: idel   
                                    

      real, intent(in)  :: b(ih,jh) 
      real, intent(in)  :: missing  
      
      real, intent(out) :: a(im,jm) 
      integer, intent(out) :: rc    
      
!
!     Averages high res array *b* producing low res array *a*.
!

!                           ---

      integer :: idel_, jdel_, i, j, i_, j_, i1, i2, j1, j2, n
      real*8  :: accum

#define VALID(x) (abs(x-missing) > 0.01*abs(missing))

      rc = 0

!     Consistency check
!     -----------------
      idel_ = nint(float(ih)/float(im))
      jdel_ = nint(float(jh)/float(jm))
      if ( idel_ /= idel .OR. jdel_ /= idel ) then
         rc = 1
         return
      end if
      if ( mod(idel,2) == 0 ) then
         rc = 2
         return
      end if
      
!     Average *good* values of high res data
!     --------------------------------------
      do j = 1, jm

       j1 = (j-1)*idel + 1
       j2 = min(jh,j1+idel-1)

       do i = 1, im

          i1 = (i-1)*idel + 1
          i2 = min(ih,i1+idel-1)

          accum = 0.0
          n = 0
          do j_ = j1, j2
             do i_ = i1, i2
                if ( VALID(b(i_,j_)) ) then
                   accum = accum + b(i_,j_)
                   n = n + 1
                end if
             end do
          end do

          if ( n == 0 ) then
             a(i,j) = missing
          else
             a(i,j) = accum / n
          end if

       end do
      end do
      return
      end 

!..................................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Interpxy3d --- Interpolates 3D field to observation locations
!
! !INTERFACE:
!
      subroutine InterpXY3d (gField,glon,glat,lon,lat,
     &                       periodic,nobs,im,jm,km,oField)
     
! !USES:
!
      Implicit NONE

! !INPUT PARAMETERS:
!

      integer, intent(in)        :: nobs      ! Number of observations

                                              ! longitude in degrees [-180,+180]
      real,    intent(in)        :: lon(nobs) 
                                              ! latitude  in degrees [-90,90]
      real,    intent(in)        :: lat(nobs) 


      integer, intent(in)        :: im        ! zonal dimension
      integer, intent(in)        :: jm        ! meridional dimension
      integer, intent(in)        :: km        ! vertical dimension: 
                                              ! = 1 for 2D fields
                                              ! = km for mid-layer fields
                                              ! = km+1 for edge fields

      integer, intent(in)        :: periodic
!f2py integer  intent(in), optional :: periodic = 0
      real,    intent(in)        :: glon(im)
      real,    intent(in)        :: glat(jm)

      real,    intent(in)        :: gField(im,jm,km) 

! !OUTPUT PARAMETERS:
!
                                              ! Interpolated profile
      real,    intent(out)       :: oField(km,nobs)

! !DESCRIPTION: This routine interpolates gridded model fields to observation
!               locations. This routine implements only the horizontal
!  interpolation.
!
!  IMPORTANT:   The input lon coordinates must be in [lon_0, lon_0+360].
!               The input field cannot have any UNDEFs.
!               For now, grid assumed to be global and zonally periodic
!               Latitudes must be in the range [-90,90]
!               Observation and grid must have the same reference, that is,
!               if input lon is in the range [-180,180] then grid must be 
!               in the same range.
!
! !SEE ALSO:
!
!              Module m_insitu which uses the same linear interpolation algorithm.
!
!
! !REVISION HISTORY:
!
!  10feb2010  da Silva  Simplified m_interp routine for profile interpolation.
!
!EOP
!-------------------------------------------------------------------------
 
      character(len=*), parameter :: myname = 'InterpXY_'

! Local
      integer i, j, k, nob
      real    o_lon, o_lat
      real    m_dlon, m_dlat
      real    alfa, beta

      real a11(km)         !W-S
      real a12(km)         !W-N
      real a21(km)         !E-S
      real a22(km)         !E-N
      real a00(km)         !temp

      integer i1, i2

!                         ------

#if 0
      print *, 'Dimensions: ', im, jm, km, nobs
      print *, '       Lon: ', minval(glon), maxval(glon)
      print *, '       Lat: ', minval(glat), maxval(glat)
      print *, '       lon: ', minval(lon), maxval(lon)
      print *, '       lat: ', minval(lat), maxval(lat)
#endif

      if ( nobs .eq. 0 ) return    ! nothing to do (keep this)

      if ( periodic == 1 ) then
           m_dlon = float(im) / 360.     ! assume global & periodic
           m_dlat = float(jm-1) / 180.
      else                              ! regional 
           m_dlon = float(im-1) / ( glon(im) - glon(1) )
           m_dlat = float(jm-1) / ( glat(jm) - glat(1) )
      endif

!     Loop over observations
!     ----------------------
      do nob = 1, nobs

!       Longitude
!       ---------
        o_lon = 1. + (lon(nob)-glon(1)) * m_dlon
        if ( periodic == 1 ) then
           i   = min(im, int( o_lon ))
           alfa  = o_lon - i
           if(i .eq. im) then
              i1 = im
              i2 = 1
           else
              i1 = i
              i2 = i + 1
           endif
        else
           i   = min(im-1, int( o_lon ))
           alfa  = o_lon - i
           i1 = i
           i2 = i + 1
        end if

!       Latitude
!       --------
        o_lat = 1. + (lat(nob) - glat(1)) * m_dlat
        j   = min( jm-1, int( o_lat ) )
        beta  = o_lat - j

        a11 = gField(i1,j,  :)
        a21 = gField(i2,j,  :)
        a12 = gField(i1,j+1,:)
        a22 = gField(i2,j+1,:)
        a00 = a11 + alfa * ( a21 - a11 )

        oField(:,nob) = a00 + beta * ( a12 + alfa * ( a22 - a12 ) - a00)
     
      end do ! loop over obs

      return
  
      end subroutine InterpXY3d
