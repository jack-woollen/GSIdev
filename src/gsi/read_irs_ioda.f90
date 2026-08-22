!---------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------
      module read_irs_data

      integer(8),allocatable,dimension(:)  :: dateTime
      integer(4),allocatable,dimension(:)  :: dwellNumber
      integer(4),allocatable,dimension(:)  :: sensorChannelNumber
      integer(4),allocatable,dimension(:)  :: spatialSampleQualityLw
      integer(4),allocatable,dimension(:)  :: spatialSampleQualityMw

      real(4),allocatable,dimension(:)  :: latitude
      real(4),allocatable,dimension(:)  :: longitude
      real(4),allocatable,dimension(:)  :: sensorAzimuthAngle
      real(4),allocatable,dimension(:)  :: sensorZenithAngle
      real(4),allocatable,dimension(:)  :: solarAzimuthAngle
      real(4),allocatable,dimension(:)  :: solarZenithAngle
      real(4),allocatable,dimension(:)  :: cloudSignal
      real(4),allocatable,dimension(:)  :: cloudFraction   
      real(4),allocatable,dimension(:)  :: sensorCentralWavenumber
      real(4),allocatable,dimension(:,:):: spectralRadiance       

      character(255) gfname

      end module
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine read_irs_ioda(nobs,channels)

      use read_irs_data

      integer nobs,channels,pscores

      character(50)  group,var
      integer        nf_open
      integer        nf_inq_ncid
      integer        nf_inq_varid
      integer        nf_get_var
      integer        ncid,grpid,varid

      nobs = -1 ! for bad return

! get the array dimensions

      write(6,*) 'read_ioda: ioda file=',trim(gfname)

      call check(nf_open(trim(gfname),0,ncid) )
      call getdim(ncid,"Location",nobs)
      call getdim(ncid,"Channel",channels)
      print*,"nobs=",nobs

! allocate the data arrays

      allocate(dateTime(nobs))
      allocate(dwellNumber(nobs))
      allocate(latitude(nobs))
      allocate(longitude(nobs))
      allocate(sensorAzimuthAngle(nobs))
      allocate(sensorZenithAngle(nobs))
      allocate(solarAzimuthAngle(nobs))
      allocate(solarZenithAngle(nobs))
      allocate(cloudSignal(nobs))
      allocate(cloudFraction(nobs))
      allocate(spatialSampleQualityLw(nobs))
      allocate(spatialSampleQualityMw(nobs))
      allocate(sensorChannelNumber(channels))
      allocate(sensorCentralWavenumber(channels))
      allocate(spectralRadiance(channels,nobs))

! read the data

      group = "MetaData"; call check(NF_INQ_NCID(ncid, group, grpid))

      var="dateTime"
      call check( nf_inq_varid(grpid,var,varid) )
      call check( nf_get_var(grpid,varid,dateTime) )

      var="dwellNumber"
      call check( nf_inq_varid(grpid,var,varid) )
      call check( nf_get_var(grpid,varid,dwellNumber) )

      var="latitude"
      call check( nf_inq_varid(grpid,var,varid) )
      call check( nf_get_var(grpid,varid,latitude) )

      var="longitude"
      call check( nf_inq_varid(grpid,var,varid) )
      call check( nf_get_var(grpid,varid,longitude) )

      var="sensorAzimuthAngle"
      call check( nf_inq_varid(grpid,var,varid) )
      call check( nf_get_var(grpid,varid,sensorAzimuthAngle) )

      var="sensorZenithAngle"
      call check( nf_inq_varid(grpid,var,varid) )
      call check( nf_get_var(grpid,varid,sensorZenithAngle) )

      var="solarAzimuthAngle"
      call check( nf_inq_varid(grpid,var,varid) )
      call check( nf_get_var(grpid,varid,solarAzimuthAngle) )

      var="solarZenithAngle"
      call check( nf_inq_varid(grpid,var,varid) )
      call check( nf_get_var(grpid,varid,solarZenithAngle) )

      var="cloudSignal"
      call check( nf_inq_varid(grpid,var,varid) )
      call check( nf_get_var(grpid,varid,cloudSignal) )

      var="cloudFraction"
      call check( nf_inq_varid(grpid,var,varid) )
      call check( nf_get_var(grpid,varid,cloudFraction) )

      var="sensorChannelNumber"
      call check( nf_inq_varid(grpid,var,varid) )
      call check( nf_get_var(grpid,varid,sensorChannelNumber) )

      var="spatialSampleQualityLw"
      call check( nf_inq_varid(grpid,var,varid) )
      call check( nf_get_var(grpid,varid,spatialSampleQualityLw) )
 
      var="spatialSampleQualityMw"
      call check( nf_inq_varid(grpid,var,varid) )
      call check( nf_get_var(grpid,varid,spatialSampleQualityMw) )
 
      var="sensorCentralWavenumber"
      call check( nf_inq_varid(grpid,var,varid) )
      call check( nf_get_var(grpid,varid,sensorCentralWavenumber) )

      group = "ObsValue"; call check(NF_INQ_NCID(ncid, group, grpid))

      # read the spectral radiance in SI meters units
      var="spectralRadiance"
      call check( nf_inq_varid(grpid,var,varid) )
      call check( nf_get_var(grpid,varid,spectralRadiance) )

      # scale the radiance to CRTM units
      do chan=1,channels
      spectralRadiance(chan,:) = spectralRadiance(chan,:) * 1.e5   
      enddo

      end subroutine 
!---------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------
 subroutine check(status)

 use netcdf

 integer, intent(in) :: status

 if(status /= NF90_NOERR) then
    write (6,*) status
    write (6,*) NF90_STRERROR(status)
    call tracebackqq()
    call stop2(99)
 end if

 end subroutine check
!---------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------
 subroutine getdim(ncid,dimname,dim)
 use netcdf
 character(*)  dimname
 integer       ncid,dimid,dim  
 call check( nf90_inq_dimid(ncid,dimname,dimid))         
 call check( nf90_inquire_dimension(ncid,dimid,len=dim))
 end subroutine getdim
!---------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------
