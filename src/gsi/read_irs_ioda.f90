!---------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------
      module read_irs_data

      real(8),allocatable,dimension(:)  :: dateTime
      real(8),allocatable,dimension(:)  :: dwellNumber
      real(8),allocatable,dimension(:)  :: strokeDirection
      real(8),allocatable,dimension(:)  :: latitude
      real(8),allocatable,dimension(:)  :: longitude
      real(8),allocatable,dimension(:)  :: sensorAzimuthAngle
      real(8),allocatable,dimension(:)  :: sensorZenithAngle
      real(8),allocatable,dimension(:)  :: solarAzimuthAngle
      real(8),allocatable,dimension(:)  :: solarZenithAngle
      real(8),allocatable,dimension(:)  :: cloudSignal
      real(8),allocatable,dimension(:)  :: cloudCoverTotal
      real(8),allocatable,dimension(:,:):: mwGlobalPcScores
      real(8),allocatable,dimension(:)  :: mwGlobalPcrScores
      real(8),allocatable,dimension(:)  :: mwGlobalPcrsQuality
      real(8),allocatable,dimension(:)  :: mwSpatialSampleQuality
      real(8),allocatable,dimension(:)  :: mwResidualEnergy
      real(8),allocatable,dimension(:,:):: lwGlobalPcScores
      real(8),allocatable,dimension(:)  :: lwGlobalPcrScores
      real(8),allocatable,dimension(:)  :: lwGlobalPcrsQuality
      real(8),allocatable,dimension(:)  :: lwSpatialSampleQuality
      real(8),allocatable,dimension(:)  :: lwResidualEnergy
      real(8),allocatable,dimension(:,:):: sensorChannelNumber
      real(8),allocatable,dimension(:,:):: sensorChannelradiance 

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
      call getdim(ncid,"channel",channels)
      print*,"nobs=",nobs

! allocate the data arrays

      allocate(dateTime(nobs))
      allocate(dwellNumber(nobs))
      allocate(strokeDirection(nobs))
      allocate(latitude(nobs))
      allocate(longitude(nobs))
      allocate(sensorAzimuthAngle(nobs))
      allocate(sensorZenithAngle(nobs))
      allocate(solarAzimuthAngle(nobs))
      allocate(solarZenithAngle(nobs))
      allocate(cloudSignal(nobs))
      allocate(cloudCoverTotal(nobs))
      allocate(mwSpatialSampleQuality(nobs))
      allocate(mwResidualEnergy(nobs))
      allocate(lwSpatialSampleQuality(nobs))
      allocate(lwResidualEnergy(nobs))
      allocate(sensorChannelNumber(channels,nobs))
      allocate(sensorChannelradiance(channels,nobs))

! read the data

      group = "MetaData"; call check(NF_INQ_NCID(ncid, group, grpid))

      var="dateTime"
      call check( nf_inq_varid(grpid,var,varid) )
      call check( nf_get_var(grpid,varid,dateTime) )

      var="dwellNumber"
      call check( nf_inq_varid(grpid,var,varid) )
      call check( nf_get_var(grpid,varid,dwellNumber) )

      var="strokeDirection"
      call check( nf_inq_varid(grpid,var,varid) )
      call check( nf_get_var(grpid,varid,strokeDirection) )

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

      var="cloudCoverTotal"
      call check( nf_inq_varid(grpid,var,varid) )
      call check( nf_get_var(grpid,varid,cloudCoverTotal) )

      var="sensorChannelNumber"
      call check( nf_inq_varid(grpid,var,varid) )
      call check( nf_get_var(grpid,varid,sensorChannelNumber) )

      group = "ObsValues"; call check(NF_INQ_NCID(ncid, group, grpid))

      var="mwSpatialSampleQuality"
      call check( nf_inq_varid(grpid,var,varid) )
      call check( nf_get_var(grpid,varid,mwSpatialSampleQuality) )

      var="mwResidualEnergy"
      call check( nf_inq_varid(grpid,var,varid) )
      call check( nf_get_var(grpid,varid,mwResidualEnergy) )

      var="lwSpatialSampleQuality"
      call check( nf_inq_varid(grpid,var,varid) )
      call check( nf_get_var(grpid,varid,lwSpatialSampleQuality) )

      var="lwResidualEnergy"
      call check( nf_inq_varid(grpid,var,varid) )
      call check( nf_get_var(grpid,varid,lwResidualEnergy) )

      var="radiance"
      call check( nf_inq_varid(grpid,var,varid) )
      call check( nf_get_var(grpid,varid,sensorChannelradiance) )

      end subroutine 
!---------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------
 subroutine check(status)

 use netcdf

 integer, intent(in) :: status

 if(status /= NF90_NOERR) then
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
