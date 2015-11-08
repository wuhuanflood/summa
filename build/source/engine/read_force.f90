! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module read_force_module
USE netcdf
USE netcdf_util_module,only:check
USE get_ixname_module,only:get_ixForce
implicit none
private
public::read_force
contains


 ! ************************************************************************************************
 ! public subroutine read_force: read in forcing data
 ! ************************************************************************************************
 subroutine read_force(step_in_Forcfile,ncid, err, message)
 USE nrtype                                             ! variable types, etc.
 USE data_struc,only:gru_struc, nGRU, nHRU              ! GRU data structures
 !USE summaFileManager,only:INPUT_PATH                  ! path of the forcing data file
 !USE time_utils_module,only:extractTime,compJulday     ! extract time info from units string
 !USE multiconst,only:secprday                          ! number of seconds in a day
 !USE data_struc,only:forcFileInfo                      ! forcing file info
 !USE data_struc,only:data_step                         ! length of the data step (s)
 !USE data_struc,only:dJulianStart                      ! julian day of start time of simulation
 !USE data_struc,only:refTime,refJulday                 ! reference time
 !USE data_struc,only:fracJulDay                        ! fractional julian days since the start of year
 !USE data_struc,only:yearLength                        ! number of days in the current year
 USE data_struc,only:forc_meta                          ! metadata structures
 !USE data_struc,only:time_data,time_hru                ! time information
 USE data_struc,only:forc_gru                           ! forcing data
 !USE var_lookup,only:iLookTIME,iLookFORCE              ! named variables to define structure elements
 implicit none
 ! define dummy variables
 integer(i4b),intent(out)          :: ncid             ! the netCDF file id
 integer(i4b),intent(in)           :: step_in_Forcfile ! the step (or time offset) within the current opened forcing netCDF file
 integer(i4b),intent(out)          :: err              ! error code
 character(*),intent(out)          :: message          ! error message
 ! define local variables
 integer(i4b),parameter            :: missingInteger= -9999     ! missing integer
 real(dp),parameter                :: amiss= -1.d+30   ! missing real
 integer(i4b)                      :: varIndx          ! index of variable in its data structure
 integer(i4b)                      :: nVar_forc        ! number of variables in forc_meta struct
 character(len=256)                :: infile           ! filename
 integer(i4b)                      :: iForc            ! index of variable in forc_meta struct
 character(len=256)                :: cmessage         ! error message for downwind routine
 integer(i4b)                      :: iHRU             ! index of HRU
 integer(i4b)                      :: iGRU             ! index of HRU
 integer(i4b)                      :: hruCount         ! total number of HRUs in a GRU
 integer(i4b)                      :: hru_ix           ! index of HRU in the vector for entire domain
 ! real(dp)                          :: dsec             ! double precision seconds (not used)
 ! real(dp)                          :: juldayFirst      ! julian day of the first time step in the data file
 ! real(dp)                          :: startJulDay      ! julian day at the start of the year
 ! real(dp)                          :: currentJulday    ! Julian day of current time step
 ! logical(lgt),parameter            :: checkTime=.false.  ! flag to check the time

 ! define variables for netCDF operations
 integer(i4b), allocatable         :: start(:), count(:) ! specifications for netCDF reading
 real(dp),allocatable              :: buf_flt_1d(:)      ! temporal buffer for reading data
 integer(i4b)                      :: mode               ! netCDF file open mode
 integer(i4b)                      :: varid, grp_ncid    ! variable ids for netCDF file 


 ! Start procedure here
 err=0; message="read_force/"


 ! **********************************************************************************************
 ! ***** read forcing data to data structure, i.e., forc_hru 
 ! **********************************************************************************************

 ! **********************************************************************************************
 ! (0) get the number of variables in forc_meta or forc_hru  structure
 ! **********************************************************************************************
 ! check that metadata structures are initialized
 if(.not.associated(forc_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 nVar_forc = size(forc_meta)
 ! **********************************************************************************************
 ! (1) open HRU_ATTRIBUTES file.
 ! **********************************************************************************************

 ! get hru_id dimension length
 !call check(nf90_inq_dimid(ncid, "hru_id", hruDimID),message)
 !call check(nf90_inquire_dimension(ncid, hruDimID, len = nHRU),message)

 ! allocate space
 !call alloc_type(nHRU,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 !call alloc_attr(nHRU,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 allocate(buf_flt_1d(nHRU), stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for buffer reading'; return; endif

 allocate ( start(2), count(2), stat=err )
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for netcdf reading specification'; return; endif
 ! **********************************************************************************************
 ! (2) read local attributes (types)  from netCDF file group "local_types"
 ! (ie., 'hruIndex','vegTypeIndex','soilTypeIndex','slopeTypeIndex','downHRUindex')
 ! **********************************************************************************************
 varIndx= missingInteger
 
 ! get group id
 call check(nf90_inq_ncid(ncid, 'forcings_input', grp_ncid), message) 
 
 ! check if all desired forcing variables are available in the netCDF file 
 ! it is optional now and it should be in read_meta.f90
 do iForc=2,nVar_forc
  call check(nf90_inq_varid(grp_ncid,  trim(forc_meta(iForc)%varname),varid), message)
 end do

 ! read type variables and put into type_hru structures for all HRUs
 start = (/ 1, step_in_Forcfile/)
 count = (/ nHRU, 1 /)
 
 do iForc=2,nVar_forc
   
  varIndx = get_ixForce(forc_meta(iForc)%varname)
  ! check that the variable could be identified in the data structure
  if(varIndx < 1)then; err=20; message=trim(message)//'unable to find variable ['//trim(forc_meta(iForc)%varname)//'] in data structure'; return; endif
   
  call check(nf90_inq_varid(grp_ncid, forc_meta(iForc)%varname,varid), message)
  call check(nf90_get_var(grp_ncid, varid, buf_flt_1d, start = start, count = count), message)
 
  do iGRU=1, nGRU  
   hruCount = gru_struc(iGRU)%hruCount
   do iHRU=1, hruCount
    hru_ix = gru_struc(iGRU)%hru(iHRU)%hru_ix
    forc_gru(iGRU)%hru(iHRU)%var(varIndx)=buf_flt_1d(hru_ix)
    ! print*, forc_meta(iForc)%varname, forc_hru(iHRU)%var(varIndx)
   end do
  enddo 
 
 end do ! end of "do iForc=1,nVar_forc" 
 
 ! close the HRU_ATTRIBUTES netCDF file
 !call check(nf90_close(ncid), message)

 ! **********************************************************************************************
 ! (4) deallocate space
 ! **********************************************************************************************
 deallocate(buf_flt_1d, stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem deallocating space'; return; endif

 !  ! define the reference time for the model simulation
 !  call extractTime(forc_meta(iLookFORCE%time)%varunit,    & ! input  = units string for time data
 !                   refTime%var(iLookTIME%iyyy),           & ! output = year
 !                   refTime%var(iLookTIME%im),             & ! output = month
 !                   refTime%var(iLookTIME%id),             & ! output = day
 !                   refTime%var(iLookTIME%ih),             & ! output = hour
 !                   refTime%var(iLookTIME%imin),dsec,      & ! output = minute/second
 !                   err,cmessage)                            ! output = error code and error message
 !  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 !  ! convert the reference time to days since the beginning of time
 !  call compjulday(refTime%var(iLookTIME%iyyy),            & ! input  = year
 !                  refTime%var(iLookTIME%im),              & ! input  = month
 !                  refTime%var(iLookTIME%id),              & ! input  = day
 !                  refTime%var(iLookTIME%ih),              & ! input  = hour
 !                  refTime%var(iLookTIME%imin),dsec,       & ! input  = minute/second
 !                  refJulday,err,cmessage)                   ! output = julian day (fraction of day) + error control
 !  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 !  ! identify the start index
 !  time_data%var(:) = imiss
 !  ! put data in time structure
 !   read(cline(time_ix(iline)),*,iostat=err) time_data%var(iline)
 !   !print*,trim(time_meta(iline)%varname),time_data%var(iline)
 !  
 !  ! compute the julian date of the first time index
 !  call compjulday(time_data%var(iLookTIME%iyyy),            & ! input  = year
 !                  time_data%var(iLookTIME%im),              & ! input  = month
 !                  time_data%var(iLookTIME%id),              & ! input  = day
 !                  time_data%var(iLookTIME%ih),              & ! input  = hour
 !                  time_data%var(iLookTIME%imin),dsec,       & ! input  = minute/second
 !                  juldayFirst,err,cmessage)                   ! output = julian day (fraction of day) + error control
 ! 
 ! ! compute the julian day at the start of the year
 ! call compjulday(time_data%var(iLookTIME%iyyy),          & ! input  = year
 !                 1, 1, 1, 1, 0._dp,                      & ! input  = month, day, hour, minute, second
 !                 startJulDay,err,cmessage)                 ! output = julian day (fraction of day) + error control
 ! if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! ! compute the fractional julian day for the current time step
 ! call compjulday(time_data%var(iLookTIME%iyyy),           & ! input  = year
 !                 time_data%var(iLookTIME%im),             & ! input  = month
 !                 time_data%var(iLookTIME%id),             & ! input  = day
 !                 time_data%var(iLookTIME%ih),             & ! input  = hour
 !                 time_data%var(iLookTIME%imin),0._dp,     & ! input  = minute/second
 !                 currentJulday,err,cmessage)                ! output = julian day (fraction of day) + error control
 ! if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! ! compute the time since the start of the year (in fractional days)
 ! fracJulday = currentJulday - startJulDay
 ! ! compute time since the reference time (in seconds)
 ! forc_data%var(iLookFORCE%time) = (currentJulday-refJulday)*secprday
 
 ! test
 ! if(checkTime)then
 !  write(*,'(i4,1x,4(i2,1x),f9.3,1x,i4)') time_data%var(iLookTIME%iyyy),           & ! year
 !                                         time_data%var(iLookTIME%im),             & ! month
 !                                         time_data%var(iLookTIME%id),             & ! day
 !                                         time_data%var(iLookTIME%ih),             & ! hour
 !                                         time_data%var(iLookTIME%imin),           & ! minute
 !                                         fracJulday,                              & ! fractional julian day for the current time step
 !                                         yearLength                                 ! number of days in the current year
 !  !pause ' checking time'
 ! endif
 
 end subroutine read_force


end module read_force_module
