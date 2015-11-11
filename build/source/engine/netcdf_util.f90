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

module netcdf_util_module
USE nrtype
USE netcdf
implicit none
private
public::check
public::file_open
public::next_timestep
contains


 ! *********************************************************************************************************
 ! public subroutine file_open: open file
 ! *********************************************************************************************************
 subroutine file_open(infile,mode,ncid,err,message)
 implicit none
 ! declare dummy variables
 character(*),intent(in)              :: infile      ! filename
 integer(i4b),intent(in)              :: mode        ! file open mode
 integer(i4b),intent(out)             :: ncid        ! file unit
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! declare local variables
 logical(lgt)                         :: xist        ! .TRUE. if the file exists
 logical(lgt)                         :: xopn        ! .TRUE. if the file is already open
 ! initialize errors
 err=0; message="f-file_open/"
 ! check if the file exists
 inquire(file=trim(infile),exist=xist) ! Check for existence of file
 if(.not.xist)then
   message=trim(message)//"FileNotFound[file='"//trim(infile)//"']"
   err=10; return
 endif
 ! check if the file is already open
 inquire(file=trim(infile),opened=xopn) ! Check if the file is open
 if(xopn)then
  message=trim(message)//"FileAlreadyOpen['"//trim(infile)//"']"
  err=20; return
 endif

 ! open file
 err=nf90_open(infile, mode, ncid) 
 if(err/=nf90_noerr) then
   message=trim(message)//"OpenError['"//trim(infile)//"']"//trim(nf90_strerror(err))
   err=20; return
 endif

 end subroutine file_open

! ***********************************************************************************************
! Public subroutine next_timestep:  Given the current time and timestep, to get the next time stamp
! ***********************************************************************************************
subroutine next_timestep(step_currentyear,step_currentmonth, step_currentday, &
   step_currenthour, step_currentmin,step_currentsec, timestep)
 
 implicit none
 integer(i4b), intent(inout) :: step_currentyear, step_currentmonth, step_currentday 
 integer(i4b), intent(inout) :: step_currenthour, step_currentmin
 integer(i4b), intent(inout)    :: step_currentsec
 integer(i4b), intent(inout)    :: timestep  ! in seconds
 integer(i4b)                :: leap
 integer(i4b)                :: month(12)
 month = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
 leap=0

 step_currentsec = step_currentsec + timestep
 if(step_currentsec>=60) then
   step_currentmin = step_currentmin + INT(step_currentsec/60)
   step_currentsec = mod(step_currentsec,60)
 endif 

 if(step_currentmin>=60) then
   step_currenthour= step_currenthour + INT(step_currentmin/60)
   step_currentmin = mod(step_currentmin, 60) 
 endif

 if(step_currenthour>=24) then
   step_currentday = step_currentday + INT(step_currenthour/24)
   step_currenthour = mod(step_currenthour,24)
 endif

 if( mod(step_currentyear,4) ==0  .and.  &
   mod(step_currentyear,100) /=0 .or.    &
   mod(step_currentyear,400) ==0 ) then
   leap =1
 endif

 if (leap ==1) then
   month(2) = 29
 endif
 
 if(step_currentday > month(step_currentmonth) ) then
   step_currentday= 1
   step_currentmonth=step_currentmonth+1
 endif 

 if(step_currentmonth>12) then
   step_currentmonth= 1
   step_currentyear= step_currentyear + 1
 endif

end subroutine next_timestep


! ***********************************************************************************************
! check the status of netCDF file operation and return error message 
! ***********************************************************************************************
 subroutine check(status,message)
   integer, intent ( in)       :: status
   character(*),intent(out)    :: message     ! error message
   
   if(status /= nf90_noerr) then 
     message = trim(message)//'netCDF operation error:'//trim(nf90_strerror(status))
     print *, message
     stop "Stopped"
   end if
 end subroutine check

end module netcdf_util_module
