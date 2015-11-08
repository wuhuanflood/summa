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

module read_metad_module
USE nrtype
USE netcdf
implicit none
private
public::read_metad
! define named variables to define the type of data structures used
integer(i4b),parameter :: ix_time =1001
integer(i4b),parameter :: ix_force=1002
integer(i4b),parameter :: ix_attr =1003
integer(i4b),parameter :: ix_type =1004
integer(i4b),parameter :: ix_param=1005
integer(i4b),parameter :: ix_mvar =1006
integer(i4b),parameter :: ix_index=1007
integer(i4b),parameter :: ix_bpar =1008
integer(i4b),parameter :: ix_bvar =1009
contains


 ! ********************************************************************************************************
 ! public subroutine read_metad: populate metadata structures
 ! ********************************************************************************************************
 subroutine read_metad(err,message)
 USE summaFileManager,only:SETNGS_PATH                        ! path for metadata files
 USE summaFileManager,only:META_VAR                           ! name of metadata file
 USE data_struc,only:time_meta,forc_meta,attr_meta,type_meta  ! metadata structures
 USE data_struc,only:mpar_meta,mvar_meta,indx_meta            ! metadata structures
 USE data_struc,only:bpar_meta,bvar_meta                      ! metadata structures
 USE netcdf_util_module,only:file_open                        ! open netCDF file
 USE netcdf_util_module, only:check                           ! check status of netcdf file operation     
 implicit none
 ! declare variables
 integer(i4b),intent(out)             :: err                  ! error code
 character(*),intent(out)             :: message              ! error message
 ! local variables
 character(len=1024)                  :: cmessage             ! error message for downstream routine
 integer(i4b)                         :: mode                 ! netCDF file open mod
 integer(i4b)                         :: ncid                 ! netCDF file handle
 ! initialize errors
 err=0; message="read_metad/"

 ! open file
 mode=nf90_NoWrite
 call file_open(trim(SETNGS_PATH)//trim(META_VAR), mode, ncid, err, cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! populate time structure with metadata
 call v_metadata(ncid,ix_time,time_meta,err,cmessage)
 if(err/=0)then; err=40; message=trim(message)//'time/'//trim(cmessage); return; endif
 ! populate local attributes structure with metadata
 call v_metadata(ncid,ix_attr,attr_meta,err,cmessage)
 if(err/=0)then; err=40; message=trim(message)//'attr/'//trim(cmessage); return; endif
 ! populate local category structure with metadata
 call v_metadata(ncid,ix_type,type_meta,err,cmessage)
 if(err/=0)then; err=40; message=trim(message)//'type/'//trim(cmessage); return; endif
 ! populate forcing structure with metadata
 call v_metadata(ncid,ix_force,forc_meta,err,cmessage)
 if(err/=0)then; err=40; message=trim(message)//'forc/'//trim(cmessage); return; endif
 ! populate local parameter structure with metadata
 call v_metadata(ncid,ix_param,mpar_meta,err,cmessage)
 if(err/=0)then; err=40; message=trim(message)//'param/'//trim(cmessage); return; endif
 ! populate local model variable structure with metadata
 call v_metadata(ncid,ix_mvar,mvar_meta,err,cmessage)
 if(err/=0)then; err=40; message=trim(message)//'mvar/'//trim(cmessage); return; endif
 ! populate local model variable structure with metadata
 call v_metadata(ncid,ix_index,indx_meta,err,cmessage)
 if(err/=0)then; err=40; message=trim(message)//'indx/'//trim(cmessage); return; endif
 ! populate basin parameter structure with metadata
 call v_metadata(ncid,ix_bpar,bpar_meta,err,cmessage)
 if(err/=0)then; err=40; message=trim(message)//'bpar/'//trim(cmessage); return; endif
 ! populate basin model variable structure with metadata
 call v_metadata(ncid,ix_bvar,bvar_meta,err,cmessage)
 if(err/=0)then; err=40; message=trim(message)//'bvar/'//trim(cmessage); return; endif

 ! close the HRU_ATTRIBUTES netCDF file
 call check(nf90_close(ncid), message)
 
 end subroutine read_metad


 ! ********************************************************************************************************
 ! private subroutine v_metadata: read metadata from summa_zLocalMeta.nc file and populate the
 !  appropriate metadata structure
 ! ********************************************************************************************************
 subroutine v_metadata(ncid,ivar_lookup,meta_vec,err,message)
 USE netcdf_util_module, only: check                          ! check status of netcdf file operation     
 USE data_struc,only:var_info                                 ! metadata structure
 USE data_struc,only:time_meta,forc_meta,attr_meta,type_meta  ! metadata structures
 USE data_struc,only:mpar_meta,mvar_meta,indx_meta            ! metadata structures
 USE data_struc,only:bpar_meta,bvar_meta                      ! metadata structures

 ! identify indices of named variables
 USE get_ixname_module,only:get_ixTime
 USE get_ixname_module,only:get_ixAttr
 USE get_ixname_module,only:get_ixType
 USE get_ixname_module,only:get_ixForce
 USE get_ixname_module,only:get_ixParam
 USE get_ixname_module,only:get_ixMvar
 USE get_ixname_module,only:get_ixIndex
 USE get_ixname_module,only:get_ixBpar
 USE get_ixname_module,only:get_ixBvar
 implicit none
 ! define input
 integer(i4b),intent(in)              :: ncid           ! input filename
 integer(i4b),intent(in)              :: ivar_lookup    ! used to identify variable type
 ! define input/output
 type(var_info),intent(inout),pointer :: meta_vec(:)    ! vector of metadata
 ! define output
 integer(i4b),intent(out)             :: err            ! error code
 character(*),intent(out)             :: message        ! error message
 ! define local variables
 character(len=256)                   :: cmessage       ! error message for downwind routine
 character(LEN=80)                    :: groupName      ! single lime of information
 type(var_info)                       :: metaTemp       ! temporary metadata structure
 integer(i4b)                         :: ivar,jvar      ! index of model variable
 integer(i4b)                         :: nvar           ! the number of variables in corresponding data structure
 integer(i4b)                         :: nvar_f         ! the number of variables in a group from netcdf file
 integer(i4b)                         :: grp_ncid       ! group id in a netcdf file
 integer(i4b),allocatable             :: varids(:)      ! vector of variable ids in a group

 ! Start procedure here
 err=0; message="v_metadata/"

 !(0) get group name and the number of variables defined in the corresponding meta data structures
 select case(ivar_lookup)
   case(ix_time);  groupName = "TimeMeta";          nvar=size(time_meta);
   case(ix_attr);  groupName = "hru_AttributeMeta"; nvar=size(attr_meta);
   case(ix_type);  groupName = "hru_TypeMeta";      nvar=size(type_meta);
   case(ix_force); groupName = "hru_ForceMeta";     nvar=size(forc_meta);
   case(ix_param); groupName = "hru_ParamMeta";     nvar=size(mpar_meta);
   case(ix_mvar);  groupName = "hru_ModelVarMeta";  nvar=size(mvar_meta);
   case(ix_index); groupName = "hru_ModelIndexMeta";nvar=size(indx_meta);
   case(ix_bpar);  groupName = "BasinParamMeta";    nvar=size(bpar_meta);
   case(ix_bvar);  groupName = "BasinModelVarMeta"; nvar=size(bvar_meta);
   case default; err=35; message=trim(message)//"caseNotFound"; return
 end select
 
 !(1) get group id according to group name and allocate varids
 call check(nf90_inq_ncid(ncid, groupName, grp_ncid), message)

 allocate(varids(nvar), stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problem allocating space for varids"; return; endif
 
 !(2) get variable ids in the group
 call check(nf90_inq_varids(grp_ncid, nvar_f, varids), message)
 if(nvar/=nvar_f) then; err=20; message=trim(message)//"the number of variable in input file is not correct"; return; endif

 !(3) get atrributes of each variable
 do ivar=1, nvar
   call check(nf90_inquire_variable(grp_ncid, varids(ivar), metaTemp%varname), message) 
   call check(nf90_get_att(grp_ncid, varids(ivar), "long_name", metaTemp%vardesc), message)
   call check(nf90_get_att(grp_ncid, varids(ivar), "units", metaTemp%varunit), message)
   call check(nf90_get_att(grp_ncid, varids(ivar), "v_type", metaTemp%vartype), message)

  ! identify the index of the named variable
  select case(ivar_lookup)
   case(ix_time);  jvar = get_ixTime(metaTemp%varname)
   case(ix_attr);  jvar = get_ixAttr(metaTemp%varname)
   case(ix_type);  jvar = get_ixType(metaTemp%varname)
   case(ix_force); jvar = get_ixForce(metaTemp%varname)
   case(ix_param); jvar = get_ixParam(metaTemp%varname)
   case(ix_mvar);  jvar = get_ixMvar(metaTemp%varname)
   case(ix_index); jvar = get_ixIndex(metaTemp%varname)
   case(ix_bpar);  jvar = get_ixBpar(metaTemp%varname)
   case(ix_bvar);  jvar = get_ixBvar(metaTemp%varname)
   case default; err=35; message=trim(message)//"caseNotFound"; return
  end select

  if(jvar<=0)then; err=40; message=trim(message)//"variableNotFound[var="//trim(metaTemp%varname)//"]"; return; endif
  ! check if index is within range
  if(jvar>size(meta_vec))then; err=50; message=trim(message)//"variableExceedsVectorSize[var="//trim(metaTemp%varname)//"]"; return; endif
  ! put data into the metadata vector
  meta_vec(jvar) = metaTemp
 
 enddo  ! looping through variables in the group
 ! check that all elements are populated
 if(any(meta_vec(:)%varname==''))then
  do jvar=1,size(meta_vec)
   print*,jvar,' -> ',trim(meta_vec(jvar)%varname)
  end do
  err=40; message=trim(message)//"someVariablesNotPopulated"; return
 endif

 deallocate(varids, stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem deallocating space'; return; endif

 end subroutine v_metadata


end module read_metad_module
