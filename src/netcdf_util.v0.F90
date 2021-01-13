   module netcdf_util

   use nc4fortran, only : netcdf_file
   use musica_constants, only : musica_ik, musica_dk

   implicit none

   private
   public :: data_t, read_netcdf_file

   type data_t
     real(musica_dk) :: scaling_factor
     real(musica_dk), allocatable :: data(:,:)
     character(len=:), allocatable :: dataType
   end type data_t

   contains

   subroutine read_netcdf_file( filespec, dataTray )

   character(len=*), intent(in) :: filespec
   type(data_t), allocatable, intent(out) :: dataTray(:)

   integer(musica_ik), parameter :: noErr = 0
   integer(musica_ik), parameter :: nDataVarMax = 10     ! maximum data variables per netcdf file
   character(len=*), parameter :: Iam = "read_netcdf_file: "

   integer(musica_ik) :: nHdrs
   integer(musica_ik) :: stat
   integer(musica_ik) :: hdrNdx, ndx, dataNdx
   integer(musica_ik) :: nRows, nCols
   integer(musica_ik), allocatable :: nxsqyVars(:)
   integer(musica_ik), allocatable :: nDataRows(:), nDataCols(:)
   integer(musica_ik), allocatable :: dims(:)
   character(:), allocatable :: Hdr(:)
   character(:), allocatable :: varName, attrName
   type(data_t) :: wrk_data_t
   type(netcdf_file) :: ncObj

!-----------------------------------------------------
!  open the netcdf file
!-----------------------------------------------------
   call ncObj%initialize(filespec, ierr=stat, status='old', action='r')
   if( stat /= noErr ) then
     write(*,*) Iam,'retcode from initialize = ',stat
     stop 'FileOpenError'
   endif

!-----------------------------------------------------
!  for now hardwire header strings
!-----------------------------------------------------
   Hdr   = (/'XS_','QY_'/)
   allocate(nxsqyVars(0))
!-----------------------------------------------------
!  the dimensions
!-----------------------------------------------------
!  first get count of xs, qy data arrays
!  (no more than 10 allowed)
!-----------------------------------------------------
   varName = '            '
   nHdrs = size(Hdr)
   do hdrNdx = 1,nHdrs
     do ndx = 1,nDataVarMax
       write(varName,'(a,i1)') Hdr(hdrNdx),ndx-1
       if( .not. ncObj%exist(varName) ) then
         exit
       endif
     enddo
     nxsqyVars = [nxsqyVars,ndx - 1]
   enddo

   allocate(nDataRows(0),nDataCols(0))
!-----------------------------------------------------
!  get dimension values
!-----------------------------------------------------
   do hdrNdx = 1,nHdrs
     do ndx = 1,nxsqyVars(hdrNdx)
       call ncObj%shape( varName, dims )
       nDataRows = [nDataRows,dims(1)]
       nDataCols = [nDataCols,dims(2)]
     end do
   end do

   allocate(dataTray(0))
   varName = '                    '
   dataNdx = 0
!-----------------------------------------------------
!  get input data, transfer to "dataTray"
!-----------------------------------------------------
   do hdrNdx = 1,nHdrs
     do ndx = 1,nxsqyVars(hdrNdx)
       write(varName,'(a,i1)') Hdr(hdrNdx),ndx-1
       dataNdx = dataNdx + 1
       nRows = nDataRows(dataNdx) ; nCols = nDataCols(dataNdx)
       allocate(wrk_data_t%data(nRows,nCols))
       call ncObj%read( varName, wrk_data_t%data )

       if( Hdr(hdrNdx)(:2) == 'XS' ) then
         wrk_data_t%dataType = 'cross-section'
       elseif( Hdr(hdrNdx)(:2) == 'QY' ) then
         wrk_data_t%dataType = 'quantum-yield'
       endif
       attrName = 'scaling_factor'
       call ncObj%read_attribute( VarName,attrName,wrk_data_t%scaling_factor,stat)
       if( stat == noErr ) then
         write(*,*) Iam,': ',varName,' ',attrName,' = ',wrk_data_t%scaling_factor
       else
         write(*,*) Iam,': ',varName,' ',attrName,' does not exist'
       endif

       dataTray = [dataTray,wrk_data_t]
       if( allocated( wrk_data_t%data ) ) then
         deallocate(wrk_data_t%data)
       endif
     end do
   end do

!-----------------------------------------------------
!  close the netcdf file
!-----------------------------------------------------
   call ncObj%finalize()

   end subroutine read_netcdf_file

   end module netcdf_util
