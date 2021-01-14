! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> This hno3->oh+no2 cross_section module

!> The hno3->oh+no2_cross_section type and related functions
module micm_hno3_oh_no2_cross_section_type

  use micm_abs_cross_section_type,     only : abs_cross_section_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: hno3_oh_no2_cross_section_t, cross_section_t

  type cross_section_t
    real(musica_dk), allocatable :: array(:,:)
  end type cross_section_t

  !> Calculator for hno3_oh_no2_cross_section
  type, extends(abs_cross_section_t) :: hno3_oh_no2_cross_section_t
    !> The cross section array
    type(cross_section_t), allocatable :: cross_section(:)
  contains
    !> Initialize the cross section
    procedure :: initialize
    !> Calculate the cross section
    procedure :: calculate
    !> clean up
    final     :: finalize
  end type hno3_oh_no2_cross_section_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize hno3_oh_no2_cross_section_t object
  subroutine initialize( this, config )

    use musica_config,                   only : config_t
    use musica_string,                   only : string_t
    use micm_photolysis_wavelength_grid, only : wavelength_grid
    use netcdf_util,                     only : read_netcdf_file, data_t
    use photo_utils,                     only : inter2
    use musica_assert,                   only : die_msg

    !> cross section configuration object
    type(config_t), intent(inout) :: config
    !> hno3->oh_no2 cross section type
    class(hno3_oh_no2_cross_section_t), intent(inout) :: this

!   local variables
    integer(musica_ik), parameter :: iONE = 1_musica_ik
    real(musica_dk), parameter :: rZERO = 0.0_musica_dk
    real(musica_dk), parameter :: rONE  = 1.0_musica_dk
    character(len=*), parameter :: Iam = 'hno3->oh+no2 cross section initialize: '
    character(len=*), parameter :: Hdr = 'XS_'

    type(string_t) :: netcdfFileSpec
    type(string_t) :: addpntVal
    type(data_t), allocatable :: dataTray(:)
    type(config_t) :: tmp_config

    integer(musica_ik) :: retcode
    integer(musica_ik) :: colNdx, trayNdx, nCols
    real(musica_dk), allocatable :: data_lambda(:)
    real(musica_dk), allocatable :: data_xsect(:)
    real(musica_dk), allocatable :: mdl_lambda_edge(:)
    logical :: found
    character(len=:), allocatable :: msg
    character(len=:), allocatable :: addpntKey

    write(*,*) Iam,'entering'
    !> get cross section netcdf filespec
    call config%get( 'netcdf file', netcdfFileSpec, Iam, found=found )
has_data_file: &
    if( found ) then
    !> read netcdf file cross section data
      call read_netcdf_file( filespec=netcdfFileSpec%to_char(), Hdr=Hdr, dataTray=dataTray )

    !> check dataTray for validity
      if( size(dataTray) < iONE ) then
        call die_msg( 400000001, Iam//'something wrong with read_netcdf_file for '//trim(netcdfFileSpec%to_char()) )
      endif

      allocate( this%cross_section(size(dataTray)) )

tray_loop: &
      do trayNdx = 1,size(dataTray)
        nCols = size(dataTray(trayNdx)%data,dim=2)
        if( nCols < 2 ) then
          write(msg,*) Iam//'dataTray(',trayNdx,') data array has < 2 columns'
          call die_msg( 400000002, msg )
        endif

        if( .not. allocated(this%cross_section(trayNdx)%array) ) then
          allocate(this%cross_section(trayNdx)%array(wavelength_grid%nwave,nCols-1))
        endif

        do colNdx = 2,nCols
          if( colNdx == 2 ) then
            dataTray(trayNdx)%data(:,colNdx) = 1.e-20_musica_dk*dataTray(trayNdx)%data(:,colNdx)
          elseif( colNdx == 3 ) then
            dataTray(trayNdx)%data(:,colNdx) = 1.e-3_musica_dk*dataTray(trayNdx)%data(:,colNdx)
          endif

          data_lambda = dataTray(trayNdx)%data(:,1)
          data_xsect  = dataTray(trayNdx)%data(:,colNdx)

          if( colNdx == 2 ) then
            call this%addpnts( config, data_lambda, data_xsect )
          elseif( colNdx == 3 ) then
            tmp_config = config
            addpntKey = 'lower addpnt type'
            addpntVal = 'boundary'
            call tmp_config%add( addpntKey, addpntVal, Iam )
            addpntKey = 'upper addpnt type'
            addpntVal = 'boundary'
            call tmp_config%add( addpntKey, addpntVal, Iam )
            call this%addpnts( tmp_config, data_lambda, data_xsect )
          endif
          mdl_lambda_edge = wavelength_grid%wedge
          call inter2(xto=mdl_lambda_edge,yto=this%cross_section(trayNdx)%array(:,colNdx-1), &
                      xfrom=data_lambda,yfrom=data_xsect,ierr=retcode)
        enddo
      enddo tray_loop
    endif has_data_file

    write(*,*) Iam,'exiting'

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function calculate( this, environment ) result( cross_section )

    use micm_environment,                only : environment_t
    use micm_photolysis_wavelength_grid, only : wavelength_grid

    !> Calculated cross section
    real(kind=musica_dk)              :: cross_section(wavelength_grid%nwave)
    !> hno3->oh+no2 cross section
    class(hno3_oh_no2_cross_section_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in) :: environment

    character(len=*), parameter :: Iam = 'hno3->oh+no2 cross section calculate: '
    real(musica_dk) :: Temp

    write(*,*) Iam,'entering'

    Temp = environment%temperature - 298._musica_dk 
    cross_section = this%cross_section(1)%array(:,1)*exp( this%cross_section(1)%array(:,2)*Temp )

    write(*,*) Iam,'exiting'

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> finalize the cross section type
   subroutine finalize( this )

   type(hno3_oh_no2_cross_section_t), intent(inout) :: this

   character(len=*), parameter :: Iam = 'hno3->oh+no2 cross section finalize: '
   integer(musica_ik) :: ndx

   write(*,*) Iam,'entering'
   do ndx = 1,size(this%cross_section)
     if( allocated(this%cross_section(ndx)%array ) ) then
       deallocate(this%cross_section(ndx)%array )
     endif
   enddo
   write(*,*) Iam,'exiting'
   
   end subroutine finalize

end module micm_hno3_oh_no2_cross_section_type
