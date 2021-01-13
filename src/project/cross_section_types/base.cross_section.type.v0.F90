! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> This base_cross_section module

!> The base_cross_section type and related functions
module micm_base_cross_section_type

  use micm_abs_cross_section_type,     only : abs_cross_section_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: base_cross_section_t

  !> Calculator for base_cross_section
  type, extends(abs_cross_section_t) :: base_cross_section_t
    !> The cross section array
    real(musica_dk), allocatable :: cross_section(:)
  contains
    !> Calculate the cross section
    procedure :: calculate
    !> clean up
    final     :: finalize
  end type base_cross_section_t

  !> Constructor of base cross section objects
  interface base_cross_section_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of base_cross_section_t objects
  function constructor( config ) result( base_cross_section_obj )

    use musica_config,                   only : config_t
    use musica_string,                   only : string_t
    use micm_photolysis_wavelength_grid, only : wavelength_grid
    use netcdf_util,                     only : read_netcdf_file, data_t
    use photo_utils,                     only : addpnt, inter2
    use musica_assert,                   only : die_msg

    !> New base cross section calculator
    class(abs_cross_section_t), pointer :: base_cross_section_obj
    !> cross section configuration data
    type(config_t), intent(inout) :: config

!   local variables
    integer(musica_ik), parameter :: iONE = 1_musica_ik
    real(musica_dk), parameter :: rZERO = 0.0_musica_dk
    real(musica_dk), parameter :: rONE  = 1.0_musica_dk
    real(musica_dk), parameter :: deltax = 1.e-5_musica_dk
    character(len=*), parameter :: Iam = 'base cross_section constructor: '
    character(len=*), parameter :: Hdr = 'XS_'

    type(string_t) :: netcdfFileSpec
    type(string_t) :: addpnt_type_
    type(data_t), allocatable :: dataTray(:)

    integer(musica_ik) :: retcode
    integer :: ndx, m
    integer :: nRows,nCols
    real(musica_dk) :: lowerVal, upperVal
    real(musica_dk) :: addpnt_val_
    real(musica_dk), allocatable :: data_lambda(:)
    real(musica_dk), allocatable :: data_xsect(:)
    real(musica_dk), allocatable :: mdl_lambda_edge(:)
    logical :: found

    write(*,*) Iam,'entering'
    allocate( base_cross_section_t :: base_cross_section_obj )
    !> get cross section netcdf filespec
    call config%get( 'netcdf file', netcdfFileSpec, Iam )
    !> read netcdf file cross section data
    call read_netcdf_file( filespec=netcdfFileSpec%to_char(), Hdr=Hdr, dataTray=dataTray )

    if( size(DataTray) < iONE ) then
      call die_msg( 400000001, Iam//'something wrong with read_netcdf_file for '//trim(netcdfFileSpec%to_char()) )
    endif
    nRows = size(dataTray(1)%data,dim=1)
    nCols = size(dataTray(1)%data,dim=2)
    if( nCols < 2 ) then
      call die_msg( 400000002, Iam//'dataTray(1) data array has < 2 columns' )
    endif
    if( dataTray(1)%scaling_factor /= rONE ) then
      dataTray(1)%data(:,2) = dataTray(1)%scaling_factor*dataTray(1)%data(:,2)
    endif

    data_lambda = dataTray(1)%data(:,1)
    lowerVal = data_lambda(1) ; upperVal = data_lambda(nRows)
    data_xsect = dataTray(1)%data(:,2)

    call base_cross_section_obj%addpnts( config, data_lambda, data_xsect )

    !> add endpoints to data arrays; first the lower bound
    call config%get( 'lower addpnt type', addpnt_type_, Iam, found=found )
    if( .not. found ) then
      call config%get( 'lower addpnt value', addpnt_val_, Iam, found=found )
      if( .not. found ) then
        addpnt_val_ = rZERO
      endif
    else
        addpnt_val_ = data_xsect(1)
    endif

    call addpnt(x=data_lambda,y=data_xsect,xnew=rZERO,ynew=addpnt_val_) 
    call addpnt(x=data_lambda,y=data_xsect,xnew=(rONE-deltax)*lowerVal,ynew=addpnt_val_) 
    !> add endpoints to data arrays; now the upper bound
    call config%get( 'upper addpnt type', addpnt_type_, Iam, found=found )
    if( .not. found ) then
      call config%get( 'upper addpnt value', addpnt_val_, Iam, found=found )
      if( .not. found ) then
        addpnt_val_ = rZERO
      endif
    else
        addpnt_val_ = data_xsect(nRows)
    endif

    call addpnt(x=data_lambda,y=data_xsect,xnew=(rONE+deltax)*upperVal,ynew=addpnt_val_) 
    call addpnt(x=data_lambda,y=data_xsect,xnew=1.e38_musica_dk,ynew=addpnt_val_) 

    mdl_lambda_edge = wavelength_grid%wedge
    select type ( base_cross_section_obj )
      class is (base_cross_section_t) 
        if( .not. allocated(base_cross_section_obj%cross_section) ) then
          allocate(base_cross_section_obj%cross_section(wavelength_grid%nwave))
        endif
        call inter2(xto=mdl_lambda_edge,yto=base_cross_section_obj%cross_section,xfrom=data_lambda,yfrom=data_xsect,ierr=retcode)
      class default
        call die_msg( 400000003, Iam//"Mismatched cross section type" )
    end select

    write(*,*) Iam,'exiting'

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function calculate( this, environment ) result( cross_section )

    use micm_environment,                only : environment_t
    use micm_photolysis_wavelength_grid, only : wavelength_grid

    !> Calculated cross section
    real(kind=musica_dk)              :: cross_section(wavelength_grid%nwave)
    !> base cross section
    class(base_cross_section_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in) :: environment

    character(len=*), parameter :: Iam = 'base cross section calculate: '

    write(*,*) Iam,'entering'

    cross_section = this%cross_section

    write(*,*) Iam,'exiting'

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> finalize the cross section type
   subroutine finalize( this )

   type(base_cross_section_t), intent(inout) :: this

   character(len=*), parameter :: Iam = 'base cross section finalize: '

   write(*,*) Iam,'entering'
   if( allocated(this%cross_section ) ) then
     deallocate(this%cross_section )
   endif
   write(*,*) Iam,'exiting'
   
   end subroutine finalize

end module micm_base_cross_section_type
