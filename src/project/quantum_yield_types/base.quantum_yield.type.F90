! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> This base quantum yield module

!> The base quantum yield type and related functions
module micm_base_quantum_yield_type

  use micm_abs_quantum_yield_type,     only : abs_quantum_yield_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: base_quantum_yield_t

  integer(musica_ik), parameter :: iONE = 1_musica_ik
  real(musica_dk), parameter :: rZERO = 0.0_musica_dk
  real(musica_dk), parameter :: rONE  = 1.0_musica_dk

  !> Calculator for base quantum yield
  type, extends(abs_quantum_yield_t) :: base_quantum_yield_t
    real(musica_dk), allocatable :: quantum_yield(:)
  contains
    !> Calculate the quantum yield
    procedure :: calculate
  end type base_quantum_yield_t

  !> Constructor of base quantum yield objects
  interface base_quantum_yield_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of base quantum yield objects
  function constructor( config ) result( base_quantum_yield_obj )

    use musica_config,                   only : config_t
    use musica_string,                   only : string_t
    use micm_photolysis_wavelength_grid, only : wavelength_grid
    use netcdf_util,                     only : read_netcdf_file, data_t
    use photo_utils,                     only : inter2
    use musica_assert,                   only : die_msg

    !> New base quantum yield calculator
    class(abs_quantum_yield_t), pointer :: base_quantum_yield_obj
    !> quantum yield configuration data
    type(config_t), intent(inout) :: config

!   local variables
    character(len=*), parameter :: Iam = 'base quantum yield constructor: '
    character(len=*), parameter :: Hdr = 'QY_'

    type(string_t) :: netcdfFileSpec
    type(data_t), allocatable :: dataTray(:)

    integer(musica_ik) :: retcode
    real(musica_dk) :: quantum_yield_constant
    real(musica_dk), allocatable :: data_lambda(:)
    real(musica_dk), allocatable :: data_qyld(:)
    real(musica_dk), allocatable :: mdl_lambda_edge(:)
    logical :: found

    write(*,*) Iam,'entering'
    allocate( base_quantum_yield_t :: base_quantum_yield_obj )

    select type ( base_quantum_yield_obj )
      class is ( base_quantum_yield_t )
    !> get quantum yield netcdf filespec
        call config%get( 'netcdf file', netcdfFileSpec, Iam, found=found )
    !> read netcdf file quantum yield data
has_netcdf_file: &
        if( found ) then
          call read_netcdf_file( filespec=netcdfFileSpec%to_char(), Hdr=Hdr, dataTray=dataTray )
          if( size(DataTray) < iONE ) then
            call die_msg( 500000001, Iam//'something wrong with read_netcdf_file for '//trim(netcdfFileSpec%to_char()) )
          endif
          if( size(dataTray(1)%data,dim=2) < 2 ) then
            call die_msg( 500000002, Iam//'dataTray(1) data array has < 2 columns' )
          endif
          if( dataTray(1)%scaling_factor /= rONE ) then
            dataTray(1)%data(:,2) = dataTray(1)%scaling_factor*dataTray(1)%data(:,2)
          endif

          if( .not. allocated(base_quantum_yield_obj%quantum_yield) ) then
            allocate(base_quantum_yield_obj%quantum_yield(wavelength_grid%nwave))
          endif

          data_lambda = dataTray(1)%data(:,1)
          data_qyld   = dataTray(1)%data(:,2)
          call base_quantum_yield_obj%addpnts( config, data_lambda, data_qyld )

          mdl_lambda_edge = wavelength_grid%wedge
          call inter2(xto=mdl_lambda_edge,yto=base_quantum_yield_obj%quantum_yield, &
                      xfrom=data_lambda,yfrom=data_qyld,ierr=retcode)
        else has_netcdf_file
          call config%get( 'quantum yield constant', quantum_yield_constant, Iam, found=found )
          if( .not. found ) then
            quantum_yield_constant = rONE
          endif
          if( .not. allocated(base_quantum_yield_obj%quantum_yield) ) then
            allocate(base_quantum_yield_obj%quantum_yield(wavelength_grid%nwave))
          endif
          base_quantum_yield_obj%quantum_yield(:) = quantum_yield_constant
        endif has_netcdf_file
      class default
        call die_msg( 400000003, Iam//"Mismatched quantum yield type" )
    end select

    write(*,*) Iam,'exiting'
   
  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate quantum yield for the environmental conditions
  function calculate( this, environment ) result( quantum_yield )

    use micm_environment,                only : environment_t
    use micm_photolysis_wavelength_grid, only : wavelength_grid

    !> Calculated quantum yield
    real(kind=musica_dk)              :: quantum_yield(wavelength_grid%nwave)
    !> base quantum yield
    class(base_quantum_yield_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in) :: environment

    character(len=*), parameter :: Iam = 'base quantum yield calculate: '

    write(*,*) Iam,'entering'

    quantum_yield = this%quantum_yield

    write(*,*) Iam,'exiting'

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> finalize the quantum yield type
   subroutine finalize( this )

   type(base_quantum_yield_t), intent(inout) :: this

   character(len=*), parameter :: Iam = 'base quantum yield finalize: '

   write(*,*) Iam,'entering'
   if( allocated(this%quantum_yield ) ) then
     deallocate(this%quantum_yield )
   endif
   write(*,*) Iam,'exiting'
   
   end subroutine finalize

end module micm_base_quantum_yield_type
