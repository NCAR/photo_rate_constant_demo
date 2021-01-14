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

  type quantum_yield_t
    real(musica_dk), allocatable :: array(:,:)
  end type quantum_yield_t

  !> Calculator for base quantum yield
  type, extends(abs_quantum_yield_t) :: base_quantum_yield_t
    type(quantum_yield_t), allocatable :: quantum_yield(:)
  contains
    !> Initialize the quantum yield
    procedure :: initialize
    !> Calculate the quantum yield
    procedure :: calculate
    !> clean up
    final     :: finalize
  end type base_quantum_yield_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize base quantum yield_t object
  subroutine initialize( this, config )

    use musica_config,                   only : config_t
    use musica_string,                   only : string_t
    use micm_photolysis_wavelength_grid, only : wavelength_grid
    use netcdf_util,                     only : read_netcdf_file, data_t
    use photo_utils,                     only : inter2
    use musica_assert,                   only : die_msg

    !> quantum yield configuration data
    type(config_t), intent(inout) :: config
    !> New base quantum yield calculator
    class(base_quantum_yield_t), intent(inout) :: this

!   local variables
    character(len=*), parameter :: Iam = 'base quantum yield constructor: '
    character(len=*), parameter :: Hdr = 'QY_'

    type(string_t) :: netcdfFileSpec
    type(data_t), allocatable :: dataTray(:)

    integer(musica_ik) :: retcode
    integer(musica_ik) :: colNdx, trayNdx, nCols
    real(musica_dk) :: quantum_yield_constant
    real(musica_dk), allocatable :: data_lambda(:)
    real(musica_dk), allocatable :: data_qyld(:)
    real(musica_dk), allocatable :: mdl_lambda_edge(:)
    logical :: found
    character(len=:), allocatable :: msg

    write(*,*) Iam,'entering'

    !> check for quantum yield netcdf filespec
    call config%get( 'netcdf file', netcdfFileSpec, Iam, found=found )
    !> read netcdf file quantum yield data
has_netcdf_file: &
    if( found ) then
      call read_netcdf_file( filespec=netcdfFileSpec%to_char(), Hdr=Hdr, dataTray=dataTray )
      if( size(DataTray) < iONE ) then
        call die_msg( 500000001, Iam//'something wrong with read_netcdf_file for '//trim(netcdfFileSpec%to_char()) )
      endif

      allocate( this%quantum_yield(size(dataTray)) )

tray_loop: &
      do trayNdx = 1,size(dataTray) 
        nCols = size(dataTray(trayNdx)%data,dim=2)
        if( nCols < 2 ) then
          call die_msg( 500000002, Iam//'dataTray(1) data array has < 2 columns' )
        endif

        if( .not. allocated(this%quantum_yield(trayNdx)%array) ) then
          allocate(this%quantum_yield(trayNdx)%array(wavelength_grid%nwave,nCols-1))
        endif
        do colNdx = 2,nCols
          if( dataTray(trayNdx)%scaling_factor /= rONE ) then
            dataTray(trayNdx)%data(:,colNdx) = dataTray(trayNdx)%scaling_factor*dataTray(trayNdx)%data(:,colNdx)
          endif

          data_lambda = dataTray(trayNdx)%data(:,1)
          data_qyld   = dataTray(trayNdx)%data(:,colNdx)

          call this%addpnts( config, data_lambda, data_qyld )
          mdl_lambda_edge = wavelength_grid%wedge
          call inter2(xto=mdl_lambda_edge,yto=this%quantum_yield(trayNdx)%array(:,colNdx-1), &
                      xfrom=data_lambda,yfrom=data_qyld,ierr=retcode)
        enddo
      enddo tray_loop
    else has_netcdf_file
    !> check for quantum yield constant
      call config%get( 'quantum yield constant', quantum_yield_constant, Iam, found=found )
      if( found ) then
        allocate( this%quantum_yield(1) )
        allocate(this%quantum_yield(1)%array(wavelength_grid%nwave,1))
        this%quantum_yield(1)%array(:,1) = quantum_yield_constant
      endif
    endif has_netcdf_file

    write(*,*) Iam,'exiting'
   
  end subroutine initialize

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

    quantum_yield = this%quantum_yield(1)%array(:,1)

    write(*,*) Iam,'exiting'

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> finalize the quantum yield type
   subroutine finalize( this )

   type(base_quantum_yield_t), intent(inout) :: this

   character(len=*), parameter :: Iam = 'base quantum yield finalize: '
   integer(musica_dk) :: ndx

   write(*,*) Iam,'entering'
   do ndx = 1,size(this%quantum_yield)
     if( allocated(this%quantum_yield(ndx)%array) ) then
       deallocate( this%quantum_yield(ndx)%array )
     endif
   enddo
   write(*,*) Iam,'exiting'
   
   end subroutine finalize

end module micm_base_quantum_yield_type
