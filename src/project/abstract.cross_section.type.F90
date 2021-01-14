! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The photolysis cross section module

!> The abstract cross section type and related functions
module micm_abs_cross_section_type

  use musica_constants,                only : musica_dk, musica_ik
  use micm_environment,                only : environment_t

  implicit none
  private

  public :: abs_cross_section_t, abs_cross_section_ptr

  !> Photo rate cross section abstract type
  type, abstract :: abs_cross_section_t
  contains
    !> Calculate the photo rate cross section
    procedure(initial),   deferred :: initialize
    procedure(calculate), deferred :: calculate
    procedure                      :: addpnts
  end type abs_cross_section_t

  !> Pointer type for building sets of photo rate constants
  type :: abs_cross_section_ptr
    class(abs_cross_section_t), pointer :: val_ => null( )
  end type abs_cross_section_ptr

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the rate constant for a given set of conditions
  function calculate( this, environment ) result( cross_section )
    use micm_environment,              only : environment_t
    use musica_constants,              only : musica_dk
    use micm_photolysis_wavelength_grid,  only : wavelength_grid
    import abs_cross_section_t

    !> cross section on model photo grid
    real(kind=musica_dk)              :: cross_section(wavelength_grid%nwave)
    !> Cross section calculator
    class(abs_cross_section_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)  :: environment
 end function calculate

  !> Initialize the base cross section type
  subroutine initial( this, config )
    use musica_config,                    only : config_t
    import abs_cross_section_t

    !> Cross section calculator
    class(abs_cross_section_t), intent(inout) :: this
    !> Environmental conditions
    type(config_t), intent(inout) :: config
 end subroutine initial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

contains

  subroutine addpnts( this, config, data_lambda, data_xsect )
    use musica_config, only : config_t
    use musica_string, only : string_t
    use photo_utils,   only : addpnt

    class(abs_cross_section_t)    :: this
    type(config_t), intent(inout) :: config
    real(musica_dk), allocatable, intent(inout) :: data_lambda(:)
    real(musica_dk), allocatable, intent(inout) :: data_xsect(:)

    real(musica_dk), parameter :: rZERO = 0.0_musica_dk
    real(musica_dk), parameter :: rONE  = 1.0_musica_dk
    real(musica_dk), parameter :: deltax = 1.e-5_musica_dk
    character(len=*), parameter :: Iam = 'cross_section; addpnts: '

    integer(musica_ik) :: nRows
    real(musica_dk) :: lowerLambda, upperLambda
    real(musica_dk) :: lowerVal, upperVal
    real(musica_dk) :: addpnt_val_
    type(string_t)  :: addpnt_type_
    logical         :: found

    write(*,*) Iam,'entering'

    !> add endpoints to data arrays; first the lower bound
    nRows = size(data_lambda)
    lowerLambda = data_lambda(1) ; upperLambda = data_lambda(nRows)
    lowerVal    = data_xsect(1)  ; upperVal    = data_xsect(nRows)
    call config%get( 'lower addpnt type', addpnt_type_, Iam, found=found )
    if( .not. found ) then
      call config%get( 'lower addpnt value', addpnt_val_, Iam, found=found )
      if( .not. found ) then
        addpnt_val_ = rZERO
      endif
    else
      addpnt_val_ = lowerVal
    endif

    call addpnt(x=data_lambda,y=data_xsect,xnew=rZERO,ynew=addpnt_val_) 
    call addpnt(x=data_lambda,y=data_xsect,xnew=(rONE-deltax)*lowerLambda,ynew=addpnt_val_) 
    !> add endpoints to data arrays; now the upper bound
    call config%get( 'upper addpnt type', addpnt_type_, Iam, found=found )
    if( .not. found ) then
      call config%get( 'upper addpnt value', addpnt_val_, Iam, found=found )
      if( .not. found ) then
        addpnt_val_ = rZERO
      endif
    else
      addpnt_val_ = upperVal
    endif

    call addpnt(x=data_lambda,y=data_xsect,xnew=(rONE+deltax)*upperLambda,ynew=addpnt_val_) 
    call addpnt(x=data_lambda,y=data_xsect,xnew=1.e38_musica_dk,ynew=addpnt_val_) 

    write(*,*) Iam,'exiting'

  end subroutine addpnts

end module micm_abs_cross_section_type
