! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The micm_cross_section_factory module

!> Builder of cross section calculators
module micm_cross_section_factory

  use micm_abs_cross_section_type,     only : abs_cross_section_t
  use micm_base_cross_section_type,    only : base_cross_section_t

  implicit none

  private
  public :: cross_section_builder

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Builder of cross section calculators
  !!
  !! At minimum, the \c config argument must include a top-level key-value
  !! pair "type" whose value is a valid rate constant type. Currently, these
  !! are:
  !! - "foo"
  !! - "bar"
  !!
  function cross_section_builder( config ) result( new_cross_section_t )

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    !> New rate constant calculator
    class(abs_cross_section_t), pointer :: new_cross_section_t
    !> cross section configuration data
    type(config_t), intent(inout) :: config

    type(string_t) :: cross_section_type
    character(len=*), parameter :: Iam = 'cross section builder: '

    write(*,*) Iam,'entering'
    new_cross_section_t => null( )
    call config%get( 'cross section type', cross_section_type, Iam )

    select case( cross_section_type%to_char() )
      case( 'base cross section' )
        new_cross_section_t => base_cross_section_t( config )
      case default
        call die_msg( 450768214, "Invalid cross section type: '"//              &
                                 cross_section_type%to_char( )//"'" )
    end select
    write(*,*) Iam,'exiting'

  end function cross_section_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_cross_section_factory
