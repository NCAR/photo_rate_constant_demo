! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The micm_cross_section_factory module

!> Builder of cross section calculators
module micm_cross_section_factory

  use micm_abs_cross_section_type,     only : abs_cross_section_t
  use micm_base_cross_section_type,    only : base_cross_section_t
  use micm_n2o5_no2_no3_cross_section_type, only : n2o5_no2_no3_cross_section_t
  use micm_cl2_cl_cl_cross_section_type, only : cl2_cl_cl_cross_section_t
  use micm_hno3_oh_no2_cross_section_type, only : hno3_oh_no2_cross_section_t

  implicit none

  private
  public :: cross_section_builder

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
        allocate( base_cross_section_t :: new_cross_section_t )
      case( 'N2O5+hv->NO2+NO3 cross section' )
        allocate( n2o5_no2_no3_cross_section_t :: new_cross_section_t )
      case( 'Cl2+hv->Cl+Cl cross section' )
        allocate( cl2_cl_cl_cross_section_t :: new_cross_section_t )
      case( 'HNO3+hv->OH+NO2 cross section' )
        allocate( hno3_oh_no2_cross_section_t :: new_cross_section_t )
      case default
        call die_msg( 450768214, "Invalid cross section type: '"//              &
                                 cross_section_type%to_char( )//"'" )
    end select
    call new_cross_section_t%initialize( config )
    write(*,*) Iam,'exiting'

  end function cross_section_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_cross_section_factory
