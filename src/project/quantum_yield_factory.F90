! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The micm_quantum_yield_factory module

!> Builder of quantum yield calculators
module micm_quantum_yield_factory

  use micm_abs_quantum_yield_type,              only : abs_quantum_yield_t
  use micm_base_quantum_yield_type,             only : base_quantum_yield_t
  use micm_ch3cho_ch3_hco_quantum_yield_type,   only : ch3cho_ch3_hco_quantum_yield_t

  implicit none

  private
  public :: quantum_yield_builder

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Builder of quantum yield calculators
  function quantum_yield_builder( config ) result( new_quantum_yield_t )

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    !> New rate constant calculator
    class(abs_quantum_yield_t), pointer :: new_quantum_yield_t
    !> cross section configuration data
    type(config_t), intent(inout) :: config

    type(string_t) :: quantum_yield_type
    character(len=*), parameter :: Iam = 'quantum yield builder: '

    write(*,*) Iam,'entering'
    new_quantum_yield_t => null()
    call config%get( 'quantum yield type', quantum_yield_type, Iam )

    select case( quantum_yield_type%to_char() )
      case( 'base quantum yield' )
        allocate( base_quantum_yield_t :: new_quantum_yield_t )
      case( 'CH3CHO+hv->CH3+HCO_qy_t' )
        allocate( ch3cho_ch3_hco_quantum_yield_t :: new_quantum_yield_t )
      case default
        call die_msg( 450768214, "Invalid quantum yield type: '"//              &
                                 quantum_yield_type%to_char( )//"'" )
    end select
    call new_quantum_yield_t%initialize( config )
    write(*,*) Iam,'exiting'

  end function quantum_yield_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_quantum_yield_factory
