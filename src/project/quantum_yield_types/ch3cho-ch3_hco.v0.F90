! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> This ch3cho+hv->ch3+hco quantum yield module

!> The ch3cho+hv->ch3+hco quantum yield type and related functions
module micm_ch3cho_ch3_hco_quantum_yield_type

  use micm_abs_quantum_yield_type,     only : abs_quantum_yield_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: ch3cho_ch3_hco_quantum_yield_t

  integer(musica_ik), parameter :: iONE = 1_musica_ik
  real(musica_dk), parameter :: rZERO = 0.0_musica_dk
  real(musica_dk), parameter :: rONE  = 1.0_musica_dk

  !> Calculator for ch3cho+hv->ch3+hco quantum yield
  type, extends(abs_quantum_yield_t) :: ch3cho_ch3_hco_quantum_yield_t
    real(musica_dk), allocatable :: quantum_yield(:,:)
  contains
    !> Calculate the quantum yield
    procedure :: calculate
    !> clean up
    final     :: finalize
  end type ch3cho_ch3_hco_quantum_yield_t

  !> Constructor of ch3cho+hv->ch3+hco quantum yield objects
  interface ch3cho_ch3_hco_quantum_yield_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of ch3cho+hv->ch3+hco quantum yield objects
  function constructor( config ) result( ch3cho_ch3_hco_quantum_yield_obj )

    use musica_config,                   only : config_t
    use musica_string,                   only : string_t
    use micm_photolysis_wavelength_grid, only : wavelength_grid
    use netcdf_util,                     only : read_netcdf_file, data_t
    use photo_utils,                     only : addpnt, inter2
    use musica_assert,                   only : die_msg

    !> New ch3cho+hv->ch3+hco quantum yield calculator
    class(abs_quantum_yield_t), pointer :: ch3cho_ch3_hco_quantum_yield_obj
    !> quantum yield configuration data
    type(config_t), intent(inout) :: config

!   local variables
    real(musica_dk), parameter :: deltax = 1.e-5_musica_dk
    character(len=*), parameter :: Iam = 'ch3cho+hv->ch3_hco quantum yield constructor: '
    character(len=*), parameter :: Hdr = 'QY_'

    type(string_t) :: netcdfFileSpec
    type(string_t) :: addpnt_type_
    type(data_t), allocatable :: dataTray(:)

    integer(musica_ik) :: retcode
    integer :: colNdx, chnl
    integer :: nRows,nCols
    integer :: m
    integer :: nZeroCnt
    real(musica_dk) :: lowerVal, upperVal
    real(musica_dk) :: addpnt_val_
    real(musica_dk), allocatable :: data_lambda(:)
    real(musica_dk), allocatable :: data_qyld(:)
    real(musica_dk), allocatable :: mdl_lambda_edge(:)
    logical :: found

    write(*,*) Iam,'entering'
    allocate( ch3cho_ch3_hco_quantum_yield_t :: ch3cho_ch3_hco_quantum_yield_obj )

    select type ( ch3cho_ch3_hco_quantum_yield_obj )
      class is ( ch3cho_ch3_hco_quantum_yield_t )
    !> get quantum yield netcdf filespec
        call config%get( 'netcdf file', netcdfFileSpec, Iam, found=found )
    !> read netcdf file quantum yield data
        if( found ) then
          call read_netcdf_file( filespec=netcdfFileSpec%to_char(), Hdr=Hdr, dataTray=dataTray )
        else
          call die_msg( 400000004, Iam//'expects to read a netCDF data file' )
        endif

        if( size(DataTray) < iONE ) then
          call die_msg( 400000001, Iam//'something wrong with read_netcdf_file for '//trim(netcdfFileSpec%to_char()) )
        endif
        nRows = size(dataTray(1)%data,dim=1)
        nCols = size(dataTray(1)%data,dim=2)
        if( nCols < 3 ) then
          call die_msg( 400000002, Iam//'dataTray(1) data array has < 3 columns' )
        endif

        if( .not. allocated(ch3cho_ch3_hco_quantum_yield_obj%quantum_yield) ) then
          allocate(ch3cho_ch3_hco_quantum_yield_obj%quantum_yield(wavelength_grid%nwave,2))
        endif

        do colNdx = 2,3
          data_lambda = dataTray(1)%data(:,1)
          lowerVal = data_lambda(1) ; upperVal = data_lambda(nRows)
          data_qyld = dataTray(1)%data(:,colNdx)

    !> add endpoints to data arrays; first the lower bound
          call config%get( 'lower addpnt type', addpnt_type_, Iam, found=found )
          if( .not. found ) then
            call config%get( 'lower addpnt value', addpnt_val_, Iam, found=found )
            if( .not. found ) then
              addpnt_val_ = rZERO
            endif
          else
              addpnt_val_ = data_qyld(1)
          endif

          call addpnt(x=data_lambda,y=data_qyld,xnew=rZERO,ynew=addpnt_val_) 
          call addpnt(x=data_lambda,y=data_qyld,xnew=(rONE-deltax)*lowerVal,ynew=addpnt_val_) 
    !> add endpoints to data arrays; now the upper bound
          call config%get( 'upper addpnt type', addpnt_type_, Iam, found=found )
          if( .not. found ) then
            call config%get( 'upper addpnt value', addpnt_val_, Iam, found=found )
            if( .not. found ) then
              addpnt_val_ = rZERO
            endif
          else
              addpnt_val_ = data_qyld(nRows)
          endif

          call addpnt(x=data_lambda,y=data_qyld,xnew=(rONE+deltax)*upperVal,ynew=addpnt_val_) 
          call addpnt(x=data_lambda,y=data_qyld,xnew=1.e38_musica_dk,ynew=addpnt_val_) 

          mdl_lambda_edge = wavelength_grid%wedge
          if( colNdx == 2 ) then
            chnl = 2
          else
            chnl = 1
          endif
          call inter2(xto=mdl_lambda_edge,yto=ch3cho_ch3_hco_quantum_yield_obj%quantum_yield(:,chnl), &
                      xfrom=data_lambda,yfrom=data_qyld,ierr=retcode)
        enddo
      class default
        call die_msg( 400000003, Iam//"Mismatched quantum yield type" )
    end select

    select type ( ch3cho_ch3_hco_quantum_yield_obj )
      class is ( ch3cho_ch3_hco_quantum_yield_t )
        do m = 1,size(ch3cho_ch3_hco_quantum_yield_obj%quantum_yield,dim=1)
          if( ch3cho_ch3_hco_quantum_yield_obj%quantum_yield(m,1) /= rZERO ) then
            write(*,*) wavelength_grid%wcenter(m), ch3cho_ch3_hco_quantum_yield_obj%quantum_yield(m,1)
          endif
        enddo
        write(*,*) ' '
        do m = 1,size(ch3cho_ch3_hco_quantum_yield_obj%quantum_yield,dim=1)
          if( ch3cho_ch3_hco_quantum_yield_obj%quantum_yield(m,2) /= rZERO ) then
            write(*,*) wavelength_grid%wcenter(m), ch3cho_ch3_hco_quantum_yield_obj%quantum_yield(m,2)
          endif
        enddo
        nZeroCnt = count(ch3cho_ch3_hco_quantum_yield_obj%quantum_yield(:,1) /= rZERO)
        write(*,*) Iam,'there are ',nZeroCnt,' non-zero elements in col #1'
        nZeroCnt = count(ch3cho_ch3_hco_quantum_yield_obj%quantum_yield(:,2) /= rZERO)
        write(*,*) Iam,'there are ',nZeroCnt,' non-zero elements in col #2'
    end select

    write(*,*) Iam,'exiting'
   
  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function calculate( this, environment ) result( quantum_yield )

    use micm_environment,                only : environment_t
    use micm_photolysis_wavelength_grid, only : wavelength_grid

    !> Calculated cross section
    real(kind=musica_dk)              :: quantum_yield(wavelength_grid%nwave)
    !> base cross section
    class(ch3cho_ch3_hco_quantum_yield_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in) :: environment

    character(len=*), parameter :: Iam = 'ch3cho+hv->ch3_hco quantum yield calculate: '
    integer(musica_ik) :: m
    real(musica_dk) :: air_den_factor
    real(musica_dk), allocatable :: quantum_yield_chnl1(:)
    real(musica_dk), allocatable :: quantum_yield_chnl2(:)
    real(musica_dk), allocatable :: quantum_yield_wrk(:)

    write(*,*) Iam,'entering'

    quantum_yield_chnl1 = this%quantum_yield(:,1)
    quantum_yield_chnl2 = rONE - this%quantum_yield(:,2)
    quantum_yield_wrk = (/ (rZERO,m=1,size(this%quantum_yield,dim=1)) /)
    where( quantum_yield_chnl1(:) > rZERO )
      quantum_yield_wrk = quantum_yield_chnl2/quantum_yield_chnl1 - rONE
    elsewhere
      quantum_yield_wrk = rZERO
    endwhere
    air_den_factor = environment%number_density_air/2.46e19_musica_dk
    quantum_yield = quantum_yield_chnl1*(rONE + quantum_yield_wrk) &
                    /(rONE + quantum_yield_wrk*air_den_factor)
    quantum_yield = min( rONE,max(rZERO,quantum_yield) )

    write(*,*) Iam,'exiting'

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> finalize the quantum yield type
   subroutine finalize( this )

   type(ch3cho_ch3_hco_quantum_yield_t), intent(inout) :: this

   character(len=*), parameter :: Iam = 'ch3cho+hv->ch3_hco quantum yield finalize: '

   write(*,*) Iam,'entering'
   if( allocated(this%quantum_yield ) ) then
     deallocate(this%quantum_yield )
   endif
   write(*,*) Iam,'exiting'
   
   end subroutine finalize

end module micm_ch3cho_ch3_hco_quantum_yield_type
