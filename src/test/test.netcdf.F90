  program test_netcdf

  use musica_constants, only : musica_ik, musica_dk
  use netcdf_util, only : data_t, read_netcdf_file
  use photo_utils, only : addpnt, inter2
  use micm_photolysis_wavelength_grid, only : wavelength_grid_initialize, wavelength_grid

  implicit none

  real(musica_dk), parameter :: ZERO = 0.0_musica_dk
  real(musica_dk), parameter :: ONE  = 1.0_musica_dk
  real(musica_dk), parameter :: deltax = 1.e-5_musica_dk

  integer(musica_ik) :: retcode
  integer :: ndx, m
  integer :: nRows,nCols
  real(musica_dk) :: lowerVal, upperVal
  real(musica_dk), allocatable :: data_lambda(:)
  real(musica_dk), allocatable :: data_xsect(:)
  real(musica_dk), allocatable :: mdl_lambda_edge(:)
  real(musica_dk), allocatable :: mdl_lambda_center(:)
  real(musica_dk), allocatable :: mdl_xsect(:)
  character(len=*), parameter :: Iam = 'test_netcdf: '
  type(data_t), allocatable   :: dataTray(:)

! allocate(dataTray(0))
  call read_netcdf_file( '/Users/stacy/Documents/Python_dev/Sandbox/XSQY/Data/XSQY/CH3CHO_cross_section.nc', 'XS_', dataTray )

  write(*,*) ' '
  write(*,*) Iam,'size of dataTray = ',size(dataTray)
  write(*,*) Iam,'dataTray diagnostics'
  do ndx = 1,size(dataTray)
    nRows = size(dataTray(ndx)%data,dim=1)
    nCols = size(dataTray(ndx)%data,dim=2)
    write(*,*) Iam,'array ',ndx,' is (',nRows,' x ',nCols,')'
    write(*,*) Iam,'array ',ndx,' has scaling factor = ',dataTray(ndx)%scaling_factor
    write(*,*) Iam,'array ',ndx,' is type ',trim(dataTray(ndx)%dataType)
    write(*,*) ' '
    data_lambda = dataTray(ndx)%data(:,1)
    lowerVal = data_lambda(1) ; upperVal = data_lambda(nRows)
    write(*,*) data_lambda
    write(*,*) ' '
    data_xsect = dataTray(ndx)%data(:,2)*dataTray(ndx)%scaling_factor
    write(*,*) Iam,'data_xsect'
    do m = 1,nRows
      if( data_xsect(m) /= ZERO ) then
        write(*,*) data_lambda(m),data_xsect(m)
      endif
    enddo

    write(*,*) ' '
    call addpnt(x=data_lambda,y=data_xsect,xnew=ZERO,ynew=ZERO) 
    call addpnt(x=data_lambda,y=data_xsect,xnew=(ONE-deltax)*lowerVal,ynew=ZERO) 
    call addpnt(x=data_lambda,y=data_xsect,xnew=(ONE+deltax)*upperVal,ynew=ZERO) 
    call addpnt(x=data_lambda,y=data_xsect,xnew=1.e38_musica_dk,ynew=ZERO) 
    write(*,*) ' '
    write(*,*) Iam,'modified data_lambda'
    write(*,*) data_lambda
    write(*,*) ' '
    write(*,*) Iam,'modified data_xsect'
    write(*,*) data_xsect

    retcode = wavelength_grid_initialize( '../../wavelength_grid.nc' )
    if( retcode /= 0_musica_ik ) then
      write(*,*) Iam,'failed to read model wavelength grid'
    endif
    mdl_lambda_edge   = wavelength_grid%wedge
    mdl_lambda_center = wavelength_grid%wcenter
    if( .not. allocated(mdl_xsect) ) then
      allocate(mdl_xsect(wavelength_grid%nwave))
    endif
    mdl_xsect = -99._musica_dk
    call inter2(xto=mdl_lambda_edge,yto=mdl_xsect,xfrom=data_lambda,yfrom=data_xsect,ierr=retcode)
    write(*,*) ' '
    write(*,*) Iam,'interpolated data_xsect'
    do m = 1,wavelength_grid%nwave
      if( mdl_xsect(m) /= ZERO ) then
        write(*,*) mdl_lambda_center(m),mdl_xsect(m)
      endif
    enddo
  enddo

  end program test_netcdf
