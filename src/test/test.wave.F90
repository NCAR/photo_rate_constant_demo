
  program test_wave

  use musica_constants, only : musica_ik
  use micm_photolysis_wavelength_grid, only : wavelength_grid_initialize, wavelength_grid

  implicit none

  character(len=*), parameter :: Iam ='test_wave: '
  integer(musica_ik) :: retcode

  retcode =  wavelength_grid_initialize( '../../wavelength_grid.nc' )
  write(*,*) Iam,'retcode = ',retcode
  write(*,*) Iam,'size of wave center array = ',size(wavelength_grid%wcenter)
  write(*,*) Iam,'size of wave edge   array = ',size(wavelength_grid%wedge)
  write(*,*) Iam,'size of etf         array = ',size(wavelength_grid%etf)
  write(*,*) ' '
  write(*,*) Iam,'wcenter array'
  write(*,*) wavelength_grid%wcenter
  
  end program test_wave
