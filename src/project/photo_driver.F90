program photo_rate_demo

  use micm_photo_kinetics,             only : photo_kinetics_t
  use micm_environment,                only : environment_t
  use musica_config,                   only : config_t
  use musica_assert,                   only : die_msg
  use micm_photolysis_wavelength_grid, only : wavelength_grid
  use musica_constants,                only : musica_ik, musica_dk
  use micm_photolysis_wavelength_grid, only : wavelength_grid_initialize, wavelength_grid

  !> Kinetics
  class(photo_kinetics_t), pointer :: photo_kinetics
  !> Environmental state
  type(environment_t) :: environment
  !> Path to the configuration file
  character(len=256) :: config_file_name
  !> Configuration data
  type(config_t) :: config

  real(musica_dk), parameter :: hc = 6.626068e-34_musica_dk * 2.99792458e8_musica_dk
  character(len=*), parameter :: Iam = 'photo_rate_demo: '
  character(len=*), parameter :: wavelength_grid_filespec = '/photo-demo/wavelength_grid.nc'
  integer(musica_ik), parameter :: noErr = 0
  integer :: lambda, ndx
  real(musica_dk)              :: jaccum
  real(musica_dk), allocatable :: jvals_(:), dlambda(:), wrk(:)

  ! Get the model configuration file from the command line
  ! This is a test to see how this comes out in docker
  if( command_argument_count( ) .ne. 1 ) then
    call die_msg( 100000001, 'Usage: ./photo_rate_demo configuration_file.json' )
  end if
  call get_command_argument( 1, config_file_name )
  call config%from_file( config_file_name )

  write(*,*) Iam,'Opened and read config file ',trim(config_file_name)

  ! do all the initialization stuff
    if( wavelength_grid_initialize(wavelength_grid_filespec) /= noErr ) then
      call die_msg( 100000001,Iam//'failed to initialize the photolysis wavelength grid' )
    endif

  ! initialize the kinetics_t object
  photo_kinetics => photo_kinetics_t( config )

  ! do time looping

    ! get new environmental conditions from the host model
    environment%temperature = 274.5_musica_dk
    environment%pressure    = 100428.4_musica_dk
    environment%number_density_air = 1.e20_musica_dk

    ! calculate photo rate cross section,quantum yield
    call photo_kinetics%update_for_new_environmental_state( environment )
    write(*,*) ' '
    do ndx = 1,size(photo_kinetics%cross_section_objs_)
      write(*,*) ' '
      write(*,*) Iam,'cross section values for photorate ',ndx
      do lambda = 1,size(photo_kinetics%cross_section_values_,dim=1)
        if( photo_kinetics%cross_section_values_(lambda,ndx) /= 0.0 ) then
          write(*,*) lambda,wavelength_grid%wcenter(lambda),photo_kinetics%cross_section_values_(lambda,ndx)
        endif
      enddo
    enddo
    do ndx = 1,size(photo_kinetics%quantum_yield_objs_)
      write(*,*) ' '
      write(*,*) Iam,'quantum yield values for photorate ',ndx
      do lambda = 1,size(photo_kinetics%quantum_yield_values_,dim=1)
        if( photo_kinetics%quantum_yield_values_(lambda,ndx) /= 0.0 ) then
          write(*,*) lambda,wavelength_grid%wcenter(lambda),photo_kinetics%quantum_yield_values_(lambda,ndx)
        endif
      enddo
    enddo

  ! end time loop

  ! finalize MICM


    associate( nSize => wavelength_grid%nwave )
      allocate(dlambda(nSize),wrk(nSize))
      dlambda(:) = wavelength_grid%wedge(2:nSize+1) - wavelength_grid%wedge(1:nSize)
    end associate
    wrk(:) = dlambda(:)*wavelength_grid%etf(:)*1.e-13_musica_dk*wavelength_grid%wcenter(:)/hc
    allocate(jvals_(0))
    do ndx = 1,size(photo_kinetics%cross_section_objs_)
      jaccum = sum( photo_kinetics%cross_section_values_(:,ndx) &
                    *photo_kinetics%quantum_yield_values_(:,ndx) &
                    *wrk(:) )
      jvals_ = [jvals_,jaccum]
    enddo

    write(*,*) Iam,'Photorate constants at the top of the atmosphere (1/s)'
    write(*,*) jvals_(:)

    write(*,*) Iam,'Done'

end program photo_rate_demo
