
  program test_add_pnt

  use photo_utils, only : addpnt
  use musica_constants, only : musica_dk

  implicit none

  character(len=*), parameter :: Iam = 'test_add_pnt: '
  integer :: n
  real(musica_dk), allocatable :: x(:)
  real(musica_dk), allocatable :: y(:)

  x = (/ (real(n),n=1,10) /)
  y = (/ (20.+real(n),n=1,10) /)

  write(*,*) Iam,'original array'
  write(*,*) Iam,'size x = ',size(x)
  write(*,*) x
  write(*,*) Iam,'size y = ',size(y)
  write(*,*) y

  call addpnt( x, y, 9.5_musica_dk, -10._musica_dk )

  write(*,*) ' '
  write(*,*) Iam,'modified array'
  write(*,*) Iam,'size x = ',size(x)
  write(*,*) x
  write(*,*) Iam,'size y = ',size(y)
  write(*,*) y

  end program test_add_pnt
