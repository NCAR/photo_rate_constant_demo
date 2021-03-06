
  program test_inter2

  use photo_utils, only : inter2
  use musica_constants, only : musica_dk

  implicit none

  character(len=*), parameter :: Iam = 'test_inter2: '
  integer :: n, stat
  real(musica_dk), allocatable :: xfrom(:)
  real(musica_dk), allocatable :: yfrom(:)
  real(musica_dk), allocatable :: xto(:)
  real(musica_dk), allocatable :: yto(:)

  xfrom = (/ (real(n),n=1,10) /)
  yfrom = (/ (20.+real(n),n=1,10) /)

  xto = (/ (real(n),n=2,9) /)
  yto = (/ (0._musica_dk,n=2,9) /)

  write(*,*) Iam,'size xfrom = ',size(xfrom)
  write(*,*) xfrom
  write(*,*) Iam,'size yfrom = ',size(yfrom)
  write(*,*) yfrom
  write(*,*) ' '
  write(*,*) Iam,'size xto = ',size(xto)
  write(*,*) xto
  write(*,*) Iam,'size yto = ',size(yto)
  write(*,*) yto

  call inter2(xto=xto,yto=yto,xfrom=xfrom,yfrom=yfrom,ierr=stat)
  write(*,*) Iam,'return code = ',stat

  write(*,*) ' '
  write(*,*) Iam,'size yto = ',size(yto)
  write(*,*) yto

  end program test_inter2
