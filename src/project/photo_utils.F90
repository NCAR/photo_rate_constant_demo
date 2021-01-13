
  module photo_utils

  use musica_constants, only : musica_ik, musica_dk

  implicit none

  private
  public :: inter2, addpnt

  contains
  
    subroutine inter2(xto,yto,xfrom,yfrom,ierr)
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Map input data given on single, discrete points onto a set of target     =*
!=  bins.                                                                    =*
!=  The original input data are given on single, discrete points of an       =*
!=  arbitrary grid and are being linearly interpolated onto a specified set  =*
!=  of target bins.  In general, this is the case for most of the weighting  =*
!=  functions (action spectra, molecular cross section, and quantum yield    =*
!=  data), which have to be matched onto the specified wavelength intervals. =*
!=  The average value in each target bin is found by averaging the trapezoi- =*
!=  dal area underneath the input data curve (constructed by linearly connec-=*
!=  ting the discrete input values).                                         =*
!=  Some caution should be used near the endpoints of the grids.  If the     =*
!=  input data set does not span the range of the target grid, an error      =*
!=  message is printed and the execution is stopped, as extrapolation of the =*
!=  data is not permitted.                                                   =*
!=  If the input data does not encompass the target grid, use ADDPNT to      =*
!=  expand the input array.                                                  =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NG  - INTEGER, number of bins + 1 in the target grid                  (I)=*
!=  XG  - REAL, target grid (e.g., wavelength grid);  bin i is defined    (I)=*
!=        as [XG(i),XG(i+1)] (i = 1..NG-1)                                   =*
!=  YG  - REAL, y-data re-gridded onto XG, YG(i) specifies the value for  (O)=*
!=        bin i (i = 1..NG-1)                                                =*
!=  N   - INTEGER, number of points in input grid                         (I)=*
!=  X   - REAL, grid on which input data are defined                      (I)=*
!=  Y   - REAL, input y-data                                              (I)=*
!-----------------------------------------------------------------------------*

      integer(musica_ik), intent(out) :: ierr
      real(musica_dk), intent(in)  :: xfrom(:), yfrom(:)
      real(musica_dk), intent(in)  :: xto(:)
      real(musica_dk), intent(out) :: yto(:)

! local:
      integer(musica_ik), parameter :: ONE = 1_musica_ik
      integer(musica_ik), parameter :: TWO = 2_musica_ik
      real(musica_ik), parameter    :: ZERO = 0.0_musica_dk

      character(len=*), parameter :: Iam = 'inter2: '

      integer(musica_ik) :: nto, nfrom
      integer(musica_ik) :: ntom1, nfromm1
      integer(musica_ik) :: i, k, jstart
      real(musica_dk)    :: area, xtol, xtou
      real(musica_dk)    :: darea, slope
      real(musica_dk)    :: a1, a2, b1, b2

      ierr   = 0_musica_ik
      nfrom  = size(xfrom)
      nto    = size(xto)
      nfromm1 = nfrom - ONE
      ntom1   = nto - ONE

!-----------------------------------------------------------------------------*
!  check data grid for monotonicity
!-----------------------------------------------------------------------------*
      if( any( xfrom(2:nfrom) <= xfrom(1:nfrom-1) ) ) then
        write(*,*) Iam,'data grid not monotonically increasing'
        stop 'GridErr'
      endif
!-----------------------------------------------------------------------------*
!  do model grid x values lie completley inside data grid x values?
!-----------------------------------------------------------------------------*
      IF ( (xfrom(1) > xto(1)) .or. (xfrom(nfrom) < xto(nto)) ) THEN
        write(*,*)  Iam,'Data do not span grid; Use ADDPNT to expand data and re-run.'
        stop 'GridErr'
      ENDIF
!-----------------------------------------------------------------------------*
!  find the integral of each grid interval and use this to 
!  calculate the average y value for the interval      
!  xtol and xtou are the lower and upper limits of the to grid interval
!-----------------------------------------------------------------------------*
      jstart = ONE
to_interval_loop: &
      do i = ONE,ntom1
        area = ZERO
        xtol = xto(i)
        xtou = xto(i+1)
!-----------------------------------------------------------------------------*
!  discard data before the first grid interval and after the last grid interval
!  for internal grid intervals, start calculating area by interpolating
!  between the last point which lies in the previous interval and the
!  first point inside the current interval
!-----------------------------------------------------------------------------*
        k = jstart
        if (k <= nfromm1) then
!  if both points are before the first grid, go to the next point
          do while( k <= nfromm1 )
            if (xfrom(k+1) <= xtol) then
              jstart = k - 1
              k = k+1
            else
              exit
            endif
          enddo
!-----------------------------------------------------------------------------*
!  if the last point is beyond the end of the grid,
!  then complete and go to the next grid
!-----------------------------------------------------------------------------*
          do while( (k <= nfromm1) .and. (xfrom(k) < xtou) )
            jstart = k-1
! compute x-coordinates of increment
            a1 = MAX(xfrom(k),xtol)
            a2 = MIN(xfrom(k+1),xtou)
!  if points coincide, contribution is zero
            if (xfrom(k+1) == xfrom(k)) then
              darea = ZERO
            else
              slope = (yfrom(k+1) - yfrom(k))/(xfrom(k+1) - xfrom(k))
              b1 = yfrom(k) + slope*(a1 - xfrom(k))
              b2 = yfrom(k) + slope*(a2 - xfrom(k))
              darea = .5_musica_dk*(a2 - a1)*(b2 + b1)
            endif
!  find the area under the trapezoid from a1 to a2
            area = area + darea
! go to next point
            k = k+1
          enddo
        endif
!-----------------------------------------------------------------------------*
!  calculate the average y after summing the areas in the interval
!-----------------------------------------------------------------------------*
        yto(i) = area/(xtou - xtol)
      enddo to_interval_loop

      end subroutine inter2

      SUBROUTINE addpnt ( x, y, xnew, ynew )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Add a point <xnew,ynew> to a set of data pairs <x,y>.  x must be in      =*
!=  ascending order                                                          =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  X    - REAL vector of length LD, x-coordinates                       (IO)=*
!=  Y    - REAL vector of length LD, y-values                            (IO)=*
!=  XNEW - REAL, x-coordinate at which point is to be added               (I)=*
!=  YNEW - REAL, y-value of point to be added                             (I)=*
!-----------------------------------------------------------------------------*

! arguments
      REAL(musica_dk), allocatable, intent(inout) :: x(:), y(:)
      REAL(musica_dk), intent(in)  :: xnew, ynew

! local variables
      character(len=*), parameter :: Iam = 'addpnt: '

      INTEGER(musica_ik) :: n
      INTEGER(musica_ik) :: insertNdx
      logical            :: found
      

      n = size(x)
!-----------------------------------------------------------------------------*
!  check data grid for monotonicity
!-----------------------------------------------------------------------------*
      if( any( x(2:n) <= x(1:n-1) ) ) then
        write(*,*) Iam,'grid not monotonically increasing'
        stop 'GridErr'
      endif
!-----------------------------------------------------------------------------*
!  does xnew == any x value?
!-----------------------------------------------------------------------------*
      if( any( x(:) == xnew ) ) then
        write(*,*) Iam,'xnew exactly matches a grid x value'
        stop 'GridErr'
      endif
!-----------------------------------------------------------------------------*
! find the index at which xnew needs to be inserted into x
!-----------------------------------------------------------------------------*
      found = .true.
      if( xnew < x(1) ) then
        insertNdx = 1
      else if( xnew > x(n) ) then
        insertNdx = n
      else
        found = .false.
        do insertNdx = 2,n
          if (x(insertNdx) > xnew ) then
            found = .true.
            exit
          endif
        enddo
      endif
      if( .not. found ) then
        write(*,*) Iam,'something really wrong; all stop'
        stop 'codeErr'
      endif
!-----------------------------------------------------------------------------*
! increment x,y arrays, then insert xnew,ynew
!-----------------------------------------------------------------------------*
      x = [x,xnew]
      y = [y,ynew]
      if( xnew < x(n) ) then
        x(insertNdx:) = eoshift( x(insertNdx:),shift=-1,boundary=xnew )
        y(insertNdx:) = eoshift( y(insertNdx:),shift=-1,boundary=ynew )
      endif

      end subroutine addpnt

  end module photo_utils
