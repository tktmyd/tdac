!! ------------------------------------------------------------------------- !!
module m_matrix

  implicit none
  private

  integer, parameter :: DP = selected_real_kind(13)
  integer, parameter :: STDERR  = 0

  public matrix__gs

contains

  !! ----------------------------------------------------------------------- !!
  !>
  !! Solve the well-conditioned linear system
  !! Sum( m(i,j)*a(j) ) == v(i) for a(:)
  !! by Gauss-Seidel method.
  !<
  !! --
  subroutine matrix__gs( n, m, v, a )

    integer,  intent(in)  :: n       !< model size
    real(DP), intent(in)  :: m(n,n)  !< coefficient matrix
    real(DP), intent(in)  :: v(n)    !< values
    real(DP), intent(out) :: a(n)    !< answers
    real(DP), parameter   :: TOL = 1e-6
    real(DP) :: a0(n)
    integer :: iter
    integer :: i, j
    real(DP) :: m_max, wk
    !! --

    a(:) = 0
    a0(:) = 0
    iter = 0
    m_max = maxval(abs(m))

    do
      do i=1, n
        wk = 0
        do j=1, i-1
          wk = wk + m(i,j) * a(j)
        end do
        do j=i+1, n
          wk = wk + m(i,j) * a(j)
        end do
        a(i) = (v(i) - wk) / m(i,i)
      end do
      if( maxval(abs(a0 - a))/m_max < TOL ) exit
      a0(:) = a(:)
      iter = iter + 1
      if(iter > 10000) then
        write(STDERR,*) "does not converge"
        exit
      end if
    end do

  end subroutine matrix__gs
  !! ----------------------------------------------------------------------- !!

end module m_matrix
!! ------------------------------------------------------------------------- !!
