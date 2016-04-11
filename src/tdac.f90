!! ------------------------------------------------------------------------- !!
!>
!! TDAC: Tsunami Data Assimilation Code
!!
!! @author
!!   Takuto Maeda
!!
!! @see
!!   Maeda  et al. (2015) GRL doi:10.1002/2015GL065588
!!   Gusman et al. (2016) submitted
!<
!! --
program tdac

  use m_params
  use m_llw2d
  use m_matrix

  implicit none
  character(256) :: title_da  = "da"   !< output file title
  character(256) :: title_syn = "syn"  !< output file title for synthetic data

  !! control parameters for optimum interpolations: See document
  real(DP), parameter :: RHO = 1.     !< Ratio between obs and bg error
  real(DP), parameter :: RR  = 20000. !< Cutoff distance of error covariance (m)

  !! model sizes
  integer, parameter :: NN = 3*Nx*Ny
  integer, parameter :: NTMAX = 3000  !< #time step

  !! visualization
  integer, parameter :: ntdec = 10 !< decimation factor for visualization

  !! Model vector for data assimilation
  !!   m*(        1:  Nx*Ny): tsunami height eta(nx,ny)
  !!   m*(  Nx*Ny+1:2*Nx*Ny): vertically integrated velocity Mx(nx,ny)
  !!   m*(2*Nx*Ny+1:3*Nx*Ny): vertically integrated velocity Mx(nx,ny)
  real(DP), allocatable :: mf(:)  !< model vector: forecasted
  real(DP), allocatable :: ma(:)  !< model vector: data-assimilated
  real(DP), allocatable :: mt(:)  !< model vector: true wavefield (observation)

  real(DP), allocatable :: ww(:,:)         !< weight matrix
  real(DP), allocatable :: yt(:)           !< observed tsunami height
  real(DP), allocatable :: yf(:)           !< forecasted tsunami height
  real(DP), allocatable :: d_obs(:)        !< Forecast-observation residual
  integer,  allocatable :: ist(:), jst(:)  !< station location in digital grids
  integer  :: it
  !! ----

  !!
  !! memory allocation
  !!
  allocate( mf(NN), ma(NN), mt(NN) )
  allocate( ww(NN,No) )
  allocate( yt(No), yf(No) )
  allocate( d_obs(No) )
  allocate( ist(No), jst(No))

  mf = 0.0
  mt = 0.0
  ma = 0.0

  !!
  !! Set up tsunami model
  !!
  call llw2d__setup() ! obtain initial tsunami height
  call llw2d__initheight(mt(1:nx*ny))
  call llw2d__set_stations(ist, jst)

  !!
  !! Calculate weight matrix used in the data assimilation
  !!
  call set_weight()


  !!
  !! Time-marching of tsunami wavefiled & data assimilation
  !!
  do it=1, NTMAX

    if( mod(it-1, ntdec) == 0 ) then
      write(STDERR,* ) "timestep = ", it
    end if

    !!
    !! save tsunami wavefield snapshot for visualization
    !!   note that m*(1:Nx*Ny) corresponds to the tsunami height
    !!
    ! assimilation
    if ( mod(it-1, ntdec) == 0 ) then
      call llw2d__output_snap(ma(1:nx*ny), (it-1)/ntdec, title_da)
      call llw2d__output_snap(mt(1:nx*ny), (it-1)/ntdec, title_syn)
    endif

    !! Retrieve "observation" data from synthetic true wavefield
    !!   This part is specialized for synthetic test.
    !!   This example code calculates synthetic data parallel to
    !!   the data assimilation and use it as "observed" data.
    !!   This part can be substituted by real observation data.

    call tsunami_update( mt, mt )  !! integrate true synthetic wavefield

    call get_obs( mt, yt )             !! generate observed data

    !!
    !! Forecast-Assimilate
    !!

    !! tsunami forecast
    call tsunami_update( ma, mf ) !< ma-->mf forecasting

    !! Retrieve forecasted tsunami waveform data at stations
    call get_obs( mf, yf )

    !! Residual
    d_obs(:) = yt(:) - yf(:)

    !! Assimilation
    ma = mf + matmul( ww, d_obs )

  end do


contains

  !! ----------------------------------------------------------------------- !!
  !>
  !! Weight matrix based on the Optimum Interpolation method
  !<
  !! --
  subroutine set_weight

    real(DP), allocatable :: mu_bgo(:,:)
    real(DP), allocatable :: mu_boo(:,:)
    real(DP) :: mat(No,No)
    integer  :: ii, jj, i, j, ig
    real(DP) :: dist

    allocate( mu_bgo(nn,no))
    allocate( mu_boo(no,no ))
    mu_bgo(:,:) = 0.0
    mu_boo(:,:) = 0.0

    !! Estimate background error between numerical grid and station
    do ig=1, Nx*Ny
      do i=1, No

        !!
        !! Background error depend on spatial distance
        !!

        !! Gaussian correlation function
        ii = mod( ig, Nx )
        jj = (ig - ii)/Nx + 1
        call get_distance(ist(i), jst(i), ii, jj, dist)
        mu_bgo(ig,i) = exp( - (dist/RR)**2 )
      end do
    end do


    !! Estimate background error between stations
    do j=1, No
      do i = 1, No
        !! Gaussian correlation function
        call get_distance( ist(i), jst(i), ist(j), jst(j), dist )
        mu_boo(i,j) = exp( - (dist/RR)**2 )
      end do
    end do


    !! Calculate inverse matrix for obtaining weight matrix
    do j=1, No
      do i=1, No
        mat(i,j) = mu_boo(i,j)
        if(i==j) mat(i,j) = mat(i,j) + rho**2
      end do
    end do

    !! invert weight vector
    do ig=1, 3*nx*ny
      call matrix__gs(no, mat, mu_bgo(ig,:), ww(ig,:))
    end do

    deallocate( mu_bgo, mu_boo )

  end subroutine set_weight
  !! ----------------------------------------------------------------------- !!

  !! ----------------------------------------------------------------------- !!
  !>
  !! grid-to-grid distance
  !<
  !! --
  subroutine get_distance( i0, j0, i1, j1, dist )

    integer, intent(in) :: i0, j0, i1, j1
    real(DP), intent(out) :: dist

    dist = sqrt( ( (i0-i1)*dx )**2 + ( (j0-j1)*dy )**2  )

  end subroutine get_distance
  !! ----------------------------------------------------------------------- !!

  !! ----------------------------------------------------------------------- !!
  !>
  !! Return waveform data at stations from given model parameter matrix
  !<
  !! --
  subroutine get_obs( mm, obs )
    real(DP), intent(in) :: mm(NN)
    real(DP), intent(out) :: obs(No)
    integer :: i, ii, jj, iptr

    do i=1,No
      ii = ist(i)
      jj = jst(i)
      iptr = (jj-1) * Nx + ii
      obs(i) = mm( iptr )
    end do

  end subroutine get_obs
  !! ----------------------------------------------------------------------- !!

  !! ----------------------------------------------------------------------- !!
  !>
  !! Integrate linear long-wave equation for updating tsunami model
  !<
  !! --
  subroutine tsunami_update( xa, xf )

    real(DP), target, intent(inout) :: xa(nn)  !< input: assimilated
    real(DP), target, intent(inout) :: xf(nn)  !< output: forecasted

    real(DP), pointer :: eta_a(:)
    real(DP), pointer :: mm_a(:)
    real(DP), pointer :: nn_a(:)
    real(DP), pointer :: eta_f(:)
    real(DP), pointer :: mm_f(:)
    real(DP), pointer :: nn_f(:)
    integer :: nn = nx*ny

    eta_a => xa(1:nn)
    mm_a  => xa(nn+1:2*nn)
    nn_a  => xa(2*nn+1:3*nn)
    eta_f => xf(1:nn)
    mm_f  => xf(nn+1:2*nn)
    nn_f  => xf(2*nn+1:3*nn)

    !! Parts of model vector are aliased to tsunami heiht and velocities
    call llw2d__timestep( eta_a, mm_a, nn_a, eta_f, mm_f, nn_f )

  end subroutine tsunami_update
  !! ----------------------------------------------------------------------- !!

end program tdac
!! ------------------------------------------------------------------------- !!
