!! -------------------------------------------------------------------------- !!
!>
!! Linear Long Wave (LLW) tsunami in 2D Cartesian Coordinate
!<
!! --
module m_llw2d

  use m_params

  implicit none
  private

  public :: llw2d__setup
  public :: llw2d__timestep
  public :: llw2d__initheight
  public :: llw2d__set_stations
  public :: llw2d__output_snap

  integer,  parameter :: NXA = 20               !< absorber thickness
  integer,  parameter :: NYA = 20               !< absorber thickness
  real(DP), parameter :: APARA = 0.015_DP       !< damping for boundaries
  real(DP)            :: gg (Nx,Ny)             !< absorbing boundary
  real(DP), parameter :: CUTOFF_DEPTH = 10.0_DP  !< shallowest water depth

  integer, allocatable :: ist(:), jst(:)  !! station locations

  real(DP), allocatable :: hh(:,:)  !< ocean depth
  real(DP), allocatable :: hm(:,:)  !< x-averaged depth
  real(DP), allocatable :: hn(:,:)  !< y-averaged depth
  real(DP), allocatable :: fn(:,:), fm(:,:), fe(:,:)  !< land filters

  integer, parameter :: IO_SNP = 31   !< fortran I/O number

contains


  !! ----------------------------------------------------------------------- !!
  !>
  !! Set station locations. Users may need to modify it
  !<
  !! --
  subroutine llw2d__set_stations( ist_ret, jst_ret )

    integer, intent(out) :: ist_ret(No)
    integer, intent(out) :: jst_ret(No)

    integer :: i, j, nst

    !! synthetic station locations
    nst = 0
    do i=1, 6
      do j=1, 6
        nst = nst + 1
        ist(nst) = floor( ( (i-1) * 20 + 150 ) * 1000 / dx + 0.5 )
        jst(nst) = floor( ( (j-1) * 20 + 150 ) * 1000 / dy + 0.5 )
      end do
    end do

    !! output for visualization
    open(10, file='./out/stloc.dat')
    do i=1, nst
      write(10,*) (ist(i)-1)*dx/1000., (jst(i)-1)*dy/1000
    end do
    close(10)

    !! return copy of station location
    ist_ret(1:No) = ist(1:No)
    jst_ret(1:No) = jst(1:No)

  end subroutine llw2d__set_stations
  !! ----------------------------------------------------------------------- !!


  !! ----------------------------------------------------------------------- !!
  !>
  !! export snapshot data for visualization
  !<
  !! --
  subroutine llw2d__output_snap( eta, isnap, title )

    real(DP), intent(in) :: eta(Nx,Ny)
    integer,  intent(in) :: isnap
    character(*), intent(in) :: title
    integer :: i, j
    character(256) :: fn_out
    character(6) :: csnap

    write(csnap,'(I6.6)') isnap
    fn_out = './out/' // trim(title)//'__'//csnap//'__.dat'

    open(IO_SNP, file=trim(fn_out)  )
    do j=1, ny
          do i=1, nx
        write(io_snp,*) real((i-1)*dx/1000), real((j-1)*dy/1000), real(eta(i,j))
      end do
      write(io_snp,*)
    end do
    close(IO_SNP)

  end subroutine llw2d__output_snap
  !! ----------------------------------------------------------------------- !!


  !! ----------------------------------------------------------------------- !!
  !>
  !! Numerically integrate linear long-wave equation with one time step.
  !! From input of tsunami height (eta0) and velocities (mm0, nn0),
  !! Returns updated height (eta1) and velocities(mm1, nn1)
  !<
  !! --
  subroutine llw2d__timestep( eta0, mm0, nn0,  eta1, mm1, nn1 )

    real(DP), intent(in),  dimension(Nx,Ny) :: eta0, mm0, nn0
    real(DP), intent(out), dimension(Nx,Ny) :: eta1, mm1, nn1

    integer :: i, j
    real(DP) :: dxeta(Nx,Ny), dyeta(Nx,Ny)
    real(DP) :: dxM  (Nx,Ny), dyN  (Nx,Ny)

    ! diffs
     do j=1, Ny
        do i=2, Nx
           dxeta(i,j) = ( eta0(i,j) - eta0(i-1,j)) / dx
        end do
           dxeta(1,j) = ( eta0(1,j) -        0.0 ) / dx
     end do
     do i=1, Nx
        do j=2, Ny
           dyeta(i,j) = ( eta0(i,j) - eta0(i,j-1)) / dy
        end do
           dyeta(i,1) = ( eta0(i,1) -        0.0 ) / dy
     end do

    ! Update Velocity
    do j=1, Ny
       do i=1, Nx
          mm1(i,j) = mm0(i,j) - g0*hm(i,j)*dxeta(i,j)*dt
          nn1(i,j) = nn0(i,j) - g0*hn(i,j)*dyeta(i,j)*dt
       end do
    end do

    !! boundary condition
    do j=1, Ny
       do i=1, Nx
          mm1(i,j) = mm1(i,j) * fm(i,j) * gg(i,j)
          nn1(i,j) = nn1(i,j) * fn(i,j) * gg(i,j)
       end do
    end do

    ! diffs
    do j=1, Ny
       dxM(Nx,j) = ( 0.0 - mm1(Nx,j)  ) / dx
       do i=1, Nx-1
          dxM(i,j) = ( mm1(i+1,j) - mm1(i,j) ) / dx
       end do
    end do
    do i=1, Nx
       dyN(i,Ny) = ( 0.0 - nn1(i,Ny) ) / dy
       do j=1, Ny-1
          dyN(i,j) = ( nn1(i,j+1) - nn1(i,j) ) / dy
       end do
    end do

    ! Update Wave Heigt
    do j=1, Ny
       do i=1, Nx
          eta1(i,j) = eta0(i,j) - ( dxM(i,j) +  dyN(i,j) )*dt
       end do
    end do

    !! boundary condition
    do j=1, Ny
       do i=1, Nx
          eta1(i,j) = eta1(i,j) * fe(i,j) * gg(i,j)
       end do
    end do


  end subroutine llw2d__timestep
  !! ----------------------------------------------------------------------- !!


  !! ----------------------------------------------------------------------- !!
  subroutine llw2d__setup(  )

    integer :: i, j

    !!
    !! Memory allocation
    !!
    allocate(hh(nx,ny), hm(nx,ny), hn(nx,ny))
    allocate(fn(nx,ny), fm(nx,ny), fe(nx,ny))
    allocate(ist(no), jst(no))
    !!
    !! Bathymetry set-up. Users may need to modify it
    !!
    hh(:,:) =  3000.0_DP

    do j=1, ny
      do i=1, nx
        if( hh(i,j) < 0.0_DP ) then
          hh(i,j) = 0.0_DP
        else if ( hh(i,j) < CUTOFF_DEPTH ) then
          hh(i,j) = CUTOFF_DEPTH
        end if
      end do
    end do

    !!
    !! average bathymetry for staggered-grid computation
    !!
    do j=1, Ny
       do i=2, Nx
          hm(i,j) = ( hh(i,j)+hh(i-1,j) ) / 2
          if( hh(i,j) <= 0.0_DP .or. hh(i-1,j) <= 0.0_DP ) hm(i,j) = 0.0_DP
       end do
       hm(1,j) = hh(1,j)
    end do
    do i=1, Nx
       do j=2,Ny
          hn(i,j) = ( hh(i,j) + hh(i,j-1) ) / 2
          if( hh(i,j) <= 0.0_DP .or. hh(i,j-1) <= 0.0_DP ) hn(i,j) = 0.0_DP
       end do
       hn(i,1) = hh(i,1)
    end do


    !!
    !! Land filter
    !!
    fm(:,:) = 1.0_DP
    fn(:,:) = 1.0_DP
    fe(:,:) = 1.0_DP
    do j=1, ny
       do i=1, nx
          if( hm(i,j) < 0 ) fm(i,j) = 0.0_DP
          if( hn(i,j) < 0 ) fn(i,j) = 0.0_DP
          if( hh(i,j) < 0 ) fe(i,j) = 0.0_DP
       end do
    end do

    !!
    !! Sponge absorbing boundary condition by Cerjan (1985)
    !!
    do j = 1,ny
       do i = 1,nx
          if( i <= nxa ) then
             gg(i,j) = exp( -( ( apara * (nxa-i     ) )**2 ) )
          else if ( i >= ( nx-nxa+1 ) ) then
             gg(i,j) = exp( -( ( apara * (i-nx+nxa-1) )**2 ) )
          else if( j <= nya ) then
             gg(i,j) = exp( -( ( apara * (nya-j     ) )**2 ) )
          else if ( j >= ( ny-nya+1 ) ) then
             gg(i,j) = exp( -( ( apara  * (j-ny+nya-1) )**2 ) )
          else
             gg(i,j) = 1.0
          end if
       end do
    end do

  end subroutine llw2d__setup
  !! ----------------------------------------------------------------------- !!


  !! ----------------------------------------------------------------------- !!
  !>
  !! initial tsunami height of synthetic tsunamis
  !<
  !! --
  subroutine llw2d__initheight( eta )

    real(DP), intent(inout) :: eta(Nx,Ny)

    !! source size
    real(DP), parameter :: aa = 30000.
    real(DP), parameter :: bb = 30000.
    integer :: i, j, i0, j0
    real(DP) :: hx, hy
    real(DP) :: pi = atan(1.0) * 4


    !! bathymetry setting
    eta(:,:) = 0.0

    i0 = nx / 4
    j0 = ny / 4
    do j=1, ny
       if( -bb <= (j-j0)*dy  .and. (j-j0) * dy <= bb ) then
          hy = ( 1 + cos( pi * ( j-j0 ) * dy / bb ) ) / 2.0
       else
          hy = 0.0
       end if

       do i=1, nx

          if( -aa <= (i-i0)*dx  .and. (i-i0) * dx <= aa ) then
             hx = ( 1 + cos( pi * ( i-i0 ) * dx / aa ) ) / 2.0
          else
             hx = 0.0
          end if

          eta(i,j) = hx * hy

       end do

    end do

    !! force zero amplitude on land
    do j=1, Ny
       do i=1, Nx
          if( hh(i,j) < epsilon(1.0) ) eta(i,j) = 0.0_DP
       end do
    end do


  end subroutine llw2d__initheight
  !! ----------------------------------------------------------------------- !!


end module m_llw2d
!! ------------------------------------------------------------------------- !!
