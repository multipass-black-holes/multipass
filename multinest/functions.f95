                          !!!!!!!!!!!!!!!!!!!!!!
                              MODULE functions
                          !!!!!!!!!!!!!!!!!!!!!!

  implicit none
  integer, parameter :: prec = selected_real_kind(15,32)
  real (kind=prec), parameter :: pi = 3.14159265358979323846_prec

contains


  ! Implements the normalised power-law with cut-offs
  PURE FUNCTION POWERLAW(M, SPEC, MI, MA)
  implicit none
  real(kind=prec), intent(in) :: m(:)
  real(kind=prec), intent(in) :: spec
  real(kind=prec), intent(in), optional :: mi, ma
  real(kind=prec) :: powerlaw(size(m))

  powerlaw = m**spec * (1+spec)
  if(present(mi) .and. present(ma)) then
    powerlaw = powerlaw / (ma**(1+spec) - mi**(1+spec))
  endif

  END FUNCTION POWERLAW


  PURE FUNCTION GAUSS(M, MU, SIG)
  implicit none
  real(kind=prec), intent(in) :: m(:)
  real(kind=prec), intent(in) :: mu, sig
  real(kind=prec) :: gauss(size(m))

  gauss = exp( -0.5*(m-mu)*(m-mu) / sig / sig ) / sig / sqrt(2*pi)

  END FUNCTION GAUSS


  ! This implements the smooth function (B6)
  PURE FUNCTION SMOOTH(M, MI, DM)
  implicit none
  real(kind=prec), intent(in) :: m(:)
  real(kind=prec), intent(in) :: mi, dm
  real(kind=prec) :: smooth(size(m))

  smooth = 0.

  where ((mi < m) .and. (m < mi + dm)) &
    smooth = exp( dm / (m-mi) + dm / (m - mi - dm)) + 1

  where (m > mi + dm) &
    smooth = 1.

  END FUNCTION SMOOTH


  ! This implements the alternative smooth function using tanh
  PURE FUNCTION SMOOTH_TANH(M, MI, DM)
  implicit none
  real(kind=prec), intent(in) :: m(:)
  real(kind=prec), intent(in) :: mi, dm
  real(kind=prec) :: smooth_tanh(size(m))
  real(kind=prec) :: stargx2(size(m))

  smooth_tanh = 1.

  stargx2 = 2 * (m-mi) / dm - 1

  where ((mi < m) .and. (m < mi + dm)) &
    smooth_tanh = 0.5 * (1+tanh(2*stargx2 / (1-stargx2*stargx2)))

  END FUNCTION SMOOTH_TANH


                          !!!!!!!!!!!!!!!!!!!!!!
                            END MODULE functions
                          !!!!!!!!!!!!!!!!!!!!!!
