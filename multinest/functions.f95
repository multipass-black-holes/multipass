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

  PURE FUNCTION BTILDE(A, B, Z)
  real(kind=prec), intent(in) :: a, b
  real(kind=prec), intent(in) :: z(:)
  real(kind=prec) :: btilde(size(z))
  integer, parameter :: N = 10
  integer k, kfac
  real(kind=prec) :: gb, zk(size(z))

  gb = gamma(1-b)

  kfac = 1
  zk = 1
  Btilde = 1._prec / a
  do k=1,n
    kfac = kfac * k
    zk = zk * z
    Btilde = Btilde + gamma(1+k-b) * zk / (a+k) / gb / kfac
  enddo

  Btilde = z**a * Btilde - gamma(a)*gamma(b)/gamma(a+b)
  END FUNCTION BTILDE


  ! This implements the smooth function (B6) of [2010.14533] or (4) of
  ! [2104.02685]
  PURE FUNCTION SMOOTH_EXP(M, MI, DM)
  implicit none
  real(kind=prec), intent(in) :: m(:)
  real(kind=prec), intent(in) :: mi, dm
  real(kind=prec) :: smooth_exp(size(m))

  smooth_exp = 0.

  where ((mi < m) .and. (m < mi + dm)) &
    smooth_exp = 1/(exp( dm / (m-mi) + dm / (m - mi - dm)) + 1)

  where (m > mi + dm) &
    smooth_exp = 1.

  END FUNCTION SMOOTH_EXP


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

  ! This implements (1.17)
  PURE FUNCTION SMOOTH_ERF(m, mi, dm)
  implicit none
  real(kind=prec), intent(in) :: m(:)
  real(kind=prec), intent(in) :: mi, dm
  real(kind=prec) :: smooth_erf(size(m))
  smooth_erf = 0.5 * erf( 4 * (m-mi-dm/2) / dm ) + 0.5
  END FUNCTION SMOOTH_ERF

                          !!!!!!!!!!!!!!!!!!!!!!
                            END MODULE functions
                          !!!!!!!!!!!!!!!!!!!!!!
