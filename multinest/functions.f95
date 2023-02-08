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
  real(kind=prec) :: btilde(size(z)), kfac
  integer, parameter :: N = 20
  integer k
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

  smooth_exp = smooth_exp * 2 / dm

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


  PURE FUNCTION LVC_INT(X)
  ! This returns
  !  /\ y                         /\ y                1
  !  |    dx  Slvc(mmin+x dm) =   |    dx -------------------------
  ! \/  0                        \/  0     1 - exp((1-2x)/(x-x^2))
  !
  real(kind=prec), intent(in) :: x(:)
  real(kind=prec) :: lvc_int(size(x))

  real(kind=prec) :: y(size(x))
  integer, parameter :: n = 31
  real(kind=prec) :: d(size(x), n)
  integer j
  real(kind=prec), parameter :: tab(n) = (/ &
    0.16790870845090712667_prec,  &
    0.24999999999999999714_prec,  &
    0.090025770283039018322_prec,  &
    0._prec,  &
    -0.0086900366867500027511_prec,  &
    0._prec,  &
    0.00060381235790050824993_prec,  &
    0._prec,  &
    0.00022998828036085831279_prec,  &
    0._prec,  &
    -0.000078972992026240712111_prec,  &
    0._prec,  &
    -5.3118908386568023622e-6_prec,  &
    0._prec,  &
    6.0803996989299871997e-6_prec,  &
    0._prec,  &
    5.0847562954713294469e-7_prec,  &
    0._prec,  &
    -4.4499957697981702844e-7_prec,  &
    0._prec,  &
    -1.4537309204451550916e-7_prec,  &
    0._prec,  &
    1.175946802213709722e-8_prec,  &
    0._prec,  &
    2.8726825067562003393e-8_prec,  &
    0._prec,  &
    8.3072230051403210069e-9_prec,  &
    0._prec,  &
    -2.9478257498685722721e-9_prec,  &
    0._prec,  &
    -2.2149737143110270241e-9_prec /)

  do j=1,size(x)
    d(j, :) = tab
  enddo
  y = -1+2*x

  do j = n, 3, -1
    d(:, j-1) = d(:, j-1) + 2*y(:)*d(:, j)
    d(:, j-2) = d(:, j-2) - d(:, j)
  end do
  lvc_int = d(:, 1)+y(:)*d(:, 2)

  END FUNCTION LVC_INT

  PURE FUNCTION ERF_INT(X)
  ! This returns
  !  /\ y                         /\ y     1  +-              -+
  !  |    dx  Serf(mmin+x dm) =   |    dx --- |  1 - Erf(2-4x) |
  ! \/  0                        \/  0     2  +-              -+
  !
  real(kind=prec), intent(in) :: x(:)
  real(kind=prec) :: erf_int(size(x))
  real(kind=prec), parameter :: e4 = exp(-4.)
  real(kind=prec), parameter :: erf2 = 2*erf(2.)
  real(kind=prec), parameter :: sqpi = 1/sqrt(pi)


  erf_int = (-e4 + exp(-4 * (1-2*x)**2 )) * sqpi &
          + 4 * x - erf2 + (2-4*x) * erf(2-4*x)
  erf_int = erf_int / 8

  END FUNCTION ERF_INT


  PURE FUNCTION SMOOTH_EXPint(M, MI, DM)
  implicit none
  real(kind=prec), intent(in) :: m(:)
  real(kind=prec), intent(in) :: mi, dm
  real(kind=prec) :: smooth_expint(size(m))
  smooth_expint = dm * lvc_int((m - mi)/dm)
  END FUNCTION SMOOTH_EXPINT


  PURE FUNCTION SMOOTH_ERFint(M, MI, DM)
  implicit none
  real(kind=prec), intent(in) :: m(:)
  real(kind=prec), intent(in) :: mi, dm
  real(kind=prec) :: smooth_erfint(size(m))
  smooth_erfint = dm * erf_int((m - mi)/dm)
  END FUNCTION SMOOTH_ERFINT


                          !!!!!!!!!!!!!!!!!!!!!!
                            END MODULE functions
                          !!!!!!!!!!!!!!!!!!!!!!
