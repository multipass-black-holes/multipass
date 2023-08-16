                          !!!!!!!!!!!!!!!!!!!!!!
                              MODULE MODELS
                          !!!!!!!!!!!!!!!!!!!!!!
  use functions
  implicit none
  ! This just controls the overall normalisation
  real(kind=prec), parameter :: mpiv   = 30

  ABSTRACT INTERFACE
    PURE FUNCTION SMOOTHFN(M, MI, DM)
      use functions, only: prec
      real(kind=prec), intent(in) :: m(:), mi, dm
      real(kind=prec) :: smoothfn(size(m))
    END FUNCTION SMOOTHFN
  END INTERFACE

  TYPE PARA
    ! Needed for all
    real(kind=prec) :: mmin=0, dm=0
    ! Needed for PLP
    real(kind=prec) :: mmax=0, mum=0, sm=0, alpha=0, lp=0
    ! Needed for PLP_POW
    real(kind=prec) :: k=0
    ! Needed for PPISN
    real(kind=prec) :: mgap=0, a=0, b=0, d=0

    ! Needed for spin
    real(kind=prec) :: alpha1=0, alpha2=0, beta1=0, beta2=0

    ! Smooth function
    procedure(smoothfn), pointer, nopass :: sf, sfint
    character(len=3) :: sf_c

    real(kind=prec) :: lam21, lam12
  END TYPE PARA


  ABSTRACT INTERFACE
    PURE FUNCTION MASS1FN(M, P)
      use functions, only: prec
      import para
      real(kind=prec), intent(in) :: m(:)
      type(para), intent(in) :: p
      real(kind=prec) :: mass1fn(size(m))
    END FUNCTION MASS1FN

    PURE FUNCTION MASS2FN(M1, M2, P)
      use functions, only: prec
      import para
      real(kind=prec), intent(in) :: m1(:), m2(:)
      type(para), intent(in) :: p
      real(kind=prec) :: mass2fn(size(m1))
    END FUNCTION MASS2FN

    PURE FUNCTION REDSHIFTFN(Z, P)
      use functions, only: prec
      import para
      real(kind=prec), intent(in) :: z(:)
      type(para), intent(in) :: p
      real(kind=prec) :: redshiftfn(size(z))
    END FUNCTION REDSHIFTFN

    PURE FUNCTION SPINFN(chi1, chi2, P)
      use functions, only: prec
      import para
      real(kind=prec), intent(in) :: chi1(:), chi2(:)
      type(para), intent(in) :: p
      real(kind=prec) :: spinfn(size(chi1))
    END FUNCTION SPINFN

    PURE FUNCTION PARAFN(V)
      use functions, only: prec
      import para
      real(kind=prec), intent(in) :: v(:)
      type(para) :: parafn
    END FUNCTION PARAFN
  END INTERFACE

  TYPE SPINFNS
    procedure(spinFn), pointer, nopass :: f
  END TYPE SPINFNS

  TYPE MODEL
    integer ndim
    character(len=100) :: model_name
    procedure(mass1fn), pointer, nopass :: primary
    procedure(mass1fn), pointer, nopass :: primaryM2
    procedure(mass2fn), pointer, nopass :: secondary
    procedure(redshiftFn), pointer, nopass :: redshift
    type(spinFns), dimension(2,2) :: spin
    procedure(paraFn), pointer, nopass :: r2p
    procedure(smoothFn), pointer, nopass :: smooth, smoothInt
    character(len=3) :: smooth_c
    character(len=4) :: secondary_c
    logical :: norms
  END TYPE MODEL

contains

  PURE FUNCTION TRIVIAL(M, P)
  real(kind=prec), intent(in) :: m(:)
  type(para), intent(in) :: p
  real(kind=prec) :: trivial(size(m))
  trivial = 1.
  END FUNCTION TRIVIAL

  PURE FUNCTION TRIVIAL_SPIN(chi1, chi2, P)
  real(kind=prec), intent(in) :: chi1(:), chi2(:)
  type(para), intent(in) :: p
  real(kind=prec) :: trivial_spin(size(chi1))
  trivial_spin = 1.
  END FUNCTION TRIVIAL_SPIN

  PURE FUNCTION BETA_SPIN(chi, alpha, beta)
  real(kind=prec), intent(in) :: chi(:), alpha, beta
  real(kind=prec) :: beta_spin(size(chi))

  beta_spin = chi**(alpha-1) * (1-chi)**(beta-1) * gamma(alpha+beta) &
    / gamma(alpha) / gamma(beta)

  END FUNCTION BETA_SPIN

  PURE FUNCTION BETA_SPIN_11(chi1, chi2, P)
  real(kind=prec), intent(in) :: chi1(:), chi2(:)
  type(para), intent(in) :: p
  real(kind=prec) :: beta_spin_11(size(chi1))
  beta_spin_11 = beta_spin(chi1, p%alpha1, p%beta1) &
               * beta_spin(chi2, p%alpha1, p%beta1)
  END FUNCTION BETA_SPIN_11

  PURE FUNCTION BETA_SPIN_12(chi1, chi2, P)
  real(kind=prec), intent(in) :: chi1(:), chi2(:)
  type(para), intent(in) :: p
  real(kind=prec) :: beta_spin_12(size(chi1))
  beta_spin_12 = beta_spin(chi1, p%alpha1, p%beta1) &
               * beta_spin(chi2, p%alpha2, p%beta2)
  END FUNCTION BETA_SPIN_12

  PURE FUNCTION BETA_SPIN_21(chi1, chi2, P)
  real(kind=prec), intent(in) :: chi1(:), chi2(:)
  type(para), intent(in) :: p
  real(kind=prec) :: beta_spin_21(size(chi1))
  beta_spin_21 = beta_spin(chi1, p%alpha2, p%beta2) &
               * beta_spin(chi2, p%alpha1, p%beta1)
  END FUNCTION BETA_SPIN_21

                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   !!                             !!
                   !!     POWER LAW + PEAK        !!
                   !!                             !!
                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  PURE FUNCTION FLATM(M1, M2, P)
  real(kind=prec), intent(in) :: m1(:), m2(:)
  type(para), intent(in) :: p
  real(kind=prec) :: flatm(size(m1))
  flatm = 1/(m1-p%mmin)
  flatm = flatm * p%sf(m2, p%mmin, p%dm)
  END FUNCTION FLATM

  PURE FUNCTION POWM(M1, M2, P)
  real(kind=prec), intent(in) :: m1(:), m2(:)
  type(para), intent(in) :: p
  real(kind=prec) :: powm(size(m1))

  powm = (m2 / m1) ** p%k
  !powm = powm / (p%k + 1)
  powm = powm * p%sf(m2, p%mmin, p%dm)
  END FUNCTION POWM

  ! This implements the integral over S[m]*PLP[m] from zero to infinity
  PURE FUNCTION PLP_INT(p)
  type(para), intent(in) :: p
  integer, parameter :: Nsamples = 100
  integer i
  real(kind=prec), parameter :: linspace(1:Nsamples) = [(i / real(Nsamples), i=1,Nsamples)]
  real(kind=prec), dimension(1:Nsamples) :: Msample, int
  real(kind=prec) :: plp_int
  ! We write
  !
  !  /\inf                   /\mmin+dm                   /\inf
  !  |  dm S[m] * PLP[m] =   |       dm S[m] * PLP[m] +  |       dm PLP[m]
  ! \/ 0                    \/ mmin                     \/ mmin+dm
  !
  ! The first integral needs to be evaluated numerically

  Msample = p%mmin + linspace * p%dm
  int = ( (1-p%lp) * powerlaw(Msample, -p%alpha, p%mmin, p%mmax) &
              +   p%lp  * gauss(Msample, p%mum, p%sm) ) &
            * p%sf(Msample, p%mmin, p%dm)

  plp_int = sum(int) / Nsamples * p%dm

  plp_int = plp_int &
      +(1-p%lp)*(p%mmax*p%mmin**p%alpha - (p%mmax*p%mmin/(p%dm + p%mmin))**p%alpha*(p%dm + p%mmin)) &
            /(-p%mmax**p%alpha*p%mmin + p%mmax*p%mmin**p%alpha) &
      +p%lp * (Erf((p%mmax - p%mum)/Sqrt(2.)/p%sm) - Erf((p%dm + p%mmin - p%mum)/Sqrt(2.)/p%sm))/2.


  END FUNCTION PLP_INT

  ! This implements the power-law + peak model (PLP) for the primary mass of the
  ! 1st generation. The model is specified in B2 of [2010.14533]
  PURE FUNCTION PLP_MF(MBH, p)
  real(kind=prec), intent(in) :: mBH(:)
  type(para), intent(in) :: p
  real(kind=prec) :: plp_mf(size(mBH))

  plp_mf = 0.
  where ((p%mmin < mBH) .and. (mBH < p%mmax)) &
    plp_mf = ( (1-p%lp) * powerlaw(mbh, -p%alpha, p%mmin, p%mmax) &
              +   p%lp  * gauss(mBH, p%mum, p%sm) ) &
            * p%sf(mBH, p%mmin, p%dm)

  plp_mf = plp_mf / plp_int(p)

  END FUNCTION PLP_MF


  PURE FUNCTION R2P_PLP_FLAT(V) result(p)
  real(kind=prec), intent(in) :: v(:)
  type(para) :: p
  p%mmin = v(1)
  p%dm   = v(2)
  p%mmax = v(3)
  p%mum  = v(4)
  p%sm   = v(5)
  p%alpha= v(6)
  p%lp   = v(7)
  END FUNCTION R2P_PLP_FLAT

  PURE FUNCTION R2P_PLP_POW(V) result(p)
  real(kind=prec), intent(in) :: v(:)
  type(para) :: p
  p%mmin = v(1)
  p%dm   = v(2)
  p%mmax = v(3)
  p%mum  = v(4)
  p%sm   = v(5)
  p%alpha= v(6)
  p%lp   = v(7)
  p%k    = v(8)
  END FUNCTION R2P_PLP_POW

  PURE FUNCTION R2P_PLP_FLAT_BETA(V) result(p)
  real(kind=prec), intent(in) :: v(:)
  type(para) :: p
  p%mmin = v(1)
  p%dm   = v(2)
  p%mmax = v(3)
  p%mum  = v(4)
  p%sm   = v(5)
  p%alpha= v(6)
  p%lp   = v(7)

  ! Spin
  p%alpha1 = v(8)
  p%beta1 = v(9)
  END FUNCTION R2P_PLP_FLAT_BETA

  PURE FUNCTION R2P_PLP_POW_BETA(V) result(p)
  real(kind=prec), intent(in) :: v(:)
  type(para) :: p
  p%mmin = v(1)
  p%dm   = v(2)
  p%mmax = v(3)
  p%mum  = v(4)
  p%sm   = v(5)
  p%alpha= v(6)
  p%lp   = v(7)
  p%k    = v(8)

  ! Spin
  p%alpha1 = v(9)
  p%beta1 = v(10)
  END FUNCTION R2P_PLP_POW_BETA

                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   !!                             !!
                   !!     PULS. PAIR. INST. SN    !!
                   !!                             !!
                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! This implements (3) of [2104.02685]
  PURE FUNCTION PPISN_MF1G(M, P)
  real(kind=prec), intent(in) :: m(:)
  type(para), intent(in) :: p
  real(kind=prec), parameter :: mfudge = 0.99
  real(kind=prec) :: ppisn_mf1g(size(m))

  ppisn_mf1g = 0.

  where ((p%mmin < m) .and. (m < p%mgap * mfudge)) &
    ppisn_mf1g = m**p%b*&
      (1+2*p%a*p%a*sqrt(m/p%mgap)*(1-m/p%mgap)**(p%a-1))

  where ((p%mmin < m) .and. (m < p%mmin + p%dm)) &
    ppisn_mf1g = ppisn_mf1g * p%sf(m, p%mmin, p%dm)

  END FUNCTION PPISN_MF1G


  ! This implements (6) of [2104.02685]
  PURE FUNCTION PPISN_MF2G(M, P)
  real(kind=prec), intent(in) :: m(:)
  type(para), intent(in) :: p
  real(kind=prec), parameter :: mfudge = 0.99
  real(kind=prec) :: ppisn_mf2g(size(m))

  ppisn_mf2g = 1.
  where (m < p%mmin) ppisn_mf2g = 0.
  where ((p%mmin < m) .and. (m < p%mgap+ p%mmin + p%dm/2)) &
    ppisn_mf2g = p%sf(m, p%mmin, p%dm)
  where (p%mgap+ p%mmin + p%dm/2 < m) &
    ppisn_mf2g = (m / (p%mgap + p%mmin + p%dm/2)) ** p%d

  ppisn_mf2g = ppisn_mf2g * mpiv**p%b
  END FUNCTION PPISN_MF2G


  ! returns the integral from mmin to m1 of mf_1g for all m1 in lm1;
  ! this is the denominator of the conditional probability for m2
  ! given the assumption that m1 is in 1g, i.e. (1.12)
  !
  !                      (1g)
  !   /\ m            dN
  !   |      dm      -----
  !  \/ mmin   1      dm
  !                     1
  !
  PURE FUNCTION PPISN_PM2M1DEN_M11G(M, P)
  real(kind=prec), intent(in) :: m(:)
  type(para), intent(in) :: p
  real(kind=prec) :: ppisn_pm2m1den_m11g(size(m))

  integer, parameter :: Nsamples = 100
  real(kind=prec) :: BT(0:size(m))
  real(kind=prec) :: sLVC(2)
  integer i, cut
  real(kind=prec), parameter :: linspace(0:Nsamples) = [(i / real(Nsamples), i=0,Nsamples)]
  real(kind=prec), dimension(0:Nsamples) :: Msample, int
  real(kind=prec), parameter :: erf2 = erf(2._prec)
  real(kind=prec), parameter :: ep = 1e-8

  Msample = p%mmin + linspace * p%dm
  int = ppisn_mf1g(Msample, p) * p%dm / Nsamples

  do i=1, size(m)
    if(m(i) > p%mgap) cycle
    cut = findloc(Msample > m(i), .false., dim=1, back=.true.)
    ppisn_pm2m1den_m11g(i) =  sum(int(0:cut))
  enddo

  ! the part from Mmin + dm -- m1
  BT = Btilde(1.5_prec+p%b, p%a, [p%mmin+p%dm, m]/p%mgap)

  where((p%mgap > m) .and. (m > p%mmin + p%dm)) &
    ppisn_pm2m1den_m11g = ppisn_pm2m1den_m11g + &
        (m**(1+p%b) - (p%mmin+p%dm)**(1+p%b)) / (1+p%b) &
        + 2*p%a**2*p%mgap**(1+p%b) * (BT(1:size(m)) - BT(0))

  END FUNCTION PPISN_PM2M1DEN_M11G

  ! returns the integral from mmin to m1 of mf_2g for all m1 in lm1;
  ! this is the denominator of the conditional probability for m2
  ! given the assumption that m1 is in 2g
  !
  !                      (2g)
  !   /\ m            dN
  !   |      dm      -----
  !  \/ mmin   2      dm
  !                     2
  !
  PURE FUNCTION PPISN_PM2M1DEN_M12G(M, P)
  real(kind=prec), intent(in) :: m(:)
  type(para), intent(in) :: p
  real(kind=prec) :: ppisn_pm2m1den_m12g(size(m))
  real(kind=prec) :: mmax

  mmax = p%mgap + p%mmin + p%dm/2

  ppisn_pm2m1den_m12g = 0._prec

  where( (p%mmin < m).and.(m < p%mmin+p%dm) ) &
    ppisn_pm2m1den_m12g = p%sfInt(m, p%mmin, p%dm)
  where( (p%mmin+p%dm < m).and.(m < mmax) ) &
    ppisn_pm2m1den_m12g = -0.5*p%dm - p%mmin + m
  where( (mmax < m) ) &
    ppisn_pm2m1den_m12g = -0.5*p%dm - p%mmin + mmax + &
        ((m/mmax)**(1+p%d) - 1)*mmax/(1+p%d)


  ppisn_pm2m1den_m12g = ppisn_pm2m1den_m12g * mpiv**p%b

  END FUNCTION PPISN_PM2M1DEN_M12G


  PURE FUNCTION PPISN_M2_PHYS(M1, M2, P)
  real(kind=prec), intent(in) :: m1(:), m2(:)
  type(para), intent(in) :: p
  real(kind=prec) :: ppisn_m2_phys(size(m1))

  ppisn_m2_phys = ppisn_mf1g(m2, p) / ppisn_pm2m1den_m11g(m1, p)

  ! if m1 < mmin, pm2m1den_m11g = 0 and hence ppisn_m2_phys = nan
  where(isnan(ppisn_m2_phys)) &
    ppisn_m2_phys = 0.

  END FUNCTION PPISN_M2_PHYS


  ! We have
  !
  !        (1g)                  (1g)
  !      dN      +-       -+   dN      +-          -+
  ! P ~ -------- | M  | th |  -------- | M  | th,M  |
  !        dM    +- 1     -+     dM    +- 2       1-+
  !
  !                     (1g)                  (2g)
  !                   dN      +-       -+   dN      +-       -+
  !    + lam    N    -------- | M  | th |  -------- | M  | th |
  !         12   1      dM    +- 1     -+     dM    +- 2     -+
  !
  !                     (2g)                  (1g)
  !                   dN      +-       -+   dN      +-       -+
  !    + lam    N    -------- | M  | th |  -------- | M  | th |
  !         21   2      dM    +- 1     -+     dM    +- 2     -+
  !
  ! This function calculates the sum of the second and third term. We
  ! define the lam_{ij} as
  !
  ! lam       = N( M     in 2g  ) / N
  !    12(21)       2(1)             tot
  !
  ! and now need to solve for N_1(2)
  !
  !                            inf            (1g)    (2g)
  !                            /\           dN      dN
  ! N( M   in 2g) = lam   N    |  dM  dM   -----   -----     theta(M  > M )
  !     2              21  1  \/    1   2   dM      dM              1    2
  !                            0              1       2
  !
  !                            inf            (2g)    (1g)
  !                            /\           dN      dN
  ! N( M   in 2g) = lam   N    |  dM  dM   -----   -----     theta(M  > M )
  !     1              12  2  \/    1   2   dM      dM              1    2
  !                            0              1       2
  !
  ! We also need an N_0 for normalisation purposes
  !
  !       inf            (1g)    (1g)
  !       /\           dN      dN
  !  N    |  dM  dM   -----   -----     theta(M  > M )
  !   0  \/    1   2   dM      dM              1    2
  !       0              1       2
  !
  ! The integrals
  !
  !          inf        (1g)     M        (1g)
  !          /\       dN       /\ 1     dN
  ! int  =   |  dM   -----     |  dM   -----
  !    0    \/    1   dM      \/    2   dM
  !          0          1      0          2
  !
  !          inf        (1g)     M        (2g)
  !          /\       dN       /\ 1     dN
  ! int  =   |  dM   -----     |  dM   -----
  !    1    \/    1   dM      \/    2   dM
  !          0          1      0          2
  !
  !
  !          inf        (2g)     M        (1g)
  !          /\       dN       /\ 1     dN
  ! int  =   |  dM   -----     |  dM   -----
  !    2    \/    1   dM      \/    2   dM
  !          0          1      0          2
  !
  ! are returned by this function

  PURE FUNCTION PPISN_NINT(p)
  type(para), intent(in) :: p
  integer, parameter :: Nsamples = 1000
  integer, parameter :: Nit = 120
  integer i
  real(kind=prec), parameter :: linspace(1:Nsamples) = [(i / real(Nsamples), i=1,Nsamples)]
  real(kind=prec), dimension(1:NSamples) :: M1, M2, mf1g, subInt
  real(kind=prec) :: ppisn_nint(0:2)
  real(kind=prec) :: bt3,btt3(1),i1,i2,i3,i4,i5

  real(kind=prec) :: gammaFrac, gamma1ma

  btt3 = btilde(1.5 + p%b, p%a, (/(p%mmin+p%dm)/p%mgap/))
  bt3 = btt3(1)

  i1 = 0.
  gamma1ma = gamma(1-p%a)
  gammaFrac = gamma1ma
  do i=0,Nit
    btt3 = Btilde(3+2*p%b+i, p%a, (/(p%mmin+p%dm)/p%mgap/))
    if(isnan(btt3(1))) exit
    i1 = i1 + btt3(1) / (3+2*p%b+2*i) * gammaFrac
    gammaFrac = (1. - p%a/(1. + i)) * gammaFrac
  enddo
  i1 = -2*p%mgap/gamma1ma * i1 + p%mgap * bt3 * gamma(p%a) * gamma(1.5+p%b) / gamma(1.5+p%a+p%b)

  M1 = linspace * p%dm + p%mmin
  mf1g = ppisn_mf1g(m1, p)
  i2 = sum(mf1g) * p%dm / Nsamples
  i4 = sum(mf1g * lvc_int((m1 - p%mmin)/p%dm)) * p%dm / Nsamples

  do i=1,Nsamples
    subInt(i) = sum(mf1g(1:i)) * p%dm / Nsamples
  enddo
  i3 = sum(subInt * mf1g) * p%dm / Nsamples
  i5 = sum(subInt * smooth_exp(M1, p%mmin, p%dm)) * p%dm / Nsamples


  ppisn_nint(0) = (p%mgap**(1 + p%b) - (p%dm + p%mmin)**(1 + p%b))**2/(2.*(1 + p%b)**2) &
    + 4*p%a**4*p%mgap**(2 + 2*p%b)*bt3**2 + 4*p%a**4*p%mgap**(1 + 2*p%b)*i1 &
    + ((p%mgap**(1 + p%b) - (p%dm + p%mmin)**(1 + p%b))*i2)/(1 + p%b) &
    + bt3*( - 2*p%a**2*p%mgap**(1 + p%b)*i2 &
      + (-2*p%a**2*p%mgap**(1 + p%b)*(p%mgap**(1 + p%b) - (p%dm + p%mmin)**(1 + p%b))&
        )/(1 + p%b)) &
    + i3

  ppisn_nint(1) = + (p%mgap**(2 + p%b) - (p%dm + p%mmin)**(2 + p%b))/(2 + p%b) + p%dm*i4 &
    + (4*p%a**2*Sqrt(p%mgap)*(p%dm + p%mmin)**(1.5 + p%b)*(1 - (p%dm + p%mmin)/p%mgap)**p%a)/(3 + 2*p%a + 2*p%b) &
    + ((p%dm + 2*p%mmin)*(-p%mgap**(1 + p%b) + (p%dm + p%mmin)**(1 + p%b)))/(2.*(1 + p%b)) &
    + p%a**2*p%mgap**(1 + p%b)*(p%dm - (2*(3 + 2*p%b)*p%mgap)/(3 + 2*p%a + 2*p%b) + 2*p%mmin)*bt3

  ppisn_nint(2) = ((p%dm + p%mmin)**(1 + p%b)*(p%dm - p%mgap + p%mmin))/(1 + p%b) &
    - (4*p%a**2*Sqrt(p%mgap)*(p%dm + p%mmin)**(1.5 + p%b)*(1 - (p%dm + p%mmin)/p%mgap)**p%a)/(3 + 2*p%a + 2*p%b) &
    + (p%mgap**(2 + p%b) - (p%dm + p%mmin)**(2 + p%b))/((1 + p%b)*(2 + p%b)) &
    + ((p%mgap**(1 + p%b) - (p%dm + p%mmin)**(1 + p%b))*(-2*p%mgap + p%d*(p%dm + 2*p%mmin)))/(2.*(1 + p%b)*(1 + p%d)) &
    + p%a**2*p%mgap**(1 + p%b)*((-2*p%a*p%mgap)/(1.5 + p%a + p%b) + (2*p%mgap - p%d*(p%dm + 2*p%mmin))/(1 + p%d))*bt3 &
    + ((-2*p%dm - p%d*p%dm + 2*p%d*p%mgap - 2*p%mmin)*i2)/(2 + 2*p%d) + i5

  END FUNCTION PPISN_NINT

  ! norms = (/lam21, lam12 /)
  PURE FUNCTION PPISN_NORMS(D11, M1, M2, CHI1, CHI2, Z, M, P)
  real(kind=prec), intent(in) :: D11(:), m1(:), m2(:), chi1(:), chi2(:), z(:)
  type(model), intent(in) :: m
  type(para), intent(in) :: p
  real(kind=prec) :: ppisn_norms(size(m1)), N(0:2)

  real(kind=prec), dimension(size(m1)) :: D21, D12

  D21 = m%primaryM2(m1,p)*m%spin(2,1)%f(chi1, chi2,p)*m%redshift(z,p)
  D12 = m%primary  (m1,p)*m%spin(1,2)%f(chi1, chi2,p)*m%redshift(z,p)

  select case(m%secondary_c)
    case('phys')
      D21 = D21 * m%primary  (m2,p)
      D12 = D12 * m%primaryM2(m2,p)
    case('flat')
      D21 = D21 / (m1 - p%mmin)
      D12 = D12 / (m1 - p%mmin)
  end select

  ! We want
  !
  !     D                  D                  D
  !      11      \          21       \         12
  !  --------- + /\     ---------  + /\    ---------
  !   N   /N       21    N   /N        12   N   /N
  !    tot  0             tot  1             tot  2
  !
  !               \                \
  ! ~ D     N   + /\   D     N   + /\   D     N
  !    11    0      21  21    1      12  12    2

  N = ppisn_nint(p)

  ppisn_norms = D11 * N(0) + p%lam21 * D21 * N(1) + p%lam12 * D12 * N(2)
  where(isnan(ppisn_norms)) &
    ppisn_norms = 0.

  END FUNCTION PPISN_NORMS


  PURE FUNCTION R2P_PPISN(V) result(p)
  real(kind=prec), intent(in) :: v(:)
  type(para) :: p
  p%mmin = v(1)
  p%dm   = v(2)
  p%mgap = v(3)
  p%a    = v(4)
  p%b    = v(5)
  p%d    = v(6)
  p%lam21= 10**v(7)
  p%lam12= 10**v(8)
  END FUNCTION R2P_PPISN


  PURE FUNCTION R2P_PPISN_BETA(V) result(p)
  real(kind=prec), intent(in) :: v(:)
  type(para) :: p
  p%mmin = v(1)
  p%dm   = v(2)
  p%mgap = v(3)
  p%a    = v(4)
  p%b    = v(5)
  p%d    = v(6)
  p%lam21= 10**v(7)
  p%lam12= 10**v(8)

  ! Spin
  p%alpha1 = v(9)
  p%beta1 = v(10)
  p%alpha2 = v(11)
  p%beta2 = v(12)
  END FUNCTION R2P_PPISN_BETA


                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   !!                             !!
                   !!            TOOLS            !!
                   !!                             !!
                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION GETMODEL(mods) result(m)
  character(len=*), intent(in) :: mods
  type(model) :: m

  m%model_name = trim(mods)
  select case(mods)
    case('plp+flat+trivial+trivial')
      m%ndim = 7
      m%primary => plp_mf
      m%secondary => flatm
      m%redshift => trivial
      m%spin(1,1)%f => trivial_spin
      m%spin(1,2)%f => trivial_spin
      m%spin(2,1)%f => trivial_spin
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_plp_flat
      m%smooth => smooth_tanh
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
      m%norms = .false.
    case('plp+pow+trivial+trivial')
      m%ndim = 8
      m%primary => plp_mf
      m%secondary => powm
      m%redshift => trivial
      m%spin(1,1)%f => trivial_spin
      m%spin(1,2)%f => trivial_spin
      m%spin(2,1)%f => trivial_spin
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_plp_pow
      m%smooth => smooth_exp
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
      m%norms = .false.

    case('plp+flat+trivial+beta')
      m%ndim = 7
      m%primary => plp_mf
      m%secondary => flatm
      m%redshift => trivial
      m%spin(1,1)%f => beta_spin_11
      m%spin(1,2)%f => trivial_spin
      m%spin(2,1)%f => trivial_spin
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_plp_flat_beta
      m%smooth => smooth_tanh
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
      m%norms = .false.
    case('plp+pow+trivial+beta')
      m%ndim = 8
      m%primary => plp_mf
      m%secondary => powm
      m%redshift => trivial
      m%spin(1,1)%f => beta_spin_11
      m%spin(1,2)%f => trivial_spin
      m%spin(2,1)%f => trivial_spin
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_plp_pow_beta
      m%smooth => smooth_exp
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
      m%norms = .false.

    case("ppisn+flat+trivial+trivial")
      m%ndim = 8
      m%primary => ppisn_mf1g
      m%primaryM2 => ppisn_mf2g
      m%secondary => flatm
      m%secondary_c = "flat"
      m%redshift => trivial
      m%spin(1,1)%f => trivial_spin
      m%spin(1,2)%f => trivial_spin
      m%spin(2,1)%f => trivial_spin
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_ppisn
      m%smooth => smooth_exp
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
      m%norms = .true.
    case("ppisn+trivial+trivial")
      m%ndim = 8
      m%primary => ppisn_mf1g
      m%primaryM2 => ppisn_mf2g
      m%secondary => ppisn_m2_phys
      m%secondary_c = "phys"
      m%redshift => trivial
      m%spin(1,1)%f => trivial_spin
      m%spin(1,2)%f => trivial_spin
      m%spin(2,1)%f => trivial_spin
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_ppisn
      m%smooth => smooth_exp
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
      m%norms = .true.

    case("ppisn+trivial+beta")
      m%ndim = 8
      m%primary => ppisn_mf1g
      m%primaryM2 => ppisn_mf2g
      m%secondary => ppisn_m2_phys
      m%secondary_c = "phys"
      m%redshift => trivial
      m%spin(1,1)%f => beta_spin_11
      m%spin(1,2)%f => beta_spin_12
      m%spin(2,1)%f => beta_spin_21
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_ppisn_beta
      m%smooth => smooth_exp
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
      m%norms = .true.
    case default
      stop 9
  end select
  END FUNCTION GETMODEL

  SUBROUTINE TEST
  real(kind=prec), dimension(6) :: mtest, ans
  real(kind=prec) :: diff
  type(para) :: p
  real(kind=prec) :: BT(0:size(mtest))
  type(model) :: the_model

  mtest = (/ 1.,11., 20., 30., 40., 100. /)
  ans = (/0._prec,  0.010459028944395847_prec, &
        0.0035206207261596584_prec,  0.0017915249285508827_prec, &
        0.0012500000000000002_prec,  0._prec /)

  diff = sum(abs(ans - ppisn_mf1g(mtest, &
                para(mgap = 50._prec, &
                     a    =  .5_prec, &
                     b    = -2._prec, &
                     d    = -3._prec, &
                     mmin =  5._prec, &
                     dm   =  5._prec, &
                     lam12=  0._prec, &
                     lam21=  0._prec, &
                     sf_c = 'exp'   , &
                     sfint= smooth_expint, &
                     sf   = smooth_exp)))) / 6.
  print*, "ppisn_mf1g", diff

  ans = (/ 0._prec, 1/900._prec, 1/900._prec, &
           1/900._prec, 1/900._prec, 0.00021123263888888892_prec/)
  diff = sum(abs(ans - ppisn_mf2g(mtest, &
                para(mgap = 50._prec, &
                     a    =  .5_prec, &
                     b    = -2._prec, &
                     d    = -3._prec, &
                     mmin =  5._prec, &
                     dm   =  5._prec, &
                     lam12=  0._prec, &
                     lam21=  0._prec, &
                     sf_c = 'exp'   , &
                     sfint= smooth_expint, &
                     sf   = smooth_exp)))) / 6.
  print*, "ppisn_mf2g", diff
  print*,plp_int(para(mmax = 50._prec, &
                     mum  = 30._prec, &
                     sm   = 5._prec, &
                     alpha= 2._prec, &
                     lp   = 0.2_prec, &
                     mmin = 5._prec, &
                     dm   = 5._prec, &
                     lam12=  0._prec, &
                     lam21=  0._prec, &
                     sf_c = 'tan'   , &
                     sfint= smooth_expint, &
                     sf   = smooth_tanh))

  ans = (/ 0., 0.03674262, 0.01327075, 0.02089596, 0.00493742, 0. /)

  diff = sum(abs(ans - plp_mf(mtest, &
                para(mmax = 50._prec, &
                     mum  = 30._prec, &
                     sm   = 5._prec, &
                     alpha= 2._prec, &
                     lp   = 0.2_prec, &
                     mmin = 5._prec, &
                     dm   = 5._prec, &
                     lam12=  0._prec, &
                     lam21=  0._prec, &
                     sf_c = 'tan'   , &
                     sfint= smooth_expint, &
                     sf   = smooth_tanh)))) / 6.
  print*,"plp_mf", diff

  mtest = (/ 2., 3., 4., 5., 10., 20. /)
  ans = (/ 0., 3.14309363e-13, 5.66146660e-06, &
           1.21599290e-04, 3.97884590e-03, 1.30775022e-02 /)

  diff = sum(abs(ans - plp_mf(mtest, &
                para(mmax = 39.3856_prec,&
                     mum  = 44.7844_prec,&
                     sm   = 6.06597_prec,&
                     alpha= -1.71307_prec,&
                     lp   = 0.394535_prec,&
                     mmin = 2.61592_prec,&
                     dm   = 8.5451_prec,&
                     lam12=  0._prec, &
                     lam21=  0._prec, &
                     sf_c = 'tan'   , &
                     sfint= smooth_expint, &
                     sf   = smooth_tanh)))) / 6.

  print*,"plp_mf", diff

  mtest = (/ 2., 4., 10., 20., 25., 35. /)
  p = para(mmin = 2.61592_prec,&
           dm   = 8.5451_prec,&
           mgap = 34._prec, &
           a    = 0.23_prec, &
           b    = -2.34_prec, &
           d    = -7.63_prec, &
           lam21 = 0._prec, &
           lam12 = 0._prec, &
           sf   = smooth_erf, &
           sfint= smooth_erfint, &
           sf_c = 'erf')

  ans = (/0._prec, 0.0011329932907377792_prec, &
          0.004816986990970342_prec, 0.001047864239495318_prec,&
          0.0006707850448342742_prec, 0._prec /)
  diff = sum(abs(ans-ppisn_mf1g(mtest, p))) / 6.
  print*, "ppisn_mf1g", diff

  ans = (/ 0._prec, 0.0008369386215617741_prec, &
        0.028895302280797192_prec, 0.05158740716028034_prec, &
        0.0557564000458822_prec, 0.0_prec /)
  diff = sum((ppisn_pm2m1den_m11g(mtest, p) - ans) / (ans + 1e-15))/ 6.
  print*, "ppisn_pm2m1den_m11g", diff

  mtest = (/ 3., 30., 35., 37., 40., 45. /)
  ans = (/ 0.017562827691927917e-4_prec, 3.4957169791193885e-4_prec, &
          3.4957169791193885e-4_prec, 3.4957169791193885e-4_prec, &
          3.4957169791193885e-4_prec, 1.6828265142804126e-4_prec /)

  diff = sum(abs(ans-ppisn_mf2g(mtest, p))) / 6.
  print*, "ppisn_mf2g", diff

  mtest = (/ 2._prec, 4._prec, 10._prec, 20._prec, 25._prec, 35._prec /)
  ans = (/0._prec, 5.290164443318091e-6_prec, &
          0.0010911514507551488_prec, 0.004583419804323329_prec, &
          0.006331278293883026_prec, 0.009826995273002417_prec/)

  diff = sum(abs((ppisn_pm2m1den_m12g(mtest, p)+1e-15) / (ans+1e-15) - 1)) / 6.
  print*, "ppisn_pm2m1den_m12g", diff

  print*, "Comparison with python code"
  mtest = (/ 2., 4., 10., 20., 25., 35. /)
  p%sf => smooth_exp
  p%sfint => smooth_expint
  p%sf_c = 'exp'

  ans = (/0._prec, 0.000276814962771_prec, 0.004903900506881_prec, &
         0.001047864239495_prec, 0.000670785044834_prec, 0._prec /)
  diff = sum(abs(ans-ppisn_mf1g(mtest, p))) / 6.
  print*, "ppisn_mf1g", diff

  ans = (/2.167636805288472e-13, 3.495716979119389e-04, &
          3.495716979119389e-04, 3.495716979119389e-04, &
          3.495716979119389e-04, 1.682826514280412e-04 /)
  mtest = (/ 3., 30., 35., 37., 40., 45. /)
  diff = sum(abs(ans-ppisn_mf2g(mtest, p))) / 6.
  print*, "ppisn_mf2g", diff

  ans = (/ -7.036244987357108,  3.366450242210134,  4.699203519611515, &
            5.757592699486201,  9.013019845715352, 29.836942133255924/)
  BT = Btilde(1.5_prec+p%b, p%a, [p%mmin+p%dm, mtest]/p%mgap)
  diff = sum(abs(ans-(BT(1:size(mtest)) - BT(0)))) / 6.
  print*, 'Btilde', diff

  p%sf => smooth_exp
  ans(1:3) = ppisn_nint(p)
  diff = sum(ans(1:3) / (/0.0018863614613710922,0.4113173426155685,2.05585086404557/) - 1.) / 3
  print*, 'norm int', diff

  ! comparison with gwpopulation
  the_model = getmodel("plp+pow+trivial+trivial")
  p = the_model%r2p((/3.5_prec, 5._prec, 85._prec, 33._prec, 3._prec, 3.0_prec, 0.03_prec, 4._prec/))
  p%sf => the_model%smooth
  p%sfint => the_model%smoothint
  p%sf_c = the_model%smooth_c
  ans = plp_mf(mtest, p)
  print*, ans(2:) / (/0.00013537551253735636_prec, 0.06283120273729201_prec, &
              0.007854781123313786_prec,   0.004321979423810812_prec, &
              0.009896871567398238_prec /)

  END SUBROUTINE TEST

                          !!!!!!!!!!!!!!!!!!!!!!
                            END MODULE MODELS
                          !!!!!!!!!!!!!!!!!!!!!!
