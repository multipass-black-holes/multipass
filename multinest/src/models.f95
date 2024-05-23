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
    real(kind=prec) :: mgap=0, a=0, b=0, d=0, bq0=0, bq1=0

    ! Needed for spin
    real(kind=prec) :: alpha1=0, alpha2=0, beta1=0, beta2=0

    ! Smooth function
    procedure(smoothfn), pointer, nopass :: sf, sfint
    character(len=3) :: sf_c

    ! Needed for H0
    real(kind=prec) :: H0 = 100.
    real(kind=prec) :: gamma = 3

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

    PURE SUBROUTINE REDSHIFTFN(M1D, M2D, M1S, M2S, Z, D, P, LL)
      use functions, only: prec
      import para
      real(kind=prec), intent(in) :: m1d(:), m2d(:), Z(:), D(:)
      real(kind=prec), intent(out) :: m1s(:), m2s(:), ll(:)
      type(para), intent(in) :: p
    END SUBROUTINE REDSHIFTFN

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

    PURE FUNCTION CUTSFN(P)
      use functions, only: prec
      import para
      type(para), intent(in) :: p
      logical cutsfn
    END FUNCTION CUTSFN
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
    procedure(cutsFn), pointer, nopass :: cuts
    procedure(smoothFn), pointer, nopass :: smooth, smoothInt
    character(len=3) :: smooth_c
    character(len=4) :: secondary_c
    logical :: norms
  END TYPE MODEL

contains

  PURE FUNCTION NOCUTS(P)
  type(para), intent(in) :: p
  logical nocuts
  nocuts = .true.
  END FUNCTION NOCUTS

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

  beta_spin = 0.
  where ((chi > 0.) .and. (chi < 1.)) &
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



  PURE FUNCTION GAUSS_SPIN(chi, alpha, beta)
  real(kind=prec), intent(in) :: chi(:), alpha, beta
  real(kind=prec) :: gauss_spin(size(chi))

  gauss_spin = 0.
  where ((chi > 0.) .and. (chi < 1.)) &
    gauss_spin = gauss(chi, alpha, beta)

  END FUNCTION GAUSS_SPIN

  PURE FUNCTION GAUSS_SPIN_11(chi1, chi2, P)
  real(kind=prec), intent(in) :: chi1(:), chi2(:)
  type(para), intent(in) :: p
  real(kind=prec) :: gauss_spin_11(size(chi1))
  gauss_spin_11 = gauss_spin(chi1, p%alpha1, p%beta1) &
               * gauss_spin(chi2, p%alpha1, p%beta1)
  END FUNCTION GAUSS_SPIN_11

  PURE FUNCTION GAUSS_SPIN_12(chi1, chi2, P)
  real(kind=prec), intent(in) :: chi1(:), chi2(:)
  type(para), intent(in) :: p
  real(kind=prec) :: gauss_spin_12(size(chi1))
  gauss_spin_12 = gauss_spin(chi1, p%alpha1, p%beta1) &
               * gauss_spin(chi2, p%alpha2, p%beta2)
  END FUNCTION GAUSS_SPIN_12

  PURE FUNCTION GAUSS_SPIN_21(chi1, chi2, P)
  real(kind=prec), intent(in) :: chi1(:), chi2(:)
  type(para), intent(in) :: p
  real(kind=prec) :: gauss_spin_21(size(chi1))
  gauss_spin_21 = gauss_spin(chi1, p%alpha2, p%beta2) &
               * gauss_spin(chi2, p%alpha1, p%beta1)
  END FUNCTION GAUSS_SPIN_21


                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   !!                             !!
                   !!          RED SHIFT          !!
                   !!                             !!
                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE REDSHIFT_PLANCK(M1D, M2D, M1S, M2S, Z, D, P, LL)
  ! This implements (4.1) of 2205.11278 to calculate the source-frame
  ! masses from the detector-frame ones.
  !    m_source = m_det / (1+z)
  ! We infer the redshift z from the luminosity distance d rather than
  ! taking the LVC value (which assumed a value of H0).
  use functions, only: prec
  real(kind=prec), intent(in) :: m1d(:), m2d(:), Z(:), D(:)
  real(kind=prec), intent(out) :: m1s(:), m2s(:), LL(:)
  real(kind=prec) :: zz(size(d))
  type(para), intent(in) :: p
  real(kind=prec), parameter :: Om = 0.3111
  real(kind=prec), parameter :: conv = 299792.458 ! speed of light in km/s

#if 0
  integer, parameter :: Nnewton = 10
  real(kind=prec) :: Btild0(1), pref
  real(kind=prec) :: Dl(size(d)), dDl(size(d))
  integer :: i

  Btild0 = Btilde(1./6._prec, 1./3._prec, (/1. - Om/))
  pref = 1/3./p%H0/Om**(1/3.)/(1-Om)**(1/6.)
  zz = 1.
  do i=1,Nnewton
    Dl = pref*(1 + zz)*(Btild0(1) - Btilde(1./6._prec,1./3._prec,(1 - Om)/(1 - Om + Om*(1 + zz)**3)))
    dDl = Dl/(1 + zz) + (1 + zz)/(p%H0*Sqrt(1 + Om*(-1 + (1 + zz)**3)))

    zz = zz - (Dl-d)/dDL
  enddo

#else
  integer, parameter :: n = 101
  real(kind=prec), parameter :: tab(n) =  (/&
      +4.531037182232878e+0_prec, +4.153820161637109e+0_prec, -1.948963677061516e-1_prec, &
      +7.725543272874448e-2_prec, -3.972005516372123e-2_prec, +2.301900172743636e-2_prec, &
      -1.420795885288088e-2_prec, +9.101653707715140e-3_prec, -5.975505612039882e-3_prec, &
      +3.994371132444885e-3_prec, -2.708612173574250e-3_prec, +1.859027363293421e-3_prec, &
      -1.289395364094574e-3_prec, +9.026654802070240e-4_prec, -6.372000230371466e-4_prec, &
      +4.531663079915926e-4_prec, -3.244439107824475e-4_prec, +2.336820390773077e-4_prec, &
      -1.692197344832828e-4_prec, +1.231354550873643e-4_prec, -8.999446114561720e-5_prec, &
      +6.603407259755666e-5_prec, -4.862743197982891e-5_prec, +3.592661972759330e-5_prec, &
      -2.662262832557941e-5_prec, +1.978233141559833e-5_prec, -1.473672584251300e-5_prec, &
      +1.100366469317950e-5_prec, -8.234014388133980e-6_prec, +6.173862318468393e-6_prec, &
      -4.637814010026728e-6_prec, +3.490011036218634e-6_prec, -2.630562280508878e-6_prec, &
      +1.985795341703480e-6_prec, -1.501219857745157e-6_prec, +1.136426537214736e-6_prec, &
      -8.613747896730620e-7_prec, +6.536810991460079e-7_prec, -4.966316801954343e-7_prec, &
      +3.777214048161198e-7_prec, -2.875765163501911e-7_prec, +2.191581588077745e-7_prec, &
      -1.671720455257976e-7_prec, +1.276298638628562e-7_prec, -9.752271490446970e-8_prec, &
      +7.457744024174252e-8_prec, -5.707447072350718e-8_prec, +4.371138264098474e-8_prec, &
      -3.350054500423107e-8_prec, +2.569219200565986e-8_prec, -1.971652894195755e-8_prec, &
      +1.514009343857410e-8_prec, -1.163281742225450e-8_prec, +8.943129815480940e-9_prec, &
      -6.879122608329922e-9_prec, +5.294275898965282e-9_prec, -4.076634110589088e-9_prec, &
      +3.140585549805561e-9_prec, -2.420615107567923e-9_prec, +1.866550992949891e-9_prec, &
      -1.439946291738698e-9_prec, +1.111318426075016e-9_prec, -8.580464271909640e-10_prec, &
      +6.627613849546026e-10_prec, -5.121187903805117e-10_prec, +3.958658563584725e-10_prec, &
      -3.061149374741611e-10_prec, +2.367952375351237e-10_prec, -1.832359700021200e-10_prec, &
      +1.418383431486512e-10_prec, -1.098283066813663e-10_prec, +8.506960733731040e-11_prec, &
      -6.591235040098587e-11_prec, +5.108356185650173e-11_prec, -3.960264451755628e-11_prec, &
      +3.071085019365777e-11_prec, -2.382202861839671e-11_prec, +1.848315746583550e-11_prec, &
      -1.434464700231436e-11_prec, +1.113535279679639e-11_prec, -8.646271680583760e-12_prec, &
      +6.715701449325460e-12_prec, -5.216446657850302e-12_prec, +4.052702132802613e-12_prec, &
      -3.149132685785899e-12_prec, +2.448446784511657e-12_prec, -1.903530061642417e-12_prec, &
      +1.479270071207392e-12_prec, -1.149785281681476e-12_prec, +8.943695918415400e-13_prec, &
      -6.945628773341674e-13_prec, +5.385960156805607e-13_prec, -4.179446773537620e-13_prec, &
      +3.226805113955488e-13_prec, -2.485542557473826e-13_prec, +1.891098111244261e-13_prec, &
      -1.429951121380522e-13_prec, +1.047123118977761e-13_prec, -7.319415388656947e-14_prec, &
      +4.658437470735631e-14_prec, -2.220473529848197e-14_prec /)
  real(kind=prec) :: y(size(d)), dmat(size(d),n)
  integer j

  do j=1,size(d)
    dmat(j, :) = tab
  enddo
  y = (p%H0 * d / conv) / 10 - 1

  do j = n, 3, -1
    dmat(:, j-1) = dmat(:, j-1) + 2*y(:)*dmat(:, j)
    dmat(:, j-2) = dmat(:, j-2) - dmat(:, j)
  end do
  zz = dmat(:, 1)+y(:)*dmat(:, 2)
#endif

  ll = d**2 * (1+zz)**(p%gamma-3) / p%H0 / sqrt((1+zz)**3*Om + (1-Om))

  m1s = m1d / (1+zz)
  m2s = m2d / (1+zz)

  END SUBROUTINE REDSHIFT_PLANCK


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

  PURE FUNCTION PLP_M2F(m1,m2, p)
  real(kind=prec), intent(in) :: m1(:),m2(:)
  type(para), intent(in) :: p
  real(kind=prec) :: plp_m2f(size(m2))

  plp_m2f = 0.
  where ((p%mmin < m2) .and. (m2 < p%mmax)) &
    plp_m2f = ( (1-p%lp) * powerlaw(m2, -p%alpha, p%mmin, p%mmax) &
              +   p%lp  * gauss(m2, p%mum, p%sm) ) &
            * p%sf(m2, p%mmin, p%dm) &
            * (m2/m1) ** p%k

  plp_m2f = plp_m2f / plp_int(p)

  END FUNCTION PLP_M2F


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

  PURE FUNCTION R2P_PLP_FLAT_PLANCK(V) result(p)
  real(kind=prec), intent(in) :: v(:)
  type(para) :: p
  p%mmin = v(1)
  p%dm   = v(2)
  p%mmax = v(3)
  p%mum  = v(4)
  p%sm   = v(5)
  p%alpha= v(6)
  p%lp   = v(7)
  p%h0   = v(8)
  END FUNCTION R2P_PLP_FLAT_PLANCK

  PURE FUNCTION R2P_PLP_POW_PLANCK(V) result(p)
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
  p%h0   = v(9)
  END FUNCTION R2P_PLP_POW_PLANCK

                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   !!                             !!
                   !!     PULS. PAIR. INST. SN    !!
                   !!                             !!
                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE FUNCTION BD_CUT(P)
  type(para), intent(in) :: p
  logical bd_cut
  bd_cut = p%d < p%b
  END FUNCTION BD_CUT

  PURE FUNCTION BD_LAM_CUT(P)
  type(para), intent(in) :: p
  logical bd_lam_cut
  bd_lam_cut = (p%d < p%b) .and. (p%lam21 > p%lam12)
  END FUNCTION BD_LAM_CUT

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


  PURE FUNCTION PPISN_M2_PHYS(M1, M2, P)
  real(kind=prec), intent(in) :: m1(:), m2(:)
  type(para), intent(in) :: p
  real(kind=prec) :: ppisn_m2_phys(size(m1))

  ppisn_m2_phys = ppisn_mf1g(m2, p)

  ! if m1 < mmin, pm2m1den_m11g = 0 and hence ppisn_m2_phys = nan
  where(isnan(ppisn_m2_phys)) &
    ppisn_m2_phys = 0.

  END FUNCTION PPISN_M2_PHYS


  ! We have
  !
  !        (1g)                  (1g)
  !      dN      +-       -+   dN      +-          -+    bq0
  ! P ~ -------- | M  | th |  -------- | M  | th,M  |   q
  !        dM    +- 1     -+     dM    +- 2       1-+
  !
  !                     (1g)                  (2g)
  !                   dN      +-       -+   dN      +-       -+   bq1
  !    + lam    N    -------- | M  | th |  -------- | M  | th |  q
  !         12   1      dM    +- 1     -+     dM    +- 2     -+
  !
  !                     (2g)                  (1g)
  !                   dN      +-       -+   dN      +-       -+   bq1
  !    + lam    N    -------- | M  | th |  -------- | M  | th |  q
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
  !                            /\           dN      dN                        bq1
  ! N( M   in 2g) = lam   N    |  dM  dM   -----   -----     theta(M  > M )  q
  !     2              21  1  \/    1   2   dM      dM              1    2
  !                            0              1       2
  !
  !                            inf            (2g)    (1g)
  !                            /\           dN      dN                        bq1
  ! N( M   in 2g) = lam   N    |  dM  dM   -----   -----     theta(M  > M )  q
  !     1              12  2  \/    1   2   dM      dM              1    2
  !                            0              1       2
  !
  ! We also need an N_0 for normalisation purposes
  !
  !       inf            (1g)    (1g)
  !       /\           dN      dN                        bq0
  !  N    |  dM  dM   -----   -----     theta(M  > M )  q
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
  real(kind=prec) :: bt3,bt3m0,bt3p0,bt3m1,bt3p1,bt5,btt(1)
  real(kind=prec) :: i1,i2,i3,i4,i5

  real(kind=prec) :: gammaFrac, gamma1ma

  btt = btilde(1.5 + p%b, p%a, (/(p%mmin+p%dm)/p%mgap/)) ; bt3 = btt(1)
  btt = btilde(2.5 + 2*p%b, p%a, (/(p%mmin+p%dm)/p%mgap/)) ; bt5 = btt(1)
  btt = btilde(1.5 + p%b - p%bq0, p%a, (/(p%mmin+p%dm)/p%mgap/)) ; bt3m0 = btt(1)
  btt = btilde(1.5 + p%b - p%bq1, p%a, (/(p%mmin+p%dm)/p%mgap/)) ; bt3m1 = btt(1)
  btt = btilde(1.5 + p%b + p%bq0, p%a, (/(p%mmin+p%dm)/p%mgap/)) ; bt3p0 = btt(1)
  btt = btilde(1.5 + p%b + p%bq1, p%a, (/(p%mmin+p%dm)/p%mgap/)) ; bt3p1 = btt(1)

  i1 = 0.
  gamma1ma = gamma(1-p%a)
  gammaFrac = gamma1ma
  do i=0,Nit
    btt = Btilde(3+2*p%b+i, p%a, (/(p%mmin+p%dm)/p%mgap/))
    if(isnan(btt(1))) exit
    i1 = i1 + btt(1) / (3+2*p%b+2*i+2*p%bq0) * gammaFrac
    gammaFrac = (1. - p%a/(1. + i)) * gammaFrac
  enddo
  i1 = -2*p%mgap/gamma1ma * i1 + p%mgap * bt3m0 * gamma(p%a) * gamma(1.5+p%b+p%bq0) / gamma(1.5+p%a+p%b+p%bq0)

  M1 = linspace * p%dm + p%mmin
  mf1g = ppisn_mf1g(m1, p)
  i2 = sum(mf1g) * p%dm / Nsamples
  i4 = sum(mf1g * lvc_int((m1 - p%mmin)/p%dm)) * p%dm / Nsamples

  do i=1,Nsamples
    subInt(i) = sum(mf1g(1:i)) * p%dm / Nsamples
  enddo
  i3 = sum(subInt * mf1g) * p%dm / Nsamples
  i5 = sum(subInt * smooth_exp(M1, p%mmin, p%dm)) * p%dm / Nsamples


  ppisn_nint(0) = ((1 + p%b - p%bq0)*p%mgap**(2 + 2*p%b) &
      + (1 + p%b + p%bq0)*p%mgap**(2 + 2*p%b)*((p%dm + p%mmin)/p%mgap)**(2 + 2*p%b) &
      - 2*(1 + p%b)*p%mgap**(2 + 2*p%b)*((p%dm + p%mmin)/p%mgap)**(1 + p%b + p%bq0)&
      )/(2.*(1 + p%b)*(1 + p%b - p%bq0)*(1 + p%b + p%bq0)) &
    - (2*p%a**2*p%mgap**(2 + 2*p%b)*bt3p0)/(1 + p%b - p%bq0) &
    + bt3m0*((2*p%a**2*p%mgap**(2 + 2*p%b)*((p%dm + p%mmin)/p%mgap)**(1 + p%b + p%bq0))&
      /(1 + p%b + p%bq0) + 4*p%a**4*p%mgap**(2 + 2*p%b)*bt3p0) &
    + (4*p%a**2*p%bq0*p%mgap**(2 + 2*p%b)*bt5)/(1 + 2*p%b + p%b**2 - p%bq0**2) &
    + 4*p%a**4*p%mgap**(1 + 2*p%b)*i1 &
    + ((p%mgap**(1 + p%b) - p%mgap**(1 + p%b)*((p%dm + p%mmin)/p%mgap)**(1 + p%b))/(1 + p%b) &
      - 2*p%a**2*p%mgap**(1 + p%b)*bt3)*i2 &
    + i3

  ppisn_nint(1) = (4*p%a**2*p%mgap**(2 + p%b)*((p%dm + p%mmin)/p%mgap)**(1.5 + p%b) &
      *(1 - (p%dm + p%mmin)/p%mgap)**p%a)/((3 + 2*p%a + 2*p%b)*(1 + p%bq1)) &
    + (p%dm*(p%mgap**(1 + p%b) - p%mgap**(1 + p%b)*((p%dm + p%mmin)/p%mgap)**(1 + p%b)))/(2.*(1 + p%b)) &
    + (p%mgap**(2 + p%b) - p%mgap**(2 + p%b)*((p%dm + p%mmin)/p%mgap)**(2 + p%b))/((2 + p%b)*(1 + p%bq1)) &
    - (p%mgap**(1 + p%bq1)*((p%dm + p%mmin)/p%mgap)**(1 + p%bq1) &
      *(p%mgap**(1 + p%b - p%bq1) - p%mgap**(1 + p%b - p%bq1)*((p%dm + p%mmin)/p%mgap)**(1 + p%b - p%bq1))&
      )/((1 + p%b - p%bq1)*(1 + p%bq1)) &
    + p%a**2*p%mgap**(1 + p%b)*(-p%dm - (2*(3 + 2*p%b)*p%mgap)/((3 + 2*p%a + 2*p%b)*(1 + p%bq1)))*bt3 &
    + (2*p%a**2*p%mgap**(2 + p%b)*((p%dm + p%mmin)/p%mgap)**(1 + p%bq1)*bt3m1)/(1 + p%bq1) &
    + p%dm*i4

  ppisn_nint(2) = (4*p%a**2*p%mgap**(2 + p%b)*((p%dm + p%mmin)/p%mgap)**(1.5 + p%b) &
      *(1 - (p%dm + p%mmin)/p%mgap)**p%a)/((3 + 2*p%a + 2*p%b)*(-1 + p%bq1)) &
    + ((2*(p%mgap**(2 + p%b) - p%mgap**(2 + p%b)*((p%dm + p%mmin)/p%mgap)**(2 + p%b)))/(2 + p%b) &
      + 2*(-p%mgap**(2 + p%b + p%bq1)*((p%dm + p%mmin)/p%mgap)**(2 + p%b) &
        + p%mgap**(2 + p%b + p%bq1)*((p%dm + p%mmin)/p%mgap)**(1 + p%b + p%bq1)) &
      /((-1 + p%bq1)*p%mgap**p%bq1) &
      + ((p%mgap**(1 + p%b + p%bq1) &
          - p%mgap**(1 + p%b + p%bq1)*((p%dm + p%mmin)/p%mgap)**(1 + p%b + p%bq1))&
        *(2**(2*p%bq1)*p%d*p%mgap**(2*p%bq1)*(p%dm + 2*p%mmin)*(1 + (p%dm + 2*p%mmin)/(2.*p%mgap))**p%bq1 &
          + 2*p%mgap*(2**(2*p%bq1)*p%d*p%mgap**(2*p%bq1)*(1 + (p%dm + 2*p%mmin)/(2.*p%mgap))**p%bq1 &
          + 2**(2*p%bq1)*(-1 + p%bq1)*p%mgap**(2*p%bq1)*(1 + (p%dm + 2*p%mmin)/(2.*p%mgap))**(2*p%bq1) &
          - 2**(2*p%bq1)*p%d*p%mgap**(2*p%bq1)*(1 + (p%dm + 2*p%mmin)/(2.*p%mgap))**(2*p%bq1)))&
        )/(4**p%bq1*(-1 + p%bq1)*(-1 + p%bq1 - p%d)*p%mgap**(3*p%bq1)*(1 + (p%dm + 2*p%mmin)/(2.*p%mgap))**(2*p%bq1))&
      )/(2.*(1 + p%b + p%bq1)) &
    - (2*p%a**2*(3 + 2*p%b)*p%mgap**(2 + p%b)*bt3)/((3 + 2*p%a + 2*p%b)*(-1 + p%bq1)) &
    - (2*p%a**2*p%d*p%mgap**(2 + p%b)*(1 + (p%dm + 2*p%mmin)/(2.*p%mgap))**(1 - p%bq1)*bt3p1)&
      /(1 - 2*p%bq1 + p%bq1**2 + p%d - p%bq1*p%d) &
    + ((-2*p%dm - p%d*p%dm + 2*p%d*p%mgap - 2*p%mmin)*i2)/(2 + 2*p%d) &
    + i5

  END FUNCTION PPISN_NINT

  ! norms = (/lam21, lam12 /)
  PURE FUNCTION PPISN_NORMS(D11, M1, M2, CHI1, CHI2, M, P)
  real(kind=prec), intent(in) :: D11(:), m1(:), m2(:), chi1(:), chi2(:)
  type(model), intent(in) :: m
  type(para), intent(in) :: p
  real(kind=prec) :: ppisn_norms(size(m1)), N(0:2)

  real(kind=prec), dimension(size(m1)) :: D21, D12

  D21 = m%primaryM2(m1,p)*m%spin(2,1)%f(chi1, chi2,p)
  D12 = m%primary  (m1,p)*m%spin(1,2)%f(chi1, chi2,p)

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

  ppisn_norms = D11 / N(0) * (m2/m1)**p%bq0 &
    + (p%lam21 * D21 / N(1) + p%lam12 * D12 / N(2)) * (m2/m1)**p%bq1
  where(isnan(ppisn_norms)) &
    ppisn_norms = 0.

  END FUNCTION PPISN_NORMS


  PURE FUNCTION R2P_PPISN_NOLAM(V) result(p)
  real(kind=prec), intent(in) :: v(:)
  type(para) :: p
  p%mmin = v(1)
  p%dm   = v(2)
  p%mgap = v(3)
  p%a    = v(4)
  p%b    = v(5)
  p%d    = v(6)
  p%lam21= 0
  p%lam12= 0
  p%bq0  = v(7)
  p%bq1  = v(7)
  END FUNCTION R2P_PPISN_NOLAM

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
  p%bq0  = v(9)
  p%bq1  = v(9)
  END FUNCTION R2P_PPISN

  PURE FUNCTION R2P_PPISN_TWOPAIR(V) result(p)
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
  p%bq0  = v(9)
  p%bq1  = v(10)
  END FUNCTION R2P_PPISN_TWOPAIR


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
  p%bq0  = v(9)
  p%bq1  = v(9)

  ! Spin
  p%alpha1 = v(10)
  p%beta1 = v(11)
  p%alpha2 = v(12)
  p%beta2 = v(13)
  END FUNCTION R2P_PPISN_BETA

  PURE FUNCTION R2P_PPISN_BETA_NO_TURN_ON(V) result(p)
  real(kind=prec), intent(in) :: v(:)
  type(para) :: p
  p%mmin = 4.02_prec
  p%dm   = 5.27_prec

  p%mgap =     v(1)
  p%a    =     v(2)
  p%b    =     v(3)
  p%d    =     v(4)
  p%lam21= 10**v(5)
  p%lam12= 10**v(6)
  p%bq0  =     v(7)
  p%bq1  =     v(7)
  ! Spin
  p%alpha1 = v(8)
  p%beta1 =  v(9)
  p%alpha2 = v(10)
  p%beta2 =  v(11)
  END FUNCTION R2P_PPISN_BETA_NO_TURN_ON

  PURE FUNCTION R2P_PPISN_BETA_NO_MASS(V) result(p)
  real(kind=prec), intent(in) :: v(:)
  type(para) :: p
  p%mmin = 4.02_prec
  p%dm   = 5.27_prec
  p%mgap = 84.09
  p%a    =  0.12
  p%b    = -1.53
  p%d    = -5.86
  p%lam21= 10**(-2.49)
  p%lam12= 10**(-4.74)
  p%bq0  =  6.78
  p%bq1  =  6.78
  ! Spin
  p%alpha1 = v(1)
  p%beta1 =  v(2)
  p%alpha2 = v(3)
  p%beta2 =  v(4)
  END FUNCTION R2P_PPISN_BETA_NO_MASS

  PURE FUNCTION R2P_PPISN_1BETA_NO_TURN_ON(V) result(p)
  real(kind=prec), intent(in) :: v(:)
  type(para) :: p
  p%mmin = 4.02_prec
  p%dm   = 5.27_prec

  p%mgap =     v(1)
  p%a    =     v(2)
  p%b    =     v(3)
  p%d    =     v(4)
  p%lam21= 10**v(5)
  p%lam12= 10**v(6)
  p%bq0  =     v(7)
  p%bq1  =     v(7)
  ! Spin
  p%alpha1 = v(8)
  p%beta1 =  v(9)
  p%alpha2 = v(8)
  p%beta2 =  v(9)
  END FUNCTION R2P_PPISN_1BETA_NO_TURN_ON

  PURE FUNCTION R2P_PPISN_PLANCK(V) result(p)
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
  p%bq0  = v(9)
  p%bq1  = v(9)

  p%h0   = v(10)
  END FUNCTION R2P_PPISN_PLANCK

                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   !!                             !!
                   !!            TOOLS            !!
                   !!                             !!
                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION GETMODEL(mods) result(m)
  character(len=*), intent(in) :: mods
  type(model) :: m

  m%model_name = trim(mods)
  m%cuts => nocuts
  select case(mods)
    case('plp+flat+trivial+trivial')
      m%ndim = 7
      m%primary => plp_mf
      m%secondary => flatm
      m%redshift => null()
      m%spin(1,1)%f => trivial_spin
      m%spin(1,2)%f => trivial_spin
      m%spin(2,1)%f => trivial_spin
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_plp_flat
      m%smooth => smooth_tanh
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
      m%norms = .false.
    case('plp+plp+trivial+trivial')
      m%ndim = 8
      m%primary => plp_mf
      m%secondary => plp_m2f
      m%redshift => null()
      m%spin(1,1)%f => trivial_spin
      m%spin(1,2)%f => trivial_spin
      m%spin(2,1)%f => trivial_spin
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_plp_pow
      m%smooth => smooth_tanh
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
      m%norms = .false.
    case('plp+pow+trivial+trivial')
      m%ndim = 8
      m%primary => plp_mf
      m%secondary => powm
      m%redshift => null()
      m%spin(1,1)%f => trivial_spin
      m%spin(1,2)%f => trivial_spin
      m%spin(2,1)%f => trivial_spin
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_plp_pow
      m%smooth => smooth_exp
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
      m%norms = .false.

    case('plp+flat+planck+trivial')
      m%ndim = 8
      m%primary => plp_mf
      m%secondary => flatm
      m%redshift => redshift_planck
      m%spin(1,1)%f => trivial_spin
      m%spin(1,2)%f => trivial_spin
      m%spin(2,1)%f => trivial_spin
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_plp_flat_planck
      m%smooth => smooth_tanh
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
      m%norms = .false.
    case('plp+pow+planck+trivial')
      m%ndim = 9
      m%primary => plp_mf
      m%secondary => powm
      m%redshift => redshift_planck
      m%spin(1,1)%f => trivial_spin
      m%spin(1,2)%f => trivial_spin
      m%spin(2,1)%f => trivial_spin
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_plp_pow_planck
      m%smooth => smooth_exp
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
      m%norms = .false.
    case('plp+plp+planck+trivial')
      m%ndim = 9
      m%primary => plp_mf
      m%secondary => plp_m2f
      m%redshift => redshift_planck
      m%spin(1,1)%f => trivial_spin
      m%spin(1,2)%f => trivial_spin
      m%spin(2,1)%f => trivial_spin
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_plp_pow_planck
      m%smooth => smooth_tanh
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
      m%norms = .false.

    case('plp+flat+trivial+beta')
      m%ndim = 9
      m%primary => plp_mf
      m%secondary => flatm
      m%redshift => null()
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
      m%ndim = 10
      m%primary => plp_mf
      m%secondary => powm
      m%redshift => null()
      m%spin(1,1)%f => beta_spin_11
      m%spin(1,2)%f => trivial_spin
      m%spin(2,1)%f => trivial_spin
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_plp_pow_beta
      m%smooth => smooth_exp
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
      m%norms = .false.
    case('plp+plp+trivial+beta')
      m%ndim = 10
      m%primary => plp_mf
      m%secondary => plp_m2f
      m%redshift => null()
      m%spin(1,1)%f => beta_spin_11
      m%spin(1,2)%f => trivial_spin
      m%spin(2,1)%f => trivial_spin
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_plp_pow_beta
      m%smooth => smooth_tanh
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
      m%norms = .false.

    case("ppisn+flat+trivial+trivial")
      m%ndim = 9
      m%primary => ppisn_mf1g
      m%primaryM2 => ppisn_mf2g
      m%secondary => flatm
      m%secondary_c = "flat"
      m%redshift => null()
      m%spin(1,1)%f => trivial_spin
      m%spin(1,2)%f => trivial_spin
      m%spin(2,1)%f => trivial_spin
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_ppisn
      m%smooth => smooth_exp
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
      m%norms = .true.
      m%cuts => bd_lam_cut
    case("ppisn-lam+trivial+trivial")
      m%ndim = 7
      m%primary => ppisn_mf1g
      m%primaryM2 => ppisn_mf2g
      m%secondary => ppisn_m2_phys
      m%secondary_c = "phys"
      m%redshift => null()
      m%spin(1,1)%f => trivial_spin
      m%spin(1,2)%f => trivial_spin
      m%spin(2,1)%f => trivial_spin
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_ppisn_nolam
      m%smooth => smooth_exp
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
      m%norms = .true.
      m%cuts => bd_lam_cut
    case("ppisn+trivial+trivial")
      m%ndim = 9
      m%primary => ppisn_mf1g
      m%primaryM2 => ppisn_mf2g
      m%secondary => ppisn_m2_phys
      m%secondary_c = "phys"
      m%redshift => null()
      m%spin(1,1)%f => trivial_spin
      m%spin(1,2)%f => trivial_spin
      m%spin(2,1)%f => trivial_spin
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_ppisn
      m%smooth => smooth_exp
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
      m%norms = .true.
      m%cuts => bd_lam_cut
    case("ppisn2P+trivial+trivial")
      m%ndim = 10
      m%primary => ppisn_mf1g
      m%primaryM2 => ppisn_mf2g
      m%secondary => ppisn_m2_phys
      m%secondary_c = "phys"
      m%redshift => null()
      m%spin(1,1)%f => trivial_spin
      m%spin(1,2)%f => trivial_spin
      m%spin(2,1)%f => trivial_spin
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_ppisn_twopair
      m%smooth => smooth_exp
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
      m%norms = .true.
      m%cuts => bd_lam_cut
    case("ppisn+planck+trivial")
      m%ndim = 10
      m%primary => ppisn_mf1g
      m%primaryM2 => ppisn_mf2g
      m%secondary => ppisn_m2_phys
      m%secondary_c = "phys"
      m%redshift => redshift_planck
      m%spin(1,1)%f => trivial_spin
      m%spin(1,2)%f => trivial_spin
      m%spin(2,1)%f => trivial_spin
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_ppisn_planck
      m%smooth => smooth_exp
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
      m%norms = .true.
      m%cuts => bd_lam_cut

    case("ppisn+trivial+beta")
      m%ndim = 13
      m%primary => ppisn_mf1g
      m%primaryM2 => ppisn_mf2g
      m%secondary => ppisn_m2_phys
      m%secondary_c = "phys"
      m%redshift => null()
      m%spin(1,1)%f => beta_spin_11
      m%spin(1,2)%f => beta_spin_12
      m%spin(2,1)%f => beta_spin_21
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_ppisn_beta
      m%smooth => smooth_exp
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
      m%norms = .true.
      m%cuts => bd_lam_cut
    case("ppisn+trivial+beta-turnon")
      m%ndim = 11
      m%primary => ppisn_mf1g
      m%primaryM2 => ppisn_mf2g
      m%secondary => ppisn_m2_phys
      m%secondary_c = "phys"
      m%redshift => null()
      m%spin(1,1)%f => beta_spin_11
      m%spin(1,2)%f => beta_spin_12
      m%spin(2,1)%f => beta_spin_21
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_ppisn_beta_no_turn_on
      m%smooth => smooth_exp
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
      m%norms = .true.
      m%cuts => bd_lam_cut
    case("ppisn+trivial+1beta-turnon")
      m%ndim = 9
      m%primary => ppisn_mf1g
      m%primaryM2 => ppisn_mf2g
      m%secondary => ppisn_m2_phys
      m%secondary_c = "phys"
      m%redshift => null()
      m%spin(1,1)%f => beta_spin_11
      m%spin(1,2)%f => beta_spin_12
      m%spin(2,1)%f => beta_spin_21
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_ppisn_1beta_no_turn_on
      m%smooth => smooth_exp
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
      m%norms = .true.
      m%cuts => bd_lam_cut
    case("ppisn+trivial+gauss-turnon")
      m%ndim = 11
      m%primary => ppisn_mf1g
      m%primaryM2 => ppisn_mf2g
      m%secondary => ppisn_m2_phys
      m%secondary_c = "phys"
      m%redshift => null()
      m%spin(1,1)%f => gauss_spin_11
      m%spin(1,2)%f => gauss_spin_12
      m%spin(2,1)%f => gauss_spin_21
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_ppisn_beta_no_turn_on
      ! m%cuts => bd_cut
      m%cuts => bd_lam_cut
      m%smooth => smooth_exp
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
      m%norms = .true.
    case("ppisn+trivial+1gauss-turnon")
      m%ndim = 9
      m%primary => ppisn_mf1g
      m%primaryM2 => ppisn_mf2g
      m%secondary => ppisn_m2_phys
      m%secondary_c = "phys"
      m%redshift => null()
      m%spin(1,1)%f => gauss_spin_11
      m%spin(1,2)%f => gauss_spin_12
      m%spin(2,1)%f => gauss_spin_21
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_ppisn_1beta_no_turn_on
      m%cuts => bd_lam_cut
      ! m%cuts => bd_cut
      m%smooth => smooth_exp
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
      m%norms = .true.

    case("ppisn+trivial+gauss-mass")
      m%ndim = 4
      m%primary => ppisn_mf1g
      m%primaryM2 => ppisn_mf2g
      m%secondary => ppisn_m2_phys
      m%secondary_c = "phys"
      m%redshift => null()
      m%spin(1,1)%f => gauss_spin_11
      m%spin(1,2)%f => gauss_spin_12
      m%spin(2,1)%f => gauss_spin_21
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_ppisn_beta_no_mass
      m%smooth => smooth_exp
      m%smoothint => smooth_expint
      m%smooth_c = "tan"
    case("ppisn+trivial+beta-mass")
      m%ndim = 4
      m%primary => ppisn_mf1g
      m%primaryM2 => ppisn_mf2g
      m%secondary => ppisn_m2_phys
      m%secondary_c = "phys"
      m%redshift => null()
      m%spin(1,1)%f => beta_spin_11
      m%spin(1,2)%f => beta_spin_12
      m%spin(2,1)%f => beta_spin_21
      m%spin(2,2)%f => trivial_spin
      m%r2p => r2p_ppisn_beta_no_mass
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
  real(kind=prec), dimension(25) :: chi1test, chi2test, chians
  real(kind=prec), dimension(6) :: m1dtest, m2dtest, m1stest, m2stest, ztest, dtest
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
  ! print*,plp_int(para(mmax = 50._prec, &
  !                    mum  = 30._prec, &
  !                    sm   = 5._prec, &
  !                    alpha= 2._prec, &
  !                    lp   = 0.2_prec, &
  !                    mmin = 5._prec, &
  !                    dm   = 5._prec, &
  !                    lam12=  0._prec, &
  !                    lam21=  0._prec, &
  !                    sf_c = 'tan'   , &
  !                    sfint= smooth_expint, &
  !                    sf   = smooth_tanh))

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

  mtest = (/ 3., 30., 35., 37., 40., 45. /)
  ans = (/ 0.017562827691927917e-4_prec, 3.4957169791193885e-4_prec, &
          3.4957169791193885e-4_prec, 3.4957169791193885e-4_prec, &
          3.4957169791193885e-4_prec, 1.6828265142804126e-4_prec /)

  diff = sum(abs(ans-ppisn_mf2g(mtest, p))) / 6.
  print*, "ppisn_mf2g", diff !FIXME

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

  p%bq0 = 5 ; p%bq1 = 4
  ans(1:3) = ppisn_nint(p)
  diff = sum(ans(1:3) / (/0.0015955885789676152, 0.26672702783952185, 1.424327104942902/) - 1.) / 3
  print*, 'norm int', diff

#if 0
! I have no idea what this does...
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
#endif

  ! test spin
  p%alpha1 = 1.2
  p%beta1 = 4.2
  p%alpha2 = 3.2
  p%beta2 = 2.2

  chi1test=(/-0.4,-0.4,-0.4,-0.4,-0.4,0.1,0.1,0.1,0.1,0.1,&
    0.4,0.4,0.4,0.4,0.4,0.8,0.8,0.8,0.8,0.8,1.1,1.1,1.1,1.1,1.1/)
  chi2test = (/-0.4,0.1,0.4,0.8,1.1,-0.4,0.1,0.4,0.8,1.1,&
    -0.4,0.1,0.4,0.8,1.1,-0.4,0.1,0.4,0.8,1.1,-0.4,0.1,0.4,0.8,1.1/)

  chians = (/0.,0.,0.,0.,0.,0.,7.954291810492719,2.867619260750767, &
       0.09793534842428382,0.,0.,2.867619260750767,1.033811735921138, &
       0.03530688817316813,0.,0.,0.09793534842428382,0.03530688817316813, &
       0.0012058059597881091,0.,0.,0.,0.,0.,0./)
  print*, 'spin_11', sum(abs(chians - beta_spin_11(chi1test, chi2test, p)))

  chians = (/ 0.,0.,0.,0.,0.,0.,0.2618696892343281,3.3986598645935016, &
       4.178574866102779,0.,0.,0.09440721846093056,1.225258353684312, &
       1.5064272287218412,0.,0.,0.0032242090016225116,0.04184520079782889, &
       0.05144772095095936,0.,0.,0.,0.,0.,0. /)
  print*, 'spin_12', sum(abs(chians - beta_spin_12(chi1test, chi2test, p)))

  chians = (/ 0.,0.,0.,0.,0.,0.,0.2618696892343281,0.09440721846093056, &
       0.0032242090016225116,0.,0.,3.3986598645935016,1.225258353684312, &
       0.04184520079782889,0.,0.,4.178574866102779,1.5064272287218412, &
       0.05144772095095936,0.,0.,0.,0.,0.,0. /)
  print*, 'spin_21', sum(abs(chians - beta_spin_21(chi1test, chi2test, p)))

  ! test redshift
  p%H0 = 70
  p%gamma = 2.5
  m1dtest = (/ 1.,11., 20., 30., 40., 100. /)
  m2dtest = (/ 0.2,10., 15., 3., 33., 80. /)
  dtest = (/ 459.9190751959335, 2165.352153582382, 86.95879617595057, &
        18249.0848891852, 32310.127331891996, 45016.503681422204 /)
  ztest = (/ 0.1, 0.4, 0.02, 2.3, 3.7, 4.9 /)
  ! real(kind=prec), dimension(6) :: m1dtest, m2dtest, m1stest, m2stest, z, d
  call redshift_planck(m1dtest, m2dtest, m1stest, m2stest, ztest, dtest, p, ans)
  print*, 'm1D -> m1S', sum(abs(m1stest / (/ 0.909091, 7.85714, 19.6078, 9.09091, 8.51064, 16.9492 /) - 1)) / 6.
  print*, 'm2D -> m2S', sum(abs(m2stest / (/ 0.181818, 7.14286, 14.7059, 0.909091, 7.02128, 13.5593 /) - 1)) / 6.
  print*, 'prob', sum(abs(ans / (/ 2743.38, 45580., 105.958, 760190., 1.19771e6, 1.48307e6 /)-1))/6

  END SUBROUTINE TEST

                          !!!!!!!!!!!!!!!!!!!!!!
                            END MODULE MODELS
                          !!!!!!!!!!!!!!!!!!!!!!
