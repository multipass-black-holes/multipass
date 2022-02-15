                          !!!!!!!!!!!!!!!!!!!!!!
                              MODULE MODELS
                          !!!!!!!!!!!!!!!!!!!!!!
  use functions
  implicit none

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

    ! Smooth function
    procedure(smoothfn), pointer, nopass :: sf
    character(len=3) :: sf_c

    real(kind=prec) :: lam21, lam12, lam22
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

    PURE FUNCTION SPINFN(chieff, P)
      use functions, only: prec
      import para
      real(kind=prec), intent(in) :: chieff(:)
      type(para), intent(in) :: p
      real(kind=prec) :: spinfn(size(chieff))
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
    procedure(mass1fn), pointer, nopass :: primary
    procedure(mass1fn), pointer, nopass :: primaryM2
    procedure(mass2fn), pointer, nopass :: secondary
    procedure(redshiftFn), pointer, nopass :: redshift
    type(spinFns), dimension(2,2) :: spin
    procedure(paraFn), pointer, nopass :: r2p
    procedure(smoothFn), pointer, nopass :: smooth
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
  real(kind=prec) :: powm(size(m1)), stargx2(size(m1))

  powm = (m2 / m1) ** p%k
  !powm = powm * p%sf(m2, p%mmin, p%dm)

  stargx2 = 2 * (m2 - p%mmin) / p%dm - 1
  where ((p%mmin < m2).and.(m2 < p%mmin+p%dm)) &
    powm = powm * 0.5 * (1+tanh(2*stargx2 / (1-stargx2*stargx2)))

  END FUNCTION POWM

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
  real(kind=prec), parameter :: mpiv   = 30
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
  real(kind=prec), parameter :: mpiv   = 30
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
  PURE FUNCTION PPISN_PM2M1DEN_M11G(M, P)
  real(kind=prec), intent(in) :: m(:)
  type(para), intent(in) :: p
  real(kind=prec) :: ppisn_pm2m1den_m11g(size(m))

  integer, parameter :: Nsamples = 10000
  real(kind=prec) :: BT(0:size(m))
  real(kind=prec) :: sLVC(2)
  integer i, cut
  real(kind=prec), parameter :: linspace(0:Nsamples) = [(i/Nsamples, i=0,Nsamples)]
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
  PURE FUNCTION PPISN_PM2M1DEN_M12G(M, P)
  real(kind=prec), intent(in) :: m(:)
  type(para), intent(in) :: p
  real(kind=prec) :: ppisn_pm2m1den_m12g(size(m))

  ppisn_pm2m1den_m12g = 1._prec

  where( (p%mmin < m).and.(m < p%mmin+p%dm) ) &
    ppisn_pm2m1den_m12g = 0.5 * (m - p%mmin)  ! TODO
  where( (p%mmin+p%dm < m).and.(m < p%mmax) ) &
    ppisn_pm2m1den_m12g = (0.5*p%dm) + m
  where( (p%mmax < m) ) &
    ppisn_pm2m1den_m12g = (0.5*p%dm) + p%mmax + &
        ((m/p%mmax)**(1+p%d) - 1)*p%mmax/(1+p%d)

  END FUNCTION PPISN_PM2M1DEN_M12G


  PURE FUNCTION PPISN_M2_PHYS(M1, M2, P)
  real(kind=prec), intent(in) :: m1(:), m2(:)
  type(para), intent(in) :: p
  real(kind=prec) :: ppisn_m2_phys(size(m1))

  ppisn_m2_phys = ppisn_mf1g(m2, p) / ppisn_pm2m1den_m11g(m1, p)

  END FUNCTION PPISN_M2_PHYS


  ! This calculates the correction term (1.9) with
  ! norms = (/lam21, lam12, lam22)
  PURE FUNCTION PPISN_NORMS(M1, M2, CHI, Z, M, P)
  real(kind=prec), intent(in) :: m1(:), m2(:), chi(:), z(:)
  type(model), intent(in) :: m
  type(para), intent(in) :: p
  real(kind=prec) :: ppisn_norms(size(m1))

  real(kind=prec), dimension(size(m1)) :: L21, L12, L22
  real(kind=prec) :: lam21, lam12, lam22

  L21 = p%lam21 * m%primaryM2(m1,p)*m%spin(2,1)%f(chi,p)*m%redshift(z,p)
  L12 = p%lam12 * m%primary  (m1,p)*m%spin(1,2)%f(chi,p)*m%redshift(z,p)
  L22 = p%lam22 * m%primaryM2(m1,p)*m%spin(2,2)%f(chi,p)*m%redshift(z,p)

  select case(m%secondary_c)
    case('phys')
      L21 = L21 * m%primary  (m2,p) / ppisn_pm2m1den_m12g(m1, p)
      L12 = L12 * m%primaryM2(m2,p) / ppisn_pm2m1den_m11g(m1, p)
      L22 = L22 * m%primaryM2(m2,p) / ppisn_pm2m1den_m12g(m1, p)
    case('flat')
      L21 = L21 / (m1 - p%mmin)
      L12 = L12 / (m1 - p%mmin)
      L22 = L22 / (m1 - p%mmin)
  end select

  ppisn_norms = L21+L12+L22

  END FUNCTION PPISN_NORMS

                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   !!                             !!
                   !!            TOOLS            !!
                   !!                             !!
                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION GETMODEL(mods) result(m)
  character(len=*), intent(in) :: mods
  type(model) :: m

  select case(mods)
    case('plp+flat+trivial+trivial')
      m%primary => plp_mf
      m%secondary => flatm
      m%redshift => trivial
      m%spin(1,1)%f => trivial
      m%spin(1,2)%f => trivial
      m%spin(2,1)%f => trivial
      m%spin(2,2)%f => trivial
      m%r2p => r2p_plp_flat
      m%smooth => smooth_tanh
      m%smooth_c = "tan"
      m%norms = .false.
    case('plp+pow+trivial+trivial')
      m%primary => plp_mf
      m%secondary => powm
      m%redshift => trivial
      m%spin(1,1)%f => trivial
      m%spin(1,2)%f => trivial
      m%spin(2,1)%f => trivial
      m%spin(2,2)%f => trivial
      m%r2p => r2p_plp_pow
      m%smooth => smooth_exp
      m%smooth_c = "tan"
      m%norms = .false.
    case default
      stop 9
  end select
  END FUNCTION GETMODEL

  SUBROUTINE TEST
  real(kind=prec), dimension(6) :: mtest, ans
  real(kind=prec) :: diff
  type(para) :: p

  mtest = (/ 1.,11., 20., 30., 40., 100. /)
  ans = (/0._prec,  0.010459028944395847_prec, &
        0.0035206207261596584_prec,  0.0017915249285508827_prec, &
        0.0012500000000000002_prec,  0._prec /)

  print*,sum(abs(ans - ppisn_mf1g(mtest, &
                para(mgap = 50._prec, &
                     a    =  .5_prec, &
                     b    = -2._prec, &
                     d    = -3._prec, &
                     mmin =  5._prec, &
                     dm   =  5._prec, &
                     lam12=  0._prec, &
                     lam21=  0._prec, &
                     lam22=  0._prec, &
                     sf_c = 'exp'   , &
                     sf   = smooth_exp)))) / 6.

  ans = (/ 0._prec, 1/900._prec, 1/900._prec, &
           1/900._prec, 1/900._prec, 0.00021123263888888892_prec/)
  print*,sum(abs(ans - ppisn_mf2g(mtest, &
                para(mgap = 50._prec, &
                     a    =  .5_prec, &
                     b    = -2._prec, &
                     d    = -3._prec, &
                     mmin =  5._prec, &
                     dm   =  5._prec, &
                     lam12=  0._prec, &
                     lam21=  0._prec, &
                     lam22=  0._prec, &
                     sf_c = 'exp'   , &
                     sf   = smooth_exp)))) / 6.

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
                     lam22=  0._prec, &
                     sf_c = 'tan'   , &
                     sf   = smooth_tanh)))) / 6.
  print*,diff

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
                     lam22=  0._prec, &
                     sf_c = 'tan'   , &
                     sf   = smooth_tanh)))) / 6.

  print*,diff

  mtest = (/ 2., 4., 10., 20., 25., 35. /)
  p = para(mmin = 2.61592_prec,&
           dm   = 8.5451_prec,&
           mgap = 34._prec, &
           a    = 0.23_prec, &
           b    = -2.34_prec, &
           d    = -7.63_prec, &
           lam21 = 0._prec, &
           lam22 = 0._prec, &
           lam12 = 0._prec, &
           sf   = smooth_erf, &
           sf_c = 'erf')

  ans = (/0._prec, 0.0011329932907377792_prec, &
          0.004816986990970342_prec, 0.001047864239495318_prec,&
          0.0006707850448342742_prec, 0._prec /)
  print*, sum(abs(ans-ppisn_mf1g(mtest, p)))

  ans = (/ 0._prec, 0.0008369386215617741_prec, &
        0.028895302280797192_prec, 0.05158740716028034_prec, &
        0.0557564000458822_prec, 0.0_prec /)
  print*, sum((ppisn_pm2m1den_m11g(mtest, p) - ans) / (ans + 1e-15))

  mtest = (/ 3., 30., 35., 37., 40., 45. /)
  ans = (/ 0.017562827691927917e-4_prec, 3.4957169791193885e-4_prec, &
          3.4957169791193885e-4_prec, 3.4957169791193885e-4_prec, &
          3.4957169791193885e-4_prec, 1.6828265142804126e-4_prec /)

  print*, sum(abs(ans-ppisn_mf2g(mtest, p)))
  END SUBROUTINE TEST

                          !!!!!!!!!!!!!!!!!!!!!!
                            END MODULE MODELS
                          !!!!!!!!!!!!!!!!!!!!!!
