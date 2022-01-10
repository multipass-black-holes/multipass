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


  TYPE MODEL
    procedure(mass1fn), pointer, nopass :: primary
    procedure(mass2fn), pointer, nopass :: secondary
    procedure(redshiftFn), pointer, nopass :: redshift
    procedure(spinFn), pointer, nopass :: spin
    procedure(paraFn), pointer, nopass :: r2p
    procedure(smoothFn), pointer, nopass :: smooth

  END TYPE MODEL

contains

  PURE FUNCTION TRIVIAL(M, P)
  real(kind=prec), intent(in) :: m(:)
  type(para), intent(in) :: p
  real(kind=prec) :: trivial(size(m))
  trivial = 1.
  END FUNCTION TRIVIAL

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
  powm = powm * p%sf(m2, p%mmin, p%dm)
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

  FUNCTION GETMODEL(mods) result(m)
  character(len=*), intent(in) :: mods
  type(model) :: m

  select case(mods)
    case('plp+flat+trivial+trivial')
      m%primary => plp_mf
      m%secondary => flatm
      m%redshift => trivial
      m%spin => trivial
      m%r2p => r2p_plp_flat
      m%smooth => smooth_tanh
    case('plp+pow+trivial+trivial')
      m%primary => plp_mf
      m%secondary => powm
      m%redshift => trivial
      m%spin => trivial
      m%r2p => r2p_plp_pow
      m%smooth => smooth_tanh
    case default
      stop 9
  end select
  END FUNCTION GETMODEL

  SUBROUTINE TEST
  real(kind=prec), dimension(6) :: mtest, ans
  real(kind=prec) :: diff

  mtest = (/ 1.,11., 20., 30., 40., 100. /)
  ans = (/ 0., 0.03674262, 0.01327075, 0.02089596, 0.00493742, 0. /)

  diff = sum(abs(ans - plp_mf(mtest, &
                para(mmax = 50._prec, &
                     mum  = 30._prec, &
                     sm   = 5._prec, &
                     alpha= 2._prec, &
                     lp   = 0.2_prec, &
                     mmin = 5._prec, &
                     dm   = 5._prec, &
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
                     sf   = smooth_tanh)))) / 6.

  print*,diff
  END SUBROUTINE TEST

                          !!!!!!!!!!!!!!!!!!!!!!
                            END MODULE MODELS
                          !!!!!!!!!!!!!!!!!!!!!!
