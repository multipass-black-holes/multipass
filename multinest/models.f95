                          !!!!!!!!!!!!!!!!!!!!!!
                              MODULE MODELS
                          !!!!!!!!!!!!!!!!!!!!!!
  use functions
  implicit none

  ABSTRACT INTERFACE
    PURE FUNCTION SMOOTH_FUNC(M, MI, DM)
      use functions, only: prec
      real(kind=prec), intent(in) :: m(:), mi, dm
      real(kind=prec) :: smooth_func(size(m))
    END FUNCTION SMOOTH_FUNC
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
    procedure(smooth_func), pointer, nopass :: sf
  END TYPE PARA


contains

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
