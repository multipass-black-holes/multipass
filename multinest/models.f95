                          !!!!!!!!!!!!!!!!!!!!!!
                              MODULE MODELS
                          !!!!!!!!!!!!!!!!!!!!!!
  use functions
  implicit none

contains

  ! This implements the power-law + peak model (PLP) for the primary mass of the
  ! 1st generation. The model is specified in B2 of [2010.14533]
  PURE FUNCTION PLP_MF(MBH, MMAX, MUM, SM, ALPHA, LP, MMIN, DM)
  real(kind=prec), intent(in) :: mBH(:)
  real(kind=prec), intent(in) :: mmax, mum, sm, alpha, lp, mmin, dm
  real(kind=prec) :: plp_mf(size(mBH))

  plp_mf = 0.
  where ((mmin < mBH) .and. (mBH < mmax)) &
    plp_mf = ( (1-lp) * powerlaw(mbh, -alpha, mmin, mmax) &
              +   lp  * gauss(mBH, mum, sm) ) &
            * smooth_tanh(mBH, mmin, dm)

  END FUNCTION PLP_MF


  SUBROUTINE TEST
  real(kind=prec), dimension(6) :: mtest, ans
  real(kind=prec) :: diff

  mtest = (/ 1.,11., 20., 30., 40., 100. /)
  ans = (/ 0., 0.03674262, 0.01327075, 0.02089596, 0.00493742, 0. /)

  diff = sum(abs(ans - plp_mf(mtest, &
                50._prec, 30._prec, 5._prec, &
                2._prec, 0.2_prec, 5._prec, 5._prec))) / 6.
  print*,diff

  mtest = (/ 2., 3., 4., 5., 10., 20. /)
  ans = (/ 0., 3.14309363e-13, 5.66146660e-06, &
           1.21599290e-04, 3.97884590e-03, 1.30775022e-02 /)

  diff = sum(abs(ans - plp_mf(mtest, &
                39.3856_prec,44.7844_prec,6.06597_prec,-1.71307_prec,&
                0.394535_prec,2.61592_prec,8.5451_prec ))) / 6.

  print*,diff
  END SUBROUTINE TEST

                          !!!!!!!!!!!!!!!!!!!!!!
                            END MODULE MODELS
                          !!!!!!!!!!!!!!!!!!!!!!
