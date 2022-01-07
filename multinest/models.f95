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
              + lp * gauss(mBH, mum, sm) ) &
            * smooth(mBH, mmin, dm)

  END FUNCTION PLP_MF


                          !!!!!!!!!!!!!!!!!!!!!!
                            END MODULE MODELS
                          !!!!!!!!!!!!!!!!!!!!!!
