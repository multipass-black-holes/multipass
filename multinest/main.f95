  PROGRAM MAIN
  use loglikelihood
  use nested, only: nestrun
  implicit none
  logical, parameter :: IS        = .true. ! do nested importance sampling?
  logical, parameter :: modal     = .false.! do mode separation?
  logical, parameter :: ceff      = .false.! run in constant efficiency mode?
  logical, parameter :: fstdout   = .true. ! feedback to stdout
  logical, parameter :: resume    = .false.! resume?
  logical, parameter :: outfile   = .true. ! write output files
  logical, parameter :: MPIinit   = .true. ! do MPI init

  integer, parameter :: np        = 1000   ! number of live points
  integer, parameter :: ndim      = 8      ! ndim = 8
  integer, parameter :: npara     = ndim   ! no. of parameters
  integer, parameter :: nparaMode = ndim   ! no. of parameters with mode
  integer, parameter :: maxmode   = 100    ! max modes
  integer, parameter :: feedbackn = 10     ! require feedback and update after no. it
  integer, parameter :: seed      = 1234   ! seed
  integer, parameter :: maxiter   = 0      ! maxiter
  integer, parameter :: pWrap(ndim) = 0
  real(kind=prec), parameter :: tol     = 5.5   ! tol, defines stopping
  real(kind=prec), parameter :: efr     = 0.8   ! efr, requried effiecny
  real(kind=prec), parameter :: Ztol    = -1e30 ! all the modes with logZ < Ztol are ignored
  real(kind=prec), parameter :: logZero = -1e30 ! logZero

  integer :: context

  real(kind=prec), parameter, dimension(npara) :: &
    min_val = (/30._prec, 20._prec, 1._prec, -4._prec, 0._prec, 2._prec, 0._prec, 0._prec /)
  real(kind=prec), parameter, dimension(npara) :: &
    max_val = (/100._prec, 50._prec, 10._prec, 12._prec, 1._prec, 10._prec, 10._prec, 10._prec/)
  character(len=1000), parameter :: root = "testrun/"

  call load_inj("inj.rec")
  call load_data("data.rec")

  call nestrun(IS, modal, ceff, np, tol, efr, &
      ndim, npara, nparaMode, maxmode, feedbackn, &
      Ztol, root, seed, pWrap, fstdout, resume, &
      outfile, MPIinit, logZero, maxiter, &
      slikelihood, dumper, context)
contains

  SUBROUTINE SLIKELIHOOD(CUBE, NDIM, NPAR, LNEW, CONTEXT)
  integer :: ndim, npar, context
  real(kind=prec) :: cube(ndim)
  real(kind=prec) :: lnew
  real(kind=prec) :: para(npara)
  para = min_val + cube * (max_val - min_val)
  lnew = ll(real2model(mmax = para(1), &
                       mum  = para(2), &
                       sm   = para(3), &
                       alpha= para(4), &
                       lp   = para(5), &
                       mmin = para(6), &
                       dm   = para(7), &
                       k    = para(8)))
  END SUBROUTINE SLIKELIHOOD

  SUBROUTINE DUMPER(nSamples, nlive, nPar, physLive, posterior, &
                    paramConstr, maxLogLike, logZ, INSlogZ, &
                    logZerr, context)
  integer :: nSamples, nLive, nPar, context
  real(kind=prec), pointer, dimension(:,:) :: physLive, posterior
  real(kind=prec), pointer, dimension(:) :: paramConstr
  real(kind=prec) :: maxLogLike, logZ, INSLogZ, logZerr

  print*,nSamples

  END SUBROUTINE DUMPER

  END PROGRAM MAIN
