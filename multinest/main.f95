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
  integer, parameter :: ndim      = 9      ! ndim = 9
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

  real(kind=prec) :: start, finish

  integer :: context

  !real(kind=prec), parameter, dimension(npara) :: &
  !  min_val = (/ 2._prec,  0._prec,  30._prec, 20._prec,  1._prec, -4._prec, 0._prec,  0._prec /)
  !real(kind=prec), parameter, dimension(npara) :: &
  !  max_val = (/10._prec, 10._prec, 100._prec, 50._prec, 10._prec, 12._prec, 1._prec, 10._prec/)
  real(kind=prec), parameter, dimension(npara) :: &
    min_val = (/ 2._prec,  0._prec,  20._prec,  0._prec,- 4._prec,-10._prec, -7._prec,  -7._prec, -7._prec /)
  real(kind=prec), parameter, dimension(npara) :: &
    max_val = (/10._prec, 10._prec, 120._prec, 0.5_prec,  0._prec,  0._prec, -0.3_prec, -0.3_prec, -0.3_prec /)
  character(len=1000), parameter :: root = "testrun/test-"

  type(model) :: the_model

  call system("date > " // trim(root) // ".timing")
  call system("git diff >> " // trim(root) // ".timing")

  open(unit=9, file=trim(root) // ".timing", action="write", access="append")


  call load_inj("inj.rec")
  call load_data("data.rec")

  !the_model = getmodel('plp+pow+trivial+trivial')
  the_model = getmodel('ppisn+flat+trivial+trivial')

  call cpu_time(start)
  call nestrun(IS, modal, ceff, np, tol, efr, &
      ndim, npara, nparaMode, maxmode, feedbackn, &
      Ztol, root, seed, pWrap, fstdout, resume, &
      outfile, MPIinit, logZero, maxiter, &
      slikelihood, dumper, context)

  call cpu_time(finish)
  write(9, *) "Complete, took", finish-start, "s"
  close(unit=9)


contains

  SUBROUTINE SLIKELIHOOD(CUBE, NDIM, NPAR, LNEW, CONTEXT)
  integer :: ndim, npar, context
  real(kind=prec) :: cube(ndim)
  real(kind=prec) :: lnew
  type(para) :: p

  cube = min_val + cube * (max_val - min_val)
  p = the_model%r2p(cube)
  p%sf => the_model%smooth
  p%sfint => the_model%smoothint
  p%sf_c = the_model%smooth_c

  lnew = ll(the_model, p)
  END SUBROUTINE SLIKELIHOOD

  SUBROUTINE DUMPER(nSamples, nlive, nPar, physLive, posterior, &
                    paramConstr, maxLogLike, logZ, INSlogZ, &
                    logZerr, context)
  integer :: nSamples, nLive, nPar, context
  real(kind=prec), pointer, dimension(:,:) :: physLive, posterior
  real(kind=prec), pointer, dimension(:) :: paramConstr
  real(kind=prec) :: maxLogLike, logZ, INSLogZ, logZerr

  print*, "nSamples = ",nSamples, " time = ", finish - start, "s"
  call cpu_time(finish)
  write(9, *) "nSamples = ",nSamples, " time = ", finish - start, "s"
  call flush(9)

  END SUBROUTINE DUMPER

  END PROGRAM MAIN
