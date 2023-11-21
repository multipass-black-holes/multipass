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

  integer            :: np        = 1000   ! number of live points
  integer, parameter :: ndimmax   = 20
  integer            :: npara              ! no. of parameters
  integer            :: nparaMode          ! no. of parameters with mode
  integer, parameter :: maxmode   = 100    ! max modes
  integer, parameter :: feedbackn = 10     ! require feedback and update after no. it
  integer, parameter :: seed      = 1234   ! seed
  integer, parameter :: maxiter   = 0      ! maxiter
  integer, parameter :: pWrap(ndimMax) = 0
  real(kind=prec)            :: tol     = 5.5   ! tol, defines stopping
  real(kind=prec)            :: efr     = 0.8   ! efr, requried effiecny
  real(kind=prec), parameter :: Ztol    = -1e30 ! all the modes with logZ < Ztol are ignored
  real(kind=prec), parameter :: logZero = -1e30 ! logZero

  real(kind=prec) :: start, finish

  integer :: context

  real(kind=prec), dimension(ndimMax) :: min_val, max_val
  character(len=1000) :: arg, root = "testrun/test-"

  type(model) :: the_model

  call get_command_argument(1, arg)
  if(len_trim(arg) == 0) then
    call userinterface(5, 9)
  else
    open(unit=15, file=trim(arg), action="read")
    call userinterface(15, 9)
    close(unit=15)
  endif

  call load_inj("inj.rec")
  call load_data("data.rec")

  call cpu_time(start)
  call nestrun(IS, modal, ceff, np, tol, efr, &
      the_model%ndim, npara, nparaMode, maxmode, feedbackn, &
      Ztol, root, seed, pWrap(1:the_model%ndim), fstdout, resume, &
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

  cube = min_val(1:ndim) + cube * (max_val(1:ndim) - min_val(1:ndim))
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

  call cpu_time(finish)
  print*,     "nSamples = ",nSamples, "nlive = ", nlive, "logZ = ", logZ, "+-", logZerr, " time = ", finish - start, "s"
  write(9, *) "nSamples = ",nSamples, "nlive = ", nlive, "logZ = ", logZ, "+-", logZerr, " time = ", finish - start, "s"
  call flush(9)

  END SUBROUTINE DUMPER


  SUBROUTINE USERINTERFACE(Uin, Uout)
  implicit none
  integer, intent(in) :: Uin, Uout
  character(len=1) :: cmd, space
  character(len=1000) :: buf

  do while(cmd .ne. 's')
    if(Uin==5) then
      print*, "Available commands:"
      print*, " * [m]odel: ", trim(the_model%model_name)
      ! generate using
      ! cat models.f95| perl -ne 'print if s/ *case\(.(.*\+.*).\)/      print*, "    - \1"/g'
      print*, "    - plp+flat+trivial+trivial"
      print*, "    - plp+pow+trivial+trivial"
      print*, "    - plp+flat+planck+trivial"
      print*, "    - plp+pow+planck+trivial"
      print*, "    - plp+flat+trivial+beta"
      print*, "    - plp+pow+trivial+beta"
      print*, "    - plp+plp+trivial+beta"
      print*, "    - ppisn+flat+trivial+trivial"
      print*, "    - ppisn+trivial+trivial"
      print*, "    - ppisn+planck+trivial"
      print*, "    - ppisn+trivial+beta"
      print*, "    - ppisn+trivial+beta-turnon"
      print*, " * [n]umber of live points: ", np
      print*, " * [t]olerance (defines stopping)", tol
      print*, " * [e]fr, require efficency", efr
      print*, " * [p]rior range, lower bounds", min_val(1:the_model%ndim)
      print*, " * [P]rior range, upper bounds", max_val(1:the_model%ndim)
      print*, " * [o]utput prefix ", trim(root)
      print*, " * [s]tart"
    endif

    read(unit=Uin, fmt='(A1, A1, A)') cmd, space, buf
    if(space.ne. ' ') cycle

    select case(cmd)
      case('m')
        the_model = getmodel(trim(buf))
        npara     = the_model%ndim
        nparaMode = the_model%ndim

      case('p')
        read(unit=buf, fmt=*) min_val(1:the_model%ndim)
      case('P')
        read(unit=buf, fmt=*) max_val(1:the_model%ndim)

      case('n')
        read(unit=buf, fmt=*) np
      case('t')
        read(unit=buf, fmt=*) tol
      case('e')
        read(unit=buf, fmt=*) efr
      case('o')
        root = trim(buf)
    end select
  enddo

  call system("date > " // trim(root) // ".timing")
  call system("git diff >> " // trim(root) // ".timing")

  open(unit=Uout, file=trim(root) // ".timing", action="write", access="append")

  write(Uout, *) " * model: ", trim(the_model%model_name)
  write(Uout, *) " * number of live points: ", np
  write(Uout, *) " * tolerance (defines stopping)", tol
  write(Uout, *) " * efr, require efficency", efr
  write(Uout, *) " * prior range, lower bounds", min_val(1:the_model%ndim)
  write(Uout, *) " * Prior range, upper bounds", max_val(1:the_model%ndim)
  write(Uout, *) " * output prefix ", trim(root)

  flush(Uout)

  END SUBROUTINE USERINTERFACE

  END PROGRAM MAIN
