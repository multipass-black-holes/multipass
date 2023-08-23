  subroutine pyinterface(mods, what, r, np, dat, ndat, outp)
  use loglikelihood, only: av_likelihood
  use models
  character(len=*), intent(in) :: mods
  character(len=*), intent(in) :: what
  integer, intent(in) :: np, ndat
  real(kind=8), intent(in) :: r(np), dat(ndat)
  real(kind=8), intent(out) :: outp(ndat)

  type(model) :: the_model
  type(para) :: p

  the_model = getmodel(mods)
  p = the_model%r2p(r)
  p%sf => the_model%smooth
  p%sf_c = the_model%smooth_c

  select case(what)
    case('m1')
      outp = the_model%primary(dat, p)
  end select
  end subroutine

  SUBROUTINE LOGLIKE(outp, r, mods, datafile, injfile)
  use loglikelihood, only: ll, load_inj, load_data
  use models
  real(kind=8), intent(in) :: r(:)
  character(len=*), intent(in), optional :: mods
  character(len=*), intent(in), optional :: datafile, injfile
  real(kind=8), intent(out) :: outp
  type(model), save :: the_model
  type(para) :: p

  if(present(datafile) .and. len_trim(datafile) > 0) then
    call load_data(datafile)
  endif
  if(present(injfile) .and. len_trim(injfile) > 0) then
    call load_inj(injfile)
  endif
  if(present(mods) .and. len_trim(mods) > 0) then
    the_model = getmodel(mods)
  endif

  p = the_model%r2p(r)
  p%sf => the_model%smooth
  p%sfint => the_model%smoothint

  outp = ll(the_model, p)
  END SUBROUTINE LOGLIKE
