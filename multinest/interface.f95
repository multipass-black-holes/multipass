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
