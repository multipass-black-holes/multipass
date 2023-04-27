                          !!!!!!!!!!!!!!!!!!!!!!
                           MODULE LOGLIKELIHOOD
                          !!!!!!!!!!!!!!!!!!!!!!
  use functions
  use models, only: plp_mf, model, para, getmodel, ppisn_norms
  implicit none
  real(kind=prec), allocatable, dimension(:,:) :: injections, dat
  integer, allocatable, dimension(:) :: offsets

contains

  PURE FUNCTION MEAN(dat)
  real(kind=prec), intent(in) :: dat(:)
  real(kind=prec) :: mean
  mean = sum(dat) / size(dat)
  END FUNCTION MEAN

  ! This builds the likelihood array for an event but *does not*
  ! average it
  PURE FUNCTION AV_LIKELIHOOD(dat, m, p)
  real(kind=prec), intent(in) :: dat(:, :)
  type(model), intent(in) :: m
  type(para), intent(in) :: p
  real(kind=prec) :: av_likelihood(size(dat,1))
  integer i

  av_likelihood = m%primary(dat(:, 1), p)  &
                * m%secondary(dat(:,1), dat(:,2), p) &
                * m%redshift(dat(:,3), p) &
                * m%spin(1,1)%f(dat(:,4), p)

  if(m%norms) then
    av_likelihood = av_likelihood + &
        ppisn_norms(dat(:,1), dat(:,2), dat(:,3), dat(:,4), m, p)
  endif

  END FUNCTION AV_LIKELIHOOD


  PURE FUNCTION LL(M, P)
  type(model), intent(in) :: M
  type(para), intent(in) :: p
  real(kind=prec) :: ll
  real(kind=prec) :: avg(size(dat,1)), inj(size(injections,1))
  real(kind=prec) :: tmp, acc, Nlxi
  integer i

  avg = av_likelihood(dat, m, p)

  inj = av_likelihood(injections, m, p)
  inj = inj / injections(:,1)**(-4.35)
  inj = inj / injections(:,2)**2

  ! We need to average the avg for each event file as delimited by offsets
  acc = 0.
  do i=1,size(offsets)-1
    tmp = mean(avg(offsets(i):offsets(i+1)-1))
    if(tmp <= 0) then
      tmp = -1e6
    else
      tmp = log(tmp)
    endif
    acc = acc + tmp
  enddo

  Nlxi = (size(offsets) - 1) * log(mean(inj))

  ll = -Nlxi + acc
  END FUNCTION LL


  SUBROUTINE LOAD_INJ(FN)
  character(len=*), intent(in) :: fn
  integer :: n

  open(unit=8, action='read', form='unformatted', file=trim(fn))

  read(8) n
  allocate(injections(n,4))
  read(8) injections

  close(unit=8)
  END SUBROUTINE LOAD_INJ

  SUBROUTINE LOAD_DATA(FN)
  character(len=*), intent(in) :: fn
  integer :: nO, nD

  open(unit=8, action='read', form='unformatted', file=trim(fn))

  read(8) nD, nO
  allocate(dat(nD,4))
  allocate(offsets(nO))
  read(8) offsets
  read(8) dat

  offsets = offsets+1

  close(unit=8)
  END SUBROUTINE LOAD_DATA

  SUBROUTINE TEST
  real(kind=prec) :: ans
  type(model) :: m
  type(para) :: p
  type(model) :: the_model
  real(kind=prec), parameter, dimension(8) :: &
    min_val = (/ 2._prec,  0._prec,  30._prec, 20._prec,  1._prec, -4._prec, 0._prec,  0._prec /)
  real(kind=prec), parameter, dimension(8) :: &
    max_val = (/10._prec, 10._prec, 100._prec, 50._prec, 10._prec, 12._prec, 1._prec, 10._prec/)
  real(kind=prec) :: cube(8)

  call load_inj("inj.rec")
  call load_data("data.rec")

  p = para(mmax = 39.3856_prec,&
           mum  = 44.7844_prec,&
           sm   = 6.06597_prec,&
           alpha= -1.71307_prec,&
           lp   = 0.394535_prec,&
           mmin = 2.61592_prec,&
           dm   = 8.54510_prec,&
           k    = 4.40876_prec,&
           sf_c = 'tan',&
           lam12=  0._prec, &
           lam21=  0._prec, &
           sfint= smooth_expint, &
           sf   = smooth_tanh)

  m = getmodel('plp+pow+trivial+trivial')

  ans = mean(av_likelihood(dat(offsets(47):offsets(48)-1,:), m, p))
  print*, ans / 0.0048440957269826265 - 1

  ans = mean(av_likelihood(injections, m, p))
  print*, ans / 0.005397037082072583 - 1

  ans = ll(m, p)
  print*,ans / -3000396.124014442 - 1

  cube = (/0.6720379591,0.0260769129,0.3415273428,0.1477759480,0.4869252443,0.7336550355,0.9522889853,0.4295404553/)
  the_model = getmodel('plp+pow+trivial+trivial')
  p = the_model%r2p(min_val + cube((/6, 7, 1, 2, 3, 4, 5, 8/)) * (max_val - min_val))
  p%sf => the_model%smooth
  p%sf_c = the_model%smooth_c

  ans = ll(the_model, p) / -1000421.514
  print*,'ans',ans

  cube = (/0.6293965578,0.4330928922,0.3437212706,0.4545528889,0.8079833388,0.0935372710,0.2308301330,0.9758238792/)
  p = the_model%r2p(min_val + cube((/6, 7, 1, 2, 3, 4, 5, 8/)) * (max_val - min_val))
  p%sf => the_model%smooth
  p%sf_c = the_model%smooth_c

  ans = ll(the_model, p) / -449.972
  print*,'ans',ans

  the_model = getmodel('ppisn+trivial+trivial')
  p = the_model%r2p((/2.53325764e+00_prec,& ! mmax
                      2.21172743e-02_prec,& ! dm
                      4.76023054e+01_prec,& ! mgap
                      4.35632844e-01_prec,& ! a
                      3.22931846e-02_prec,& ! b
                      4.09565979e+00_prec,& ! d
                     -5.39496671e+00_prec,&
                     -5.33128390e-01_prec,&
                     -5.02651228e+00_prec/))
  p%sf => the_model%smooth
  p%sfint => the_model%smoothint
  ans = ll(the_model, p)
  print*, 'x',ans

  END SUBROUTINE TEST

                          !!!!!!!!!!!!!!!!!!!!!!
                         END MODULE LOGLIKELIHOOD
                          !!!!!!!!!!!!!!!!!!!!!!
