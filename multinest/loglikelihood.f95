                          !!!!!!!!!!!!!!!!!!!!!!
                           MODULE LOGLIKELIHOOD
                          !!!!!!!!!!!!!!!!!!!!!!
  use functions
  use models, only: plp_mf
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
  PURE FUNCTION AV_LIKELIHOOD(dat, args)
  real(kind=prec), intent(in) :: dat(:, :), args(:)
  real(kind=prec) :: av_likelihood(size(dat,1))
  integer i

  av_likelihood = plp_mf(dat(:, 1), &
    args(1), args(2), args(3), args(4), &
    args(5), args(6), args(7))

  av_likelihood = av_likelihood * ( dat(:, 2) / dat(:,1) ) ** args(8)
  av_likelihood = av_likelihood * smooth_tanh(dat(:,2), args(6), args(7))

  END FUNCTION AV_LIKELIHOOD


  PURE FUNCTION LL(ARGS)
  real(kind=prec), intent(in) :: args(:)
  real(kind=prec) :: ll
  real(kind=prec) :: avg(size(dat,1)), inj(size(injections,1))
  real(kind=prec) :: tmp, acc, Nlxi
  integer i

  avg = av_likelihood(dat, args)

  inj = av_likelihood(injections, args)
  inj = inj / injections(:,1)**(-4.35)
  inj = inj / injections(:,2)**2

  ! We need to average the avg for each event file as delimited by offsets
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
  real(kind=prec) :: para(8), ans
  call load_inj("inj.rec")
  call load_data("data.rec")

  para = (/ 39.3856_prec, 44.7844_prec, 6.06597_prec, -1.71307_prec, &
           0.394535_prec, 2.61592_prec, 8.54510_prec,  4.40876_prec /)

  ans = mean(av_likelihood(dat(offsets(37):offsets(38)-1,:), para))
  print*, ans / 0.0048440957269826265 - 1

  ans = mean(av_likelihood(injections, para))
  print*, ans / 0.005397037082072583 - 1

  ans = ll(para)
  print*,ans / -3000396.124014442 - 1

  END SUBROUTINE TEST

                          !!!!!!!!!!!!!!!!!!!!!!
                         END MODULE LOGLIKELIHOOD
                          !!!!!!!!!!!!!!!!!!!!!!
