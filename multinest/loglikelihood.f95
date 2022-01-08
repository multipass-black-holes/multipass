                          !!!!!!!!!!!!!!!!!!!!!!
                           MODULE LOGLIKELIHOOD
                          !!!!!!!!!!!!!!!!!!!!!!
  use models, only: plp_mf, prec
  implicit none
  real(kind=prec), allocatable :: injections(:,:)

contains

  SUBROUTINE LOAD_INJ(FN)
  character(len=*), intent(in) :: fn
  integer :: n

  open(unit=8, action='read', form='unformatted', file=trim(fn))

  read(8) n
  print*,n
  allocate(injections(n,4))
  read(8) injections

  close(unit=8)
  END SUBROUTINE LOAD_INJ

  SUBROUTINE TEST
  END SUBROUTINE TEST

                          !!!!!!!!!!!!!!!!!!!!!!
                         END MODULE LOGLIKELIHOOD
                          !!!!!!!!!!!!!!!!!!!!!!
