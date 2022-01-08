                          !!!!!!!!!!!!!!!!!!!!!!
                           MODULE LOGLIKELIHOOD
                          !!!!!!!!!!!!!!!!!!!!!!
  use models, only: plp_mf, prec
  implicit none
  real(kind=prec), allocatable, dimension(:,:) :: injections, dat
  integer, allocatable, dimension(:) :: offsets

contains

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

  close(unit=8)
  END SUBROUTINE LOAD_DATA

  SUBROUTINE TEST
  END SUBROUTINE TEST

                          !!!!!!!!!!!!!!!!!!!!!!
                         END MODULE LOGLIKELIHOOD
                          !!!!!!!!!!!!!!!!!!!!!!
