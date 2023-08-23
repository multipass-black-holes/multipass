  PROGRAM TEST
  use models, only: test_models => test
  use loglikelihood, only: test_ll => test
  implicit none

  call test_models
  call test_ll
  END PROGRAM TEST
