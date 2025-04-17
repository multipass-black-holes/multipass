project(
  'multinest',
  'fortran'
)

fc = meson.get_compiler('fortran')
add_project_arguments('-ffree-line-length-none', language : 'fortran')

blas = dependency('blas')
lapack = dependency('lapack')
deps = [blas, lapack]


multinest_prefix = 'MultiNest/MultiNest_v3.12_CMake/multinest/src/'
multinest_src = [
  multinest_prefix + 'cwrapper.f90',
  multinest_prefix + 'kmeans_clstr.f90',
  multinest_prefix + 'nested.F90',
  multinest_prefix + 'posterior.F90',
  multinest_prefix + 'priors.f90',
  multinest_prefix + 'utils1.f90',
  multinest_prefix + 'utils.f90',
  multinest_prefix + 'xmeans_clstr.f90',
]
multinest_args = ['-std=legacy']

mpi = dependency('mpi', language: 'fortran', required : false)
if mpi.found()
  deps += mpi
  multinest_args += '-DMPI'
endif
multinest = static_library('multinest', multinest_src, dependencies : deps, fortran_args : multinest_args)

libmultipass = library(
  'multipass',
  [
    'src/functions.f95',
    'src/models.f95',
    'src/loglikelihood.f95',
  ],
  fortran_args : [
    '-cpp'
  ]
)

multipass = executable(
  'main', 'src/main.f95',
  link_with : [
    libmultipass,
    multinest
  ],
  fortran_args : [
    '-cpp'
  ],
  dependencies : deps
)

test('tests', executable(
  'tests', 'src/test.f95',
  link_with : [
    libmultipass,
    multinest
  ],
  fortran_args : [
    '-cpp'
  ]
))

f2py = find_program('f2py', required : false)
if f2py.found()
  add_languages('c')

  py_mod = import('python')
  py = py_mod.find_installation(pure: false)
  py_dep = py.dependency()

  incdir_numpy = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
    check : true
  ).stdout().strip()

  incdir_f2py = run_command(py,
      ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
      check : true
  ).stdout().strip()

  inc_np = include_directories(incdir_numpy, incdir_f2py)

  interface_module = custom_target(
    'interfacemodule',
    input : ['src/interface.f95'],
    output : ['interfacemodule.c', 'interface-f2pywrappers.f', 'interface-f2pywrappers2.f90'],
    command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'interface', '--lower']
  )

  py.extension_module('interface',
    [
      'src/interface.f95',
      interface_module
    ],
    incdir_f2py / 'fortranobject.c',
    include_directories: inc_np,
    dependencies : py_dep,
    link_with : [
      libmultipass,
    ],
  )
endif
