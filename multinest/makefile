FC=gfortran
FFLAGS=
LD=gfortran
LFLAGS=

lib: MultiNest/modules/nested.mod MultiNest/lib/libmultinest.a

FFLAGS+=-IMultiNest/modules/
LFLAGS+=-LMultiNest/lib/ -lmultinest

MultiNest/build/Makefile:
	cd MultiNest/build && cmake ..

MultiNest/modules/nested.mod: MultiNest/build/Makefile
	make -C MultiNest/build

MultiNest/lib/libmultinest.a: MultiNest/build/Makefile
	make -C MultiNest/build
