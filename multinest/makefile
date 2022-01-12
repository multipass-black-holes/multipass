FC=gfortran
FFLAGS=
LD=gfortran
LFLAGS=

MPI=/usr/lib64/openmpi3
CMAKE=cmake3

OBJ=models.o functions.o loglikelihood.o

FFLAGS+=-IMultiNest/modules/
LFLAGS+=-llapack -lblas

ifeq ($(MPI),)
MPISUFFIX=""
else
MPISUFFIX=_mpi
LFLAGS+= -L$(MPI)/lib/ -lmpi_mpifh -lmpi
endif

all: test main

lib: MultiNest/modules/nested.mod MultiNest/lib/libmultinest$(MPISUFFIX).a


MultiNest/build/Makefile:
	cd MultiNest/build && $(CMAKE) .. -DCMAKE_Fortran_COMPILER=$(FC)

MultiNest/modules/nested.mod: MultiNest/build/Makefile
	make -C MultiNest/build

MultiNest/lib/libmultinest$(MPISUFFIX).a: MultiNest/build/Makefile
	make -C MultiNest/build

%.mod %.o: %.f95
	$(FC) $(FFLAGS) -c $<

models.o: functions.mod
loglikelihood.o: models.mod
test.o: models.mod loglikelihood.mod
main.o: loglikelihood.mod MultiNest/modules/nested.mod

test: test.o $(OBJ)
	$(LD) $(LFLAGS) -o $@ $^
main: main.o $(OBJ) MultiNest/lib/libmultinest$(MPISUFFIX).a
	$(LD)  -o $@ $^ $(LFLAGS)


clean:
	rm -fr *.o *.mod test main MultiNest/build/{CMake*,cmake*,Makefile,src}
