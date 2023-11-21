include makefile.sys

FFLAGS += -g -fPIC -cpp -J .obj -I .obj

EXT_SUFFIX := $(shell python3-config --extension-suffix)

OBJ=.obj/models.o .obj/functions.o .obj/loglikelihood.o

FFLAGS+=-IMultiNest/modules/

ifeq ($(MPI),)
MPISUFFIX=
else
MPISUFFIX=_mpi
LFLAGS+= -L$(MPI)/lib/ -lmpi_mpifh -lmpi
endif

all: test main interface$(EXT_SUFFIX)

lib: MultiNest/modules/nested.mod MultiNest/lib/libmultinest$(MPISUFFIX).a


MultiNest/build/Makefile:
	cd MultiNest/build && $(CMAKE) .. -DCMAKE_Fortran_COMPILER=$(FC)

MultiNest/modules/nested.mod: MultiNest/build/Makefile
	make -C MultiNest/build

MultiNest/lib/libmultinest$(MPISUFFIX).a: MultiNest/build/Makefile
	make -C MultiNest/build

.obj/%.mod .obj/%.o: src/%.f95
	$(FC) $(FFLAGS) -o $(@:.mod=.o) -c $<

.obj/models.o: .obj/functions.mod
.obj/loglikelihood.o: .obj/models.mod
.obj/test.o: .obj/models.mod .obj/loglikelihood.mod
.obj/main.o: .obj/loglikelihood.mod MultiNest/modules/nested.mod

test: .obj/test.o $(OBJ)
	$(LD) $(LFLAGS) -o $@ $^
main: .obj/main.o $(OBJ) MultiNest/lib/libmultinest$(MPISUFFIX).a
	$(LD)  -o $@ $^ $(LFLAGS)

interface$(EXT_SUFFIX): src/interface.f95 .obj/loglikelihood.o
	$(F2PY) -m interface -I.obj -c .obj/models.o .obj/functions.o .obj/loglikelihood.o $<

clean:
	rm -fr .obj/*.o .obj/*.mod test main MultiNest/build/{CMake*,cmake*,Makefile,src} *.so
