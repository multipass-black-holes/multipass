include makefile.sys

FFLAGS += -g -fPIC -cpp -J .obj -I .obj

EXT_SUFFIX := $(shell python3-config --extension-suffix)

OBJ=.obj/models.o .obj/functions.o .obj/loglikelihood.o

MN_PREFIX=MultiNest/MultiNest_v3.12_CMake/multinest/

FFLAGS+=-I$(MN_PREFIX)/modules/

ifeq ($(MPI),)
MPISUFFIX=
else
MPISUFFIX=_mpi
LFLAGS+= -L$(MPI)/lib/ -lmpi_mpifh -lmpi
endif

all: test main interface$(EXT_SUFFIX)

lib: $(MN_PREFIX)/modules/nested.mod $(MN_PREFIX)/lib/libmultinest$(MPISUFFIX).a

$(MN_PREFIX)/build/Makefile:
	mkdir -p $(MN_PREFIX)/build
	cd $(MN_PREFIX)/build && $(CMAKE) .. -DCMAKE_Fortran_COMPILER=$(FC)

$(MN_PREFIX)/modules/nested.mod: $(MN_PREFIX)/build/Makefile
	make -C $(MN_PREFIX)/build/ multinest$(MPISUFFIX)_static
	cp $(MN_PREFIX)/build/src/nested.mod $@

$(MN_PREFIX)/lib/libmultinest$(MPISUFFIX).a: $(MN_PREFIX)/build/Makefile
	make -C $(MN_PREFIX)/build/ multinest$(MPISUFFIX)_static

.obj/%.mod .obj/%.o: src/%.f95
	$(FC) $(FFLAGS) -o $(@:.mod=.o) -c $<

.obj/models.o: .obj/functions.mod
.obj/loglikelihood.o: .obj/models.mod
.obj/test.o: .obj/models.mod .obj/loglikelihood.mod
.obj/main.o: .obj/loglikelihood.mod $(MN_PREFIX)/modules/nested.mod

test: .obj/test.o $(OBJ)
	$(LD) $(LFLAGS) -o $@ $^
main: .obj/main.o $(OBJ) $(MN_PREFIX)/lib/libmultinest$(MPISUFFIX).a
	$(LD)  -o $@ $^ $(LFLAGS)

interface$(EXT_SUFFIX): src/interface.f95 .obj/loglikelihood.o
	$(F2PY) -m interface -I.obj -c .obj/models.o .obj/functions.o .obj/loglikelihood.o $<

clean:
	make -C $(MN_PREFIX)/build clean || true
	rm -fr .obj/*.o .obj/*.mod test main $(MN_PREFIX)/build *.so
