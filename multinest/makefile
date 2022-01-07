FC=gfortran
FFLAGS=
LD=gfortran
LFLAGS=

OBJ=models.o functions.o

lib: MultiNest/modules/nested.mod MultiNest/lib/libmultinest.a

FFLAGS+=-IMultiNest/modules/
LFLAGS+= MultiNest/lib/libmultinest.a

MultiNest/build/Makefile:
	cd MultiNest/build && cmake ..

MultiNest/modules/nested.mod: MultiNest/build/Makefile
	make -C MultiNest/build

MultiNest/lib/libmultinest.a: MultiNest/build/Makefile
	make -C MultiNest/build

%.mod %.o: %.f95
	$(FC) $(FFLAGS) -c $<

models.o: functions.mod

clean:
	rm -f *.o *.mod
