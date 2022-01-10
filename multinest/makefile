FC=gfortran
FFLAGS=
LD=gfortran
LFLAGS=

OBJ=models.o functions.o loglikelihood.o

lib: MultiNest/modules/nested.mod MultiNest/lib/libmultinest.a

FFLAGS+=-IMultiNest/modules/
LFLAGS+=-llapack -lblas

MultiNest/build/Makefile:
	cd MultiNest/build && cmake ..

MultiNest/modules/nested.mod: MultiNest/build/Makefile
	make -C MultiNest/build

MultiNest/lib/libmultinest.a: MultiNest/build/Makefile
	make -C MultiNest/build

%.mod %.o: %.f95
	$(FC) $(FFLAGS) -c $<

models.o: functions.mod
loglikelihood.o: models.mod
test.o: models.mod loglikelihood.mod
main.o: loglikelihood.mod

test: test.o $(OBJ)
	$(LD) $(LFLAGS) -o $@ $^
main: main.o $(OBJ) MultiNest/lib/libmultinest.a
	$(LD) $(LFLAGS) -o $@ $^


clean:
	rm -f *.o *.mod test main
