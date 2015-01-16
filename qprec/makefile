Q77=../qprec/q77 -I../qprec

all: f77290 qcomplex.o ../libooptools.a ../libooptools_check.a
	-mkdir ../$(HOSTTYPE)
	-mkdir ../$(HOSTTYPE)/lib ../$(HOSTTYPE)/include
	ar r ../libooptools.a qcomplex.o
	mv ../libooptools.a ../$(HOSTTYPE)/lib/libooptoolsQ.a
	ar r ../libooptools_check.a qcomplex.o
	mv ../libooptools_check.a ../$(HOSTTYPE)/lib/libooptoolsQ_check.a
	f77290 ../include/looptools.h ../$(HOSTTYPE)/include/looptools.h90
	-cp qcomplex.mod ../$(HOSTTYPE)/include/

../libooptools.a: f77290 qcomplex.o
	cd .. ; $(MAKE) FC="$(Q77)" FFC="$(Q77)" CC="$(Q77)" libooptools.a

../libooptools_check.a: f77290 qcomplex.o
	cd .. ; $(MAKE) FC="$(Q77)" FFC="$(Q77)" CC="$(Q77)" libooptools_check.a

f77290: f77290.c
	gcc -O -o f77290 f77290.c

qcomplex.o: qcomplex.f90
	f90 -O -c qcomplex.f90

clean:
	cd .. ; $(MAKE) clean
	rm -f f77290 qcomplex.mod qcomplex.o

