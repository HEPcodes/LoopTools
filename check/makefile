OBJS = A0_check.o B0_check.o B1_check.o B00_check.o B11_check.o \
  DB0_check.o DB1_check.o DB00_check.o DB11_check.o \
  C0_check.o D0_check.o aux_B.o aux_CD.o ini.o

all: $(OBJS)

.F.o: defs.h
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJS)

