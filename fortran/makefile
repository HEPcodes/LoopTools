OBJS = A0.o B0.o B1.o B00.o B11.o DB0.o DB00.o DB1.o DB11.o \
  C0.o C0i.o Cget.o D0.o D0i.o Dget.o ini.o cache.o

all: compiled_with_$(CHECK) $(OBJS)

.F.o: defs.h
	$(FC) $(FFLAGS) -DCHECK=$(CHECK) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

compiled_with_$(CHECK):
	rm -f $(OBJS) compiled_with_*
	touch compiled_with_$(CHECK)

clean:
	rm -f $(OBJS) compiled_with_*

