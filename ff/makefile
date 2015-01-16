OBJS = ff2dl2.o ffabcd.o ffcc0p.o ffcel2.o ffcli2.o ffcrr.o ffcxr.o \
  ffcxs3.o ffcxs4.o ffcxyz.o ffdcc0.o ffdcxs.o ffdel2.o ffdel3.o \
  ffdel4.o ffdxc0.o ffinit.o fftran.o ffxa0.o ffxb0.o ffxb1.o \
  ffxb2p.o ffxc0.o ffxc0i.o ffxc0p.o ffxc0p0.o ffxd0.o ffxd0h.o \
  ffxd0i.o ffxd0m0.o ffxd0p.o ffxdb0.o ffxdb1.o ffxdbd.o ffxli2.o \
  ffxxyz.o

all: $(OBJS)

.f.o: ff.h
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJS)

