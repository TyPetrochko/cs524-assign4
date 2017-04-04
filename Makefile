CC = icc
CFLAGS = -g -O3 -xHost -fno-alias -std=c99
FC = ifort
FFLAGS = -g -O3 -xHost -fno-alias
MPICC = mpicc

all: parallel serial

parallel: parallel.o /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(MPICC) -o $@ $(CFLAGS) $^

serial:	serial.o /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(CC) -o $@ $(CFLAGS) $^

fserial: fserial.o /home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(FC) -o $@ $(FFLAGS) $^

.c.o:
	$(CC) $(CFLAGS) -c $<

.f.o:
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f serial fserial *.o

run:
	./build.sh
	./run.sh

