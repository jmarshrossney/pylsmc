SHELL = /bin/sh
ENV = $(CONDAHOME)/envs/lsmc

# Compilers and flags
F90       = gfortran
CC        = gcc
SWIG      = swig
LD        = ld
FFLAGS    = -fPIC -O3
CFLAGS    = -fPIC -O3 -I$(ENV)/include/python3.9/ -I$(ENV)/lib/python3.9/site-packages/numpy/core/include/

# Objects to build into a library
OBJECTS    = util.o energy.o energy_wrap.o

# Main build target - builds the Python-wrapped Fortran library
pot: $(OBJECTS) 

	$(F90) -shared -o _energy.so $(OBJECTS) 

.PRECIOUS: %.o %_wrap.c
.PHONY:  clean

%: %.o
%.o: %.f90
	$(F90) $(FFLAGS) -c -o $@ $<

%.o: %.F90
	$(F90) $(FFLAGS) -c -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

%_wrap.c: %.i %.h
	$(SWIG) -python $<

clean : 
	rm -f _pot.so *.o energy_wrap.* energy.py energy2.py _energy2.so _energy.so *.mod *.pyc


