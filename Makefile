#!!! set the directory where lib/libfftw3.a and include/fftw3.h are included
FFTWDIR=/usr/local

CC	= icc
CFLAGS 	= -O2
INCDIRS_FFTW	= -I$(FFTWDIR)/include
LIBS = -L$(FFTWDIR)/lib -lfftw3

PROG	= pk_recon

.SUFFIXES: .cpp

.cpp.o	:
	$(CC) $(INCDIRS_FFTW) -c $<

default	: $(PROG)

pk_recon	: pk_recon.o
	$(CC) $(CFLAGS) $^ $(LIBS) -o $@

pk_recon.o : recon.h fourier_ver3.h

clean:
	-rm -f $(PROG) *.o core
tidy:
	-rm -f $(PROG) 


