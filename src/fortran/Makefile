PATH_LIBHADES=../
SOURCES = ${PATH_LIBHADES}/libhades.c  ${PATH_LIBHADES}/expm.c  ${PATH_LIBHADES}/odeint.c  ${PATH_LIBHADES}/optimize.c ${PATH_LIBHADES}/parse_npy_dict.c
CFLAGS += -std=c99 -malign-double -Wall -Wmissing-prototypes -Wstrict-prototypes -I${PATH_LIBHADES}

all:
	rm -f *.o
	$(CC) $(OPT) $(CFLAGS) $(LDFLAGS) -c ${SOURCES} expm_wrapper.c
	ld -r *.o -o expm_fortran.o
	rm -r `ls *.o | grep -v expm_fortran.o`

test:
	f95 -O3 -Wall test.f expm_fortran.o -lm -llapack -lblas -o test

clean:
	rm -f *.o
