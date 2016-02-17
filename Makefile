.PHONY: all clean tests doc

TESTS   = tests/

CC      = gcc
OPT     = -O3 -march=native -flto
CFLAGS  = -std=c99
CFLAGS += -Wall -Wextra -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wwrite-strings -Wno-unused-parameter
CFLAGS += -I src/

all:
	$(CC) $(CFLAGS) $(OPT) -c -fpic src/libhades.c -o libhades.o
	$(CC) $(CFLAGS) $(OPT) -c -fpic src/optimize.c -o optimize.o
	$(CC) $(CFLAGS) $(OPT) -c -fpic src/odeint.c   -o odeint.o
	$(CC) $(CFLAGS) $(OPT) -c -fpic src/expm.c     -o expm.o
	$(CC) $(CFLAGS) $(OPT) -c -fpic src/parse_npy_dict.c -o parse_npy_dict.o
	$(CC) $(OPT) -shared -Wl,-soname,libhades.so -o libhades.so libhades.o optimize.o odeint.o expm.o parse_npy_dict.o

install:
	cp libhades.so /usr/lib
	cp src/libhades*.h /usr/include
	cp -R src/libhades /usr/include

tests:
	${CC} -std=c99 -Wall -O3 ${TESTS}/tests.c -I ${SRC} -L . -I /:/home/physik/theo1/hartmmic/libs -L /:/home/physik/theo1/hartmmic/libs -lm -lhades -lunittest -llapacke -llapack -lblas -lgfortran -o libhades_tests

clean:
	rm -rf *.so *.o

doc:
	doxygen Doxyfile
