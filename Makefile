.PHONY: all clean tests doc

TESTS   = tests/
OPT     = -O3
CC      = gcc

all:
	${CC} -std=c99 -Wall ${OPT} -I src/ -c -fpic src/libhades.c -o libhades.o
	${CC} -std=c99 -Wall ${OPT} -I src/ -c -fpic src/optimize.c -o optimize.o
	${CC} -std=c99 -Wall ${OPT} -I src/ -c -fpic src/odeint.c   -o odeint.o
	${CC} -std=c99 -Wall ${OPT} -I src/ -c -fpic src/expm.c     -o expm.o
	${CC} -std=c99 -Wall ${OPT} -I src/ -c -fpic src/parse_npy_dict.c -o parse_npy_dict.o
	${CC} -shared -Wl,-soname,libhades.so -o libhades.so libhades.o optimize.o odeint.o expm.o parse_npy_dict.o

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
