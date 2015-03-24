.PHONY: all clean tests doc

SOURCES = src/libhades.c  src/odeint.c  src/optimize.c
TESTS   = tests/
OPT     = -O3
CC      = gcc

all:
	${CC} -std=c99 -Wall -Werror ${OPT} -I src/ -c -fpic src/libhades.c -o libhades.o
	${CC} -std=c99 -Wall -Werror ${OPT} -I src/ -c -fpic src/optimize.c -o optimize.o
	${CC} -std=c99 -Wall -Werror ${OPT} -I src/ -c -fpic src/odeint.c   -o odeint.o
	${CC} -shared -Wl,-soname,libhades.so -o libhades.so libhades.o optimize.o odeint.o

install:
	cp libhades.so ~/libs

tests:
	${CC} -std=c99 -Wall -O3 ${TESTS}/tests.c -I ${SRC} -L . -I /:/home/physik/theo1/hartmmic/libs -L /:/home/physik/theo1/hartmmic/libs -lhades -lunittest -llapacke -llapack -lblas -lgfortran -o libhades_tests

clean:
	rm -rf libhades.so libhades.o

doc:
	doxygen Doxyfile
