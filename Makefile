.PHONY: all clean tests doc

SRC=src/
TESTS=tests/
CC=gcc

all:
	${CC} -std=c99 -I src/ -c -fpic -Wall -O3 ${SRC}/libhades.c -o libhades.o
	${CC} -shared -Wl,-soname,libhades.so -o libhades.so libhades.o

install:
	cp libhades.so ${SRC}/libhades.h ~/libs

tests:
	${CC} -std=c99 -Wall -O3 ${TESTS}/tests.c -I ${SRC} -L . -I /:/home/physik/theo1/hartmmic/libs -L /:/home/physik/theo1/hartmmic/libs -lhades -lunittest -llapacke -llapack -lblas -lgfortran -o libhades_tests

clean:
	rm -rf libhades.so libhades.o

doc:
	doxygen Doxyfile
