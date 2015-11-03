all:
	gcc -O2 -Wall invert.c    -lm -lhades -llapack -lblas -o invert
	gcc -O2 -Wall solve.c     -lm -lhades -llapack -lblas -o solve
	gcc -O2 -Wall expm.c      -lm -lhades -llapack -lblas -o expm
	gcc -O2 -Wall lu_lapack.c -lm -lhades -llapack -lblas -o lu_lapack

clean:
	rm -f invert solve expm lu_lapack
