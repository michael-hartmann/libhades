all:
	gcc -O2 -Wall invert.c -lm -lhades -llapack -lblas -o invert
	gcc -O2 -Wall solve.c  -lm -lhades -llapack -lblas -o solve
	gcc -O2 -Wall expm.c   -lm -lhades -llapack -lblas -o expm

clean:
	rm -f invert solve expm
