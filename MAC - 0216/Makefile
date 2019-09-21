all: write

run: sw1.c
	gcc sw1.c -o sw1
	valgrind --leak-check=yes ./sw1 < e1

write: sw1.c plot results/nave1
	gcc sw1.c -o sw1 -lm
	./sw1 < e1 > results/nave1
	gcc plot.c -o plot	
	./plot < results/nave1

