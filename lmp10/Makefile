aproksymator: main.o splines.o points.o aproksymator_na_bazie.o libge
	$(CC) -o aproksymator  main.o splines.o points.o aproksymator_na_bazie.o -L gaus -l ge

aprox: main.o splines.o points.o laguerr.o libge
	$(CC) -o aprox  main.o splines.o points.o laguerr.o -L gaus -l ge -lm

intrp: main.o splines.o points.o interpolator.o libge
	$(CC) -o intrp  main.o splines.o points.o interpolator.o -L gaus -l ge

prosta: main.o splines.o points.o prosta.o
	$(CC) -o prosta  main.o splines.o points.o prosta.o	

laguerr.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c laguerr.c

interpolator.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c interpolator.c

aproksymator_na_bazie.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c aproksymator_na_bazie.c

main.o: main.c
	$(CC) -c -o main.o main.c

points.o: points.c
	$(CC) -c -o points.o points.c

splines.o: splines.c 
	$(CC) -c -o splines.o splines.c

libge: gaus/matrix.o gaus/pivot.o gaus/piv_ge_solver.o
	ar rvs gaus/libge.a gaus/matrix.o gaus/pivot.o gaus/piv_ge_solver.o
	
gaus/matrix.o: 
	$(CC) -c -o gaus/matrix.o gaus/matrix.c
gaus/pivot.o: 
	$(CC) -c -o gaus/pivot.o gaus/pivot.c 
gaus/piv_ge_solver.o: 
	$(CC) -c -o gaus/pivot.o gaus/piv_ge_solver.c

all: aprox intrp prosta aproksymator

.PHONY: clean

clean:
	-rm *.o gaus/*.o gaus/libge.a aprox aproksymator intrp prosta
