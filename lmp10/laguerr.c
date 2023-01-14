#include "makespl.h"
#include "piv_ge_solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

/* UWAGA: liczbę używanych f. bazowych można ustawić przez wartość
          zmiennej środowiskowej APPROX_BASE_SIZE
*/

/*
 * Funkcje bazowe: n - liczba funkcji a,b - granice przedzialu aproksymacji i
 * - numer funkcji x - wspolrzedna dla ktorej obliczana jest wartosc funkcji
 */

double fi( double a, double b, int n, int i, double x )
{
	switch (i)
    {
		case 0:
		  return 1;
		case 1:
		  return -x + 1;
		case 2:
		  return ( ( pow( x, 2 ) - 4 * x + 2 ) / 2 );
		case 3:
		  return ( ( -1 * pow( x, 3 ) + 9 * pow( x, 2 ) - 18 * x + 6 ) / 6 );
		case 4:
		  return ( ( pow( x, 4 ) - 16 * pow( x, 3 ) + 72 * pow( x, 2 ) - 96 * x + 24 ) / 24 );
		case 5:
		  return ( ( -1 * pow( x, 5 ) + 25 * pow( x, 4 ) - 200 * pow( x, 3 ) + 600 * pow( x, 2 ) - 600 * x + 120 ) / 120 );
		case 6:
		  return ( ( pow( x, 6 ) - 36 * pow( x, 5 ) + 450 * pow( x, 4 ) - 2400 * pow( x, 3 ) + 5400 * pow( x, 2 ) - 4320 * x + 720 ) / 720 );
		case 7:
		  return ( ( -1 * pow( x, 7 ) + 49 * pow( x, 6 ) - 882 * pow( x, 5 ) + 7350 * pow( x, 4 ) - 29400 * pow( x, 3 ) + 52920 * pow( x, 2 ) - 35280 * x + 5040 ) / 5040 );
		case 8:
		  return ( ( pow( x, 8 ) - 64 * pow( x, 7 ) + 1568 * pow( x, 6 ) - 18816 * pow( x, 5 ) + 117600 * pow( x, 4 ) - 376320 * pow( x, 3 ) + 564480 * pow( x, 2 ) - 322560 * x + 40320 ) / 40320 );
		case 9:
		  return ( ( -1 * pow( x, 9 ) + 81 * pow( x, 8 ) - 2592 * pow( x, 7 ) + 42336 * pow( x, 6 ) - 381024 * pow( x, 5 ) + 1905120 * pow( x, 4 ) - 5080320 * pow( x, 3 ) + 6531840 * pow( x, 2 ) - 3265920 * x + 362880 ) / 362880 );
		default:
		  return 0;
		}
}

/* Pierwsza pochodna fi */
double
dfi(double a, double b, int n, int i, double x)
{
	switch (i)
    {
		case 0:
		  return 0;
		case 1:
		  return -1;
		case 2:
		  return x - 2;
		case 3:
		  return ( ( -3 * pow( x, 2 ) + 18 * x - 18 ) / 6 );
		case 4:
		  return ( ( 4 * pow( x, 3 ) - 48 * pow( x, 2 ) + 144 * x - 96 ) / 24 );
		case 5:
		  return ( ( -5 * pow( x, 4 ) + 100 * pow( x, 3 ) - 600 * pow( x, 2 ) + 1200 * x - 600 ) / 120 );
		case 6:
		  return ( ( 6 * pow( x, 5 ) - 180 * pow( x, 4 ) + 1800 * pow( x, 3 ) - 7200 * pow( x, 2 ) + 10800 * x - 4320 ) / 720 );
		case 7:
		  return ( ( -7 * pow( x, 6 ) + 294 * pow( x, 5 ) - 4410 * pow( x, 4 ) + 29400 * pow( x, 3 ) - 88200 * pow( x, 2 ) + 105840 * x - 35280 ) / 5040 );
		case 8:
		  return ( ( 8 * pow( x, 7 ) - 64 * 7 * pow( x, 6 ) + 1568 * 6 * pow( x, 5 ) - 18816 * 5 * pow( x, 4 ) + 117600 * 4 * pow( x, 3 ) - 376320 * 3 * pow( x, 2 ) + 564480 * 2 * x - 322560 ) / 40320 );
		case 9:
		  return ( ( -9 * pow( x, 8 ) + 81 * 8 * pow( x, 7 ) - 2592 * 7 * pow( x, 6 ) + 42336 * 6 * pow( x, 5 ) - 381024 * 5 * pow( x, 4 ) + 1905120 * 4 * pow( x, 3 ) - 5080320 * 3 * pow( x, 2 ) + 6531840 * 2 * x - 3265920 ) / 362880 );
		default:
		  return 0;
    }
}

/* Druga pochodna fi */
double
d2fi(double a, double b, int n, int i, double x)
{
	switch (i)
    {
		case 0:
		  return 0;
		case 1:
		  return 0;
		case 2:
		  return 1;
		case 3:
		  return -x + 3;
		case 4:
		  return ( ( 12 * pow( x, 2 ) - 96 * x + 144 ) / 24 );
		case 5:
		  return ( ( -20 * pow( x, 3 ) + 300 * pow( x, 2 ) - 1200 * x + 1200 ) / 120 );
		case 6:
		  return ( ( 30 * pow( x, 4 ) - 720 * pow( x, 3 ) + 5400 * pow( x, 2 ) - 14400 * x + 10800 ) / 720 );
		case 7:
		  return ( ( -42 * pow( x, 5 ) + 294 * 5 * pow( x, 4 ) - 4410 * 4 * pow( x, 3 ) + 29400 * 3 * pow( x, 2 ) - 88200 * 2 * x + 105840 ) / 5040 );
		case 8:
		  return ( ( 56 * pow( x, 6 ) - 64 * 7 * 6 * pow( x, 5 ) + 1568 * 6 * 5 * pow( x, 4 ) - 18816 * 5 * 4 * pow( x, 3 ) + 117600 * 4 * 3 * pow( x, 2 ) - 376320 * 3 * 2 * x + 564480 * 2 ) / 40320 );
		case 9:
		  return ( ( -72 * pow( x, 7 ) + 81 * 8 * 7 * pow( x, 6 ) - 2592 * 7 * 6 * pow( x, 5 ) + 42336 * 6 * 5 * pow( x, 4 ) - 381024 * 5 * 4 * pow( x, 3 ) + 1905120 * 4 * 3 * pow( x, 2 ) - 5080320 * 3 * 2 * x + 6531840 * 2 ) / 362880 );
		default:
		  return 0;
    }
}

/* Trzecia pochodna fi */
double
d3fi(double a, double b, int n, int i, double x)
{
	switch (i)
    {
		case 0:
		  return 0;
		case 1:
		  return 0;
		case 2:
		  return 0;
		case 3:
		  return -1;
		case 4:
		  return x - 4;
		case 5:
		  return ( ( -60 * pow( x, 2 ) + 600 * x - 1200 ) / 120 );
		case 6:
		  return ( ( 120 * pow( x, 3 ) - 2160 * pow( x, 2 ) + 10800 * x - 14400 ) / 720 );
		case 7:
		  return ( ( -210 * pow( x, 4 ) + 294 * 5 * 4 * pow( x, 3 ) - 4410 * 4 * 3 * pow( x, 2 ) + 29400 * 3 * 2 * x - 88200 * 2 ) / 5040 );
		case 8:
		  return ( ( 56 * 6 * pow( x, 5 ) - 64 * 7 * 6 * 5 * pow( x, 4 ) + 1568 * 6 * 5 * 4 * pow( x, 3 ) - 18816 * 5 * 4 * 3 * pow( x, 2 ) + 117600 * 4 * 3 * 2 * x - 376320 * 3 * 2 ) / 40320 );
		case 9:
		  return ( ( -72 * 7 * pow( x, 6 ) + 81 * 8 * 7 * 6 * pow( x, 5 ) - 2592 * 7 * 6 * 5 * pow( x, 4 ) + 42336 * 6 * 5 * 4 * pow( x, 3 ) - 381024 * 5 * 4 * 3 * pow( x, 2 ) + 1905120 * 4 * 3 * 2 * x - 5080320 * 3 * 2 ) / 362880 );
		default:
		  return 0;
    }
}

/* Pomocnicza f. do rysowania bazy */
void
xfi(double a, double b, int n, int i, FILE *out)
{
	double		h = (b - a) / (n - 1);
	double		h3 = h * h * h;
	int		hi         [5] = {i - 2, i - 1, i, i + 1, i + 2};
	double		hx      [5];
	int		j;

	for (j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];

	fprintf( out, "# nb=%d, i=%d: hi=[", n, i );
	for( j= 0; j < 5; j++ )
		fprintf( out, " %d", hi[j] );
	fprintf( out, "] hx=[" );
	for( j= 0; j < 5; j++ )
		fprintf( out, " %g", hx[j] );
	fprintf( out, "]\n" );
}

void
make_spl(points_t * pts, spline_t * spl)
{

	matrix_t       *eqs= NULL;
	double         *x = pts->x;
	double         *y = pts->y;
	double		a = x[0];
	double		b = x[pts->n - 1];
	int		i, j, k;
	int		nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
  char *nbEnv= getenv( "APPROX_BASE_SIZE" );

	if( nbEnv != NULL && atoi( nbEnv ) > 0 )
		nb = atoi( nbEnv );

	eqs = make_matrix(nb, nb + 1);

#ifdef DEBUG
#define TESTBASE 500
	{
		FILE           *tst = fopen("debug_base_plot.txt", "w");
		double		dx = (b - a) / (TESTBASE - 1);
		for( j= 0; j < nb; j++ )
			xfi( a, b, nb, j, tst );
		for (i = 0; i < TESTBASE; i++) {
			fprintf(tst, "%g", a + i * dx);
			for (j = 0; j < nb; j++) {
				fprintf(tst, " %g", fi  (a, b, nb, j, a + i * dx));
				fprintf(tst, " %g", dfi (a, b, nb, j, a + i * dx));
				fprintf(tst, " %g", d2fi(a, b, nb, j, a + i * dx));
				fprintf(tst, " %g", d3fi(a, b, nb, j, a + i * dx));
			}
			fprintf(tst, "\n");
		}
		fclose(tst);
	}
#endif

	for (j = 0; j < nb; j++) {
		for (i = 0; i < nb; i++)
			for (k = 0; k < pts->n; k++)
				add_to_entry_matrix(eqs, j, i, fi(a, b, nb, i, x[k]) * fi(a, b, nb, j, x[k]));

		for (k = 0; k < pts->n; k++)
			add_to_entry_matrix(eqs, j, nb, y[k] * fi(a, b, nb, j, x[k]));
	}

#ifdef DEBUG
	write_matrix(eqs, stdout);
#endif

	if (piv_ge_solver(eqs)) {
		spl->n = 0;
		return;
	}
#ifdef DEBUG
	write_matrix(eqs, stdout);
#endif

	if (alloc_spl(spl, 1) == 0) {
		for (i = 0; i < spl->n; i++) {
			double xx = spl->x[i] = ( a + b ) / 2;
			xx+= 10.0*DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (k = 0; k < nb; k++) {
				double		ck = get_entry_matrix(eqs, k, nb);
				spl->f[i]  += ck * fi  (a, b, nb, k, xx);
				spl->f1[i] += ck * dfi (a, b, nb, k, xx);
				spl->f2[i] += ck * d2fi(a, b, nb, k, xx);
				spl->f3[i] += ck * d3fi(a, b, nb, k, xx);
			}
		}
	}

#ifdef DEBUG
	{
		FILE           *tst = fopen("debug_spline_plot.txt", "w");
		double		dx = (b - a) / (TESTBASE - 1);
		for (i = 0; i < TESTBASE; i++) {
			double yi= 0;
			double dyi= 0;
			double d2yi= 0;
			double d3yi= 0;
			double xi= a + i * dx;
			for( k= 0; k < nb; k++ ) {
							yi += get_entry_matrix(eqs, k, nb) * fi(a, b, nb, k, xi);
							dyi += get_entry_matrix(eqs, k, nb) * dfi(a, b, nb, k, xi);
							d2yi += get_entry_matrix(eqs, k, nb) * d2fi(a, b, nb, k, xi);
							d3yi += get_entry_matrix(eqs, k, nb) * d3fi(a, b, nb, k, xi);
			}
			fprintf(tst, "%g %g %g %g %g\n", xi, yi, dyi, d2yi, d3yi );
		}
		fclose(tst);
	}
#endif
	
	free_matrix( eqs );
}
