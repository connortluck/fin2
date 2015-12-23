#include <stdio.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_math.h>
#include "timer.c"
#include "timer.h"

extern double g (double *t, size_t dim, void *params);

double dipole_approx (double da);

double gaussian (double *x, int dim);

int main (void)
{
    double res, err;
	//initialize variables
    double x1[] = { 0., 0., 0., 0., 0., 0., };
    double xu[] = { 1., 1., 1., 1., 1., 1., };
    double dmin = 1.001;
    double dmax = 4.;
    double dist;
	size_t size = 6;
    int nt = 20;
    double nz = (dmax - dmin) / (nt - 1);
    double vegas[20], di_p[20], distance[20];
	
	//establish random number generator
    gsl_rng *da = gsl_rng_alloc (gsl_rng_taus2);
    unsigned long seed = 1UL;
	gsl_rng_set (da, seed);
	size_t calls = 1000000;
   
    dist = dmin;

    gsl_monte_function G = { &g, size, &dist };

    gsl_monte_vegas_state *sa = gsl_monte_vegas_alloc (size);

    gsl_monte_vegas_init (sa);


    // Monte Carlo Calculations
	//initialize variables
    double sum;
    double x[6];
	long i, j, nn;
	nn = 1000000;
	
   	//start timer and calculation
	timer_start ();
    double home[20];
  
    dist = dmin;
    for (j = 0; j < nt; j++)
    {
        sum = 0.;
        for (i = 0; i < nn; i++)
        {
            // random point in 6 dimensions of the cube
            for (int k = 0; k < (int) size; k++)
            {
                x[k] = gsl_rng_uniform (da);
            }
            sum += g (x, size, &dist);
        }
        res = sum/nn;
        dist += nz;
        home[j] = res;
    }   

    timer_stop ();   

	// Vegas Calculations
    timer_start (); //timer 

    for (int i = 0; i < nt; i++)
    {
        gsl_monte_vegas_integrate (&G, x1, xu, size, calls / 5, da, sa, &res,
            &err);
        do
        {
            gsl_monte_vegas_integrate (&G, x1, xu, size, calls, da, sa, &res,
                &err);
        }
        while (fabs (gsl_monte_vegas_chisq (sa) - 1.0) > 0.2);
       
        
        fflush (stdout);
        dist += nt;
        vegas[i] = res;
    distance[i] = dist;
    di_p[i] = -2. / pow (dist, 3.);
    }

    timer_stop();
   
    gsl_monte_vegas_free (sa);


  
    /*
    printf ("Time for Vegas:   %6f\n", t1);
    printf ("Time for HomeMC:  %6f\n", t2);
    */

    gsl_rng_free (da);

    //Begin calculating error
    double homeerr = 0.0;

    for (int i = 0; i < nt; i++)
    {
    homeerr += fabs(home[i] - vegas[i]);
    }

    printf("#    Dist               Vegas           Home       Dipolapprox\n");
    for( int l = 0; l < nt; l++)
    {
    double xx = distance[l];
    double yy = fabs(vegas[l]);
    double zz = fabs(home[l]);
    double zx = fabs(di_p[l]);
    printf("   %.6f         %.6f       %.6f      %.6f\n", xx, yy, zz, zx);
    }

    /*
    printf("\nAverage error for Home integration: %.6f\n", homeerr_avg);
    printf("Average error for Vegas integration: %.6f\n", err);
    printf("Ratio of the two errors: %.6f\n", differ);
    */


    return 0;
}
