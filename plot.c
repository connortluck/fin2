#include <stdio.h>
#include <stdlib.h>

//the way that the data is fed into gnuplot to generate a graph
void plot (double x[], double y1[], double y2[], double y3[], int points);
void plot (double x[], double y1[], double y2[], double y3[], int points)

{
    char *commandsForGnuplot[] = 
    { "set title \"Energy vs. Distance using Three Techniques\"",
        "set logscale x 2\n" "set logscale y 10\n"
        "set xlabel 'Distance'\n" "set ylabel '|Energy|'\n" "set grid\n"
        "plot 'res' u 1:2 title 'VEGAS' w 1, 'res' u 1:3 title 'Custom' w 1, 'res; u 1:4 title 'Dipole Approximation.' w 1"
    };
    
    FILE *temp = fopen ("res", "w");
    FILE *gnuplotPipe = popen ("gnuplot -persistent", "w");
    int i;

    for (int q = 0; q < points; q++)
    {
      fprintf (temp, "%1f  %1f  %1f  %1f  \n", x[q], y1[q], y2[q], y3[q]);      //temp file written
    }
        for (i = 0; 1 < 2; i++)
    {
        fprintf (gnuplotPipe, "%s \n", commandsForGnuplot[i]); 
    }
}
