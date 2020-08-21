//
// Created by horacio on 25/06/2020.
//

#include <stdio.h>

#include "core/time_eval.h"
#include "core/gnuplot_i.h"


int main() {
    int beg = 10;
    int end = 50;
    int n = end - beg + 1;
    double limit = 10;
    int reps = 50;

    // I'm going to rescale following the number of repetitions
    double dreps = (double)reps;

    printf("Settings:\n    N range: [%d, %d]\n    Bound for values in matrice: %f\n    Repetitions for single execution stats: %d\n", beg, end, limit, reps);

    // compute the total time taken
    printf("Test begins now, this make take some time...\n");

    // first method zero
    printf("Full vector computation method\n");
    clock_t t = clock();
    double * stats0 = execution_time_array(beg, end, limit, reps, FULL_V);
    t = clock() - t;
    double tt = ((double) t) / CLOCKS_PER_SEC;
    printf("\ntotal time taken with method full_v: %f\n\n", tt);
    printf("\ntotal time taken with method full_v after descaling: %f\n\n", tt / dreps);

    /*
    // second, method 1
    printf("Single value computation method\n");
    t = clock();
    double * stats1 = execution_time_array(beg, end, limit, reps, SINGLE_V);
    t = clock() - t;
    tt = ((double) t) / CLOCKS_PER_SEC;
    printf("\ntotal time taken with method single_v: %f\n\n", tt);
    printf("\ntotal time taken with method single_v after descaling: %f\n\n", tt / dreps);
    */

    /*
    for(int i=0; i<n; i++) {
        printf("N=%d => min=%f, avg=%f, max=%f\n", beg+i, stats[i*3], stats[i*3 + 1], stats[i*3 + 2]);
    }
     */

    // for plot
    gnuplot_ctrl * h1;
    h1 = gnuplot_init();
    gnuplot_setstyle(h1, "lines" );

    // created the arrays
    double x[n];
    double tmin0[n], dtmin0[n];
    double tavg0[n], dtavg0[n];
    double tmax0[n], dtmax0[n];

    /*
    double tmin1[n], dtmin1[n];
    double tavg1[n], dtavg1[n];
    double tmax1[n], dtmax1[n];
     */

    for(int i=0; i<n; i++) {
        x[i] = (double)(beg + i);
        tmin0[i] = stats0[i*3];
        tavg0[i] = stats0[i*3 + 1];
        tmax0[i] = stats0[i*3 + 2];

        /*
        tmin1[i] = stats1[i*3];
        tavg1[i] = stats1[i*3 + 1];
        tmax1[i] = stats1[i*3 + 2];
        */

        dtmin0[i] = stats0[i*3] / dreps;
        dtavg0[i] = stats0[i*3 + 1] / dreps;
        dtmax0[i] = stats0[i*3 + 2] / dreps;

        /*
        dtmin1[i] = stats1[i*3] / dreps;
        dtavg1[i] = stats1[i*3 + 1] / dreps;
        dtmax1[i] = stats1[i*3 + 2] / dreps;
         */
    }

    gnuplot_plot_xy(h1, x, tmin0, n, "minimal time full_v");
    gnuplot_plot_xy(h1, x, tavg0, n, "average time full_v");
    gnuplot_plot_xy(h1, x, tmax0, n, "maximal time full_v");

    /*
    gnuplot_plot_xy(h1, x, tmin1, n, "minimal time single_v");
    gnuplot_plot_xy(h1, x, tavg1, n, "average time single_v");
    gnuplot_plot_xy(h1, x, tmax1, n, "maximal time single_v");
    */

    // after descaling
    gnuplot_plot_xy(h1, x, dtmin0, n, "afdes: minimal time full_v");
    gnuplot_plot_xy(h1, x, dtavg0, n, "afdes: average time full_v");
    gnuplot_plot_xy(h1, x, dtmax0, n, "afdes: maximal time full_v");

    /*
    gnuplot_plot_xy(h1, x, dtmin1, n, "afdes: minimal time single_v");
    gnuplot_plot_xy(h1, x, dtavg1, n, "afdes: average time single_v");
    gnuplot_plot_xy(h1, x, dtmax1, n, "afdes: maximal time single_v");
     */

    // wait for key
    getchar();

    // clean and return
    free(stats0);
    // free(stats1);
    gnuplot_close(h1);

    return 0;
}
