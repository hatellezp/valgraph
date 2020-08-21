#include <stdio.h>

#include "core/time_eval.h"
#include "core/algrank.h"

void print_d_array(double * v, int n) {
    for(int i=0; i<n; i++) {
        printf("%f ", v[i]);
    }
}

void print_d_array_as_p(double * v, int n, double tolerance) {
    int d,i;
    for(i=n-1; i>=0; i--) {
        if (fabs(v[i]) > tolerance) {
            d = i;
            break;
        } else {
            continue;
        }
    }

    for(i=d; i>=1; i--) {
        printf("%f*x^*%d + ", v[i], i);
    }
    printf("%f", v[0]);
}
void print_d_array_as_p2(double * v, int n, double tolerance) {
    printf("%f*x^%d", v[n-1], n-1);
    for(int i=n-2; i>=0; i--) {
        if (fabsl(v[i]) > tolerance) {
            if (v[i] > 0) {
                printf(" + %f*x^*%d", v[i], i);
            } else {
                printf(" %f*x^*%d", v[i], i);
            }
        } else {
            continue;
        }
    }
}

void generate_time_execution_polynomial(int beg, int end, double limit, int reps, double coeff[], int typ) {
    double * stats0 = execution_time_array(beg, end, limit, reps, FULL_V);

    int n = end - beg + 1;
    int i;
    double * x = malloc(sizeof(double)*n);
    double * y = malloc(sizeof(double)*n);
    int offs = typ % 3;;
    double scale = 1000.;
    for(i=0; i<n; i++) {
        x[i] = (double)(beg + i);
        y[i] = scale * ((double)(stats0[i + offs])); // now it depends of typ
    }

    printf("-------------------\n");
    printf("this is x array\n");
    print_d_array(x, n);
    printf("\n");
    printf("this is y array\n");
    print_d_array(y, n);
    printf("\n");
    printf("---------------------\n");

    polynomial_coefficient_computation(x, y, coeff, n);
    free(x); free(y);
}


int main() {
    double tolerance = 0.0001;
    int beg = 2;
    int end = 10;
    int n = end - beg + 1;
    double limit = 100;
    int reps = 1000;
    double * pcoeff = malloc(sizeof(double)*n);

    generate_time_execution_polynomial(beg, end, limit, reps, pcoeff, 1);

    print_d_array(pcoeff, n);
    printf("\n");
    print_d_array_as_p2(pcoeff, n, tolerance);

    free(pcoeff);
    return 0;
}
