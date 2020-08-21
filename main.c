//
// Created by horacio on 30/06/2020.
//

#include <stdio.h>

#include "core/polynomial.h"
#include "core/polynomial_vector.h"
#include "core/polynomial_matrix.h"



int main() {
    double * ar1 = malloc(sizeof(double) * 4);
    double * ar2 = malloc(sizeof(double) * 4);

    ar1[0] = 0.;
    ar2[1] = 0.;
    ar1[1] = 1.;
    ar2[0] = 1.;

    ar1[2] = ar2[3] = 4.;
    ar1[3] = ar2[2] = 5.;

    vgraph_poly * p1 = vgraph_poly_from_array(ar1, 4);
    vgraph_poly * p2 = vgraph_poly_from_array(ar2, 4);

    vgraph_poly_print(p1); // 1*x + 0
    printf("\n"); // 0*x + 1

    vgraph_poly_print(p2);
    printf("\n");

    vgraph_poly_swap(&p1, &p2);

    printf("----------------------\n");

    vgraph_poly_print(p1);
    printf("\n");

    vgraph_poly_print(p2);
    printf("\n");


    vgraph_poly_free(p1);
    vgraph_poly_free(p2);

    return 0;
}
