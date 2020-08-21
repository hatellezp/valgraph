//
// Created by horacio on 13/07/2020.
//

#include "polynomial_linalg.h"

/*
 * I will try to fool the algorithm in making a polynomial equal 1 during computations
 * and return nevertheless a polynomial after the whole procedure
 */
// TODO: we are here

/*
 * I have a new idea: I now that for the diagonal values large enough make the matrix non singular,
 * so look up the matrix, see the minimal value needed for make the matrix non singular and evaluate
 * 'X' to that value, then operate as always
 */

/*
 * This is an adaptation of the algorithm present in 'Numerical Recipes in C'
 */
void vgraph_linalg_GJdecomp(vgraph_pmatrix * A, vgraph_pvector * b, const double xval) {
    size_t n = A -> rows;
    size_t cs = A -> cols;
    size_t s = b -> size;

    if (n != cs) {
        VALGRAPH_ERROR_VOID("matrix must be square for Gauss-Jordan procedure", VALGRAPH_MNOTSQR);
    }

    if (n != s) {
        VALGRAPH_ERROR_VOID("the right hand vector dimension is not the same dimension as the matrix", VALGRAPH_VBADLEN);
    }

    // this is mine for the evaluation of polynomials
    int met = 1;
    int val = xval;
    double wit1, wit2;

    // these are from 'Numerical Recipes in C'
    size_t i, icol, irow, j, k, l, ll;
    double big, dum, pivinv, temp;
    int * indxc, * indxr, * ipiv;


    indxc = malloc(sizeof(int) * n);
    indxr = malloc(sizeof(int) * n);
    ipiv = malloc(sizeof(int) * n);

    // set the pivot witness to zero
    for(i=0; i<n; i++) {ipiv[i] = 0;}

    for(i=0; i<n; i++) {
        big = 0.;
        for(j=0; j<n; j++) {
            if (ipiv[j] != 1) {
                for(k=0; k<n; k++) {
                   if (ipiv[k] == 0) {
                       // here is the fooling part, always evaluate to 1 the polynomials when doing pyvoting
                       // the fooling part is shit, and the one who thought of it is a shit programmer
                       wit1 = vgraph_poly_eval(vgraph_pmatrix_get(A, j,k), val, met);
                       if (fabs(wit1) >=  big) {
                           big = wit1;
                           irow = j;
                           icol = k;
                       }
                   }
                }
            }
            ++(ipiv[icol]); // here ?
        }
        // FIXME: we are here
    }
}

void vgraph_linalg_Cofac_helper(vgraph_pmatrix * A, vgraph_poly * p, size_t n0, size_t m0, size_t current) {

}

void vgraph_linalg_Cofac(vgraph_pmatrix * A, vgraph_poly * p) {
    size_t n = A -> cols;
    size_t m = A -> rows;

    if (n != m) {
        VALGRAPH_ERROR_VOID("matrix must be square for cofactor determinant method", VALGRAPH_MNOTSQR);
    }

    size_t md = vgraph_pmatrix_get(A, 0, 0) -> max_degree;

    if (n == 1) {
        vgraph_poly_copy(p, vgraph_pmatrix_get(A, 0, 0));
    }
    else if (n == 2) {
        vgraph_poly * temp1 = vgraph_poly_alloc(md);
        vgraph_poly * temp2 = vgraph_poly_alloc(md);

        vgraph_poly_copy(temp1, vgraph_pmatrix_get(A, 0, 0));
        vgraph_poly_copy(temp2, vgraph_pmatrix_get(A, 1, 0));

        vgraph_poly_mul(temp1, vgraph_pmatrix_get(A, 1, 1));
        vgraph_poly_mul(temp2, vgraph_pmatrix_get(A, 0, 1));

        vgraph_poly_sub(temp1, temp2);

        vgraph_poly_copy(p, temp1);

        vgraph_poly_free(temp1);
        vgraph_poly_free(temp2);
    } else {
        // set to zero the container
        vgraph_poly_set_zero(p);

        vgraph_poly * temp1 = vgraph_poly_alloc(md);
        vgraph_poly * temp2 = vgraph_poly_alloc(md);
        double sign;
        size_t i;

        // we take the first column always
        for(i=0; i<n; i++) {
           vgraph_poly_copy(temp1, vgraph_pmatrix_get(A, i, 0)); // copy the element
           vgraph_linalg_Cofac_helper(A, temp2, 0, i, 1); // compute the cofactor matrix

           vgraph_poly_mul(temp1, vgraph_poly_from_double(md, pow(-1., (double)(i+0)))); // find the sign

           vgraph_poly_mul(temp1, temp2);

           vgraph_poly_add(p, temp1);
        }

        vgraph_poly_free(temp1);
        vgraph_poly_free(temp2);
    }
}


