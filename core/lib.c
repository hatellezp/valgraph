//
// Created by horacio on 20/08/2020.
//

#include "lib.h"

// inner and helper functions
void from_real_to_complex_matrix(gsl_matrix * from, gsl_matrix_complex * to) {
    size_t fromc, fromr, toc, tor;
    fromr = from -> size1;
    fromc = from -> size2;
    tor = to -> size1;
    toc = to -> size2;

    if (fromr != tor || fromc != toc) {
        VALGRAPH_ERROR_VOID("matrix dimensions mismathc", VALGRAPH_MBADSIZE);
    } else {
        size_t i,j;
        for(i=0; i<fromr; i++) {
            for(j=0; j<fromc; j++) {
                gsl_matrix_complex_set(to, i, j, gsl_complex_rect(gsl_matrix_get(from, i, j), 0.));
            }
        }
    }
}

// the matrix must be square
gsl_matrix * create_matrix_dyn(size_t n, double * dyn_array) {
    gsl_matrix  * m = gsl_matrix_alloc(n , n);
    size_t i,j;
    for(i=0; i<n; i++) {
        for(j=0; i<n; i++) {
            gsl_matrix_set(m, i, j, dyn_array[i*n + j]);
        }
    }

    return m;
}

gsl_matrix * create_matrix_stat(size_t n, double stat_array[n][n]) {
    gsl_matrix * m = gsl_matrix_alloc(n, n);
    size_t i,j;

    for(i=0; i<n; i++) {
        for(j=0; j<n; j++) {
            gsl_matrix_set(m, i, j, stat_array[i][j]);
        }
    }

    return m;
}

gsl_matrix_complex * create_matrix_complex_dyn(size_t n, double * dyn_array) {
    return from_array_dynamic(n, n, dyn_array);
}

gsl_matrix_complex * create_matrix_complex_stat(size_t n, double stat_array[n][n]) {
    return from_array_static(n, n, stat_array);
}

enum State solve(gsl_matrix * m, gsl_vector * x, double a, double b, double c) {
    /*
     * tranforms (a*I - b*A) X = c*(1,...1) on (I - (b/a)*A)X = (c/a) * (1,...,1)
     * and then solve
     *
     * UPDATE:
     * we solve for ((a/b)*I - A)X = (c/b) * (1,...., 1)
     */
    size_t mr, mc, xs;
    mr = m -> size1;
    mc = m -> size2;
    xs = x -> size;

    if (mr != mc) {
        VALGRAPH_ERROR_VAL("matrix must be square", VALGRAPH_MNOTSQR, 0);
    } else if (mr != xs) {
        VALGRAPH_ERROR_VAL("matrix and vector must have same dimension", VALGRAPH_VBADLEN, 0);
    } else if (a == 0.) {
        VALGRAPH_ERROR_VAL("a parameter must be non zero", VALGRAPH_BADPAR, 0);
    } else if (b == 0.) {
        for(size_t i=0; i<xs; i++) {
            gsl_vector_set(x, i, c / a);

            return SUCCESS;
        }
    } else {
        size_t i,j; // the only index that I need

        // transform your parameters to complex numbers
        gsl_complex identity_mod = gsl_complex_rect(a / b, 0.); // matrix modifier
        gsl_complex vector_mod = gsl_complex_rect(c / b, 0.); // vector modifier

        // create complex matrix and complex vector
        gsl_matrix_complex * mcomplex = gsl_matrix_complex_alloc(mr, mc);
        gsl_vector_complex * xcomplex = gsl_vector_complex_alloc(xs);

        // populate the complex matrix
        from_real_to_complex_matrix(m, mcomplex);

        // no need to populate the vector, it will store the solution
        // call compute_full_v
        // TODO: what to do if matrix are not invertible ??
        compute_full_v(mcomplex, xcomplex, identity_mod, vector_mod);

        for(i=0; i<xs; i++) {
            gsl_vector_set(x, i, GSL_REAL(gsl_vector_complex_get(xcomplex, i)));
        }

        // now you can free your gsl stuff, and yours only
        gsl_matrix_complex_free(mcomplex);
        gsl_vector_complex_free(xcomplex);

        return SUCCESS;
    }
}

enum State is_equivalent(gsl_matrix * m, double a, double b, double a2, double b2) {
    size_t mr, mc;
    mr = m -> size1;
    mc = m -> size2;
    if (mr != mc) {
        VALGRAPH_ERROR_VAL("matrix must be square", VALGRAPH_MNOTSQR, 0);
    } else if (a == 0 || a2 == 0) {
        VALGRAPH_ERROR_VAL("a parameter must be non zero", VALGRAPH_BADPAR, 0);
    } else {
        gsl_vector * sol1 = gsl_vector_alloc(mr);
        gsl_vector * sol2 = gsl_vector_alloc(mr);

        // TODO: what to do if matrix are not invertible ??

        enum State state1 = solve(m, sol1, a, b, 1.);
        enum State state2 = solve(m, sol2, a2, b2, 1.);

        if (state1 == SUCCESS && state2 == SUCCESS) {
           size_t i, j;
           double sol11, sol12, sol21, sol22;

           for(i=0; i<(mr-1); i++) {
               for(j=i+1; j<mr; j++) {
                   sol11 = gsl_vector_get(sol1, i);
                   sol12 = gsl_vector_get(sol1, j);
                   sol21 = gsl_vector_get(sol2, i);
                   sol22 = gsl_vector_get(sol2, j);

                   if (sol11 == sol12 && sol21 != sol22) {
                       return UNEQUIVALENT;
                   } else if (sol11 != sol12 && sol21 == sol22) {
                       return UNEQUIVALENT;
                   } else if (((sol11 - sol12) * (sol21 - sol22)) < 0.) {
                       return UNEQUIVALENT;
                   } else {
                       continue;
                   }
               }
           }

           // free the vectors
           gsl_vector_free(sol1);
           gsl_vector_free(sol2);

           // everything is good
           return EQUIVALENT;
        } else {
            // free the vectors
            gsl_vector_free(sol1);
            gsl_vector_free(sol2);

            // TODO: fix this!
            VALGRAPH_ERROR_VAL("solving went wrong", VALGRAPH_FAILURE, 0);
        }
    }
}

vgraph_result find_bound_A(gsl_matrix * m, double b) {
    size_t mc, mr;
    mr = m -> size1;
    mc = m -> size2;

    if (mc != mr) {
        VALGRAPH_ERROR_VAL("Matrix must be square", VALGRAPH_MNOTSQR, BAD_RESULT);
    } else if (b == 0.) {
       // set the default bound to 1
       vgraph_result res = {SUCCESS, 1.};
        return res;
    } else {
        // allocate the complex counterpart
        gsl_matrix_complex * mcomplex = gsl_matrix_complex_alloc(mr, mc);

        // load the matrix in its complex counterpart
        from_real_to_complex_matrix(m, mcomplex);

        double bound = bound_only_computation(mcomplex, VALGRAPH_TOLERANCE, FULL_V);
        vgraph_result res = {SUCCESS, (1. + fabs(b))*bound}; // TODO: is this correct? verify

        return res;
    }
}
vgraph_result find_bound_B(gsl_matrix * m, double a) {
    size_t mc, mr;
    mr = m -> size1;
    mc = m -> size2;

    // TODO: now idea of what I'm doing...
    //     : with the error handlign

    if (mc != mr) {
        VALGRAPH_ERROR_VAL("Matrix must be square", VALGRAPH_MNOTSQR, BAD_RESULT);
    } else if (a == 0) {
        VALGRAPH_ERROR_VAL("'a' parameter must be non zero", VALGRAPH_BADPAR, BAD_RESULT);
    } else {
        // allocate the matrix counterpart
        gsl_matrix_complex * mcomplex = gsl_matrix_complex_alloc(mr, mc);

        // load the matrix in its complex counterpart
        from_real_to_complex_matrix(m, mcomplex);

        // this value can't never be zero, it is initialized to one, and can only increase
        double bound = bound_only_computation(mcomplex, VALGRAPH_TOLERANCE, FULL_V);

        vgraph_result res = { SUCCESS, a/bound}; // TODO: as before, veirfy this

        return res;
    }
}
