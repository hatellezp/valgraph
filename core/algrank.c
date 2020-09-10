//
// Created by horacio on 29/06/2020.
//

#include "algrank.h"
#include "c_errno.h"

gsl_complex compute_v(gsl_matrix_complex * m, gsl_complex mmod, gsl_complex vmod, size_t ind) {
    size_t r, c;
    r = m -> size1;
    c = m -> size2;

    if (r != c) { // we need a square matrix
        VALGRAPH_ERROR_VAL("matrix must be square", VALGRAPH_MNOTSQR, gsl_complex_rect(0., 0.));
    } else {
        // allocate for matrices
        gsl_matrix_complex * mcopy = gsl_matrix_complex_alloc(r, c);
        gsl_matrix_complex * miden = gsl_matrix_complex_alloc(r, c);

        // set identity and copy to the other matrix
        gsl_matrix_complex_set_identity(miden) ;
        gsl_matrix_complex_memcpy(mcopy, m);

        // multiply iden by v and add mcopy
        gsl_matrix_complex_scale(mcopy, mmod);
        gsl_matrix_complex_sub(miden, mcopy);

        // change the i-columnt with ones
        gsl_vector_complex * v = gsl_vector_complex_ones(r);
        gsl_vector_complex_scale(v, vmod);

        gsl_matrix_complex_chcol(mcopy, v, ind);


        // compute the determinant
        gsl_permutation * p = gsl_permutation_alloc(r); // same size as c
        int signum;
        gsl_linalg_complex_LU_decomp(mcopy, p, &signum);
        gsl_complex det = gsl_linalg_complex_LU_det(mcopy, signum);

        // free the memory
        gsl_matrix_complex_free(mcopy);
        gsl_matrix_complex_free(miden);
        gsl_vector_complex_free(v);
        gsl_permutation_free(p);

        // return the value
        return det;
    }
}

// this function copies the matrix m and erases from the copy the diagonal values, setting them to 0 !!!
// this function does the same work as compute_v but for a full vector this time
void compute_full_v(gsl_matrix_complex * m, gsl_vector_complex * x, gsl_complex identity_mod, gsl_complex vector_mod) {
    size_t r, c, xs;
    r = m -> size1;
    c = m -> size2;
    xs = x -> size;

    if (r != c) { // we need a square matrix
        VALGRAPH_ERROR_VOID("The matrix must be square", VALGRAPH_MNOTSQR);
    } else if (xs != r) {
        VALGRAPH_ERROR_VOID("the solution vector must have same dimension as the matrix", VALGRAPH_VBADLEN);
    } else {
        // allocate for matrices
        gsl_matrix_complex * mcopy = gsl_matrix_complex_alloc(r, c);
        gsl_matrix_complex * miden = gsl_matrix_complex_alloc(r, c);

        // set identity and copy to the other matrix
        gsl_matrix_complex_set_identity(miden) ;
        gsl_matrix_complex_memcpy(mcopy, m);

        // set the diagonal of mcopy to zero, as a measure of security
        for(size_t i=0; i<r; i++) {
            gsl_matrix_complex_set(mcopy, i, i, gsl_complex_rect(0., 0.));
        }

        // UPDATE: the new computation is: (identity_mod*1 - A)
        gsl_matrix_complex_scale(miden, identity_mod); // perform identity_mod * I
        gsl_matrix_complex_sub(miden, mcopy); // perform identity_mod*1 - A

        // BE CAREFUL: the important matrix lives in miden now

        // create right-hand vector
        gsl_vector_complex * v = gsl_vector_complex_ones(r);
        // scale by vmod
        gsl_vector_complex_scale(v, vector_mod); // perform vmod*(1,...,1)

        // no determinant needed, compute directly the solution to (identity_mod*1 - A)X = vector_mod*(1,...,1)
        gsl_permutation * p = gsl_permutation_alloc(r); // same size as c
        int signum;
        gsl_linalg_complex_LU_decomp(miden, p, &signum);

        gsl_linalg_complex_LU_solve(miden, p, v, x); // here we retrieve the solution and store it in x

        // free the memory
        gsl_matrix_complex_free(mcopy);
        gsl_matrix_complex_free(miden);
        gsl_vector_complex_free(v);
        gsl_permutation_free(p);

        // no need of return
        // the solution is stored in x
    }
}

void create_unity_roots(gsl_vector_complex * units, size_t n, int inverse) { // I don't like the inverse arg as an int
    double t = (inverse) ? -1. : 1.; // be sure of this
    size_t k, usize;
    usize = units -> size;
    if (usize != n) {
        VALGRAPH_ERROR_VOID("vector length and number of roots are different", VALGRAPH_VBADLEN);
    } else {
        for(k=0; k<n; k++) {
            double arg = t*((2.0*M_PI*k) / ((double)n));
            double real_part = cos(arg);
            double imag_part = sin(arg);
            gsl_complex z = gsl_complex_rect(real_part, imag_part);
            gsl_vector_complex_set(units, k, z);
        }
    }
}

double bound_only_computation(gsl_matrix_complex * m, double tolerance, int method) {
    size_t r,c;
    r = m -> size1;
    c = m -> size2;

    if (r != c) { // verify is an square matrix
        printf("only square matrix!\n");
        VALGRAPH_ERROR_VAL("matrix must be square", VALGRAPH_MNOTSQR, 0);
    } else {
        size_t N = r - 1; // we need only r-1 points, because we know the polynomial has at most degree r-2
        size_t i, indv, indp, indpi, indpj;

        // initialize the fftw necessary parts
        fftw_complex * in, * out;
        fftw_plan p;
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_MEASURE); // TODO: test FFTW_ESTIMATE too

        // get unity roots
        gsl_vector_complex * rts = gsl_vector_complex_alloc(N);
        create_unity_roots(rts, N, UR_BACKWARDS);

        // create the matrix that will keep the values
        // this keeps array p_i(r) where i is the number of the corresponding polynomial and r is a root of the unity
        // allocate the memory to a 1D array
        double * pvalues = malloc(sizeof(double)*(N+1)*N*2); // number of poly times number of values times 2

        // the two different methods
        // the first computing a value at a time, the second a full vector for each loop
        if (method == SINGLE_V) { // single value method
            // populate the array

            /*
            for(indp=0; indp<N+1; indp++) { // indp is the polynomial index
                for(indv=0; indv<N; indv++) { // indv is the value index (the unity root)
                    gsl_complex root = gsl_vector_complex_get(rts, indv);
                    gsl_complex value = compute_v(m, root, indp); // computes det((root*id + m)|indp)

                    // now put the value in the array
                    pvalues[(N+1)*indp + N*indv] = GSL_REAL(value);
                    pvalues[(N+1)*indp + N*indv + 1] = GSL_IMAG(value);
                }
            }
             */
        } else if (method == FULL_V) { // full vector method
            gsl_vector_complex * values = gsl_vector_complex_alloc(N); // create vector
            for(indv=0; indv<N; indv++) {
                gsl_complex root = gsl_vector_complex_get(rts, indv);
                compute_full_v(m, values, root, gsl_complex_rect(1., 0.)); // here vmod should be (1,0) and mmod should be the root

                for(indp=0; indp<N+1; indp++) {
                    gsl_complex z = gsl_vector_complex_get(values, indp);
                    pvalues[(N+1)*indp + N*indv] = GSL_REAL(z);
                    pvalues[(N+1)*indp + N*indv + 1] = GSL_IMAG(z);
                }
            }
            // free the vector
            gsl_vector_complex_free(values);
        }

        // this values are common for all internal loops
        double current_max = 1.; // I prefer one as current max, this way we avoid some negative unpleasant errors
        size_t real_degree;
        double possible_coeff_real, max_coeff;
        double Ndouble = (double)N;

        // from the pvalues array populate the in array and compute the fft inverse
        for(indpi=0; indpi<N; indpi++) { // be careful with the indexes ! 0 < indpi < indpj < N+1 (number of polys)
            for(indpj=indpi+1; indpj<N+1; indpj++) {
                for(indv=0; indv<N; indv++) {
                    // the in array has the correct values
                    // remember that we are interested in the differences p_i(a) - p_j(a)
                    in[indv][0] = pvalues[(N+1)*indpi + N*indv] - pvalues[(N+1)*indpj + N*indv];
                    in[indv][1] = pvalues[(N+1)*indpi + N*indv + 1] - pvalues[(N+1)*indpj + N*indv + 1];

                    // after execution of p capture the values in the out array
                    fftw_execute(p);

                    // find the real degree of the polynomial
                    // I think these updates are not necessary ...
                    real_degree = -1;
                    possible_coeff_real = tolerance / 2.;
                    max_coeff = tolerance / 2.;

                    // loop to find the real degree
                    for(i=0; i<N; i++) {
                        possible_coeff_real = out[N-1-i][0] / Ndouble; // I only need the real coefficient
                        if (fabs(possible_coeff_real) > tolerance) { // is a credible coefficient
                            real_degree = N-1-i;
                            max_coeff = possible_coeff_real;
                            // get out of loop
                            break;
                        } else { // continue the search
                            continue;
                        }
                    }

                    // now that we have the real degree we compute the bound
                    // where the max comprehend all polynomials
                    for(size_t deg=0; deg<real_degree; deg++) {
                        current_max = fmax(current_max, (-out[deg][0])/max_coeff);
                    }
                }
            }
        }

        // at the end don't forget to free the memory
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
        free(pvalues); // TODO: make yourself sure this is safe
        gsl_vector_complex_free(rts);

        // return the bound value
        current_max = 2. + current_max; // add 2 to be sure be outside the roots
        return current_max;
    }
}

void compute_ranking(gsl_matrix_complex * m, double a, double b, double c, int * error, double * ranking) {
    int rows, cols;
    rows = m -> size1;
    cols = m -> size2;

    if (rows == cols) { // verify is a square matrix

        // allocate matrices and structs for LU solving
        gsl_matrix_complex * mcopy = gsl_matrix_complex_alloc(rows, cols);
        gsl_matrix_complex * miden = gsl_matrix_complex_alloc(rows, cols);
        gsl_permutation * p = gsl_permutation_alloc(rows);
        gsl_vector_complex * v = gsl_vector_complex_ones(rows);
        gsl_vector_complex * sol = gsl_vector_complex_ones(rows);
        int i, signum;

        // put m in the mcopy matrix
        gsl_matrix_complex_memcpy(mcopy, m);

        // create the a*id - b*m matrix
        gsl_matrix_complex_set_identity(miden);
        gsl_matrix_complex_scale(miden, gsl_complex_rect(a, 0.));

        gsl_matrix_complex_scale(mcopy, gsl_complex_rect(b, 0.));

        // this is really important, we work with (a*Id - b*M)X = c
        gsl_matrix_complex_sub(miden, mcopy);
        // gsl_matrix_complex_add(miden, mcopy); TODO: there is also this possibility

        // don't forget to create the vector c and also vector solution
        gsl_vector_complex_scale(v, gsl_complex_rect(c, 0.));

        gsl_linalg_complex_LU_decomp(miden, p, &signum);

        // know if you can solve the system
        gsl_complex det = gsl_linalg_complex_LU_det(miden, signum);
        if (GSL_REAL(det) != 0.) {
            gsl_linalg_complex_LU_solve(miden, p, v, sol);

            // take the solution to the ranking vector
            for(i=0; i<rows; i++) {
                ranking[i] = GSL_REAL(gsl_vector_complex_get(sol, i));
            }

            (* error) = 0; // no error
        } else {
            (* error) = 1; // an error is sent
        }

        // free all the memory
        gsl_matrix_complex_free(miden);
        gsl_matrix_complex_free(mcopy);
        gsl_permutation_free(p);
        gsl_vector_complex_free(v);
        gsl_vector_complex_free(sol);
    }

}

void compute_ranking_stable(gsl_matrix_complex * m, double tolerance, int method, int * error, double * ranking) {
    double bound = bound_only_computation(m, tolerance, method);
    int n = m -> size1;

    compute_ranking(m, bound, 1., 1., error, ranking);
}

/*
 * thanks to 'Numerical Recipes in C' for the O(n^2) algorithm
 */
void polynomial_coefficient_computation(double * x, double * y, double * cof, size_t n) {
    int k, j, i;
    double phi, ff, b;
    double * s = malloc(sizeof(double)*n);

    for(i=0; i<n; i++) s[i] = cof[i] = 0.;

    s[n-1] = -x[0];

    for(i=1; i<n; i++) {
        for(j=n-i; j<n-1; j++) {
            s[j] -= x[i] * s[j+1];
        }
        s[n-1] -= x[i];
    }

    for(j=0; j<n; j++) {
        phi = n; // or n+1 ??
        for(k=n-1; k>=1; k--) {
            phi = ((double)k) * s[k] + x[j]*phi;
        }
        ff = y[j] / phi;
        b = 1.0;
        for(k=n-1; k>=0; k--) {
            cof[k] += b*ff;
            b = s[k] + x[j]*b;
        }
    }
    free(s);
}
