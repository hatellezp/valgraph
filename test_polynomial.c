//
// Created by horacio on 13/07/2020.
//

#include <stdio.h>

#include "core/polynomial.h"
#include <gsl/gsl_vector.h>

void test_allocation() {
    printf("---------------------------------------\n");
    printf("--test allocation and freed of memory--\n");
    printf("---------------------------------------\n");
    printf("\n");

    printf("allocating polynomial of zero degree\n");
    vgraph_poly * pzero = vgraph_poly_alloc(0);
    printf("-success\n");

    printf("freeing polynomial of zero degree\n");
    vgraph_poly_free(pzero);
    printf("-success\n");

    printf("allocating polynomial of degree one\n");
    vgraph_poly * pone = vgraph_poly_alloc(1);
    printf("-success\n");

    printf("freeing polynomial of degree one\n");
    vgraph_poly_free(pone);
    printf("-success\n");


    printf("allocating polynomial of degree two\n");
    vgraph_poly * ptwo = vgraph_poly_alloc(2);
    printf("-success\n");

    printf("freeing polynomial of degree two\n");
    vgraph_poly_free(ptwo);
    printf("-success\n");

    printf("allocating two polynomials at the same time");
    vgraph_poly * p1 = vgraph_poly_alloc(3);
    vgraph_poly * p2 = vgraph_poly_alloc(4);
    printf("-success\n");

    printf("freeing both polynomials\n");
    vgraph_poly_free(p1);
    vgraph_poly_free(p2);
    printf("-success\n");

    printf("\n");
    printf("---------------------------------------\n");
    printf("------all allocation test success------\n");
    printf("-----------------------------------------");
}

void test_print() {
    printf("-----------------------------------------\n");
    printf("--testing pretty printer of polynomials--\n");
    printf("-----------------------------------------\n");
    printf("\n");

    vgraph_poly * p0 = vgraph_poly_alloc(0);
    vgraph_poly * p1 = vgraph_poly_alloc(1);
    vgraph_poly * p2 = vgraph_poly_alloc(2);

    printf("should print: 0.\n");
    printf("prints      : ");
    vgraph_poly_print(p0);
    printf("\n");

    gsl_vector_set((p0 -> coeffs), 0, 1.);

    printf("\n");
    printf("should print: 1.000000\n");
    printf("prints      : ");
    vgraph_poly_print(p0);
    printf("\n");

    gsl_vector_set((p1 -> coeffs), 1, 2.);

    printf("\n");
    printf("should print: 2.000000*x\n");
    printf("prints      : ");
    vgraph_poly_print(p1);
    printf("\n");

    gsl_vector_set((p2 -> coeffs), 2, 3.);
    gsl_vector_set((p2 -> coeffs), 0, 1.);
    p2 -> real_degree = 2;

    printf("\n");
    printf("should print: 3.000000*x^2 + 1.000000\n");
    printf("prints      : ");
    vgraph_poly_print(p2);
    printf("\n");

    vgraph_poly_free(p0); vgraph_poly_free(p1); vgraph_poly_free(p2);

    printf("\n");
    printf("-----------------------------------------\n");
    printf("-----user should evaluate the output-----\n");
    printf("-----------------------------------------\n");
}

void test_set_get() {
    printf("-----------------------------------------\n");
    printf("----testing setter and getter methods----\n");
    printf("-----------------------------------------\n");
    printf("\n");

    vgraph_poly * p = vgraph_poly_alloc(5);
    vgraph_poly * p1 = vgraph_poly_alloc(5);

    int s1, s2, s3, s4, s5, g1, g2;

    vgraph_poly_set(p, 0, 5.);
    if (gsl_vector_get((p -> coeffs), 0) == 5.) {
        printf("set test 1 success\n");
        s1 = 1;
    } else {
        printf("set test 1 failure\n");
        s1 = 0;
    }

    vgraph_poly_set(p, 5, 3.);
    if (gsl_vector_get((p -> coeffs), 5) == 3.) {
        printf("set test 2 success\n");
        s2 = 1;
    } else {
        printf("set test 2 failure\n");
        s2 = 0;
    }

    double val = vgraph_poly_get(p, 0);
    if (gsl_vector_get((p -> coeffs), 0) == val) {
        printf("get test 1 success\n");
        g1 = 1;
    } else {
        printf("get test 1 failure\n");
        g1 = 0;
    }

    val = vgraph_poly_get(p, 5);
    if (gsl_vector_get((p -> coeffs), 5) == val) {
        printf("get test 2 success\n");
        g2 = 1;
    } else {
        printf("get test 2 failure\n");
        g2 = 0;
    }

    vgraph_poly_set(p1, 2, 2.);
    vgraph_poly_set(p1, 4, 2.);

    if (p1 -> real_degree == 4) {
        s3 = 1;
        printf("set test 3 success\n");
    } else {
        s3 = 0;
        printf("set test 3 failure\n");
    }

    vgraph_poly_set(p1, 4, 0.);

    if (p1 -> real_degree == 2) {
        s4 = 1;
        printf("set test 4 success\n");
    } else {
        s4 = 0;
        printf("set test 5 failure\n");
    }

    vgraph_poly_set(p1, 2, 0.);

    if (p1 -> real_degree == 0) {
        s5 = 1;
        printf("set test 5 success\n");
    } else {
        s5 = 0;
        printf("set test 5 failure\n");
    }

    if (s1 * s2 * s3 * s4 * s5 * g1 * g2) {
        printf("\n");
        printf("-----------------------------------------\n");
        printf("--------all modifier test success--------\n");
        printf("-----------------------------------------\n");
    } else {
        printf("\n");
        printf("-----------------------------------------\n");
        printf("----!!! some modifier tests failed !!!---\n");
        printf("-----------------------------------------\n");
    }
}

void test_equality() {
    printf("-----------------------------------------\n");
    printf("---------testing equality method---------\n");
    printf("-----------------------------------------\n");
    printf("\n");

    vgraph_poly * p0 = vgraph_poly_alloc(0);
    vgraph_poly * p1 = vgraph_poly_alloc(0);
    vgraph_poly * p2 = vgraph_poly_alloc(1);
    vgraph_poly * p3 = vgraph_poly_alloc(1);
    vgraph_poly * p4 = vgraph_poly_alloc(3);
    vgraph_poly * p5 = vgraph_poly_alloc(3);
    vgraph_poly * p6 = vgraph_poly_alloc(5);

    int * t = malloc(sizeof(int) * 10);

    if (vgraph_poly_equal(p0, p1)) {
        t[0] = 1;
        printf("1st  equality test success\n");
    } else {
        t[0] = 0;
        printf("1st  equality test failure\n");
    }

    vgraph_poly_set(p0, 0, 2.);

    if (!vgraph_poly_equal(p0, p1)) {
        t[1] = 1;
        printf("2nd  equality test success\n");
    } else {
        t[1] = 0;
        printf("2nd  equality test failure\n");
    }

    if (!vgraph_poly_equal(p0, p2)) {
        t[2] = 1;
        printf("3rd  equality test success\n");
    } else {
        t[2] = 0;
        printf("3rd  equality test failure\n");
    }

    vgraph_poly_set(p1, 0, 1.);

    if (!vgraph_poly_equal(p1, p4)) {
        t[3] = 1;
        printf("4th  equality test success\n");
    } else {
        t[3] = 0;
        printf("4th  equality test failure\n");
    }

    if (vgraph_poly_equal(p2, p3)) {
        t[4] = 1;
        printf("5th  equality test success\n");
    } else {
        t[4] = 0;
        printf("5th  equality test failure\n");
    }

    vgraph_poly_set(p2, 1, 4.);

    if (!vgraph_poly_equal(p2, p3)) {
        t[5] = 1;
        printf("6th  equality test success\n");
    } else {
        t[5] = 0;
        printf("6th  equality test failure\n");
    }

    if (vgraph_poly_equal(p4, p5)) {
        t[6] = 1;
        printf("7th  equality test success\n");
    } else {
        t[6] = 0;
        printf("7th  equality test failure\n");
    }

    vgraph_poly_set(p4, 1, 4.);

    if (!vgraph_poly_equal(p2, p4)) {
        t[7] = 1;
        printf("8ht  equality test success\n");
    } else {
        t[7] = 0;
        printf("8ht  equality test failure\n");
    }

    vgraph_poly_set(p5, 1, 4.);

    if (vgraph_poly_equal(p4, p5)) {
        t[8] = 1;
        printf("9th  equality test success\n");
    } else {
        t[8] = 0;
        printf("9th  equality test failure\n");
    }


    if (!vgraph_poly_equal(p5, p6)) {
        t[9] = 1;
        printf("10th equality test success\n");
    } else {
        t[9] = 0;
        printf("10th equality test failure\n");
    }

    int tt=1;
    for(int i=0; i<10; i++) {
       tt *= t[i];
    }

    free(t);
    vgraph_poly_free(p0);
    vgraph_poly_free(p1);
    vgraph_poly_free(p2);
    vgraph_poly_free(p3);
    vgraph_poly_free(p4);
    vgraph_poly_free(p5);
    vgraph_poly_free(p6);

    if (tt) {
        printf("\n");
        printf("-----------------------------------------\n");
        printf("--------all equality test success--------\n");
        printf("-----------------------------------------\n");
    } else {
        printf("\n");
        printf("-----------------------------------------\n");
        printf("--------some equality test failed--------\n");
        printf("-----------------------------------------\n");
    }
}

void test_creation() {
    printf("-----------------------------------------\n");
    printf("--------testing creation methods---------\n");
    printf("-----------------------------------------\n");
    printf("\n");

    printf("testing creating polynomial from array\n");
    double * ar = malloc(5);
    ar[0] = ar[3] = 0.;
    ar[1] = ar[2] = 3.;
    ar[4] = 2.;

    vgraph_poly * p0 = vgraph_poly_alloc(4);
    vgraph_poly * p1 = vgraph_poly_from_array(ar, 5);

    vgraph_poly_set(p0, 1, 3.);
    vgraph_poly_set(p0, 2, 3.);
    vgraph_poly_set(p0, 4, 2.);

    if (vgraph_poly_equal(p0, p1)) {
        printf("first test success\n");
    }


    vgraph_poly_free(p0);
    vgraph_poly_free(p1);
    free(ar);
}

int main() {

    // test allocation
    test_allocation();

    printf("\n");

    // test pretty print
    test_print();

    printf("\n");

    // test setter and getter methods
    test_set_get();

    printf("\n");

    // test equality method
    test_equality();

    printf("\n");

    // test creation method
    test_creation();

    return 0;
}