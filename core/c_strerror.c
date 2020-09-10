//
// Created by horacio on 10/07/2020.
//

#include "c_errno.h"

const char *
valgraph_strerror (const int valgraph_errno) {
    switch (valgraph_errno) {
        case VALGRAPH_SUCCESS:
            return "success";
        case VALGRAPH_FAILURE:
            return "failure";
        case VALGRAPH_PBADDEG:
            return "polynomials of different degree";
        case VALGRAPH_NOMEM:
            return "failed to allocate memory";
        case VALGRAPH_NOIMPL:
            return "method or functionality to be implemented yet";
        case VALGRAPH_OUTOFB:
            return "out of bounds operation";
        case VALGRAPH_VBADLEN:
            return "pvector(s) of different dimension";
        case VALGRAPH_MBADSIZE:
            return "pmatrix(es) of different dimension";
        case VALGRAPH_MNOTSQR:
            return "pmatrix is not square";
        case VALGRAPH_MBADMUL:
            return "mismatched dimension for pmatrix multiplication";
        case VALGRAPH_PCORRP:
            return "badly constructed or corrupted polynomial";
        case VALGRAPH_BADPAR:
            return "bad parameter value input";
        default:
            return "unknown error code";
    }
}
