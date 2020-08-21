//
// Created by horacio on 10/07/2020.
//

#include <stdio.h>
#include <errno.h>

#ifndef VALGRAPHCORE_C_ERRNO_H
#define VALGRAPHCORE_C_ERRNO_H

enum {
    VALGRAPH_SUCCESS = 0,
    VALGRAPH_FAILURE = 1,
    VALGRAPH_PBADDEG = 2, // polynomial bad length
    VALGRAPH_NOMEM = 3, // malloc failed
    VALGRAPH_NOIMPL = 4, // method not implemented yet
    VALGRAPH_OUTOFB = 5, // out of bound error
    VALGRAPH_VBADLEN = 6, // pvector of mismatched length
    VALGRAPH_MBADSIZE = 7, // pmatrix of mismatched length
    VALGRAPH_MNOTSQR = 8, // pmatrix not square
    VALGRAPH_MBADMUL = 9, // pmatrix can't be multiplied
    VALGRAPH_PCORRP = 10, // badly constructed polynomial
};


void valgraph_error (const char * reason, const char * file, int line,
                int valgraph_errno);

void valgraph_stream_printf (const char *label, const char *file,
                        int line, const char *reason);

const char * valgraph_strerror (const int valgraph_errno);

typedef void valgraph_error_handler_t (const char * reason, const char * file,
                                  int line, int valgraph_errno);

typedef void valgraph_stream_handler_t (const char * label, const char * file,
                                   int line, const char * reason);

valgraph_error_handler_t *
valgraph_set_error_handler (valgraph_error_handler_t * new_handler);

valgraph_error_handler_t *
valgraph_set_error_handler_off (void);

valgraph_stream_handler_t *
valgraph_set_stream_handler (valgraph_stream_handler_t * new_handler);

FILE * valgraph_set_stream (FILE * new_stream);

/* VALGRAPH_ERROR: call the error handler, and return the error code */

#define VALGRAPH_ERROR(reason, valgraph_errno) \
       do { \
       valgraph_error (reason, __FILE__, __LINE__, valgraph_errno) ; \
       return valgraph_errno ; \
       } while (0)

/* VALGRAPH_ERROR_VAL: call the error handler, and return the given value */

#define VALGRAPH_ERROR_VAL(reason, valgraph_errno, value) \
       do { \
       valgraph_error (reason, __FILE__, __LINE__, valgraph_errno) ; \
       return value ; \
       } while (0)

/* valgraph_ERROR_VOID: call the error handler, and then return
   (for void functions which still need to generate an error) */

#define VALGRAPH_ERROR_VOID(reason, valgraph_errno) \
       do { \
       valgraph_error (reason, __FILE__, __LINE__, valgraph_errno) ; \
       return ; \
       } while (0)

/* VALGRAPH_ERROR_NULL suitable for out-of-memory conditions */

#define VALGRAPH_ERROR_NULL(reason, valgraph_errno) VALGRAPH_ERROR_VAL(reason, valgraph_errno, 0)

#endif //VALGRAPHCORE_C_ERRNO_H
