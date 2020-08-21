//
// Created by horacio on 10/07/2020.
//

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include "c_errno.h"

valgraph_error_handler_t * valgraph_error_handler = NULL;

static void no_error_handler (const char *reason, const char *file, int line, int valgraph_errno);

void
valgraph_error (const char * reason, const char * file, int line, int valgraph_errno)
{
    if (valgraph_error_handler)
    {
        (*valgraph_error_handler) (reason, file, line, valgraph_errno);
        return ;
    }

    valgraph_stream_printf ("ERROR", file, line, reason);

    fflush (stdout);
    fprintf (stderr, "Default VALGRAPH_CORE error handler invoked.\n");
    fflush (stderr);

    abort ();
}

valgraph_error_handler_t *
valgraph_set_error_handler (valgraph_error_handler_t * new_handler)
{
    valgraph_error_handler_t * previous_handler = valgraph_error_handler;
    valgraph_error_handler = new_handler;
    return previous_handler;
}


valgraph_error_handler_t *
valgraph_set_error_handler_off (void)
{
    valgraph_error_handler_t * previous_handler = valgraph_error_handler;
    valgraph_error_handler = no_error_handler;
    return previous_handler;
}

static void
no_error_handler (const char *reason, const char *file, int line, int valgraph_errno)
{
    /* do nothing */
    reason = 0;
    file = 0;
    line = 0;
    valgraph_errno = 0;
    return;
}




