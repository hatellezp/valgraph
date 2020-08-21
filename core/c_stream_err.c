//
// Created by horacio on 10/07/2020.
//

#include <stdio.h>

#include "c_errno.h"

FILE * valgraph_stream = NULL ;
valgraph_stream_handler_t * valgraph_stream_handler = NULL;

void
valgraph_stream_printf (const char *label, const char *file, int line,
                   const char *reason)
{
    if (valgraph_stream == NULL)
    {
        valgraph_stream = stderr;
    }
    if (valgraph_stream_handler)
    {
        (*valgraph_stream_handler) (label, file, line, reason);
        return;
    }
    fprintf (valgraph_stream, "core: %s:%d: %s: %s\n", file, line, label, reason);

}

valgraph_stream_handler_t *
valgraph_set_stream_handler (valgraph_stream_handler_t * new_handler)
{
    valgraph_stream_handler_t * previous_handler = valgraph_stream_handler;
    valgraph_stream_handler = new_handler;
    return previous_handler;
}

FILE *
valgraph_set_stream (FILE * new_stream)
{
    FILE * previous_stream;
    if (valgraph_stream == NULL) {
        valgraph_stream = stderr;
    }
    previous_stream = valgraph_stream;
    valgraph_stream = new_stream;
    return previous_stream;
}
