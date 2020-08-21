//
// Created by horacio on 30/06/2020.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef VALGRAPHCORE_CSV_READER_H
#define VALGRAPHCORE_CSV_READER_H

double * read_csv(char * filename, char delimiter);
void write_csv(char * filename, double * array, int rows, int cols, char delimiter);

#endif //VALGRAPHCORE_CSV_READER_H
