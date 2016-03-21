#ifndef MAIN_DIRECTION_H
#define MAIN_DIRECTION_H

#include "fixedpoint.h"

#define DATA_MAXLEN 256

typedef struct _DataVec {
	int len;
	int data[DATA_MAXLEN][3];
} DataVec;

typedef struct _StVals {
	fixed avg[3];
	fixed cov[3][3];
} StVals;

void calcStatistics(const DataVec* dv, StVals* st);
void eigenSystem(fixed mat[3][3], fixed e[3], fixed v[3][3]);

#endif
