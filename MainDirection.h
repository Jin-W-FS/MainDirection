#ifndef MAIN_DIRECTION_H
#define MAIN_DIRECTION_H

#include "fixedpoint.h"

#define DATA_MAXLEN 256

typedef struct _DataVec {
	int len;
	int data[DATA_MAXLEN][3];
} DataVec;

#define STVAL_SHIFT 4	// StVals = real * 16
typedef struct _StVals {
	long avg[3];	// * 16
	long cov[3][3]; // * 16
} StVals;

void calcStatistics(const DataVec* dv, StVals* st);
void stcov2fixed(const StVals* st, fixed mat[3][3]);
void eigenSystem(fixed mat[3][3], fixed e[3], fixed v[3][3]);

#endif
