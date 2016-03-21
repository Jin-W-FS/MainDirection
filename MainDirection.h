#ifndef MAIN_DIRECTION_H
#define MAIN_DIRECTION_H

#define REALNUM_USE_FLOAT	// REALNUM_USE_FIXED
#include "realnum.h"

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

#define EPS (int2real(1) / 10)	// 0.1

void calcStatistics(const DataVec* dv, StVals* st);
void stcov2real(const StVals* st, real mat[3][3]);
void eigenSystem(real mat[3][3], real e[3], real v[3][3]);

#endif
