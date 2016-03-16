#ifndef MAIN_DIRECTION_H
#define MAIN_DIRECTION_H

#define DATA_MAXLEN 200

typedef struct _StVals {
	long avg[3];
	long cov[3][3];
} StVals;

typedef struct _DataVec {
	int len;
	int data[DATA_MAXLEN][3];
} DataVec;

void calcStatistics(const DataVec* dv, StVals* st);
void eigenSystem(long mat[3][3], float e[3], float v[3][3]);

#endif
