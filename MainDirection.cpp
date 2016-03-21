#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "MainDirection.h"

/// sum up statistics
#define lsftst(x)	((x) << STVAL_SHIFT)
#define rsftst(x)	((x) >> STVAL_SHIFT)
void calcAverage(const DataVec* dv, long avg[3])
{
	const int (*p)[3];
	avg[0] = avg[1] = avg[2] = 0;
	for (p = dv->data; p < dv->data + dv->len; ++p) {
		avg[0] += (*p)[0];
		avg[1] += (*p)[1];
		avg[2] += (*p)[2];
	}
	avg[0] = lsftst(avg[0]) / dv->len;
	avg[1] = lsftst(avg[1]) / dv->len;
	avg[2] = lsftst(avg[2]) / dv->len;
}

void calcCovariance(const DataVec* dv, const long avg[3], long cov[3][3])
{
	const int (*p)[3];
	long len = dv->len;
	long x, y, z;
	long xx = 0, xy = 0, xz = 0, yy = 0, yz = 0, zz = 0;
	for (p = dv->data; p < dv->data + len; ++p) {
		x = lsftst((*p)[0]) - avg[0];
		y = lsftst((*p)[1]) - avg[1];
		z = lsftst((*p)[2]) - avg[2];
		xx += rsftst(x * x); xy += rsftst(x * y); xz += rsftst(x * z);
							 yy += rsftst(y * y); yz += rsftst(y * z);
												  zz += rsftst(z * z);
	}
	cov[0][0] = xx / len;
	cov[0][1] = cov[1][0] = xy / len;
	cov[0][2] = cov[2][0] = xz / len;
	cov[1][1] = yy / len;
	cov[1][2] = cov[2][1] = yz / len;
	cov[2][2] = zz / len;
}

void calcStatistics(const DataVec* dv, StVals* st)
{
	calcAverage(dv, st->avg);
	calcCovariance(dv, st->avg, st->cov);
}

/// covert
void stcov2fixed(const StVals* st, fixed mat[3][3])
{
	int i;
	const long* p;
	fixed* q;
	for (i = 0, p = st->cov[0], q = mat[0]; i < 9; i++, p++, q++) {
#if FIXEDLEN >= STVAL_SHIFT
		(*q) = (*p) << (FIXEDLEN - STVAL_SHIFT);
#else
		(*q) = (*p) >> (STVAL_SHIFT - FIXEDLEN);
#endif
	}
}

/// fixedpoint math
#define sign(x) ((x < 0) ? -1 : 1)
#define SIZE_MAT (9 * sizeof(fixed))

void eyes(fixed M[3][3])
{
	memset(M, 0, SIZE_MAT);
	M[0][0] = M[1][1] = M[2][2] = int2fixed(1);
}

fixed offmax(fixed M[3][3], int* idx)
{
	fixed v[3] = { fixabs(M[0][1]), fixabs(M[1][2]), fixabs(M[2][0]) };
	*idx = 0;
	if (v[1] > v[*idx]) *idx = 1;
	if (v[2] > v[*idx]) *idx = 2;
	return v[*idx];
}

// M = M1 * M2, safe if M == M1 or M == M2
void matMul(fixed M[3][3], fixed M1[3][3], fixed M2[3][3])
{
	fixed O[3][3];
	int i, j;
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			O[i][j] = fixmul(M1[i][0], M2[0][j]) + fixmul(M1[i][1], M2[1][j]) + fixmul(M1[i][2], M2[2][j]);
		}
	}
	memcpy(M, O, SIZE_MAT);
}

#include <stdio.h>
#define EPS (int2fixed(1) / 10)	// 0.1
int givensRotation(fixed M[3][3], fixed J[3][3], fixed R[3][3], int i, int j)
{
	//printf("M:" "\t{{ %.3f, %.3f, %.3f },\n"
	//		"\t { %.3f, %.3f, %.3f },\n"
	//		"\t { %.3f, %.3f, %.3f }}\n\n", 
	//		fixed2float(M[0][0]), fixed2float(M[0][1]), fixed2float(M[0][2]),
	//		fixed2float(M[1][0]), fixed2float(M[1][1]), fixed2float(M[1][2]),
	//		fixed2float(M[2][0]), fixed2float(M[2][1]), fixed2float(M[2][2]));
	fixed tao, t, c, s;
	// check
	eyes(R);
	if (fixabs(M[i][j]) < EPS) return 0;
	if (i > j) { int t = i; i = j; j = t; }	// assert i < j

	// calculate cos/sin
	tao = fixdiv((M[i][i] - M[j][j]) / 2, M[i][j]);
	t = fixdiv(sign(tao) * int2fixed(1), fixabs(tao) + fixsqrt(int2fixed(1) + fixmul(tao, tao)));	// tan
	c = fixdiv(int2fixed(1), fixsqrt(int2fixed(1) + fixmul(t, t)));	// cos
	s = fixmul(c, t);				// sin
	if (s == 0) return 0;

	// rotate matrix
	R[i][i] = c; R[i][j] = -s;
	R[j][i] = s; R[j][j] = c;
	matMul(M, M, R);	// M = M * R'
	R[i][j] = s; R[j][i] = -s;
	matMul(M, R, M);	// M = R * M	
	matMul(J, R, J);	// J = R * J
	return 1;
}

// sort by eigen value(e) (ascending)
void _swap_eigens(fixed e[3], fixed v[3][3], int i, int j) {
	fixed temp;
	fixed buff[3];
	temp = e[i];
	e[i] = e[j];
	e[j] = temp;
	memcpy(buff, v[i], sizeof(buff));
	memcpy(v[i], v[j], sizeof(buff));
	memcpy(v[j], buff, sizeof(buff));
}

void sortEigens(fixed e[3], fixed v[3][3])
{
	int s = 0;	// swap times (even or odd)
	int i = 0;	// idx of min(e)
	if (e[1] < e[0]) i = 1;
	if (e[2] < e[i]) i = 2;
	if (i != 0) { _swap_eigens(e, v, 0, i); s = !s; }
	if (e[2] < e[1]) { _swap_eigens(e, v, 1, 2); s = !s; }
	if (s) { for (i = 0; i < 3; i++) v[2][i] = -v[2][i]; }	// keep Vz = Vx \cross Vy
}

// calculate eigen value and eigen vector
void eigenSystem(fixed mat[3][3], fixed e[3], fixed v[3][3])
{
	int i;
	fixed rot[3][3];
	// init
	eyes(v);
	// loop
	int idx;
	while (offmax(mat, &idx) > EPS) {
		if (!givensRotation(mat, v, rot, idx, (idx+1)%3)) break;
	}
	// get result
	for (i = 0; i < 3; i++) {
		e[i] = mat[i][i];
	}
	sortEigens(e, v);
}

#ifdef TEST
#include <stdio.h>
int main()
{
	float M1[3][3] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	float M2[3][3] = { 4, 5, 6, 1, 2, 3, 7, 8, 9 };
	matMul(M1, M1, M2);
	printf(	"{{ %f, %f, %f },\n"
			" { %f, %f, %f },\n"
			" { %f, %f, %f }}\n",
			M1[0][0], M1[0][1], M1[0][2],
			M1[1][0], M1[1][1], M1[1][2],
			M1[2][0], M1[2][1], M1[2][2]);
	system("pause");
	return 0; 
}
#endif
