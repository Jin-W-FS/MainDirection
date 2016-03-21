#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "MainDirection.h"

void calcAverage(const DataVec* dv, fixed avg[3])
{
	const int (*p)[3];
	avg[0] = avg[1] = avg[2] = 0;
	for (p = dv->data; p < dv->data + dv->len; ++p) {
		avg[0] += (*p)[0];
		avg[1] += (*p)[1];
		avg[2] += (*p)[2];
	}
	avg[0] = fixdiv(avg[0], dv->len); // == fixdiv(fixed(avg[0]), fixed(dv->len))
	avg[1] = fixdiv(avg[1], dv->len);
	avg[2] = fixdiv(avg[2], dv->len);
}

void calcCovariance(const DataVec* dv, const fixed avg[3], fixed cov[3][3])
{
	const int (*p)[3];
	long len = dv->len;
	fixed x, y, z;
	fixed xx = 0, xy = 0, xz = 0, yy = 0, yz = 0, zz = 0;
	for (p = dv->data; p < dv->data + len; ++p) {
		x = int2fixed((*p)[0]) - avg[0];
		y = int2fixed((*p)[1]) - avg[1];
		z = int2fixed((*p)[2]) - avg[2];
		xx += fixmul(x, x); xy += fixmul(x, y); xz += fixmul(x, z);
		                    yy += fixmul(y, y); yz += fixmul(y, z);
		                                        zz += fixmul(z, z);
	}
	cov[0][0] = xx / len;	// fixed / long -> fixed
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

#define sign(x) ((x < 0) ? -1 : 1)
#define SIZE_MAT (9 * sizeof(fixed))

void eyes(fixed M[3][3])
{
	memset(M, 0, SIZE_MAT);
	M[0][0] = M[1][1] = M[2][2] = int2fixed(1);
}

fixed offsum(fixed M[3][3])
{
	return fixabs(M[0][1]) + fixabs(M[0][2]) + fixabs(M[1][2]);
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

#define EPS (int2fixed(1) / 10)	// 0.1
int givensRotation(fixed M[3][3], fixed J[3][3], fixed R[3][3], int i, int j)
{
	fixed tao, t, c, s;
	// check
	eyes(R);
	if (fixabs(M[i][j]) < EPS) return 0;
	if (i > j) { int t = i; i = j; j = t; }	// assert i < j

	// calculate cos/sin
	tao = fixdiv(M[i][i] - M[j][j], 2 * M[i][j]);
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
void eigenSystem(fixed mat_[3][3], fixed e[3], fixed v[3][3])
{
	int i, j;
	fixed mat[3][3];
	fixed rot[3][3];
	// init
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			mat[i][j] = mat_[i][j];
		}
	}
	eyes(v);
	// loop
	while (offsum(mat) > 5 * EPS) {
		int c = 0;
		c |= givensRotation(mat, v, rot, 0, 1);
		c |= givensRotation(mat, v, rot, 0, 2);
		c |= givensRotation(mat, v, rot, 1, 2);
		if (!c) break;
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
