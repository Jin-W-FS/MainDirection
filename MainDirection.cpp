#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "MainDirection.h"

static inline long qround(long x, long y)
{
	ldiv_t qr = ldiv(x, y);
	return qr.quot + (2 * qr.rem >= y);
}

void calcAverage(const DataVec* dv, long avg[3])
{
	const int (*p)[3];
	avg[0] = avg[1] = avg[2] = 0;
	for (p = dv->data; p < dv->data + dv->len; ++p) {
		avg[0] += (*p)[0];
		avg[1] += (*p)[1];
		avg[2] += (*p)[2];
	}
	avg[0] = qround(avg[0], dv->len);
	avg[1] = qround(avg[1], dv->len);
	avg[2] = qround(avg[2], dv->len);
}

void calcCovariance(const DataVec* dv, const long avg[3], long cov[3][3])
{
	const int (*p)[3];
	long len = dv->len;
	long x, y, z;
	long xx = 0, xy = 0, xz = 0, yy = 0, yz = 0, zz = 0;
	for (p = dv->data; p < dv->data + len; ++p) {
		x = (*p)[0] - avg[0];
		y = (*p)[1] - avg[1];
		z = (*p)[2] - avg[2];
		xx += x * x; xy += x * y; xz += x * z;
				     yy += y * y; yz += y * z;
		                          zz += z * z;
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

#define EPS 0.001f
#define sign(x) ((x < 0) ? -1 : 1)
#define SIZE_MAT (9 * sizeof(float))

void eyes(float M[3][3])
{
	memset(M, 0, SIZE_MAT);
	M[0][0] = M[1][1] = M[2][2] = 1;
}

float offsum(float M[3][3])
{
	return fabsf(M[0][1]) + fabsf(M[0][2]) + fabsf(M[1][2]);
}

// M = M1 * M2, safe if M == M1 or M == M2
void matMul(float M[3][3], float M1[3][3], float M2[3][3])
{
	float O[3][3];
	int i, j;
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			O[i][j] = M1[i][0] * M2[0][j] + M1[i][1] * M2[1][j] + M1[i][2] * M2[2][j];
		}
	}
	memcpy(M, O, SIZE_MAT);
}

int givensRotation(float M[3][3], float J[3][3], float R[3][3], int i, int j)
{
	float tao, t, c, s;
	// check
	eyes(R);
	if (fabsf(M[i][j]) < EPS) return 0;
	if (i > j) { int t = i; i = j; j = t; }	// assert i < j

	// calculate cos/sin
	tao = (M[i][i] - M[j][j]) / (2 * M[i][j]);
	t = sign(tao) / (fabsf(tao) + sqrtf(1 + tao*tao));	// tan
	c = 1 / sqrtf(1 + t*t);	// cos
	 s = c * t;				// sin

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
void _swap_eigens(float e[3], float v[3][3], int i, int j) {
	float temp;
	float buff[3];
	temp = e[i];
	e[i] = e[j];
	e[j] = temp;
	memcpy(buff, v[i], sizeof(buff));
	memcpy(v[i], v[j], sizeof(buff));
	memcpy(v[j], buff, sizeof(buff));
}

void sortEigens(float e[3], float v[3][3])
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
void eigenSystem(long mat_[3][3], float e[3], float v[3][3])
{
	int i, j;
	float mat[3][3];
	float rot[3][3];
	// init
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			mat[i][j] = mat_[i][j];
		}
	}
	eyes(v);
	// loop
	while (offsum(mat) > 5 * EPS) {
		givensRotation(mat, v, rot, 0, 1);
		givensRotation(mat, v, rot, 0, 2);
		givensRotation(mat, v, rot, 1, 2);
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